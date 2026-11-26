#include "H5DLPEventFile.h"
#include <cstdio>
#include <iostream>
#include <stdexcept>

namespace {
    const H5Z_filter_t kBlosc = 32001;
    void Check(hid_t id, const std::string& what) {
        if (id < 0) throw std::runtime_error("H5DLP::EventFile: " + what + " failed");
    }
}

namespace H5DLP {

    hid_t EventFile::Half() {
        hid_t t = H5Tcopy(H5T_IEEE_F32LE);
        H5Tset_fields(t, 15, 10, 5, 0, 10);
        H5Tset_size(t, 2);
        H5Tset_ebias(t, 15);
        return t;
    }

    void EventFile::Attr(hid_t obj, const std::string& name, const std::string& v) {
        hid_t sp = H5Screate(H5S_SCALAR);
        hid_t ty = H5Tcopy(H5T_C_S1);
        H5Tset_size(ty, H5T_VARIABLE);
        H5Tset_cset(ty, H5T_CSET_UTF8);
        hid_t a = H5Acreate2(obj, name.c_str(), ty, sp, H5P_DEFAULT, H5P_DEFAULT);
        Check(a, "attribute " + name);
        const char* text = v.c_str();
        Check(H5Awrite(a, ty, &text), "attribute " + name);
        H5Aclose(a); H5Tclose(ty); H5Sclose(sp);
    }

    void EventFile::Attr(hid_t obj, const std::string& name, double v) {
        hid_t sp = H5Screate(H5S_SCALAR);
        bool integer = (v == (double)(long long) v);
        hid_t a = H5Acreate2(obj, name.c_str(), integer ? H5T_STD_I64LE : H5T_IEEE_F64LE,
                             sp, H5P_DEFAULT, H5P_DEFAULT);
        Check(a, "attribute " + name);
        Check(H5Awrite(a, H5T_NATIVE_DOUBLE, &v), "attribute " + name);
        H5Aclose(a); H5Sclose(sp);
    }

    void EventFile::Open(const std::string& fname, const std::string& role) {
        _file = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        if (_file < 0) { _file = 0; throw std::runtime_error("cannot create " + fname); }
        _num_events = 0;
        _blosc = H5Zfilter_avail(kBlosc) > 0;
        if (!_blosc) std::cout << "H5DLP::EventFile: blosc filter not available, using deflate" << std::endl;
        hid_t s = H5Gcreate2(_file, "schema", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        Check(s, "schema");
        Attr(s, "role", role);
        Attr(s, "format_version", 3.0);
        Attr(s, "compression", std::string(_blosc ? "blosc zstd level 4, byte shuffle"
                                                 : "deflate level 4, byte shuffle"));
        H5Gclose(s);
        _event_group = H5Gcreate2(_file, "event", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        _event = new EventDataset();
        _event->Prepare(_event_group, "geant4");
    }

    void EventFile::Close() {
        if (_entry) { H5Gclose(_entry); _entry = 0; }
        if (_event) { _event->Close(); delete _event; _event = nullptr; }
        if (_event_group) { H5Gclose(_event_group); _event_group = 0; }
        H5Fclose(_file);
        _file = 0;
    }

    void EventFile::BeginEntry() {
        char name[32];
        std::snprintf(name, sizeof(name), "entry_%06zu", _num_events);
        _entry = H5Gcreate2(_file, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        Check(_entry, name);
    }

    void EventFile::EndEntry(const H5DLP::Event& header) {
        Attr(_entry, "run_id", (double) header.run_id);
        Attr(_entry, "event_id", (double) header.event_id);
        H5Gclose(_entry);
        _entry = 0;
        _event->Set(header);
        _event->Write();
        ++_num_events;
    }

    void EventFile::Group(const std::string& name) {
        hid_t g = H5Gcreate2(_entry, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        Check(g, name);
        H5Gclose(g);
    }

    void EventFile::SchemaAttr(const std::string& name, const std::string& value) {
        hid_t s = H5Gopen2(_file, "schema", H5P_DEFAULT); Attr(s, name, value); H5Gclose(s);
    }
    void EventFile::SchemaAttr(const std::string& name, double value) {
        hid_t s = H5Gopen2(_file, "schema", H5P_DEFAULT); Attr(s, name, value); H5Gclose(s);
    }

    void EventFile::StringTable(const std::string& name, const std::vector<std::string>& v) {
        hsize_t n = v.size();
        hid_t sp = H5Screate_simple(1, &n, NULL);
        hid_t ty = H5Tcopy(H5T_C_S1);
        H5Tset_size(ty, H5T_VARIABLE);
        hid_t d = H5Dcreate2(_file, name.c_str(), ty, sp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        Check(d, name);
        std::vector<const char*> p;
        for (const auto& s : v) p.push_back(s.c_str());
        if (n) Check(H5Dwrite(d, ty, H5S_ALL, H5S_ALL, H5P_DEFAULT, p.data()), name);
        H5Dclose(d); H5Tclose(ty); H5Sclose(sp);
    }

    void EventFile::Write(const std::string& name, hid_t mem_type, hid_t file_type,
                          const void* data, size_t n, size_t ncol,
                          const std::string& units, const std::string& description,
                          const NumAttrs& attrs, bool plugin) {
        hsize_t dims[2] = {n, ncol};
        int rank = ncol ? 2 : 1;
        hid_t sp = H5Screate_simple(rank, dims, NULL);
        hid_t pl = H5Pcreate(H5P_DATASET_CREATE);
        if (n) {
            H5Pset_chunk(pl, rank, dims);
            if (_blosc && plugin) {
                unsigned cd[7] = {0, 0, 0, 0, 4, 1, 5};  // level 4, byte shuffle, zstd
                Check(H5Pset_filter(pl, kBlosc, H5Z_FLAG_MANDATORY, 7, cd), "blosc filter");
            }
            else {
                H5Pset_shuffle(pl);
                H5Pset_deflate(pl, 4);
            }
        }
        hid_t d = H5Dcreate2(_entry, name.c_str(), file_type, sp, H5P_DEFAULT, pl, H5P_DEFAULT);
        Check(d, name);
        if (n) Check(H5Dwrite(d, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data), name);
        Attr(d, "units", units);
        Attr(d, "description", description);
        for (const auto& a : attrs) Attr(d, a.first, a.second);
        H5Dclose(d); H5Pclose(pl); H5Sclose(sp);
    }
}
