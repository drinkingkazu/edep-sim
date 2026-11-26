#ifndef H5DLPEVENTFILE_H
#define H5DLPEVENTFILE_H
#include "hdf5.h"
#include "H5DLPDatasets.h"
#include <string>
#include <utility>
#include <vector>

namespace H5DLP {

    /// One HDF5 file of the step/mctruth layout: one group per event
    /// (entry_NNNNNN) holding plain fixed-type datasets, plus the event table
    /// event/geant4 of the older layout (used to identify the events).
    ///
    /// Datasets are one chunk each, compressed with blosc-zstd when the HDF5
    /// filter plugin 32001 can be loaded (HDF5_PLUGIN_PATH), otherwise with
    /// shuffle + deflate.
    class EventFile {
    public:
        typedef std::vector<std::pair<std::string, double> > NumAttrs;

        EventFile() : _file(0), _entry(0), _num_events(0), _blosc(false) {}
        ~EventFile() { if (_file) Close(); }

        void Open(const std::string& fname, const std::string& role);
        void Close();
        bool IsOpen() const { return _file != 0; }

        void BeginEntry();
        /// Finish the entry and add its row to event/geant4.
        void EndEntry(const H5DLP::Event& header);

        /// Create a group in the current entry.
        void Group(const std::string& name);

        /// Attributes of the file's schema group.
        void SchemaAttr(const std::string& name, const std::string& value);
        void SchemaAttr(const std::string& name, double value);
        /// A root level list of strings (code tables).
        void StringTable(const std::string& name, const std::vector<std::string>& v);

        /// Write n x ncol values (ncol 0: one dimension) to the current entry.
        /// plugin false: deflate, readable without the blosc plugin.
        void Write(const std::string& name, hid_t mem_type, hid_t file_type,
                   const void* data, size_t n, size_t ncol,
                   const std::string& units, const std::string& description,
                   const NumAttrs& attrs = NumAttrs(), bool plugin = true);

        static hid_t Half();  ///< IEEE half precision (the caller closes it)

    private:
        static void Attr(hid_t obj, const std::string& name, const std::string& v);
        static void Attr(hid_t obj, const std::string& name, double v);
        hid_t _file, _entry;
        size_t _num_events;
        bool _blosc;
        H5DLP::EventDataset* _event = nullptr;
        hid_t _event_group = 0;
    };
}
#endif
