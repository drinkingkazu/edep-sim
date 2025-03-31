#include "H5DLPFileStorage.h"
#include <stdexcept>
namespace H5DLP {

    FileStorage::FileStorage() : _file(0), _num_events(0)
    {}

    void FileStorage::OpenFile(std::string fname)
    {
        if(_file)
            throw std::runtime_error("file is already open");
        
        _file = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

        // Create the default set of groups
        CreateGroup("ass");
        CreateGroup("event");
    }

    template <typename DatasetType>
    void FileStorage::CloseDataset() {
        std::cout<<"    Closing dataset type: " 
                 << H5DLP::get_name<DatasetType>() 
                 << " (empty datasets will be deleted)" 
                 << std::endl;
        auto& dataset_map = GetDatasetMap<DatasetType>();
        size_t ctr_empty=0;
        size_t ctr_saved=0;
        for(auto &it : dataset_map) {
            for(auto &it2 : it.second) {
                if(it2.second.GetNumWritten() < 1) {
                    // delete empty datasets
                    hid_t target = _file;
                    if(!it.first.empty()) target = _gmap[it.first];
                    if(H5Ldelete(target, it2.first.c_str(), H5P_DEFAULT) < 0) {
                        std::cerr << "    Failed to delete empty dataset " 
                                  << it.first << "/" << it2.first 
                                  << " from file." << std::endl;
                    }
                    ctr_empty++; // Count empty datasets
                }else{
                    std::cout << "        Saved dataset " 
                                  << it.first << "/" << it2.first 
                                  << std::endl;
                    ctr_saved++; // Count saved datasets
                }
                it2.second.Close(); // Close each dataset
            }
        }
        std::cout << "        Saved " << ctr_saved << ", deleted " 
        << ctr_empty << " empty datasets." << std::endl;
    }
        

    void FileStorage::CloseFile()
    {   
        if(!_file)
            return;
        std::cout<<"Closing HDF5 file... (empty datasets will be deleted and not stored)"<<std::endl;
        CloseDataset<H5DLP::Particle>();
        CloseDataset<H5DLP::PStep>();
        CloseDataset<H5DLP::Primary>();
        CloseDataset<H5DLP::Vertex>();
        for(auto &it : _dmap_ass) {
            it.second.Close();
        }
        if(_event.IsOpen())
            _event.Close();

        for(auto &it : _gmap) {
            H5Gclose(it.second);
        }
        if(_file)
            H5Fclose(_file);

        _file = 0;
        _num_events=0;
        _gmap.clear();
        _dmap_name.clear();
        _dmap_ass.clear();
        _dmap_step.clear();
        _dmap_part.clear();
        _dmap_prim.clear();
        _dmap_vertex.clear();
    }

    bool FileStorage::ExistsGroup(const std::string &group_name)
    {
        auto it = _gmap.find(group_name);
        if(it == _gmap.end()) {
            return false;
        }
        return true;
    }

    bool FileStorage::ExistsDataset(const std::string &name, std::string group_name)
    {
        auto it = _dmap_name.find(group_name);
        if(it == _dmap_name.end()) {
            return false;
        }
        auto it2 = it->second.find(name);
        if(it2 == it->second.end()) {
            return false;
        }
        return true;
    }

    template <typename DatasetType>
    std::map<std::string, std::map<std::string, H5DLP::VLArrayDataset<DatasetType> > >&
    FileStorage::GetDatasetMap() {
        static_assert(sizeof(DatasetType) == 0, "Unsupported dataset type");
    }
    
    // Specialization for Particle
    template <>
    std::map<std::string, std::map<std::string, H5DLP::VLArrayDataset<H5DLP::Particle> > >& 
    FileStorage::GetDatasetMap<H5DLP::Particle>() {
        return _dmap_part;
    }
    
    // Specialization for PStep
    template <>
    std::map<std::string, std::map<std::string, H5DLP::VLArrayDataset<H5DLP::PStep> > >& 
    FileStorage::GetDatasetMap<H5DLP::PStep>() {
        return _dmap_step;
    }

    // Specialization for Primary
    template <>
    std::map<std::string, std::map<std::string, H5DLP::VLArrayDataset<H5DLP::Primary> > >&
    FileStorage::GetDatasetMap<H5DLP::Primary>() {
        return _dmap_prim;
    }
    // Specialization for Vertex
    template <>
    std::map<std::string, std::map<std::string, H5DLP::VLArrayDataset<H5DLP::Vertex> > >&
    FileStorage::GetDatasetMap<H5DLP::Vertex>() {
        return _dmap_vertex;
    }

    template <typename DatasetType>
    H5DLP::VLArrayDataset<DatasetType>& FileStorage::GetVLArray(const std::string &name, std::string group_name) {
        auto& dataset_map = GetDatasetMap<DatasetType>();
    
        if (!ExistsGroup(group_name)) {
            if (_num_events) {
                throw std::runtime_error("Error: File already has data, cannot create new group.");
            }
            this->CreateGroup(group_name);
        }
    
        // Check if the dataset already exists
        auto it = _dmap_name[group_name].find(name);
        if (it != _dmap_name[group_name].end()) {
            // The dataset already exists. Check if this is the right type.
            auto iit = dataset_map.find(group_name);
            if (iit != dataset_map.end()) {
                auto iit2 = iit->second.find(name);
                if (iit2 != iit->second.end()) {
                    return iit2->second;
                } else {
                    std::cerr << "Error: Dataset group " << group_name.c_str()
                    << " and name " << name.c_str() << " already exists but with different type."
                    << std::endl;
                    throw std::runtime_error("dataset not found");
                }
            } else {
                std::cerr << "Error: Dataset group " << group_name.c_str()
                          << " and name " << name.c_str() << " already exists but with different type."
                          << std::endl;
                throw std::runtime_error("dataset not found");
            }
        }
    
        // Create the dataset if it doesn't exist
        if (_num_events) {
            std::cerr << "Error: File already has data, cannot create new dataset." << std::endl;
            throw std::runtime_error("cannot create new dataset instance");
        }
    
        auto &ds = dataset_map[group_name][name];
        if (group_name.empty()) {
            ds.Prepare(_file, name);
        } else {
            ds.Prepare(_gmap[group_name], name);
        }
    
        _dmap_name[group_name].insert(name);
        return ds;
    }


    H5DLP::Event& FileStorage::GetEvent(const std::string &name)
    {
        static const std::string group_name="event";

        // Check if the dataset already exists
        auto it = _dmap_name[group_name].find(name);
        if (it != _dmap_name[group_name].end()) {
            if(_event.IsOpen()) {
                return _event.Get();
            }else{
                std::cerr << "Error: the group/name " << group_name.c_str() << "/" << name.c_str()
                << "already exists for a different type dataset." << std::endl;
                throw std::runtime_error("event dataset not found");
            }
        }

        if(_num_events) {
            std::cerr << "Error: File already has data, cannot create new event." << std::endl;
            throw std::runtime_error("cannot create new event instance");
        }
        _event.Prepare(_gmap["event"], name);
        _dmap_name[group_name].insert(name);
        return _event.Get();
    }

    H5DLP::VLArrayDataset<H5DLP::Ass>& FileStorage::GetAss(const std::string &name)
    {
        static const std::string group_name="ass";

        if(_dmap_ass.find(name) != _dmap_ass.end())
            return _dmap_ass[name];

        if(_num_events) {
            std::cerr << "Error: File already has data, cannot create new ass." << std::endl;
            throw std::runtime_error("cannot create new association instance");
        }

        auto &ds = _dmap_ass[name];
        ds.Prepare(_gmap["ass"], name);
    
        return ds;
    }


    template <typename DatasetType>
    void FileStorage::StoreDataset() {
        auto& dataset_map = GetDatasetMap<DatasetType>();
        for(auto &it : dataset_map) {
            for(auto &it2 : it.second) {
                if(it2.second.GetNumEvents() != _num_events) {
                    std::cerr << "Error: Number of events mismatch in " << it.first.c_str() 
                    << "/" << it2.first.c_str() << " (expected " << _num_events 
                    << " but found " << it2.second.GetNumEvents() << ")" << std::endl;
                    throw std::logic_error("dataset's event count inconsistency detected");
                }
                if(!it2.second.IsOpen()){
                    std::cerr << "Error: " << it2.first.c_str() <<" dataset is not open." << std::endl;
                    throw std::logic_error("dataset's open state inconsistency detected");
                }
                it2.second.Write();
            }
        }
    }

    void FileStorage::StoreAss() {
        for(auto &it : _dmap_ass) {
            if(it.second.GetNumEvents() != _num_events) {
                std::cerr << "Error: Number of events mismatch in ass/" << it.first.c_str() 
                << " (expected " << _num_events << " but found " << it.second.GetNumEvents() 
                << ")" << std::endl;
                throw std::runtime_error("association's event count inconsistency detected");
            }
            if(!it.second.IsOpen()){
                std::cerr << "Error: ass/" << it.first.c_str() <<" dataset is not open." << std::endl;
                throw std::runtime_error("association's open state inconsistency detected");
            }
            it.second.Write();
        }
    }

    void FileStorage::Write() {

        StoreDataset<H5DLP::Particle>();
        StoreDataset<H5DLP::PStep>();
        StoreDataset<H5DLP::Primary>();
        StoreDataset<H5DLP::Vertex>();
        StoreAss();

        if(_event.IsOpen()) _event.Write();
        
        _num_events++;
    }

    void FileStorage::CreateGroup(const std::string &group_name)
    {
        if(!_file) {
            std::cerr << "Error: File not opened." << std::endl;
            throw std::runtime_error("failed to create a new group");
        }
        if(_num_events) {
            std::cerr << "Error: File already has data, cannot create new group." << std::endl;
            throw std::runtime_error("failed to create a new group");
        }
        if(group_name.empty()) {
            std::cerr << "Error: Group name is empty." << std::endl;
            throw std::runtime_error("failed to create a new group");
        }
        auto it = _gmap.find(group_name);
        if(it == _gmap.end()) {
            hid_t group = H5Gcreate2(_file, group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            _gmap[group_name] = group;
            _dmap_name[group_name] = std::set<std::string>();
        }
    }
}

namespace H5DLP {
    // Explicit instantiations for GetVLArray
    template H5DLP::VLArrayDataset<H5DLP::Primary>& FileStorage::GetVLArray<H5DLP::Primary>(const std::string&, std::string);
    template H5DLP::VLArrayDataset<H5DLP::Particle>& FileStorage::GetVLArray<H5DLP::Particle>(const std::string&, std::string);
    template H5DLP::VLArrayDataset<H5DLP::PStep>& FileStorage::GetVLArray<H5DLP::PStep>(const std::string&, std::string);
    template H5DLP::VLArrayDataset<H5DLP::Vertex>& FileStorage::GetVLArray<H5DLP::Vertex>(const std::string&, std::string);
}