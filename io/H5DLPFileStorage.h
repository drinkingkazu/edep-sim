#ifndef H5DLPFILESTORAGE_H
#define H5DLPFILESTORAGE_H
#include "H5DLPDatasets.h"
#include <iostream>
#include <map>
#include <set>

namespace H5DLP {

    /**
     * @brief Class to manage HDF5 file operations.
     */
    class FileStorage {

        public:
            FileStorage();
            ~FileStorage() {if(_file) this->CloseFile(); }
    
            /**
             * @brief Open an HDF5 file.
             * 
             * @param fname The name of the file to open.
             */
            void OpenFile(std::string fname);
    
            /**
             * @brief Close the HDF5 file.
             */
            void CloseFile();
    
            /**
             * @brief Check if the file is open.
             */
            bool IsOpen() { return _file != 0; }
    
            /**
             * @brief Check if a group exists.
             * 
             * @param name The name of the group to check.
             */
            bool ExistsGroup(const std::string &group_name);
    
            /**
             * @brief Check if a dataset exists.
             * 
             * @param name The name of the dataset to check.
             * @param group_name The name of the group to check in.
             */
            bool ExistsDataset(const std::string &name, std::string group_name = "");
    
            /**
             * @brief Get the particle storage of DatasetType for a given name.
             * 
             * @param name The name of the dataset.
             * @return Reference to the VLArrayDataset<H5DLP::DatasetType> object.
             */            
            template<typename DatasetType>
            H5DLP::VLArrayDataset<DatasetType>& GetVLArray(const std::string& name, std::string group_name = "");

            /**
             * @brief Get the event dataset.
             * 
             * @param name The name of the event dataset.
             * @return Reference to the event object.
             */
            H5DLP::Event& GetEvent(const std::string& name);

            /**
             * @brief Get the association dataset.
             * 
             * @param name The name of the association dataset.
             * @return Reference to the association dataset.
             */
            H5DLP::VLArrayDataset<H5DLP::Ass>& GetAss(const std::string& name);

            /**
             * @brief Write all datasets to the HDF5 file.
             */
            void Write();
    
        private:
    
            /**
             * @brief Create a new group in the HDF5 file.
             * 
             * @param group_name The name of the group to create.
             */
            void CreateGroup(const std::string &group_name);

            /**
             * @brief Get the dataset map for a given DatasetType.
             * 
             * @tparam DatasetType The type of the dataset.
             * @return Reference to the dataset map.
             */
            template <typename DatasetType>
            std::map<std::string, std::map<std::string, H5DLP::VLArrayDataset<DatasetType> > >& GetDatasetMap();

            /**
             * @brief Store the dataset for a given DatasetType.
             * 
             * @tparam DatasetType The type of the dataset.
             */
            template <typename DatasetType>
            void StoreDataset();

            /**
             * @brief Store the association datasets.
             */
            void StoreAss();
    
        private:
            hid_t _file;
            size_t _num_events;
            std::map<std::string, hid_t> _gmap;
            std::map<std::string, std::set<std::string> > _dmap_name;
            H5DLP::EventDataset _event;
            std::map<std::string, H5DLP::VLArrayDataset<H5DLP::Ass> > _dmap_ass;
            std::map<std::string, std::map<std::string, H5DLP::VLArrayDataset<H5DLP::Particle> > > _dmap_part;
            std::map<std::string, std::map<std::string, H5DLP::VLArrayDataset<H5DLP::PStep   > > > _dmap_step;
            std::map<std::string, std::map<std::string, H5DLP::VLArrayDataset<H5DLP::Primary > > > _dmap_prim;
            std::map<std::string, std::map<std::string, H5DLP::VLArrayDataset<H5DLP::Vertex  > > > _dmap_vertex;
    
    };
    
}
#endif