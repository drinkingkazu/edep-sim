#ifndef H5DLPDATASETS_H
#define H5DLPDATASETS_H
#include "H5DLPDatatypes.h"
#include <vector>
#include <exception>

namespace H5DLP {

    class EventDataset {

        public:
            /**
             * @brief Constructor for the EventDataset class.
             * 
             * Initializes the dataset and its properties.
             */
            EventDataset();
            /**
             * @brief Destructor for the EventDataset class.
             * 
             * Closes the HDF5 dataset and releases resources.
             */
            ~EventDataset() { if(_compound_type) { this->Close(); } }
    
            /**
             * @brief Close the dataset and release resources.
             */
            void Close();
    
            /**
             * @brief Prepare the dataset for writing.
             * 
             * @param file The HDF5 file identifier.
             * @param name The name of the dataset.
             */
            void Prepare(hid_t &file, const std::string name = "steps");

            /**
             * @brief Check if the dataset is open.
             * 
             * @return bool True if the dataset is open, false otherwise.
             */
            bool IsOpen() const { return _dataset != 0; }
    
            /**
             * @brief Get the number of events in the dataset.
             * 
             * @return size_t The number of events.
             */
            size_t GetNumEvents() const { return _num_events; }

            /**
             * @brief Set the event data for the dataset.
             * 
             * @param event The event data to set.
             */
            void Set(const H5DLP::Event &event)
            { _event[0] = event;}

            /**
             * @brief Get the event data from the dataset.
             * 
             * @return H5DLP::Event& Reference to the event data.
             */
            H5DLP::Event& Get() 
            { return _event[0]; }

            /**
             * @brief Write the event data to the dataset.
             * 
             * This function writes the event data to the HDF5 dataset.
             */
            void Write();
    
        private:
            H5DLP::Event _event[1];
            hid_t _compound_type, _space, _plist, _dataset;
            size_t _num_events;
        };




    /**
     * @brief Template class to manage HDF5 variable length dataset storage for supported data types.
     * 
     * @tparam T The type of subject data type.
     */
    template<typename T>
    class VLArrayDataset {

    public:
        VLArrayDataset() : _compound_type(0), _plist(0), _dataset(0), _vlen_type(0), _num_events(0) {}
        ~VLArrayDataset() { if(_compound_type) { this->Close(); } }

        /**
         * @brief Close the storage.
         */
        void Close()
        {
            if(_compound_type) {
                H5Tclose(_compound_type); 
                H5Pclose(_plist);
                H5Dclose(_dataset);
            }
            _compound_type = _plist = _dataset = _vlen_type = 0;
            _num_events = 0;
        }

        /**
         * @brief Prepare the storage HDF5 dataset for writing.
         * 
         * @param file The HDF5 file identifier.
         * @param name The name of the dataset.
         */
        void Prepare(hid_t &file, const std::string name)
        {
            if(_compound_type)
                throw std::exception();

            // Define the compound datatype for the memory
            _compound_type = H5DLP::get_h5type<T>();

            _vlen_type = H5Tvlen_create(_compound_type);

            // Create a simple dataspace with initial size 0 and unlimited maximum size
            hsize_t dims[1] = {0};
            hsize_t maxdims[1] = {H5S_UNLIMITED};
            hid_t space = H5Screate_simple(1, dims, maxdims);

            // Create the dataset creation property list and set the chunk size
            _plist = H5Pcreate(H5P_DATASET_CREATE);
            hsize_t chunk_dims[1] = {1};
            H5Pset_chunk(_plist, 1, chunk_dims);

            // Create the dataset
            _dataset = H5Dcreate2(file,
                name.c_str(),
                _vlen_type, 
                space,
                H5P_DEFAULT,
                _plist,
                H5P_DEFAULT);

            _num_events = 0;

            H5Sclose(space);
        }

        /**
         * @brief Check if the dataset is open.
         * 
         * @return bool True if the dataset is open, false otherwise.
         */
        bool IsOpen() const { return _dataset != 0; }

        /**
         * @brief Get the number of events in the dataset.
         * 
         * @return size_t The number of events.
         */
        size_t GetNumEvents() const { return _num_events; }

        /**
         * @brief Add a unit data entry to the dataset.
         * 
         * @param data The unit data to add.
         */
        void Add(const T &data)
        {
            _data_v.push_back(data);
        }

        /**
         * @brief Add a vector of data to the dataset.
         * 
         * @param data The vector of data to add.
         */
        void Add(const std::vector<T> &data_v)
        {
            _data_v.reserve(_data_v.size() + data_v.size());
            _data_v.insert(_data_v.end(), data_v.begin(), data_v.end());
        }

        void Clear() { _data_v.clear();}

        /**
         * @brief Get the vector of data.
         * 
         * @return std::vector<T>& Reference to the vector of data.
         */
        size_t Size() const { return _data_v.size(); }

        /**
         * @brief Get the data at the specified index.
         * 
         * @return Reference to a unit of data at the specified index of the vector.
         */
        T& At(size_t i) { return _data_v.at(i); }

        /**
         * @brief Reserve the vector memory space.
         * 
         * @param n The size of the vector to reserve.
         */
        void Reserve(size_t n) { _data_v.reserve(n); }

        /**
         * @brief Write the dataset to the HDF5 file.
         */
        void Write()
        {
            // Initialize data to write
            _vlen_data[0].len = _data_v.size();
            _vlen_data[0].p = _data_v.data();

            // Extend the dataset
            hsize_t new_dims[1] = {_num_events + 1};
            H5Dset_extent(_dataset, new_dims);

            // Select the hyperslab
            hid_t dataspace = H5Dget_space(_dataset);
            hsize_t start[1] = {_num_events};
            hsize_t count[1] = {1};
            H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, NULL, count, NULL);

            // Create a memory dataspace
            auto memspace = H5Screate_simple(1, count, NULL);

            // Write the data
            H5Dwrite(_dataset, _vlen_type, memspace, dataspace, H5P_DEFAULT, _vlen_data);

            H5Sclose(memspace);
            H5Sclose(dataspace);
            H5Dflush(_dataset);

            _data_v.clear();
            _num_events += 1;
        }

    private:
        hid_t _compound_type, _plist, _dataset, _vlen_type;
        hvl_t _vlen_data[1];
        size_t _num_events;
        std::vector<T> _data_v;

    };
}

#endif