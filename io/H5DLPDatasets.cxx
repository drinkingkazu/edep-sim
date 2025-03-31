#ifndef H5DLPDATASETS_CXX
#define H5DLPDATASETS_CXX
#include "H5DLPDatasets.h"
#include <stdexcept>

namespace H5DLP {
    
    EventDataset::EventDataset() : _compound_type(0), 
    _space(0),
    _plist(0),
    _dataset(0),
    _num_events(0) 
    {}

    void EventDataset::Close()
    {
        if(_compound_type) {
            H5Tclose(_compound_type); 
            H5Sclose(_space);
            H5Pclose(_plist);
            H5Dclose(_dataset);
        }
        _compound_type = _space = _plist = _dataset = 0;
        _num_events = 0;
    }

    void EventDataset::Prepare(hid_t &file, const std::string name)
    {
        if(_compound_type)
            throw std::runtime_error("dataset preparation already executed (cannot be done more than once)");

        // Define the compound datatype for the memory
        _compound_type = get_h5type<H5DLP::Event>();

        // Define the dataspace for the dataset.
        hsize_t dims[1] = {0}; // Initial size 
        hsize_t maxdims[1] = {H5S_UNLIMITED}; // Maximum size

        // Create the space
        _space = H5Screate_simple(1, dims, maxdims);
        // Create the property list
        _plist = H5Pcreate(H5P_DATASET_CREATE);
        hsize_t chunk_dims[1] = {1}; // Chunk dimensions
        H5Pset_chunk(_plist, 1, chunk_dims);

        // Create the dataset
        _dataset = H5Dcreate(file,
            name.c_str(),
            _compound_type, 
            _space,
            H5P_DEFAULT,
            _plist,
            H5P_DEFAULT);

        _num_events = 0;
    }

    void EventDataset::Write()
    {    
        hsize_t new_dims[1] = {_num_events+1};
        H5Dset_extent(_dataset, new_dims);

        // Write additional data to extended parts of the dataset
        hsize_t start[1] = {_num_events};
        hsize_t count[1] = {1};
        hid_t filespace = H5Dget_space(_dataset);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, NULL);

        hid_t memspace = H5Screate_simple(1, count, NULL);
        H5Dwrite(_dataset, _compound_type, memspace, filespace, H5P_DEFAULT, _event);

        H5Sclose(memspace);
        H5Sclose(filespace);

        _event[0] = H5DLP::Event(); // Reset the event
        H5Dflush(_dataset);
        _num_events += 1;
    }
}
#endif