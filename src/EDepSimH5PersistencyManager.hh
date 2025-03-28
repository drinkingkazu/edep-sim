////////////////////////////////////////////////////////////
// 
#ifndef EDepSim_H5PersistencyManager_hh_seen
#define EDepSim_H5PersistencyManager_hh_seen

namespace H5DLP {
    class Particle;
    class PStep;
    class FileStorage;
}

class TFile;
class TTree;
class TGeoManager;

#include "EDepSimPersistencyManager.hh"

/// Provide a HDF5 output for the geant 4 events.  This just takes the summary
/// from EDepSim::PersistencyManager and dumps it as HDF5 datasets.
namespace EDepSim {class H5PersistencyManager;}
class EDepSim::H5PersistencyManager : public EDepSim::PersistencyManager {
public:
    /// Creates a HDF5 persistency manager.  Through the "magic" of
    /// G4VPersistencyManager the ultimate base class, this declared to the G4
    /// persistency management system.  You can only have one active
    /// persistency class at any give moment.
    H5PersistencyManager();
    virtual ~H5PersistencyManager();

    /// Return true if the ROOT output file is active.  This means that the
    /// output file is open and ready to accept data.
    bool IsOpen();

    /// Stores an event to the output file.
    virtual G4bool Store(const G4Event* anEvent);
    virtual G4bool Store(const G4Run* aRun);
    virtual G4bool Store(const G4VPhysicalVolume* aWorld);

    /// Retrieve information from a file.  These are not implemented.
    virtual G4bool Retrieve(G4Event *&e) {e=NULL; return false;}
    virtual G4bool Retrieve(G4Run* &r) {r=NULL; return false;}
    virtual G4bool Retrieve(G4VPhysicalVolume* &w) {w=NULL; return false;}

    /// Interface with PersistencyMessenger (open and close the
    /// database).
    virtual G4bool Open(G4String dbname);
    virtual G4bool Close(void);

private:
    /// Make the MC Header and add it to truth.
    void MakeMCHeader(const G4Event* src);

private:

    /// The number of events saved to the output file since the last write.
    int fEventsNotSaved;

};
#endif