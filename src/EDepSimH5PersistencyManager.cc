////////////////////////////////////////////////////////////
//

#include "EDepSimLog.hh"

#include "EDepSimH5PersistencyManager.hh"
#include "EDepSimRootGeometryManager.hh"
#include "EDepSimUserPrimaryGeneratorAction.hh"
#include "kinem/EDepSimKinemPassThrough.hh"
#include "kinem/EDepSimPrimaryGenerator.hh"
#include "kinem/EDepSimRooTrackerKinematicsGenerator.hh"
#include "kinem/EDepSimVKinematicsGenerator.hh"

#include <globals.hh>

#include <G4Event.hh>
#include <G4Run.hh>
#include <G4RunManager.hh>

#include <TGeoManager.h>


EDepSim::H5PersistencyManager::H5PersistencyManager() 
    : EDepSim::PersistencyManager() {}

EDepSim::H5PersistencyManager::~H5PersistencyManager() {
    if (fH5Summary.IsOpen()) {
        fH5Summary.CloseFile();
    }
}

bool EDepSim::H5PersistencyManager::IsOpen() { if (fH5Summary.IsOpen()) fH5Summary.CloseFile(); }

bool EDepSim::H5PersistencyManager::Open(G4String filename) {
    if (fH5Summary.IsOpen()) {
        EDepSimError("EDepSim::H5PersistencyManager::Open "
                 << "-- Close current file" );
        return false;
    }
    if (fROOTOutput) {
        EDepSimError("EDepSim::H5PersistencyManager::Open "
                   << "-- Cannot open HDF5 file when fROOTOutput is true");
        return false;
    }
    
    SetFilename(filename);

    EDepSimLog("EDepSim::H5PersistencyManager::Open " << GetFilename());

    fH5Summary.OpenFile(GetFilename());
       
    fEventsNotSaved = 0;

    return true;
}

bool EDepSim::H5PersistencyManager::Close() {
    if (!fH5Summary.IsOpen()) {
        EDepSimError("EDepSim::H5PersistencyManager::Close "
                   << "-- No Output File");
        return false;
    }
    fH5Summary.CloseFile();

    return true;
}

bool EDepSim::H5PersistencyManager::Store(const G4Event* anEvent) {
    if (!fH5Summary.IsOpen()) {
        EDepSimError("EDepSim::H5PersistencyManager::Store "
                   << "-- No Output File");
        return false;
    }

    bool save = UpdateSummaries(anEvent);

    if (save) {
        fH5Summary.Write();

        auto genAction = dynamic_cast<const EDepSim::UserPrimaryGeneratorAction*>(
            G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());

        for (int i = 0; i < genAction->GetGeneratorCount(); ++i) {
            auto gen = dynamic_cast<const EDepSim::PrimaryGenerator*>(
                genAction->GetGenerator(i));
            auto vkinGen = const_cast<EDepSim::VKinematicsGenerator*>(
                gen->GetKinematicsGenerator());
            auto kinGen = dynamic_cast<EDepSim::RooTrackerKinematicsGenerator*>(
                vkinGen);

            if (kinGen) {
                EDepSim::KinemPassThrough::GetInstance()->AddEntry(kinGen);
                break;
            }
        }
    }

    return true;
}

bool EDepSim::H5PersistencyManager::Store(const G4Run*) {
    return false;
}

bool EDepSim::H5PersistencyManager::Store(const G4VPhysicalVolume*) {
    if (!fH5Summary.IsOpen()) {
        EDepSimError("EDepSim::H5PersistencyManager::Store "
                   << "-- No Output File");
        return false;
    }
    if (!gGeoManager) {
        EDepSimError("EDepSim::H5PersistencyManage::Store(world)"
                  << " -- Cannot be run before /edep/update");
        return false; 
    }
    //gGeoManager->Write();
    return true;
}

