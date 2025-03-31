//
// Implementation for concrete G4PersistencyManager.
//

#include <TPRegexp.h>

#include "EDepSimPersistencyManager.hh"
#include "EDepSimPersistencyMessenger.hh"
#include "EDepSimVertexInfo.hh"
#include "EDepSimTrajectory.hh"
#include "EDepSimTrajectoryPoint.hh"
#include "EDepSimTrajectoryMap.hh"
#include "EDepSimHitSegment.hh"
#include "EDepSimException.hh"
#include "EDepSimUserRunAction.hh"
#include "EDepSimLog.hh"

#include <G4ios.hh>
#include <G4RunManager.hh>
#include <G4Event.hh>
#include <G4Run.hh>
#include <G4PrimaryVertex.hh>
#include <G4PrimaryParticle.hh>
#include <G4StepStatus.hh>
#include <G4ProcessType.hh>
#include <G4EmProcessSubType.hh>
#include <G4UnitsTable.hh>
#include <G4ParticleTable.hh>
#include <G4SDManager.hh>
#include <G4HCtable.hh>

#include <G4SystemOfUnits.hh>
#include <G4PhysicalConstants.hh>

#include <memory>
#include <chrono>

// Handle foraging the edep-sim information into convenience classes that are
// independent of geant4 and edep-sim internal dependencies.  The classes are
// suitable for saving directly to root (and will be by default if this class
// isn't inherited).  This class can be inherited, and the convenience objects
// can be accessed using specific getters.  If this class is used directly,
// then it will create a simple tree that can be analyzed using root.  The
// best way to access the tree is with the root interface to python (no
// headers needed), or by using ROOT TFile::MakeProject.
EDepSim::PersistencyManager::PersistencyManager()
    : G4VPersistencyManager(), fFilename("/dev/null"),
      fLengthThreshold(10*mm),
      fGammaThreshold(5*MeV), fNeutronThreshold(50*MeV),
      fTrajectoryPointAccuracy(1.*mm), fTrajectoryPointDeposit(0*MeV),
      fSaveAllPrimaryTrajectories(true),
      fROOTOutput(false), fHDF5Output(false),
      fMergeHitSegments(false){
    fPersistencyMessenger = new EDepSim::PersistencyMessenger(this);
}

// This is called by the G4RunManager during AnalyzeEvent.
EDepSim::PersistencyManager::~PersistencyManager() {
    ClearTrajectoryBoundaries();
    ClearTrajectoryBulks();
    delete fPersistencyMessenger;
}

G4bool EDepSim::PersistencyManager::Open(G4String filename) {
    EDepSimSevere(" -- Open is not implimented for " << filename);
    SetFilename(filename);
    return false;
}

/// Make sure the output file is closed.
G4bool EDepSim::PersistencyManager::Close(void) {
    EDepSimSevere(" -- Close is not implimented.");
    return false;
}

// This is called by the G4RunManager during AnalyzeEvent.
G4bool EDepSim::PersistencyManager::Store(const G4Event* anEvent) {
    UpdateSummaries(anEvent);
    return false;
}

// This is called by the G4RunManager during AnalyzeEvent.
G4bool EDepSim::PersistencyManager::Store(const G4Run* aRun) {
    if (!aRun) return false;
    EDepSimSevere(" -- Run store called without a save method for "
               << "GEANT4 run " << aRun->GetRunID());
    return false;
}

// This is called by the G4RunManager during AnalyzeEvent.
G4bool EDepSim::PersistencyManager::Store(const G4VPhysicalVolume* aWorld) {
    if (!aWorld) return false;
    EDepSimSevere(" -- Geometry store called without a save method for "
               << aWorld->GetName());
    return false;
}

void EDepSim::PersistencyManager::AddTrajectoryBoundary(const G4String& b) {
    TPRegexp* bound = new TPRegexp(b.c_str());
    fTrajectoryBoundaries.push_back(bound);
}

void EDepSim::PersistencyManager::AddTrajectoryBulk(const G4String& b) {
    //TPRegexp* bulk = new TPRegexp(b.c_str());
    fTrajectoryBulks.push_back(b);
}

void EDepSim::PersistencyManager::ClearTrajectoryBoundaries() {
    for (std::vector<TPRegexp*>::iterator r = fTrajectoryBoundaries.begin();
         r != fTrajectoryBoundaries.end();
         ++r) {
        delete (*r);
    }
    fTrajectoryBoundaries.clear();
}

void EDepSim::PersistencyManager::ClearTrajectoryBulks() {
    /*for (auto * r fTrajectoryBulks) {
        delete r;
    }*/
    fTrajectoryBulks.clear();
}


//TODO -- Define which Boundaries to save on?
bool EDepSim::PersistencyManager::SaveTrajectoryBoundary(G4VTrajectory* g4Traj,
                                                    G4StepStatus status,
                                                    G4String currentVolume,
                                                    G4String prevVolume) {
    if (status != fGeomBoundary) return false;
    std::string particleInfo = ":" + g4Traj->GetParticleName();
    if (std::abs(g4Traj->GetCharge())<0.1) particleInfo += ":neutral";
    else particleInfo += ":charged";
    std::string current = particleInfo + ":" + currentVolume;
    std::string previous = particleInfo + ":" + prevVolume;
    for (std::vector<TPRegexp*>::iterator r = fTrajectoryBoundaries.begin();
         r != fTrajectoryBoundaries.end();
         ++r) {
        // Check if a watched volume is being entered.
        if ((*r)->Match(current)>0 && (*r)->Match(previous)<1) {
            EDepSimNamedDebug("boundary","Entering " << current);
            return true;
        }
        // Check if a watched volume is being exited.
        if ((*r)->Match(current)<1 && (*r)->Match(previous)>0) {
            EDepSimNamedDebug("boundary","Exiting " << current);
            return true;
        }
    }

    //New: Check Bulk
    /*for (std::vector<TPRegexp*>::iterator r = fTrajectoryBulks.begin();
         r != fTrajectoryBulks.end();
         ++r) {*/
    /*for (auto * r : fTrajectoryBulks) {
         //std::cout << "Checking " << r->GetPattern() << " " << current << std::endl;
        // Check if the particle is traveling through a relevant volume.
        if (r->Match(current)>0 && r->Match(previous)<1) {
            EDepSimNamedDebug("boundary","Bulk " << current);
            std::cout << "Saving " << current << std::endl;
            return true;
        }
    }*/
    return false;
}

bool EDepSim::PersistencyManager::SaveTrajectoryBulk(G4VTrajectory* g4Traj,
                                                    //G4StepStatus status,
                                                    G4String material) {
    //if (status != fGeomBoundary) return false;
    std::string particleInfo = ":" + g4Traj->GetParticleName();
    if (std::abs(g4Traj->GetCharge())<0.1) particleInfo += ":neutral";
    else particleInfo += ":charged";
    std::string current = particleInfo + ":" + material;
    for (const auto r : fTrajectoryBulks) {
         //std::cout << "Checking " << r->GetPattern() << " " << current << " " << r->Match(current) <<std::endl;
         //std::cout << "Checking " << r << " " << current <<std::endl;
        // Check if a watched volume is being traveled through.
        //if (r->Match(current)>0) {
        if (r == current) {
            //std::cout << "\tMatched" << std::endl;
            EDepSimNamedDebug("bulk","bulk " << current);
            return true;
        }
    }

    return false;
}

static bool EventHasHits(const G4Event* event)
{
    G4HCofThisEvent* hitCollections = event->GetHCofThisEvent();
    if (!hitCollections) return false;
    for (int i=0; i < hitCollections->GetNumberOfCollections(); ++i) {
        G4VHitsCollection* g4Hits = hitCollections->GetHC(i);
        if (g4Hits->GetSize() > 0)
            return true;
    }
    return false;
}

bool EDepSim::PersistencyManager::UpdateSummaries(const G4Event* event) {

    const G4Run* runInfo = G4RunManager::GetRunManager()->GetCurrentRun();

    // Summarize the trajectories first so that fTrackIdMap is filled.
    EDepSimDebug("   "<<__PRETTY_FUNCTION__);

    MarkTrajectories(event);

    BuildIndexMap(event);

    if(fROOTOutput) {
        EDepSimDebug("   Constructing the outputs for ROOT");

        fEventSummary.RunId = runInfo->GetRunID();
        fEventSummary.EventId = event->GetEventID();
        EDepSimLog("Event Summary for run " << fEventSummary.RunId
                   << " event " << fEventSummary.EventId);
    
        if (GetRequireEventsWithHits() && not EventHasHits(event)) {
            EDepSimLog("   No hits and /edep/db/set/requireEventsWithHits is true");
            return false;
        }

        SummarizePrimaries(fEventSummary.Primaries,event->GetPrimaryVertex());
        EDepSimDebug("      Primaries " << fEventSummary.Primaries.size());

        SummarizeTrajectories(fEventSummary.Trajectories,event);
        EDepSimDebug("      Trajectories " << fEventSummary.Trajectories.size());

        SummarizeSegmentDetectors(fEventSummary.SegmentDetectors, event);
        EDepSimDebug("      Segment Detectors "
                << fEventSummary.SegmentDetectors.size());

    }else{
        EDepSimDebug("   Constructing the outputs for HDF5");

        if(GetRequireEventsWithHits()) {
            EDepSimError("   HDF5 mode does not support /edep/db/set/requireEventsWithHits");
            return false;
        }

        auto& header = fH5Summary.GetEvent("geant4");
        header.run_id = runInfo->GetRunID();
        header.event_id = event->GetEventID();

        auto& dest_prim   = fH5Summary.GetVLArray<H5DLP::Primary>("geant4","primary");
        auto& dest_vertex = fH5Summary.GetVLArray<H5DLP::Vertex>("geant4","vertex");
        SummarizePrimariesH5(dest_prim, dest_vertex, event->GetPrimaryVertex());
        EDepSimDebug("      Primaries " << dest_prim.Size());

        header.num_vertices = dest_vertex.Size();
        header.num_primaries = dest_prim.Size();

        auto& dest_part = fH5Summary.GetVLArray<H5DLP::Particle>("geant4","particle");
        SummarizeTrajectoriesH5(dest_part, event);
        EDepSimDebug("      Trajectories " << dest_part.Size());

        header.num_particles = dest_part.Size();

        SummarizeSegmentDetectorsH5(fH5Summary, header, event);
    }
    EDepSimDebug("   "<<__PRETTY_FUNCTION__<<" done");

    return true;
}


void EDepSim::PersistencyManager::SummarizePrimariesH5(
    H5DLP::VLArrayDataset<H5DLP::Primary>& dest_part,
    H5DLP::VLArrayDataset<H5DLP::Vertex>& dest,
    const G4PrimaryVertex* src) {
    EDepSimDebug("   "<<__PRETTY_FUNCTION__);

    if (!src) return;
    auto start = std::chrono::high_resolution_clock::now();

    int interaction_id = 0;

    while (src) {

        H5DLP::Vertex vtx;

        vtx.interaction_id = interaction_id;
        vtx.ke_sum = 0;
        vtx.energy_sum = 0;
        vtx.num_particles = 0;
        vtx.x = src->GetX0();
        vtx.y = src->GetY0();
        vtx.z = src->GetZ0();
        vtx.t = src->GetT0();

        // Add the particles associated with the vertex to the summary.
        for (int i=0; i< src->GetNumberOfParticle(); ++i) {

            H5DLP::Primary prim;
            prim.mass = 0.;

            G4PrimaryParticle *g4Prim = src->GetPrimary(i);
            if (g4Prim->GetG4code()) {
                prim.name = g4Prim->GetG4code()->GetParticleName();
                prim.mass = g4Prim->GetG4code()->GetPDGMass();
            }
            prim.pdg = g4Prim->GetPDGcode();
            prim.track_id = g4Prim->GetTrackID() - 1;
            prim.px=g4Prim->GetPx();
            prim.py=g4Prim->GetPy();
            prim.pz=g4Prim->GetPz();
            prim.ke= sqrt(pow(prim.px,2) + pow(prim.py,2) + pow(prim.pz,2));
            prim.interaction_id = interaction_id;

            dest_part.Add(prim);

            vtx.num_particles++;
            vtx.ke_sum += prim.ke;
            vtx.energy_sum += sqrt(pow(prim.ke,2) + pow(prim.mass,2));
        }

        // Check to see if there is anyu user information associated with the
        // vertex.
        EDepSim::VertexInfo* srcInfo
            = dynamic_cast<EDepSim::VertexInfo*>(src->GetUserInformation());
        if (srcInfo) {
            
            vtx.generator = srcInfo->GetName();
            vtx.reaction  = srcInfo->GetReaction();
            vtx.filename  = srcInfo->GetFilename();
            vtx.generator_vertex_id = srcInfo->GetInteractionNumber();
            vtx.xs        = srcInfo->GetCrossSection();
            vtx.diff_xs   = srcInfo->GetDiffCrossSection();
            vtx.weight    = srcInfo->GetWeight();
            vtx.probability = srcInfo->GetProbability();
            const G4PrimaryVertex* infoVertex
                = srcInfo->GetInformationalVertex();
            if (infoVertex) {

                auto& info_dest_part = fH5Summary.GetVLArray<H5DLP::Primary>(vtx.generator.c_str(),"primary");
                auto& info_dest = fH5Summary.GetVLArray<H5DLP::Vertex>(vtx.generator.c_str(),"vertex");
                SummarizePrimariesH5(info_dest_part,info_dest,srcInfo->GetInformationalVertex());

                for(size_t i=0; i<info_dest_part.Size(); ++i) {
                    info_dest_part.At(i).interaction_id = interaction_id;
                }
                for(size_t i=0; i<info_dest.Size(); ++i) {
                    info_dest.At(i).interaction_id = interaction_id;
                }
                EDepSimWarn("InformationalVertex found but not saving in H5 ");
            }
        }

        interaction_id++;
        dest.Add(vtx);
        src = src->GetNext();
    }
    auto end = std::chrono::high_resolution_clock::now();
    EDepSimLog("   Time report - SummarizePrimariesH5: " 
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() 
        << " [ms]");
}

void EDepSim::PersistencyManager::SummarizePrimaries(
    std::vector<TG4PrimaryVertex>& dest,
    const G4PrimaryVertex* src) {
    dest.clear();
    EDepSimDebug("   "<<__PRETTY_FUNCTION__);

    if (!src) return;
    auto start = std::chrono::high_resolution_clock::now();

    while (src) {
        TG4PrimaryVertex vtx;

        vtx.Position.SetX(src->GetX0());
        vtx.Position.SetY(src->GetY0());
        vtx.Position.SetZ(src->GetZ0());
        vtx.Position.SetT(src->GetT0());

        // Add the particles associated with the vertex to the summary.
        for (int i=0; i< src->GetNumberOfParticle(); ++i) {
            TG4PrimaryParticle prim;
            G4PrimaryParticle *g4Prim = src->GetPrimary(i);
            double E = 0.0;
            if (g4Prim->GetG4code()) {
                prim.Name = g4Prim->GetG4code()->GetParticleName();
                E = pow(g4Prim->GetG4code()->GetPDGMass(),2);
            }
            else {
                prim.Name = "unknown";
            }
            prim.PDGCode = g4Prim->GetPDGcode();
            prim.TrackId = g4Prim->GetTrackID() - 1;
            prim.Momentum.SetX(g4Prim->GetPx());
            prim.Momentum.SetY(g4Prim->GetPy());
            prim.Momentum.SetZ(g4Prim->GetPz());
            E += pow(prim.Momentum.P(),2);
            if (E>0) E = std::sqrt(E);
            else E = 0;
            prim.Momentum.SetE(E);
            vtx.Particles.push_back(prim);
        }

        // Check to see if there is anyu user information associated with the
        // vertex.
        EDepSim::VertexInfo* srcInfo
            = dynamic_cast<EDepSim::VertexInfo*>(src->GetUserInformation());
        if (srcInfo) {
            vtx.GeneratorName = srcInfo->GetName();
            vtx.Reaction = srcInfo->GetReaction();
            vtx.Filename = srcInfo->GetFilename();
            vtx.InteractionNumber = srcInfo->GetInteractionNumber();
            vtx.CrossSection = srcInfo->GetCrossSection();
            vtx.DiffCrossSection = srcInfo->GetDiffCrossSection();
            vtx.Weight = srcInfo->GetWeight();
            vtx.Probability = srcInfo->GetProbability();

            const G4PrimaryVertex* infoVertex
                = srcInfo->GetInformationalVertex();
            if (infoVertex) {
                SummarizePrimaries(vtx.Informational,
                                   srcInfo->GetInformationalVertex());
            }
        }

        dest.push_back(vtx);
        src = src->GetNext();
    }
    auto end = std::chrono::high_resolution_clock::now();
    EDepSimLog("   Time report - SummarizePrimaries: " 
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() 
        << " [ms]");
}

void EDepSim::PersistencyManager::SummarizeTrajectories(
    TG4TrajectoryContainer& dest,
    const G4Event* event) {
    dest.clear();
    EDepSimDebug("   "<<__PRETTY_FUNCTION__);

    auto start = std::chrono::high_resolution_clock::now();

    //MarkTrajectories(event);

    const G4TrajectoryContainer* trajectories = event->GetTrajectoryContainer();
    if (!trajectories) {
        EDepSimError("No trajectory container");
        return;
    }
    if (!trajectories->GetVector()) {
        EDepSimError("No trajectory vector in trajectory container");
        return;
    }

    int max_track_id=-1;
    TG4TrajectoryContainer tempContainer;
    for (TrajectoryVector::iterator t = trajectories->GetVector()->begin();
         t != trajectories->GetVector()->end();
         ++t) {
        EDepSim::Trajectory* ndTraj = dynamic_cast<EDepSim::Trajectory*>(*t);

        int track_id = ndTraj->GetTrackID();

        // Check if the trajectory should be saved.
        if(fTrack2OutputIndex.size() <= track_id) 
            continue;
        if(fTrack2OutputIndex[track_id]==-1)
            continue;
        
        if(max_track_id < track_id)
            max_track_id = track_id;

        // Set the particle type information.
        G4ParticleDefinition* part
            = G4ParticleTable::GetParticleTable()->FindParticle(
                ndTraj->GetParticleName());
        if (!part) {
            EDepSimError(std::string("EDepSim::RootPersistencyManager::")
                      + "No particle information for "
                      + ndTraj->GetParticleName());
        }

        TG4Trajectory traj;
        traj.TrackId = track_id;
        traj.Name = ndTraj->GetParticleName();
        traj.PDGCode = ndTraj->GetPDGEncoding();
        traj.ParentId = ndTraj->GetParentID();
        traj.InitialMomentum.SetXYZM(ndTraj->GetInitialMomentum().x(),
                                     ndTraj->GetInitialMomentum().y(),
                                     ndTraj->GetInitialMomentum().z(),
                                     part->GetPDGMass());
        CopyTrajectoryPoints(traj,ndTraj);
        int throttle = 999999;
        do {
            if (traj.ParentId == 0) break;
            EDepSim::Trajectory* pTraj
                = dynamic_cast<EDepSim::Trajectory*>(
                    EDepSim::TrajectoryMap::Get(traj.ParentId));
            if (!pTraj) {
                EDepSimError("Trajectory " << traj.ParentId << " does not exist");
                throw;
                break;
            }
            if (pTraj->SaveTrajectory()) break;
            traj.ParentId = pTraj->GetParentID();
        } while (--throttle > 0);
        tempContainer.push_back(traj);
    }

    // Reorder the trajectories and store in dest container
    dest.resize(tempContainer.size());
    for(size_t i=0; i<tempContainer.size(); ++i) {
        auto& traj = tempContainer[i];
        int output_index = fTrack2OutputIndex[traj.TrackId];
        dest[output_index] = traj;
    }

    /// Rewrite the track ids so that they are consecutive.
    for (TG4TrajectoryContainer::iterator
             t = dest.begin();
         t != dest.end(); ++t) {
        t->TrackId  = fTrack2OutputIndex[t->TrackId];
        t->ParentId = fTrack2OutputIndex[t->ParentId];
    }
    auto end = std::chrono::high_resolution_clock::now();
    EDepSimLog("   Time report - SummarizeTrajectories: " 
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() 
        << " [ms]");
}

void EDepSim::PersistencyManager::SummarizeTrajectoriesH5(
    H5DLP::VLArrayDataset<H5DLP::Particle>& dest,
    const G4Event* event) {

    EDepSimDebug("   "<<__PRETTY_FUNCTION__);

    auto start = std::chrono::high_resolution_clock::now();

    const G4TrajectoryContainer* trajectories = event->GetTrajectoryContainer();
    if (!trajectories) {
        EDepSimError("No trajectory container");
        return;
    }
    if (!trajectories->GetVector()) {
        EDepSimError("No trajectory vector in trajectory container");
        return;
    }

    int max_track_id=-1;
    std::vector<H5DLP::Particle> array1, array2;
    array1.reserve(trajectories->GetVector()->size());
    array2.reserve(trajectories->GetVector()->size());
    for (TrajectoryVector::iterator t = trajectories->GetVector()->begin();
         t != trajectories->GetVector()->end();
         ++t) {
        EDepSim::Trajectory* ndTraj = dynamic_cast<EDepSim::Trajectory*>(*t);

        // Check if the trajectory should be saved.
        int track_id = ndTraj->GetTrackID();
        if(fTrack2OutputIndex.size() <= track_id) 
        {   /*
            EDepSim::TrajectoryPoint* edepPoint
            = dynamic_cast<EDepSim::TrajectoryPoint*>(ndTraj->GetPoint(0));
            G4String prevVolumeName = edepPoint->GetPhysVolName();
            std::cout<<"Found a track ID " << track_id << " > " << fTrack2OutputIndex.size() 
                     << " (size of fTrack2OutputIndex). This means that the track id was not mapped properly. " 
                     << std::endl;
            std::cout<<"  PDG "<<ndTraj->GetPDGEncoding()
            << " Px " << ndTraj->GetInitialMomentum().x()
            << " Py " << ndTraj->GetInitialMomentum().y()
            << " Pz " << ndTraj->GetInitialMomentum().z()
            << " volume: "  << prevVolumeName
            << std::endl;
            */
            continue;
        }
        if(fTrack2OutputIndex[track_id]==-1)
            continue;
        
        if(max_track_id < track_id)
            max_track_id = track_id;

        // Set the particle type information.
        G4ParticleDefinition* g4part
            = G4ParticleTable::GetParticleTable()->FindParticle(
                ndTraj->GetParticleName());
        if (!g4part) {
            EDepSimError(std::string("EDepSim::RootPersistencyManager::")
                      + "No particle information for "
                      + ndTraj->GetParticleName());
        }

        H5DLP::Particle part;

        part.track_id = track_id;
        part.mass = g4part->GetPDGMass();
        part.pdg = ndTraj->GetPDGEncoding();
        part.parent_track_id = ndTraj->GetParentID();
        part.ancestor_track_id = EDepSim::TrajectoryMap::FindPrimaryId(part.track_id);
        part.px = ndTraj->GetInitialMomentum().x();
        part.py = ndTraj->GetInitialMomentum().y();
        part.pz = ndTraj->GetInitialMomentum().z();
        part.ke = std::sqrt(pow(part.px,2) + pow(part.py,2) + pow(part.pz,2));

        // Get the first point
        auto initPoint = dynamic_cast<EDepSim::TrajectoryPoint*>(ndTraj->GetPoint(0));
        auto lastPoint = dynamic_cast<EDepSim::TrajectoryPoint*>(ndTraj->GetPoint(ndTraj->GetPointEntries()-1));

        part.x = initPoint->GetPosition().x();
        part.y = initPoint->GetPosition().y();
        part.z = initPoint->GetPosition().z();
        part.t = initPoint->GetTime();
        part.end_x = lastPoint->GetPosition().x();
        part.end_y = lastPoint->GetPosition().y();
        part.end_z = lastPoint->GetPosition().z();
        part.end_t = lastPoint->GetTime();

        part.end_px = lastPoint->GetMomentum().x();
        part.end_py = lastPoint->GetMomentum().y();
        part.end_pz = lastPoint->GetMomentum().z();
        part.end_ke = std::sqrt(pow(part.end_px,2) + pow(part.end_py,2) + pow(part.end_pz,2));
        part.proc_start = initPoint->GetProcessType();
        part.subproc_start = initPoint->GetProcessSubType();
        part.proc_name_start = initPoint->GetProcessName();
        part.proc_end = lastPoint->GetProcessType();
        part.subproc_end = lastPoint->GetProcessSubType();
        part.proc_name_end = lastPoint->GetProcessName();

        int throttle = 999999;
        do {
            if (part.parent_track_id == 0) break;
            EDepSim::Trajectory* pTraj
                = dynamic_cast<EDepSim::Trajectory*>(
                    EDepSim::TrajectoryMap::Get(part.parent_track_id));
            if (!pTraj) {
                EDepSimError("Trajectory " << part.parent_track_id << " does not exist");
                throw;
                break;
            }
            if (pTraj->SaveTrajectory()) break;
            part.parent_track_id = pTraj->GetParentID();
        } while (--throttle > 0);
        array1.push_back(part);
    }

    // Reorder the trajectories and store in dest container
    array2.resize(array1.size());
    for(size_t i=0; i<array1.size(); ++i) {
        auto& part = array1[i];
        int output_index = fTrack2OutputIndex.at(part.track_id);
        array2.at(output_index) = part;
    }

    /// Rewrite the track ids so that they are consecutive.
    for (auto& part : array2) {
        part.track_id = fTrack2OutputIndex.at(part.track_id);
        part.parent_track_id = fTrack2OutputIndex.at(part.parent_track_id);
        part.ancestor_track_id = fTrack2OutputIndex.at(part.ancestor_track_id);
        dest.Add(part);
    }

    auto end = std::chrono::high_resolution_clock::now();
    EDepSimLog("   Time report - SummarizeTrajectoriesH5: " 
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() 
        << " [ms]");

}

void EDepSim::PersistencyManager::BuildIndexMap(const G4Event* event) {
    EDepSimDebug("   "<<__PRETTY_FUNCTION__);

    auto start = std::chrono::high_resolution_clock::now();
    std::fill(fTrack2OutputIndex.begin(), fTrack2OutputIndex.end(), -1);
    std::fill(fTrack2InputIndex.begin(),  fTrack2InputIndex.end(),  -1);

    const G4TrajectoryContainer* trajectories = event->GetTrajectoryContainer();
    if (!trajectories) {
        EDepSimError("No trajectory container");
        return;
    }
    if (!trajectories->GetVector()) {
        EDepSimError("No trajectory vector in trajectory container");
        return;
    }

    size_t array_size = std::max(fTrack2OutputIndex.size(),trajectories->GetVector()->size());
    fTrack2OutputIndex.resize(array_size, -1);
    fTrack2InputIndex.resize (array_size, -1);

    int input_index = 0;
    int output_count = 0;
    int max_track_id = -1;
    for (TrajectoryVector::iterator t = trajectories->GetVector()->begin();
         t != trajectories->GetVector()->end();
         ++t) {
        EDepSim::Trajectory* ndTraj = dynamic_cast<EDepSim::Trajectory*>(*t);

        // Check if the trajectory should be saved.
        if (!ndTraj->SaveTrajectory()) {
            input_index++;
            continue;
        }
        // Set the particle type information.
        G4ParticleDefinition* part
            = G4ParticleTable::GetParticleTable()->FindParticle(
                ndTraj->GetParticleName());
        if (!part) {
            EDepSimError(std::string("EDepSim::RootPersistencyManager::")
                      + "No particle information for "
                      + ndTraj->GetParticleName());
        }
        int track_id = ndTraj->GetTrackID();
        if(track_id > max_track_id)
            max_track_id = track_id;

        if(track_id >= (int)fTrack2OutputIndex.size()) {
            fTrack2OutputIndex.resize(track_id+1, -1);
            fTrack2InputIndex.resize(track_id+1, -1);
        }
        fTrack2InputIndex[track_id] = input_index++;
        output_count++;
    }

    fOutputIndex2Track.resize(std::max((int)(fOutputIndex2Track.size()),output_count),-1);
    std::fill(fOutputIndex2Track.begin(), fOutputIndex2Track.end(), -1);

    int output_index = 0;
    for (int track_id=0; track_id<max_track_id; ++track_id) {
        if (fTrack2InputIndex[track_id] <0) continue;
        fTrack2OutputIndex[track_id] = output_index;
        fOutputIndex2Track[output_index] = track_id;
        output_index++;
        //std::cout<<track_id<<std::endl;
    }
    //std::cout<<"Size: "<<fTrack2OutputIndex.size()<<" "<<fOutputIndex2Track.size()<<" track id: "<<max_track_id<<std::endl;
    auto end = std::chrono::high_resolution_clock::now();
    EDepSimLog("   Time report - BuildIndexMap: " 
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() 
        << " [ms]");
}

void EDepSim::PersistencyManager::MarkTrajectories(const G4Event* event) {
    const G4TrajectoryContainer* trajectories = event->GetTrajectoryContainer();
    if (!trajectories) {
        EDepSimVerbose("No Trajectories ");
        return;
    }
    EDepSimDebug("   "<<__PRETTY_FUNCTION__);
    auto start = std::chrono::high_resolution_clock::now();

    // Go through all of the trajectories and save:
    //
    //   ** Trajectories from primary particles if
    //       1) a daughter deposited energy in a sensitive detector
    //       2) or, SaveAllPrimaryTrajectories() is true
    //
    //   ** Trajectories created by a particle decay if
    //       1) a daughter deposited energy in a sensitve detector
    //       2) or, SaveAllPrimaryTrajectories() is true.
    //
    //   ** Charged particle trajectories that pass through a sensitive
    //         detector.
    //
    //   ** Gamma-rays and neutrons above a threshold which have daughters
    //         depositing energy in a sensitive detector.
    //
    for (TrajectoryVector::iterator t = trajectories->GetVector()->begin();
         t != trajectories->GetVector()->end();
         ++t) {
        EDepSim::Trajectory* ndTraj = dynamic_cast<EDepSim::Trajectory*>(*t);
        std::string particleName = ndTraj->GetParticleName();
        std::string processName = ndTraj->GetProcessName();
        double initialMomentum = ndTraj->GetInitialMomentum().mag();

        // Check if all primary particle trajectories should be saved.  The
        // primary particle should always be saved if it, or any of it's
        // children, deposited energy in a sensitive detector.  If the primary
        // didn't deposit energy in a sensitive detector, then it will only be
        // saved if SaveAllPrimaryTrajectories is true.
        if (ndTraj->GetParentID() == 0) {
            if (ndTraj->GetSDTotalEnergyDeposit()>1*eV
                || GetSaveAllPrimaryTrajectories()) {
                ndTraj->MarkTrajectory();
                continue;
            }
        }

        int pid = ndTraj->GetParentID();
        std::string parentName = "";
        while (pid != 0) {
          EDepSim::Trajectory * pTraj = dynamic_cast<EDepSim::Trajectory*>(
              EDepSim::TrajectoryMap::Get(pid));
          parentName = pTraj->GetParticleName();
          //std::cout << pid << " " << parentName << std::endl;
          pid = pTraj->GetParentID();
        }
        //std::cout << "progenitor " << parentName << std::endl;

        bool progenitor_is_nu = (std::find(fNuNames.begin(),
                                           fNuNames.end(),
                                           parentName) != fNuNames.end());
        bool is_had = (std::find(fHadNames.begin(),
                                 fHadNames.end(),
                                 particleName) != fHadNames.end());
        if (progenitor_is_nu && is_had) {
          //std::cout << "Saving progeny " << particleName << std::endl;
	  ndTraj->MarkTrajectory(false);
        }


        // Don't save the neutrinos
        if (particleName == "anti_nu_e") continue;
        if (particleName == "anti_nu_mu") continue;
        if (particleName == "anti_nu_tau") continue;
        if (particleName == "nu_e") continue;
        if (particleName == "nu_mu") continue;
        if (particleName == "nu_tau") continue;

        // Save any decay product if it caused any energy deposit.
        if (processName == "Decay") {
            if (ndTraj->GetSDTotalEnergyDeposit()>1*eV
                || GetSaveAllPrimaryTrajectories()) {
                ndTraj->MarkTrajectory(false);
                continue;
            }
        }

        // Save particles that produce charged track inside a sensitive
        // detector.  This doesn't automatically save, but the parents will be
        // automatically considered for saving by the next bit of code.
        if (ndTraj->GetSDLength() > GetLengthThreshold()) {
            ndTraj->MarkTrajectory(false);
            continue;
        }

        // For the next part, only consider particles where the children have
        // deposited energy in a sensitive detector.
        if (ndTraj->GetSDTotalEnergyDeposit()<1*eV) {
            continue;
        }

        // Save higher energy gamma rays that have descendents depositing
        // energy in a sensitive detector.  This only affects secondary
        // photons since primary photons are handled above.
        if (particleName == "gamma" && initialMomentum > GetGammaThreshold()) {
            ndTraj->MarkTrajectory(false);
            continue;
        }

        // Save higher energy neutrons that have descendents depositing energy
        // in a sensitive detector.  This only affects secondary neutrons
        // since primary neutrons are controlled above.
        if (particleName == "neutron"
            && initialMomentum > GetNeutronThreshold()) {
            ndTraj->MarkTrajectory(false);
            continue;
        }

        // Save all pi0s that pass the GetSDTotalEnergyDeposit cut above
        if (particleName == "pi0") {
            ndTraj->MarkTrajectory(false);
            continue;
        }
    }

    // Go through all of the event hit collections and make sure that all
    // primary trajectories and trajectories contributing to a hit are saved.
    // These are mostly a sub-set of the trajectories marked in the previous
    // step, but there are a few corner cases where trajectories are not saved
    // because of theshold issues.
    G4HCofThisEvent* hitCollections = event->GetHCofThisEvent();
    if (!hitCollections) return;
    for (int i=0; i < hitCollections->GetNumberOfCollections(); ++i) {
        G4VHitsCollection* g4Hits = hitCollections->GetHC(i);
        if (g4Hits->GetSize()<1) continue;
        for (unsigned int h=0; h<g4Hits->GetSize(); ++h) {
            EDepSim::HitSegment* g4Hit
                = dynamic_cast<EDepSim::HitSegment*>(g4Hits->GetHit(h));
            if (!g4Hit) {
                EDepSimError("Not a hit segment");
                continue;
            }

            // Explicitly save the primary.  It will probably be marked again
            // with the contributors, but that's OK.  This catches some corner
            // cases where the primary isn't what you would expect.
            int primaryId = g4Hit->GetPrimaryTrajectoryId();
            EDepSim::Trajectory* ndTraj
                = dynamic_cast<EDepSim::Trajectory*>(
                    EDepSim::TrajectoryMap::Get(primaryId));
            if (ndTraj) {
                ndTraj->MarkTrajectory(false);
            }
            else {
                EDepSimError("Primary trajectory not found");
            }

            // Make sure that all the contributors associated with this hit
            // are saved.
            for (int j = 0; j < g4Hit->GetContributorCount(); ++j) {
                int contribId = g4Hit->GetContributor(j);
                EDepSim::Trajectory* contribTraj
                    = dynamic_cast<EDepSim::Trajectory*>(
                        EDepSim::TrajectoryMap::Get(contribId));
                if (contribTraj) {
                    contribTraj->MarkTrajectory(false);
                }
                else {
                    EDepSimError("Contributor trajectory not found");
                }
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    EDepSimLog("   Time report - MarkTrajectories: " 
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() 
        << " [ms]");
}

void EDepSim::PersistencyManager::CopyTrajectoryPoints(TG4Trajectory& traj,
                                                  G4VTrajectory* g4Traj) {
    std::vector<int> selected;

    // Choose the trajectory points that are going to be saved.
    SelectTrajectoryPoints(selected, g4Traj);

    // Make sure the selected trajectory points are in order and unique.
    std::sort(selected.begin(),selected.end());
    selected.erase(std::unique(selected.begin(), selected.end()),
                   selected.end());
    //std::cout << "Selected: " << selected.size() << std::endl;
    ////////////////////////////////////
    // Save the trajectories.
    ////////////////////////////////////
    for (std::vector<int>::iterator tp = selected.begin();
         tp != selected.end(); ++tp) {
        EDepSim::TrajectoryPoint* edepPoint
            = dynamic_cast<EDepSim::TrajectoryPoint*>(g4Traj->GetPoint(*tp));
        TG4TrajectoryPoint point;
        point.Position.SetXYZT(edepPoint->GetPosition().x(),
                               edepPoint->GetPosition().y(),
                               edepPoint->GetPosition().z(),
                               edepPoint->GetTime());
        point.Momentum.SetXYZ(edepPoint->GetMomentum().x(),
                              edepPoint->GetMomentum().y(),
                              edepPoint->GetMomentum().z());
        point.Process = edepPoint->GetProcessType();
        point.Subprocess = edepPoint->GetProcessSubType();
        point.Material = edepPoint->GetMaterial();
        traj.Points.push_back(point);
    }
}

void
EDepSim::PersistencyManager::SummarizeSegmentDetectors(
    TG4HitSegmentDetectors& dest,
    const G4Event* event) {
    EDepSimDebug("   "<<__PRETTY_FUNCTION__);

    auto start = std::chrono::high_resolution_clock::now();
    dest.clear();

    G4HCofThisEvent* HCofEvent = event->GetHCofThisEvent();
    if (!HCofEvent) return;
    G4SDManager *sdM = G4SDManager::GetSDMpointer();
    G4HCtable *hcT = sdM->GetHCtable();
    // Copy each of the hit categories into the output event.
    for (int i=0; i<hcT->entries(); ++i) {
        G4String SDname = hcT->GetSDname(i);
        G4String HCname = hcT->GetHCname(i);
        int HCId = sdM->GetCollectionID(SDname+"/"+HCname);
        G4VHitsCollection* g4Hits = HCofEvent->GetHC(HCId);
        if (g4Hits->GetSize()<1) continue;
        EDepSim::HitSegment* hitSeg
            = dynamic_cast<EDepSim::HitSegment*>(g4Hits->GetHit(0));
        if (!hitSeg) continue;
        SummarizeHitSegments(dest[SDname],g4Hits);
    }
    auto end = std::chrono::high_resolution_clock::now();
    EDepSimLog("   Time report - SummarizeSegmentDetectors: " 
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() 
        << " [ms]");
}

void
EDepSim::PersistencyManager::SummarizeSegmentDetectorsH5(
    H5DLP::FileStorage& dest,
    H5DLP::Event& header,
    const G4Event* event) {
    EDepSimDebug("   "<<__PRETTY_FUNCTION__);

    auto start = std::chrono::high_resolution_clock::now();

    header.num_steps = 0;
    auto& part_v = dest.GetVLArray<H5DLP::Particle>("geant4","particle");
    std::string ass_name_base = "particle_pstep_";

    G4HCofThisEvent* HCofEvent = event->GetHCofThisEvent();
    if (!HCofEvent) return;
    G4SDManager *sdM = G4SDManager::GetSDMpointer();
    G4HCtable *hcT = sdM->GetHCtable();
    // Copy each of the hit categories into the output event.
    for (int i=0; i<hcT->entries(); ++i) {
        G4String SDname = hcT->GetSDname(i);
        G4String HCname = hcT->GetHCname(i);
        int HCId = sdM->GetCollectionID(SDname+"/"+HCname);

        std::string ass_name = ass_name_base + SDname.c_str();
        auto& ass_v = dest.GetAss(ass_name);
        auto& dest_sd = dest.GetVLArray<H5DLP::PStep>(SDname.c_str(),"pstep");

        G4VHitsCollection* g4Hits = HCofEvent->GetHC(HCId);
        if (g4Hits->GetSize()<1) continue;
        EDepSim::HitSegment* hitSeg
            = dynamic_cast<EDepSim::HitSegment*>(g4Hits->GetHit(0));
        if (!hitSeg) continue;

        SummarizeHitSegmentsH5(dest_sd,ass_v,part_v,g4Hits);
        EDepSimDebug("      Segment Detector " << SDname.c_str() << " got " << dest_sd.Size() << " steps");
        header.num_steps += dest_sd.Size();
    }
    auto end = std::chrono::high_resolution_clock::now();
    EDepSimLog("   Time report - SummarizeSegmentDetectorsH5: " 
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() 
        << " [ms]");
}

void
EDepSim::PersistencyManager::SummarizeHitSegments(TG4HitSegmentContainer& dest,
    G4VHitsCollection* g4Hits) {
    EDepSimDebug("   "<<__PRETTY_FUNCTION__);

    dest.clear();

    EDepSim::HitSegment* g4Hit = dynamic_cast<EDepSim::HitSegment*>(g4Hits->GetHit(0));
    if (!g4Hit) return;

    for (std::size_t h=0; h<g4Hits->GetSize(); ++h) {
        g4Hit = dynamic_cast<EDepSim::HitSegment*>(g4Hits->GetHit(h));
        TG4HitSegment hit;
        hit.PrimaryId = fTrackIdMap[g4Hit->GetPrimaryTrajectoryId()];
        hit.EnergyDeposit = g4Hit->GetEnergyDeposit();
        hit.SecondaryDeposit = g4Hit->GetSecondaryDeposit();
        hit.TrackLength = g4Hit->GetTrackLength();
        CopyHitContributors(hit.Contrib,g4Hit->GetContributors());
        hit.Start.SetXYZT(g4Hit->GetStart().x(),
                          g4Hit->GetStart().y(),
                          g4Hit->GetStart().z(),
                          g4Hit->GetStart().t());
        hit.Stop.SetXYZT(g4Hit->GetStop().x(),
                          g4Hit->GetStop().y(),
                          g4Hit->GetStop().z(),
                          g4Hit->GetStop().t());
        hit.StartMomentum.SetXYZT(g4Hit->GetStartMomentum().x(),
            g4Hit->GetStartMomentum().y(),
            g4Hit->GetStartMomentum().z(),
            g4Hit->GetStartMomentum().t());
        hit.StopMomentum.SetXYZT(g4Hit->GetStopMomentum().x(),
            g4Hit->GetStopMomentum().y(),
            g4Hit->GetStopMomentum().z(),
            g4Hit->GetStopMomentum().t());
        hit.StartStatus = g4Hit->GetStartStatus();
        hit.StartProcessType = g4Hit->GetStartProcessType();
        hit.StartProcessSubType = g4Hit->GetStartProcessSubType();
        //hit.StartProcessName = g4Hit->GetStartProcessName();
        hit.StopStatus = g4Hit->GetStopStatus();
        hit.StopProcessType = g4Hit->GetStopProcessType();
        hit.StopProcessSubType = g4Hit->GetStopProcessSubType();
        //hit.StopProcessName = g4Hit->GetStopProcessName();
        dest.push_back(hit);
    }
}

void
EDepSim::PersistencyManager::SummarizeHitSegmentsH5(H5DLP::VLArrayDataset<H5DLP::PStep>& dest,
    H5DLP::VLArrayDataset<H5DLP::Ass>& ass_v,
    H5DLP::VLArrayDataset<H5DLP::Particle>& part_v,
    G4VHitsCollection* g4Hits) {
    EDepSimDebug("   "<<__PRETTY_FUNCTION__);

    EDepSim::HitSegment* g4Hit = dynamic_cast<EDepSim::HitSegment*>(g4Hits->GetHit(0));
    if (!g4Hit) return;

    // Sort the order of an outer array
    std::vector<size_t> reservation;
    reservation.resize(fOutputIndex2Track.size()+1,0);

    for(size_t i=0; i<g4Hits->GetSize(); ++i) {
        g4Hit = dynamic_cast<EDepSim::HitSegment*>(g4Hits->GetHit(i));
        for(auto const& track_id : g4Hit->GetContributors()) {
            if(track_id >= (int)(fTrack2OutputIndex.size()) || fTrack2OutputIndex[track_id] == -1) {
                reservation[fOutputIndex2Track.size()]++;
            }else{
                reservation[fTrack2OutputIndex[track_id]]++;
            }
        }
    }
    
    std::vector<std::vector<H5DLP::PStep> > steps_v;
    steps_v.resize(reservation.size());
    ass_v.Clear();
    H5DLP::Ass ass;
    ass_v.Reserve(reservation.size());
    size_t index_offset=0;
    for(size_t oindex=0; oindex<reservation.size(); ++oindex) {
        if(reservation[oindex] > 0)
            steps_v[oindex].reserve(reservation[oindex]);
        ass.start = index_offset;
        ass.end = index_offset+reservation[oindex];
        index_offset += reservation[oindex];
        ass_v.Add(ass);
    }

    size_t total_array_size=0;
    for (std::size_t h=0; h<g4Hits->GetSize(); ++h) {
        g4Hit = dynamic_cast<EDepSim::HitSegment*>(g4Hits->GetHit(h));
        H5DLP::PStep step;
        //step.primary_id = fTrackIdMap[g4Hit->GetPrimaryTrajectoryId()];
        auto mid_point = g4Hit->GetStart() + (g4Hit->GetStop() - g4Hit->GetStart())/2.;
        step.x = mid_point.x();
        step.y = mid_point.y();
        step.z = mid_point.z();
        step.t = mid_point.t();
        step.theta = g4Hit->GetStartMomentum().theta();
        step.phi = g4Hit->GetStartMomentum().phi();
        step.p = g4Hit->GetStartMomentum().mag();
        step.de = g4Hit->GetEnergyDeposit();
        step.dx = g4Hit->GetTrackLength();
        step.proc_start = g4Hit->GetStartProcessType();
        step.subproc_start = g4Hit->GetStartProcessSubType();
        step.proc_stop = g4Hit->GetStopProcessType();
        step.subproc_stop = g4Hit->GetStopProcessSubType();

        for(size_t cindex=0; cindex < g4Hit->GetContributors().size(); cindex++) {
            step.track_id = g4Hit->GetContributors()[cindex];
            step.de = g4Hit->GetContributions()[cindex];
            size_t output_index;
            if(step.track_id >= (int)(fTrack2OutputIndex.size()) || fTrack2OutputIndex[step.track_id] == -1) {
                output_index = fOutputIndex2Track.size();
                step.track_id = H5DLP::kINVALID_INT;
            }else{
                output_index = fTrack2OutputIndex[step.track_id];
                step.track_id = output_index;
            }

            step.ancestor_track_id = H5DLP::kINVALID_INT;
            step.pdg = H5DLP::kINVALID_INT;
            if(output_index < fOutputIndex2Track.size()) {
                step.ancestor_track_id = part_v.At(output_index).ancestor_track_id;
                step.pdg = part_v.At(output_index).pdg;
            }
            steps_v[output_index].push_back(step);
            total_array_size++;
        }    
    }

    dest.Reserve(total_array_size);
    for(auto const& steps : steps_v) {
        dest.Add(steps);
    }
}

void EDepSim::PersistencyManager::CopyHitContributors(std::vector<int>& dest,
                                                 const std::vector<int>& src) {

    dest.clear();

    for (std::vector<int>::const_iterator c = src.begin();
         c != src.end(); ++c) {
        // Check each contributor to make sure that it is a valid
        // trajectory.  If it isn't in the trajectory map, then set it
        // to a parent that is.
        EDepSim::Trajectory* ndTraj
            = dynamic_cast<EDepSim::Trajectory*>(
                EDepSim::TrajectoryMap::Get(*c));
        while (ndTraj && !ndTraj->SaveTrajectory()) {
            ndTraj = dynamic_cast<EDepSim::Trajectory*>(
                EDepSim::TrajectoryMap::Get(ndTraj->GetParentID()));
        }
        if (!ndTraj) {
            // This trajectory does not exist OR this and none of its parent is to be saved
            dest.push_back(-1 * (*c));
            continue;
        }
        if (fTrackIdMap.find(ndTraj->GetTrackID()) != fTrackIdMap.end()) {
            dest.push_back(fTrackIdMap[ndTraj->GetTrackID()]);
        }
        else {
            // Unexpected: logic error
            EDepSimError("Contributor with unknown trajectory: "
                      << ndTraj->GetTrackID());
            //dest.push_back(-2);
        }
    }

    // Remove the duplicate entries.
    //std::sort(dest.begin(),dest.end());
    //dest.erase(std::unique(dest.begin(),dest.end()),dest.end());
    
    // Duplicate entries should not exist 2025-03-28: assert
    for(std::size_t i=0; i<dest.size(); ++i) {
        for(std::size_t j=i+1; j<dest.size(); ++j) {
            if(dest[i] == dest[j]) {
                EDepSimError("Duplicate contributor: " << dest[i]);
            }
        }
    }

}

double EDepSim::PersistencyManager::FindTrajectoryAccuracy(
    G4VTrajectory* g4Traj, int point1, int point2) {
    if ((point2-point1) < 2) return 0;

    G4ThreeVector p1 = g4Traj->GetPoint(point1)->GetPosition();
    G4ThreeVector p2 = g4Traj->GetPoint(point2)->GetPosition();

    if ( (p2-p1).mag() < fTrajectoryPointAccuracy) return 0;

    G4ThreeVector dir = (p2-p1).unit();

    int step = (point2-point1)/10 + 1;

    double approach = 0.0;
    for (int p = point1+1; p<point2; p = p + step) {
        p2 = g4Traj->GetPoint(p)->GetPosition() - p1;
        approach = std::max((p2 - (dir*p2)*dir).mag(), approach);
    }

    return approach;
}

int EDepSim::PersistencyManager::SplitTrajectory(G4VTrajectory* g4Traj,
                                            int point1, int point2) {

    int point3 = 0.5*(point1 + point2);
    if (point3 <= point1) EDepSimThrow("Points too close to split");
    if (point2 <= point3) EDepSimThrow("Points too close to split");
    double bestAccuracy = FindTrajectoryAccuracy(g4Traj, point1, point3);

    bestAccuracy = std::max(FindTrajectoryAccuracy(g4Traj, point3, point2),
                            bestAccuracy);
    for (int p = point1+1; p<point2-1; ++p) {
        double a1 = FindTrajectoryAccuracy(g4Traj,point1, p);
        double a2 = FindTrajectoryAccuracy(g4Traj,p, point2);
        double accuracy = std::max(a1,a2);
        if (accuracy < bestAccuracy) {
            point3 = p;
            bestAccuracy  = accuracy;
        }
    }

    return point3;
}

void
EDepSim::PersistencyManager::SelectTrajectoryPoints(std::vector<int>& selected,
                                                    G4VTrajectory* g4Traj) {

    selected.clear();
    if (g4Traj->GetPointEntries() < 1) {
        EDepSimError("Trajectory with no points"
                     << " " << g4Traj->GetTrackID()
                     << " " << g4Traj->GetParentID()
                     << " " << g4Traj->GetParticleName());
        return;
    }

    ////////////////////////////////////
    // Save the first point of the trajectory.
    ////////////////////////////////////
    selected.push_back(0);

    /////////////////////////////////////
    // Save the last point of the trajectory.
    /////////////////////////////////////
    int lastIndex = g4Traj->GetPointEntries()-1;
    if (lastIndex < 1) {
        EDepSimError("Trajectory with one point"
                     << " " << g4Traj->GetTrackID()
                     << " " << g4Traj->GetParentID()
                     << " " << g4Traj->GetParticleName());
        return;
    }
    selected.push_back(lastIndex);

    //////////////////////////////////////////////
    // Short out trajectories that don't create any energy deposit in a
    // sensitive detector.  These are trajectories that disappear from the
    // detector, so we don't need to record the extra trajectory points.  The
    // starting and stopping point of the particle are recorded in the
    // trajectory.
    //////////////////////////////////////////////
    EDepSim::Trajectory* ndTraj = dynamic_cast<EDepSim::Trajectory*>(g4Traj);
    //std::cout << "\tSelecting " << lastIndex << " " <<
    //             ndTraj->GetSDTotalEnergyDeposit() << std::endl;
    //if (ndTraj->GetSDTotalEnergyDeposit() < 1*eV) return;

    // Find the trajectory points where particles are entering and leaving the
    // detectors.
    EDepSim::TrajectoryPoint* edepPoint
        = dynamic_cast<EDepSim::TrajectoryPoint*>(g4Traj->GetPoint(0));
    G4String prevVolumeName = edepPoint->GetPhysVolName();
        //std::cout << 0 << " " << edepPoint->GetMaterial() << " " <<
        //             edepPoint->GetPosition().x() << " " <<
        //             edepPoint->GetPosition().y() << " " <<
        //             edepPoint->GetPosition().z() << " " <<
        //             edepPoint->GetProcessType() << " _ _ _" << std::endl;
    //std::cout << 0 << " " << edepPoint->GetMaterial() << " " <<
    //             edepPoint->GetPosition().z() << std::endl;
    bool saved_prev = true;
    for (int tp = 1; tp < lastIndex; ++tp) {
        edepPoint
            = dynamic_cast<EDepSim::TrajectoryPoint*>(g4Traj->GetPoint(tp));
        G4String volumeName = edepPoint->GetPhysVolName();
        // Save the point on a boundary crossing for volumes where we are
        // saving the entry and exit points.
        bool save_bulk = SaveTrajectoryBulk(g4Traj, edepPoint->GetMaterial());
        bool save_bound = SaveTrajectoryBoundary(
            g4Traj, edepPoint->GetStepStatus(), volumeName,prevVolumeName);
        //std::cout << tp << " " << edepPoint->GetMaterial() << " " <<
        //             edepPoint->GetPosition().x() << " " <<
        //             edepPoint->GetPosition().y() << " " <<
        //             edepPoint->GetPosition().z() << " " <<
        //             edepPoint->GetProcessType() << " " << save_bulk << " " <<
        //             save_bound << " " << saved_prev << std::endl;
        //Saving previous point if skipped last time, and this is a good material
        //Last point was a Transportation step
        if (save_bulk && !saved_prev) {
          selected.push_back(tp-1);
        }

        if (save_bound || save_bulk) {
          selected.push_back(tp);
          saved_prev = true;
        }
        else {
          saved_prev = false;
        }
        //TODO -- Add volume where we save all points?
        //        for a given particle -- maybe define by material
        prevVolumeName = volumeName;
    }
    //edepPoint
    //    = dynamic_cast<EDepSim::TrajectoryPoint*>(g4Traj->GetPoint(lastIndex));
    //    std::cout << lastIndex << " " << edepPoint->GetMaterial() << " " <<
    //                 edepPoint->GetPosition().x() << " " <<
    //                 edepPoint->GetPosition().y() << " " <<
    //                 edepPoint->GetPosition().z() << " " <<
    //                 edepPoint->GetProcessType() << " _ _ _" << std::endl;

    // Save trajectory points where there is a "big" interaction.
    for (int tp = 1; tp < lastIndex; ++tp) {
      //std::cout << tp << " processA: " << edepPoint->GetProcessType() << std::endl;
        edepPoint
            = dynamic_cast<EDepSim::TrajectoryPoint*>(g4Traj->GetPoint(tp));
        // Just navigation....
        if (edepPoint->GetProcessType() == fTransportation) continue;
        // Not much energy deposit...
        if (edepPoint->GetProcessDeposit() < GetTrajectoryPointDeposit())
            continue;
        // Don't save optical photons...
        if (edepPoint->GetProcessType() == fOptical) continue;
        // Not a physics step...
        if (edepPoint->GetProcessType() == fGeneral) continue;
        if (edepPoint->GetProcessType() == fUserDefined) continue;
        // Don't save continuous ionization steps.
        if (edepPoint->GetProcessType() == fElectromagnetic
            && edepPoint->GetProcessSubType() == fIonisation) continue;
        // Don't save multiple scattering.
        if (edepPoint->GetProcessType() == fElectromagnetic
            && edepPoint->GetProcessSubType() == fMultipleScattering) continue;
        //std::cout << "Pushing back" << std::endl;
        selected.push_back(tp);
    }

    // Make sure there aren't any duplicates in the selected trajectory points.
    std::sort(selected.begin(), selected.end());
    selected.erase(std::unique(selected.begin(), selected.end()),
                   selected.end());

    if (ndTraj->GetSDEnergyDeposit() < 1*eV) return;

    double desiredAccuracy = GetTrajectoryPointAccuracy();
    // Make sure that the trajectory accuracy stays in tolerance.
    for (int throttle = 0; throttle < 1000; ++throttle) {
        bool addPoint = false;
        for (std::vector<int>::iterator p1 = selected.begin();
             p1 != selected.end();
             ++p1) {
            std::vector<int>::iterator p2 = p1+1;
            if (p2==selected.end()) break;
            double trajectoryAccuracy = FindTrajectoryAccuracy(g4Traj,*p1,*p2);
            if (trajectoryAccuracy <= desiredAccuracy) continue;
            int split = SplitTrajectory(g4Traj,*p1,*p2);
            if (split < 0) continue;
            selected.push_back(split);
            addPoint = true;
            break;
        }
        std::sort(selected.begin(), selected.end());
        selected.erase(std::unique(selected.begin(), selected.end()),
                       selected.end());
        if (!addPoint) break;
    }
}
