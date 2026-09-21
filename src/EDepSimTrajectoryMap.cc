#include "EDepSimTrajectoryMap.hh"
#include "EDepSimTrajectory.hh"
#include "EDepSimException.hh"

#include <G4VTrajectory.hh>
#include <G4VTrajectoryPoint.hh>
#include <G4ThreeVector.hh>

#include <EDepSimLog.hh>

std::map<int, G4VTrajectory*> EDepSim::TrajectoryMap::fMap;

void EDepSim::TrajectoryMap::Clear() {
    fMap.clear();
}

void EDepSim::TrajectoryMap::Add(G4VTrajectory* traj) {
        int trackId = traj->GetTrackID();
        fMap[trackId] = traj;
}

int EDepSim::TrajectoryMap::FindPrimaryId(int trackId) {
    int currentId = trackId;
    int parentId = trackId;
    int loopCount=0;
    for (loopCount=0;loopCount<10000;++loopCount) {
        G4VTrajectory* t = Get(currentId);
        if (!t) break;
        parentId = t->GetParentID();
        // Check to see the search loop should terminate.
        G4VTrajectory* p = Get(parentId);
        // There is no parent so break so this is a primary trajectory.
        if (!p) break;
        // Decay products are primary trajectories since they should be
        // independently reconstructed
        EDepSim::Trajectory * edepTraj = dynamic_cast<EDepSim::Trajectory*>(t);
        if (!edepTraj) EDepSimThrow("Invalid Trajectory");
        if (edepTraj->GetProcessName() == "Decay") break;
        // A parent ID of zero means that this particle is a primary particle,
        // so that makes it a primary trajectory too.
        if (parentId == 0) break;
        currentId = parentId;
    }
    if (loopCount>9999) {
        EDepSimLog("Infinite Loop in EDepSim::TrajectoryMap::FindPrimaryId(): "
                 << "Track Id: " << trackId);
    }
    
    return currentId;
}

int EDepSim::TrajectoryMap::FindAncestorId(int trackId) {
    // Unlike FindPrimaryId(), this walks all the way to the top of the parent
    // chain: a decay vertex does not start a new ancestor.  The returned track
    // is the highest track in the chain that is flagged to be saved, so that
    // the answer stays consistent with the parent ids that are written out
    // (those are also remapped onto the closest saved parent).
    int ancestorId = trackId;
    int currentId = trackId;
    int loopCount=0;
    for (loopCount=0;loopCount<10000;++loopCount) {
        G4VTrajectory* t = Get(currentId);
        // The chain leaves the trajectory map, so the last saved track found
        // is the best ancestor available.
        if (!t) break;
        EDepSim::Trajectory* edepTraj = dynamic_cast<EDepSim::Trajectory*>(t);
        if (!edepTraj) EDepSimThrow("Invalid Trajectory");
        if (edepTraj->SaveTrajectory()) ancestorId = currentId;
        int parentId = edepTraj->GetParentID();
        // A parent id of zero means this is a primary particle, and a self
        // referencing parent would be an infinite loop.
        if (parentId == 0 || parentId == currentId) break;
        currentId = parentId;
    }
    if (loopCount>9999) {
        EDepSimLog("Infinite Loop in EDepSim::TrajectoryMap::FindAncestorId(): "
                 << "Track Id: " << trackId);
    }

    return ancestorId;
}

G4VTrajectory* EDepSim::TrajectoryMap::Get(int trackId) {
    std::map<int,G4VTrajectory*>::iterator t = fMap.find(trackId);
    if (t == fMap.end()) {
        return NULL;
    }
    return t->second;
}

