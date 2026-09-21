#ifndef EDepSim_TrajectoryMap_hh_seen
#define EDepSim_TrajectoryMap_hh_seen
////////////////////////////////////////////////////////////
// $Id: EDepSim::TrajectoryMap.hh,v 1.1 2007/01/01 05:32:49 mcgrew Exp $
//

#include <map>

class G4VTrajectory;

/// Maintain a singleton map of track Id to the trajectory in the trajectory
/// container.  THIS IS NOT THREAD SAFE AND CAN NOT BE USED WITH MULTITHREAD
/// GEANT4.
namespace EDepSim {class TrajectoryMap;}
class EDepSim::TrajectoryMap {
public:
    ~TrajectoryMap() {}

    /// Provide a map between the track id and the trajectory object.
    static G4VTrajectory* Get(int trackId);

    /// Add a trajectory to the map.
    static void Add(G4VTrajectory* traj);

    /// Clear the trajectory map.  This must be done in the
    /// EDepSim::UserEventAction::BeginOfEventAction() method.
    static void Clear();

    /// Find the primary track ID for the current track.  This is the primary
    /// that is the ultimate parent of the current track.  A decay product
    /// counts as a primary, since it should be independently reconstructed,
    /// so this stops at the first decay vertex going up the chain.  Use
    /// FindAncestorId() if you want the true root of the parent chain.
    static int FindPrimaryId(int trackId);

    /// Find the ancestor track ID for the current track.  This is the top of
    /// the parent chain (the track descending from the generator level), and
    /// unlike FindPrimaryId() it is not reset at decay vertices.  Tracks that
    /// are not being saved are skipped, so the result is the highest saved
    /// track in the chain.
    static int FindAncestorId(int trackId);

private:
    /// A map to the trajectories information indexed the the track id. Be
    /// careful since the trajectory information is owned by the event, so if
    /// you try to use this after a trajectory has been deleted... bad things
    /// will happen.
    static std::map<int,G4VTrajectory*> fMap;

    /// The constructor is private so that it can only be created using the
    /// static get method.
    TrajectoryMap() {}
};
#endif
