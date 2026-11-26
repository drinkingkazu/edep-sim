//
// The step and mctruth files: the HDF5 output with one group per event.
//
// step     what a detector simulation reads: the Geant4 steps and their groups.
// mctruth  everything else edep-sim has: per-step truth (same rows as step),
//          particles, trajectory points, primaries, vertices.
//
// Positions are cell indices k on a fixed grid, stored as row to row
// differences; a position is (k + cell_offset)*scale_factor.  Group run lengths
// are computed from these stored positions, so a reader gets the same groups.
//
#include "EDepSimPersistencyManager.hh"
#include "EDepSimTrajectory.hh"
#include "EDepSimTrajectoryPoint.hh"
#include "EDepSimTrajectoryMap.hh"
#include "EDepSimHitSegment.hh"
#include "EDepSimException.hh"
#include "EDepSimLog.hh"

#include <G4Event.hh>
#include <G4HCtable.hh>
#include <G4Run.hh>
#include <G4RunManager.hh>
#include <G4SDManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4TrajectoryContainer.hh>
#include <G4Version.hh>

#include <cmath>
#include <cstring>
#include <limits>

namespace {
    const double kGrid = 0.01*mm;
    const double kCellOffset = 0.5;

    // Fixed size rows of the tables.
    const size_t kName = 32;
    const size_t kLongName = 256;

    void CopyName(char* dest, size_t size, const std::string& src) {
        if (src.size() >= size) EDepSimError("Name truncated: " << src);
        std::strncpy(dest, src.c_str(), size-1);
    }

    struct ParticleRow {
        float x, y, z, t, px, py, pz, ke;
        float end_x, end_y, end_z, end_t, end_px, end_py, end_pz, end_ke;
        int track_id, parent_track_id, ancestor_track_id, pdg;
        float mass;
        int proc_start, subproc_start, proc_end, subproc_end;
        char proc_name_start[kName], proc_name_end[kName], name[kName];
        int g4_track_id, g4_parent_id, primary_track_id, num_points;
    };
    struct PrimaryRow {
        float px, py, pz, ke;
        int interaction_id, track_id, pdg;
        float mass;
        char name[kName];
    };
    struct VertexRow {
        float x, y, z, t, energy_sum, ke_sum;
        int num_particles;
        char generator[kName], reaction[kLongName], filename[kLongName];
        int generator_vertex_id;
        float xs, diff_xs, weight, probability;
        int interaction_id;
    };

    /// A compound type built member by member.
    struct Compound {
        hid_t id;
        explicit Compound(size_t size) : id(H5Tcreate(H5T_COMPOUND, size)) {}
        ~Compound() { H5Tclose(id); }
        void Add(const char* name, size_t offset, hid_t type) {
            H5Tinsert(id, name, offset, type);
        }
        void Str(const char* name, size_t offset, size_t size) {
            hid_t t = H5Tcopy(H5T_C_S1);
            H5Tset_size(t, size);
            H5Tinsert(id, name, offset, t);
            H5Tclose(t);
        }
    };
#define MEMBER(C, S, M, T) C.Add(#M, HOFFSET(S, M), T)
#define STRING(C, S, M) C.Str(#M, HOFFSET(S, M), sizeof(S::M))

    void WriteTable(H5DLP::EventFile& f, const std::string& name, const Compound& c,
                    const void* data, size_t n, const std::string& description) {
        // the small tables are read by production tools without the plugin
        f.Write(name, c.id, c.id, data, n, 0, "mm, ns, MeV", description,
                H5DLP::EventFile::NumAttrs(), name == "particle");
    }

    void WritePrimaries(H5DLP::EventFile& f, const std::string& name,
                        const std::vector<H5DLP::Primary>& in) {
        std::vector<PrimaryRow> rows(in.size());
        for (size_t i = 0; i < in.size(); ++i) {
            PrimaryRow& r = rows[i];
            std::memset(&r, 0, sizeof(r));
            r.px = in[i].px; r.py = in[i].py; r.pz = in[i].pz; r.ke = in[i].ke;
            r.interaction_id = in[i].interaction_id; r.track_id = in[i].track_id;
            r.pdg = in[i].pdg; r.mass = in[i].mass;
            CopyName(r.name, sizeof(r.name), in[i].name);
        }
        Compound c(sizeof(PrimaryRow));
        MEMBER(c, PrimaryRow, px, H5T_NATIVE_FLOAT); MEMBER(c, PrimaryRow, py, H5T_NATIVE_FLOAT);
        MEMBER(c, PrimaryRow, pz, H5T_NATIVE_FLOAT); MEMBER(c, PrimaryRow, ke, H5T_NATIVE_FLOAT);
        MEMBER(c, PrimaryRow, interaction_id, H5T_NATIVE_INT);
        MEMBER(c, PrimaryRow, track_id, H5T_NATIVE_INT); MEMBER(c, PrimaryRow, pdg, H5T_NATIVE_INT);
        MEMBER(c, PrimaryRow, mass, H5T_NATIVE_FLOAT); STRING(c, PrimaryRow, name);
        WriteTable(f, name, c, rows.data(), rows.size(), "primary particles");
    }

    void WriteVertices(H5DLP::EventFile& f, const std::string& name,
                       const std::vector<H5DLP::Vertex>& in) {
        std::vector<VertexRow> rows(in.size());
        for (size_t i = 0; i < in.size(); ++i) {
            VertexRow& r = rows[i];
            std::memset(&r, 0, sizeof(r));
            r.x = in[i].x; r.y = in[i].y; r.z = in[i].z; r.t = in[i].t;
            r.energy_sum = in[i].energy_sum; r.ke_sum = in[i].ke_sum;
            r.num_particles = in[i].num_particles;
            CopyName(r.generator, sizeof(r.generator), in[i].generator);
            CopyName(r.reaction, sizeof(r.reaction), in[i].reaction);
            CopyName(r.filename, sizeof(r.filename), in[i].filename);
            r.generator_vertex_id = in[i].generator_vertex_id;
            r.xs = in[i].xs; r.diff_xs = in[i].diff_xs;
            r.weight = in[i].weight; r.probability = in[i].probability;
            r.interaction_id = in[i].interaction_id;
        }
        Compound c(sizeof(VertexRow));
        MEMBER(c, VertexRow, x, H5T_NATIVE_FLOAT); MEMBER(c, VertexRow, y, H5T_NATIVE_FLOAT);
        MEMBER(c, VertexRow, z, H5T_NATIVE_FLOAT); MEMBER(c, VertexRow, t, H5T_NATIVE_FLOAT);
        MEMBER(c, VertexRow, energy_sum, H5T_NATIVE_FLOAT);
        MEMBER(c, VertexRow, ke_sum, H5T_NATIVE_FLOAT);
        MEMBER(c, VertexRow, num_particles, H5T_NATIVE_INT);
        STRING(c, VertexRow, generator); STRING(c, VertexRow, reaction);
        STRING(c, VertexRow, filename);
        MEMBER(c, VertexRow, generator_vertex_id, H5T_NATIVE_INT);
        MEMBER(c, VertexRow, xs, H5T_NATIVE_FLOAT); MEMBER(c, VertexRow, diff_xs, H5T_NATIVE_FLOAT);
        MEMBER(c, VertexRow, weight, H5T_NATIVE_FLOAT);
        MEMBER(c, VertexRow, probability, H5T_NATIVE_FLOAT);
        MEMBER(c, VertexRow, interaction_id, H5T_NATIVE_INT);
        WriteTable(f, name, c, rows.data(), rows.size(), "generator vertices");
    }

    int Cell(double x) {
        double k = std::floor(x/kGrid);
        if (std::abs(k) > std::numeric_limits<int>::max()/2)
            EDepSimThrow("Position outside of the range of the position grid");
        return (int) k;
    }

    /// Absolute cells (n x 3) to row to row differences.
    std::vector<int> Differences(const std::vector<int>& cell) {
        std::vector<int> d(cell);
        for (size_t i = cell.size(); i-- > 3; ) d[i] = cell[i] - cell[i-3];
        return d;
    }

    template <class T> struct H5Type;
    template <> struct H5Type<float> { static hid_t Id() { return H5T_NATIVE_FLOAT; } };
    template <> struct H5Type<int> { static hid_t Id() { return H5T_NATIVE_INT; } };
    template <> struct H5Type<short> { static hid_t Id() { return H5T_NATIVE_SHORT; } };
    template <> struct H5Type<unsigned short> { static hid_t Id() { return H5T_NATIVE_USHORT; } };
    template <> struct H5Type<signed char> { static hid_t Id() { return H5T_NATIVE_SCHAR; } };
    template <> struct H5Type<unsigned char> { static hid_t Id() { return H5T_NATIVE_UCHAR; } };

    template <class T>
    void WriteColumn(H5DLP::EventFile& f, const std::string& name, const std::vector<T>& v,
                     const std::string& units, const std::string& description,
                     size_t ncol = 0,
                     const H5DLP::EventFile::NumAttrs& attrs = H5DLP::EventFile::NumAttrs()) {
        f.Write(name, H5Type<T>::Id(), H5Type<T>::Id(), v.data(),
                ncol ? v.size()/ncol : v.size(), ncol, units, description, attrs);
    }
    /// IEEE half precision, rounded to nearest even.  (The conversion in HDF5
    /// 1.10 gives wrong values when the rounding carries into the exponent.)
    unsigned short ToHalf(float value) {
        unsigned int f;
        std::memcpy(&f, &value, sizeof(f));
        unsigned int sign = (f >> 16) & 0x8000u;
        unsigned int a = f & 0x7fffffffu;
        if (a > 0x7f800000u) return sign | 0x7e00u;          // NaN
        if (a >= 0x477ff000u) return sign | 0x7c00u;         // rounds to infinity
        unsigned int man, shift, half;
        if (a >= 0x38800000u) {                              // normal
            half = (a - 0x38000000u) >> 13;
            man = a; shift = 13;
        }
        else {                                               // subnormal or zero
            if (a <= 0x33000000u) return sign;               // at most half the smallest
            shift = 126 - (a >> 23);                         // 14..24
            man = (a & 0x7fffffu) | 0x800000u;
            half = man >> shift;
        }
        unsigned int rest = man & ((1u << shift) - 1u), mid = 1u << (shift - 1u);
        if (rest > mid || (rest == mid && (half & 1u))) ++half;
        return sign | half;
    }

    void WriteHalf(H5DLP::EventFile& f, const std::string& name, const std::vector<float>& v,
                   const std::string& units, const std::string& description) {
        std::vector<unsigned short> h(v.size());
        for (size_t i = 0; i < v.size(); ++i) h[i] = ToHalf(v[i]);
        hid_t half = H5DLP::EventFile::Half();
        f.Write(name, half, half, h.data(), h.size(), 0, units, description);
        H5Tclose(half);
    }
    void WriteCells(H5DLP::EventFile& f, const std::string& name, const std::vector<int>& cell,
                    const std::string& description) {
        H5DLP::EventFile::NumAttrs a;
        a.push_back(std::make_pair("scale_factor", kGrid/mm));
        a.push_back(std::make_pair("cell_offset", kCellOffset));
        WriteColumn(f, name, Differences(cell), "mm",
                    description + " (x,y,z): k = cumsum over rows, "
                    "position = (k + cell_offset)*scale_factor", 3, a);
    }
}

int EDepSim::PersistencyManager::Code(StringTable& table, const std::string& value) {
    for (size_t i = 0; i < table.size(); ++i) if (table[i] == value) return i;
    table.push_back(value);
    return table.size()-1;
}

void EDepSim::PersistencyManager::OpenStepTruth() {
    std::string name = GetFilename();
    size_t dot = name.rfind('.');
    size_t slash = name.rfind('/');
    if (dot == std::string::npos || (slash != std::string::npos && dot < slash)) dot = name.size();
    std::string base = name.substr(0, dot), ext = name.substr(dot);
    fStepFile.Open(base + "_step" + ext, "step");
    fTruthFile.Open(base + "_mctruth" + ext, "mctruth");
    fMaterialNames.clear(); fProcessNames.clear(); fVolumeNames.clear();
    for (H5DLP::EventFile* f : {&fStepFile, &fTruthFile}) {
        f->SchemaAttr("geant4_version", G4Version);
        f->SchemaAttr("units", "mm, ns, MeV, MeV/c, rad");
        f->SchemaAttr("group_size", fGroupSize);
        f->SchemaAttr("group_gap_mm", fGroupGap/mm);
        f->SchemaAttr("length_threshold_mm", GetLengthThreshold()/mm);
        f->SchemaAttr("gamma_threshold_MeV", GetGammaThreshold()/MeV);
        f->SchemaAttr("neutron_threshold_MeV", GetNeutronThreshold()/MeV);
        f->SchemaAttr("trajectory_accuracy_mm", GetTrajectoryPointAccuracy()/mm);
        f->SchemaAttr("trajectory_deposit_MeV", GetTrajectoryPointDeposit()/MeV);
        f->SchemaAttr("save_all_primary_trajectories", GetSaveAllPrimaryTrajectories());
    }
}

void EDepSim::PersistencyManager::CloseStepTruth() {
    if (fStepFile.IsOpen()) fStepFile.Close();
    if (!fTruthFile.IsOpen()) return;
    fTruthFile.StringTable("material_names", fMaterialNames);
    fTruthFile.StringTable("process_names", fProcessNames);
    fTruthFile.StringTable("volume_names", fVolumeNames);
    fTruthFile.Close();
}

bool EDepSim::PersistencyManager::WriteStepTruth(const G4Event* event) {
    H5DLP::Event header;
    header.run_id = G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();
    header.event_id = event->GetEventID();

    fInfoVertices.clear(); fInfoPrimaries.clear();
    H5DLP::VLArrayDataset<H5DLP::Primary> prim;
    H5DLP::VLArrayDataset<H5DLP::Vertex> vertex;
    H5DLP::VLArrayDataset<H5DLP::Particle> part;
    SummarizePrimariesH5(prim, vertex, event->GetPrimaryVertex());
    SummarizeTrajectoriesH5(part, event);
    header.num_vertices = vertex.Size();
    header.num_primaries = prim.Size();
    header.num_particles = part.Size();

    // The hits, one bucket per saved particle.  Segments are not merged, so
    // each is one Geant4 step of one track.
    std::vector<std::vector<EDepSim::HitSegment*> > bucket(part.Size());
    std::string sdName;
    if (G4HCofThisEvent* hce = event->GetHCofThisEvent()) {
        G4SDManager* sdM = G4SDManager::GetSDMpointer();
        G4HCtable* hcT = sdM->GetHCtable();
        for (int i = 0; i < hcT->entries(); ++i) {
            G4String name = hcT->GetSDname(i);
            G4VHitsCollection* hits = hce->GetHC(
                sdM->GetCollectionID(name + "/" + hcT->GetHCname(i)));
            if (!hits || hits->GetSize() < 1) continue;
            if (!dynamic_cast<EDepSim::HitSegment*>(hits->GetHit(0))) continue;
            if (!sdName.empty())
                EDepSimThrow("The step file supports one sensitive detector with hits");
            sdName = name;
            for (size_t h = 0; h < hits->GetSize(); ++h) {
                EDepSim::HitSegment* hit
                    = dynamic_cast<EDepSim::HitSegment*>(hits->GetHit(h));
                if (hit->GetContributors().size() != 1)
                    EDepSimThrow("The step file needs unmerged hit segments");
                int id = hit->GetContributors()[0];
                if (id >= (int) fTrack2OutputIndex.size() || fTrack2OutputIndex[id] < 0)
                    EDepSimThrow("Step of a track that is not saved");
                bucket[fTrack2OutputIndex[id]].push_back(hit);
            }
        }
    }

    std::vector<int> cell, startCell, stopCell, trackId;
    std::vector<float> t, de, dx, theta, phi, p;
    std::vector<float> startT, stopT, stopTheta, stopPhi, stopP, secondary;
    std::vector<short> procStart, procStop, subStart, subStop;
    std::vector<signed char> statStart, statStop;
    std::vector<unsigned short> runs;
    double prev[3] = {0, 0, 0};
    for (size_t track = 0; track < bucket.size(); ++track) {
        bool first = true;
        for (EDepSim::HitSegment* hit : bucket[track]) {
            const G4LorentzVector& a = hit->GetStart();
            const G4LorentzVector& b = hit->GetStop();
            // float first: the same value as in the older layout
            G4LorentzVector mid = a + (b - a)/2.;
            float x[3] = {(float) mid.x(), (float) mid.y(), (float) mid.z()};
            double now[3], d2 = 0;
            for (int k = 0; k < 3; ++k) {
                cell.push_back(Cell(x[k]));
                startCell.push_back(Cell(a.vect()[k]));
                stopCell.push_back(Cell(b.vect()[k]));
                now[k] = (cell.back() + kCellOffset)*kGrid;
                d2 += (now[k] - prev[k])*(now[k] - prev[k]);
                prev[k] = now[k];
            }
            if (first || std::sqrt(d2) > fGroupGap || runs.back() >= fGroupSize)
                runs.push_back(0);
            ++runs.back();
            first = false;
            trackId.push_back(track);
            t.push_back(mid.t());
            de.push_back(hit->GetContributions()[0]);
            dx.push_back(hit->GetTrackLength());
            theta.push_back(hit->GetStartMomentum().theta());
            phi.push_back(hit->GetStartMomentum().phi());
            p.push_back(hit->GetStartMomentum().vect().mag());
            startT.push_back(a.t());
            stopT.push_back(b.t());
            stopTheta.push_back(hit->GetStopMomentum().theta());
            stopPhi.push_back(hit->GetStopMomentum().phi());
            stopP.push_back(hit->GetStopMomentum().vect().mag());
            secondary.push_back(hit->GetSecondaryDeposit());
            procStart.push_back(hit->GetStartProcessType());
            procStop.push_back(hit->GetStopProcessType());
            subStart.push_back(hit->GetStartProcessSubType());
            subStop.push_back(hit->GetStopProcessSubType());
            statStart.push_back(hit->GetStartStatus());
            statStop.push_back(hit->GetStopStatus());
        }
    }
    header.num_steps = trackId.size();

    // ------------------------------------------------------------- step
    H5DLP::EventFile& S = fStepFile;
    S.BeginEntry();
    WriteCells(S, "positions", cell, "middle of the step");
    WriteColumn(S, "t", t, "ns", "time at the middle of the step");
    WriteColumn(S, "de", de, "MeV", "energy deposit");
    WriteHalf(S, "dx", dx, "mm", "step length; 0 for neutral steps stored as a point");
    WriteHalf(S, "theta", theta, "rad", "direction at the start of the step, from the z axis");
    WriteHalf(S, "phi", phi, "rad", "direction at the start of the step");
    WriteColumn(S, "p", p, "MeV/c", "momentum at the start of the step");
    WriteColumn(S, "track_id", trackId, "", "row of the particle table in mctruth");
    H5DLP::EventFile::NumAttrs ga;
    ga.push_back(std::make_pair("id_base", 1.));
    ga.push_back(std::make_pair("group_size", (double) fGroupSize));
    ga.push_back(std::make_pair("gap_threshold_mm", fGroupGap/mm));
    WriteColumn(S, "group_run_lengths", runs, "steps",
                "group k (1 = first) is this many consecutive rows: up to group_size steps "
                "of one track, not more than gap_threshold_mm apart", 0, ga);
    S.EndEntry(header);

    // ------------------------------------------------------------- mctruth
    H5DLP::EventFile& M = fTruthFile;
    M.BeginEntry();
    M.Group("steps");
    WriteCells(M, "steps/start", startCell, "start of the step");
    for (size_t i = 0; i < stopCell.size(); ++i) stopCell[i] -= startCell[i];
    H5DLP::EventFile::NumAttrs sa;
    sa.push_back(std::make_pair("scale_factor", kGrid/mm));
    WriteColumn(M, "steps/stop_minus_start", stopCell, "mm",
                "cell of the stop minus cell of the start (x,y,z), same grid as start", 3, sa);
    WriteColumn(M, "steps/start_t", startT, "ns", "time at the start of the step");
    WriteColumn(M, "steps/stop_t", stopT, "ns", "time at the stop of the step");
    WriteHalf(M, "steps/stop_theta", stopTheta, "rad", "direction at the stop of the step");
    WriteHalf(M, "steps/stop_phi", stopPhi, "rad", "direction at the stop of the step");
    WriteColumn(M, "steps/stop_p", stopP, "MeV/c", "momentum at the stop of the step");
    WriteColumn(M, "steps/secondary_deposit", secondary, "MeV",
                "sum of the non-ionizing deposit set by the Doke-Birks process; not a usable "
                "light yield (equals de for nearly all steps at 500 V/cm)");
    WriteColumn(M, "steps/proc_start", procStart, "", "G4ProcessType before the step");
    WriteColumn(M, "steps/proc_stop", procStop, "", "G4ProcessType that ended the step");
    WriteColumn(M, "steps/subproc_start", subStart, "", "process sub type, -1 none");
    WriteColumn(M, "steps/subproc_stop", subStop, "", "process sub type, -1 none");
    WriteColumn(M, "steps/start_status", statStart, "", "G4StepStatus at the start");
    WriteColumn(M, "steps/stop_status", statStop, "", "G4StepStatus at the stop");

    // Particles and their selected trajectory points.
    std::vector<ParticleRow> rows(part.Size());
    std::vector<float> px, pt, pm, pd;
    std::vector<short> pproc, psub, pname;
    std::vector<signed char> pstat;
    std::vector<unsigned char> pmat, pvol;
    for (size_t i = 0; i < part.Size(); ++i) {
        const H5DLP::Particle& in = part.At(i);
        ParticleRow& r = rows[i];
        std::memset(&r, 0, sizeof(r));
        r.x = in.x; r.y = in.y; r.z = in.z; r.t = in.t;
        r.px = in.px; r.py = in.py; r.pz = in.pz; r.ke = in.ke;
        r.end_x = in.end_x; r.end_y = in.end_y; r.end_z = in.end_z; r.end_t = in.end_t;
        r.end_px = in.end_px; r.end_py = in.end_py; r.end_pz = in.end_pz; r.end_ke = in.end_ke;
        r.track_id = in.track_id; r.parent_track_id = in.parent_track_id;
        r.ancestor_track_id = in.ancestor_track_id; r.pdg = in.pdg; r.mass = in.mass;
        r.proc_start = in.proc_start; r.subproc_start = in.subproc_start;
        r.proc_end = in.proc_end; r.subproc_end = in.subproc_end;
        CopyName(r.proc_name_start, kName, in.proc_name_start);
        CopyName(r.proc_name_end, kName, in.proc_name_end);
        int g4id = fOutputIndex2Track.at(i);
        EDepSim::Trajectory* traj
            = dynamic_cast<EDepSim::Trajectory*>(EDepSim::TrajectoryMap::Get(g4id));
        if (!traj) EDepSimThrow("Particle without a trajectory");
        CopyName(r.name, kName, traj->GetParticleName());
        r.g4_track_id = g4id;
        r.g4_parent_id = traj->GetParentID();
        int primary = EDepSim::TrajectoryMap::FindPrimaryId(g4id);
        r.primary_track_id = primary < (int) fTrack2OutputIndex.size()
            ? fTrack2OutputIndex[primary] : -1;
        std::vector<int> selected;
        SelectTrajectoryPoints(selected, traj);
        r.num_points = selected.size();
        for (int k : selected) {
            const EDepSim::TrajectoryPoint* q
                = dynamic_cast<EDepSim::TrajectoryPoint*>(traj->GetPoint(k));
            for (int c = 0; c < 3; ++c) {
                px.push_back(q->GetPosition()[c]);
                pm.push_back(q->GetMomentum()[c]);
            }
            pt.push_back(q->GetTime());
            pd.push_back(q->GetProcessDeposit());
            pproc.push_back(q->GetProcessType());
            psub.push_back(q->GetProcessSubType());
            pstat.push_back(q->GetStepStatus());
            int material = Code(fMaterialNames, q->GetMaterial());
            int volume = Code(fVolumeNames, q->GetPhysVolName());
            if (material > 255 || volume > 255) EDepSimThrow("Too many materials or volumes");
            pmat.push_back(material);
            pvol.push_back(volume);
            pname.push_back(Code(fProcessNames, q->GetProcessName()));
        }
    }
    Compound c(sizeof(ParticleRow));
    MEMBER(c, ParticleRow, x, H5T_NATIVE_FLOAT); MEMBER(c, ParticleRow, y, H5T_NATIVE_FLOAT);
    MEMBER(c, ParticleRow, z, H5T_NATIVE_FLOAT); MEMBER(c, ParticleRow, t, H5T_NATIVE_FLOAT);
    MEMBER(c, ParticleRow, px, H5T_NATIVE_FLOAT); MEMBER(c, ParticleRow, py, H5T_NATIVE_FLOAT);
    MEMBER(c, ParticleRow, pz, H5T_NATIVE_FLOAT); MEMBER(c, ParticleRow, ke, H5T_NATIVE_FLOAT);
    MEMBER(c, ParticleRow, end_x, H5T_NATIVE_FLOAT); MEMBER(c, ParticleRow, end_y, H5T_NATIVE_FLOAT);
    MEMBER(c, ParticleRow, end_z, H5T_NATIVE_FLOAT); MEMBER(c, ParticleRow, end_t, H5T_NATIVE_FLOAT);
    MEMBER(c, ParticleRow, end_px, H5T_NATIVE_FLOAT); MEMBER(c, ParticleRow, end_py, H5T_NATIVE_FLOAT);
    MEMBER(c, ParticleRow, end_pz, H5T_NATIVE_FLOAT); MEMBER(c, ParticleRow, end_ke, H5T_NATIVE_FLOAT);
    MEMBER(c, ParticleRow, track_id, H5T_NATIVE_INT);
    MEMBER(c, ParticleRow, parent_track_id, H5T_NATIVE_INT);
    MEMBER(c, ParticleRow, ancestor_track_id, H5T_NATIVE_INT);
    MEMBER(c, ParticleRow, pdg, H5T_NATIVE_INT); MEMBER(c, ParticleRow, mass, H5T_NATIVE_FLOAT);
    MEMBER(c, ParticleRow, proc_start, H5T_NATIVE_INT);
    MEMBER(c, ParticleRow, subproc_start, H5T_NATIVE_INT);
    STRING(c, ParticleRow, proc_name_start);
    MEMBER(c, ParticleRow, proc_end, H5T_NATIVE_INT);
    MEMBER(c, ParticleRow, subproc_end, H5T_NATIVE_INT);
    STRING(c, ParticleRow, proc_name_end);
    STRING(c, ParticleRow, name);
    MEMBER(c, ParticleRow, g4_track_id, H5T_NATIVE_INT);
    MEMBER(c, ParticleRow, g4_parent_id, H5T_NATIVE_INT);
    MEMBER(c, ParticleRow, primary_track_id, H5T_NATIVE_INT);
    MEMBER(c, ParticleRow, num_points, H5T_NATIVE_INT);
    WriteTable(M, "particle", c, rows.data(), rows.size(),
               "saved particles; parent_track_id is the closest saved ancestor, g4_* are the "
               "Geant4 ids (parent_track_id was remapped if the parent row's g4_track_id is not "
               "g4_parent_id); primary_track_id stops at decays");
    WritePrimaries(M, "primary", prim.Data());
    WriteVertices(M, "vertex", vertex.Data());
    if (!fInfoVertices.empty()) {
        WriteVertices(M, "info_vertex", fInfoVertices);
        WritePrimaries(M, "info_primary", fInfoPrimaries);
    }
    M.Group("points");
    WriteColumn(M, "points/position", px, "mm",
                "trajectory points, particle by particle (particle.num_points rows each)", 3);
    WriteColumn(M, "points/t", pt, "ns", "time");
    WriteColumn(M, "points/momentum", pm, "MeV/c", "momentum (x,y,z)", 3);
    WriteColumn(M, "points/process_deposit", pd, "MeV", "energy deposit of the step to this point");
    WriteColumn(M, "points/proc", pproc, "", "G4ProcessType");
    WriteColumn(M, "points/subproc", psub, "", "process sub type, -1 none");
    WriteColumn(M, "points/process_name", pname, "", "row of /process_names");
    WriteColumn(M, "points/status", pstat, "", "G4StepStatus");
    WriteColumn(M, "points/material", pmat, "", "row of /material_names");
    WriteColumn(M, "points/volume", pvol, "", "row of /volume_names");
    M.EndEntry(header);
    return true;
}
