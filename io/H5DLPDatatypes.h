#ifndef H5DLPDATATYPES_H
#define H5DLPDATATYPES_H

#include <string>
#include <limits>
#include "hdf5.h"

namespace H5DLP {
    /**
     * @brief Constants for invalid values.
     * 
     * These constants are used to represent invalid or uninitialized values
     * in the data structures.
     */
    const float  kINVALID_FLOAT = std::numeric_limits<float>::max();
    const int    kINVALID_INT   = std::numeric_limits<int  >::max();
    const short  kINVALID_SHORT = std::numeric_limits<short>::max();
    const unsigned long kINVALID_SIZE = std::numeric_limits<unsigned long>::max();

    /**
     * @brief Get the HDF5 type for a given C++ type.
     * 
     * This function is specialized for different types to return the
     * corresponding HDF5 type identifier.
     * 
     * @tparam T The C++ type to get the HDF5 type for.
     * @return hid_t The HDF5 type identifier.
     */
    template <typename T>
    hid_t get_h5type();

    /**
     * @brief Structure representing an association as a slice of another array.
     */
    struct Ass {

        unsigned long start; ///< the first index
        unsigned long end;   ///< the last index

        Ass() : start(kINVALID_SIZE), end(kINVALID_SIZE) {}
    };

    /**
     * @brief Get the HDF5 compound type for the Ass structure.
     * 
     * @return hid_t The HDF5 compound type identifier.
     */
    template<>
    inline hid_t get_h5type<Ass>() {
        hid_t compound_type = H5Tcreate(H5T_COMPOUND, sizeof(H5DLP::Ass));
        H5Tinsert(compound_type, "start", HOFFSET(H5DLP::Ass, start), H5T_NATIVE_ULONG);
        H5Tinsert(compound_type, "end",   HOFFSET(H5DLP::Ass, end  ), H5T_NATIVE_ULONG);
        return compound_type;
    };

    /**
     * @brief Structure representing a step in a particle's path.
     */
    struct PStep {
        float x;            ///< X coordinate
        float y;            ///< Y coordinate
        float z;            ///< Z coordinate
        float t;            ///< Time

        float theta;        ///< Theta angle
        float phi;          ///< Phi angle
        float p;            ///< Momentum
        int   pdg;          ///< PDG code
        int   track_id;     ///< Track ID
        int   ancestor_track_id; ///< Track ID of the ancestor (primary) particle

        float de;           ///< Energy deposition
        float dx;           ///< Step length

        int   proc_start;   ///< Process start
        short subproc_start;///< Subprocess start
        int   proc_stop;    ///< Process stop
        short subproc_stop; ///< Subprocess stop

        PStep() : x(kINVALID_FLOAT), y(kINVALID_FLOAT), z(kINVALID_FLOAT),
                  t(kINVALID_FLOAT), theta(kINVALID_FLOAT), phi(kINVALID_FLOAT),
                  p(kINVALID_FLOAT), pdg(kINVALID_INT), track_id(kINVALID_INT),
                  ancestor_track_id(kINVALID_INT), de(kINVALID_FLOAT), dx(kINVALID_FLOAT),
                  proc_start(kINVALID_INT), subproc_start(kINVALID_SHORT),
                  proc_stop(kINVALID_INT), subproc_stop(kINVALID_SHORT) {} ///< Default constructor
    };

    /**
     * @brief Get the HDF5 compound type for the PStep structure.
     * 
     * @return hid_t The HDF5 compound type identifier.
     */
    template<>
    inline hid_t get_h5type<PStep>() {

        hid_t compound_type = H5Tcreate(H5T_COMPOUND, sizeof(H5DLP::PStep));

        // Set individual types
        H5Tinsert(compound_type, "x",                 HOFFSET(H5DLP::PStep, x                ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "y",                 HOFFSET(H5DLP::PStep, y                ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "z",                 HOFFSET(H5DLP::PStep, z                ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "t",                 HOFFSET(H5DLP::PStep, t                ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "theta",             HOFFSET(H5DLP::PStep, theta            ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "phi",               HOFFSET(H5DLP::PStep, phi              ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "p",                 HOFFSET(H5DLP::PStep, p                ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "pdg",               HOFFSET(H5DLP::PStep, pdg              ), H5T_NATIVE_INT  );
        H5Tinsert(compound_type, "track_id",          HOFFSET(H5DLP::PStep, track_id         ), H5T_NATIVE_INT  );
        H5Tinsert(compound_type, "ancestor_track_id", HOFFSET(H5DLP::PStep, ancestor_track_id), H5T_NATIVE_INT  );
        H5Tinsert(compound_type, "de",                HOFFSET(H5DLP::PStep, de               ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "dx",                HOFFSET(H5DLP::PStep, dx               ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "proc_start",        HOFFSET(H5DLP::PStep, proc_start       ), H5T_NATIVE_INT  );
        H5Tinsert(compound_type, "subproc_start",     HOFFSET(H5DLP::PStep, subproc_start    ), H5T_NATIVE_SHORT);
        H5Tinsert(compound_type, "proc_stop",         HOFFSET(H5DLP::PStep, proc_stop        ), H5T_NATIVE_INT  );
        H5Tinsert(compound_type, "subproc_stop",      HOFFSET(H5DLP::PStep, subproc_stop     ), H5T_NATIVE_SHORT);

        return compound_type;
    };

    /**
     * @brief Structure representing a particle.
     */
    struct Particle {
        float x;                ///< X coordinate
        float y;                ///< Y coordinate
        float z;                ///< Z coordinate
        float t;                ///< Time

        float px;               ///< Momentum in X direction
        float py;               ///< Momentum in Y direction
        float pz;               ///< Momentum in Z direction
        float ke;               ///< Kinetic energy

        float end_x;            ///< End X coordinate
        float end_y;            ///< End Y coordinate
        float end_z;            ///< End Z coordinate
        float end_t;            ///< End time
        float end_px;           ///< End momentum in X direction
        float end_py;           ///< End momentum in Y direction
        float end_pz;           ///< End momentum in Z direction
        float end_ke;           ///< End kinetic energy

        int   track_id;         ///< Track ID
        int   parent_track_id;  ///< Parent track ID
        int   ancestor_track_id;///< Ancestor track ID
        int   pdg;              ///< PDG code
        float mass;             ///< Mass

        int   proc_start;            ///< Process type for creation
        int   subproc_start;         ///< Subprocess type for creation
        std::string proc_name_start; ///< Process name for creation
        int   proc_end;              ///< Process type for the end
        int   subproc_end;           ///< Subprocess type for the end
        std::string proc_name_end;   ///< Process name for the end

        Particle() : x(kINVALID_FLOAT), y(kINVALID_FLOAT), z(kINVALID_FLOAT),
                  t(kINVALID_FLOAT), px(kINVALID_FLOAT), py(kINVALID_FLOAT),
                  pz(kINVALID_FLOAT), ke(kINVALID_FLOAT), end_x(kINVALID_FLOAT),
                  end_y(kINVALID_FLOAT), end_z(kINVALID_FLOAT), end_t(kINVALID_FLOAT),
                  end_px(kINVALID_FLOAT), end_py(kINVALID_FLOAT), end_pz(kINVALID_FLOAT),
                  end_ke(kINVALID_FLOAT), track_id(kINVALID_INT), parent_track_id(kINVALID_INT),
                  ancestor_track_id(kINVALID_INT), pdg(kINVALID_INT), mass(kINVALID_FLOAT),
                  proc_start(kINVALID_INT), subproc_start(kINVALID_SHORT),
                  proc_end(kINVALID_INT), subproc_end(kINVALID_SHORT) {} ///< Default constructor
    };

    /**
     * @brief Get the HDF5 compound type for the Particle structure.
     * 
     * @return hid_t The HDF5 compound type identifier.
     */
    template<>
    inline hid_t get_h5type<Particle>() {

        hid_t compound_type = H5Tcreate(H5T_COMPOUND, sizeof(H5DLP::Particle));

        // Set individual types
        H5Tinsert(compound_type, "x",                HOFFSET(H5DLP::Particle, x                ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "y",                HOFFSET(H5DLP::Particle, y                ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "z",                HOFFSET(H5DLP::Particle, z                ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "t",                HOFFSET(H5DLP::Particle, t                ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "px",               HOFFSET(H5DLP::Particle, px               ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "py",               HOFFSET(H5DLP::Particle, py               ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "pz",               HOFFSET(H5DLP::Particle, pz               ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "ke",               HOFFSET(H5DLP::Particle, ke               ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "end_x",            HOFFSET(H5DLP::Particle, end_x            ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "end_y",            HOFFSET(H5DLP::Particle, end_y            ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "end_z",            HOFFSET(H5DLP::Particle, end_z            ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "end_t",            HOFFSET(H5DLP::Particle, end_t            ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "end_px",           HOFFSET(H5DLP::Particle, end_px           ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "end_py",           HOFFSET(H5DLP::Particle, end_py           ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "end_pz",           HOFFSET(H5DLP::Particle, end_pz           ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "end_ke",           HOFFSET(H5DLP::Particle, end_ke           ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "track_id",         HOFFSET(H5DLP::Particle, track_id         ), H5T_NATIVE_INT  );
        H5Tinsert(compound_type, "parent_track_id",  HOFFSET(H5DLP::Particle, parent_track_id  ), H5T_NATIVE_INT  );
        H5Tinsert(compound_type, "ancestor_track_id",HOFFSET(H5DLP::Particle, ancestor_track_id), H5T_NATIVE_INT  );
        H5Tinsert(compound_type, "pdg",              HOFFSET(H5DLP::Particle, pdg              ), H5T_NATIVE_INT  );
        H5Tinsert(compound_type, "mass",             HOFFSET(H5DLP::Particle, mass             ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "proc_start",       HOFFSET(H5DLP::Particle, proc_start       ), H5T_NATIVE_INT  );
        H5Tinsert(compound_type, "subproc_start",    HOFFSET(H5DLP::Particle, subproc_start    ), H5T_NATIVE_INT  );

        hid_t str_type = H5Tcopy(H5T_C_S1);
        H5Tset_size(str_type, H5T_VARIABLE);
        H5Tinsert(compound_type, "proc_name_start", HOFFSET(H5DLP::Particle, proc_name_start), str_type        );
        H5Tinsert(compound_type, "proc_end",        HOFFSET(H5DLP::Particle, proc_end       ), H5T_NATIVE_INT  );
        H5Tinsert(compound_type, "subproc_end",     HOFFSET(H5DLP::Particle, subproc_end    ), H5T_NATIVE_INT  );
        H5Tinsert(compound_type, "proc_name_end",   HOFFSET(H5DLP::Particle, proc_name_end  ), str_type        );

        H5Tclose(str_type);
        return compound_type;
    };


    /**
     * @brief Structure representing a Primary (particle).
     */
    struct Primary {

        float px;               ///< Momentum in X direction
        float py;               ///< Momentum in Y direction
        float pz;               ///< Momentum in Z direction
        float ke;               ///< Kinetic energy

        int   interaction_id;   ///< Interaction number defined by the generator
        int   track_id;         ///< Track ID
        int   pdg;              ///< PDG code
        float mass;             ///< Mass

        std::string name;       ///< Particle name

        Primary() : px(kINVALID_FLOAT), py(kINVALID_FLOAT), pz(kINVALID_FLOAT),
                  ke(kINVALID_FLOAT), interaction_id(kINVALID_INT), track_id(kINVALID_INT),
                  pdg(kINVALID_INT), mass(kINVALID_FLOAT), name("") {} ///< Default constructor
    };

    /**
     * @brief Get the HDF5 compound type for the Primary (particle) structure.
     * 
     * @return hid_t The HDF5 compound type identifier.
     */
    template<>
    inline hid_t get_h5type<Primary>() {

        hid_t compound_type = H5Tcreate(H5T_COMPOUND, sizeof(H5DLP::Primary));

        // Set individual types
        H5Tinsert(compound_type, "px",               HOFFSET(H5DLP::Primary, px               ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "py",               HOFFSET(H5DLP::Primary, py               ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "pz",               HOFFSET(H5DLP::Primary, pz               ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "ke",               HOFFSET(H5DLP::Primary, ke               ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "interaction_id",   HOFFSET(H5DLP::Primary, interaction_id   ), H5T_NATIVE_INT  );
        H5Tinsert(compound_type, "track_id",         HOFFSET(H5DLP::Primary, track_id         ), H5T_NATIVE_INT  );
        H5Tinsert(compound_type, "pdg",              HOFFSET(H5DLP::Primary, pdg              ), H5T_NATIVE_INT  );
        H5Tinsert(compound_type, "mass",             HOFFSET(H5DLP::Primary, mass             ), H5T_NATIVE_FLOAT);
        hid_t str_type = H5Tcopy(H5T_C_S1);
        H5Tset_size(str_type, H5T_VARIABLE);
        H5Tinsert(compound_type, "name",             HOFFSET(H5DLP::Primary, name             ), str_type);
        H5Tclose(str_type);
        return compound_type;
    };


    /**
     * @brief Structure representing a vertex.
     */
    /// This structure is used to store the information about the primary
    /// vertices in the event. 
    struct Vertex {
        float x;            ///< X coordinate
        float y;            ///< Y coordinate
        float z;            ///< Z coordinate
        float t;            ///< Time
        float energy_sum;   ///< Total energy sum
        float ke_sum;       ///< Total kinetic energy sum
        int num_particles;  ///< Number of particles in the vertex

        std::string generator;     ///< Name of the generator
        std::string reaction;      ///< Reaction type
        std::string filename;      ///< Input generator file name
        int   generator_vertex_id; ///< Interaction number defined by the generator
        float xs;                  ///< Cross section
        float diff_xs;             ///< Differential cross section
        float weight;              ///< Weight of the interaction
        float probability;         ///< Probability of the interaction
        short interaction_id;      ///< Unique instance identifier shared with Primary/Particle instances

        Vertex() : x(kINVALID_FLOAT), y(kINVALID_FLOAT), z(kINVALID_FLOAT),
                  t(kINVALID_FLOAT), energy_sum(kINVALID_FLOAT), ke_sum(kINVALID_FLOAT),
                  num_particles(kINVALID_INT), generator(""), reaction(""), filename(""),
                  generator_vertex_id(kINVALID_INT), xs(kINVALID_FLOAT), diff_xs(kINVALID_FLOAT),
                  weight(kINVALID_FLOAT), probability(kINVALID_FLOAT), interaction_id(kINVALID_SHORT) {} ///< Default constructor
    };

    /**
     * @brief Get the HDF5 compound type for the Vertex structure.
     * 
     * @return hid_t The HDF5 compound type identifier.
     */
    template<>
    inline hid_t get_h5type<Vertex>() {

        hid_t compound_type = H5Tcreate(H5T_COMPOUND, sizeof(H5DLP::Vertex));

        // Set individual types
        H5Tinsert(compound_type, "x",                HOFFSET(H5DLP::Vertex, x             ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "y",                HOFFSET(H5DLP::Vertex, y             ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "z",                HOFFSET(H5DLP::Vertex, z             ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "t",                HOFFSET(H5DLP::Vertex, t             ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "energy_sum",       HOFFSET(H5DLP::Vertex, energy_sum    ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "ke_sum",           HOFFSET(H5DLP::Vertex, ke_sum        ), H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "num_particles",    HOFFSET(H5DLP::Vertex, num_particles ), H5T_NATIVE_INT  );

        hid_t str_type = H5Tcopy(H5T_C_S1);
        H5Tset_size(str_type, H5T_VARIABLE);
        H5Tinsert(compound_type, "generator",       HOFFSET(H5DLP::Vertex, generator      ), str_type);
        H5Tinsert(compound_type, "reaction",        HOFFSET(H5DLP::Vertex, reaction       ), str_type);
        H5Tinsert(compound_type, "filename",        HOFFSET(H5DLP::Vertex, filename       ), str_type);
        H5Tinsert(compound_type, "generator_vertex_id", 
                  HOFFSET(H5DLP::Vertex, generator_vertex_id),H5T_NATIVE_INT  );
        H5Tinsert(compound_type, "xs",              HOFFSET(H5DLP::Vertex, xs             ),H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "diff_xs",         HOFFSET(H5DLP::Vertex,diff_xs         ),H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "weight",          HOFFSET(H5DLP::Vertex, weight         ),H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "probability",     HOFFSET(H5DLP::Vertex, probability    ),H5T_NATIVE_FLOAT);
        H5Tinsert(compound_type, "interaction_id",  HOFFSET(H5DLP::Vertex, interaction_id ),H5T_NATIVE_SHORT);
        H5Tclose(str_type);
        return compound_type;
    };

    /**
     * @brief Structure representing an event.
     */
    struct Event {
        int run_id;         ///< Run ID
        int event_id;       ///< Event ID
        int num_vertices;   ///< Number of vertices
        int num_primaries;  ///< Number of primaries
        int num_particles;  ///< Number of particles
        int num_steps;      ///< Number of steps

        Event() : run_id(kINVALID_INT), event_id(kINVALID_INT),
        num_vertices(kINVALID_INT), num_primaries(kINVALID_INT),
        num_particles(kINVALID_INT), num_steps(kINVALID_INT) {} ///< Default constructor
    };

    /**
     * @brief Get the HDF5 compound type for the Event structure.
     * 
     * @return hid_t The HDF5 compound type identifier.
     */
    template<>
    inline hid_t get_h5type<Event>() {

        hid_t compound_type = H5Tcreate(H5T_COMPOUND, sizeof(H5DLP::Event));

        // Set individual types
        H5Tinsert(compound_type, "run_id",         HOFFSET(H5DLP::Event, run_id         ), H5T_NATIVE_INT  );
        H5Tinsert(compound_type, "event_id",       HOFFSET(H5DLP::Event, event_id       ), H5T_NATIVE_INT  );
        H5Tinsert(compound_type, "num_vertices",   HOFFSET(H5DLP::Event, num_vertices   ), H5T_NATIVE_INT  );
        H5Tinsert(compound_type, "num_primaries",  HOFFSET(H5DLP::Event, num_primaries  ), H5T_NATIVE_INT  );
        H5Tinsert(compound_type, "num_particles",  HOFFSET(H5DLP::Event, num_particles  ), H5T_NATIVE_INT  );
        H5Tinsert(compound_type, "num_steps",      HOFFSET(H5DLP::Event, num_steps      ), H5T_NATIVE_INT  );

        return compound_type;
    };
}
#endif