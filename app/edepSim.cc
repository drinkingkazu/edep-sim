#include "EDepSimCreateRunManager.hh"
#include "EDepSimPersistencyManager.hh"
#include "EDepSimRootPersistencyManager.hh"
#include "EDepSimH5PersistencyManager.hh"
#include "EDepSimLog.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"

#include "globals.hh"

#include <iostream>
#include <fstream>
#include <csignal>
#include <cstdlib>

void usage () {
    std::cout << "Usage: edep-sim [options] [macros]" << std::endl;
    std::cout << "    -C      -- Toggle validating the geometry" << std::endl;
    std::cout << "    -d      -- Increase the debug level" << std::endl;
    std::cout << "    -D <name>=[error,severe,warn,debug,trace]"
              << std::endl
              << "            -- Change the named debug level"
              << std::endl;
    std::cout << "    -e <n>  -- Add /run/beamOn <n> after last macro."
              << std::endl;
    std::cout << "    -g      -- Set a GDML file" << std::endl;
    std::cout << "    -o      -- Set the output file" << std::endl;
    std::cout << "    -p      -- Select the physics list" << std::endl;
    std::cout << "    -s      -- Set the seed from the time" << std::endl;
    std::cout << "    -u      -- Do update before running the macros"
              << std::endl;
    std::cout << "    -U      -- Start an interactive run after the macros"
              << std::endl
              << "               are processed."
              << std::endl;
    std::cout << "    -v      -- Increase the verbosity" << std::endl;
    std::cout << "    -V <name>=[quiet,log,info,verbose]" << std::endl
              << "            -- Change the named log level"
              << std::endl;
    std::cout << "    -h      -- This help message." << std::endl;

    exit(1);
}

bool hasSuffix(const std::string& str, const std::string& suffix) {
    if (str.length() < suffix.length()) return false;
    return str.compare(str.length() - suffix.length(), suffix.length(), suffix) == 0;
}

int main(int argc,char** argv) {
    EDepSim::LogManager::Configure();

    std::string outputFilename;
    std::string physicsList = "";
    std::string gdmlFilename = "";

    int errflg = 0;
    int c = 0;
    bool useUI = false;
    bool setSeed=false;
    bool doUpdate=false;
    bool validateGeometry=true;
    int debugLevel = 0;
    std::map<std::string, EDepSim::LogManager::ErrorPriority> namedDebugLevel;

    int logLevel = 1; // Will choose default logging level...
    std::map<std::string, EDepSim::LogManager::LogPriority> namedLogLevel;

    // If filled, run this many events.  The command line argument must be an
    // integer or G4 will give an error.  There could be error checking here,
    // but assume the user is "intelligent", and won't try to run "-e ten"
    // instead of "-e 10".
    std::string beamOnCount = "";

    if (argc<2) usage();

    while (!errflg && ((c=getopt(argc,argv,"CdD:e:g:o:p:qsuUvV:h:H")) != -1)) {
        switch (c) {
        case 'C': {
            // Toggle the validateGeometry flag.  The default value is set
            // above.
            validateGeometry = !validateGeometry;
            break;
        }
        case 'd':
        {
            // increase the debugging level.
            ++debugLevel;
            break;
        }
        case 'D':
        {
            // Set the debug level for a named trace.
            std::string arg(optarg);
            std::size_t sep = arg.find("=");
            if (sep != std::string::npos) {
                std::string name = arg.substr(0,sep);
                std::string levelName = arg.substr(sep+1);
                switch (levelName[0]) {
                case 'e': case 'E':
                    namedDebugLevel[name.c_str()]
                        = EDepSim::LogManager::ErrorLevel;
                    break;
                case 's': case 'S':
                    namedDebugLevel[name.c_str()]
                        = EDepSim::LogManager::SevereLevel;
                    break;
                case 'w': case 'W':
                    namedDebugLevel[name.c_str()]
                        = EDepSim::LogManager::WarnLevel;
                    break;
                case 'd': case 'D':
                    namedDebugLevel[name.c_str()]
                        = EDepSim::LogManager::DebugLevel;
                    break;
                case 't': case 'T':
                    namedDebugLevel[name.c_str()]
                        = EDepSim::LogManager::TraceLevel;
                    break;
                default:
                    usage();
                }
            }
            break;
        }
        case 'e': {
            // Add a /run/beamOn command after the last macro has been
            // processed.
            beamOnCount = optarg;
            break;
        }
        case 'g': {
            gdmlFilename = optarg;
            break;
        }
        case 'o': {
            outputFilename = optarg;
            break;
        }
        case 'p': {
            physicsList = optarg;
            break;
        }
        case 's': {
            // Force a '/edep/random/timeRandomSeed'
            setSeed = true;
            break;
        }
        case 'u': {
            // Force a '/edep/update' before running the first macro.
            doUpdate = true;
            break;
        }
        case 'U': {
            // Use a tcsh-style command line interface
            useUI = true;
            break;
        }
        case 'q':
        {
            // decrease the verbosity level.
            if (logLevel>0) --logLevel;
            else logLevel = 0;
            break;
        }
        case 'v':
        {
            // increase the verbosity level.
            if (logLevel>0) ++logLevel;
            else logLevel = 2;
            break;
        }
        case 'V':
        {
            // Set the debug level for a named trace.
            std::string arg(optarg);
            std::size_t sep = arg.find("=");
            if (sep != std::string::npos) {
                std::string name = arg.substr(0,sep);
                std::string levelName = arg.substr(sep+1);
                switch (levelName[0]) {
                case 'q': case 'Q':
                    namedLogLevel[name.c_str()]
                        = EDepSim::LogManager::QuietLevel;
                    break;
                case 'l': case 'L':
                    namedLogLevel[name.c_str()]
                        = EDepSim::LogManager::LogLevel;
                    break;
                case 'i': case 'I':
                    namedLogLevel[name.c_str()]
                        = EDepSim::LogManager::InfoLevel;
                    break;
                case 'v': case 'V':
                    namedLogLevel[name.c_str()]
                        = EDepSim::LogManager::VerboseLevel;
                    break;
                default:
                    usage();
                }
            }
            break;
        }
        case 'h':
        default:
            usage();
        }
    }

    if (logLevel < 1) {
        EDepSim::LogManager::SetLogLevel(EDepSim::LogManager::QuietLevel);
    }
    else if (logLevel == 1) {
        EDepSim::LogManager::SetLogLevel(EDepSim::LogManager::LogLevel);
        EDepSimLog("Set log level to LogLevel");
    }
    else if (logLevel == 2) {
        EDepSim::LogManager::SetLogLevel(EDepSim::LogManager::InfoLevel);
        EDepSimInfo("Set log level to InfoLevel");
    }
    else if (logLevel >= 3) {
        EDepSim::LogManager::SetLogLevel(EDepSim::LogManager::VerboseLevel);
        EDepSimVerbose("Set log level to VerboseLevel");
    }

    for (std::map<std::string,EDepSim::LogManager::LogPriority>::iterator i
             = namedLogLevel.begin();
         i != namedLogLevel.end();
         ++i) {
        EDepSim::LogManager::SetLogLevel(i->first.c_str(), i->second);
    }

    if (debugLevel == 1) {
        EDepSim::LogManager::SetDebugLevel(EDepSim::LogManager::WarnLevel);
        EDepSimWarn("Set debug level to WarnLevel");
    }
    else if (debugLevel == 2) {
        EDepSim::LogManager::SetDebugLevel(EDepSim::LogManager::DebugLevel);
        EDepSimDebug("Set debug level to DebugLevel");
    }
    else if (debugLevel >= 2) {
        EDepSim::LogManager::SetDebugLevel(EDepSim::LogManager::TraceLevel);
        EDepSimTrace("Set debug level to TraceLevel");
    }

    for (std::map<std::string,EDepSim::LogManager::ErrorPriority>::iterator i
             = namedDebugLevel.begin();
         i != namedDebugLevel.end();
         ++i) {
        EDepSim::LogManager::SetDebugLevel(i->first.c_str(), i->second);
    }

    // Set the mandatory initialization classes
    // Construct the default run manager
    G4RunManager* runManager = EDepSim::CreateRunManager(physicsList);

    // Create the persistency manager.  The persistency manager must derive
    // from G4VPersistencyManager which will make this object available to the
    // G4RunManager as a singleton.  There can only be one persistency manager
    // at a time.  The Store methods will then be called by the run managers
    // Analyze methods.  The persistency manager doesn't *have* to be derived
    // from EDepSim::RootPersistencyManager, but it's probably the best way to
    // save internal data structures.
    EDepSim::PersistencyManager* persistencyManager = NULL;

    if(hasSuffix(outputFilename,".root")) {
        persistencyManager = new EDepSim::RootPersistencyManager();
        persistencyManager->SetStoreROOT(true);
    }
    else if(hasSuffix(outputFilename,".h5") ||
            hasSuffix(outputFilename,".hdf5")) {
        persistencyManager = new EDepSim::H5PersistencyManager();
        persistencyManager->SetStoreHDF5(true);
    }
    else {
        persistencyManager = new EDepSim::PersistencyManager();
        EDepSimWarn("No file format specified.  No output will be saved.");
        outputFilename = "";
    }

    // Get the pointer to the UI manager.  This is used to control how
    // edep-sim will run.  All arguments are passed in using GEANT4 macro
    // commands.
    G4UImanager* UI = G4UImanager::GetUIpointer();

    // Open the file if one was declared on the command line.
    if (gdmlFilename != "") {
        UI->ApplyCommand("/edep/gdml/read "+gdmlFilename);
    }

    // Set the defaults for the simulation.
    UI->ApplyCommand("/edep/control edepsim-defaults 1.0");

    if(!outputFilename.empty())
        UI->ApplyCommand("/edep/db/open "+outputFilename);

    std::signal(SIGILL,  SIG_DFL);
    std::signal(SIGBUS,  SIG_DFL);
    std::signal(SIGSEGV, SIG_DFL);

    // Signal that the geometry should be validated for overlaps before
    // generating the first event.  This causes the executable to throw an
    // exception if the geometry has overlaps.
    if (validateGeometry) UI->ApplyCommand("/edep/validateGeometry");

    // Set the random seed from the time.
    if (setSeed) UI->ApplyCommand("/edep/random/timeRandomSeed");

    // Set the defaults for the simulation and get ready to run.  This needs
    // to be done before the first event is generated, but can also be done in
    // the users macro file.  It's executed here if the "-u" option was
    // provided on the command line.
    if (doUpdate) UI->ApplyCommand("/edep/update");

    if (useUI) {
        // Make a command line available if one was requested.
        G4UIExecutive* ui = new G4UIExecutive(argc, argv);
        for (int i=optind; i<argc; ++i) {
            std::string macroFilename = argv[i];
            std::cout << "## Run macro: " << macroFilename << std::endl;
            UI->ApplyCommand("/control/execute "+macroFilename);
        }
        ui->SessionStart();
        delete ui;
    }
    else if (optind < argc) {
        // Run a macro from the command line.
        for (int i=optind; i<argc; ++i) {
            std::string macroFilename = argv[i];
            UI->ApplyCommand("/control/execute " +macroFilename);
        }
        // If a event count was provided with the -e format, then add a beamOn
        // command.  The beamOn command can also be in the users macro, but
        // this makes it a little more convenient to write generic macros and
        // then specify the number of interactions to simulate from the
        // command line.
        if (!beamOnCount.empty()) {
            UI->ApplyCommand("/run/beamOn " + beamOnCount);
        }
    }

    // If we have the persistency manager, then make sure it's closed.
    if (persistencyManager) {
        persistencyManager->Close();
        delete persistencyManager;
    }
    delete runManager;

    return 0;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/edep-build.sh force"
// End:
