#ifndef EDepSim_ParticleBombGenerator_hh
#define EDepSim_ParticleBombGenerator_hh

#include <string>
#include <fstream>

#include <globals.hh>
#include <G4LorentzVector.hh>

#include "kinem/EDepSimVKinematicsGenerator.hh"

#include "yaml-cpp/yaml.h"
#include "DLPGenerator/ParticleBomb/ParticleBomb.h"


class G4Event;
class G4VPrimaryGenerator;

namespace EDepSim
{
  class ParticleBombGenerator : public VKinematicsGenerator 
  {
    public:
      ParticleBombGenerator(const G4String&, const G4String&, bool); 
      ~ParticleBombGenerator();
      
      GeneratorStatus
      GeneratePrimaryVertex(G4Event*, const G4LorentzVector&);

    private:
      /// Strength requested by "ShootInward: True".  This is deliberately not
      /// DLPGenerator::kDEFAULT_SHOOT_INWARD_POWER, which is 0 and is what an
      /// absent key means.  True asks for the balanced default strength the
      /// generator documents, which is 1.
      static constexpr double kShootInwardOnPower = 1.;

      DLPGenerator::ParticleBomb _generator;
      DLPGenerator::GenParamParticle _parse_particle(const YAML::Node&);
      DLPGenerator::GenParamInteraction _parse_interaction(const YAML::Node&);
      double _parse_shoot_inward(const YAML::Node&);
  };
}
#endif
