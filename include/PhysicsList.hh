#ifndef PHYSICSLIST_HH
#define PHYSICSLIST_HH
 
#include "FTFP_BERT.hh"
#include "globals.hh"
 
// FTFP_BERT: Geant4's standard hadronic+EM physics list, the same family
// FairShip itself uses (see gconfig/g4Config.C: "FTFP_BERT_HP_EMZ"). On top
// of it, this registers one extra particle: the HNL. It is stable (no
// G4DecayTable), so no G4Decay process is attached to it, and being an
// otherwise-unrecognised exotic neutral particle, the standard physics
// list attaches no EM/hadronic processes either -- only G4Transportation
// acts on it, i.e. it is simply carried through the geometry in a straight
// line, exactly like the original script's Pythia "mayDecay = off" choice.
class PhysicsList : public FTFP_BERT {
 public:
  explicit PhysicsList(G4double hnlMassGeV);
  ~PhysicsList() override = default;
 
  void ConstructParticle() override;
 
 private:
  G4double fHnlMassGeV;
};
 
#endif  // PHYSICSLIST_HH
