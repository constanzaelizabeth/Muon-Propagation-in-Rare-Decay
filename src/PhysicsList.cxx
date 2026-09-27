#include "PhysicsList.hh"
 
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
 
namespace {
constexpr G4int kHnlPdg = 9900015;
}
 
PhysicsList::PhysicsList(G4double hnlMassGeV) : FTFP_BERT(), fHnlMassGeV(hnlMassGeV) {}
 
void PhysicsList::ConstructParticle() {
  FTFP_BERT::ConstructParticle();
 
  G4ParticleTable* table = G4ParticleTable::GetParticleTable();
  if (table->FindParticle(kHnlPdg)) return;  // already defined
 
  // NOTE: the G4ParticleDefinition constructor argument order below
  // (name, mass, width, charge, 2*spin, parity, C-conjugation, isospin,
  // isospin3, G-parity, particleType, lepton#, baryon#, PDG encoding,
  // stable, lifetime, decay table) is the classic/stable signature used
  // across Geant4 releases. If your installed version's
  // G4ParticleDefinition.hh differs, adjust the argument order to match.
  new G4ParticleDefinition("N2",              // name
                           fHnlMassGeV * GeV,  // mass
                           0.0,                // decay width
                           0.0,                // charge
                           1,                  // 2*spin (cosmetic only)
                           0,                  // parity
                           0,                  // C-conjugation
                           0,                  // isospin
                           0,                  // isospin3
                           0,                  // G-parity
                           "lepton",           // particle type (cosmetic)
                           0,                  // lepton number
                           0,                  // baryon number
                           kHnlPdg,            // PDG encoding
                           true,               // stable: no G4Decay process
                           -1.0,               // lifetime (unused, stable)
                           nullptr);           // no decay table
}
 
