#ifndef PRIMARYGENERATORACTION_HH
#define PRIMARYGENERATORACTION_HH
 
#include <random>
#include <string>
#include <vector>
 
#include "G4VUserPrimaryGeneratorAction.hh"
#include "Pythia8/Pythia.h"
#include "globals.hh"
 
#include "Rtypes.h"  // ROOT typedefs: Long64_t, Double_t
 
class G4Event;
class G4ParticleDefinition;
 
// Shared with TrackingAction so the output ROOT tree can carry the
// original CSV event number and per-row weight for every G4Track.
struct SharedEventInfo {
  Long64_t csvEvento = 0;
  Double_t weight = 1.0;
};
 
// Reads the same CSV as the original main_B_HNL_mu.cc (columns: evento,
// indice, id, m, px, py, pz, e, xProd, yProd, zProd, tProd, peso), and uses
// Pythia8 (decay-only) to force the B+/B- -> mu+/- N2(HNL) decay, sampling
// the target-interaction depth exactly as before.
//
// From that decay vertex onward, propagation is handed to Geant4: the two
// daughters are injected as G4PrimaryParticles on a single G4PrimaryVertex,
// and it is Geant4's own physics list + the global BellMagneticField that
// propagate them (energy loss, multiple scattering, magnetic deflection,
// the muon's decay-in-flight). The HNL is registered as a stable particle
// (see PhysicsList) so Geant4 simply carries it in a straight line.
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
 public:
  PrimaryGeneratorAction(const std::string& csvPath, G4double hnlMassGeV, G4int seed,
                         G4double targetThicknessCm, G4double interactionLengthCm,
                         SharedEventInfo* shared);
  ~PrimaryGeneratorAction() override = default;
 
  void GeneratePrimaries(G4Event* event) override;
 
  Long64_t GetNEvents() const { return static_cast<Long64_t>(fRows.size()); }
  bool IsReady() const { return fReady; }
 
 private:
  struct BRow {
    Long64_t evento;
    int id;
    double m, px, py, pz, e;
    double xProd_mm, yProd_mm, tProd_mmc;
    double peso;
  };
 
  bool LoadCsv(const std::string& csvPath);
 
  Pythia8::Pythia fPythia;
  bool fBPlusValido{false};
  double fCtau0B_mm{0.};
  bool fReady{false};
 
  std::vector<BRow> fRows;
  size_t fNextRow{0};
 
  std::mt19937_64 fVertexRng;
  std::uniform_real_distribution<double> fUniform01;
 
  double fTargetThickness_cm, fInteractionLength_cm;
 
  G4ParticleDefinition* fMuMinus{nullptr};
  G4ParticleDefinition* fMuPlus{nullptr};
  G4ParticleDefinition* fHnl{nullptr};
 
  SharedEventInfo* fShared;
};
 
#endif  // PRIMARYGENERATORACTION_HH
 
