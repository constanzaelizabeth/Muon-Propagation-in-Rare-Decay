#ifndef TRACKINGACTION_HH
#define TRACKINGACTION_HH
 
#include <string>
 
#include "G4UserTrackingAction.hh"
#include "globals.hh"
 
#include "Rtypes.h"  // ROOT typedefs: Int_t, Long64_t, Double_t
 
class TTree;
class G4Track;
struct SharedEventInfo;
 
// Fills one row per recorded point into a flat ROOT TTree -- a generic
// stand-in for the 4 fixed TTrees the original main_B_HNL_mu.cc wrote by
// hand.
//
// Every G4Track gets at least one row, written at creation time by
// PreUserTrackingAction: this is the track's *creation vertex* (for a
// primary -- e.g. the muon injected by PrimaryGeneratorAction -- that's
// the B-decay vertex; for a secondary -- e.g. an electron from a muon
// decaying in flight -- that's the parent's decay vertex).
//
// SteppingAction additionally calls the public RecordPoint() below once
// per step for the tracks it cares about (currently: muons), appending
// one more row per step. Reading all rows that share a trackId, in order,
// reconstructs that track's complete trajectory; the last such row, if
// its "process" is "Decay", is by definition that track's decay vertex.
class TrackingAction : public G4UserTrackingAction {
 public:
  TrackingAction(TTree* tree, SharedEventInfo* shared);
  ~TrackingAction() override = default;
 
  void PreUserTrackingAction(const G4Track* track) override;
 
  // Escribe una fila con el estado ACTUAL de 'track' (posicion, momento,
  // tiempo, volumen) y el nombre de proceso dado. PreUserTrackingAction la
  // usa internamente para la fila de creacion (processName = proceso
  // creador, o "primary"); SteppingAction la llama una vez por step para
  // las filas de trayectoria (processName = proceso que definio ese step).
  void RecordPoint(const G4Track* track, const std::string& processName);
 
 private:
  TTree* fTree;
  SharedEventInfo* fShared;
 
  Int_t fEvent{0};
  Long64_t fCsvEvento{0};
  Int_t fTrackId{0}, fParentId{0}, fPdgId{0};
  Double_t fPx{0}, fPy{0}, fPz{0}, fE{0};
  Double_t fX{0}, fY{0}, fZ{0}, fT{0};
  Double_t fWeight{0};
  std::string fVolume, fProcess;
};
 
#endif  // TRACKINGACTION_HH
 
