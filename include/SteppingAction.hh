#ifndef STEPPINGACTION_HH
#define STEPPINGACTION_HH
 
#include <set>
#include <utility>
 
#include "G4UserSteppingAction.hh"
#include "globals.hh"
 
class TrackingAction;
class DetectorConstruction;
class G4Step;
struct SharedEventInfo;
 
// Guarda, para todo track de muon (mu-/mu+), un punto por step en el
// mismo MCTrack tree que llena TrackingAction: la union de esas filas
// (mismo trackId, ordenadas por tiempo) reconstruye su trayectoria
// completa, curvada por el campo. El vertice de creacion del muon ya
// queda escrito aparte por TrackingAction::PreUserTrackingAction (fila
// con process=="primary"); si el muon decae en vuelo, la ULTIMA fila que
// agrega esta clase queda con process=="Decay" y es, por definicion, su
// vertice de decaimiento.
//
// Ademas cuenta, durante la corrida (no interpolado despues), cuantos
// muones cruzan hacia adelante la cara final del decay volume dentro de
// su apertura transversal: conteo crudo (numero de muones) y conteo
// pesado (suma de shared->weight de esos muones).
class SteppingAction : public G4UserSteppingAction {
 public:
  SteppingAction(TrackingAction* trackingAction, const DetectorConstruction* detector,
                SharedEventInfo* shared);
  ~SteppingAction() override = default;
 
  void UserSteppingAction(const G4Step* step) override;
 
  G4double GetRawMuonCount() const { return fRawCount; }
  G4double GetWeightedMuonCount() const { return fWeightedCount; }
 
 private:
  TrackingAction* fTrackingAction;
  const DetectorConstruction* fDetector;
  SharedEventInfo* fShared;
 
  G4double fRawCount{0.};
  G4double fWeightedCount{0.};
 
  // (eventId, trackId) ya contados, para no contar dos veces el mismo
  // muon si su trayectoria cruzara el plano final mas de una vez.
  std::set<std::pair<G4int, G4int>> fCountedTracks;
};
 
#endif  // STEPPINGACTION_HH
