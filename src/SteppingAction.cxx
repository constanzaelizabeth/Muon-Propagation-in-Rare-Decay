/*#include "SteppingAction.hh"
 
#include <cmath>
#include <cstdlib>
 
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"  // SharedEventInfo
#include "TrackingAction.hh"
 
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
 
SteppingAction::SteppingAction(TrackingAction* trackingAction,
                               const DetectorConstruction* detector, SharedEventInfo* shared)
    : fTrackingAction(trackingAction), fDetector(detector), fShared(shared) {}
 
void SteppingAction::UserSteppingAction(const G4Step* step) {
  const G4Track* track = step->GetTrack();
  if (std::abs(track->GetParticleDefinition()->GetPDGEncoding()) != 13) return;  // solo muones
 
  // Guarda este step para que la trayectoria real (curvada por el campo)
  // del muon quede en el arbol, no solo su vertice de creacion.
  if (fTrackingAction) fTrackingAction->RecordStep(step);
 
  // Conteo en vivo de muones que llegan al final del decay volume dentro
  // de su apertura transversal. El set de "ya contados" se reinicia al
  // empezar cada evento (los trackId se reutilizan evento a evento).
  const G4int eventId = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  if (eventId != fCurrentEventId) {
    fCountedTrackIds.clear();
    fCurrentEventId = eventId;
  }
 
  const G4int trackId = track->GetTrackID();
  if (fCountedTrackIds.count(trackId)) return;  // este track ya fue contado
 
  if (!fDetector) return;
  const G4double zEnd = fDetector->GetDecayVolumeEndZ();
  const G4double zPre = step->GetPreStepPoint()->GetPosition().z();
  const G4double zPost = step->GetPostStepPoint()->GetPosition().z();
  if (!(zPre < zEnd && zPost >= zEnd)) return;  // no acaba de cruzar hacia adelante
 
  // Interpolacion lineal entre los dos puntos del step para hallar la
  // posicion transversal en el cruce. Los steps son cortos -- sobre todo
  // dentro de la region de campo, gracias a la precision de 1mm del chord
  // finder -- asi que esto es una buena aproximacion del cruce real.
  const G4double xPre = step->GetPreStepPoint()->GetPosition().x();
  const G4double yPre = step->GetPreStepPoint()->GetPosition().y();
  const G4double xPost = step->GetPostStepPoint()->GetPosition().x();
  const G4double yPost = step->GetPostStepPoint()->GetPosition().y();
  const G4double dz = zPost - zPre;
  G4double xCross = xPost;
  G4double yCross = yPost;
  if (dz > 1e-9 * mm) {
    const G4double frac = (zEnd - zPre) / dz;
    xCross = xPre + frac * (xPost - xPre);
    yCross = yPre + frac * (yPost - yPre);
  }
 
  const G4double half = fDetector->GetApertureHalfWidth();
  if (std::fabs(xCross) > half || std::fabs(yCross) > half) return;  // salio de la apertura
 
  fCountedTrackIds.insert(trackId);
  fRawCount += 1.0;
  fWeightedCount += fShared ? fShared->weight : 1.0;
}*/
 
#include "SteppingAction.hh"
 
#include <cmath>
 
#include "DetectorConstruction.hh"
#include "TrackingAction.hh"
#include "PrimaryGeneratorAction.hh"  // SharedEventInfo
 
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
 
namespace {
constexpr G4int kMuMinusPdg = 13;
constexpr G4int kMuPlusPdg = -13;
}  // namespace
 
SteppingAction::SteppingAction(TrackingAction* trackingAction,
                              const DetectorConstruction* detector, SharedEventInfo* shared)
    : fTrackingAction(trackingAction), fDetector(detector), fShared(shared) {}
 
void SteppingAction::UserSteppingAction(const G4Step* step) {
  const G4Track* track = step->GetTrack();
  const G4int pdg = track->GetParticleDefinition()->GetPDGEncoding();
  if (pdg != kMuMinusPdg && pdg != kMuPlusPdg) return;
 
  // ---   Trayectoria: una fila mas por step, en el mismo MCTrack tree   --
  if (fTrackingAction) {
    const G4VProcess* proc = step->GetPostStepPoint()->GetProcessDefinedStep();
    const std::string processName = proc ? proc->GetProcessName() : "Unknown";
    fTrackingAction->RecordPoint(track, processName);
  }
 
  // ---   Conteo de muones que llegan al final del decay volume   ---------
  if (!fDetector) return;
 
  const G4double zEnd = fDetector->GetDecayVolumeEndZ();
  const G4double halfAperture = fDetector->GetApertureHalfWidth();
 
  const G4ThreeVector& prePos = step->GetPreStepPoint()->GetPosition();
  const G4ThreeVector& postPos = step->GetPostStepPoint()->GetPosition();
 
  // Solo contar el cruce hacia adelante del plano z = zEnd (evita doble
  // conteo si un step ya empezara justo sobre el plano).
  if (prePos.z() < zEnd && postPos.z() >= zEnd && std::fabs(postPos.x()) <= halfAperture &&
      std::fabs(postPos.y()) <= halfAperture) {
    const G4int eventId = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    const G4int trackId = track->GetTrackID();
    if (fCountedTracks.insert({eventId, trackId}).second) {
      fRawCount += 1.0;
      fWeightedCount += fShared ? fShared->weight : 1.0;
    }
  }
}
