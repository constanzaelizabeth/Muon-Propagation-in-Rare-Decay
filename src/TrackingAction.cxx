/*#include "TrackingAction.hh"
 
#include "PrimaryGeneratorAction.hh"  // SharedEventInfo
 
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4VPhysicalVolume.hh"
 
#include "TTree.h"
 
TrackingAction::TrackingAction(TTree* tree, SharedEventInfo* shared)
    : fTree(tree), fShared(shared) {
  fTree->Branch("event", &fEvent);
  fTree->Branch("csvEvento", &fCsvEvento);
  fTree->Branch("trackId", &fTrackId);
  fTree->Branch("parentId", &fParentId);
  fTree->Branch("pdgId", &fPdgId);
  fTree->Branch("px", &fPx);
  fTree->Branch("py", &fPy);
  fTree->Branch("pz", &fPz);
  fTree->Branch("e", &fE);
  fTree->Branch("x", &fX);
  fTree->Branch("y", &fY);
  fTree->Branch("z", &fZ);
  fTree->Branch("t", &fT);
  fTree->Branch("weight", &fWeight);
  fTree->Branch("volume", &fVolume);
  fTree->Branch("process", &fProcess);
}
 
void TrackingAction::PreUserTrackingAction(const G4Track* track) {
  fEvent = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  fCsvEvento = fShared ? fShared->csvEvento : 0;
  fWeight = fShared ? fShared->weight : 1.0;
  fTrackId = track->GetTrackID();
  fParentId = track->GetParentID();
  fPdgId = track->GetParticleDefinition()->GetPDGEncoding();
 
  const G4ThreeVector& p = track->GetMomentum();
  fPx = p.x() / GeV;
  fPy = p.y() / GeV;
  fPz = p.z() / GeV;
  fE = track->GetTotalEnergy() / GeV;
 
  const G4ThreeVector& x = track->GetPosition();
  fX = x.x() / cm;
  fY = x.y() / cm;
  fZ = x.z() / cm;
  fT = track->GetGlobalTime() / ns;
 
  fVolume = track->GetVolume() ? track->GetVolume()->GetName() : "outOfWorld";
  fProcess = track->GetCreatorProcess() ? track->GetCreatorProcess()->GetProcessName()
                                       : "primary";
  fTree->Fill();
}
 
void TrackingAction::RecordStep(const G4Step* step) {
  const G4Track* track = step->GetTrack();
  fEvent = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  fCsvEvento = fShared ? fShared->csvEvento : 0;
  fWeight = fShared ? fShared->weight : 1.0;
  fTrackId = track->GetTrackID();
  fParentId = track->GetParentID();
  fPdgId = track->GetParticleDefinition()->GetPDGEncoding();
 
  const G4StepPoint* post = step->GetPostStepPoint();
 
  const G4ThreeVector& p = post->GetMomentum();
  fPx = p.x() / GeV;
  fPy = p.y() / GeV;
  fPz = p.z() / GeV;
  fE = post->GetTotalEnergy() / GeV;
 
  const G4ThreeVector& x = post->GetPosition();
  fX = x.x() / cm;
  fY = x.y() / cm;
  fZ = x.z() / cm;
  fT = post->GetGlobalTime() / ns;
 
  fVolume = post->GetPhysicalVolume() ? post->GetPhysicalVolume()->GetName() : "outOfWorld";
  fProcess = post->GetProcessDefinedStep() ? post->GetProcessDefinedStep()->GetProcessName()
                                          : "step";
  fTree->Fill();
}*/
 
 
#include "TrackingAction.hh"
 
#include "PrimaryGeneratorAction.hh"  // SharedEventInfo
 
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4VPhysicalVolume.hh"
 
#include "TTree.h"
 
TrackingAction::TrackingAction(TTree* tree, SharedEventInfo* shared)
    : fTree(tree), fShared(shared) {
  fTree->Branch("event", &fEvent);
  fTree->Branch("csvEvento", &fCsvEvento);
  fTree->Branch("trackId", &fTrackId);
  fTree->Branch("parentId", &fParentId);
  fTree->Branch("pdgId", &fPdgId);
  fTree->Branch("px", &fPx);
  fTree->Branch("py", &fPy);
  fTree->Branch("pz", &fPz);
  fTree->Branch("e", &fE);
  fTree->Branch("x", &fX);
  fTree->Branch("y", &fY);
  fTree->Branch("z", &fZ);
  fTree->Branch("t", &fT);
  fTree->Branch("weight", &fWeight);
  fTree->Branch("volume", &fVolume);
  fTree->Branch("process", &fProcess);
}
 
void TrackingAction::PreUserTrackingAction(const G4Track* track) {
  // Vertice de creacion del track: "primary" si nace directamente de
  // PrimaryGeneratorAction (caso del muon y del HNL), o el nombre del
  // proceso que lo produjo (p.ej. "Decay" para los hijos de un muon que
  // decae en vuelo).
  const std::string processName =
      track->GetCreatorProcess() ? track->GetCreatorProcess()->GetProcessName() : "primary";
  RecordPoint(track, processName);
}
 
void TrackingAction::RecordPoint(const G4Track* track, const std::string& processName) {
  fEvent = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  fCsvEvento = fShared ? fShared->csvEvento : 0;
  fWeight = fShared ? fShared->weight : 1.0;
  fTrackId = track->GetTrackID();
  fParentId = track->GetParentID();
  fPdgId = track->GetParticleDefinition()->GetPDGEncoding();
 
  const G4ThreeVector& p = track->GetMomentum();
  fPx = p.x() / GeV;
  fPy = p.y() / GeV;
  fPz = p.z() / GeV;
  fE = track->GetTotalEnergy() / GeV;
 
  const G4ThreeVector& x = track->GetPosition();
  fX = x.x() / cm;
  fY = x.y() / cm;
  fZ = x.z() / cm;
  fT = track->GetGlobalTime() / ns;
 
  fVolume = track->GetVolume() ? track->GetVolume()->GetName() : "outOfWorld";
  fProcess = processName;
 
  fTree->Fill();
}
 
