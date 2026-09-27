/*// SPDX-License-Identifier: LGPL-3.0-or-later
//
// bToMuHnlSim.cc
// ---------------------------------------------------------------------
// Standalone Geant4 + Pythia8 + ROOT application (no FairRoot/FairShip).
// Reads the same B+/B- CSV as the original main_B_HNL_mu.cc and uses
// Pythia8 to compute the B+/- -> mu+/- N2(HNL) decay vertex exactly as
// before. From that vertex onward, Geant4 propagates the muon and the HNL
// through:
//
//   [ tungsten target ] -> [ SHiP magnet geometry, ported from
//   ShipMagnet.cxx ] -> [ decay volume, vacuum ]
//
// target and decay volume share a square apertureCm x apertureCm cross
// section. The magnet's field is a port of FairShip's own ShipBellField.
//
// SteppingAction records every step of every muon track into the MCTrack
// tree (its full propagated, field-curved trajectory -- not just its
// creation vertex), and counts live, during the run, how many muons reach
// the decay volume's downstream face while still inside its transverse
// aperture. Those counts are printed here and also stored in the output
// ROOT file, so no post-hoc counting is needed in plot_mu_fluence.C.
// ---------------------------------------------------------------------
 
#include <iostream>
#include <string>
 
#include "G4RunManagerFactory.hh"
#include "G4SystemOfUnits.hh"
 
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "SteppingAction.hh"
#include "TrackingAction.hh"
 
#include "TFile.h"
#include "TParameter.h"
#include "TTree.h"
 
int main(int argc, char** argv) {
  if (argc < 3 || argc > 10) {
    std::cerr << "Uso: bToMuHnlSim <input.csv> <output.root> [hnlMass=2.0] "
                "[seed=12345] [targetThickness_cm=150.0] "
                "[interactionLength_cm=9.9] [decayVolumeLength_cm=5000.0] "
                "[aperture_cm=25.0] [peakField_T=0.14361]\n";
    return 1;
  }
 
  const std::string inputCsv = argv[1];
  const std::string outputName = argv[2];
  const double hnlMass = (argc >= 4) ? std::stod(argv[3]) : 2.0;
  const int seed = (argc >= 5) ? std::stoi(argv[4]) : 12345;
  const double targetThickness_cm = (argc >= 6) ? std::stod(argv[5]) : 150.0;
  const double interactionLength_cm = (argc >= 7) ? std::stod(argv[6]) : 9.9;
  const double decayVolumeLength_cm = (argc >= 8) ? std::stod(argv[7]) : 5000.0;
  const double aperture_cm = (argc >= 9) ? std::stod(argv[8]) : 25.0;
  const double peakField_T = (argc >= 10) ? std::stod(argv[9]) : 0.14361;
 
  // -----   ROOT output   -------------------------------------------------
  auto* outFile = new TFile(outputName.c_str(), "RECREATE");
  auto* tree = new TTree("MCTrack", "Geant4 track truth (primaries + secondaries)");
 
  SharedEventInfo shared;
 
  // -----   Geant4 run manager (serial: matches the original single-
  //         threaded, deterministic-seed script)   -------------------------
  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
 
  // coilThick_cm debe cumplir coilThick_cm < aperture_cm/2 (ver
  // DetectorConstruction::BuildMagnet, G4Box "CVert": su semi-altura es
  // yap - fCoilThick, que debe ser > 0). Con aperture_cm=25.0 (default,
  // consistente con el rango x,y de -12.5 a 12.5 cm del CSV de entrada),
  // yap = 12.5 cm, así que coilThick_cm=5.0 deja margen (yap - coilThick =
  // 7.5 cm > 0) sin inflar la apertura real del imán/target.
  const double yokeWidth_cm = 200.0, yokeDepth_cm = 200.0, coilThick_cm = 5.0;
  auto* detector =
      new DetectorConstruction(aperture_cm, targetThickness_cm, decayVolumeLength_cm,
                               yokeWidth_cm, yokeDepth_cm, coilThick_cm, peakField_T);
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new PhysicsList(hnlMass));
 
  auto* generator = new PrimaryGeneratorAction(inputCsv, hnlMass, seed, targetThickness_cm,
                                               interactionLength_cm, &shared);
  runManager->SetUserAction(generator);
 
  auto* trackingAction = new TrackingAction(tree, &shared);
  runManager->SetUserAction(trackingAction);
 
  // SteppingAction necesita el TrackingAction (para escribir cada step de
  // muon en el arbol) y el DetectorConstruction (para saber donde termina
  // el decay volume y cual es su apertura), asi que se crea despues de
  // ambos.
  auto* steppingAction = new SteppingAction(trackingAction, detector, &shared);
  runManager->SetUserAction(steppingAction);
 
  runManager->Initialize();
 
  if (!generator->IsReady()) {
    std::cerr << "Error: no se pudo inicializar el generador (revisa el CSV y mHNL).\n";
    delete runManager;
    outFile->Close();
    return 1;
  }
 
  const G4int nEvents = static_cast<G4int>(generator->GetNEvents());
  std::cout << "Corriendo " << nEvents << " eventos...\n";
  runManager->BeamOn(nEvents);
 
  // -----   Conteo de muones al final del decay volume, hecho durante la
  //         corrida por SteppingAction (no interpolado despues) -----------
  std::cout << "----------------------------------------------------------\n";
  std::cout << "Muones que llegan al final del decay volume (Z="
            << detector->GetDecayVolumeEndZ() / m << " m, dentro de la apertura |x|,|y|<="
            << detector->GetApertureHalfWidth() / m << " m):\n";
  std::cout << "  Conteo crudo (sin peso) : " << steppingAction->GetRawMuonCount() << "\n";
  std::cout << "  Conteo pesado (muon/pr) : " << steppingAction->GetWeightedMuonCount()
            << "\n";
  std::cout << "----------------------------------------------------------\n";
 
  tree->Write();
  auto* rawCountParam =
      new TParameter<double>("nMuonesDecayVolumeRaw", steppingAction->GetRawMuonCount());
  auto* weightedCountParam = new TParameter<double>("nMuonesDecayVolumePesado",
                                                    steppingAction->GetWeightedMuonCount());
  rawCountParam->Write();
  weightedCountParam->Write();
  outFile->Close();
  std::cout << "Salida: " << outputName << '\n';
 
  delete runManager;
  return 0;
}
*/
 
// SPDX-License-Identifier: LGPL-3.0-or-later
//
// bToMuHnlSim.cc
// ---------------------------------------------------------------------
// Standalone Geant4 + Pythia8 + ROOT application (no FairRoot/FairShip).
// Reads the same B+/B- CSV as the original main_B_HNL_mu.cc and uses
// Pythia8 to compute the B+/- -> mu+/- N2(HNL) decay vertex exactly as
// before. From that vertex onward, Geant4 propagates the muon and the HNL
// through:
//
//   [ tungsten target ] -> [ SHiP magnet geometry, ported from
//   ShipMagnet.cxx ] -> [ decay volume, vacuum ]
//
// target and decay volume share a square apertureCm x apertureCm cross
// section. The magnet's field is a port of FairShip's own ShipBellField.
//
// SteppingAction records every step of every muon track into the MCTrack
// tree (its full propagated, field-curved trajectory -- not just its
// creation vertex), and counts live, during the run, how many muons reach
// the decay volume's downstream face while still inside its transverse
// aperture. Those counts are printed here and also stored in the output
// ROOT file, so no post-hoc counting is needed in plot_mu_fluence.C.
// ---------------------------------------------------------------------
// SPDX-License-Identifier: LGPL-3.0-or-later
//
// bToMuHnlSim.cc
// ---------------------------------------------------------------------
// Standalone Geant4 + Pythia8 + ROOT application (no FairRoot/FairShip).
// Reads the same B+/B- CSV as the original main_B_HNL_mu.cc and uses
// Pythia8 to compute the B+/- -> mu+/- N2(HNL) decay vertex exactly as
// before. From that vertex onward, Geant4 propagates the muon and the HNL
// through:
//
//   [ tungsten target ] -> [ series of SHiP magnet units, ported from
//   ShipMagnet.cxx, alternating field polarity ] -> [ decay volume, vacuum ]
//
// target and decay volume share a square apertureCm x apertureCm cross
// section. Each magnet's field is a port of FairShip's own ShipBellField;
// consecutive magnets carry opposite polarity (see BellMagneticField) so a
// charged particle is swept further out, rather than bent back, as it
// crosses the series.
// ---------------------------------------------------------------------
 
// SPDX-License-Identifier: LGPL-3.0-or-later
//
// bToMuHnlSim.cc
// ---------------------------------------------------------------------
// Standalone Geant4 + Pythia8 + ROOT application (no FairRoot/FairShip).
// Reads the same B+/B- CSV as the original main_B_HNL_mu.cc and uses
// Pythia8 to compute the B+/- -> mu+/- N2(HNL) decay vertex exactly as
// before. From that vertex onward, Geant4 propagates the muon and the HNL
// through:
//
//   [ tungsten target ] -> [ series of SHiP magnet units, ported from
//   ShipMagnet.cxx, alternating field polarity ] -> [ decay volume, vacuum ]
//
// target and decay volume share a square apertureCm x apertureCm cross
// section. Each magnet's field is a port of FairShip's own ShipBellField;
// consecutive magnets carry opposite polarity (see BellMagneticField) so a
// charged particle is swept further out, rather than bent back, as it
// crosses the series.
// ---------------------------------------------------------------------
 
// SPDX-License-Identifier: LGPL-3.0-or-later
//
// bToMuHnlSim.cc
// ---------------------------------------------------------------------
// Standalone Geant4 + Pythia8 + ROOT application (no FairRoot/FairShip).
// Reads the same B+/B- CSV as the original main_B_HNL_mu.cc and uses
// Pythia8 to compute the B+/- -> mu+/- N2(HNL) decay vertex exactly as
// before. From that vertex onward, Geant4 propagates the muon and the HNL
// through:
//
//   [ tungsten target ] -> [ series of magnet units whose field logic is
//   ported from FairShip's ShipMuonShield::ConstructGeometry() -- see
//   BellMagneticField's header comment ] -> [ decay volume, vacuum ]
//
// target and decay volume share a square apertureCm x apertureCm cross
// section. This remains an approximation of SHiP's real, 41-parameter-
// optimized, 7/8-magnet baseline shield: the field ZONE LAYOUT and SIGN
// LOGIC are ported directly from ShipMuonShield.cxx, but per-magnet
// dimensions are shared (not individually tapered/tuned) since the actual
// optimized shield_params were not available when this was written.
// ---------------------------------------------------------------------
 
// SPDX-License-Identifier: LGPL-3.0-or-later
//
// bToMuHnlSim.cc
// ---------------------------------------------------------------------
// Standalone Geant4 + Pythia8 + ROOT application (no FairRoot/FairShip).
// Reads the same B+/B- CSV as the original main_B_HNL_mu.cc and uses
// Pythia8 to compute the B+/- -> mu+/- N2(HNL) decay vertex exactly as
// before. From that vertex onward, Geant4 propagates the muon and the HNL
// through:
//
//   [ tungsten target ] -> [ SHiP muon shield, "TRY_2026" preset, ported
//   1:1 from FairShip's ShipMuonShield.cxx (8 magnets, real geometry) ]
//   -> [ decay volume, vacuum ]
//
// target and decay volume share a square apertureCm x apertureCm cross
// section. The shield's field is the real TRY_2026 field map (a port of
// FairShip's ShipBFieldMap, see MuonShieldFieldMap/BellMagneticField) --
// not an idealized analytic formula -- read from the ROOT file whose path
// is given on the command line.
// ---------------------------------------------------------------------
 
// SPDX-License-Identifier: LGPL-3.0-or-later
//
// bToMuHnlSim.cc
// ---------------------------------------------------------------------
// Standalone Geant4 + Pythia8 + ROOT application (no FairRoot/FairShip).
// Reads the same B+/B- CSV as the original main_B_HNL_mu.cc and uses
// Pythia8 to compute the B+/- -> mu+/- N2(HNL) decay vertex exactly as
// before. From that vertex onward, Geant4 propagates the muon and the HNL
// through:
//
//   [ tungsten target ] -> [ SHiP muon shield, "TRY_2026" preset, ported
//   1:1 from FairShip's ShipMuonShield.cxx (8 magnets, real geometry) ]
//   -> [ decay volume, vacuum ]
//
// target and decay volume share a square apertureCm x apertureCm cross
// section. The shield's field is the real TRY_2026 field map (a port of
// FairShip's ShipBFieldMap, see MuonShieldFieldMap/BellMagneticField) --
// not an idealized analytic formula -- read from the ROOT file whose path
// is given on the command line.
// ---------------------------------------------------------------------
 
// SPDX-License-Identifier: LGPL-3.0-or-later
//
// bToMuHnlSim.cc
// ---------------------------------------------------------------------
// Standalone Geant4 + Pythia8 + ROOT application (no FairRoot/FairShip).
// Reads the same B+/B- CSV as the original main_B_HNL_mu.cc and uses
// Pythia8 to compute the B+/- -> mu+/- N2(HNL) decay vertex exactly as
// before. From that vertex onward, Geant4 propagates the muon and the HNL
// through:
//
//   [ tungsten target ] -> [ SHiP muon shield, "TRY_2026" preset, ported
//   1:1 from FairShip's ShipMuonShield.cxx (8 magnets, real geometry) ]
//   -> [ decay volume, vacuum ]
//
// target and decay volume share a square apertureCm x apertureCm cross
// section. The shield's field is the real TRY_2026 field map (a port of
// FairShip's ShipBFieldMap, see MuonShieldFieldMap/BellMagneticField) --
// not an idealized analytic formula -- read from the ROOT file whose path
// is given on the command line.
// ---------------------------------------------------------------------
 
#include <iostream>
#include <stdexcept>
#include <string>
 
#include "G4RunManagerFactory.hh"
#include "G4SystemOfUnits.hh"
 
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "TrackingAction.hh"
 
#include "TFile.h"
#include "TTree.h"
 
int main(int argc, char** argv) {
  if (argc < 3 || argc > 10) {
    std::cerr << "Uso: bToMuHnlSim <input.csv> <output.root> [hnlMass=3.0] "
                "[seed=12345] [targetThickness_cm=150.0] "
                "[interactionLength_cm=9.9] [decayVolumeLength_cm=5000.0] "
                "[aperture_cm=25.0] [fieldMapPath=files/TRY_2026.root]\n"
                "Notas:\n"
                " - El blindaje muonico (geometria y campo) ya no es "
                "configurable: es el preset real TRY_2026 de FairShip, "
                "8 imanes (ver MuonShieldData.hh).\n"
                " - fieldMapPath debe apuntar al archivo TRY_2026.root "
                "real (el mapa de campo -- ver MuonShieldFieldMap).\n";
    return 1;
  }
 
  const std::string inputCsv = argv[1];
  const std::string outputName = argv[2];
  const double hnlMass = (argc >= 4) ? std::stod(argv[3]) : 3.0;
  const int seed = (argc >= 5) ? std::stoi(argv[4]) : 12345;
  const double targetThickness_cm = (argc >= 6) ? std::stod(argv[5]) : 150.0;
  const double interactionLength_cm = (argc >= 7) ? std::stod(argv[6]) : 9.9;
  const double decayVolumeLength_cm = (argc >= 8) ? std::stod(argv[7]) : 5000.0;
  const double aperture_cm = (argc >= 9) ? std::stod(argv[8]) : 25.0;
  const std::string fieldMapPath =
      (argc >= 10) ? argv[9] : "files/TRY_2026.root";
 
  // -----   ROOT output   -------------------------------------------------
  auto* outFile = new TFile(outputName.c_str(), "RECREATE");
  auto* tree = new TTree("MCTrack", "Geant4 track truth (primaries + secondaries)");
 
  SharedEventInfo shared;
 
  // -----   Geant4 run manager (serial: matches the original single-
  //         threaded, deterministic-seed script)   -------------------------
  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
 
  auto* detector = new DetectorConstruction(aperture_cm, targetThickness_cm,
                                            decayVolumeLength_cm, fieldMapPath);
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new PhysicsList(hnlMass));
 
  auto* generator = new PrimaryGeneratorAction(inputCsv, hnlMass, seed, targetThickness_cm,
                                               interactionLength_cm, &shared);
  runManager->SetUserAction(generator);
  runManager->SetUserAction(new TrackingAction(tree, &shared));
 
  try {
    runManager->Initialize();
  } catch (const std::exception& e) {
    std::cerr << "Error al inicializar la geometria/campo: " << e.what() << "\n"
              << "(revisa que fieldMapPath apunte al archivo TRY_2026.root real)\n";
    delete runManager;
    outFile->Close();
    return 1;
  }
 
  if (!generator->IsReady()) {
    std::cerr << "Error: no se pudo inicializar el generador (revisa el CSV y mHNL).\n";
    delete runManager;
    outFile->Close();
    return 1;
  }
 
  const G4int nEvents = static_cast<G4int>(generator->GetNEvents());
  std::cout << "Corriendo " << nEvents << " eventos...\n";
  runManager->BeamOn(nEvents);
 
  outFile->cd();
  tree->Write();
  outFile->Close();
  std::cout << "Salida: " << outputName << '\n';
 
  delete runManager;
  return 0;
}
