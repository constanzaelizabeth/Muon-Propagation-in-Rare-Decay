#include "DetectorConstruction.hh"
#include "BellMagneticField.hh"
#include "MuonShieldGeometry.hh"
 
#include <algorithm>
#include <array>
#include <string>
 
#include "G4Box.hh"
#include "G4ChordFinder.hh"
#include "G4ClassicalRK4.hh"
#include "G4Colour.hh"
#include "G4FieldManager.hh"
#include "G4GenericTrap.hh"
#include "G4LogicalVolume.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4TransportationManager.hh"
#include "G4VisAttributes.hh"
 
DetectorConstruction::DetectorConstruction(G4double apertureCm, G4double targetLengthCm,
                                          G4double decayVolumeLengthCm,
                                          const std::string& fieldMapPath)
    : fAperture(apertureCm * cm),
      fTargetLength(targetLengthCm * cm),
      fDecayVolumeLength(decayVolumeLengthCm * cm),
      fFieldMapPath(fieldMapPath) {}
 
G4VPhysicalVolume* DetectorConstruction::Construct() {
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* worldMat = nist->FindOrBuildMaterial("G4_Galactic");
  G4Material* tungsten = nist->FindOrBuildMaterial("G4_W");
  G4Material* vacuum = nist->FindOrBuildMaterial("G4_Galactic");
  G4Material* iron = nist->FindOrBuildMaterial("G4_Fe");
 
  // z convention: target's upstream (front) face at z = 0. The shield
  // starts right at the target's back face -- ShipMuonShield's
  // "z_end_of_proximity_shielding", here simply the target's back face
  // since this simulation has no separate proximity-shielding volume.
  const G4double shieldStartZNoOffset = fTargetLength;
  const std::array<double, 8> centersCmNoOffset =
      MuonShieldGeometry::MagnetCentersCm(shieldStartZNoOffset / cm);
  const G4double shieldEndZNoOffset =
      (centersCmNoOffset[7] + kTry2026Magnets[7].halfLenZ) * cm;
  const G4double decayVolumeFrontZ = shieldEndZNoOffset;
  const G4double decayVolumeBackZ = decayVolumeFrontZ + fDecayVolumeLength;
 
  const G4double worldHalfLength = decayVolumeBackZ / 2. + 1. * m;
  // Transverse World size must fit the shield's own field zones, not just
  // the aperture.
  const G4double shieldTransverseHalf =
      MuonShieldGeometry::MaxTransverseExtentCm() * cm;
  const G4double worldTransverseHalf =
      std::max(fAperture / 2., shieldTransverseHalf) + 1. * m;
 
  // World is centered on the beamline; fZOffset re-expresses our
  // "target front face at z=0" convention inside World's local frame.
  fZOffset = -worldHalfLength;
  fShieldStartZ = fZOffset + shieldStartZNoOffset;
  fDecayVolumeEndZ = fZOffset + decayVolumeBackZ;
 
  auto* worldSolid =
      new G4Box("World", worldTransverseHalf, worldTransverseHalf, worldHalfLength);
  auto* worldLog = new G4LogicalVolume(worldSolid, worldMat, "World");
  worldLog->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto* worldPhys =
      new G4PVPlacement(nullptr, G4ThreeVector(), worldLog, "World", nullptr, false, 0);
 
  // ---   Target: tungsten box, front face at z = 0   ---------------------
  {
    auto* solid = new G4Box("Target", fAperture / 2., fAperture / 2., fTargetLength / 2.);
    auto* logic = new G4LogicalVolume(solid, tungsten, "Target");
    logic->SetVisAttributes(new G4VisAttributes(G4Colour(0.6, 0.6, 0.6)));
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, fZOffset + fTargetLength / 2.), logic,
                      "Target", worldLog, false, 0);
  }
 
  BuildShield(worldLog, iron);
 
  // ---   Decay volume: vacuum box, front face right after the shield   ---
  {
    auto* solid =
        new G4Box("DecayVolume", fAperture / 2., fAperture / 2., fDecayVolumeLength / 2.);
    auto* logic = new G4LogicalVolume(solid, vacuum, "DecayVolume");
    logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    new G4PVPlacement(
        nullptr, G4ThreeVector(0, 0, fZOffset + decayVolumeFrontZ + fDecayVolumeLength / 2.),
        logic, "DecayVolume", worldLog, false, 0);
  }
 
  return worldPhys;
}
 
void DetectorConstruction::BuildShield(G4LogicalVolume* worldLog, G4Material* iron) {
  // Single source of truth for magnet centers: BellMagneticField (built in
  // ConstructSDandField(), below) recomputes the exact same centers from
  // the exact same fShieldStartZ, so the passive geometry and the field
  // zones line up exactly.
  const std::array<double, 8> centersCm =
      MuonShieldGeometry::MagnetCentersCm(fShieldStartZ / cm);
  for (std::size_t i = 0; i < kTry2026Magnets.size(); ++i) {
    BuildOneMagnet(worldLog, iron, i, centersCm[i] * cm);
  }
}
 
void DetectorConstruction::BuildOneMagnet(G4LogicalVolume* worldLog, G4Material* iron,
                                          std::size_t index, G4double zCenterWorld) {
  // Ported 1:1 from ShipMuonShield::CreateMagnet()'s 8-piece decomposition
  // (see MuonShieldGeometry.{hh,cc} for the corner formulas and mirror
  // helpers). Placing all 8 real field zones as iron volumes -- not a
  // simplified "picture frame" placeholder -- so the passive material
  // budget (energy loss, multiple scattering) matches FairShip's.
  //
  // NOT ported: the outer "AbsorberVol" box FairShip subtracts magnet 0's
  // shapes from (extra solid iron filling the gaps *between* magnet 0's 8
  // field pieces specifically, making it double as the hadron absorber).
  // So magnet 0 here has somewhat less iron than the real hadron absorber
  // in those gaps. Flagged deliberately -- add it later if that material
  // budget matters for your analysis.
  const auto& p = kTry2026Magnets[index];
  const auto shapes = MuonShieldGeometry::BuildMagnetShapes(p);
  const std::string suffix = "_" + std::to_string(index);
  auto* vis = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));
 
  struct Piece {
    const char* name;
    MuonShieldGeometry::Corners corners;
  };
  const Piece pieces[8] = {
      {"MiddleMagL", shapes.mainL},
      {"MiddleMagR", MuonShieldGeometry::MirrorXY(shapes.mainL)},
      {"MagRetL", shapes.mainSideL},
      {"MagRetR", MuonShieldGeometry::MirrorXY(shapes.mainSideL)},
      {"MagTopLeft", shapes.topLeft},
      {"MagTopRight", MuonShieldGeometry::MirrorX(shapes.topLeft)},
      {"MagBotLeft", MuonShieldGeometry::MirrorY(shapes.topLeft)},
      {"MagBotRight", MuonShieldGeometry::MirrorXY(shapes.topLeft)},
  };
 
  for (int i = 0; i < 8; ++i) {
    const G4String name = std::string(pieces[i].name) + suffix;
    auto* solid = MuonShieldGeometry::MakeTrap(name, pieces[i].corners, p.halfLenZ);
    auto* logic = new G4LogicalVolume(solid, iron, name);
    logic->SetVisAttributes(vis);
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zCenterWorld), logic, name, worldLog,
                      false, static_cast<G4int>(index));
  }
}
 
void DetectorConstruction::ConstructSDandField() {
  // Global field: see BellMagneticField for the exact 8-magnet, per-zone
  // layout, ported from ShipMuonShield::ConstructGeometry(). Built from
  // the same fShieldStartZ that BuildShield() (above) used, so the field
  // zones line up exactly with the passive geometry.
  auto* field = new BellMagneticField(fShieldStartZ, fFieldMapPath);
  auto* equation = new G4Mag_UsualEqRhs(field);
  auto* stepper = new G4ClassicalRK4(equation);
  auto* chordFinder = new G4ChordFinder(field, 1.0 * mm, stepper);
 
  auto* fieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  fieldManager->SetDetectorField(field);
  fieldManager->SetChordFinder(chordFinder);
}
