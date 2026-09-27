#ifndef DETECTORCONSTRUCTION_HH
#define DETECTORCONSTRUCTION_HH
 
#include <cstddef>
#include <string>
 
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
 
class G4LogicalVolume;
class G4Material;
 
// Beamline, along z:
//   [ tungsten target ] -> [ SHiP muon shield, "TRY_2026" preset, ported
//   1:1 from FairShip's ShipMuonShield.cxx (8 magnets, real field zones) ]
//   -> [ decay volume, vacuum ]
//
// The shield itself (geometry AND field) is no longer parametrized from
// outside: it's the real, hardcoded TRY_2026 configuration (see
// MuonShieldData.hh). Only the target and decay volume, which aren't part
// of FairShip's shield model, stay as constructor arguments.
class DetectorConstruction : public G4VUserDetectorConstruction {
 public:
  // apertureCm: half-width x2 of the (square) Target and DecayVolume
  // cross-section -- independent of the shield's own (much larger,
  // per-magnet-varying) transverse extent.
  // fieldMapPath: path to the TRY_2026 field map ROOT file (see
  // MuonShieldFieldMap and BellMagneticField) -- e.g. "files/TRY_2026.root",
  // wherever you've placed a copy of it.
  DetectorConstruction(G4double apertureCm, G4double targetLengthCm,
                       G4double decayVolumeLengthCm,
                       const std::string& fieldMapPath);
  ~DetectorConstruction() override = default;
 
  G4VPhysicalVolume* Construct() override;
  void ConstructSDandField() override;
 
  // World-local z of the decay volume's downstream face; used by
  // SteppingAction to count muons reaching it without duplicating this
  // geometry logic.
  G4double GetDecayVolumeEndZ() const { return fDecayVolumeEndZ; }
  // Half-width of the Target/DecayVolume aperture; used by SteppingAction
  // for the same muon-counting check.
  G4double GetApertureHalfWidth() const { return fAperture / 2.; }
 
 private:
  void BuildShield(G4LogicalVolume* worldLog, G4Material* iron);
  void BuildOneMagnet(G4LogicalVolume* worldLog, G4Material* iron,
                      std::size_t index, G4double zCenterWorld);
 
  G4double fAperture, fTargetLength, fDecayVolumeLength;
  std::string fFieldMapPath;
 
  // World-local z of the "target front face = 0" origin.
  G4double fZOffset{0.};
  // World-local z where the shield begins (= the target's back face; this
  // simulation has no separate "proximity shielding" volume between the
  // two, unlike full FairShip).
  G4double fShieldStartZ{0.};
  G4double fDecayVolumeEndZ{0.};
};
 
#endif  // DETECTORCONSTRUCTION_HH
 

