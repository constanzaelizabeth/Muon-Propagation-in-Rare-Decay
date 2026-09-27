#ifndef BELLMAGNETICFIELD_HH
#define BELLMAGNETICFIELD_HH
 
#include <memory>
 
#include "G4MagneticField.hh"
#include "G4String.hh"
 
class MuonShieldFieldMap;
 
// Wraps the real TRY_2026 field map -- a faithful Geant4+ROOT port of
// FairShip's ShipBFieldMap, reading the exact same files/TRY_2026.root
// grid, see MuonShieldFieldMap.{hh,cc} -- as a G4MagneticField.
//
// This is what FairShip actually uses for TRY_2026: WithConstField=False
// routes the shield's field through ShipFieldMaker::defineFieldMap(), not
// through the idealized per-region uniform-box field a plain read of
// ShipMuonShield::CreateMagnet() would suggest (see addVMCFields() in
// python/geomGeant4.py). The earlier box-model version of this class is
// gone; MuonShieldGeometry.{hh,cc} still builds the same box shapes, but
// now only for the passive iron volumes in DetectorConstruction -- the
// field itself no longer uses them.
class BellMagneticField : public G4MagneticField {
 public:
  // zStartWorld: world-local z (G4-internal units) of the shield's
  // entrance (magnet 0's front face) -- ShipMuonShield's Entrance[0],
  // which the map file's own local z=0 is defined relative to.
  // mapFilePath: path to TRY_2026.root (or whichever preset's map you're
  // using).
  BellMagneticField(G4double zStartWorld, const G4String& mapFilePath);
  ~BellMagneticField() override;
 
  void GetFieldValue(const G4double point[4], G4double* bField) const override;
 
 private:
  std::unique_ptr<MuonShieldFieldMap> fFieldMap;
};
 
#endif  // BELLMAGNETICFIELD_HH
 
