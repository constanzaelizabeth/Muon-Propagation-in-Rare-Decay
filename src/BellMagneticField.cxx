#include "BellMagneticField.hh"
 
#include "G4SystemOfUnits.hh"
#include "MuonShieldFieldMap.hh"
 
namespace {
// ShipBFieldMap's own convention: tesla = 10 in its "kilogauss" unit
// system (same as ShipBellField.cxx / ShipConstField.cxx). Converting one
// of its field values to G4-internal units is therefore (value/10)*tesla.
constexpr G4double kKilogaussToG4 = tesla / 10.0;
}  // namespace
 
BellMagneticField::BellMagneticField(G4double zStartWorld,
                                     const G4String& mapFilePath)
    : fFieldMap(std::make_unique<MuonShieldFieldMap>(
          std::string(mapFilePath), zStartWorld / cm, /*useSymmetry=*/true)) {}
 
BellMagneticField::~BellMagneticField() = default;
 
void BellMagneticField::GetFieldValue(const G4double point[4],
                                      G4double* bField) const {
  const double xCm = point[0] / cm;
  const double yCm = point[1] / cm;
  const double zCm = point[2] / cm;
 
  double BkG[3];
  fFieldMap->Field(xCm, yCm, zCm, BkG);
 
  bField[0] = BkG[0] * kKilogaussToG4;
  bField[1] = BkG[1] * kKilogaussToG4;
  bField[2] = BkG[2] * kKilogaussToG4;
}
 
