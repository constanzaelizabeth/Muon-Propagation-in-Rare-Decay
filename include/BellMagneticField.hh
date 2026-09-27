#ifndef BELLMAGNETICFIELD_HH
#define BELLMAGNETICFIELD_HH
 
#include <memory>
 
#include "G4MagneticField.hh"
#include "G4String.hh"
 
class MuonShieldFieldMap;
 
class BellMagneticField : public G4MagneticField {
 public:

  BellMagneticField(G4double zStartWorld, const G4String& mapFilePath);
  ~BellMagneticField() override;
 
  void GetFieldValue(const G4double point[4], G4double* bField) const override;
 
 private:
  std::unique_ptr<MuonShieldFieldMap> fFieldMap;
};
 
#endif  // BELLMAGNETICFIELD_HH
 
