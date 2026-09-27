#ifndef MUONSHIELDFIELDMAP_HH
#define MUONSHIELDFIELDMAP_HH
 
#include <array>
#include <string>
#include <vector>
 
// Faithful, standalone (no FairLogger/TGeoMatrix/FairModule/TVirtualMagField
// dependency) Geant4+ROOT port of FairShip's field/ShipBFieldMap.{h,cxx}.
// Reads the exact same ROOT-file layout: a "Range" TTree with one entry
// (branches xMin/xMax/dx/yMin/yMax/dy/zMin/zMax/dz) and a "Data" TTree with
// one entry per grid node -- in ascending z, then y, then x order --
// giving Bx/By/Bz in Tesla. Same quadrant-symmetry sign rules (Bx
// antisymmetric in x and y, By symmetric, Bz antisymmetric in y only) and
// the same trilinear interpolation, transcribed 1:1 from ShipBFieldMap.cxx.
//
// Units: exactly ShipBFieldMap's own convention. Positions in and out are
// raw cm (ShipUnit's native unit), NOT G4-internal units; the field is
// returned in raw "ShipUnit kilogauss" (i.e. tesla = 10 in this
// convention -- the same one ShipBellField.cxx and ShipConstField.cxx use).
// BellMagneticField.cc is the only place that converts to/from G4's own
// units.
//
// Not ported: the text-file path (readTextFile()) and the rotation
// (phi/theta/psi) support -- addVMCFields() calls defineFieldMap with a
// .root file and zero rotation angles for the muon shield, so neither is
// exercised here.
class MuonShieldFieldMap {
 public:
  // mapFileName: path to the .root file (e.g. "files/TRY_2026.root").
  // zOffsetCm: z shift so the map's own local z=0 lands at this global z
  // (cm) -- ShipMuonShield's Entrance[0], the front face of magnet 0 (x
  // and y offsets are always 0 for this shield, matching addVMCFields()).
  // useSymmetry: quadrant-symmetry flag (true for TRY_2026 -- addVMCFields
  // hardcodes quadSymm = True for the muon shield map regardless of
  // preset).
  MuonShieldFieldMap(const std::string& mapFileName, double zOffsetCm,
                     bool useSymmetry);
 
  // xCm, yCm, zCm: global position, raw cm.
  // BkG[3]: (Bx, By, Bz) out, raw ShipUnit kilogauss (tesla = 10) -- NOT
  // yet G4 units.
  void Field(double xCm, double yCm, double zCm, double BkG[3]) const;
 
 private:
  void ReadRootFile(const std::string& mapFileName);
  void SetLimits();
  bool InsideRange(double x, double y, double z) const;
 
  struct BinInfo {
    int bin;
    double frac;
  };
  enum Axis { kX, kY, kZ };
  BinInfo GetBinInfo(double u, Axis axis) const;
  int GetMapBin(int iX, int iY, int iZ) const;
  static double TriLinearInterp(double A, double B, double C, double D,
                                double E, double F, double G, double H,
                                double xFrac, double xFrac1, double yFrac,
                                double yFrac1, double zFrac, double zFrac1);
 
  // (Bx, By, Bz) per grid node, raw ShipUnit kilogauss; index order
  // matches ShipBFieldMap::getMapBin: (iX * fNy + iY) * fNz + iZ.
  std::vector<std::array<double, 3>> fFieldMap;
 
  double fZOffsetCm{0.};
  bool fSymmetric{false};
 
  double fXmin{0.}, fXmax{0.}, fDx{0.};
  double fYmin{0.}, fYmax{0.}, fDy{0.};
  double fZmin{0.}, fZmax{0.}, fDz{0.};
  int fNx{0}, fNy{0}, fNz{0}, fN{0};
};
 
#endif  // MUONSHIELDFIELDMAP_HH
 
