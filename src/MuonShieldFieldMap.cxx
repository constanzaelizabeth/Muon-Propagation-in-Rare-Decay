#include "MuonShieldFieldMap.hh"
 
#include <stdexcept>
 
#include "TFile.h"
#include "TTree.h"
 
namespace {
// ShipBFieldMap's own Tesla_ member: the file stores Bx/By/Bz in Tesla,
// multiplied by this at load time so everything downstream lives in
// ShipUnit's native "kilogauss" convention (tesla = 10), matching every
// other FairShip field class we've ported so far.
constexpr double kTeslaToKilogauss = 10.0;
}  // namespace
 
MuonShieldFieldMap::MuonShieldFieldMap(const std::string& mapFileName,
                                       double zOffsetCm, bool useSymmetry)
    : fZOffsetCm(zOffsetCm), fSymmetric(useSymmetry) {
  ReadRootFile(mapFileName);
}
 
void MuonShieldFieldMap::ReadRootFile(const std::string& mapFileName) {
  TFile* file = TFile::Open(mapFileName.c_str());
  if (!file || file->IsZombie()) {
    throw std::runtime_error("MuonShieldFieldMap: could not open " + mapFileName);
  }
 
  auto* rangeTree = dynamic_cast<TTree*>(file->Get("Range"));
  if (!rangeTree) {
    file->Close();
    throw std::runtime_error("MuonShieldFieldMap: no 'Range' tree in " + mapFileName);
  }
 
  Float_t xMin, xMax, dx, yMin, yMax, dy, zMin, zMax, dz;
  rangeTree->SetBranchAddress("xMin", &xMin);
  rangeTree->SetBranchAddress("xMax", &xMax);
  rangeTree->SetBranchAddress("dx", &dx);
  rangeTree->SetBranchAddress("yMin", &yMin);
  rangeTree->SetBranchAddress("yMax", &yMax);
  rangeTree->SetBranchAddress("dy", &dy);
  rangeTree->SetBranchAddress("zMin", &zMin);
  rangeTree->SetBranchAddress("zMax", &zMax);
  rangeTree->SetBranchAddress("dz", &dz);
  rangeTree->GetEntry(0);
 
  fXmin = xMin;
  fXmax = xMax;
  fDx = dx;
  fYmin = yMin;
  fYmax = yMax;
  fDy = dy;
  fZmin = zMin;
  fZmax = zMax;
  fDz = dz;
  SetLimits();
 
  auto* dataTree = dynamic_cast<TTree*>(file->Get("Data"));
  if (!dataTree) {
    file->Close();
    throw std::runtime_error("MuonShieldFieldMap: no 'Data' tree in " + mapFileName);
  }
 
  Float_t Bx, By, Bz;
  dataTree->SetBranchStatus("*", 0);
  dataTree->SetBranchStatus("Bx", 1);
  dataTree->SetBranchStatus("By", 1);
  dataTree->SetBranchStatus("Bz", 1);
  dataTree->SetBranchAddress("Bx", &Bx);
  dataTree->SetBranchAddress("By", &By);
  dataTree->SetBranchAddress("Bz", &Bz);
 
  const Long64_t nEntries = dataTree->GetEntries();
  if (nEntries != fN) {
    file->Close();
    throw std::runtime_error("MuonShieldFieldMap: expected " + std::to_string(fN) +
                             " entries in 'Data', found " +
                             std::to_string(nEntries));
  }
 
  fFieldMap.reserve(nEntries);
  for (Long64_t i = 0; i < nEntries; ++i) {
    dataTree->GetEntry(i);
    fFieldMap.push_back({static_cast<double>(Bx) * kTeslaToKilogauss,
                         static_cast<double>(By) * kTeslaToKilogauss,
                         static_cast<double>(Bz) * kTeslaToKilogauss});
  }
 
  file->Close();
  delete file;
}
 
void MuonShieldFieldMap::SetLimits() {
  // The number of bins includes both the minimum and maximum values; +1.5
  // (not +1.0) rounds correctly to the nearest integer, exactly as
  // ShipBFieldMap::setLimits() does.
  fNx = (fDx > 0.0) ? static_cast<int>((fXmax - fXmin) / fDx + 1.5) : 0;
  fNy = (fDy > 0.0) ? static_cast<int>((fYmax - fYmin) / fDy + 1.5) : 0;
  fNz = (fDz > 0.0) ? static_cast<int>((fZmax - fZmin) / fDz + 1.5) : 0;
 
  if (fNx < 2 || fNy < 2 || fNz < 2) {
    throw std::runtime_error(
        "MuonShieldFieldMap: field map needs at least 2 bins per axis");
  }
 
  fN = fNx * fNy * fNz;
}
 
bool MuonShieldFieldMap::InsideRange(double x, double y, double z) const {
  return x >= fXmin && x <= fXmax && y >= fYmin && y <= fYmax && z >= fZmin &&
        z <= fZmax;
}
 
MuonShieldFieldMap::BinInfo MuonShieldFieldMap::GetBinInfo(double u,
                                                            Axis axis) const {
  double du = 0.0, uMin = 0.0;
  int Nu = 0;
  switch (axis) {
    case kX:
      du = fDx;
      uMin = fXmin;
      Nu = fNx;
      break;
    case kY:
      du = fDy;
      uMin = fYmin;
      Nu = fNy;
      break;
    case kZ:
      du = fDz;
      uMin = fZmin;
      Nu = fNz;
      break;
  }
 
  int iBin = -1;
  double frac = 0.0;
  if (du > 1e-10) {
    const double dist = (u - uMin) / du;
    iBin = static_cast<int>(dist);
    frac = dist - iBin;
  }
 
  if (iBin < 0 || iBin >= Nu) {
    iBin = -1;
    frac = 0.0;
  } else if (iBin >= Nu - 1) {
    // Last bin: the +1 neighbour used by the trilinear interpolation would
    // fall outside the axis, so clamp to the last valid cell and use a
    // fractional weight of 1 (use the upper-edge slice, don't read past
    // the end).
    iBin = Nu - 2;
    frac = 1.0;
  }
 
  return {iBin, frac};
}
 
int MuonShieldFieldMap::GetMapBin(int iX, int iY, int iZ) const {
  int index = (iX * fNy + iY) * fNz + iZ;
  if (index < 0) {
    index = 0;
  } else if (index >= fN) {
    index = fN - 1;
  }
  return index;
}
 
double MuonShieldFieldMap::TriLinearInterp(double A, double B, double C,
                                          double D, double E, double F,
                                          double G, double H, double xFrac,
                                          double xFrac1, double yFrac,
                                          double yFrac1, double zFrac,
                                          double zFrac1) {
  const double F00 = A * xFrac1 + B * xFrac;
  const double F10 = C * xFrac1 + D * xFrac;
  const double F01 = E * xFrac1 + F * xFrac;
  const double F11 = G * xFrac1 + H * xFrac;
 
  const double F0 = F00 * yFrac1 + F10 * yFrac;
  const double F1 = F01 * yFrac1 + F11 * yFrac;
 
  return F0 * zFrac1 + F1 * zFrac;
}
 
void MuonShieldFieldMap::Field(double xCm, double yCm, double zCm,
                               double BkG[3]) const {
  // Local coords: z-translation only, matching addVMCFields()'s
  // defineFieldMap(..., TVector3(0,0,Entrance[0]), TVector3(0,0,0), ...)
  // -- zero rotation, zero x/y offset.
  double x = xCm;
  double y = yCm;
  const double z = zCm - fZOffsetCm;
 
  double BxSign = 1.0;
  double BzSign = 1.0;
 
  if (fSymmetric) {
    // The map stores only the x > 0, y > 0 quadrant; fold negative
    // coordinates in and track the induced signs (dipole symmetry: Bx is
    // antisymmetric in x and y, Bz antisymmetric in y only, By always
    // symmetric).
    if (x < 0.0) {
      x = -x;
      BxSign *= -1.0;
    }
    if (y < 0.0) {
      y = -y;
      BxSign *= -1.0;
      BzSign = -1.0;
    }
  }
 
  BkG[0] = 0.0;
  BkG[1] = 0.0;
  BkG[2] = 0.0;
 
  if (!InsideRange(x, y, z)) return;
 
  const BinInfo xBin = GetBinInfo(x, kX);
  const BinInfo yBin = GetBinInfo(y, kY);
  const BinInfo zBin = GetBinInfo(z, kZ);
  if (xBin.bin == -1 || yBin.bin == -1 || zBin.bin == -1) return;
 
  const int iX = xBin.bin, iY = yBin.bin, iZ = zBin.bin;
  const int iX1 = iX + 1, iY1 = iY + 1, iZ1 = iZ + 1;
 
  const auto& cA = fFieldMap[GetMapBin(iX, iY, iZ)];
  const auto& cB = fFieldMap[GetMapBin(iX1, iY, iZ)];
  const auto& cC = fFieldMap[GetMapBin(iX, iY1, iZ)];
  const auto& cD = fFieldMap[GetMapBin(iX1, iY1, iZ)];
  const auto& cE = fFieldMap[GetMapBin(iX, iY, iZ1)];
  const auto& cF = fFieldMap[GetMapBin(iX1, iY, iZ1)];
  const auto& cG = fFieldMap[GetMapBin(iX, iY1, iZ1)];
  const auto& cH = fFieldMap[GetMapBin(iX1, iY1, iZ1)];
 
  const double xFrac = xBin.frac, xFrac1 = 1.0 - xFrac;
  const double yFrac = yBin.frac, yFrac1 = 1.0 - yFrac;
  const double zFrac = zBin.frac, zFrac1 = 1.0 - zFrac;
 
  double Bc[3];
  for (int i = 0; i < 3; ++i) {
    Bc[i] = TriLinearInterp(cA[i], cB[i], cC[i], cD[i], cE[i], cF[i], cG[i],
                            cH[i], xFrac, xFrac1, yFrac, yFrac1, zFrac, zFrac1);
  }
 
  BkG[0] = Bc[0] * BxSign;
  BkG[1] = Bc[1];
  BkG[2] = Bc[2] * BzSign;
}
