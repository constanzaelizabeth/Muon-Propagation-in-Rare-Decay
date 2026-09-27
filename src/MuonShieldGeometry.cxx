#include "MuonShieldGeometry.hh"
 
#include <algorithm>
#include <cmath>
#include <vector>
 
#include "G4GenericTrap.hh"
#include "G4SystemOfUnits.hh"
 
namespace MuonShieldGeometry {
 
namespace {
 
// "anti_overlap" in ShipMuonShield::CreateMagnet(): a small gap left
// between adjacent field regions so Geant4/TGeo don't choke on exactly
// touching solids.
constexpr double kAntiOverlap = 0.1;  // cm
 
// ShipMuonShield::CreateMagnet():
//   gap = std::ceil(std::max(100. / dY, gap));
// (and the same for gap2/dY2). The result (renamed coil_gap/coil_gap2
// there) is what both cornersTL and cornersMainSideL actually use --
// despite the original using two different local names for it, it is
// numerically the exact same value in both formulas, since `gap`/`gap2`
// are overwritten in place before either shape is built.
double EffectiveGap(double gapParam, double dY) {
  return std::ceil(std::max(100.0 / dY, gapParam));
}
 
}  // namespace
 
MagnetShapes BuildMagnetShapes(const MuonShieldMagnetParams& p) {
  const double dX = p.dXIn;
  const double dY = p.dYIn;
  const double dX2 = p.dXOut;
  const double dY2 = p.dYOut;
  const double middleGap = p.midGapIn;
  const double middleGap2 = p.midGapOut;
  const double dYyoke1 = p.dYyokeIn;
  const double dYyoke2 = p.dYyokeOut;
  const double ratioYoke1 = p.ratioYokeIn;
  const double ratioYoke2 = p.ratioYokeOut;
  const double coilGap = EffectiveGap(p.gapIn, dY);
  const double coilGap2 = EffectiveGap(p.gapOut, dY2);
 
  MagnetShapes s;
 
  // ---   cornersMainL   (ShipMuonShield.cxx, verbatim)   ------------------
  s.mainL = {
      G4TwoVector(middleGap, -(dY + dYyoke1) - kAntiOverlap),
      G4TwoVector(middleGap, dY + dYyoke1 - kAntiOverlap),
      G4TwoVector(dX + middleGap, dY - kAntiOverlap),
      G4TwoVector(dX + middleGap, -(dY - kAntiOverlap)),
      G4TwoVector(middleGap2, -(dY2 + dYyoke2 - kAntiOverlap)),
      G4TwoVector(middleGap2, dY2 + dYyoke2 - kAntiOverlap),
      G4TwoVector(dX2 + middleGap2, dY2 - kAntiOverlap),
      G4TwoVector(dX2 + middleGap2, -(dY2 - kAntiOverlap)),
  };
 
  // ---   cornersTL   (verbatim; uses the effective/ceil'd gap)   ---------
  s.topLeft = {
      G4TwoVector(middleGap + dX, dY),
      G4TwoVector(middleGap, dY + dYyoke1),
      G4TwoVector(dX + ratioYoke1 * dX + middleGap + coilGap, dY + dYyoke1),
      G4TwoVector(dX + middleGap + coilGap, dY),
      G4TwoVector(middleGap2 + dX2, dY2),
      G4TwoVector(middleGap2, dY2 + dYyoke2),
      G4TwoVector(dX2 + ratioYoke2 * dX2 + middleGap2 + coilGap2,
                 dY2 + dYyoke2),
      G4TwoVector(dX2 + middleGap2 + coilGap2, dY2),
  };
 
  // ---   cornersMainSideL   (verbatim; same effective gap as above)   ----
  s.mainSideL = {
      G4TwoVector(dX + middleGap + coilGap, -dY),
      G4TwoVector(dX + middleGap + coilGap, dY),
      G4TwoVector(dX + ratioYoke1 * dX + middleGap + coilGap, dY + dYyoke1),
      G4TwoVector(dX + ratioYoke1 * dX + middleGap + coilGap,
                 -(dY + dYyoke1)),
      G4TwoVector(dX2 + middleGap2 + coilGap2, -dY2),
      G4TwoVector(dX2 + middleGap2 + coilGap2, dY2),
      G4TwoVector(dX2 + ratioYoke2 * dX2 + middleGap2 + coilGap2,
                 dY2 + dYyoke2),
      G4TwoVector(dX2 + ratioYoke2 * dX2 + middleGap2 + coilGap2,
                 -(dY2 + dYyoke2)),
  };
 
  return s;
}
 
Corners MirrorX(const Corners& c) {
  Corners out;
  for (std::size_t i = 0; i < 4; ++i) {
    out[3 - i] = G4TwoVector(-c[i].x(), c[i].y());
    out[7 - i] = G4TwoVector(-c[4 + i].x(), c[4 + i].y());
  }
  return out;
}
 
Corners MirrorY(const Corners& c) {
  Corners out;
  for (std::size_t i = 0; i < 4; ++i) {
    out[3 - i] = G4TwoVector(c[i].x(), -c[i].y());
    out[7 - i] = G4TwoVector(c[4 + i].x(), -c[4 + i].y());
  }
  return out;
}
 
Corners MirrorXY(const Corners& c) {
  Corners out;
  for (std::size_t i = 0; i < 8; ++i) {
    out[i] = G4TwoVector(-c[i].x(), -c[i].y());
  }
  return out;
}
 
G4GenericTrap* MakeTrap(const G4String& name, const Corners& corners,
                        double halfLenZcm) {
  std::vector<G4TwoVector> pts;
  pts.reserve(8);
  for (const auto& v : corners) {
    pts.emplace_back(v.x() * cm, v.y() * cm);
  }
  return new G4GenericTrap(name, halfLenZcm * cm, pts);
}
 
std::array<double, 8> MagnetCentersCm(double zStartCm) {
  std::array<double, 8> z{};
  for (std::size_t i = 0; i < kTry2026Magnets.size(); ++i) {
    const auto& p = kTry2026Magnets[i];
    if (i == 0) {
      z[i] = zStartCm + p.gapBeforeZ + p.halfLenZ;
    } else {
      const auto& prev = kTry2026Magnets[i - 1];
      z[i] = z[i - 1] + prev.halfLenZ + p.gapBeforeZ + p.halfLenZ;
    }
  }
  return z;
}
 
double MaxTransverseExtentCm() {
  double maxExtent = 0.0;
  for (const auto& p : kTry2026Magnets) {
    const MagnetShapes s = BuildMagnetShapes(p);
    const Corners allShapes[8] = {
        s.mainL,          MirrorXY(s.mainL),      s.mainSideL,
        MirrorXY(s.mainSideL), s.topLeft,          MirrorX(s.topLeft),
        MirrorY(s.topLeft),    MirrorXY(s.topLeft),
    };
    for (const auto& shape : allShapes) {
      for (const auto& v : shape) {
        maxExtent = std::max({maxExtent, std::fabs(v.x()), std::fabs(v.y())});
      }
    }
  }
  return maxExtent;
}
 
}  // namespace MuonShieldGeometry
 
