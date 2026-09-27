#ifndef MUONSHIELDGEOMETRY_HH
#define MUONSHIELDGEOMETRY_HH
 
#include <array>
 
#include "G4String.hh"
#include "G4TwoVector.hh"
#include "MuonShieldData.hh"
 
class G4GenericTrap;
 
// Ported 1:1 from ShipMuonShield::CreateMagnet()'s non-superconducting
// branch (see MuonShieldData.hh for why the SC_key branch is skipped) and
// ShipMuonShield::CreateArb8()/Initialize()'s Z-position recursion.
//
// G4GenericTrap is Geant4's equivalent of ROOT's TGeoArb8: both describe a
// prism via 8 (x,y) vertices -- indices 0..3 for the face at -halfZ,
// indices 4..7 for the face at +halfZ -- linearly interpolated in between.
// That equivalence is what makes this a faithful, not approximate, port:
// the corner formulas below are the same ones ShipMuonShield.cxx computes,
// just written once and handed to Geant4's own version of the same solid.
namespace MuonShieldGeometry {
 
// One shape's 8 vertices, in raw cm (see MuonShieldData.hh for the unit
// convention), same layout as TGeoArb8's flat 16-double `corners` array /
// G4GenericTrap's 8-point constructor argument.
using Corners = std::array<G4TwoVector, 8>;
 
struct MagnetShapes {
  Corners mainL;      // core, "left" half (ShipMuonShield's "_MiddleMagL")
  Corners mainSideL;  // return flux, "left" half ("_MagRetL")
  Corners topLeft;    // top-left yoke bar ("_MagTopLeft")
  // mainR/mainSideR/topRight/bottomLeft/bottomRight are exact mirrors (see
  // MirrorX/MirrorY/MirrorXY below) of the three shapes above, and always
  // carry the *same* field vector as their mirror partner -- see
  // BellMagneticField.cc for how that lets the field evaluate only these
  // three canonical solids per magnet instead of all eight.
};
 
// Builds the 3 canonical shapes for one magnet, in that magnet's own local
// (x, y) frame -- x_translation = y_translation = 0 in every CreateArb8()
// call in the original, so no transverse offset is needed here.
MagnetShapes BuildMagnetShapes(const MuonShieldMagnetParams& p);
 
// x -> -x, with vertex order reversed within each face to preserve
// winding (matches ShipMuonShield's cornersTL -> cornersTR transform
// exactly, reindex included). Use to get topRight from topLeft.
Corners MirrorX(const Corners& c);
 
// y -> -y, with the same reindexing as MirrorX (derived from
// ShipMuonShield's cornersTR -> cornersBL transform, which composes a
// direct negate with the TL->TR reindex). Use to get bottomLeft from
// topLeft.
Corners MirrorY(const Corners& c);
 
// (x, y) -> (-x, -y), vertex order UNCHANGED (a double flip preserves
// winding, matching ShipMuonShield's direct std::negate transforms for
// cornersMainR, cornersMainSideR and cornersBR). Use to get mainR from
// mainL, mainSideR from mainSideL, and bottomRight from topLeft.
Corners MirrorXY(const Corners& c);
 
// Builds a placeable G4GenericTrap ("Arb8" equivalent) from 8 raw-cm
// corners and a raw-cm half-length along z; converts to G4-internal units
// internally (this is the only place that conversion happens for the
// shield's transverse shapes).
G4GenericTrap* MakeTrap(const G4String& name, const Corners& corners,
                        double halfLenZcm);
 
// Z (world-local, i.e. already including whatever z-offset the caller
// wants magnet 0 measured from) of each of the 8 magnets' centers, in raw
// cm, following ShipMuonShield::Initialize()'s recursion:
//   Z[0] = zStartCm + dZ[0] + halfLenZ[0]
//   Z[i] = Z[i-1] + halfLenZ[i-1] + dZ[i] + halfLenZ[i]   (i > 0)
std::array<double, 8> MagnetCentersCm(double zStartCm);
 
// Largest |x| or |y| (raw cm) reached by any of the 8 magnets' full shapes
// (all 8 mirrored pieces, not just the 3 canonical ones) -- used to size
// the World volume so it actually contains the shield.
double MaxTransverseExtentCm();
 
}  // namespace MuonShieldGeometry
 
#endif  // MUONSHIELDGEOMETRY_HH
 
