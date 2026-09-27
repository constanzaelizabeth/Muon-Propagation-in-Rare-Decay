#ifndef MUONSHIELDDATA_HH
#define MUONSHIELDDATA_HH
 
#include <array>
 
// Ported 1:1 from FairShip's passive/ShipMuonShield.{h,cxx} and the
// "TRY_2026" preset in python/geometry_config.py
// (shield_db["TRY_2026"]["params"] / ["hybrid"] / ["WithConstField"]).
//
// IMPORTANT -- units: every number here is a raw magnitude in the same
// convention ShipMuonShield.cxx itself uses internally (ShipUnit "cm", the
// native TGeo unit -- i.e. these are plain cm, NOT yet Geant4-internal
// units). Multiply lengths by G4's `cm` and BgoalTesla by G4's `tesla`
// only at the point where a G4 solid / field vector is actually built
// (see MuonShieldGeometry.hh and BellMagneticField.cc); keep every
// intermediate computation (the corner formulas, the effective-gap
// ceiling, the Z-position recursion) in this same raw-cm domain, exactly
// as the original does, to avoid subtle unit-conversion drift relative to
// FairShip's own numbers.
//
// shield_db["TRY_2026"]["hybrid"] == False, so ShipMuonShield's SC_key is
// always false for this preset (SC_key = (nM==3 && fSC_mag), and
// fSC_mag = hybrid = False) -- the superconducting-magnet branch of
// ShipMuonShield::CreateMagnet() is therefore never exercised by TRY_2026
// and is intentionally NOT ported here.
//
// shield_db["TRY_2026"]["WithConstField"] == False: in real FairShip this
// means ShipMuonShield itself never calls TGeoVolume::SetField() on these
// volumes: the ACTIVE field is wired up separately (via ShipFieldMaker +
// an external text config we don't have, see python/geomGeant4.py). What
// we port here is the field ShipMuonShield::ConstructGeometry() computes
// for each of the 8 sub-volumes per magnet (Bgoal, direction, and the 4
// field vectors) -- i.e. the physical field the shield is *designed* to
// have, independent of which particular runtime mechanism activates it.
 
enum class MuonShieldFieldDirection { kUp, kDown };
 
struct MuonShieldMagnetParams {
  double gapBeforeZ;  // "dZ" column: gap before this magnet, along z (cm)
  double halfLenZ;    // "Z_rel" column: this magnet's own half-length (cm)
  double dXIn, dXOut;
  double dYIn, dYOut;
  double gapIn, gapOut;
  double ratioYokeIn, ratioYokeOut;
  double dYyokeIn, dYyokeOut;
  double midGapIn, midGapOut;
  double BgoalTesla;  // signed; can be 0 (magnet 5) or negative (6, 7)
  MuonShieldFieldDirection dir;
};
 
// Column order matches Initialize(): dZ, Z_rel, dXIn, dXOut, dYIn, dYOut,
// gapIn, gapOut, ratio_yokesIn, ratio_yokesOut, dY_yokeIn, dY_yokeOut,
// midGapIn, midGapOut, Bgoal. Row 0 is "MagnAbsorb" (the hadron absorber,
// which ShipMuonShield builds with the exact same CreateMagnet() field
// machinery as every other magnet). fieldDirection is ShipMuonShield's
// hardcoded {up,up,up,up,up,down,down,down} (it has exactly 8 entries, one
// per TRY_2026 magnet, so no fallback/alternation logic is needed here).
inline constexpr std::array<MuonShieldMagnetParams, 8> kTry2026Magnets = {{
    {0.0, 115.0, 40.0, 40.0, 119.0, 119.0, 61.5, 61.5, 1.5375, 1.5375, 50.0,
     50.0, 0.0, 0.0, 2.05, MuonShieldFieldDirection::kUp},
    {15.0, 150.0, 62.0, 62.0, 11.0, 11.0, 12.0, 12.0, 1.0, 1.0, 69.0, 69.0,
     0.0, 0.0, 1.87, MuonShieldFieldDirection::kUp},
    {15.0, 225.0, 70.0, 70.0, 12.0, 12.0, 9.0, 9.0, 1.0, 1.0, 77.0, 77.0, 0.0,
     0.0, 1.865, MuonShieldFieldDirection::kUp},
    {16.0, 225.0, 64.0, 64.0, 14.0, 14.0, 8.0, 8.0, 1.3, 1.3, 71.0, 71.0, 0.0,
     0.0, 1.92, MuonShieldFieldDirection::kUp},
    {15.0, 225.0, 42.0, 42.0, 13.0, 13.0, 8.0, 8.0, 2.6, 2.6, 47.0, 47.0, 0.0,
     0.0, 1.89, MuonShieldFieldDirection::kUp},
    {16.0, 131.0, 40.0, 40.0, 30.0, 30.0, 0.0, 0.0, 2.25, 2.25, 22.0, 22.0,
     0.0, 0.0, 0.0, MuonShieldFieldDirection::kDown},
    {15.0, 168.0, 20.0, 30.0, 30.0, 30.0, 22.0, 12.0, 5.4, 3.6, 33.0, 33.0,
     0.0, 0.0, -1.42, MuonShieldFieldDirection::kDown},
    {19.0, 200.0, 38.0, 69.0, 51.0, 51.0, 8.0, 8.0, 3.2894736842105265,
     1.36231884057971, 76.0, 76.0, 0.0, 0.0, -2.035,
     MuonShieldFieldDirection::kDown},
}};
 
#endif  // MUONSHIELDDATA_HH
 
