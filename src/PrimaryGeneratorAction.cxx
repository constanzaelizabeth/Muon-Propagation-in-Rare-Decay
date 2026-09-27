#include "PrimaryGeneratorAction.hh"
 
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <unordered_map>
 
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
 
namespace {
constexpr int kHnlPdg = 9900015;
constexpr int kMuPdg = 13;
constexpr double kCmToMm = 10.0;
constexpr double kC_mm_per_ns = 299.792458;  // matches G4's c_light exactly
 
std::vector<std::string> splitCsv(const std::string& line) {
  std::vector<std::string> fields;
  std::string field;
  bool quoted = false;
  for (size_t i = 0; i < line.size(); ++i) {
    const char c = line[i];
    if (c == '"') {
      if (quoted && i + 1 < line.size() && line[i + 1] == '"') {
        field.push_back('"');
        ++i;
      } else {
        quoted = !quoted;
      }
    } else if (c == ',' && !quoted) {
      fields.push_back(field);
      field.clear();
    } else {
      field.push_back(c);
    }
  }
  fields.push_back(field);
  return fields;
}
 
double asDouble(const std::vector<std::string>& row,
               const std::unordered_map<std::string, size_t>& column,
               const std::string& name) {
  const auto it = column.find(name);
  if (it == column.end() || it->second >= row.size()) {
    throw std::runtime_error("Falta la columna obligatoria '" + name + "'");
  }
  return std::stod(row[it->second]);
}
 
long long asInteger(const std::vector<std::string>& row,
                    const std::unordered_map<std::string, size_t>& column,
                    const std::string& name) {
  return static_cast<long long>(asDouble(row, column, name));
}
}  // namespace
 
PrimaryGeneratorAction::PrimaryGeneratorAction(const std::string& csvPath, G4double hnlMassGeV,
                                              G4int seed, G4double targetThicknessCm,
                                              G4double interactionLengthCm,
                                              SharedEventInfo* shared)
    : fUniform01(1e-12, 1.0),
      fTargetThickness_cm(targetThicknessCm),
      fInteractionLength_cm(interactionLengthCm),
      fShared(shared) {
  fVertexRng.seed(static_cast<uint64_t>(seed) ^ 0x9E3779B97f4A7C15ULL);
 
  fPythia.readString("ProcessLevel:all = off");
  fPythia.readString("PartonLevel:all = off");
  fPythia.readString("HadronLevel:all = off");
  fPythia.readString("Random:setSeed = on");
  fPythia.readString("Random:seed = " + std::to_string(seed));
 
  fPythia.readString(std::to_string(kHnlPdg) + ":new = N2 N2 2 0 0 " +
                     std::to_string(hnlMassGeV) + " 0.0 0.0 0.0 0.0");
  fPythia.readString(std::to_string(kHnlPdg) + ":mayDecay = off");
 
  const double mB = fPythia.particleData.m0(521);
  const double mMu = fPythia.particleData.m0(13);
  fCtau0B_mm = fPythia.particleData.tau0(521);
  fBPlusValido = (hnlMassGeV < mB - mMu);
 
  if (fBPlusValido) {
    fPythia.readString("521:onMode = off");
    fPythia.readString("521:oneChannel = 1 1.0 0 -13 9900015");
  } else {
    G4cerr << "PrimaryGeneratorAction: mHNL = " << hnlMassGeV
          << " GeV no es cinematicamente valido para B+/- -> mu N2 (mB - mMu = "
          << (mB - mMu) << " GeV). No se generaran eventos." << G4endl;
  }
 
  if (!fPythia.init()) {
    G4cerr << "PrimaryGeneratorAction: pythia.init() fallo." << G4endl;
    return;
  }
 
  if (!LoadCsv(csvPath) || fRows.empty()) {
    G4cerr << "PrimaryGeneratorAction: no se cargaron filas B+/B- validas de " << csvPath
          << G4endl;
    return;
  }
 
  fReady = true;
  G4cout << "PrimaryGeneratorAction: " << fRows.size() << " candidatos B+/B- cargados de "
        << csvPath << " (mHNL = " << hnlMassGeV << " GeV)" << G4endl;
}
 
bool PrimaryGeneratorAction::LoadCsv(const std::string& csvPath) {
  std::ifstream input(csvPath);
  if (!input) {
    G4cerr << "PrimaryGeneratorAction: no se pudo abrir " << csvPath << G4endl;
    return false;
  }
  std::string headerLine;
  if (!std::getline(input, headerLine)) return false;
  if (!headerLine.empty() && headerLine.back() == '\r') headerLine.pop_back();
 
  const auto headers = splitCsv(headerLine);
  std::unordered_map<std::string, size_t> column;
  for (size_t i = 0; i < headers.size(); ++i) column[headers[i]] = i;
 
  std::string line;
  long long inputRows = 0;
  while (std::getline(input, line)) {
    ++inputRows;
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) continue;
    try {
      const auto fields = splitCsv(line);
      const int id = static_cast<int>(asInteger(fields, column, "id"));
      if (std::abs(id) != 521) continue;
      if (!fBPlusValido) continue;
 
      BRow row;
      row.evento = asInteger(fields, column, "evento");
      row.id = id;
      row.m = asDouble(fields, column, "m");
      row.px = asDouble(fields, column, "px");
      row.py = asDouble(fields, column, "py");
      row.pz = asDouble(fields, column, "pz");
      double energy = asDouble(fields, column, "e");
      const double onShellEnergy = std::sqrt(row.px * row.px + row.py * row.py +
                                             row.pz * row.pz + row.m * row.m);
      if (!std::isfinite(energy) || energy < onShellEnergy) energy = onShellEnergy;
      row.e = energy;
      row.xProd_mm = asDouble(fields, column, "xProd") * kCmToMm;
      row.yProd_mm = asDouble(fields, column, "yProd") * kCmToMm;
      row.tProd_mmc = asDouble(fields, column, "tProd") * kC_mm_per_ns;
      row.peso = asDouble(fields, column, "peso");
      fRows.push_back(row);
    } catch (const std::exception& error) {
      G4cerr << "PrimaryGeneratorAction: se omitio la fila " << inputRows + 1 << ": "
            << error.what() << G4endl;
    }
  }
  return true;
}
 
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
  if (!fReady) {
    G4cerr << "PrimaryGeneratorAction: no inicializado, evento vacio." << G4endl;
    return;
  }
 
  if (!fMuMinus || !fMuPlus || !fHnl) {
    G4ParticleTable* table = G4ParticleTable::GetParticleTable();
    fMuMinus = table->FindParticle("mu-");
    fMuPlus = table->FindParticle("mu+");
    fHnl = table->FindParticle(kHnlPdg);
  }
 
  for (size_t attempt = 0; attempt < fRows.size(); ++attempt) {
    if (fNextRow >= fRows.size()) {
      fNextRow = 0;
      G4cout << "PrimaryGeneratorAction: fin del CSV, se reinicia." << G4endl;
    }
    const BRow& row = fRows[fNextRow++];
 
    // Depth at which the B actually interacts/is produced inside the
    // tungsten target (exponential absorption profile, same as before).
    const double U0 = fUniform01(fVertexRng);
    const double L = fTargetThickness_cm;
    const double lambda = fInteractionLength_cm;
    const double factor = 1.0 - std::exp(-L / lambda);
    const double zProd_cm = -lambda * std::log(1.0 - U0 * factor);
    const double zProd_mm = zProd_cm * kCmToMm;
 
    fPythia.event.reset();
    const int iMadre = fPythia.event.append(row.id, 23, 0, 0, 0, 0, 0, 0, row.px, row.py,
                                            row.pz, row.e, row.m);
    fPythia.event[iMadre].vProd(row.xProd_mm, row.yProd_mm, zProd_mm, row.tProd_mmc);
 
    // Sample how far the B itself flies before decaying, fixing the
    // actual B+/- -> mu + HNL vertex (same model as before).
    {
      Pythia8::Particle& p = fPythia.event[iMadre];
      const double px = p.px(), py = p.py(), pz = p.pz();
      const double pabs = std::sqrt(px * px + py * py + pz * pz);
      const double m = p.m();
      const double e = p.e();
      if (pabs > 0.0 && m > 0.0 && fCtau0B_mm > 0.0) {
        const double gammaBeta = pabs / m;
        const double beta = pabs / e;
        const double U = fUniform01(fVertexRng);
        const double properLength_mm = -fCtau0B_mm * std::log(U);
        const double labLength_mm = gammaBeta * properLength_mm;
        const double nx = px / pabs, ny = py / pabs, nz = pz / pabs;
        p.vProd(p.xProd() + labLength_mm * nx, p.yProd() + labLength_mm * ny,
               p.zProd() + labLength_mm * nz, p.tProd() + labLength_mm / beta);
      }
    }
 
    if (!fPythia.moreDecays()) continue;
 
    const int d1 = fPythia.event[iMadre].daughter1();
    const int d2 = fPythia.event[iMadre].daughter2();
    if (d1 <= 0 || d2 < d1) continue;
 
    int iMuon = -1, iHnl = -1;
    for (int i = d1; i <= d2; ++i) {
      if (std::abs(fPythia.event[i].id()) == kMuPdg) {
        iMuon = i;
      } else if (std::abs(fPythia.event[i].id()) == kHnlPdg) {
        iHnl = i;
      }
    }
    if (iMuon < 0 || iHnl < 0) continue;
 
    const Pythia8::Particle& mu = fPythia.event[iMuon];
    const Pythia8::Particle& hnl = fPythia.event[iHnl];
 
    // Pythia8 lengths are already mm and its "time" branch is c*t in mm;
    // G4's c_light is also expressed in mm/ns, so dividing by it gives ns
    // directly -- exactly the units G4PrimaryVertex wants.
    auto* vertex = new G4PrimaryVertex(mu.xProd() * mm, mu.yProd() * mm, mu.zProd() * mm,
                                       mu.tProd() / c_light);
 
    G4ParticleDefinition* muDef = (mu.id() == kMuPdg) ? fMuMinus : fMuPlus;
    auto* muPrim =
        new G4PrimaryParticle(muDef, mu.px() * GeV, mu.py() * GeV, mu.pz() * GeV);
    muPrim->SetMass(mu.m() * GeV);
    vertex->SetPrimary(muPrim);
 
    auto* hnlPrim =
        new G4PrimaryParticle(fHnl, hnl.px() * GeV, hnl.py() * GeV, hnl.pz() * GeV);
    hnlPrim->SetMass(hnl.m() * GeV);
    vertex->SetPrimary(hnlPrim);
 
    event->AddPrimaryVertex(vertex);
 
    if (fShared) {
      fShared->csvEvento = row.evento;
      fShared->weight = row.peso;
    }
    return;
  }
 
  G4cerr << "PrimaryGeneratorAction: se agotaron los reintentos sin producir un "
        << "decaimiento B -> mu HNL valido." << G4endl;
}
 
