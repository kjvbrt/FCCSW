#include "UpstreamDownstreamMaterial.h"

// FCCSW
#include "DetInterface/IGeoSvc.h"

// datamodel
#include "datamodel/CaloHitCollection.h"
#include "datamodel/MCParticleCollection.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "GaudiKernel/ITHistSvc.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TVector2.h"

// DD4hep
#include "DD4hep/Detector.h"
#include "DD4hep/Readout.h"


DECLARE_COMPONENT(UpstreamDownstreamMaterial)


UpstreamDownstreamMaterial::UpstreamDownstreamMaterial(const std::string& aName, ISvcLocator* aSvcLoc)
    : GaudiAlgorithm(aName, aSvcLoc), m_geoSvc("GeoSvc", aName) {
  declareProperty("deposits", m_deposits, "Energy deposits (input)");
  declareProperty("particle", m_particle, "Generated single-particle event (input)");
}


UpstreamDownstreamMaterial::~UpstreamDownstreamMaterial() {}


StatusCode UpstreamDownstreamMaterial::initialize() {
  if (GaudiAlgorithm::initialize().isFailure()) return StatusCode::FAILURE;

  // Check geometry service
  if (!m_geoSvc) {
    error() << "Unable to locate Geometry Service. "
            << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
    return StatusCode::FAILURE;
  }

  // Check if readouts exist
  if (m_geoSvc->lcdd()->readouts().find(m_readoutName) == m_geoSvc->lcdd()->readouts().end()) {
    error() << "Readout <<" << m_readoutName << ">> does not exist." << endmsg;
    return StatusCode::FAILURE;
  }

  // Check segmentation
  auto segmentation = m_geoSvc->lcdd()->readout(m_readoutName).segmentation().segmentation();
  if (segmentation) {
    verbose() << "Segmentation name: " << segmentation->name() << endmsg;
  } else {
    warning() << "Segmentation not found, hit position histograms won't get filled!" << endmsg;
  }

  // Check histogram service
  m_histSvc = service("THistSvc");
  if (!m_histSvc) {
    error() << "Unable to locate Histogram Service" << endmsg;
    return StatusCode::FAILURE;
  }

  // Define and register histograms
  for (uint i = 0; i < m_numLayers; i++) {
    m_cellEnergyPhi.push_back(new TH1F(("upstreamEnergy_phi" + std::to_string(i)).c_str(),
                                       ("Energy deposited in layer " + std::to_string(i)).c_str(), 1000, -m_phi,
                                       m_phi));
    m_cellEnergyPhi.back()->Sumw2();
    if (m_histSvc->regHist("/det/upstreamEnergy_phi" + std::to_string(i), m_cellEnergyPhi.back()).isFailure()) {
      error() << "Couldn't register histogram" << endmsg;
      return StatusCode::FAILURE;
    }

    m_gUpstreamEnergyCellEnergy.emplace_back(new TGraph());
    m_gUpstreamEnergyCellEnergy.back()->SetName(("upstreamEnergy_presamplerEnergy_graph_" + std::to_string(i)).c_str());

    // Downstream energy vs. energy in every layer
    std::string histName = "downstreamEnergy_energyInLayer_" + std::to_string(i);
    std::string histTitle = "Downstream energy vs energy deposited in layer id: " + std::to_string(i);
    histTitle += ";Energy in layer id: " + std::to_string(i) + " [GeV];Downstream Energy [GeV]";
    m_downstreamEnergyEnergyInLayer.emplace_back(new TH2F(histName.c_str(), histTitle.c_str(), 50, 0., 0., 50, 0., 0.));
    m_downstreamEnergyEnergyInLayer.back()->Sumw2();
    if (m_histSvc->regHist("/det/" + histName, m_downstreamEnergyEnergyInLayer.back()).isFailure()) {
      error() << "Couldn't register histogram \"" << histName << "\"!" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  m_hEnergyInLayers = new TH1F("energyInLayers", "Energy deposited in layer ",
                               m_numLayers, 0, m_numLayers);
  m_hEnergyInLayers->Sumw2();
  if (m_histSvc->regHist("/det/energyInLayers", m_hEnergyInLayers).isFailure()) {
    error() << "Couldn't register histogram \"energyInLayers\"!" << endmsg;
    return StatusCode::FAILURE;
  }

  m_hSumEinLayers = new TH1F("sumEinLayers",
                             "Sum of energy deposited in all layers;E [GeV];N_{evt}",
                             200, 0., 0.);
  m_hSumEinLayers->Sumw2();
  if (m_histSvc->regHist("/det/sumEinLayers", m_hSumEinLayers).isFailure()) {
    error() << "Couldn't register histogram \"sumEinLayers\"!" << endmsg;
    return StatusCode::FAILURE;
  }

  m_hEnergyInCaloAndCryo = new TH1F("energyInCaloAndCryo",
                                    "Energy deposited in calorimeter and in cryostat;E [GeV];N_{evt}",
                                    200, 0., 0.);
  m_hEnergyInCaloAndCryo->Sumw2();
  if (m_histSvc->regHist("/det/energyInCaloAndCryo", m_hEnergyInCaloAndCryo).isFailure()) {
    error() << "Couldn't register histogram \"EnergyInCaloAndCryo\"!" << endmsg;
    return StatusCode::FAILURE;
  }

  m_hEnergyInCryo = new TH1F("energyInCryo",
                             "Energy deposited in cryostat;E [GeV];N_{evt}",
                             200, 0., 0.);
  m_hEnergyInCryo->Sumw2();
  if (m_histSvc->regHist("/det/energyInCryo", m_hEnergyInCryo).isFailure()) {
    error() << "Couldn't register histogram \"EnergyInCryo\"!" << endmsg;
    return StatusCode::FAILURE;
  }

  m_hEnergyInCryoFront = new TH1F("energyInCryoFront",
                                  "Energy deposited in cryostat front;E [GeV];N_{evt}",
                                  200, 0., 0.);
  m_hEnergyInCryoFront->Sumw2();
  if (m_histSvc->regHist("/det/energyInCryoFront", m_hEnergyInCryoFront).isFailure()) {
    error() << "Couldn't register histogram \"EnergyInCryoFront\"!" << endmsg;
    return StatusCode::FAILURE;
  }

  m_hEnergyInCryoBack = new TH1F("energyInCryoBack",
                                  "Energy deposited in cryostat back;E [GeV];N_{evt}",
                                  200, 0., 0.);
  m_hEnergyInCryoBack->Sumw2();
  if (m_histSvc->regHist("/det/energyInCryoBack", m_hEnergyInCryoBack).isFailure()) {
    error() << "Couldn't register histogram \"EnergyInCryoBack\"!" << endmsg;
    return StatusCode::FAILURE;
  }

  m_hEnergyInCryoSides = new TH1F("energyInCryoSides",
                                  "Energy deposited in cryostat sides;E [GeV];N_{evt}",
                                  200, 0., 0.);
  m_hEnergyInCryoSides->Sumw2();
  if (m_histSvc->regHist("/det/energyInCryoSides", m_hEnergyInCryoSides).isFailure()) {
    error() << "Couldn't register histogram \"EnergyInCryoSides\"!" << endmsg;
    return StatusCode::FAILURE;
  }

  m_hEnergyInCryoLArBathFront = new TH1F("energyInCryoLArBathFront",
                                  "Energy deposited in cryostat LAr bath front;E [GeV];N_{evt}",
                                  200, 0., 0.);
  m_hEnergyInCryoLArBathFront->Sumw2();
  if (m_histSvc->regHist("/det/energyInCryoLArBathFront", m_hEnergyInCryoLArBathFront).isFailure()) {
    error() << "Couldn't register histogram \"EnergyInCryoLArBathFront\"!" << endmsg;
    return StatusCode::FAILURE;
  }

  m_hEnergyInCryoLArBathBack = new TH1F("energyInCryoLArBathBack",
                                  "Energy deposited in cryostat LAr bath back;E [GeV];N_{evt}",
                                  200, 0., 0.);
  m_hEnergyInCryoLArBathBack->Sumw2();
  if (m_histSvc->regHist("/det/energyInCryoLArBathBack", m_hEnergyInCryoLArBathBack).isFailure()) {
    error() << "Couldn't register histogram \"EnergyInCryoLArBathBack\"!" << endmsg;
    return StatusCode::FAILURE;
  }

  m_hParticleMomentumXY = new TH2F("particleMomentumXY",
                                   "Particle momentum in XY plane;p_{x};p_{y}",
                                   200, 0., 0., 200, 0., 0.);
  if (m_histSvc->regHist("/det/particleMomentumXY", m_hParticleMomentumXY).isFailure()) {
    error() << "Couldn't register histogram \"particleMomentumXY\"!" << endmsg;
    return StatusCode::FAILURE;
  }

  m_hParticleMomentumZY = new TH2F("hParticleMomentumZY",
                                   "Particle momentum in ZY plane;p_{z};p_{y}",
                                   200, 0., 0., 200, 0., 0.);
  if (m_histSvc->regHist("/det/ParticleMomentumZY", m_hParticleMomentumZY).isFailure()) {
    error() << "Couldn't register histogram \"hParticleMomentumZY\"!" << endmsg;
    return StatusCode::FAILURE;
  }

  m_hHitPositionXY = new TH2F("hHitPositionXY",
                                   "Hit position in XY plane;p_{x};p_{y}",
                                   200, -3., 3., 200, -3., 3.);
  if (m_histSvc->regHist("/det/HitPositionXY", m_hHitPositionXY).isFailure()) {
    error() << "Couldn't register histogram: hHitPositionXY" << endmsg;
    return StatusCode::FAILURE;
  }

  m_hHitPositionZY = new TH2F("hHitPositionZY",
                                   "Hit position in ZY plane;p_{z};p_{y}",
                                   200, -3., 3., 200, -3., 3.);
  if (m_histSvc->regHist("/det/HitPositionZY", m_hHitPositionZY).isFailure()) {
    error() << "Couldn't register histogram: hHitPositionZY" << endmsg;
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}


StatusCode UpstreamDownstreamMaterial::execute() {
  verbose() << "Event Number: " << m_gUpstreamEnergyCellEnergy.front()->GetN() << endmsg;

  auto decoder = m_geoSvc->lcdd()->readout(m_readoutName).idSpec().decoder();
  auto segmentation = m_geoSvc->lcdd()->readout(m_readoutName).segmentation().segmentation();

  // first check MC phi angle
  const auto particle = m_particle.get();
  double phi = 0;
  for (const auto& part : *particle) {
    phi = atan2(part.core().p4.py, part.core().p4.px);
    m_hParticleMomentumXY->Fill(part.core().p4.px, part.core().p4.py);
    m_hParticleMomentumZY->Fill(part.core().p4.pz, part.core().p4.py);
  }

  // Get the energy deposited in the cryostat(s) and in the detector (each layer)
  const auto deposits = m_deposits.get();
  std::vector<float> sumEinLayer;
  sumEinLayer.assign(m_numLayers, 0);
  float sumEinCryo = 0.;
  float sumEinCryoFront = 0.;
  float sumEinCryoBack = 0.;
  float sumEinCryoSides = 0.;
  float sumEinCryoLArBathFront = 0.;
  float sumEinCryoLArBathBack = 0.;
  for (const auto& hit: *deposits) {
    dd4hep::DDSegmentation::CellID cellId = hit.core().cellId;
    size_t cryoId = decoder->get(cellId, "cryo");
    if (cryoId == 0) {
      size_t layerId = decoder->get(cellId, "layer");
      sumEinLayer[layerId] += hit.core().energy;
    } else {
      sumEinCryo += hit.core().energy;
      size_t cryoTypeId = decoder->get(cellId, "type");
      switch(cryoTypeId) {
        case 1: sumEinCryoFront += hit.core().energy;
                break;
        case 2: sumEinCryoBack += hit.core().energy;
                break;
        case 3: sumEinCryoSides += hit.core().energy;
                break;
        case 4: sumEinCryoLArBathFront += hit.core().energy;
                break;
        case 5: sumEinCryoLArBathBack += hit.core().energy;
                break;
        default: warning() << "Wrong cryostat type ID found!\ncryoTypeId: " << cryoTypeId << endmsg;
      }
    }
  }

  // Get position of the deposits
  if (segmentation) {
    for (const auto& hit: *deposits) {
      dd4hep::DDSegmentation::CellID cellId = hit.core().cellId;
      dd4hep::DDSegmentation::Vector3D position = segmentation->position(cellId);
      m_hHitPositionXY->Fill(position.x(), position.y());
      m_hHitPositionZY->Fill(position.z(), position.y());
    }
  }

  // Calibrate energy in the calorimeter layers
  for (size_t i = 0; i < m_numLayers; ++i) {
    sumEinLayer[i] /= m_samplingFraction[i];
  }

  // Sum energy deposited in all calorimeter layers
  float sumEinCalo = 0.;
  for (size_t i = 0; i < m_numLayers; ++i) {
    sumEinCalo += sumEinLayer[i];
  }

  // Fill histograms and graphs
  for (size_t i = 0; i < m_numLayers; ++i) {
    m_cellEnergyPhi[i]->Fill(phi, sumEinLayer[i]);
    m_hEnergyInLayers->Fill(i, sumEinLayer[i]);
    m_gUpstreamEnergyCellEnergy.at(i)->SetPoint(m_gUpstreamEnergyCellEnergy.at(i)->GetN(),
                                                sumEinLayer[i], sumEinCryoFront + sumEinCryoLArBathFront);
    m_downstreamEnergyEnergyInLayer.at(i)->Fill(sumEinLayer[i], sumEinCryoBack + sumEinCryoLArBathBack);
    verbose() << "Energy deposited in layer " << i << ": " << sumEinLayer[i] << " GeV" << endmsg;
  }
  m_hSumEinLayers->Fill(sumEinCalo);
  m_hEnergyInCaloAndCryo->Fill(sumEinCalo + sumEinCryo);
  verbose() << "Sum of energy deposited in all calorimeter layers: " << sumEinCalo << " GeV" << endmsg;

  m_hEnergyInCryo->Fill(sumEinCryo);
  m_hEnergyInCryoFront->Fill(sumEinCryoFront);
  m_hEnergyInCryoBack->Fill(sumEinCryoBack);
  m_hEnergyInCryoSides->Fill(sumEinCryoSides);
  m_hEnergyInCryoLArBathFront->Fill(sumEinCryoLArBathFront);
  m_hEnergyInCryoLArBathBack->Fill(sumEinCryoLArBathBack);
  verbose() << "Energy deposited in the cryostat: " << sumEinCryo << " GeV" << endmsg;
  verbose() << "Energy deposited in the cryostat front: " << sumEinCryoFront << " GeV" << endmsg;
  verbose() << "Energy deposited in the cryostat back: " << sumEinCryoBack<< " GeV" << endmsg;
  verbose() << "Energy deposited in the cryostat sides: " << sumEinCryoSides << " GeV" << endmsg;
  verbose() << "Energy deposited in the cryostat LAr bath front: " << sumEinCryoLArBathFront << " GeV" << endmsg;
  verbose() << "Energy deposited in the cryostat LAr bath back: " << sumEinCryoLArBathBack << " GeV" << endmsg;
  verbose() << "Energy deposited in the calorimeter and in the cryostat: " << sumEinCalo + sumEinCryo << " GeV" << endmsg;

  return StatusCode::SUCCESS;
}


StatusCode UpstreamDownstreamMaterial::finalize() {
  size_t nEvt = m_gUpstreamEnergyCellEnergy.front()->GetN();
  size_t nBin = std::sqrt(nEvt) + 1;
  verbose() << "Number of generated events: " << nEvt << endmsg;
  verbose() << "Number of bins: " << nBin << endmsg;

  for (size_t i = 0; i < m_gUpstreamEnergyCellEnergy.size(); ++i) {
    double xMax = *std::max_element(m_gUpstreamEnergyCellEnergy.at(i)->GetX(),
                                    m_gUpstreamEnergyCellEnergy.at(i)->GetX() + nEvt);
    double yMax = *std::max_element(m_gUpstreamEnergyCellEnergy.at(i)->GetY(),
                                    m_gUpstreamEnergyCellEnergy.at(i)->GetY() + nEvt);
    xMax = 1.000001*xMax;
    yMax = 1.000001*yMax;

    std::string histName = "upstreamEnergy_presamplerEnergy_" + std::to_string(i);
    m_upstreamEnergyCellEnergy.emplace_back(
        new TH2F(histName.c_str(),
                 ("Upstream energy vs energy deposited in layer " + std::to_string(i)).c_str(),
                 nBin, 0., xMax, nBin, 0., yMax));
    m_upstreamEnergyCellEnergy.back()->Sumw2();
    if (m_histSvc->regHist("/det/" + histName, m_upstreamEnergyCellEnergy.back()).isFailure()) {
      error() << "Couldn't register histogram \"" << histName << "\"!" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  for (size_t i = 0; i < m_upstreamEnergyCellEnergy.size(); ++i) {
    GraphToHist(m_gUpstreamEnergyCellEnergy.at(i), m_upstreamEnergyCellEnergy.at(i));
  }

  return GaudiAlgorithm::finalize();
}


void UpstreamDownstreamMaterial::GraphToHist(TGraph* graph, TH2F* hist) {
  for (int i = 0; i < graph->GetN(); ++i) {
    hist->Fill(graph->GetPointX(i), graph->GetPointY(i));
  }
}
