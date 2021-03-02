#include "UpstreamDownstreamMaterial.h"

// FCCSW
#include "DetInterface/IGeoSvc.h"

// datamodel
#include "datamodel/CaloHitCollection.h"
#include "datamodel/MCParticleCollection.h"

// Gaudi
#include "CLHEP/Vector/ThreeVector.h"
#include "GaudiKernel/ITHistSvc.h"
#include "GaudiKernel/GaudiException.h"

// ROOT
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

  // Check histogram service
  m_histSvc = service("THistSvc");
  if (!m_histSvc) {
    error() << "Unable to locate Histogram Service" << endmsg;
    return StatusCode::FAILURE;
  }

  // Get intended number of events
  auto appMgr = service<IProperty>("ApplicationMgr");
  size_t evtMax = 100;
  if (appMgr.isValid()) {
    evtMax = std::stoi(appMgr->getProperty("EvtMax").toString());
  }

  // Define and register histogram vectors
  for (size_t i = 0; i < m_numLayers; i++) {
    // Energy in Layer vs. phi
    {
      std::string histName = "energyInLayer_" + std::to_string(i) + "_phi";
      m_hEnergyInLayerPhi.push_back(new TH1F(histName.c_str(),
                                             ("Energy deposited in layer " + std::to_string(i) + "vs. #phi").c_str(),
                                             1000, -m_phi, m_phi));
      m_hEnergyInLayerPhi.back()->Sumw2();
      if (m_histSvc->regHist("/det/" + histName, m_hEnergyInLayerPhi.back()).isFailure()) {
        error() << "Couldn't register histogram \"" << histName << "\"!"<< endmsg;
        return StatusCode::FAILURE;
      }
    }

    // Upstream energy vs. energy in every layer
    {
      size_t nBins = std::ceil(std::sqrt(evtMax));
      std::string histName = "upstreamEnergy_energyInLayer_" + std::to_string(i);
      std::string histTitle = "Upstream energy vs energy deposited in layer id: " + std::to_string(i);
      histTitle += ";Energy in layer id: " + std::to_string(i) + " [GeV];Upstream Energy [GeV]";
      m_hUpstreamEnergyEnergyInLayer.emplace_back(new TH2F(histName.c_str(), histTitle.c_str(),
                                                           nBins, 0., 0., nBins, 0., 0.));
      m_hUpstreamEnergyEnergyInLayer.back()->Sumw2();
      if (m_histSvc->regHist("/det/" + histName, m_hUpstreamEnergyEnergyInLayer.back()).isFailure()) {
        error() << "Couldn't register histogram \"" << histName << "\"!" << endmsg;
        return StatusCode::FAILURE;
      }
    }

    // Downstream energy vs. energy in every layer
    {
      size_t nBins = std::ceil(std::sqrt(evtMax));
      std::string histName = "downstreamEnergy_energyInLayer_" + std::to_string(i);
      std::string histTitle = "Downstream energy vs energy deposited in layer id: " + std::to_string(i);
      histTitle += ";Energy in layer id: " + std::to_string(i) + " [GeV];Downstream Energy [GeV]";
      m_hDownstreamEnergyEnergyInLayer.emplace_back(new TH2F(histName.c_str(), histTitle.c_str(),
                                                             nBins, 0., 0., nBins, 0., 0.));
      m_hDownstreamEnergyEnergyInLayer.back()->Sumw2();
      if (m_histSvc->regHist("/det/" + histName, m_hDownstreamEnergyEnergyInLayer.back()).isFailure()) {
        error() << "Couldn't register histogram \"" << histName << "\"!" << endmsg;
        return StatusCode::FAILURE;
      }
    }
  }

  // Define other histograms

  // Energies in cryostat parts and calorimeter
  m_hEnergyInLayers = bookHisto1D("energyInLayers", "Energy deposited in layer ", m_numLayers, 0, m_numLayers);
  m_hEnergyInCalo = bookHisto1D("energyInCalo", "Sum of energy deposited in all layers;E [GeV];N_{evt}", 200, 0., 0.);
  m_hEnergyInCryo = bookHisto1D("energyInCryo", "Energy deposited in cryostat;E [GeV];N_{evt}", 200, 0., 0.);
  m_hEnergyInCryoFront = bookHisto1D("energyInCryoFront", "Energy deposited in cryostat front;E [GeV];N_{evt}",
                                     200, 0., 0.);
  m_hEnergyInCryoBack = bookHisto1D("energyInCryoBack", "Energy deposited in cryostat back;E [GeV];N_{evt}",
                                    200, 0., 0.);
  m_hEnergyInCryoSides = bookHisto1D("energyInCryoSides", "Energy deposited in cryostat sides;E [GeV];N_{evt}",
                                     200, 0., 0.);
  m_hEnergyInCryoLArBathFront = bookHisto1D("energyInCryoLArBathFront",
                                            "Energy deposited in cryostat LAr bath front;E [GeV];N_{evt}",
                                            200, 0., 0.);
  m_hEnergyInCryoLArBathBack = bookHisto1D("energyInCryoLArBathBack",
                                           "Energy deposited in cryostat LAr bath back;E [GeV];N_{evt}",
                                           200, 0., 0.);

  // Particle position
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

  // Sums of calorimeter and cryostat energies
  m_hEnergyInCaloAndCryo = bookHisto1D("energyInCaloAndCryo",
                                       "Energy deposited in calorimeter and in cryostat;E [GeV];N_{evt}",
                                       200, 0., 0.);
  m_hEnergyInCaloAndCryoFront = bookHisto1D("energyInCaloAndCryoFront",
                                            "Energy deposited in calorimeter and in cryostat front;E [GeV];N_{evt}",
                                            200, 0., 0.);
  m_hEnergyInCaloAndCryoBack = bookHisto1D("energyInCaloAndCryoBack",
                                           "Energy deposited in calorimeter and in cryostat back;E [GeV];N_{evt}",
                                           200, 0., 0.);
  m_hEnergyInCryoFrontAndLArBathFront = bookHisto1D(
      "energyInCryoFrontAndLArBathFront",
      "Energy deposited in cryostat front and LAr bath front;E [GeV];N_{evt}",
      200, 0., 0.);
  m_hEnergyInCryoBackAndLArBathBack = bookHisto1D("energyInCryoBackAndLArBathBack",
                                                  "Energy deposited in cryostat back and LAr bath back;E [GeV];N_{evt}",
                                                  200, 0., 0.);

  // Correction histograms
  m_hEnergyFromUpCorr = bookHisto1D("energyFromUpCorr", "Energy from upstream correction;E [GeV];N_{evt}", 200, 0., 0.);
  m_hEnergyFromDownCorr = bookHisto1D("energyFromDownCorr", "Energy from downstream correction;E [GeV];N_{evt}",
                                      200, 0., 0.);
  m_hEnergyFromUpDownCorr = bookHisto1D("energyFromUpDownCorr", "Energy from up/downstream correction;E [GeV];N_{evt}",
                                        200, 0., 0.);

  m_hEnergyInCaloAndUpCorr = bookHisto1D("energyInCaloAndUpCorr",
                                         "Energy deposited in calorimeter and upstream correction;E [GeV];N_{evt}",
                                         200, 0., 0.);
  m_hEnergyInCaloAndDownCorr = bookHisto1D("energyInCaloAndDownCorr",
                                           "Energy deposited in calorimeter and downstream correction;E [GeV];N_{evt}",
                                           200, 0., 0.);
  m_hEnergyInCaloAndUpDownCorr = bookHisto1D(
      "energyInCaloAndUpDownCorr",
      "Energy deposited in calorimeter and up/downstream correction;E [GeV];N_{evt}",
      200, 0., 0.);

  return StatusCode::SUCCESS;
}


StatusCode UpstreamDownstreamMaterial::execute() {
  auto decoder = m_geoSvc->lcdd()->readout(m_readoutName).idSpec().decoder();

  // Check MC phi angle
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

  // Calibrate energy in the calorimeter layers
  for (size_t i = 0; i < m_numLayers; ++i) {
    sumEinLayer[i] /= m_samplingFraction[i];
  }

  // Sum energy deposited in all calorimeter layers
  float sumEinCalo = 0.;
  for (size_t i = 0; i < m_numLayers; ++i) {
    sumEinCalo += sumEinLayer[i];
  }

  // Fill histograms
  for (size_t i = 0; i < m_numLayers; ++i) {
    m_hEnergyInLayerPhi.at(i)->Fill(phi, sumEinLayer[i]);
    m_hUpstreamEnergyEnergyInLayer.at(i)->Fill(sumEinLayer[i], sumEinCryoFront + sumEinCryoLArBathFront);
    m_hDownstreamEnergyEnergyInLayer.at(i)->Fill(sumEinLayer[i], sumEinCryoBack + sumEinCryoLArBathBack);
    m_hEnergyInLayers->Fill(i, sumEinLayer[i]);
  }
  m_hEnergyInCalo->Fill(sumEinCalo);
  m_hEnergyInCryo->Fill(sumEinCryo);
  m_hEnergyInCryoFront->Fill(sumEinCryoFront);
  m_hEnergyInCryoBack->Fill(sumEinCryoBack);
  m_hEnergyInCryoSides->Fill(sumEinCryoSides);
  m_hEnergyInCryoLArBathFront->Fill(sumEinCryoLArBathFront);
  m_hEnergyInCryoLArBathBack->Fill(sumEinCryoLArBathBack);

  m_hEnergyInCaloAndCryo->Fill(sumEinCalo + sumEinCryo);
  m_hEnergyInCaloAndCryoFront->Fill(sumEinCalo + sumEinCryoFront + sumEinCryoLArBathFront);
  m_hEnergyInCaloAndCryoBack->Fill(sumEinCalo + sumEinCryoBack + sumEinCryoLArBathBack);
  m_hEnergyInCryoFrontAndLArBathFront->Fill(sumEinCryoFront + sumEinCryoLArBathFront);
  m_hEnergyInCryoBackAndLArBathBack->Fill(sumEinCryoBack + sumEinCryoLArBathBack);

  double downstreamEnergy = getDownstreamCorr(sumEinLayer[m_numLayers - 1], sumEinCalo);
  double upstreamEnergy = getUpstreamCorr(sumEinLayer[0], sumEinCalo + downstreamEnergy);
  m_hEnergyFromUpCorr->Fill(upstreamEnergy);
  m_hEnergyFromDownCorr->Fill(downstreamEnergy);
  m_hEnergyFromUpDownCorr->Fill(upstreamEnergy + downstreamEnergy);
  m_hEnergyInCaloAndUpCorr->Fill(sumEinCalo + upstreamEnergy);
  m_hEnergyInCaloAndDownCorr->Fill(sumEinCalo + downstreamEnergy);
  m_hEnergyInCaloAndUpDownCorr->Fill(sumEinCalo + upstreamEnergy + downstreamEnergy);

  // Printouts
  verbose() << "********************************************************************" << endmsg;
  verbose() << "Event num: " << m_hEnergyInCalo->GetEntries() << endmsg;
  verbose() << "Energy in layers:" << endmsg;
  for (size_t i = 0; i < m_numLayers; ++i) {
    verbose() << "  * layer " << i << ": " << sumEinLayer[i] << " GeV" << endmsg;
  }

  verbose() << "Energy in calorimeter: " << sumEinCalo << " GeV" << endmsg;
  verbose() << "Energy in cryostat: " << sumEinCryo << " GeV" << endmsg;
  verbose() << "Energy in cryostat front: " << sumEinCryoFront << " GeV" << endmsg;
  verbose() << "Energy in cryostat back: " << sumEinCryoBack << " GeV" << endmsg;
  verbose() << "Energy in cryostat sides: " << sumEinCryoSides << " GeV" << endmsg;
  verbose() << "Energy in cryostat LAr bath front: " << sumEinCryoLArBathFront << " GeV" << endmsg;
  verbose() << "Energy in cryostat LAr bath back: " << sumEinCryoLArBathBack << " GeV" << endmsg << endmsg;

  verbose() << "Energy in calorimeter and in the cryostat front: "
            << sumEinCalo + sumEinCryoFront + sumEinCryoLArBathFront << " GeV" << endmsg;
  verbose() << "Energy in calorimeter and in the cryostat back: "
            << sumEinCalo + sumEinCryoBack + sumEinCryoLArBathBack << " GeV" << endmsg << endmsg;
  verbose() << "Energy in calorimeter and in the cryostat: " << sumEinCalo + sumEinCryo << " GeV" << endmsg;

  verbose() << "Upstream correction: " << upstreamEnergy << " GeV" << endmsg;
  verbose() << "Downstream correction: " << downstreamEnergy << " GeV" << endmsg;
  verbose() << "Upstream and downstream correction: " << upstreamEnergy + downstreamEnergy << " GeV" << endmsg;
  verbose() << "Energy in calorimeter and upstream corr.: " << sumEinCalo + upstreamEnergy << " GeV" << endmsg;
  verbose() << "Energy in calorimeter and downstream corr.: " << sumEinCalo + downstreamEnergy << " GeV" << endmsg;
  verbose() << "Energy in calorimeter and up/downstream corr.: "
            << sumEinCalo + upstreamEnergy + downstreamEnergy << " GeV" << endmsg << endmsg;

  return StatusCode::SUCCESS;
}


StatusCode UpstreamDownstreamMaterial::finalize() {
  // Force buffer flush to calculate binning
  for (size_t i = 0; i < m_numLayers; ++i) {
    m_hUpstreamEnergyEnergyInLayer.at(i)->BufferEmpty();
    m_hDownstreamEnergyEnergyInLayer.at(i)->BufferEmpty();
  }

  m_hEnergyInCalo->BufferEmpty();
  m_hEnergyInCryo->BufferEmpty();
  m_hEnergyInCryoFront->BufferEmpty();
  m_hEnergyInCryoBack->BufferEmpty();
  m_hEnergyInCryoSides->BufferEmpty();
  m_hEnergyInCryoLArBathFront->BufferEmpty();
  m_hEnergyInCryoLArBathBack->BufferEmpty();
  m_hParticleMomentumXY->BufferEmpty();
  m_hParticleMomentumZY->BufferEmpty();

  m_hEnergyInCaloAndCryo->BufferEmpty();
  m_hEnergyInCaloAndCryoFront->BufferEmpty();
  m_hEnergyInCaloAndCryoBack->BufferEmpty();
  m_hEnergyInCryoFrontAndLArBathFront->BufferEmpty();
  m_hEnergyInCryoBackAndLArBathBack->BufferEmpty();

  m_hEnergyFromUpCorr->BufferEmpty();
  m_hEnergyFromDownCorr->BufferEmpty();
  m_hEnergyFromUpDownCorr->BufferEmpty();

  m_hEnergyInCaloAndUpCorr->BufferEmpty();
  m_hEnergyInCaloAndDownCorr->BufferEmpty();
  m_hEnergyInCaloAndUpDownCorr->BufferEmpty();

  return GaudiAlgorithm::finalize();
}


double UpstreamDownstreamMaterial::getUpstreamCorr(double energyInFirstLayer, double clusterEnergy) {
  double a = 0.03248;
  double b = 0.0001858;
  double c = 1.293;
  double d = -0.9203;

  double upsilon0 = a + b * clusterEnergy;
  double upsilon1 = c + d / std::sqrt(clusterEnergy);

  return upsilon0 + upsilon1 * energyInFirstLayer;
}


double UpstreamDownstreamMaterial::getDownstreamCorr(double energyInLastLayer, double clusterEnergy) {
  double a = 0.01324;
  double b = 0.002386;
  double c = 0.7686;
  double d = -1.385;
  double e = 0.001384;
  double f = 3.454;
  double g = 6.638;

  double delta0 = a + b * clusterEnergy;
  double delta1 = c + d / std::sqrt(clusterEnergy);
  double delta2 = e + f / (clusterEnergy - g);

  return delta0 + delta1 * energyInLastLayer + delta2 * std::pow(energyInLastLayer, 2);
}


TH1F* UpstreamDownstreamMaterial::bookHisto1D(const std::string& histName, const std::string& histTitle,
                                              size_t nBins, double xMin, double xMax) {
  TH1F* hist = new TH1F(histName.c_str(), histTitle.c_str(), nBins, xMin, xMax);
  hist->Sumw2();

  if (m_histSvc->regHist(("/det/" + histName).c_str(), hist).isFailure()) {
    throw GaudiException("Couldn't register histogram \"" + histName + "\"!",
                         "histogram booking error",
                         StatusCode::FAILURE);
  }

  return hist;
}
