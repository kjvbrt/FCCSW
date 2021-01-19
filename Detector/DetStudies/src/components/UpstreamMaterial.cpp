#include "UpstreamMaterial.h"

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

DECLARE_COMPONENT(UpstreamMaterial)

UpstreamMaterial::UpstreamMaterial(const std::string& aName, ISvcLocator* aSvcLoc) : GaudiAlgorithm(aName, aSvcLoc), m_geoSvc("GeoSvc", aName) {
  declareProperty("deposits", m_deposits, "Energy deposits (input)");
  declareProperty("particle", m_particle, "Generated single-particle event (input)");
}
UpstreamMaterial::~UpstreamMaterial() {}

StatusCode UpstreamMaterial::initialize() {
  if (GaudiAlgorithm::initialize().isFailure()) return StatusCode::FAILURE;

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

  for (uint i = 0; i < m_numLayers; i++) {
    m_cellEnergyPhi.push_back(new TH1F(("upstreamEnergy_phi" + std::to_string(i)).c_str(),
                                       ("Energy deposited in layer " + std::to_string(i)).c_str(), 1000, -m_phi,
                                       m_phi));
    if (m_histSvc->regHist("/det/upstreamEnergy_phi" + std::to_string(i), m_cellEnergyPhi.back()).isFailure()) {
      error() << "Couldn't register histogram" << endmsg;
      return StatusCode::FAILURE;
    }

    // m_upstreamEnergyCellEnergy.emplace_back(
    //     new TH2F(("upstreamEnergy_presamplerEnergy" + std::to_string(i)).c_str(),
    //              ("Upstream energy vs energy deposited in layer " + std::to_string(i)).c_str(),
    //              2000, 0., (i == 0 ? 0.1 : 0.7) * m_energy, 2000, 0., std::log(m_energy)));
    // if (m_histSvc
    //         ->regHist("/det/upstreamEnergy_presamplerEnergy" + std::to_string(i), m_upstreamEnergyCellEnergy.back())
    //         .isFailure()) {
    //   error() << "Couldn't register hist" << endmsg;
    //   return StatusCode::FAILURE;
    // }

    m_gUpstreamEnergyCellEnergy.emplace_back(new TGraph());
    m_gUpstreamEnergyCellEnergy.back()->SetName(("upstreamEnergy_presamplerEnergy_graph_" + std::to_string(i)).c_str());
  }

  m_hSumEinLayers = new TH1F("hSumEinLayers",
                             "Sum of energy deposited in all layers;E [GeV];N_{evt}",
                             400, 0.7*m_energy, 1.3*m_energy);
  if (m_histSvc->regHist("/det/sumEinLayers", m_hSumEinLayers).isFailure()) {
    error() << "Couldn't register histogram: hSumEinLayers" << endmsg;
    return StatusCode::FAILURE;
  }

  m_hEinCryo = new TH1F("hEinCryo",
                        "Energy deposited in cryostat;E [GeV];N_{evt}",
                        200, 0., 0.2*m_energy);
  if (m_histSvc->regHist("/det/EinCryo", m_hEinCryo).isFailure()) {
    error() << "Couldn't register histogram: hEinCryo" << endmsg;
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}

StatusCode UpstreamMaterial::execute() {
  verbose() << "Event Number: " << m_gUpstreamEnergyCellEnergy.front()->GetN() << endmsg;

  auto decoder = m_geoSvc->lcdd()->readout(m_readoutName).idSpec().decoder();
  double sumEupstream = 0.;
  std::vector<double> sumEcells;
  sumEcells.assign(m_numLayers, 0);

  // first check MC phi angle
  const auto particle = m_particle.get();
  double phi = 0;
  for (const auto& part : *particle) {
    phi = atan2(part.core().p4.py, part.core().p4.px);
  }

  // get the energy deposited in the cryostat and in the detector (each layer)
  const auto deposits = m_deposits.get();
  for (const auto& hit : *deposits) {
    dd4hep::DDSegmentation::CellID cellId = hit.core().cellId;
    size_t cryoId = decoder->get(cellId, "cryo");
    if (cryoId == 0) {
      size_t layerId = decoder->get(cellId, "layer");
      sumEcells[layerId - m_firstLayerId] += hit.core().energy;
    } else {
      sumEupstream += hit.core().energy;
    }
  }

  // Calibrate the energy in the detector
  for (uint i = 0; i < m_numLayers; i++) {
    sumEcells[i] /= m_samplingFraction[i];
    m_cellEnergyPhi[i]->Fill(phi, sumEcells[i]);
    // m_upstreamEnergyCellEnergy[i]->Fill(sumEcells[i], sumEupstream);
    m_gUpstreamEnergyCellEnergy.at(i)->SetPoint(m_gUpstreamEnergyCellEnergy.at(i)->GetN(), sumEcells[i], sumEupstream);
    verbose() << "Energy deposited in layer " << i << ": " << sumEcells[i] << " GeV" << endmsg;
  }
  m_hEinCryo->Fill(sumEupstream);
  verbose() << "Energy deposited in the cryostat: " << sumEupstream << " GeV" << endmsg;

  // Sum energy deposited in all calorimeter layers
  {
    double sumEinLayers = 0.;
    for (uint i = 0; i < m_numLayers; i++) {
      sumEinLayers += sumEcells[i];
    }
    m_hSumEinLayers->Fill(sumEinLayers);
    verbose() << "Sum of energy deposited in all layers: " << sumEinLayers << " GeV" << endmsg;
  }

  return StatusCode::SUCCESS;
}


StatusCode UpstreamMaterial::finalize() {
  unsigned long int nEvt = m_gUpstreamEnergyCellEnergy.front()->GetN();
  size_t nBin = std::sqrt(nEvt) + 1;
  // size_t nBin = std::log2(nEvt) + 2;
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
    if (m_histSvc->regHist("/det/" + histName, m_upstreamEnergyCellEnergy.back()).isFailure()) {
      error() << "Couldn't register histogram:\n" << histName << endmsg;
      return StatusCode::FAILURE;
    }
  }

  for (size_t i = 0; i < m_upstreamEnergyCellEnergy.size(); ++i) {
    GraphToHist(m_gUpstreamEnergyCellEnergy.at(i), m_upstreamEnergyCellEnergy.at(i));
  }

  return GaudiAlgorithm::finalize();
}


void UpstreamMaterial::GraphToHist(TGraph* graph, TH2F* hist) {
  for (int i = 0; i < graph->GetN(); ++i) {
    hist->Fill(graph->GetPointX(i), graph->GetPointY(i));
  }
}
