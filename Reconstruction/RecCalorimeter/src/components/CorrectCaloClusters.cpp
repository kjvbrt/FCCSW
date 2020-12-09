#include "CorrectCaloClusters.h"

// FCCSW
#include "DetCommon/DetUtils.h"
#include "DetInterface/IGeoSvc.h"
#include "DetSegmentation/FCCSWGridPhiEta.h"
#include "GaudiKernel/ITHistSvc.h"

// DD4hep
#include "DD4hep/Detector.h"
#include "DDSegmentation/MultiSegmentation.h"

// our EDM
#include "datamodel/CaloCluster.h"
#include "datamodel/CaloClusterCollection.h"
#include "datamodel/CaloHitCollection.h"
#include "datamodel/MCParticleCollection.h"
#include "datamodel/GenVertexCollection.h"

// Root
#include "TFile.h"
#include "TLorentzVector.h"
#include "TFitResult.h"
#include "TGraphErrors.h"

DECLARE_COMPONENT(CorrectCaloClusters)

CorrectCaloClusters::CorrectCaloClusters(const std::string& name,
                                         ISvcLocator* svcLoc)
    : GaudiAlgorithm(name, svcLoc),
      m_geoSvc("GeoSvc", "CorrectCaloClusters") {
  declareProperty("inClusters", m_inClusters,
                  "Input cluster collection");
  declareProperty("outClusters", m_outClusters,
                  "Corrected (output) cluster collection");
}

StatusCode CorrectCaloClusters::initialize() {
  StatusCode sc = GaudiAlgorithm::initialize();
  if (sc.isFailure()) {
    return sc;
  }

  // Check if readouts exist
  {
    bool readoutMissing = false;
    for (size_t i = 0; i < m_readoutNames.size(); ++i) {
      auto readouts = m_geoSvc->lcdd()->readouts();
      if (readouts.find(m_readoutNames.value().at(i)) == readouts.end()) {
        readoutMissing = true;
      }
    }
    if (readoutMissing) {
      error() << "Missing readout, exiting!" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  // Check if readout related parameters have the same size
  if (m_systemIDs.size() != m_readoutNames.size()) {
    error() << "System ID array size and readout names array size does not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_systemIDs.size() != m_numLayers.size()) {
    error() << "System ID array size and number of layers array size does not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_systemIDs.size() != m_firstLayerIDs.size()) {
    error() << "System ID array size and first layer ID array size does not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_systemIDs.size() != m_samplingFractions.size()) {
    error() << "System ID array size and sampling fractions array size does not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_systemIDs.size() != m_P00.size()) {
    error() <<  "System ID array size and P00 parameter array size does not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_systemIDs.size() != m_P01.size()) {
    error() <<  "System ID array size and P01 parameter array size does not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_systemIDs.size() != m_P10.size()) {
    error() <<  "System ID array size and P10 parameter array size does not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_systemIDs.size() != m_P11.size()) {
    error() <<  "System ID array size and P11 parameter array size does not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }
  return StatusCode::SUCCESS;
}


StatusCode CorrectCaloClusters::execute() {
  // Get the input collection with clusters
  const fcc::CaloClusterCollection* inClusters = m_inClusters.get();

  // Initialize output clusters
  fcc::CaloClusterCollection* outClusters = initializeOutputClusters(inClusters);
  if (!outClusters) {
    error() << "Something went wrong in initialization of the output cluster collection, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (inClusters->size() != outClusters->size()) {
    error() << "Sizes of input and output cluster collections does not match, exiting!" << endmsg;
    return StatusCode::FAILURE;
  }

  // Apply upstream correction
  {
    StatusCode sc = applyUpstreamCorr(inClusters, outClusters);
    if (sc.isFailure()) {
      return sc;
    }
  }


  return StatusCode::SUCCESS;
}


StatusCode CorrectCaloClusters::finalize() {

  return GaudiAlgorithm::finalize();
}


fcc::CaloClusterCollection* CorrectCaloClusters::initializeOutputClusters(
    const fcc::CaloClusterCollection* inClusters) {

  fcc::CaloClusterCollection* outClusters = m_outClusters.createAndPut();

  for (auto const& inCluster: *inClusters) {
    fcc::CaloCluster outCluster = outClusters->create();

    outCluster.core().position.x = inCluster.core().position.x;
    outCluster.core().position.y = inCluster.core().position.y;
    outCluster.core().position.z = inCluster.core().position.z;
    verbose() << "Cluster position:" << endmsg;
    verbose() << "    x: " << outCluster.core().position.x << endmsg;
    verbose() << "    y: " << outCluster.core().position.y << endmsg;
    verbose() << "    z: " << outCluster.core().position.z << endmsg;
    outCluster.core().energy = inCluster.core().energy;
    verbose() << "Cluster energy: " << outCluster.core().energy << endmsg;
  }

  return outClusters;
}


StatusCode CorrectCaloClusters::applyUpstreamCorr(const fcc::CaloClusterCollection* inClusters,
                                                  fcc::CaloClusterCollection* outClusters) {
  for (size_t i = 0; i < m_readoutNames.size(); ++i) {
    for (size_t j = 0; j < inClusters->size(); ++j) {
      double energyInFirstLayer = getEnergyInFirstLayer(inClusters->at(j),
                                                        m_readoutNames[i],
                                                        m_systemIDs[i],
                                                        m_firstLayerIDs[i]);
      verbose() << "Energy in first layer: " << energyInFirstLayer << endmsg;

      if (energyInFirstLayer == 0.) {
        continue;
      }

      double P0 = m_P00[i] + m_P01[i] * inClusters->at(i).energy();
      double P1 = m_P10[i] + m_P11[i] / std::sqrt(inClusters->at(i).energy());
      double energyCorr = P0 + P1 * energyInFirstLayer * m_samplingFractions[i][m_firstLayerIDs[i]];
      verbose() << "Energy correction: " << energyCorr << endmsg;
      outClusters->at(i).core().energy += energyCorr;
      verbose() << "Corrected cluster energy: " << outClusters->at(i).core().energy << endmsg;
    }
  }

  return StatusCode::SUCCESS;
}


double CorrectCaloClusters::getEnergyInFirstLayer(const fcc::CaloCluster& cluster,
                                                  const std::string& readoutName,
                                                  size_t systemID,
                                                  size_t firstLayerID) {
  dd4hep::DDSegmentation::BitFieldCoder* decoder = m_geoSvc->lcdd()->readout(readoutName).idSpec().decoder();

  double energy = 0;
  for (auto cell = cluster.hits_begin(); cell != cluster.hits_end(); ++cell) {
    dd4hep::DDSegmentation::CellID cellID = cell->core().cellId;
    if (decoder->get(cellID, "system") != systemID) {
      continue;
    }
    if (decoder->get(cellID, "layer") != firstLayerID) {
      continue;
    }

    energy += cell->core().energy;
  }

  return energy;
}
