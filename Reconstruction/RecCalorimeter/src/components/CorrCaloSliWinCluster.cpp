#include "CorrCaloSliWinCluster.h"

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

DECLARE_COMPONENT(CorrCaloSliWinCluster)

CorrCaloSliWinCluster::CorrCaloSliWinCluster(const std::string& name,
                                             ISvcLocator* svcLoc)
    : GaudiAlgorithm(name, svcLoc),
      m_geoSvc("GeoSvc", "CorrCaloSliWinCluster") {
  declareProperty("clusters", m_inClusters, "Input cluster collection");
  declareProperty("correctedClusters", m_correctedClusters,
                  "Corrected (output) cluster collection");
}

StatusCode CorrCaloSliWinCluster::initialize() {
  StatusCode sc = GaudiAlgorithm::initialize();
  if (sc.isFailure()) {
    return sc;
  }

  for (size_t i = 0; i < m_readoutNames.size(); ++i) {
    // check if readouts exist
    auto readouts = m_geoSvc->lcdd()->readouts();
    if (readouts.find(m_readoutNames.value().at(i)) == readouts.end()) {
      error() << "Readout <<" << m_readoutNames.value().at(i) << ">> not found!"
              << endmsg;
      return StatusCode::FAILURE;
    }
  }

  return StatusCode::SUCCESS;
}


StatusCode CorrCaloSliWinCluster::execute() {
  info() << "Execute start" <<  endmsg;
  // Get the input collection with clusters
  const fcc::CaloClusterCollection* inClusters = m_inClusters.get();
  fcc::CaloClusterCollection* correctedClusters = m_correctedClusters.createAndPut();
  info() << "Clusters in Collection: " << inClusters->size() << endmsg;


  info() << "Execute end" <<  endmsg;

  return StatusCode::SUCCESS;
}


StatusCode CorrCaloSliWinCluster::finalize() {

  return GaudiAlgorithm::finalize();
}
