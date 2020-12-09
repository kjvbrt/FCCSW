#ifndef RECCALORIMETER_CORRECTCALOCLUSTERS_H
#define RECCALORIMETER_CORRECTCALOCLUSTERS_H

// FCCSW
#include "FWCore/DataHandle.h"
#include "GaudiKernel/RndmGenerators.h"
class IGeoSvc;
class IRndmGenSvc;
class ITHistSvc;

// Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiKernel/MsgStream.h"

// EDM
#include "datamodel/CaloCluster.h"

namespace fcc {
  class CaloClusterCollection;
  class CaloHitCollection;
  class MCParticleCollection;
  class GenVertexCollection;
}

namespace dd4hep {
  namespace DDSegmentation {
    class FCCSWGridPhiEta;
    class MultiSegmentation;
    class BitFieldCoder;
  }
}

/** @class CorrCaloSliWinCluster
 *
 *  Apply corrections to clusters reconstructed in ECAL.
 *
 *  Corrections:
 *  1) Energy correction for the upstream material. The energy upstream is calculated as (P00 + P01 * E_clu) + (P10 +
 *     P11 * sqrt(E_clu)) * E_firstLayer. Parameters P00, P01, P10 and P11 can be eta-dependent and they are specified
 *     in *presamplerShiftP0*, *presamplerShiftP1*, *presamplerScaleP0* and *presamplerScaleP1*, respectively.
 *      The eta values for which the parameters are extracted are defined in *etaValues*.
 *      Energy deposited in the first layer is the uncalibrated energy, therefore sampling fraction used for calibration
 *      needs to be given in *samplingFraction*.
 *
 *  Several histograms are filled in, to monitor the upstream correction, the pileup noise, and the energy prior to the
 * corrections, as well as afterwards, and the position (eta, phi) resolution.
 *
 *  Based on similar correction by Anna Zaborowska.
 *
 *  @author Juraj Smiesko
 *
 */

class CorrectCaloClusters : public GaudiAlgorithm {

public:
  CorrectCaloClusters(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize();

  StatusCode execute();

  StatusCode finalize();

private:
  /**
   * Initialize output calorimeter cluster collection.
   *
   * @param[in] inClusters  Pointer to the input cluster collection.
   *
   * @return                Pointer to the output cluster collection.
   */
  fcc::CaloClusterCollection* initializeOutputClusters(const fcc::CaloClusterCollection* inClusters);

  /**
   * Apply upstream correction to the output clusters.
   *
   * @param[in]  inClusters   Pointer to the input cluster collection.
   * @param[out] outClusters  Pointer to the output cluster collection.
   *
   * @return                  Status code.
   */
  StatusCode applyUpstreamCorr(const fcc::CaloClusterCollection* inClusters, fcc::CaloClusterCollection* outClusters);

  /**
   * Get sum of energy from cells in the first layer.
   * This energy is not calibrated.
   *
   * @param[in]  cluster       Pointer to cluster of interest.
   * @param[in]  readoutName   Name of the readout.
   * @param[in]  firstLayerID  ID of the first layer of the readout.
   */
  double getEnergyInFirstLayer(const fcc::CaloCluster& cluster,
                               const std::string& readoutName,
                               size_t systemID,
                               size_t firstLayerID);

  /// Handle for input calorimeter clusters collection
  DataHandle<fcc::CaloClusterCollection> m_inClusters {
    "inClusters", Gaudi::DataHandle::Reader, this
  };
  /// Handle for corrected (output) calorimeter clusters collection
  DataHandle<fcc::CaloClusterCollection> m_outClusters {
      "outClusters", Gaudi::DataHandle::Writer, this
  };

  /// Pointer to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc;

  /// IDs of the detectors
  Gaudi::Property<std::vector<size_t>> m_systemIDs {
      this, "systemIDs", {4}, "IDs of systems"
  };
  /// Names of the detector readouts, corresponding to system IDs
  Gaudi::Property<std::vector<std::string>> m_readoutNames {
      this, "readoutNames", {"ECalBarrelPhiEta"},
      "Names of the detector readout, corresponding to systemId"
  };
  /// Numbers of layers of the systems
  Gaudi::Property<std::vector<size_t>> m_numLayers {
      this, "numLayers", {8}, "Numbers of layers of the systems"};
  /// IDs of the first layer in the systems
  Gaudi::Property<std::vector<size_t>> m_firstLayerIDs {
      this, "firstLayerIDs", {0}, "IDs of first layer in the systems"
  };
  /// Values of sampling fractions used for energy calibration of the systems
  Gaudi::Property<std::vector<std::vector<double>>> m_samplingFractions {
      this, "samplingFractions",
      {{0.299041341789, 0.1306220735, 0.163243999965, 0.186360269398,
        0.203778124831, 0.216211280314, 0.227140796653, 0.243315422934}},
      "Values of sampling fractions used in energy calibration of the systems"};

  /// Upstream correction parameter P00 for upstream correction
  Gaudi::Property<std::vector<double>> m_P00 {
      this, "P00", {0.1037}, "Upstream material parameter P00"};
  /// Upstream correction parameter P01 for upstream correction
  Gaudi::Property<std::vector<double>> m_P01 {
      this, "P01", {0.0007507}, "Upstream material parameter P01"};
  /// Upstream correction parameter P10 for upstream correction
  Gaudi::Property<std::vector<double>> m_P10 {
      this, "P10", {0.1382}, "Upstream material parameter P10"};
  /// Upstream correction parameter P11 for upstream correction
  Gaudi::Property<std::vector<double>> m_P11 {
      this, "P11", {1.002}, "Upstream material parameter P11"};
};

#endif /* RECCALORIMETER_CORRECTCALOCLUSTERS_H */
