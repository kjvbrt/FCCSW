#ifndef RECCALORIMETER_CORRCALOSLIWINCLUSTER_H
#define RECCALORIMETER_CORRCALOSLIWINCLUSTER_H

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
 *  Apply corrections to a reconstructed cluster in EMCal barrel.
 *  Corrections for eta position and upstream energy are based on the cells in cluster, therefore they do not
 * acknowledge the energy sharing between neighbouring clusters.
 *
 *  Corrections:
 *  1) Pseudorapidity position (correction for finite granularity of detector using log-weighting)
 *      Correction uses only cells from the first defined detector in m_systemId (default ECal Barrel). The number of
 *      layers used for this correction is specified in *numLayers*. Weights are defined for each layer in *etaRecalcWeights*.
 *  2) Energy correction for pileup noise, which is parametrised per cluster as  P0 * number_of_cells ^P1. Parameters
 *      are eta-dependent and read from file *noiseFileName* as histograms called as in *pileupHistoName*
 *      (+ 0/1 respectively).
 *  3) Energy correction for the upstream material. The energy upstream is calculated as (P00 + P01 * E_clu) + (P10 +
 *      P11 * sqrt(E_clu) ) * E_firstLayer. Parameters P00, P01, P10 and P11 are eta-dependent and specified in
 *      *presamplerShiftP0*, *presamplerShiftP1*, *presamplerScaleP0* and *presamplerScaleP1*, respectively.
 *      The eta values for which the parameters are extracted are defined in *etaValues*.
 *      Energy deposited in the first layer is the uncalibrated energy, therefore sampling fraction used for calibration
 *      needs to be given in *samplingFraction*.
 *
 *  Several histograms are filled in, to monitor the upstream correction, the pileup noise, and the energy prior to the
 * corrections, as well as afterwards, and the position (eta, phi) resolution.
 *
 *  @author Anna Zaborowska
 *
 */

class CorrCaloSliWinCluster : public GaudiAlgorithm {

public:
  CorrCaloSliWinCluster(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize();

  StatusCode execute();

  StatusCode finalize();

private:
  /** Load noise histograms to the memory
   *  @return Status code informing whether loading of histograms was successful
   */
  //StatusCode GetNoiseHistos();
  /** Find the appropriate noise constant from the histogram
   *  @param[in] aEta       Pseudorapidity value of the centre of the cluster
   *  @param[in] aNumCells  Number of cells in a cluster
   *  @return               Width of the Gaussian distribution of noise per
   *                        cluster
   */
  //double getNoiseConstantPerCluster(double aEta, uint numCells);
  /**
   * Get energy in the first layer of the provided cluster.
   * @param[in] cluster   Reference to the cluster.
   * @param[in] systemID  System ID of the cluster.
   * @return              Sum of energy deposited in cells in first layer of
   *                      the calorimeter associated with the cluster.
   */
  //double energyInFirstLayer(const fcc::CaloCluster& cluster,
  //                          const size_t systemID);

  /// Handle for input calorimeter clusters collection
  DataHandle<fcc::CaloClusterCollection> m_inClusters{
      "clusters", Gaudi::DataHandle::Reader, this};
  /// Handle for corrected (output) calorimeter clusters collection
  DataHandle<fcc::CaloClusterCollection> m_correctedClusters{
      "correctedClusters", Gaudi::DataHandle::Writer, this};
  /// Pointer to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc;
  /// number of layers in the systems as in m_systemId
  Gaudi::Property<uint> m_numLayers{
      this, "numLayers", 8,
      "Number of layers for which the eta position is calculated"};
  /// ID of the first layer
  Gaudi::Property<uint> m_firstLayerID{
      this, "firstLayerID", 0, "ID of first layer"};
  /// IDs of the detectors
  Gaudi::Property<std::vector<uint>> m_systemIDs{
      this, "systemIDs", {5}, "IDs of systems"};
  /// Names of the detector readout, corresponding to system IDs in m_systemId
  Gaudi::Property<std::vector<std::string>> m_readoutNames{
      this, "readoutNames", {"ECalBarrelPhiEta"},
      "Names of the detector readout, corresponding to systemId"};
  /// maximum eta value of generated particles (for ranges in histograms)
  Gaudi::Property<double> m_etaMax{this, "etaMaxAxis", 1.6805, "Max eta value"};
  /// maximum phi value of generated particles (for ranges in histograms)
  Gaudi::Property<double> m_phiMax{this, "phiMaxAxis", M_PI, "Max phi value"};
  /// segmentation of detetor in eta (for number of bins in histograms)
  Gaudi::Property<double> m_dEta{this, "dEta", 0.01, "Segmentation in eta"};
  /// segmentation of detetor in phi (for number of bins in histograms)
  Gaudi::Property<double> m_dPhi{this, "dPhi", 2*M_PI/704, "Segmentation in phi"};
};

#endif /* RECCALORIMETER_CORRCALOSLIWINCLUSTER_H */
