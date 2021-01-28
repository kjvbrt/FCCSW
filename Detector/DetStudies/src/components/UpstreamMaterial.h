#ifndef DETSTUDIES_UPSTREAMMATERIAL_H
#define DETSTUDIES_UPSTREAMMATERIAL_H

// GAUDI
#include "GaudiAlg/GaudiAlgorithm.h"

// FCCSW
#include "FWCore/DataHandle.h"
class IGeoSvc;

// datamodel
namespace fcc {
class CaloHitCollection;
class MCParticleCollection;
}

class TH2F;
class TH1F;
class TGraph;
class ITHistSvc;

/** @class UpstreamMaterial UpstreamMaterial.h
 *
 * Histograms of energy deposited in the dead material of the calorimeter.
 * Cryostat material needs to be marked as sensitive (and therefore ID 'cryo'==1).
 * Dependence of the energy deposited in the dead material on the energy deposited in each calorimeter layer is
 * plotted.
 * Dependence of the energy deposited in the dead material on the azimuthal angle of the incoming particle (MC truth)
 * is plotted.
 *
 *  @author Anna Zaborowska
 *  @author Juraj Smiesko
 */

class UpstreamMaterial : public GaudiAlgorithm {
public:
  explicit UpstreamMaterial(const std::string&, ISvcLocator*);
  virtual ~UpstreamMaterial();
  /**  Initialize.
   *   @return status code
   */
  virtual StatusCode initialize() final;
  /**  Fills the histograms.
   *   @return status code
   */
  virtual StatusCode execute() final;
  /**  Finalize.
   *   @return status code
   */
  virtual StatusCode finalize() final;

private:
  /**
   *  Convert graph to hist
   */
  void GraphToHist(TGraph* graph, TH2F* hist);
  // Energy range in the histogram axis
  Gaudi::Property<double> m_energy{this, "energyAxis", 100, "Max energy for the axis of plot"};
  // Phi in the histogram axis
  Gaudi::Property<double> m_phi{this, "phiAxis", M_PI, "Max azimuthal angle for the axis of plot"};
  /// Handle for the energy deposits
  DataHandle<fcc::CaloHitCollection> m_deposits{"det/caloDeposits", Gaudi::DataHandle::Reader, this};
  /// Handle for the particle
  DataHandle<fcc::MCParticleCollection> m_particle{"det/particles", Gaudi::DataHandle::Reader, this};
  /// Pointer to the interface of histogram service
  SmartIF<ITHistSvc> m_histSvc;
  /// Pointer to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc;
  /// Pointer to the graphs showing upstream energy vs energy deposited in every layer
  std::vector<TGraph*> m_gUpstreamEnergyCellEnergy;
  /// Pointer to the histograms showing upstream energy vs energy deposited in every layer
  std::vector<TH2F*> m_upstreamEnergyCellEnergy;
  /// Pointer to histograms showing the energy in Phi
  std::vector<TH1F*> m_cellEnergyPhi;
  /// Pointer to histogram with energy deposited in calorimeter layers
  TH1F* m_hEnergyInLayers;
  /// Pointer to histogram with sum of energy deposited in all calorimeter layers
  TH1F* m_hSumEinLayers;
  /// Pointer to histogram with energy deposited in cryostat
  TH1F* m_hEnergyInCryo;
  /// Pointer to histogram with energy deposited in cryostat front
  TH1F* m_hEnergyInCryoFront;
  /// Pointer to histogram with energy deposited in cryostat back
  TH1F* m_hEnergyInCryoBack;
  /// Pointer to histogram with energy deposited in cryostat sides
  TH1F* m_hEnergyInCryoSides;
  /// Pointer to histogram with energy deposited in cryostat LAr bath front
  TH1F* m_hEnergyInCryoLArBathFront;
  /// Pointer to histogram with energy deposited in cryostat LAr bath back
  TH1F* m_hEnergyInCryoLArBathBack;
  /// Pointer to histogram showing particle momentum in XY plane
  TH2F* m_hParticleMomentumXY;
  /// Pointer to histogram showing particle momentum in ZY plane
  TH2F* m_hParticleMomentumZY;
  /// Pointer to histogram showing hit position in XY plane (r is assumed to be 1)
  TH2F* m_hHitPositionXY;
  /// Pointer to histogram showing hit position in ZY plane (r is assumed to be 1)
  TH2F* m_hHitPositionZY;
  /// Name of the active field
  Gaudi::Property<std::string> m_activeFieldName{this, "activeFieldName", "active", "Name of active field"};
  /// Name of the cells/layer field
  Gaudi::Property<std::string> m_layerFieldName{this, "layerFieldName", "layer", "Name of layer"};
  /// Number of layers/cells cells
  Gaudi::Property<uint> m_numLayers{this, "numLayers", 8, "Number of layers"};
  /// Id of the first layer (current design starts layer ids at 0)
  Gaudi::Property<uint> m_firstLayerId{this, "firstLayerId", 0, "ID of first layer"};
  /// Values of sampling fraction
  Gaudi::Property<std::vector<double>> m_samplingFraction{
      this, "samplingFraction", {}, "Values of sampling fraction per layer"};
  /// Name of the detector readout
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "", "Name of the readout"};
};
#endif /* DETSTUDIES_UPSTREAMMATERIAL_H */
