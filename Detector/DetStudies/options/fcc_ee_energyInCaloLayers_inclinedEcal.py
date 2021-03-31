from Gaudi.Configuration import *

from GaudiKernel.SystemOfUnits import MeV, GeV

# Electron momentum in GeV
momentum = 50
# Theta and its spread in degrees
theta = 90.
thetaSpread = 10.

samplingFractions = [0.303451138049] * 1 + [0.111872504159] * 1 + [0.135806495306] * 1 + [0.151772636618] * 1 + \
                    [0.163397436122] * 1 + [0.172566977313] * 1 + [0.179855253903] * 1 + [0.186838417657] * 1 + \
                    [0.192865946689] * 1 + [0.197420241611] * 1 + [0.202066552306] * 1 + [0.22646764465] * 1

#---------------------------------------------------------------------------------------------------------------------

# Random string for output files
import uuid
rndstr = uuid.uuid4().hex[0:16]

# Data service
from Configurables import FCCDataSvc
podioevent = FCCDataSvc("EventDataSvc")

# Particle gun setup
_pi = 3.14159

from Configurables import  MomentumRangeParticleGun
pgun = MomentumRangeParticleGun("ParticleGun_Electron")
pgun.PdgCodes = [11]
pgun.MomentumMin = momentum * GeV
pgun.MomentumMax = momentum * GeV
pgun.PhiMin = 0
pgun.PhiMax = 2 * _pi
# theta = 90 degrees (eta = 0)
pgun.ThetaMin = (theta - thetaSpread/2) * _pi / 180.
pgun.ThetaMax = (theta + thetaSpread/2) * _pi / 180.

from Configurables import GenAlg
genalg_pgun = GenAlg()
genalg_pgun.SignalProvider = pgun
genalg_pgun.hepmc.Path = "hepmc"

from Configurables import HepMCToEDMConverter
hepmc_converter = HepMCToEDMConverter()
hepmc_converter.hepmc.Path="hepmc"
hepmc_converter.genparticles.Path="GenParticles"
hepmc_converter.genvertices.Path="GenVertices"

# DD4hep geometry service
from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc",
                    detectors=['file:Detector/DetFCCeeIDEA-LAr/compact/FCCee_DectEmptyMaster.xml',
                               'file:Detector/DetFCCeeECalInclined/compact/FCCee_ECalBarrel_upstream.xml'],
                    OutputLevel = INFO)

# Geant4 service
# Configures the Geant simulation: geometry, physics list and user actions
from Configurables import SimG4Svc
geantservice = SimG4Svc("SimG4Svc",
                        detector='SimG4DD4hepDetector',
                        physicslist="SimG4FtfpBert",
                        actions="SimG4FullSimActions")
geantservice.g4PostInitCommands += ["/run/setCut 0.1 mm"]

# Geant4 algorithm
# Translates EDM to G4Event, passes the event to G4, writes out outputs via tools
# and a tool that saves the calorimeter hits
from Configurables import SimG4Alg, SimG4SaveCalHits
saveecaltool = SimG4SaveCalHits("saveECalBarrelHits",
                                readoutNames = ["ECalBarrelEta"])
saveecaltool.positionedCaloHits.Path = "ECalBarrelPositionedHits"
saveecaltool.caloHits.Path = "ECalBarrelHits"

from Configurables import SimG4PrimariesFromEdmTool
particle_converter = SimG4PrimariesFromEdmTool("EdmConverter")
particle_converter.genParticles.Path = "GenParticles"

# next, create the G4 algorithm, giving the list of names of tools ("XX/YY")
geantsim = SimG4Alg("SimG4Alg",
                    outputs= ["SimG4SaveCalHits/saveECalBarrelHits"],
                    eventProvider = particle_converter,
                    OutputLevel = INFO)

from Configurables import CreateCaloCells
createcellsBarrel = CreateCaloCells("CreateCaloCellsBarrel",
                                    doCellCalibration=False,
                                    addCellNoise=False, filterCellNoise=False)
createcellsBarrel.hits.Path="ECalBarrelHits"
createcellsBarrel.cells.Path="ECalBarrelCells"

from Configurables import EnergyInCaloLayers
caloLayers = EnergyInCaloLayers("caloLayers",
                                readoutName="ECalBarrelEta",
                                numLayers=12,
                                # sampling fraction is given as the upstream correction will be applied on calibrated
                                # cells
                                samplingFraction=samplingFractions,
                                OutputLevel=VERBOSE)
caloLayers.deposits.Path = "ECalBarrelCells"
caloLayers.particle.Path = "GenParticles"

# Print CPU information
from Configurables import AuditorSvc, ChronoAuditor
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
geantsim.AuditExecute = True
caloLayers.AuditExecute = True

from Configurables import PodioOutput
### PODIO algorithm
out = PodioOutput("out", OutputLevel=INFO)
out.outputCommands = ["drop *", "keep energyInLayer", "keep energyInCryo", "keep particleVec"]
out.filename = "fccee_upstreamDownstreamMaterial_inclinedEcal_%ideg_%iGeV_%s.root" % (theta, momentum, rndstr)

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr(TopAlg=[genalg_pgun, hepmc_converter, geantsim, createcellsBarrel, caloLayers, out],
               EvtSel='NONE',
               EvtMax=10,
               # order is important, as GeoSvc is needed by G4SimSvc
               ExtSvc=[podioevent, geoservice, geantservice, audsvc],
               OutputLevel=INFO
)
