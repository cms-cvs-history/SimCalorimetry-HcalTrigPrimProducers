import FWCore.ParameterSet.Config as cms

process = cms.Process("HcalTPG")

# Get message logger, geometry, and conditions
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'DESIGN_36_V10::All'

# use hardcoded HCAL conditions values
process.es_hardcode.toGet.extend(['Gains', 'Pedestals', 'PedestalWidths', 'QIEData', 'ElectronicsMap','ChannelQuality','RespCorrs','ZSThresholds','L1TriggerObjects','TimeCorrs','PFCorrs','LUTCorrs'])

# Use text file for the LutMetadata
process.es_ascii = cms.ESSource("HcalTextCalibrations", 
                                input = cms.VPSet ( 
        cms.PSet ( 
            object = cms.string ('LutMetadata'),
            file = cms.FileInPath('CondFormats/HcalObjects/data/HcalLutMetadata_v1.00_mc_Run1.txt')
            )
        )
)

# ESPrefers
process.es_prefer_hcalAscii    = cms.ESPrefer("HcalTextCalibrations"    , "es_ascii") 
process.es_prefer_hcalHardcode = cms.ESPrefer("HcalHardcodeCalibrations", "es_hardcode")

#Set SLHC modes
process.HcalTopologyIdealEP.SLHCMode = cms.untracked.bool(True)
process.es_hardcode.SLHCMode = cms.untracked.bool(True)
process.es_hardcode.H2Mode = cms.untracked.bool(False)

# Event output
process.load("Configuration.EventContent.EventContent_cff")

# Configure the source
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(5)
)

process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring( 
        'file:///afs/cern.ch/user/e/eberry/SLHC_Jake/CMSSW_3_6_3_SLHC1_patch1/src/test_500GeVPion_HF.root'
        ))

# Configure the output
process.dump = cms.OutputModule("PoolOutputModule",                              
   outputCommands = cms.untracked.vstring(
        'keep *',
        ),
   fileName = cms.untracked.string('tpg_test.root')
)

# import trigger primitives
process.load('SimCalorimetry.HcalTrigPrimProducers.hcalupgradetpdigi_cff')


process.p1 = cms.Path ( 
    process.simHcalUpgradeTriggerPrimitiveDigis 
)

process.end = cms.EndPath ( 
    process.dump 
)

process.schedule = cms.Schedule ()
process.schedule.append ( process.p1  ) 
process.schedule.append ( process.end ) 
