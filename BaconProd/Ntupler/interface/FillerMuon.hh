#ifndef BACONPROD_NTUPLER_FILLERMUON_HH
#define BACONPROD_NTUPLER_FILLERMUON_HH

#include "BaconProd/Utils/interface/TriggerTools.hh"
#include "BaconAna/Utils/interface/TTrigger.hh"
#include <vector>
#include <string>
#include <bitset>

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

// forward class declarations
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
class TClonesArray;
namespace trigger {
  class TriggerEvent;
}


namespace baconhep
{
  class FillerMuon
  {
    public:
    FillerMuon(const edm::ParameterSet &iConfig,edm::ConsumesCollector && iC);
      ~FillerMuon();
      
      bool fill(std::vector<float>               *muon_pt,
		std::vector<float>               *muon_eta,
		std::vector<float>               *muon_phi,
		std::vector<float>               *muon_ptErr,
		std::vector<float>               *muon_staPt,
		std::vector<float>               *muon_staEta,
		std::vector<float>               *muon_staPhi,
		std::vector<float>               *muon_pfPt,
		std::vector<float>               *muon_pfEta,
		std::vector<float>               *muon_pfPhi,
		std::vector<int>                 *muon_q,
		std::vector<float>               *muon_trkIso,
		std::vector<float>               *muon_ecalIso,
		std::vector<float>               *muon_hcalIso,
		std::vector<float>               *muon_chHadIso,
		std::vector<float>               *muon_gammaIso,
		std::vector<float>               *muon_neuHadIso,
		std::vector<float>               *muon_puIso,
		std::vector<float>               *muon_d0,
		std::vector<float>               *muon_dz,
		std::vector<float>               *muon_sip3d,
		std::vector<float>               *muon_tkNchi2,
		std::vector<float>               *muon_muNchi2,
		std::vector<float>               *muon_trkKink,
		std::vector<float>               *muon_glbKink,
		std::vector<int>                 *muon_nValidHits,
		std::vector<unsigned int>        *muon_typeBits,
		std::vector<unsigned int>        *muon_selectorBits,
		std::vector<unsigned int>        *muon_pogIDBits,
		std::vector<unsigned int>        *muon_nTkHits,
		std::vector<unsigned int>        *muon_nPixHits,
		std::vector<unsigned int>        *muon_nTkLayers,
		std::vector<unsigned int>        *muon_nPixLayers,
		std::vector<unsigned int>        *muon_nMatchStn,
		std::vector<int>                 *trkID,
		std::vector<TriggerObjects>      *muon_hltMatchBits,
		std::vector<float>               *vf_tC,
		std::vector<float>               *vf_dOF,
		std::vector<float>               *vf_nC,
		std::vector<float>               *vf_Prob,
		std::vector<int>                 *category,
		std::vector<int>                 *vf_Valid,
		std::vector<float>               *vf_ip,
		std::vector<float>               *tri_iso,
		std::vector<int>                 *tri_isoNtrk,
		std::vector<float>               *invmass,   
                const edm::Event		 &iEvent,	   // event info
	        const edm::EventSetup		 &iSetup,	   // event setup info
	        const reco::Vertex		 &pv,	           // event primary vertex
	        const std::vector<TriggerRecord> &triggerRecords,  // list of trigger names and objects
	        const trigger::TriggerEvent	 &triggerEvent);   // event trigger objects

      
      // Muon cuts
      double fMinPt;
      
      // EDM object collection names
      std::string fMuonName;
      edm::EDGetTokenT<reco::MuonCollection> fMuonName_token;
      std::string fPFCandName;
      edm::EDGetTokenT<reco::PFCandidateCollection> fPFCandName_token;
      std::string fTrackName;
      edm::EDGetTokenT<reco::TrackCollection> fTrackName_token;

      // general tracks cuts
      bool   fSaveTracks;
      double fTrackMinPt;
  };
}
#endif
