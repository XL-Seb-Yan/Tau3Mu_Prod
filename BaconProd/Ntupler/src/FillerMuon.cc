#include "BaconProd/Ntupler/interface/FillerMuon.hh"
#include "BaconProd/Utils/interface/TriggerTools.hh"
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/Utils/interface/TTrigger.hh"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <bitset>

using namespace baconhep;

//--------------------------------------------------------------------------------------------------
FillerMuon::FillerMuon(const edm::ParameterSet &iConfig,edm::ConsumesCollector && iC):
  fMinPt     (iConfig.getUntrackedParameter<double>("minPt",0)),
  fMuonName  (iConfig.getUntrackedParameter<std::string>("edmName","muons")),
  fPFCandName(iConfig.getUntrackedParameter<std::string>("edmPFCandName","particleFlow")),
  fTrackName (iConfig.getUntrackedParameter<std::string>("edmTrackName","generalTracks")),
  fSaveTracks(iConfig.getUntrackedParameter<bool>("doSaveTracks",false)),
  fTrackMinPt(iConfig.getUntrackedParameter<double>("minTrackPt",20))
{
  fMuonName_token = iC.consumes<reco::MuonCollection>(fMuonName);
  fPFCandName_token = iC.consumes<reco::PFCandidateCollection>(fPFCandName);
  fTrackName_token = iC.consumes<reco::TrackCollection>(fTrackName);
}

//--------------------------------------------------------------------------------------------------
FillerMuon::~FillerMuon(){}
//--------------------------------------------------------------------------------------------------

// Here are the implementaions of triplet filling //
// Scenario A
void FillA(std::vector<reco::Muon> muonsel,
	   const reco::TrackCollection *trackCol,
	   std::vector<float> *muon_pt,
           std::vector<float> *muon_eta,
           std::vector<float> *muon_phi,
           std::vector<float> *muon_ptErr,
           std::vector<float> *muon_staPt,
           std::vector<float> *muon_staEta,
           std::vector<float> *muon_staPhi,
           std::vector<float> *muon_pfPt,
           std::vector<float> *muon_pfEta,
           std::vector<float> *muon_pfPhi,
           std::vector<int>   *muon_q,
           std::vector<float> *muon_trkIso,
           std::vector<float> *muon_ecalIso,
           std::vector<float> *muon_hcalIso,
           std::vector<float> *muon_chHadIso,
           std::vector<float> *muon_gammaIso,
           std::vector<float> *muon_neuHadIso,
           std::vector<float> *muon_puIso,
       	   std::vector<float> *muon_d0,
           std::vector<float> *muon_dz,
           std::vector<float> *muon_sip3d,
           std::vector<float> *muon_tkNchi2,
           std::vector<float> *muon_muNchi2,
           std::vector<float> *muon_trkKink,
           std::vector<float> *muon_glbKink,
           std::vector<int>   *muon_nValidHits,
           std::vector<unsigned int>   *muon_typeBits,
           std::vector<unsigned int>   *muon_selectorBits,
           std::vector<unsigned int>   *muon_pogIDBits,
           std::vector<unsigned int>   *muon_nTkHits,
           std::vector<unsigned int>   *muon_nPixHits,
           std::vector<unsigned int>   *muon_nTkLayers,
           std::vector<unsigned int>   *muon_nPixLayers,
           std::vector<unsigned int>   *muon_nMatchStn,
           std::vector<int>   *muon_trkID,
           std::vector<TriggerObjects> *muon_hltMatchBits,
           std::vector<float> *vf_tC,
           std::vector<float> *vf_dOF,
           std::vector<float> *vf_nC,
           std::vector<float> *vf_Prob,
	   std::vector<int>   *category,
	   std::vector<int>   *vf_Valid,
	   std::vector<float> *invmass,
	   const edm::Event &iEvent, const edm::EventSetup &iSetup, const reco::Vertex &pv, 
           const std::vector<TriggerRecord> &triggerRecords,
           const trigger::TriggerEvent &triggerEvent)
{
  // Load trigger menu
  const baconhep::TTrigger triggerMenu("/afs/cern.ch/user/x/xuyan/public/TriggerMenu/HLT_50nsGRun");

  reco::Muon m[3];
  edm::ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  // Loop all combinations
  for(int i=0; i<(int)muonsel.size()-2; i++){
    m[0] = muonsel[i];
    if(m[0].innerTrack().isNull()) continue; //Check reference validity
    for(int j=i+1; j<(int)muonsel.size()-1; j++){
      m[1] = muonsel[j];
      if(m[1].innerTrack().isNull()) continue;
      for(int k=j+1; k<(int)muonsel.size(); k++){
	m[2] = muonsel[k];
	if(m[2].innerTrack().isNull()) continue;

	// Triplet mass region
	TLorentzVector lv1, lv2, lv3;
	lv1.SetPtEtaPhiM(m[0].muonBestTrack()->pt(), m[0].muonBestTrack()->eta(), m[0].muonBestTrack()->phi(), 0.105658369);
	lv2.SetPtEtaPhiM(m[1].muonBestTrack()->pt(), m[1].muonBestTrack()->eta(), m[1].muonBestTrack()->phi(), 0.105658369);
	lv3.SetPtEtaPhiM(m[2].muonBestTrack()->pt(), m[2].muonBestTrack()->eta(), m[2].muonBestTrack()->phi(), 0.105658369);
	float mass_cut = (lv1+lv2+lv3).M();
	if(mass_cut > 2.4 || mass_cut < 1.4) continue;
	invmass->push_back(mass_cut);
	invmass->push_back(-99); //Take the empty space for normalization channel invmass

	// Build tracks
	std::vector<reco::TransientTrack> t_trks;
	reco::TrackRef trk1 = m[0].innerTrack();
	reco::TrackRef trk2 = m[1].innerTrack();
	reco::TrackRef trk3 = m[2].innerTrack();
	t_trks.push_back(theB->build(trk1));
	t_trks.push_back(theB->build(trk2));
	t_trks.push_back(theB->build(trk3));

	// Vertex fit
	KalmanVertexFitter kvf;
	TransientVertex fv = kvf.vertex(t_trks);
	if(fv.isValid()){
	  vf_Valid->push_back(1);
	  vf_tC->push_back(fv.totalChiSquared());
	  vf_dOF->push_back(fv.degreesOfFreedom());
	  vf_nC->push_back(fv.totalChiSquared()/fv.degreesOfFreedom());
	  vf_Prob->push_back(TMath::Prob(fv.totalChiSquared(),(int)fv.degreesOfFreedom()));
	}
	else{
	  vf_Valid->push_back(0);
	  // Vertex Fitting info
	  vf_tC->push_back(-99);
	  vf_dOF->push_back(-99);
	  vf_nC->push_back(-99);
	  vf_Prob->push_back(-99);
	}

	// Store all possible triples
	// Muon info
	for(int i=0; i<3; i++){
	  // Kinematics
	  muon_pt->push_back(m[i].muonBestTrack()->pt());
	  muon_eta->push_back(m[i].muonBestTrack()->eta());
	  muon_phi->push_back(m[i].muonBestTrack()->phi());
	  muon_ptErr->push_back(m[i].muonBestTrack()->ptError());
	  muon_q->push_back(m[i].muonBestTrack()->charge());
	  muon_staPt->push_back(m[i].standAloneMuon().isNonnull() ? m[i].standAloneMuon()->pt() : 0);
	  muon_staEta->push_back(m[i].standAloneMuon().isNonnull() ? m[i].standAloneMuon()->eta() : 0);
	  muon_staPhi->push_back(m[i].standAloneMuon().isNonnull() ? m[i].standAloneMuon()->phi() : 0);
	  muon_pfPt->push_back(m[i].pfP4().pt());
	  muon_pfEta->push_back(m[i].pfP4().eta());
	  muon_pfPhi->push_back(m[i].pfP4().phi());

	  // Isolation
	  muon_trkIso->push_back(m[i].isolationR03().sumPt);
	  muon_ecalIso->push_back(m[i].isolationR03().emEt);
	  muon_chHadIso->push_back(m[i].pfIsolationR04().sumChargedHadronPt);
	  muon_gammaIso->push_back(m[i].pfIsolationR04().sumPhotonEt);
	  muon_neuHadIso->push_back(m[i].pfIsolationR04().sumNeutralHadronEt);
	  muon_puIso->push_back(m[i].pfIsolationR04().sumPUPt);

	  // Impact Parameter
	  muon_d0->push_back((-1)*(m[i].muonBestTrack()->dxy(pv.position())));
	  muon_dz->push_back(m[i].muonBestTrack()->dz(pv.position()));

	  const reco::TransientTrack &tt = theB->build(m[i].muonBestTrack());
	  const double thesign = ((-1)*(m[i].muonBestTrack()->dxy(pv.position())) >= 0) ? 1. : -1.;
	  const std::pair<bool,Measurement1D> &ip3d = IPTools::absoluteImpactParameter3D(tt,pv);
	  muon_sip3d->push_back(ip3d.first ? thesign*ip3d.second.value() / ip3d.second.error() : -999.);

	  // Identification
	  muon_tkNchi2->push_back(m[i].innerTrack().isNonnull() ? m[i].innerTrack()->normalizedChi2(): -999.);
	  muon_muNchi2->push_back(m[i].isGlobalMuon() ? m[i].globalTrack()->normalizedChi2() : -999.);
	  muon_trkKink->push_back(m[i].combinedQuality().trkKink);
	  muon_glbKink->push_back(m[i].combinedQuality().glbKink);
	  muon_typeBits->push_back(m[i].type());

	  unsigned int selectorbits=0;
	  if(muon::isGoodMuon(m[i],muon::All))                                    selectorbits |= baconhep::kAll;
	  if(muon::isGoodMuon(m[i],muon::AllGlobalMuons))                         selectorbits |= baconhep::kAllGlobalMuons;
	  if(muon::isGoodMuon(m[i],muon::AllStandAloneMuons))                     selectorbits |= baconhep::kAllStandAloneMuons;
	  if(muon::isGoodMuon(m[i],muon::AllTrackerMuons))                        selectorbits |= baconhep::kAllTrackerMuons;
	  if(muon::isGoodMuon(m[i],muon::TrackerMuonArbitrated))                  selectorbits |= baconhep::kTrackerMuonArbitrated;
	  if(muon::isGoodMuon(m[i],muon::AllArbitrated))                          selectorbits |= baconhep::kAllArbitrated;
	  if(muon::isGoodMuon(m[i],muon::GlobalMuonPromptTight))                  selectorbits |= baconhep::kGlobalMuonPromptTight;
	  if(muon::isGoodMuon(m[i],muon::TMLastStationLoose))                     selectorbits |= baconhep::kTMLastStationLoose;
	  if(muon::isGoodMuon(m[i],muon::TMLastStationTight))                     selectorbits |= baconhep::kTMLastStationTight;
	  if(muon::isGoodMuon(m[i],muon::TM2DCompatibilityLoose))                 selectorbits |= baconhep::kTM2DCompatibilityLoose;
	  if(muon::isGoodMuon(m[i],muon::TM2DCompatibilityTight))                 selectorbits |= baconhep::kTM2DCompatibilityTight;
	  if(muon::isGoodMuon(m[i],muon::TMOneStationLoose))                      selectorbits |= baconhep::kTMOneStationLoose;
	  if(muon::isGoodMuon(m[i],muon::TMOneStationTight))                      selectorbits |= baconhep::kTMOneStationTight;
	  if(muon::isGoodMuon(m[i],muon::TMLastStationOptimizedLowPtLoose))       selectorbits |= baconhep::kTMLastStationOptimizedLowPtLoose;
	  if(muon::isGoodMuon(m[i],muon::TMLastStationOptimizedLowPtTight))       selectorbits |= baconhep::kTMLastStationOptimizedLowPtTight;
	  if(muon::isGoodMuon(m[i],muon::GMTkChiCompatibility))                   selectorbits |= baconhep::kGMTkChiCompatibility;
	  if(muon::isGoodMuon(m[i],muon::GMStaChiCompatibility))                  selectorbits |= baconhep::kGMStaChiCompatibility;
	  if(muon::isGoodMuon(m[i],muon::GMTkKinkTight))                          selectorbits |= baconhep::kGMTkKinkTight;
	  if(muon::isGoodMuon(m[i],muon::TMLastStationAngLoose))                  selectorbits |= baconhep::kTMLastStationAngLoose;
	  if(muon::isGoodMuon(m[i],muon::TMLastStationAngTight))                  selectorbits |= baconhep::kTMLastStationAngTight;
	  if(muon::isGoodMuon(m[i],muon::TMOneStationAngLoose))                   selectorbits |= baconhep::kTMOneStationAngLoose;
	  if(muon::isGoodMuon(m[i],muon::TMOneStationAngTight))                   selectorbits |= baconhep::kTMOneStationAngTight;
	  if(muon::isGoodMuon(m[i],muon::TMLastStationOptimizedBarrelLowPtLoose)) selectorbits |= baconhep::kTMLastStationOptimizedBarrelLowPtLoose;
	  if(muon::isGoodMuon(m[i],muon::TMLastStationOptimizedBarrelLowPtTight)) selectorbits |= baconhep::kTMLastStationOptimizedBarrelLowPtTight;
	  if(muon::isGoodMuon(m[i],muon::RPCMuLoose))                             selectorbits |= baconhep::kRPCMuLoose;
	  muon_selectorBits->push_back(selectorbits);

	  unsigned int pogidbits=0;
	  if(muon::isLooseMuon(m[i]))      pogidbits |= baconhep::kPOGLooseMuon;
	  if(muon::isTightMuon(m[i], pv))  pogidbits |= baconhep::kPOGTightMuon;
	  if(muon::isSoftMuon(m[i], pv))   pogidbits |= baconhep::kPOGSoftMuon;
	  if(muon::isHighPtMuon(m[i], pv)) pogidbits |= baconhep::kPOGHighPtMuon;
	  muon_pogIDBits->push_back(pogidbits);
	    
	  muon_nValidHits->push_back(m[i].isGlobalMuon() ? m[i].globalTrack()->hitPattern().numberOfValidMuonHits(): 0);
	  muon_nTkHits->push_back(m[i].innerTrack().isNonnull() ? m[i].innerTrack()->hitPattern().numberOfValidTrackerHits(): 0);
	  muon_nPixHits->push_back(m[i].innerTrack().isNonnull() ? m[i].innerTrack()->hitPattern().numberOfValidPixelHits() : 0);
	  muon_nTkLayers->push_back(m[i].innerTrack().isNonnull() ? m[i].innerTrack()->hitPattern().trackerLayersWithMeasurement() : 0);
	  muon_nPixLayers->push_back(m[i].innerTrack().isNonnull() ? m[i].innerTrack()->hitPattern().pixelLayersWithMeasurement() : 0); 
	  muon_nMatchStn->push_back(m[i].numberOfMatchedStations());

	  int trkid = -1;
	  if(m[i].innerTrack().isNonnull()) {
	    int trkIndex = -1;
	    for(reco::TrackCollection::const_iterator itTrk = trackCol->begin(); itTrk!=trackCol->end(); ++itTrk) {
	      trkIndex++;
	      if(m[i].innerTrack().get() == &(*itTrk)) {
		trkid = trkIndex; //Check if this is correct+++++++++++++++++++++++++++++++++++
		break;
	      }
	    }
	  }
	  muon_trkID->push_back(trkid);
	  muon_hltMatchBits->push_back(TriggerTools::matchHLT(m[i].eta(), m[i].phi(), triggerRecords, triggerEvent));
	}

	// Category info
	category->push_back(1);
      }
    }
  }
}
//--------------------------------------End of Scenario A------------------------------------------

// Scenario B
void FillB(std::vector<reco::Muon> muonsel,
	   std::vector<reco::Track> trksel,
	   const reco::PFCandidateCollection *pfCandCol,
	   const reco::TrackCollection *trackCol,
	   std::vector<float> *muon_pt,
           std::vector<float> *muon_eta,
           std::vector<float> *muon_phi,
           std::vector<float> *muon_ptErr,
           std::vector<float> *muon_staPt,
           std::vector<float> *muon_staEta,
           std::vector<float> *muon_staPhi,
           std::vector<float> *muon_pfPt,
           std::vector<float> *muon_pfEta,
           std::vector<float> *muon_pfPhi,
           std::vector<int>   *muon_q,
           std::vector<float> *muon_trkIso,
           std::vector<float> *muon_ecalIso,
           std::vector<float> *muon_hcalIso,
           std::vector<float> *muon_chHadIso,
           std::vector<float> *muon_gammaIso,
           std::vector<float> *muon_neuHadIso,
           std::vector<float> *muon_puIso,
       	   std::vector<float> *muon_d0,
           std::vector<float> *muon_dz,
           std::vector<float> *muon_sip3d,
           std::vector<float> *muon_tkNchi2,
           std::vector<float> *muon_muNchi2,
           std::vector<float> *muon_trkKink,
           std::vector<float> *muon_glbKink,
           std::vector<int>   *muon_nValidHits,
           std::vector<unsigned int>   *muon_typeBits,
           std::vector<unsigned int>   *muon_selectorBits,
           std::vector<unsigned int>   *muon_pogIDBits,
           std::vector<unsigned int>   *muon_nTkHits,
           std::vector<unsigned int>   *muon_nPixHits,
           std::vector<unsigned int>   *muon_nTkLayers,
           std::vector<unsigned int>   *muon_nPixLayers,
           std::vector<unsigned int>   *muon_nMatchStn,
           std::vector<int>   *muon_trkID,
           std::vector<TriggerObjects> *muon_hltMatchBits,
           std::vector<float> *vf_tC,
           std::vector<float> *vf_dOF,
           std::vector<float> *vf_nC,
           std::vector<float> *vf_Prob,
	   std::vector<int>   *category,
	   std::vector<int>   *vf_Valid,
	   std::vector<float> *invmass,
	   const edm::Event &iEvent, const edm::EventSetup &iSetup, const reco::Vertex &pv, 
           const std::vector<TriggerRecord> &triggerRecords,
           const trigger::TriggerEvent &triggerEvent)
{
  // Load trigger menu
  const baconhep::TTrigger triggerMenu("/afs/cern.ch/user/x/xuyan/public/TriggerMenu/HLT_50nsGRun");
  
  reco::Muon m[2];
  reco::Track t[1];
  edm::ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  // Loop all combinations
  for(int i=0; i<(int)muonsel.size()-1; i++){
    m[0] = muonsel[i];
    if(m[0].innerTrack().isNull()) continue; //Check reference validity
    for(int j=i+1; j<(int)muonsel.size(); j++){
      m[1] = muonsel[j];
      if(m[1].innerTrack().isNull()) continue;
      for(int k=0; k<(int)trksel.size(); k++){
	t[0] = trksel[k];

	// Triplet mass region
	TLorentzVector lv1, lv2, lv3;
	lv1.SetPtEtaPhiM(m[0].muonBestTrack()->pt(), m[0].muonBestTrack()->eta(), m[0].muonBestTrack()->phi(), 0.105658369);
	lv2.SetPtEtaPhiM(m[1].muonBestTrack()->pt(), m[1].muonBestTrack()->eta(), m[1].muonBestTrack()->phi(), 0.105658369);
	lv3.SetPtEtaPhiM(t[0].pt(), t[0].eta(), t[0].phi(), 0.105658369);
	float mass_cut_sig = (lv1+lv2+lv3).M();
	lv3.SetPtEtaPhiM(t[0].pt(), t[0].eta(), t[0].phi(), 0.13957);
	float mass_cut_nor = (lv1+lv2+lv3).M();
	bool isReject_sig = false;
	bool isReject_nor = false;
	if(mass_cut_sig > 2.4 || mass_cut_sig < 1.4) isReject_sig = true;
	if(mass_cut_nor > 2.4 || mass_cut_nor < 1.4) isReject_nor = true;
	if(isReject_sig && isReject_nor) continue;
	invmass->push_back(mass_cut_sig);
	invmass->push_back(mass_cut_nor);

	// Build tracks
	std::vector<reco::TransientTrack> t_trks;
	reco::TrackRef trk1 = m[0].innerTrack();
	reco::TrackRef trk2 = m[1].innerTrack();
	t_trks.push_back(theB->build(trk1));
	t_trks.push_back(theB->build(trk2));
	t_trks.push_back(theB->build(&(t[0])));

	// Vertex fit
	KalmanVertexFitter kvf;
	TransientVertex fv = kvf.vertex(t_trks);
	if(fv.isValid()){
	  vf_Valid->push_back(1);
	  vf_tC->push_back(fv.totalChiSquared());
	  vf_dOF->push_back(fv.degreesOfFreedom());
	  vf_nC->push_back(fv.totalChiSquared()/fv.degreesOfFreedom());
	  vf_Prob->push_back(TMath::Prob(fv.totalChiSquared(),(int)fv.degreesOfFreedom()));
	}
	else{
	  vf_Valid->push_back(0);
	  // Vertex Fitting info
	  vf_tC->push_back(-99);
	  vf_dOF->push_back(-99);
	  vf_nC->push_back(-99);
	  vf_Prob->push_back(-99);
	}
	
	// Store all possible triples
	// Muon info
	for(int i=0; i<2; i++){
	  // Kinematics
	  muon_pt->push_back(m[i].muonBestTrack()->pt());
	  muon_eta->push_back(m[i].muonBestTrack()->eta());
	  muon_phi->push_back(m[i].muonBestTrack()->phi());
	  muon_ptErr->push_back(m[i].muonBestTrack()->ptError());
	  muon_q->push_back(m[i].muonBestTrack()->charge());
	  muon_staPt->push_back(m[i].standAloneMuon().isNonnull() ? m[i].standAloneMuon()->pt() : 0);
	  muon_staEta->push_back(m[i].standAloneMuon().isNonnull() ? m[i].standAloneMuon()->eta() : 0);
	  muon_staPhi->push_back(m[i].standAloneMuon().isNonnull() ? m[i].standAloneMuon()->phi() : 0);
	  muon_pfPt->push_back(m[i].pfP4().pt());
	  muon_pfEta->push_back(m[i].pfP4().eta());
	  muon_pfPhi->push_back(m[i].pfP4().phi());

	  // Isolation
	  muon_trkIso->push_back(m[i].isolationR03().sumPt);
	  muon_ecalIso->push_back(m[i].isolationR03().emEt);
	  muon_chHadIso->push_back(m[i].pfIsolationR04().sumChargedHadronPt);
	  muon_gammaIso->push_back(m[i].pfIsolationR04().sumPhotonEt);
	  muon_neuHadIso->push_back(m[i].pfIsolationR04().sumNeutralHadronEt);
	  muon_puIso->push_back(m[i].pfIsolationR04().sumPUPt);

	  // Impact Parameter
	  muon_d0->push_back((-1)*(m[i].muonBestTrack()->dxy(pv.position())));
	  muon_dz->push_back(m[i].muonBestTrack()->dz(pv.position()));

	  const reco::TransientTrack &tt = theB->build(m[i].muonBestTrack());
	  const double thesign = ((-1)*(m[i].muonBestTrack()->dxy(pv.position())) >= 0) ? 1. : -1.;
	  const std::pair<bool,Measurement1D> &ip3d = IPTools::absoluteImpactParameter3D(tt,pv);
	  muon_sip3d->push_back(ip3d.first ? thesign*ip3d.second.value() / ip3d.second.error() : -999.);

	  // Identification
	  muon_tkNchi2->push_back(m[i].innerTrack().isNonnull() ? m[i].innerTrack()->normalizedChi2(): -999.);
	  muon_muNchi2->push_back(m[i].isGlobalMuon() ? m[i].globalTrack()->normalizedChi2() : -999.);
	  muon_trkKink->push_back(m[i].combinedQuality().trkKink);
	  muon_glbKink->push_back(m[i].combinedQuality().glbKink);
	  muon_typeBits->push_back(m[i].type());

	  unsigned int selectorbits=0;
	  if(muon::isGoodMuon(m[i],muon::All))                                    selectorbits |= baconhep::kAll;
	  if(muon::isGoodMuon(m[i],muon::AllGlobalMuons))                         selectorbits |= baconhep::kAllGlobalMuons;
	  if(muon::isGoodMuon(m[i],muon::AllStandAloneMuons))                     selectorbits |= baconhep::kAllStandAloneMuons;
	  if(muon::isGoodMuon(m[i],muon::AllTrackerMuons))                        selectorbits |= baconhep::kAllTrackerMuons;
	  if(muon::isGoodMuon(m[i],muon::TrackerMuonArbitrated))                  selectorbits |= baconhep::kTrackerMuonArbitrated;
	  if(muon::isGoodMuon(m[i],muon::AllArbitrated))                          selectorbits |= baconhep::kAllArbitrated;
	  if(muon::isGoodMuon(m[i],muon::GlobalMuonPromptTight))                  selectorbits |= baconhep::kGlobalMuonPromptTight;
	  if(muon::isGoodMuon(m[i],muon::TMLastStationLoose))                     selectorbits |= baconhep::kTMLastStationLoose;
	  if(muon::isGoodMuon(m[i],muon::TMLastStationTight))                     selectorbits |= baconhep::kTMLastStationTight;
	  if(muon::isGoodMuon(m[i],muon::TM2DCompatibilityLoose))                 selectorbits |= baconhep::kTM2DCompatibilityLoose;
	  if(muon::isGoodMuon(m[i],muon::TM2DCompatibilityTight))                 selectorbits |= baconhep::kTM2DCompatibilityTight;
	  if(muon::isGoodMuon(m[i],muon::TMOneStationLoose))                      selectorbits |= baconhep::kTMOneStationLoose;
	  if(muon::isGoodMuon(m[i],muon::TMOneStationTight))                      selectorbits |= baconhep::kTMOneStationTight;
	  if(muon::isGoodMuon(m[i],muon::TMLastStationOptimizedLowPtLoose))       selectorbits |= baconhep::kTMLastStationOptimizedLowPtLoose;
	  if(muon::isGoodMuon(m[i],muon::TMLastStationOptimizedLowPtTight))       selectorbits |= baconhep::kTMLastStationOptimizedLowPtTight;
	  if(muon::isGoodMuon(m[i],muon::GMTkChiCompatibility))                   selectorbits |= baconhep::kGMTkChiCompatibility;
	  if(muon::isGoodMuon(m[i],muon::GMStaChiCompatibility))                  selectorbits |= baconhep::kGMStaChiCompatibility;
	  if(muon::isGoodMuon(m[i],muon::GMTkKinkTight))                          selectorbits |= baconhep::kGMTkKinkTight;
	  if(muon::isGoodMuon(m[i],muon::TMLastStationAngLoose))                  selectorbits |= baconhep::kTMLastStationAngLoose;
	  if(muon::isGoodMuon(m[i],muon::TMLastStationAngTight))                  selectorbits |= baconhep::kTMLastStationAngTight;
	  if(muon::isGoodMuon(m[i],muon::TMOneStationAngLoose))                   selectorbits |= baconhep::kTMOneStationAngLoose;
	  if(muon::isGoodMuon(m[i],muon::TMOneStationAngTight))                   selectorbits |= baconhep::kTMOneStationAngTight;
	  if(muon::isGoodMuon(m[i],muon::TMLastStationOptimizedBarrelLowPtLoose)) selectorbits |= baconhep::kTMLastStationOptimizedBarrelLowPtLoose;
	  if(muon::isGoodMuon(m[i],muon::TMLastStationOptimizedBarrelLowPtTight)) selectorbits |= baconhep::kTMLastStationOptimizedBarrelLowPtTight;
	  if(muon::isGoodMuon(m[i],muon::RPCMuLoose))                             selectorbits |= baconhep::kRPCMuLoose;
	  muon_selectorBits->push_back(selectorbits);

	  unsigned int pogidbits=0;
	  if(muon::isLooseMuon(m[i]))      pogidbits |= baconhep::kPOGLooseMuon;
	  if(muon::isTightMuon(m[i], pv))  pogidbits |= baconhep::kPOGTightMuon;
	  if(muon::isSoftMuon(m[i], pv))   pogidbits |= baconhep::kPOGSoftMuon;
	  if(muon::isHighPtMuon(m[i], pv)) pogidbits |= baconhep::kPOGHighPtMuon;
	  muon_pogIDBits->push_back(pogidbits);
	    
	  muon_nValidHits->push_back(m[i].isGlobalMuon() ? m[i].globalTrack()->hitPattern().numberOfValidMuonHits(): 0);
	  muon_nTkHits->push_back(m[i].innerTrack().isNonnull() ? m[i].innerTrack()->hitPattern().numberOfValidTrackerHits(): 0);
	  muon_nPixHits->push_back(m[i].innerTrack().isNonnull() ? m[i].innerTrack()->hitPattern().numberOfValidPixelHits() : 0);
	  muon_nTkLayers->push_back(m[i].innerTrack().isNonnull() ? m[i].innerTrack()->hitPattern().trackerLayersWithMeasurement() : 0);
	  muon_nPixLayers->push_back(m[i].innerTrack().isNonnull() ? m[i].innerTrack()->hitPattern().pixelLayersWithMeasurement() : 0); 
	  muon_nMatchStn->push_back(m[i].numberOfMatchedStations());

	  int trkid = -1;
	  if(m[i].innerTrack().isNonnull()) {
	    int trkIndex = -1;
	    for(reco::TrackCollection::const_iterator itTrk = trackCol->begin(); itTrk!=trackCol->end(); ++itTrk) {
	      trkIndex++;
	      if(m[i].innerTrack().get() == &(*itTrk)) {
		trkid = trkIndex; //Check if this is correct+++++++++++++++++++++++++++++++++++
		break;
	      }
	    }
	  }
	  muon_trkID->push_back(trkid);
	  muon_hltMatchBits->push_back(TriggerTools::matchHLT(m[i].eta(), m[i].phi(), triggerRecords, triggerEvent));
	}

	// Track info
	for(int i=0; i<1; i++){
	  // Kinematics
	  muon_pt->push_back(t[i].pt());
	  muon_eta->push_back(t[i].eta());
	  muon_phi->push_back(t[i].phi());
	  muon_ptErr->push_back(t[i].ptError());
	  muon_q->push_back(t[i].charge());
	  muon_staPt->push_back(0);
	  muon_staEta->push_back(0);
	  muon_staPhi->push_back(0);
	  muon_pfPt->push_back(0);
	  muon_pfEta->push_back(0);
	  muon_pfPhi->push_back(0);
	  for(reco::PFCandidateCollection::const_iterator itPF = pfCandCol->begin(); itPF!=pfCandCol->end(); ++itPF){
	    if(itPF->trackRef().isNonnull() && &(t[i]) == itPF->trackRef().get()) {
	      muon_pfPt->push_back(itPF->pt());
	      muon_pfEta->push_back(itPF->eta());
	      muon_pfPhi->push_back(itPF->phi());
	    }
	  }

	  // Isolation
	  muon_trkIso->push_back(-1);
	  muon_ecalIso->push_back(-1);
	  muon_chHadIso->push_back(-1);
	  muon_gammaIso->push_back(-1);
	  muon_neuHadIso->push_back(-1);
	  muon_puIso->push_back(-1);

	  // Impact Parameter
	  muon_d0->push_back((-1)*(t[i].dxy(pv.position())));
	  muon_dz->push_back(t[i].dz(pv.position()));

	  const reco::TransientTrack &tt = theB->build(&(t[i]));
	  const double thesign = ((-1)*(t[i].dxy(pv.position())) >= 0) ? 1. : -1.;
	  const std::pair<bool,Measurement1D> &ip3d = IPTools::absoluteImpactParameter3D(tt,pv);
	  muon_sip3d->push_back(ip3d.first ? thesign*ip3d.second.value() / ip3d.second.error() : -999.);

	  // Identification
	  muon_tkNchi2->push_back(t[i].normalizedChi2());
	  muon_muNchi2->push_back(-999.);
	  muon_trkKink->push_back(0);
	  muon_glbKink->push_back(0);
	  muon_typeBits->push_back(0);
	  muon_selectorBits->push_back(0);
	  muon_pogIDBits->push_back(0);
	  muon_nValidHits->push_back(0);
	  muon_nTkHits->push_back(t[i].hitPattern().numberOfValidTrackerHits());
	  muon_nPixHits->push_back(t[i].hitPattern().numberOfValidPixelHits());
	  muon_nTkLayers->push_back(t[i].hitPattern().trackerLayersWithMeasurement());
	  muon_nPixLayers->push_back(t[i].hitPattern().pixelLayersWithMeasurement()); 
	  muon_nMatchStn->push_back(0);
	  muon_trkID->push_back(0);
	  muon_hltMatchBits->push_back(TriggerTools::matchHLT(t[i].eta(), t[i].phi(), triggerRecords, triggerEvent));
	}

	// Category info
	category->push_back(2);
      }
    }
  }
}
//--------------------------------------End of Scenario B------------------------------------------

// Main Progame ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//--------------------------------------------------------------------------------------------------
bool FillerMuon::fill(std::vector<float> *muon_pt,
		      std::vector<float> *muon_eta,
		      std::vector<float> *muon_phi,
		      std::vector<float> *muon_ptErr,
		      std::vector<float> *muon_staPt,
		      std::vector<float> *muon_staEta,
		      std::vector<float> *muon_staPhi,
		      std::vector<float> *muon_pfPt,
		      std::vector<float> *muon_pfEta,
		      std::vector<float> *muon_pfPhi,
		      std::vector<int>   *muon_q,
		      std::vector<float> *muon_trkIso,
		      std::vector<float> *muon_ecalIso,
		      std::vector<float> *muon_hcalIso,
		      std::vector<float> *muon_chHadIso,
		      std::vector<float> *muon_gammaIso,
		      std::vector<float> *muon_neuHadIso,
		      std::vector<float> *muon_puIso,
		      std::vector<float> *muon_d0,
		      std::vector<float> *muon_dz,
		      std::vector<float> *muon_sip3d,
		      std::vector<float> *muon_tkNchi2,
		      std::vector<float> *muon_muNchi2,
		      std::vector<float> *muon_trkKink,
		      std::vector<float> *muon_glbKink,
		      std::vector<int>   *muon_nValidHits,
		      std::vector<unsigned int>   *muon_typeBits,
		      std::vector<unsigned int>   *muon_selectorBits,
		      std::vector<unsigned int>   *muon_pogIDBits,
		      std::vector<unsigned int>   *muon_nTkHits,
		      std::vector<unsigned int>   *muon_nPixHits,
		      std::vector<unsigned int>   *muon_nTkLayers,
		      std::vector<unsigned int>   *muon_nPixLayers,
		      std::vector<unsigned int>   *muon_nMatchStn,
		      std::vector<int>   *muon_trkID,
		      std::vector<TriggerObjects> *muon_hltMatchBits,
		      std::vector<float> *vf_tC,
		      std::vector<float> *vf_dOF,
		      std::vector<float> *vf_nC,
		      std::vector<float> *vf_Prob,
		      std::vector<int>   *category,
	              std::vector<int>   *vf_Valid,
		      std::vector<float> *invmass,
                      const edm::Event &iEvent, const edm::EventSetup &iSetup, const reco::Vertex &pv, 
		      const std::vector<TriggerRecord> &triggerRecords,
		      const trigger::TriggerEvent &triggerEvent)
{
  
  // Get muon collection
  edm::Handle<reco::MuonCollection> hMuonProduct;
  iEvent.getByToken(fMuonName_token,hMuonProduct);
  assert(hMuonProduct.isValid());
  const reco::MuonCollection *muonCol = hMuonProduct.product();
  
  // Get PF-candidates collection
  edm::Handle<reco::PFCandidateCollection> hPFCandProduct;
  iEvent.getByToken(fPFCandName_token,hPFCandProduct);
  assert(hPFCandProduct.isValid());
  const reco::PFCandidateCollection *pfCandCol = hPFCandProduct.product();
  
  // Get track collection
  edm::Handle<reco::TrackCollection> hTrackProduct;
  iEvent.getByToken(fTrackName_token,hTrackProduct);
  assert(hTrackProduct.isValid());
  const reco::TrackCollection *trackCol = hTrackProduct.product();
  
  // Track builder for computing 3D impact parameter
  edm::ESHandle<TransientTrackBuilder> hTransientTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",hTransientTrackBuilder);
  const TransientTrackBuilder *transientTrackBuilder = hTransientTrackBuilder.product();

  // Muon simple selecion
  std::vector<reco::Muon> muonsel;
  for(std::vector<reco::Muon>::const_iterator itMu=muonCol->begin(); itMu!=muonCol->end(); ++itMu){
    
    // Kinematic cuts
    if(itMu->pt() < fMinPt) continue;
    if(abs(itMu->eta()) > 2.4) continue;

    // Muon type
    if(!(itMu->type() & baconhep::EMuType::kGlobal || itMu->type() & baconhep::EMuType::kTracker))
      continue;

    // Push to muon array
    muonsel.push_back(*itMu);
  }

  // GeneralTrack simple selection
  std::vector<reco::Track> trksel;
  int trkIndex = -1;
  for(reco::TrackCollection::const_iterator itTrk = trackCol->begin(); itTrk!=trackCol->end(); ++itTrk) {
    trkIndex++;
    // check track is not a muon
    bool isMuon = false;
    for(reco::MuonCollection::const_iterator itMu = muonCol->begin(); itMu!=muonCol->end(); ++itMu) {
      if(itMu->innerTrack().isNonnull() && itMu->innerTrack().get() == &(*itTrk)) {
	isMuon = true;
	break;
      }
    }
    if(isMuon) continue;

    // Kinematic cuts
    if(itTrk->pt() < fTrackMinPt) continue;
    if(abs(itTrk->eta()) > 2.4) continue;

    // Push to trk array
    trksel.push_back(*itTrk);
  }

  // There are several possible compositions of the triplets and we will discuss this one by one //
  // If we have more than 3 muons in muonsel, we can form triplets contain: A. 3 muons, B. 2muons C. 1muon and D. 0muon //
  // If we have 2 muons in muonsel, we can form triplets contain: B. 2 muons, C. 1muon D. 0muon //
  // If we have 1 muons in muonsel, we can form triplets contain: C. 1muon D. 0muon //
  // If we do not have muon in muonsel, we can only form triplets of D: 3 trks //
  // Throw triplets of C, D category
  // Implementations are at the beginning of this script // 

  if(muonsel.size() > 2){
    FillA(muonsel, trackCol, muon_pt, muon_eta, muon_phi, muon_ptErr, muon_staPt, muon_staEta, muon_staPhi, muon_pfPt, muon_pfEta, muon_pfPhi, muon_q,
	  muon_trkIso, muon_ecalIso, muon_hcalIso, muon_chHadIso, muon_gammaIso, muon_neuHadIso, muon_puIso, muon_d0, muon_dz, muon_sip3d,
	  muon_tkNchi2, muon_muNchi2, muon_trkKink, muon_glbKink, muon_nValidHits, muon_typeBits, muon_selectorBits, muon_pogIDBits, muon_nTkHits,
	  muon_nPixHits, muon_nTkLayers, muon_nPixLayers, muon_nMatchStn, muon_trkID, muon_hltMatchBits, vf_tC, vf_dOF, vf_nC, vf_Prob, category, vf_Valid, invmass, iEvent, iSetup, pv, triggerRecords, triggerEvent);
    FillB(muonsel, trksel, pfCandCol, trackCol, muon_pt, muon_eta, muon_phi, muon_ptErr, muon_staPt, muon_staEta, muon_staPhi, muon_pfPt, muon_pfEta, muon_pfPhi, muon_q,
	  muon_trkIso, muon_ecalIso, muon_hcalIso, muon_chHadIso, muon_gammaIso, muon_neuHadIso, muon_puIso, muon_d0, muon_dz, muon_sip3d,
	  muon_tkNchi2, muon_muNchi2, muon_trkKink, muon_glbKink, muon_nValidHits, muon_typeBits, muon_selectorBits, muon_pogIDBits, muon_nTkHits,
	  muon_nPixHits, muon_nTkLayers, muon_nPixLayers, muon_nMatchStn, muon_trkID, muon_hltMatchBits, vf_tC, vf_dOF, vf_nC, vf_Prob, category, vf_Valid, invmass, iEvent, iSetup, pv, triggerRecords, triggerEvent);
  }
  if(muonsel.size() == 2){
    FillB(muonsel, trksel, pfCandCol, trackCol, muon_pt, muon_eta, muon_phi, muon_ptErr, muon_staPt, muon_staEta, muon_staPhi, muon_pfPt, muon_pfEta, muon_pfPhi, muon_q,
	  muon_trkIso, muon_ecalIso, muon_hcalIso, muon_chHadIso, muon_gammaIso, muon_neuHadIso, muon_puIso, muon_d0, muon_dz, muon_sip3d,
	  muon_tkNchi2, muon_muNchi2, muon_trkKink, muon_glbKink, muon_nValidHits, muon_typeBits, muon_selectorBits, muon_pogIDBits, muon_nTkHits,
	  muon_nPixHits, muon_nTkLayers, muon_nPixLayers, muon_nMatchStn, muon_trkID, muon_hltMatchBits, vf_tC, vf_dOF, vf_nC, vf_Prob, category, vf_Valid, invmass, iEvent, iSetup, pv, triggerRecords, triggerEvent);
  }
   // End of events processing //

  // Check if no triplet is found in this event
  if(category->size() == 0)
    return false;
  else
    return true;
}
