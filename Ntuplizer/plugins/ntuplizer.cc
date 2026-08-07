#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
//#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include <string>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>

#include "TLorentzVector.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"

namespace MTYPE {
  const char* DSA = "DSA";
  const char* DGL = "DGL";
}

float dxy_value(const reco::GenParticle &p, const reco::Vertex &pv){
    float vx = p.vx();
    float vy = p.vy();
    float phi = p.phi();
    float pv_x = pv.x();
    float pv_y = pv.y();
  
    float dxy = -(vx-pv_x)*sin(phi) + (vy-pv_y)*cos(phi);
    return dxy;
}

bool passTagID(const reco::Track *track, const char* mtype) {
    bool passID = false;
    if (mtype==MTYPE::DSA) {
       if (track->phi() >= -0.8 || track->phi() <= -2.1)        {return passID;}
       if (abs(track->eta()) >= 0.7)                            {return passID;}
       if (track->pt() <= 12.5)                                 {return passID;}
       if (track->ptError()/track->pt() >= 0.2)                 {return passID;}
       if (track->hitPattern().numberOfValidMuonDTHits() <= 30) {return passID;}
       if (track->normalizedChi2() >= 2)                        {return passID;}
       passID = true;
    } else if (mtype==MTYPE::DGL) {
       if (track->phi() >= -0.6 || track->phi() <= -2.6)      {return passID;}
       if (abs(track->eta()) >= 0.9)                          {return passID;}
       if (track->pt() <= 20)                                 {return passID;}
       if (track->ptError()/track->pt() >= 0.3)               {return passID;}
       if (track->hitPattern().numberOfMuonHits() <= 12)      {return passID;}
       if (track->hitPattern().numberOfValidStripHits() <= 5) {return passID;}
       passID = true;
    } else {
      std::cout << "Error (in passTagID): wrong muon type" << std::endl;
    }
    return passID;
}

bool passProbeID(const reco::Track *track, const TVector3 &v_tag, const char* mtype) {
    bool passID = false;
    if (mtype==MTYPE::DSA) {
       if (track->hitPattern().numberOfValidMuonDTHits()+track->hitPattern().numberOfValidMuonCSCHits() <= 12) {return passID;}
       if (track->pt() <= 3.5) {return passID;}
       TVector3 v_probe = TVector3();
       v_probe.SetPtEtaPhi(track->pt(), track->eta(), track->phi());
       if (v_probe.Angle(v_tag) <= 2.1) {return passID;}
       passID = true;
    } else if (mtype==MTYPE::DGL) {
       if (track->pt() <= 20) {return passID;}
       TVector3 v_probe = TVector3();
       v_probe.SetPtEtaPhi(track->pt(), track->eta(), track->phi());
       if (v_probe.Angle(v_tag) <= 2.8) {return passID;}
       passID = true;
    } else {
      std::cout << "Error (in passProbeID): wrong muon type" << std::endl;
    }
    return passID;
}

class ntuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit ntuplizer(const edm::ParameterSet&);
      ~ntuplizer();

      edm::ConsumesCollector iC = consumesCollector();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      edm::ParameterSet parameters;

      bool isData = true;
      bool isAOD  = false;
      //
      // --- Tokens and Handles
      //

      // trigger bits
      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::Handle<edm::TriggerResults> triggerBits;

      // displacedGlobalMuons (reco::Track)
      edm::EDGetTokenT<edm::View<reco::Track> > dglToken;
      edm::Handle<edm::View<reco::Track> > dgls;
      // displacedStandAloneMuons (reco::Track)
      edm::EDGetTokenT<edm::View<reco::Track> > dsaToken;
      edm::Handle<edm::View<reco::Track> > dsas;
      // displacedMuons (reco::Muon // pat::Muon)
      edm::EDGetTokenT<edm::View<reco::Muon> > dmuToken;
      edm::Handle<edm::View<reco::Muon> > dmuons;
      // generator particles (MC only)
      edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticleToken;
      edm::Handle<edm::View<reco::GenParticle> > genParticles;

      // Trigger tags
      std::vector<std::string> HLTPaths_;
      bool triggerPass[200] = {false};

      // Event
      Int_t event = 0;
      Int_t lumiBlock = 0;
      Int_t run = 0;

      // Generator cosmic-muon states (MC only)
      Int_t gen_status1_nMuon = 0;
      bool gen_entry_valid = false;
      Int_t gen_entry_pdgId = 0;
      Float_t gen_entry_charge = 0.;
      Float_t gen_entry_pt = 0.;
      Float_t gen_entry_eta = 0.;
      Float_t gen_entry_phi = 0.;
      Float_t gen_entry_vx = 0.;
      Float_t gen_entry_vy = 0.;
      Float_t gen_entry_vz = 0.;
      Int_t gen_status3_nMuon = 0;
      bool gen_initial_valid = false;
      Int_t gen_initial_pdgId = 0;
      Float_t gen_initial_pt = 0.;
      Float_t gen_initial_eta = 0.;
      Float_t gen_initial_phi = 0.;
      Float_t gen_initial_vx = 0.;
      Float_t gen_initial_vy = 0.;
      Float_t gen_initial_vz = 0.;
      Float_t gen_entry_over_initial_pt = 0.;

      // ----------------------------------
      // displacedMuons
      // ----------------------------------
      Int_t ndmu = 0;
      Int_t dmu_isDSA[200] = {0};
      Int_t dmu_isDGL[200] = {0};
      Int_t dmu_isDTK[200] = {0};
      Int_t dmu_isMatchesValid[200] = {0};
      Int_t dmu_numberOfMatches[200] = {0};
      Int_t dmu_numberOfChambers[200] = {0};
      Int_t dmu_numberOfChambersCSCorDT[200] = {0};
      Int_t dmu_numberOfMatchedStations[200] = {0};
      Int_t dmu_numberOfMatchedRPCLayers[200] = {0};

      Float_t dmu_dsa_pt[200] = {0.};
      Float_t dmu_dsa_eta[200] = {0.};
      Float_t dmu_dsa_phi[200] = {0.};
      Float_t dmu_dsa_ptError[200] = {0.};
      Float_t dmu_dsa_p[200] = {0.};
      Float_t dmu_dsa_qoverp[200] = {0.};
      Float_t dmu_dsa_qoverpError[200] = {0.};
      Float_t dmu_dsa_qoverpt[200] = {0.};
      Float_t dmu_dsa_dxy[200] = {0.};
      Float_t dmu_dsa_dz[200] = {0.};
      Float_t dmu_dsa_refx[200] = {0.};
      Float_t dmu_dsa_refy[200] = {0.};
      Float_t dmu_dsa_refz[200] = {0.};
      Float_t dmu_dsa_innerx[200] = {0.};
      Float_t dmu_dsa_innery[200] = {0.};
      Float_t dmu_dsa_innerz[200] = {0.};
      Float_t dmu_dsa_outerx[200] = {0.};
      Float_t dmu_dsa_outery[200] = {0.};
      Float_t dmu_dsa_outerz[200] = {0.};
      Float_t dmu_dsa_normalizedChi2[200] = {0.};
      Float_t dmu_dsa_charge[200] = {0.};
      Int_t dmu_dsa_collectionIndex[200] = {0};
      Int_t dmu_dsa_side[200] = {0};
      Int_t dmu_dsa_nMuonHits[200] = {0};
      Int_t dmu_dsa_nValidMuonHits[200] = {0};
      Int_t dmu_dsa_nValidMuonDTHits[200] = {0};
      Int_t dmu_dsa_nValidMuonCSCHits[200] = {0};
      Int_t dmu_dsa_nValidMuonRPCHits[200] = {0};
      Int_t dmu_dsa_nValidStripHits[200] = {0};
      Int_t dmu_dsa_nhits[200] = {0};
      Int_t dmu_dsa_dtStationsWithValidHits[200] = {0};
      Int_t dmu_dsa_cscStationsWithValidHits[200] = {0};
      Int_t dmu_dsa_nsegments[200] = {0};
      // Variables for tag and probe
      bool dmu_dsa_passTagID[200] = {false};
      bool dmu_dsa_hasProbe[200] = {false};
      Int_t dmu_dsa_probeID[200] = {0};
      Float_t dmu_dsa_cosAlpha[200] = {0.};

      Int_t evt_dsa_nReco = 0;
      Float_t evt_dsa_pt_1 = 0.;
      Float_t evt_dsa_pt_2 = 0.;
      Int_t evt_dsa_index_1 = -1;
      Int_t evt_dsa_index_2 = -1;
      Int_t evt_dsa_side_1 = 0;
      Int_t evt_dsa_side_2 = 0;
      Float_t evt_dsa_absqoverpt_1 = 0.;
      Float_t evt_dsa_absqoverpt_2 = 0.;
      Float_t evt_dsa_qoverpt_1 = 0.;
      Float_t evt_dsa_qoverpt_2 = 0.;
      Float_t evt_dsa_pt_asymmetry = 0.;
      Float_t evt_dsa_abs_pt_asymmetry = 0.;
      Float_t evt_dsa_cosAlpha_12 = 0.;
      Float_t evt_dsa_residual_12 = 0.;
      Float_t evt_dsa_residual_21 = 0.;
      bool evt_dsa_oppositeSides = false;
      bool evt_dsa_passRawPair = false;
      bool evt_dsa_passQualityPair = false;
      bool evt_dsa_passResolutionPair = false;
      Int_t evt_dsa_index_upper = -1;
      Int_t evt_dsa_index_lower = -1;
      Float_t evt_dsa_pt_upper = 0.;
      Float_t evt_dsa_pt_lower = 0.;
      Float_t evt_dsa_qoverpt_upper = 0.;
      Float_t evt_dsa_qoverpt_lower = 0.;
      Float_t evt_dsa_residual_lower_upper = 0.;
      bool evt_dsa_gen_response_valid = false;
      Float_t evt_dsa_response_upper_gen_entry = 0.;
      Float_t evt_dsa_response_lower_gen_entry = 0.;
      Float_t evt_dsa_gen_cosAlpha_upper = 0.;
      Float_t evt_dsa_gen_cosAlpha_lower_reversed = 0.;

      Float_t dmu_dgl_pt[200] = {0.};
      Int_t dmu_dgl_hasOuterTrack[200] = {0};
      Float_t dmu_dgl_outer_pt[200] = {0.};
      Int_t dmu_dgl_outer_side[200] = {0};
      Float_t dmu_dgl_eta[200] = {0.};
      Float_t dmu_dgl_phi[200] = {0.};
      Float_t dmu_dgl_ptError[200] = {0.};
      Float_t dmu_dgl_dxy[200] = {0.};
      Float_t dmu_dgl_dz[200] = {0.};
      Float_t dmu_dgl_normalizedChi2[200] = {0.};
      Float_t dmu_dgl_charge[200] = {0.};
      Int_t dmu_dgl_nMuonHits[200] = {0};
      Int_t dmu_dgl_nValidMuonHits[200] = {0};
      Int_t dmu_dgl_nValidMuonDTHits[200] = {0};
      Int_t dmu_dgl_nValidMuonCSCHits[200] = {0};
      Int_t dmu_dgl_nValidMuonRPCHits[200] = {0};
      Int_t dmu_dgl_nValidStripHits[200] = {0};
      Int_t dmu_dgl_nhits[200] = {0};
      // Variables for tag and probe
      bool dmu_dgl_passTagID[200] = {false};
      bool dmu_dgl_hasProbe[200] = {false};
      Int_t dmu_dgl_probeID[200] = {0};
      Float_t dmu_dgl_cosAlpha[200] = {0.};

      Float_t dmu_dtk_pt[200] = {0.};
      Float_t dmu_dtk_eta[200] = {0.};
      Float_t dmu_dtk_phi[200] = {0.};
      Float_t dmu_dtk_ptError[200] = {0.};
      Float_t dmu_dtk_dxy[200] = {0.};
      Float_t dmu_dtk_dz[200] = {0.};
      Float_t dmu_dtk_normalizedChi2[200] = {0.};
      Float_t dmu_dtk_charge[200] = {0.};
      Int_t dmu_dtk_nMuonHits[200] = {0};
      Int_t dmu_dtk_nValidMuonHits[200] = {0};
      Int_t dmu_dtk_nValidMuonDTHits[200] = {0};
      Int_t dmu_dtk_nValidMuonCSCHits[200] = {0};
      Int_t dmu_dtk_nValidMuonRPCHits[200] = {0};
      Int_t dmu_dtk_nValidStripHits[200] = {0};
      Int_t dmu_dtk_nhits[200] = {0};

      //
      // --- Output
      //
      std::string output_filename;
      TH1F *counts;
      TFile *file_out;
      TTree *tree_out;

};

// Constructor
ntuplizer::ntuplizer(const edm::ParameterSet& iConfig) {

   usesResource("TFileService");

   parameters = iConfig;

   // Analyzer parameters
   isData = parameters.getParameter<bool>("isData");
   isAOD = parameters.getParameter<bool>("isAOD");  
 
   counts = new TH1F("counts", "", 1, 0, 1);

   dglToken = consumes<edm::View<reco::Track> >  (parameters.getParameter<edm::InputTag>("displacedGlobalCollection"));
   dsaToken = consumes<edm::View<reco::Track> >  (parameters.getParameter<edm::InputTag>("displacedStandAloneCollection"));
   dmuToken = consumes<edm::View<reco::Muon> >  (parameters.getParameter<edm::InputTag>("displacedMuonCollection"));
   if (!isData) {
     genParticleToken = consumes<edm::View<reco::GenParticle> >(
       parameters.getParameter<edm::InputTag>("genParticleCollection"));
   }

   triggerBits_ = consumes<edm::TriggerResults> (parameters.getParameter<edm::InputTag>("bits"));
}


// Destructor
ntuplizer::~ntuplizer() {
}


// beginJob (Before first event)
void ntuplizer::beginJob() {

   std::cout << "Begin Job" << std::endl;

   // Init the file and the TTree
   output_filename = parameters.getParameter<std::string>("nameOfOutput");
   file_out = new TFile(output_filename.c_str(), "RECREATE");
   tree_out = new TTree("Events", "Events");

   // Load HLT paths
   HLTPaths_.push_back("HLT_L2Mu10_NoVertex_NoBPTX3BX");
   HLTPaths_.push_back("HLT_L2Mu10_NoVertex_NoBPTX");

   // TTree branches
   tree_out->Branch("event", &event, "event/I");
   tree_out->Branch("lumiBlock", &lumiBlock, "lumiBlock/I");
   tree_out->Branch("run", &run, "run/I");
   tree_out->Branch("gen_status1_nMuon", &gen_status1_nMuon, "gen_status1_nMuon/I");
   tree_out->Branch("gen_entry_valid", &gen_entry_valid, "gen_entry_valid/O");
   tree_out->Branch("gen_entry_pdgId", &gen_entry_pdgId, "gen_entry_pdgId/I");
   tree_out->Branch("gen_entry_charge", &gen_entry_charge, "gen_entry_charge/F");
   tree_out->Branch("gen_entry_pt", &gen_entry_pt, "gen_entry_pt/F");
   tree_out->Branch("gen_entry_eta", &gen_entry_eta, "gen_entry_eta/F");
   tree_out->Branch("gen_entry_phi", &gen_entry_phi, "gen_entry_phi/F");
   tree_out->Branch("gen_entry_vx", &gen_entry_vx, "gen_entry_vx/F");
   tree_out->Branch("gen_entry_vy", &gen_entry_vy, "gen_entry_vy/F");
   tree_out->Branch("gen_entry_vz", &gen_entry_vz, "gen_entry_vz/F");
   tree_out->Branch("gen_status3_nMuon", &gen_status3_nMuon, "gen_status3_nMuon/I");
   tree_out->Branch("gen_initial_valid", &gen_initial_valid, "gen_initial_valid/O");
   tree_out->Branch("gen_initial_pdgId", &gen_initial_pdgId, "gen_initial_pdgId/I");
   tree_out->Branch("gen_initial_pt", &gen_initial_pt, "gen_initial_pt/F");
   tree_out->Branch("gen_initial_eta", &gen_initial_eta, "gen_initial_eta/F");
   tree_out->Branch("gen_initial_phi", &gen_initial_phi, "gen_initial_phi/F");
   tree_out->Branch("gen_initial_vx", &gen_initial_vx, "gen_initial_vx/F");
   tree_out->Branch("gen_initial_vy", &gen_initial_vy, "gen_initial_vy/F");
   tree_out->Branch("gen_initial_vz", &gen_initial_vz, "gen_initial_vz/F");
   tree_out->Branch("gen_entry_over_initial_pt", &gen_entry_over_initial_pt, "gen_entry_over_initial_pt/F");

   // ----------------------------------
   // displacedMuons
   // ----------------------------------
   tree_out->Branch("ndmu", &ndmu, "ndmu/I");
   tree_out->Branch("dmu_isDSA", dmu_isDSA, "dmu_isDSA[ndmu]/I");
   tree_out->Branch("dmu_isDGL", dmu_isDGL, "dmu_isDGL[ndmu]/I");
   tree_out->Branch("dmu_isDTK", dmu_isDTK, "dmu_isDTK[ndmu]/I");
   tree_out->Branch("dmu_isMatchesValid", dmu_isMatchesValid, "dmu_isMatchesValid[ndmu]/I");
   tree_out->Branch("dmu_numberOfMatches", dmu_numberOfMatches, "dmu_numberOfMatches[ndmu]/I");
   tree_out->Branch("dmu_numberOfChambers", dmu_numberOfChambers, "dmu_numberOfChambers[ndmu]/I");
   tree_out->Branch("dmu_numberOfChambersCSCorDT", dmu_numberOfChambersCSCorDT, "dmu_numberOfChambersCSCorDT[ndmu]/I");
   tree_out->Branch("dmu_numberOfMatchedStations", dmu_numberOfMatchedStations, "dmu_numberOfMatchedStations[ndmu]/I");
   tree_out->Branch("dmu_numberOfMatchedRPCLayers", dmu_numberOfMatchedRPCLayers, "dmu_numberOfMatchedRPCLayers[ndmu]/I");
   // dmu_dsa
   tree_out->Branch("dmu_dsa_pt", dmu_dsa_pt, "dmu_dsa_pt[ndmu]/F");
   tree_out->Branch("dmu_dsa_eta", dmu_dsa_eta, "dmu_dsa_eta[ndmu]/F");
   tree_out->Branch("dmu_dsa_phi", dmu_dsa_phi, "dmu_dsa_phi[ndmu]/F");
   tree_out->Branch("dmu_dsa_ptError", dmu_dsa_ptError, "dmu_dsa_ptError[ndmu]/F");
   tree_out->Branch("dmu_dsa_p", dmu_dsa_p, "dmu_dsa_p[ndmu]/F");
   tree_out->Branch("dmu_dsa_qoverp", dmu_dsa_qoverp, "dmu_dsa_qoverp[ndmu]/F");
   tree_out->Branch("dmu_dsa_qoverpError", dmu_dsa_qoverpError, "dmu_dsa_qoverpError[ndmu]/F");
   tree_out->Branch("dmu_dsa_qoverpt", dmu_dsa_qoverpt, "dmu_dsa_qoverpt[ndmu]/F");
   tree_out->Branch("dmu_dsa_dxy", dmu_dsa_dxy, "dmu_dsa_dxy[ndmu]/F");
   tree_out->Branch("dmu_dsa_dz", dmu_dsa_dz, "dmu_dsa_dz[ndmu]/F");
   tree_out->Branch("dmu_dsa_refx", dmu_dsa_refx, "dmu_dsa_refx[ndmu]/F");
   tree_out->Branch("dmu_dsa_refy", dmu_dsa_refy, "dmu_dsa_refy[ndmu]/F");
   tree_out->Branch("dmu_dsa_refz", dmu_dsa_refz, "dmu_dsa_refz[ndmu]/F");
   tree_out->Branch("dmu_dsa_innerx", dmu_dsa_innerx, "dmu_dsa_innerx[ndmu]/F");
   tree_out->Branch("dmu_dsa_innery", dmu_dsa_innery, "dmu_dsa_innery[ndmu]/F");
   tree_out->Branch("dmu_dsa_innerz", dmu_dsa_innerz, "dmu_dsa_innerz[ndmu]/F");
   tree_out->Branch("dmu_dsa_outerx", dmu_dsa_outerx, "dmu_dsa_outerx[ndmu]/F");
   tree_out->Branch("dmu_dsa_outery", dmu_dsa_outery, "dmu_dsa_outery[ndmu]/F");
   tree_out->Branch("dmu_dsa_outerz", dmu_dsa_outerz, "dmu_dsa_outerz[ndmu]/F");
   tree_out->Branch("dmu_dsa_normalizedChi2", dmu_dsa_normalizedChi2, "dmu_dsa_normalizedChi2[ndmu]/F");
   tree_out->Branch("dmu_dsa_charge", dmu_dsa_charge, "dmu_dsa_charge[ndmu]/F");
   tree_out->Branch("dmu_dsa_collectionIndex", dmu_dsa_collectionIndex, "dmu_dsa_collectionIndex[ndmu]/I");
   tree_out->Branch("dmu_dsa_side", dmu_dsa_side, "dmu_dsa_side[ndmu]/I");
   tree_out->Branch("dmu_dsa_nMuonHits", dmu_dsa_nMuonHits, "dmu_dsa_nMuonHits[ndmu]/I");
   tree_out->Branch("dmu_dsa_nValidMuonHits", dmu_dsa_nValidMuonHits, "dmu_dsa_nValidMuonHits[ndmu]/I");
   tree_out->Branch("dmu_dsa_nValidMuonDTHits", dmu_dsa_nValidMuonDTHits, "dmu_dsa_nValidMuonDTHits[ndmu]/I");
   tree_out->Branch("dmu_dsa_nValidMuonCSCHits", dmu_dsa_nValidMuonCSCHits, "dmu_dsa_nValidMuonCSCHits[ndmu]/I");
   tree_out->Branch("dmu_dsa_nValidMuonRPCHits", dmu_dsa_nValidMuonRPCHits, "dmu_dsa_nValidMuonRPCHits[ndmu]/I");
   tree_out->Branch("dmu_dsa_nValidStripHits", dmu_dsa_nValidStripHits, "dmu_dsa_nValidStripHits[ndmu]/I");
   tree_out->Branch("dmu_dsa_nhits", dmu_dsa_nhits, "dmu_dsa_nhits[ndmu]/I");
   tree_out->Branch("dmu_dsa_dtStationsWithValidHits", dmu_dsa_dtStationsWithValidHits, "dmu_dsa_dtStationsWithValidHits[ndmu]/I");
   tree_out->Branch("dmu_dsa_cscStationsWithValidHits", dmu_dsa_cscStationsWithValidHits, "dmu_dsa_cscStationsWithValidHits[ndmu]/I");
   tree_out->Branch("dmu_dsa_nsegments", dmu_dsa_nsegments, "dmu_dsa_nsegments[ndmu]/I");
   tree_out->Branch("dmu_dsa_passTagID", dmu_dsa_passTagID, "dmu_dsa_passTagID[ndmu]/O");
   tree_out->Branch("dmu_dsa_hasProbe", dmu_dsa_hasProbe, "dmu_dsa_hasProbe[ndmu]/O");
   tree_out->Branch("dmu_dsa_probeID", dmu_dsa_probeID, "dmu_dsa_probeID[ndmu]/I");
   tree_out->Branch("dmu_dsa_cosAlpha", dmu_dsa_cosAlpha, "dmu_dsa_cosAlpha[ndmu]/F");
   tree_out->Branch("evt_dsa_nReco", &evt_dsa_nReco, "evt_dsa_nReco/I");
   tree_out->Branch("evt_dsa_pt_1", &evt_dsa_pt_1, "evt_dsa_pt_1/F");
   tree_out->Branch("evt_dsa_pt_2", &evt_dsa_pt_2, "evt_dsa_pt_2/F");
   tree_out->Branch("evt_dsa_index_1", &evt_dsa_index_1, "evt_dsa_index_1/I");
   tree_out->Branch("evt_dsa_index_2", &evt_dsa_index_2, "evt_dsa_index_2/I");
   tree_out->Branch("evt_dsa_side_1", &evt_dsa_side_1, "evt_dsa_side_1/I");
   tree_out->Branch("evt_dsa_side_2", &evt_dsa_side_2, "evt_dsa_side_2/I");
   tree_out->Branch("evt_dsa_absqoverpt_1", &evt_dsa_absqoverpt_1, "evt_dsa_absqoverpt_1/F");
   tree_out->Branch("evt_dsa_absqoverpt_2", &evt_dsa_absqoverpt_2, "evt_dsa_absqoverpt_2/F");
   tree_out->Branch("evt_dsa_qoverpt_1", &evt_dsa_qoverpt_1, "evt_dsa_qoverpt_1/F");
   tree_out->Branch("evt_dsa_qoverpt_2", &evt_dsa_qoverpt_2, "evt_dsa_qoverpt_2/F");
   tree_out->Branch("evt_dsa_pt_asymmetry", &evt_dsa_pt_asymmetry, "evt_dsa_pt_asymmetry/F");
   tree_out->Branch("evt_dsa_abs_pt_asymmetry", &evt_dsa_abs_pt_asymmetry, "evt_dsa_abs_pt_asymmetry/F");
   tree_out->Branch("evt_dsa_cosAlpha_12", &evt_dsa_cosAlpha_12, "evt_dsa_cosAlpha_12/F");
   tree_out->Branch("evt_dsa_residual_12", &evt_dsa_residual_12, "evt_dsa_residual_12/F");
   tree_out->Branch("evt_dsa_residual_21", &evt_dsa_residual_21, "evt_dsa_residual_21/F");
   tree_out->Branch("evt_dsa_oppositeSides", &evt_dsa_oppositeSides, "evt_dsa_oppositeSides/O");
   tree_out->Branch("evt_dsa_passRawPair", &evt_dsa_passRawPair, "evt_dsa_passRawPair/O");
   tree_out->Branch("evt_dsa_passQualityPair", &evt_dsa_passQualityPair, "evt_dsa_passQualityPair/O");
   tree_out->Branch("evt_dsa_passResolutionPair", &evt_dsa_passResolutionPair, "evt_dsa_passResolutionPair/O");
   tree_out->Branch("evt_dsa_index_upper", &evt_dsa_index_upper, "evt_dsa_index_upper/I");
   tree_out->Branch("evt_dsa_index_lower", &evt_dsa_index_lower, "evt_dsa_index_lower/I");
   tree_out->Branch("evt_dsa_pt_upper", &evt_dsa_pt_upper, "evt_dsa_pt_upper/F");
   tree_out->Branch("evt_dsa_pt_lower", &evt_dsa_pt_lower, "evt_dsa_pt_lower/F");
   tree_out->Branch("evt_dsa_qoverpt_upper", &evt_dsa_qoverpt_upper, "evt_dsa_qoverpt_upper/F");
   tree_out->Branch("evt_dsa_qoverpt_lower", &evt_dsa_qoverpt_lower, "evt_dsa_qoverpt_lower/F");
   tree_out->Branch("evt_dsa_residual_lower_upper", &evt_dsa_residual_lower_upper, "evt_dsa_residual_lower_upper/F");
   tree_out->Branch("evt_dsa_gen_response_valid", &evt_dsa_gen_response_valid, "evt_dsa_gen_response_valid/O");
   tree_out->Branch("evt_dsa_response_upper_gen_entry", &evt_dsa_response_upper_gen_entry, "evt_dsa_response_upper_gen_entry/F");
   tree_out->Branch("evt_dsa_response_lower_gen_entry", &evt_dsa_response_lower_gen_entry, "evt_dsa_response_lower_gen_entry/F");
   tree_out->Branch("evt_dsa_gen_cosAlpha_upper", &evt_dsa_gen_cosAlpha_upper, "evt_dsa_gen_cosAlpha_upper/F");
   tree_out->Branch("evt_dsa_gen_cosAlpha_lower_reversed", &evt_dsa_gen_cosAlpha_lower_reversed, "evt_dsa_gen_cosAlpha_lower_reversed/F");
   // dmu_dgl
   tree_out->Branch("dmu_dgl_pt", dmu_dgl_pt, "dmu_dgl_pt[ndmu]/F");
   tree_out->Branch("dmu_dgl_hasOuterTrack", dmu_dgl_hasOuterTrack, "dmu_dgl_hasOuterTrack[ndmu]/I");
   tree_out->Branch("dmu_dgl_outer_pt", dmu_dgl_outer_pt, "dmu_dgl_outer_pt[ndmu]/F");
   tree_out->Branch("dmu_dgl_outer_side", dmu_dgl_outer_side, "dmu_dgl_outer_side[ndmu]/I");
   tree_out->Branch("dmu_dgl_eta", dmu_dgl_eta, "dmu_dgl_eta[ndmu]/F");
   tree_out->Branch("dmu_dgl_phi", dmu_dgl_phi, "dmu_dgl_phi[ndmu]/F");
   tree_out->Branch("dmu_dgl_ptError", dmu_dgl_ptError, "dmu_dgl_ptError[ndmu]/F");
   tree_out->Branch("dmu_dgl_dxy", dmu_dgl_dxy, "dmu_dgl_dxy[ndmu]/F");
   tree_out->Branch("dmu_dgl_dz", dmu_dgl_dz, "dmu_dgl_dz[ndmu]/F");
   tree_out->Branch("dmu_dgl_normalizedChi2", dmu_dgl_normalizedChi2, "dmu_dgl_normalizedChi2[ndmu]/F");
   tree_out->Branch("dmu_dgl_charge", dmu_dgl_charge, "dmu_dgl_charge[ndmu]/F");
   tree_out->Branch("dmu_dgl_nMuonHits", dmu_dgl_nMuonHits, "dmu_dgl_nMuonHits[ndmu]/I");
   tree_out->Branch("dmu_dgl_nValidMuonHits", dmu_dgl_nValidMuonHits, "dmu_dgl_nValidMuonHits[ndmu]/I");
   tree_out->Branch("dmu_dgl_nValidMuonDTHits", dmu_dgl_nValidMuonDTHits, "dmu_dgl_nValidMuonDTHits[ndmu]/I");
   tree_out->Branch("dmu_dgl_nValidMuonCSCHits", dmu_dgl_nValidMuonCSCHits, "dmu_dgl_nValidMuonCSCHits[ndmu]/I");
   tree_out->Branch("dmu_dgl_nValidMuonRPCHits", dmu_dgl_nValidMuonRPCHits, "dmu_dgl_nValidMuonRPCHits[ndmu]/I");
   tree_out->Branch("dmu_dgl_nValidStripHits", dmu_dgl_nValidStripHits, "dmu_dgl_nValidStripHits[ndmu]/I");
   tree_out->Branch("dmu_dgl_nhits", dmu_dgl_nhits, "dmu_dgl_nhits[ndmu]/I");
   tree_out->Branch("dmu_dgl_passTagID", dmu_dgl_passTagID, "dmu_dgl_passTagID[ndmu]/O");
   tree_out->Branch("dmu_dgl_hasProbe", dmu_dgl_hasProbe, "dmu_dgl_hasProbe[ndmu]/O");
   tree_out->Branch("dmu_dgl_probeID", dmu_dgl_probeID, "dmu_dgl_probeID[ndmu]/I");
   tree_out->Branch("dmu_dgl_cosAlpha", dmu_dgl_cosAlpha, "dmu_dgl_cosAlpha[ndmu]/F");

   // Trigger branches
   for (unsigned int ihlt = 0; ihlt < HLTPaths_.size(); ihlt++) {
     tree_out->Branch(TString(HLTPaths_[ihlt]), &triggerPass[ihlt]);
   }

}

// endJob (After event loop has finished)
void ntuplizer::endJob()
{

    std::cout << "End Job" << std::endl;
    file_out->cd();
    tree_out->Write();
    counts->Write();
    file_out->Close();

}


// fillDescriptions
void ntuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

// Analyze (per event)
void ntuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

   iEvent.getByToken(dglToken, dgls);
   iEvent.getByToken(dsaToken, dsas);
   iEvent.getByToken(dmuToken, dmuons);
   iEvent.getByToken(triggerBits_, triggerBits);
   if (!isData) {
     iEvent.getByToken(genParticleToken, genParticles);
   }

   // Count number of events read
   counts->Fill(0);


   // -> Event info
   event = iEvent.id().event();
   lumiBlock = iEvent.id().luminosityBlock();
   run = iEvent.id().run();

   gen_status1_nMuon = 0;
   gen_entry_valid = false;
   gen_entry_pdgId = 0;
   gen_entry_charge = 0.;
   gen_entry_pt = 0.;
   gen_entry_eta = 0.;
   gen_entry_phi = 0.;
   gen_entry_vx = 0.;
   gen_entry_vy = 0.;
   gen_entry_vz = 0.;
   gen_status3_nMuon = 0;
   gen_initial_valid = false;
   gen_initial_pdgId = 0;
   gen_initial_pt = 0.;
   gen_initial_eta = 0.;
   gen_initial_phi = 0.;
   gen_initial_vx = 0.;
   gen_initial_vy = 0.;
   gen_initial_vz = 0.;
   gen_entry_over_initial_pt = 0.;

   if (!isData && genParticles.isValid()) {
     const reco::GenParticle* entryMuon = nullptr;
     for (const auto& particle : *genParticles) {
       if (std::abs(particle.pdgId()) != 13 || particle.status() != 1) continue;
       gen_status1_nMuon++;
       if (entryMuon == nullptr || particle.vy() > entryMuon->vy()) {
         entryMuon = &particle;
       }
     }

     if (entryMuon != nullptr) {
       gen_entry_valid = true;
       gen_entry_pdgId = entryMuon->pdgId();
       gen_entry_charge = entryMuon->charge();
       gen_entry_pt = entryMuon->pt();
       gen_entry_eta = entryMuon->eta();
       gen_entry_phi = entryMuon->phi();
       gen_entry_vx = entryMuon->vx();
       gen_entry_vy = entryMuon->vy();
       gen_entry_vz = entryMuon->vz();
     }

     const reco::GenParticle* initialMuon = nullptr;
     for (const auto& particle : *genParticles) {
       if (std::abs(particle.pdgId()) != 13 || particle.status() != 3) continue;
       gen_status3_nMuon++;
       if (entryMuon != nullptr && particle.pdgId() != entryMuon->pdgId()) continue;
       if (initialMuon == nullptr || particle.vy() > initialMuon->vy()) {
         initialMuon = &particle;
       }
     }

     if (initialMuon != nullptr) {
       gen_initial_valid = true;
       gen_initial_pdgId = initialMuon->pdgId();
       gen_initial_pt = initialMuon->pt();
       gen_initial_eta = initialMuon->eta();
       gen_initial_phi = initialMuon->phi();
       gen_initial_vx = initialMuon->vx();
       gen_initial_vy = initialMuon->vy();
       gen_initial_vz = initialMuon->vz();
       if (gen_initial_pt > 0.f && gen_entry_valid) {
         gen_entry_over_initial_pt = gen_entry_pt / gen_initial_pt;
       }
     }
   }

   // ----------------------------------
   // displacedMuons Collection
   // ----------------------------------
   ndmu = 0;
   evt_dsa_nReco = 0;
   evt_dsa_pt_1 = 0.;
   evt_dsa_pt_2 = 0.;
   evt_dsa_index_1 = -1;
   evt_dsa_index_2 = -1;
   evt_dsa_side_1 = 0;
   evt_dsa_side_2 = 0;
   evt_dsa_absqoverpt_1 = 0.;
   evt_dsa_absqoverpt_2 = 0.;
   evt_dsa_qoverpt_1 = 0.;
   evt_dsa_qoverpt_2 = 0.;
   evt_dsa_pt_asymmetry = 0.;
   evt_dsa_abs_pt_asymmetry = 0.;
   evt_dsa_cosAlpha_12 = 0.;
   evt_dsa_residual_12 = 0.;
   evt_dsa_residual_21 = 0.;
   evt_dsa_oppositeSides = false;
   evt_dsa_passRawPair = false;
   evt_dsa_passQualityPair = false;
   evt_dsa_passResolutionPair = false;
   evt_dsa_index_upper = -1;
   evt_dsa_index_lower = -1;
   evt_dsa_pt_upper = 0.;
   evt_dsa_pt_lower = 0.;
   evt_dsa_qoverpt_upper = 0.;
   evt_dsa_qoverpt_lower = 0.;
   evt_dsa_residual_lower_upper = 0.;
   evt_dsa_gen_response_valid = false;
   evt_dsa_response_upper_gen_entry = 0.;
   evt_dsa_response_lower_gen_entry = 0.;
   evt_dsa_gen_cosAlpha_upper = 0.;
   evt_dsa_gen_cosAlpha_lower_reversed = 0.;
   for (unsigned int i = 0; i < dmuons->size(); i++) {
     if (ndmu >= 200) break;
     //std::cout << " - - ndmu: " << ndmu << std::endl;
     const reco::Muon& dmuon(dmuons->at(i));
     dmu_isDGL[ndmu] = dmuon.isGlobalMuon();
     dmu_isDSA[ndmu] = dmuon.isStandAloneMuon();
     dmu_isDTK[ndmu] = dmuon.isTrackerMuon();
     dmu_isMatchesValid[ndmu] = dmuon.isMatchesValid();
     dmu_numberOfMatches[ndmu] = dmuon.numberOfMatches();
     dmu_numberOfChambers[ndmu] = dmuon.numberOfChambers();
     dmu_numberOfChambersCSCorDT[ndmu] = dmuon.numberOfChambersCSCorDT();
     dmu_numberOfMatchedStations[ndmu] = dmuon.numberOfMatchedStations();
     dmu_numberOfMatchedRPCLayers[ndmu] = dmuon.numberOfMatchedRPCLayers();

     // Access the DGL track associated to the displacedMuon
     //std::cout << "isGlobalMuon: " << dmuon.isGlobalMuon() << std::endl;
     if ( dmuon.isGlobalMuon() ) {
       const reco::Track* globalTrack = (dmuon.combinedMuon()).get();
       dmu_dgl_pt[ndmu] = globalTrack->pt();
       dmu_dgl_hasOuterTrack[ndmu] = 0;
       dmu_dgl_outer_pt[ndmu] = 0.f;
       dmu_dgl_outer_side[ndmu] = 0;
       const reco::TrackRef outerTrackRef = dmuon.standAloneMuon();
       if (outerTrackRef.isNonnull() && outerTrackRef.isAvailable()) {
         dmu_dgl_hasOuterTrack[ndmu] = 1;
         dmu_dgl_outer_pt[ndmu] = outerTrackRef->pt();
         if (outerTrackRef->extra().isNonnull() && outerTrackRef->extra().isAvailable()) {
           const float midpointY = 0.5f * (
             outerTrackRef->innerPosition().y() + outerTrackRef->outerPosition().y());
           dmu_dgl_outer_side[ndmu] = (midpointY > 0.f ? 1 : (midpointY < 0.f ? -1 : 0));
         }
       }
       dmu_dgl_eta[ndmu] = globalTrack->eta();
       dmu_dgl_phi[ndmu] = globalTrack->phi();
       dmu_dgl_ptError[ndmu] = globalTrack->ptError();
       dmu_dgl_dxy[ndmu] = globalTrack->dxy();
       dmu_dgl_dz[ndmu] = globalTrack->dz();
       dmu_dgl_normalizedChi2[ndmu] = globalTrack->normalizedChi2();
       dmu_dgl_charge[ndmu] = globalTrack->charge();
       dmu_dgl_nMuonHits[ndmu] = globalTrack->hitPattern().numberOfMuonHits();
       dmu_dgl_nValidMuonHits[ndmu] = globalTrack->hitPattern().numberOfValidMuonHits();
       dmu_dgl_nValidMuonDTHits[ndmu] = globalTrack->hitPattern().numberOfValidMuonDTHits();
       dmu_dgl_nValidMuonCSCHits[ndmu] = globalTrack->hitPattern().numberOfValidMuonCSCHits();
       dmu_dgl_nValidMuonRPCHits[ndmu] = globalTrack->hitPattern().numberOfValidMuonRPCHits();
       dmu_dgl_nValidStripHits[ndmu] = globalTrack->hitPattern().numberOfValidStripHits();
       dmu_dgl_nhits[ndmu] = globalTrack->hitPattern().numberOfValidHits();
     } else {
       dmu_dgl_pt[ndmu] = 0;
       dmu_dgl_hasOuterTrack[ndmu] = 0;
       dmu_dgl_outer_pt[ndmu] = 0;
       dmu_dgl_outer_side[ndmu] = 0;
       dmu_dgl_eta[ndmu] = 0;
       dmu_dgl_phi[ndmu] = 0;
       dmu_dgl_ptError[ndmu] = 0;
       dmu_dgl_dxy[ndmu] = 0;
       dmu_dgl_dz[ndmu] = 0;
       dmu_dgl_normalizedChi2[ndmu] = 0;
       dmu_dgl_charge[ndmu] = 0;
       dmu_dgl_nMuonHits[ndmu] = 0;
       dmu_dgl_nValidMuonHits[ndmu] = 0;
       dmu_dgl_nValidMuonDTHits[ndmu] = 0;
       dmu_dgl_nValidMuonCSCHits[ndmu] = 0;
       dmu_dgl_nValidMuonRPCHits[ndmu] = 0;
       dmu_dgl_nValidStripHits[ndmu] = 0;
       dmu_dgl_nhits[ndmu] = 0;
     }    

     // Access the DSA track associated to the displacedMuon
     //std::cout << "isStandAloneMuon: " << dmuon.isStandAloneMuon() << std::endl;
     if ( dmuon.isStandAloneMuon() ) {
       const reco::Track* outerTrack = (dmuon.standAloneMuon()).get();
       dmu_dsa_pt[ndmu] = outerTrack->pt();
       dmu_dsa_eta[ndmu] = outerTrack->eta();
       dmu_dsa_phi[ndmu] = outerTrack->phi();
       dmu_dsa_ptError[ndmu] = outerTrack->ptError();
       dmu_dsa_p[ndmu] = outerTrack->p();
       dmu_dsa_qoverp[ndmu] = outerTrack->qoverp();
       dmu_dsa_qoverpError[ndmu] = outerTrack->qoverpError();
       dmu_dsa_qoverpt[ndmu] = outerTrack->qoverp() * std::cosh(outerTrack->eta());
       dmu_dsa_dxy[ndmu] = outerTrack->dxy();
       dmu_dsa_dz[ndmu] = outerTrack->dz();
       dmu_dsa_refx[ndmu] = outerTrack->referencePoint().x();
       dmu_dsa_refy[ndmu] = outerTrack->referencePoint().y();
       dmu_dsa_refz[ndmu] = outerTrack->referencePoint().z();
       dmu_dsa_normalizedChi2[ndmu] = outerTrack->normalizedChi2();
       dmu_dsa_charge[ndmu] = outerTrack->charge();
       dmu_dsa_collectionIndex[ndmu] = static_cast<Int_t>(i);
       dmu_dsa_innerx[ndmu] = 0.f;
       dmu_dsa_innery[ndmu] = 0.f;
       dmu_dsa_innerz[ndmu] = 0.f;
       dmu_dsa_outerx[ndmu] = 0.f;
       dmu_dsa_outery[ndmu] = 0.f;
       dmu_dsa_outerz[ndmu] = 0.f;
       dmu_dsa_side[ndmu] = 0;
       if (outerTrack->extra().isNonnull() && outerTrack->extra().isAvailable()) {
         dmu_dsa_innerx[ndmu] = outerTrack->innerPosition().x();
         dmu_dsa_innery[ndmu] = outerTrack->innerPosition().y();
         dmu_dsa_innerz[ndmu] = outerTrack->innerPosition().z();
         dmu_dsa_outerx[ndmu] = outerTrack->outerPosition().x();
         dmu_dsa_outery[ndmu] = outerTrack->outerPosition().y();
         dmu_dsa_outerz[ndmu] = outerTrack->outerPosition().z();
         const float midpointY = 0.5f * (dmu_dsa_innery[ndmu] + dmu_dsa_outery[ndmu]);
         dmu_dsa_side[ndmu] = (midpointY > 0.f ? 1 : (midpointY < 0.f ? -1 : 0));
       }
       dmu_dsa_nMuonHits[ndmu] = outerTrack->hitPattern().numberOfMuonHits();
       dmu_dsa_nValidMuonHits[ndmu] = outerTrack->hitPattern().numberOfValidMuonHits();
       dmu_dsa_nValidMuonDTHits[ndmu] = outerTrack->hitPattern().numberOfValidMuonDTHits();
       dmu_dsa_nValidMuonCSCHits[ndmu] = outerTrack->hitPattern().numberOfValidMuonCSCHits();
       dmu_dsa_nValidMuonRPCHits[ndmu] = outerTrack->hitPattern().numberOfValidMuonRPCHits();
       dmu_dsa_nValidStripHits[ndmu] = outerTrack->hitPattern().numberOfValidStripHits();
       dmu_dsa_nhits[ndmu] = outerTrack->hitPattern().numberOfValidHits();
       dmu_dsa_dtStationsWithValidHits[ndmu] = outerTrack->hitPattern().dtStationsWithValidHits();
       dmu_dsa_cscStationsWithValidHits[ndmu] = outerTrack->hitPattern().cscStationsWithValidHits();
       if (isAOD) {
         // Number of DT+CSC segments
         unsigned int nsegments = 0;
         for (trackingRecHit_iterator hit = outerTrack->recHitsBegin(); hit != outerTrack->recHitsEnd(); ++hit) {
           if (!(*hit)->isValid()) continue;
           DetId id = (*hit)->geographicalId();
           if (id.det() != DetId::Muon) continue;
           if (id.subdetId() == MuonSubdetId::DT || id.subdetId() == MuonSubdetId::CSC) {
             nsegments++;
           }
         }
         dmu_dsa_nsegments[ndmu] = nsegments;
       }
     } else {
       dmu_dsa_pt[ndmu] = 0;
       dmu_dsa_eta[ndmu] = 0;
       dmu_dsa_phi[ndmu] = 0;
       dmu_dsa_ptError[ndmu] = 0;
       dmu_dsa_p[ndmu] = 0;
       dmu_dsa_qoverp[ndmu] = 0;
       dmu_dsa_qoverpError[ndmu] = 0;
       dmu_dsa_qoverpt[ndmu] = 0;
       dmu_dsa_dxy[ndmu] = 0;
       dmu_dsa_dz[ndmu] = 0;
       dmu_dsa_refx[ndmu] = 0;
       dmu_dsa_refy[ndmu] = 0;
       dmu_dsa_refz[ndmu] = 0;
       dmu_dsa_innerx[ndmu] = 0;
       dmu_dsa_innery[ndmu] = 0;
       dmu_dsa_innerz[ndmu] = 0;
       dmu_dsa_outerx[ndmu] = 0;
       dmu_dsa_outery[ndmu] = 0;
       dmu_dsa_outerz[ndmu] = 0;
       dmu_dsa_normalizedChi2[ndmu] = 0;
       dmu_dsa_charge[ndmu] = 0;
       dmu_dsa_collectionIndex[ndmu] = -1;
       dmu_dsa_side[ndmu] = 0;
       dmu_dsa_nMuonHits[ndmu] = 0;
       dmu_dsa_nValidMuonHits[ndmu] = 0;
       dmu_dsa_nValidMuonDTHits[ndmu] = 0;
       dmu_dsa_nValidMuonCSCHits[ndmu] = 0;
       dmu_dsa_nValidMuonRPCHits[ndmu] = 0;
       dmu_dsa_nValidStripHits[ndmu] = 0;
       dmu_dsa_nhits[ndmu] = 0;
       dmu_dsa_dtStationsWithValidHits[ndmu] = 0;
       dmu_dsa_cscStationsWithValidHits[ndmu] = 0;
       dmu_dsa_nsegments[ndmu] = 0;
     }

     ndmu++;
     //std::cout << "End muon" << std::endl;
   }

   std::vector<int> dsaRecoIndices;
   for (int i = 0; i < ndmu; ++i) {
     if (dmu_isDSA[i]) {
       dsaRecoIndices.push_back(i);
     }
   }

   evt_dsa_nReco = static_cast<Int_t>(dsaRecoIndices.size());
   if (evt_dsa_nReco == 2) {
     const int idx1 = dsaRecoIndices[0];
     const int idx2 = dsaRecoIndices[1];
     evt_dsa_index_1 = dmu_dsa_collectionIndex[idx1];
     evt_dsa_index_2 = dmu_dsa_collectionIndex[idx2];
     evt_dsa_side_1 = dmu_dsa_side[idx1];
     evt_dsa_side_2 = dmu_dsa_side[idx2];
     evt_dsa_pt_1 = dmu_dsa_pt[idx1];
     evt_dsa_pt_2 = dmu_dsa_pt[idx2];
     evt_dsa_absqoverpt_1 = (dmu_dsa_pt[idx1] != 0.f ? std::abs(dmu_dsa_charge[idx1] / dmu_dsa_pt[idx1]) : 0.f);
     evt_dsa_absqoverpt_2 = (dmu_dsa_pt[idx2] != 0.f ? std::abs(dmu_dsa_charge[idx2] / dmu_dsa_pt[idx2]) : 0.f);
     evt_dsa_qoverpt_1 = dmu_dsa_qoverpt[idx1];
     evt_dsa_qoverpt_2 = dmu_dsa_qoverpt[idx2];
     const float ptSum = evt_dsa_pt_1 + evt_dsa_pt_2;
     evt_dsa_pt_asymmetry = (ptSum != 0.f ? (evt_dsa_pt_2 - evt_dsa_pt_1) / ptSum : 0.f);
     evt_dsa_abs_pt_asymmetry = std::abs(evt_dsa_pt_asymmetry);
     evt_dsa_residual_12 = (evt_dsa_qoverpt_1 != 0.f ?
       (evt_dsa_qoverpt_2 - evt_dsa_qoverpt_1) / evt_dsa_qoverpt_1 : 0.f);
     evt_dsa_residual_21 = (evt_dsa_qoverpt_2 != 0.f ?
       (evt_dsa_qoverpt_1 - evt_dsa_qoverpt_2) / evt_dsa_qoverpt_2 : 0.f);
     TVector3 v1 = TVector3();
     TVector3 v2 = TVector3();
     v1.SetPtEtaPhi(dmu_dsa_pt[idx1], dmu_dsa_eta[idx1], dmu_dsa_phi[idx1]);
     v2.SetPtEtaPhi(dmu_dsa_pt[idx2], dmu_dsa_eta[idx2], dmu_dsa_phi[idx2]);
     evt_dsa_cosAlpha_12 = cos(v1.Angle(v2));

     evt_dsa_oppositeSides = (evt_dsa_side_1 * evt_dsa_side_2 == -1);
     if (evt_dsa_oppositeSides) {
       const int upperIdx = (evt_dsa_side_1 == 1 ? idx1 : idx2);
       const int lowerIdx = (evt_dsa_side_1 == -1 ? idx1 : idx2);
       evt_dsa_index_upper = upperIdx;
       evt_dsa_index_lower = lowerIdx;
       evt_dsa_pt_upper = dmu_dsa_pt[upperIdx];
       evt_dsa_pt_lower = dmu_dsa_pt[lowerIdx];
       evt_dsa_qoverpt_upper = dmu_dsa_qoverpt[upperIdx];
       evt_dsa_qoverpt_lower = dmu_dsa_qoverpt[lowerIdx];
       evt_dsa_residual_lower_upper = (evt_dsa_qoverpt_upper != 0.f ?
         (evt_dsa_qoverpt_lower - evt_dsa_qoverpt_upper) / evt_dsa_qoverpt_upper : 0.f);

       if (gen_entry_valid && gen_entry_pt > 0.f) {
         evt_dsa_gen_response_valid = true;
         evt_dsa_response_upper_gen_entry = (evt_dsa_pt_upper - gen_entry_pt) / gen_entry_pt;
         evt_dsa_response_lower_gen_entry = (evt_dsa_pt_lower - gen_entry_pt) / gen_entry_pt;

         TVector3 genDirection;
         TVector3 upperDirection;
         TVector3 lowerDirection;
         genDirection.SetPtEtaPhi(gen_entry_pt, gen_entry_eta, gen_entry_phi);
         upperDirection.SetPtEtaPhi(
           dmu_dsa_pt[upperIdx], dmu_dsa_eta[upperIdx], dmu_dsa_phi[upperIdx]);
         lowerDirection.SetPtEtaPhi(
           dmu_dsa_pt[lowerIdx], dmu_dsa_eta[lowerIdx], dmu_dsa_phi[lowerIdx]);
         if (genDirection.Mag() > 0. && upperDirection.Mag() > 0.) {
           evt_dsa_gen_cosAlpha_upper = genDirection.Unit().Dot(upperDirection.Unit());
         }
         if (genDirection.Mag() > 0. && lowerDirection.Mag() > 0.) {
           evt_dsa_gen_cosAlpha_lower_reversed = -genDirection.Unit().Dot(lowerDirection.Unit());
         }
       }

       evt_dsa_passRawPair = (evt_dsa_cosAlpha_12 < std::cos(2.1));
       if (evt_dsa_passRawPair) {
         const bool passCommonPt = (dmu_dsa_pt[upperIdx] > 20.f && dmu_dsa_pt[lowerIdx] > 20.f);
         const bool passCommonEta =
           (std::abs(dmu_dsa_eta[upperIdx]) < 0.7f && std::abs(dmu_dsa_eta[lowerIdx]) < 0.7f);
         const bool passCommonDTHits =
           (dmu_dsa_nValidMuonDTHits[upperIdx] >= 31 && dmu_dsa_nValidMuonDTHits[lowerIdx] >= 31);
         const bool passCommonChi2 =
           (dmu_dsa_normalizedChi2[upperIdx] < 5.f && dmu_dsa_normalizedChi2[lowerIdx] < 5.f);
         evt_dsa_passQualityPair = passCommonPt && passCommonEta && passCommonDTHits && passCommonChi2;

         if (evt_dsa_passQualityPair) {
           const float upperRelPtError = dmu_dsa_ptError[upperIdx] / dmu_dsa_pt[upperIdx];
           const float lowerRelPtError = dmu_dsa_ptError[lowerIdx] / dmu_dsa_pt[lowerIdx];
           evt_dsa_passResolutionPair = (upperRelPtError < 0.5f && lowerRelPtError < 0.5f);
         }
       }
     }
   }

   // ----------------------------------
   // Tag and probe code
   // ----------------------------------
   ndmu = 0;
   for (unsigned int i = 0; i < dmuons->size(); i++) {
     if (ndmu >= 200) break;
     const reco::Muon& dmuon(dmuons->at(i));
     // Access the DGL track associated to the displacedMuon
     //std::cout << "isGlobalMuon: " << dmuon.isGlobalMuon() << std::endl;
     if ( dmuon.isGlobalMuon() ) {
       const reco::Track* globalTrack = (dmuon.combinedMuon()).get();
       // Fill tag and probe variables
       //   First, reset the variables
       dmu_dgl_passTagID[ndmu] = false;
       dmu_dgl_hasProbe[ndmu] = false;
       dmu_dgl_probeID[ndmu] = 0;
       dmu_dgl_cosAlpha[ndmu] = 0.;
       // Check if muon passes tag ID
       dmu_dgl_passTagID[ndmu] = passTagID(globalTrack, "DGL");
       if (dmu_dgl_passTagID[ndmu]) {
         // Search probe
         TVector3 v_tag = TVector3();
         v_tag.SetPtEtaPhi(globalTrack->pt(), globalTrack->eta(), globalTrack->phi());
         const reco::Muon *muonProbeTemp = nullptr; // pointer for temporal probe (initialized to nullptr)
         for (unsigned int j = 0; j < dmuons->size(); j++) { // Loop over the rest of the muons
           if (i == j) {continue;}
           const reco::Muon& muonProbeCandidate(dmuons->at(j));
           if (!muonProbeCandidate.isGlobalMuon()) {continue;} // Get only dgls
           const reco::Track *trackProbeCandidate = (muonProbeCandidate.combinedMuon()).get();
           if (passProbeID(trackProbeCandidate, v_tag, "DGL")) { 
             TVector3 v_probe = TVector3();
             v_probe.SetPtEtaPhi(trackProbeCandidate->pt(), trackProbeCandidate->eta(), trackProbeCandidate->phi());
             if (!dmu_dgl_hasProbe[ndmu]) {
               dmu_dgl_hasProbe[ndmu] = true;
               muonProbeTemp = &(dmuons->at(j));
               dmu_dgl_probeID[ndmu] = j;
               dmu_dgl_cosAlpha[ndmu] = cos(v_tag.Angle(v_probe));
             } else {
               const reco::Track *trackProbeTemp = (muonProbeTemp->combinedMuon()).get();
               if (trackProbeCandidate->pt() > trackProbeTemp->pt()) {
                 muonProbeTemp = &(dmuons->at(j));
                 dmu_dgl_probeID[ndmu] = j;
                 dmu_dgl_cosAlpha[ndmu] = cos(v_tag.Angle(v_probe));
               } else {
                 std::cout << ">> Probe candidate " << j << " has lower pt than " << dmu_dgl_probeID[ndmu] << std::endl;
               }
             }
           }
         }
       }
     } else {
       dmu_dgl_passTagID[ndmu] = false;
       dmu_dgl_hasProbe[ndmu] = false;
       dmu_dgl_probeID[ndmu] = 0;
       dmu_dgl_cosAlpha[ndmu] = 0.;
     }
     // Access the DSA track associated to the displacedMuon
     //std::cout << "isStandAloneMuon: " << dmuon.isStandAloneMuon() << std::endl;
     if ( dmuon.isStandAloneMuon() ) {
       const reco::Track* outerTrack = (dmuon.standAloneMuon()).get();
       // Fill tag and probe variables
       //   First, reset the variables
       dmu_dsa_passTagID[ndmu] = false;
       dmu_dsa_hasProbe[ndmu] = false;
       dmu_dsa_probeID[ndmu] = 0;
       dmu_dsa_cosAlpha[ndmu] = 0.;
       // Check if muon passes tag ID
       dmu_dsa_passTagID[ndmu] = passTagID(outerTrack, "DSA");
       if (dmu_dsa_passTagID[ndmu]) {
         // Search probe
         TVector3 v_tag = TVector3();
         v_tag.SetPtEtaPhi(outerTrack->pt(), outerTrack->eta(), outerTrack->phi());
         const reco::Muon *muonProbeTemp = nullptr; // pointer for temporal probe (initialized to nullptr)
         for (unsigned int j = 0; j < dmuons->size(); j++) { // Loop over the rest of the muons
           if (i == j) {continue;}
           const reco::Muon& muonProbeCandidate(dmuons->at(j));
           if (!muonProbeCandidate.isStandAloneMuon()) {continue;} // Get only dsas
           const reco::Track *trackProbeCandidate = (muonProbeCandidate.standAloneMuon()).get();
           if (passProbeID(trackProbeCandidate, v_tag, "DSA")) { 
             TVector3 v_probe = TVector3();
             v_probe.SetPtEtaPhi(trackProbeCandidate->pt(), trackProbeCandidate->eta(), trackProbeCandidate->phi());
             if (!dmu_dsa_hasProbe[ndmu]) {
               dmu_dsa_hasProbe[ndmu] = true;
               muonProbeTemp = &(dmuons->at(j));
               dmu_dsa_probeID[ndmu] = j;
               dmu_dsa_cosAlpha[ndmu] = cos(v_tag.Angle(v_probe));
             } else {
               const reco::Track *trackProbeTemp = (muonProbeTemp->standAloneMuon()).get();
               if (trackProbeCandidate->pt() > trackProbeTemp->pt()) {
                 muonProbeTemp = &(dmuons->at(j));
                 dmu_dsa_probeID[ndmu] = j;
                 dmu_dsa_cosAlpha[ndmu] = cos(v_tag.Angle(v_probe));
               } else {
                 std::cout << ">> Probe candidate " << j << " has lower pt than " << dmu_dsa_probeID[ndmu] << std::endl;
               }
             }
           }
         }
       }
     } else {
       dmu_dsa_passTagID[ndmu] = false;
       dmu_dsa_hasProbe[ndmu] = false;
       dmu_dsa_probeID[ndmu] = 0;
       dmu_dsa_cosAlpha[ndmu] = 0.;
     }
     ndmu++;
   }

   // Check if trigger fired:
   const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
   unsigned int ipath = 0;
   for (auto path : HLTPaths_) {
     std::string path_v = path + "_v";
     // std::cout << path << "\t" << std::endl;
     bool fired = false;
     for (unsigned int itrg = 0; itrg < triggerBits->size(); ++itrg) {
       TString TrigPath = names.triggerName(itrg);
       if (!triggerBits->accept(itrg))
         continue;
       if (!TrigPath.Contains(path_v)){
         continue;
       }
       fired = true;
     }
     triggerPass[ipath] = fired;
     ipath++;
   } 

   //-> Fill tree
   tree_out->Fill();

}

DEFINE_FWK_MODULE(ntuplizer);
