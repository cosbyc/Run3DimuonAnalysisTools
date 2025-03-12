// -*- C++ -*-
//
// Package:    Run3ScoutingAnalysisTools/MuMuGammaTreeMaker
// Class:      MuMuGammaTreeMaker
//
/**\class MuMuGammaTreeMaker MuMuGammaTreeMaker.cc Run3ScoutingAnalysisTools/MuMuGammaTreeMaker/plugins/MuMuGammaTreeMaker.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/

// system include files
#include <memory>
#include <TTree.h>
#include <TLorentzVector.h>
#include <Math/Vector4D.h>
#include "TMath.h"
#include <TPRegexp.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

//#include "TMtuple.hh"
#include "MMGUtils.hh"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

class EEGammaTreeMaker : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns, edm::one::WatchLuminosityBlocks> {
public:
  explicit EEGammaTreeMaker(const edm::ParameterSet&);
  ~EEGammaTreeMaker() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  bool isAncestor(const reco::GenParticle* ancestor, const reco::Candidate* particle);
  bool isMatched(const reco::Candidate* gen_particle, const ROOT::Math::PtEtaPhiMVector* reco_vector, float cand_mass);
  bool isPi0(const std::vector<float>& photonsPt, const std::vector<float>& photonsEta, const std::vector<float>& photonsPhi);

private:
  static constexpr float kMinElectronPt = 3.0f;
  static constexpr float kMinVtxProb = 0.005f; // loose minimum probability as in dimuon HLT path
  static constexpr float kMinPfCandPhotonPt = 1.0f;
  static constexpr float kMaxPfCandDieleDeltaR = 0.5f;
  enum class PDGID : int {
    ID_PHOTON = 22,
    ID_ETA = 221,
    ID_ETA_PRIME = 331,
    ID_OMEGA = 223,
    ID_RHO = 113,
    ID_PHI = 333
  };

  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  const edm::InputTag triggerResultsTag;
  const edm::EDGetTokenT<edm::TriggerResults>             triggerResultsToken;
  const edm::EDGetTokenT<std::vector<pat::Electron> >      electronsToken;
  const edm::EDGetTokenT<std::vector<reco::Vertex> >    primaryVerticesToken;
  const edm::EDGetTokenT<std::vector<reco::VertexCompositePtrCandidate> >    verticesToken;
  const edm::EDGetTokenT<std::vector<pat::Photon> >         photonsToken;
  const edm::EDGetTokenT<std::vector<pat::PackedCandidate> >         pfCandsToken;
  const edm::EDGetTokenT<std::vector<reco::GenParticle> >         prunedGenToken;
  const edm::EDGetTokenT<std::vector<pat::PackedGenParticle> >    packedGenToken;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> esToken;

  std::vector<std::string> triggerPathsVector;
  std::map<std::string, int> triggerPathsMap;

  bool doL1;
  bool doGEN;
  bool doFullGEN;
  triggerExpression::Data triggerCache_;

  edm::InputTag                algInputTag_;
  edm::InputTag                extInputTag_;
  edm::EDGetToken              algToken_;
  std::unique_ptr<l1t::L1TGlobalUtil> l1GtUtils_;
  std::vector<std::string>     l1Seeds_;
  std::vector<bool>            l1Result_;
  std::vector<bool>            hltResult_;

  TTree* tree;

  // event info
  int eventNum;
  int lumiSec;
  int runNum;

  // dielectron variables
  float ee_mass;
  float ee_pt;
  float ee_deltaR;

  // primary vertex variables
  reco::Vertex pv;
  int npv;
  float pv_pos[3];

  // fitted vertex variables
  reco::Vertex vertex;
  float vtx_prob;
  float vtx_pos[3];
  float vtx_posError[3];
  float vtx_chi2;

  // electron variables
  ROOT::Math::PtEtaPhiMVector e_v[2];
  std::vector<bool> e_ID[2];
  size_t   e_idx[2]; // let standard be 0 = e-, 1 = e+
  float e_pfIso[2];
  float e_pt[2];
  float e_eta[2];
  float e_phi[2];
  float e_dxy[2];
  float e_dz[2];
  float e_trkChi2[2];
  float e_trkNdof[2];
  float e_pdgId[2];

  // photon variables
  // slimmed photons
  std::vector<float> slimmedPhoton_pt;
  std::vector<float> slimmedPhoton_eta;
  std::vector<float> slimmedPhoton_phi;
  std::vector<float> slimmedPhoton_mass;
  std::vector<float> slimmedPhoton_sigmaIetaIeta;
  std::vector<float> slimmedPhoton_hOverE;
  std::vector<float> slimmedPhoton_ecalIso;
  std::vector<float> slimmedPhoton_hcalIso;
  std::vector<float> slimmedPhoton_trkIso;
  std::vector<float> slimmedPhoton_r9;
  // pfCand photons
  std::vector<float> pfCandPhoton_deltaR;
  std::vector<float> pfCandPhoton_iso;
  std::vector<float> pfCandPhoton_pt;
  std::vector<float> pfCandPhoton_eta;
  std::vector<float> pfCandPhoton_phi;
  std::vector<float> pfCandPhoton_energy;
  std::vector<float> pfCandPhoton_et;
  std::vector<float> pfCandPhoton_et2;

  /*
  // doGEN variables
  int motherID[2];
  int motherGenID;
  int genID[2];
  int simType[2];
  int simExtType[2];
  std::vector<int> matchedDaughtersIDs;
  std::vector<float> matchedPhotonPt;
  std::vector<float> matchedPhotonEta;
  std::vector<float> matchedPhotonPhi;
  // flags for GEN matched decays
  bool isEta2EE;
  bool isEta2EEGamma;
  bool isEtaPrime2EE;
  bool isEtaPrime2EEGamma;
  bool isOmega2EE;
  bool isOmega2Pi0EE;
  bool isRho2EE;
  // bool isRho2Pi0EE; //not observed?
  bool isPhi2EE;
  bool isPhi2KK;
  */

  int motherID[2];
  int grandMotherID[2];
  int simType[2];
  int simExtType[2];
  int matchedGenIdx;

  std::vector< std::vector<int> > gen_matchedDaughtersIDs;
  std::vector< std::vector<float> > gen_matchedPhotonPt;
  std::vector< std::vector<float> > gen_matchedPhotonEta;
  std::vector< std::vector<float> > gen_matchedPhotonPhi;

  std::vector<int> gen_motherID;
  std::vector<int> gen_grandMotherID;
  std::vector<float> gen_ee_mass;
  std::vector<float> gen_e_pt[2];
  std::vector<float> gen_e_eta[2];
  std::vector<float> gen_e_phi[2];
  std::vector<bool> gen_e_recoMatch[2];

  std::vector<bool> gen_isEta2EE;
  std::vector<bool> gen_isEta2EEGamma;
  std::vector<bool> gen_isEtaPrime2EE;
  std::vector<bool> gen_isEtaPrime2EEGamma;
  std::vector<bool> gen_isOmega2EE;
  std::vector<bool> gen_isOmega2Pi0EE;
  std::vector<bool> gen_isRho2EE;
  // std::vector<bool> isRho2Pi0EE;
  std::vector<bool> gen_isPhi2EE;
  std::vector<bool> gen_isPhi2KK;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
EEGammaTreeMaker::EEGammaTreeMaker(const edm::ParameterSet& iConfig):
    triggerResultsTag        (iConfig.getParameter<edm::InputTag>("triggerresults")),
    triggerResultsToken      (consumes<edm::TriggerResults>                    (triggerResultsTag)),
    electronsToken               (consumes<std::vector<pat::Electron> >             (iConfig.getParameter<edm::InputTag>("electrons"))),
    primaryVerticesToken     (consumes<std::vector<reco::Vertex> >           (iConfig.getParameter<edm::InputTag>("primaryVertices"))),
    verticesToken            (consumes<std::vector<reco::VertexCompositePtrCandidate> >           (iConfig.getParameter<edm::InputTag>("displacedVertices"))),
    photonsToken             (consumes<std::vector<pat::Photon> >         (iConfig.getParameter<edm::InputTag>("photons"))),
    pfCandsToken             (consumes<std::vector<pat::PackedCandidate> >         (iConfig.getParameter<edm::InputTag>("pfcands"))),
    prunedGenToken  (consumes<std::vector<reco::GenParticle> >      (iConfig.getParameter<edm::InputTag>("prunedGenParticles"))),
    packedGenToken  (consumes<std::vector<pat::PackedGenParticle> > (iConfig.getParameter<edm::InputTag>("packedGenParticles"))),
    esToken(esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"))),
    doL1                     (iConfig.existsAs<bool>("doL1")               ?    iConfig.getParameter<bool>  ("doL1")            : false),
    doGEN                    (iConfig.existsAs<bool>("doGEN")              ?    iConfig.getParameter<bool>  ("doGEN")           : false),
    doFullGEN                (iConfig.existsAs<bool>("doFullGEN")          ?    iConfig.getParameter<bool>  ("doFullGEN")       : false)
{
    usesResource("TFileService");
    if (doL1) {
        algInputTag_ = iConfig.getParameter<edm::InputTag>("AlgInputTag");
        extInputTag_ = iConfig.getParameter<edm::InputTag>("l1tExtBlkInputTag");
        algToken_ = consumes<BXVector<GlobalAlgBlk>>(algInputTag_);
        l1Seeds_ = iConfig.getParameter<std::vector<std::string> >("l1Seeds");
        l1GtUtils_ = std::make_unique<l1t::L1TGlobalUtil>(iConfig, consumesCollector(), *this, algInputTag_, extInputTag_, l1t::UseEventSetupIn::Event);
    }
    else {
        l1Seeds_ = std::vector<std::string>();
        l1GtUtils_ = 0;
    }
}

EEGammaTreeMaker::~EEGammaTreeMaker() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

//Check recursively if any ancestor of particle is the given one
bool EEGammaTreeMaker::isAncestor(const reco::GenParticle* ancestor, const reco::Candidate* particle) {
  // cast to the base class to make direct comparison
  const reco::Candidate* ancestorPtr = ancestor;
  //particle is already the ancestor
          if(ancestorPtr == particle ) return true;

  //otherwise loop on mothers, if any and return true if the ancestor is found
          for(size_t i=0;i< particle->numberOfMothers();i++)
          {
            const reco::Candidate* motherPtr = particle->mother(i);
            if(isAncestor(ancestor, motherPtr)) return true;
          }
  //if we did not return yet, then particle and ancestor are not relatives
          return false;
}

//check if invariant mass of 2 photons is close to pi0
bool EEGammaTreeMaker::isPi0(const std::vector<float>& photonsPt, const std::vector<float>& photonsEta, const std::vector<float>& photonsPhi) {
  ROOT::Math::PtEtaPhiMVector photon1( photonsPt[0], photonsEta[0], photonsPhi[0], 0.);
  ROOT::Math::PtEtaPhiMVector photon2( photonsPt[1], photonsEta[1], photonsPhi[1], 0.);
  ROOT::Math::PtEtaPhiMVector diPhoton;
  float diPhotonMass;

  diPhoton = photon1 + photon2;
  diPhotonMass = diPhoton.M();

  return ((diPhotonMass >= (PI0_MASS - PI0_MASS_SHIFT)) and (diPhotonMass <= (PI0_MASS + PI0_MASS_SHIFT)));

}

// check if a vector of gen particle and reco vector match (by dR)
bool EEGammaTreeMaker::isMatched(const reco::Candidate* gen_particle, const ROOT::Math::PtEtaPhiMVector* reco_vector, float cand_mass) {
  bool is_matched = false;
  ROOT::Math::PtEtaPhiMVector gen_vec( gen_particle->pt(), gen_particle->eta(), gen_particle->phi(), cand_mass);
  float dr_gen_reco = ROOT::Math::VectorUtil::DeltaR(gen_vec,*reco_vector);
  is_matched = dr_gen_reco <= MIN_DR_TRUTH;

  return is_matched;
}


// ------------ method called for each event  ------------
void EEGammaTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;
  using namespace reco;

  // RECO
  bool fillTree = false;
  Handle<vector<pat::Electron> > electronsH;
  iEvent.getByToken(electronsToken, electronsH);
  if ( electronsH->size() > 1 ) { // check for 2+ electrons
    Handle<vector<reco::Vertex> > primaryVerticesH;
    iEvent.getByToken(primaryVerticesToken, primaryVerticesH);
    npv = primaryVerticesH->size();
    //std::cout << "npv " << npv << std::endl;
    if ( npv > 0 ) { // check for primary vertex
      pv = *primaryVerticesH->begin();

      float bestProbVtx = 0;
      edm::ESHandle<TransientTrackBuilder> theB = iSetup.getHandle(esToken);
      KalmanVertexFitter kvf(true);
      TransientVertex tv;
      for ( unsigned int iElectron = 0; iElectron < electronsH->size(); iElectron++ ) {
        if ( (electronsH->at(iElectron)).pt() > kMinElectronPt ) {
          for ( unsigned int jElectron = iElectron+1; jElectron < electronsH->size(); jElectron++ ) {
	          if ( (!electronsH->at(jElectron).gsfTrack().isNull()) && (!electronsH->at(iElectron).gsfTrack().isNull())) {
	            if ( (electronsH->at(jElectron)).pt() > kMinElectronPt && ( (electronsH->at(iElectron)).charge() * (electronsH->at(jElectron)).charge() < 0 ) ) {
              reco::Track part_1 = *((electronsH->at(iElectron)).gsfTrack());
	            reco::Track part_2 = *((electronsH->at(jElectron)).gsfTrack());
              vector<reco::TransientTrack> transient_tracks{};
              transient_tracks.push_back(theB->build(&part_1));
              transient_tracks.push_back(theB->build(&part_2));
              tv = kvf.vertex(transient_tracks);

              /*reco::Track part_1 = *((electronsH->at(iElectron)).track());
              reco::Track part_2 = *((electronsH->at(jElectron)).track());
              vector<reco::TransientTrack> transient_tracks{};
              transient_tracks.push_back(theB->build(fix_track(&part_1)));
              transient_tracks.push_back(theB->build(fix_track(&part_2)));
              tv = kvf.vertex(transient_tracks);*/
	      
              if (!tv.isValid()) {
                std::cout << "ij " << iElectron << jElectron << "Vertex not valid." << std::endl;
              } else {
                reco::Vertex currentVtx = reco::Vertex(tv);
                float currentVtxProb = TMath::Prob( currentVtx.chi2() , currentVtx.ndof() );
                //std::cout << "  currentVtxProb " << currentVtxProb << std::endl;
                if (currentVtxProb > kMinVtxProb && currentVtxProb > bestProbVtx) {
                  vertex = currentVtx;
                  bestProbVtx = currentVtxProb;
                  bool iElectronIsPositive = (electronsH->at(iElectron)).charge() > 0;
                  e_idx[iElectronIsPositive] = iElectron; // if iElectron is negative/positive, save to e_idx[0]/e_idx[1]
                  e_idx[!iElectronIsPositive] = jElectron; // if iElectron is negative/positive, jElectron is positive/negative, save jElectron to e_idx[1]/e_idx[0]
                }
              }
	          }
	        }
        }
      }
    }      
      
      //if (bestProbVtx > kMinVtxProb) { // will fill tree if true
      if (true) { // will fill tree if true
        //std::cout<< "fill tree"<<std::endl;
        fillTree = true;
        /*vtx_prob = bestProbVtx;
        vtx_chi2 = vertex.normalizedChi2();
        vtx_pos[0] = vertex.x();
        vtx_pos[1] = vertex.y();
        vtx_pos[2] = vertex.z();
        vtx_posError[0] = vertex.xError();
        vtx_posError[1] = vertex.yError();
        vtx_posError[2] = vertex.zError();

        pv_pos[0] = pv.x();
        pv_pos[1] = pv.y();
        pv_pos[2] = pv.z();
	*/
        eventNum = iEvent.id().event();
        lumiSec = iEvent.luminosityBlock();
        runNum = iEvent.id().run();

        for ( int i : {0, 1} ) {
	  if ( e_idx[i] < electronsH->size() ) continue;
          auto iE = electronsH->at(e_idx[i]);

          e_pt[i] = iE.pt();
          e_eta[i] = iE.eta();
          e_phi[i] = iE.phi();
          e_v[i].SetCoordinates(e_pt[i],e_eta[i],e_phi[i],ele_mass);

          // Compute the pfIso for the electron. Note: PUCorr = 0.5*electrons_iter->puChargedHadronIso()                                                                              
          // -----------------------------------------------------------------------------------                                                                              
          e_pfIso[i] = (iE.chargedHadronIso() + std::max(iE.photonIso() + iE.neutralHadronIso() - 0.5 * iE.puChargedHadronIso() , 0.0)) / iE.pt();

          e_dxy[i] = iE.track()->dxy();
          e_dz[i] = iE.track()->dz();
          e_trkChi2[i] = iE.track()->chi2();
          e_trkNdof[i] = iE.track()->ndof();
	  
          e_pdgId[i] = iE.pdgId();

          e_ID[i].clear();
          e_ID[i].push_back(iE.electronID("simpleEleId95relIso"));
          e_ID[i].push_back(iE.electronID("simpleEleId90relIso"));
          e_ID[i].push_back(iE.electronID("simpleEleId85relIso"));
          e_ID[i].push_back(iE.electronID("simpleEleId80relIso"));
          e_ID[i].push_back(iE.electronID("simpleEleId70relIso"));
          e_ID[i].push_back(iE.electronID("simpleEleId60relIso"));

	}

        ROOT::Math::PtEtaPhiMVector diele = e_v[0] + e_v[1];
        ee_mass = diele.M();
        ee_pt = diele.Pt();
        ee_deltaR = ROOT::Math::VectorUtil::DeltaR(e_v[0], e_v[1]);

        Handle<vector<pat::Photon> > photonsH;
        iEvent.getByToken(photonsToken, photonsH);
        slimmedPhoton_pt.clear();
        slimmedPhoton_eta.clear();
        slimmedPhoton_phi.clear();
        slimmedPhoton_mass.clear();
        slimmedPhoton_sigmaIetaIeta.clear();
        slimmedPhoton_hOverE.clear();
        slimmedPhoton_ecalIso.clear();
        slimmedPhoton_hcalIso.clear();
        slimmedPhoton_trkIso.clear();
        slimmedPhoton_r9.clear();
        for (auto photons_iter = photonsH->begin(); photons_iter != photonsH->end(); ++photons_iter) {
          slimmedPhoton_pt.push_back(photons_iter->pt());
          slimmedPhoton_eta.push_back(photons_iter->eta());
          slimmedPhoton_phi.push_back(photons_iter->phi());
          slimmedPhoton_mass.push_back(photons_iter->mass());
          slimmedPhoton_sigmaIetaIeta.push_back(photons_iter->sigmaIetaIeta());
          slimmedPhoton_hOverE.push_back(photons_iter->hadronicOverEm());
          slimmedPhoton_ecalIso.push_back(photons_iter->ecalIso());
          slimmedPhoton_hcalIso.push_back(photons_iter->hcalIso());
          slimmedPhoton_trkIso.push_back(photons_iter->trackIso());
          slimmedPhoton_r9.push_back(photons_iter->r9());
        }

        Handle<vector<pat::PackedCandidate> > pfCandH;
        iEvent.getByToken(pfCandsToken, pfCandH);
        pfCandPhoton_deltaR.clear();
        pfCandPhoton_iso.clear();
        pfCandPhoton_pt.clear();
        pfCandPhoton_eta.clear();
        pfCandPhoton_phi.clear();
        pfCandPhoton_energy.clear();
        pfCandPhoton_et.clear();
        pfCandPhoton_et2.clear();
        for (auto pfCand_iter = pfCandH->begin(); pfCand_iter != pfCandH->end(); ++pfCand_iter) {
          if (abs(pfCand_iter->pdgId()) != 22 or pfCand_iter->pt() < kMinPfCandPhotonPt) continue;
          float pfCandDieleDr = deltaR(diele.Eta(), diele.Phi(), pfCand_iter->eta(), pfCand_iter->phi());
          if (pfCandDieleDr < kMaxPfCandDieleDeltaR) {
            pfCandPhoton_deltaR.push_back(pfCandDieleDr);
            pfCandPhoton_iso.push_back(photonPfIso03(*pfCand_iter,pfCandH)/pfCand_iter->pt());
            pfCandPhoton_pt.push_back(pfCand_iter->pt());
            pfCandPhoton_eta.push_back(pfCand_iter->eta());
            pfCandPhoton_phi.push_back(pfCand_iter->phi());
            pfCandPhoton_energy.push_back(pfCand_iter->energy());
            pfCandPhoton_et.push_back(pfCand_iter->et());
            pfCandPhoton_et2.push_back(pfCand_iter->et2());
          }
        }

        l1Result_.clear();
        if (doL1) {
            l1GtUtils_->retrieveL1(iEvent,iSetup,algToken_);
            for( unsigned int iseed = 0; iseed < l1Seeds_.size(); iseed++ ) {
                bool l1htbit = 0;
                l1GtUtils_->getFinalDecisionByName(string(l1Seeds_[iseed]), l1htbit);
                l1Result_.push_back( l1htbit );
            }
        }

        Handle<TriggerResults> triggerResultsH;
        iEvent.getByToken(triggerResultsToken, triggerResultsH);
        hltResult_.clear();
        for (size_t i = 0; i < triggerPathsVector.size(); i++) {
            hltResult_.push_back(triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]]));
        }
      }
    }
  }

  bool fillGenTree = true;
  if ( doFullGEN or ( doGEN and fillTree ) ) {

    gen_motherID.clear();
    gen_grandMotherID.clear();
    matchedGenIdx = -1;
    gen_e_recoMatch[0].clear();
    gen_e_recoMatch[1].clear();

    gen_ee_mass.clear();
    gen_e_pt[0].clear();
    gen_e_pt[1].clear();
    gen_e_eta[0].clear();
    gen_e_eta[1].clear();
    gen_e_phi[0].clear();
    gen_e_phi[1].clear();

    gen_isEta2EE.clear();
    gen_isEta2EEGamma.clear();
    gen_isEtaPrime2EE.clear();
    gen_isEtaPrime2EEGamma.clear();
    gen_isOmega2EE.clear();
    gen_isOmega2Pi0EE.clear();
    gen_isRho2EE.clear();
    gen_isPhi2EE.clear();
    gen_isPhi2KK.clear();

    gen_matchedDaughtersIDs.clear();
    gen_matchedPhotonPt.clear();
    gen_matchedPhotonEta.clear();
    gen_matchedPhotonPhi.clear();

    bool isEta2EE            = false;
    bool isEta2EEGamma       = false;
    bool isEtaPrime2EE       = false;
    bool isEtaPrime2EEGamma  = false;
    bool isOmega2EE          = false;
    bool isOmega2Pi0EE       = false;
    bool isRho2EE            = false;
    bool isPhi2EE            = false;
    bool isPhi2KK            = false;

    std::vector<int> matchedDaughtersIDs;
    std::vector<float> matchedPhotonPt;
    std::vector<float> matchedPhotonEta;
    std::vector<float> matchedPhotonPhi;

    Handle<vector<reco::GenParticle> > prunedGenParticles;
    iEvent.getByToken(prunedGenToken, prunedGenParticles);

    Handle<vector<pat::PackedGenParticle> > packedGenParticles;
    iEvent.getByToken(packedGenToken, packedGenParticles);

    for (auto genp = prunedGenParticles->begin(); genp != prunedGenParticles->end(); ++genp) {            
      if (abs(genp->pdgId())==221 or abs(genp->pdgId())==113 or abs(genp->pdgId())==223 or abs(genp->pdgId())==331 or abs(genp->pdgId()==333)) {
      // std::cout<<"genp id: "<<genp->pdgId()<<" pt: "<<genp->pt()<<" eta: "<<genp->eta()<<" status: "<<genp->status();
        int nDaughterElectrons = 0;
        int ep_idx = 99;
        int em_idx = 99;

        for (int i=0; i<(int)genp->numberOfDaughters(); ++i) {

          // check if electrons and save idx of daughters
          if (genp->daughter(i)->pdgId()==11){
            em_idx = i;
            nDaughterElectrons+=1;
          } else if  (genp->daughter(i)->pdgId()== -11){
            ep_idx = i;
            nDaughterElectrons+=1;
          }

        }
        if (nDaughterElectrons==2) {
          matchedDaughtersIDs.clear();
          // try to match pair of reco electrons with a pair of daughters
          auto daught_ep =  &(*genp->daughter(ep_idx));
          auto daught_em =  &(*genp->daughter(em_idx));
          
          // std::cout<<"Mother with 2 electrons ID == " << genp->pdgId() <<std::endl;
          for (int i=0; i<(int)genp->numberOfDaughters(); ++i) {
            matchedDaughtersIDs.push_back(genp->daughter(i)->pdgId());
            // std::cout<<"Daughter "<< i << "  ID == " << genp->daughter(i)->pdgId()<<std::endl;
          }

          bool isMatched1 = false;
          bool isMatched2 = false;
          if (fillTree) {
            isMatched1 = isMatched(daught_ep, &e_v[1], ele_mass);
            isMatched2 = isMatched(daught_em, &e_v[0], ele_mass);
          }

          if ( ( doFullGEN and daught_em->pt() > 2 and daught_ep->pt() > 2 and abs(daught_em->eta()) < 2.5 and abs(daught_ep->eta()) < 2.5 ) or ( isMatched1 and isMatched2 ) ){

            gen_matchedDaughtersIDs.push_back(matchedDaughtersIDs);

            if (isMatched1 and isMatched2) matchedGenIdx = gen_motherID.size();

            gen_e_recoMatch[0].push_back(isMatched1);
            gen_e_recoMatch[1].push_back(isMatched2);

            gen_motherID.push_back(genp->pdgId());
            //gen_grandMotherID.push_back(daught_ep->grandMotherPdgId());


            //std::cout << "ep gm: " << daught_ep->grandMotherPdgId() << " em gm: " << daught_em->grandMotherPdgId() << " grandMotherID1: " << grandMotherID[0] << " grandMotherID2: " << grandMotherID[1] << std::endl;


            ROOT::Math::PtEtaPhiMVector gen_e1_v(daught_em->pt(), daught_em->eta(), daught_em->phi(), ele_mass);
            ROOT::Math::PtEtaPhiMVector gen_e2_v(daught_ep->pt(), daught_ep->eta(), daught_ep->phi(), ele_mass);
            gen_ee_mass.push_back( (gen_e1_v + gen_e2_v).M() );
            gen_e_pt[0].push_back(daught_em->pt());
            gen_e_pt[1].push_back(daught_ep->pt());
            gen_e_eta[0].push_back(daught_em->eta());
            gen_e_eta[1].push_back(daught_ep->eta());
            gen_e_phi[0].push_back(daught_em->phi());
            gen_e_phi[1].push_back(daught_ep->phi());

            // std::cout<<"Matched "<< daught_ep->pt() << "  " << ep_reco->Pt() << "  " << daught_em->pt() << "  " << em_reco->Pt() << "for mother ID == " << motherGenID <<std::endl;
            //check for photons
            int nPhotons = 0;
            matchedPhotonPt.clear();
            matchedPhotonEta.clear();
            matchedPhotonPhi.clear();
            // look for photons
            for (auto pgp = packedGenParticles->begin(); pgp != packedGenParticles->end(); ++pgp) {
              // derefference and cast to the base class to make direct comparison
              const reco::Candidate* pgpPtr = &(*pgp);
              if (pgp->pdgId()==22 and isAncestor( &(*genp) , pgpPtr)) {
                  //std::cout<<"packed photon "<<pgp->pt()<<" mother pt: "<<pgp->motherRef()->pt()<<std::endl;
                  matchedPhotonPt.push_back( pgp->pt() );
                  matchedPhotonEta.push_back( pgp->eta() );
                  matchedPhotonPhi.push_back( pgp->phi() );
                  nPhotons++;   
              }
            }
            gen_matchedPhotonPt.push_back(matchedPhotonPt);
            gen_matchedPhotonEta.push_back(matchedPhotonEta);
            gen_matchedPhotonPhi.push_back(matchedPhotonPhi);
            bool aPi0 = false;
            // check if 2 photons make a pi0
            if (nPhotons == 2){
              aPi0 = isPi0( matchedPhotonPt, matchedPhotonEta, matchedPhotonPhi);
              // if (aPi0) std::cout<<"Matched pi0 to two photons!" <<std::endl;
            }
            // isPhi2KK
            if      ((abs(genp->pdgId()) == 221) and (nPhotons == 0)) isEta2EE = true;
            else if ((abs(genp->pdgId()) == 221) and (nPhotons == 1))    isEta2EEGamma = true;
            else if ((abs(genp->pdgId()) == 331) and (nPhotons == 0)) isEtaPrime2EE = true;
            else if ((abs(genp->pdgId()) == 331) and (nPhotons == 1))    isEtaPrime2EEGamma = true;
            else if ((abs(genp->pdgId()) == 223) and (nPhotons == 0)) isOmega2EE = true;
            else if ((abs(genp->pdgId()) == 223) and aPi0)    isOmega2Pi0EE = true;      
            else if ((abs(genp->pdgId()) == 113) and (nPhotons == 0)) isRho2EE = true;
            else if ((abs(genp->pdgId()) == 333) and (nPhotons == 0)) isPhi2EE = true;

            gen_isEta2EE.push_back(isEta2EE);
            gen_isEta2EEGamma.push_back(isEta2EEGamma);
            gen_isEtaPrime2EE.push_back(isEtaPrime2EE);
            gen_isEtaPrime2EEGamma.push_back(isEtaPrime2EEGamma);
            gen_isOmega2EE.push_back(isOmega2EE);
            gen_isOmega2Pi0EE.push_back(isOmega2Pi0EE);
            gen_isRho2EE.push_back(isRho2EE);
            gen_isPhi2EE.push_back(isPhi2EE);
            gen_isPhi2KK.push_back(isPhi2KK);

            eventNum = iEvent.id().event();
            lumiSec = iEvent.luminosityBlock();
            runNum = iEvent.id().run();
            fillGenTree = true;
          }
        }
      }
    }
  }

  if (fillTree and fillGenTree) {
    tree->Fill();
  } else if (fillTree) {
    motherID[0]=0;
    motherID[1]=0;
    grandMotherID[0]=0;
    grandMotherID[1]=0;
    simType[0]=0;
    simType[1]=0;
    simExtType[0]=0;
    simExtType[1]=0;

    gen_motherID.clear();
    gen_grandMotherID.clear();
    matchedGenIdx = -1;
    gen_e_recoMatch[0].clear();
    gen_e_recoMatch[1].clear();

    gen_ee_mass.clear();
    gen_e_pt[0].clear();
    gen_e_pt[1].clear();
    gen_e_eta[0].clear();
    gen_e_eta[1].clear();
    gen_e_phi[0].clear();
    gen_e_phi[1].clear();

    gen_isEta2EE.clear();
    gen_isEta2EEGamma.clear();
    gen_isEtaPrime2EE.clear();
    gen_isEtaPrime2EEGamma.clear();
    gen_isOmega2EE.clear();
    gen_isOmega2Pi0EE.clear();
    gen_isRho2EE.clear();
    gen_isPhi2EE.clear();
    gen_isPhi2KK.clear();
    
    gen_matchedDaughtersIDs.clear();
    gen_matchedPhotonPt.clear();
    gen_matchedPhotonEta.clear();
    gen_matchedPhotonPhi.clear();

    tree->Fill();
  } else if (fillGenTree) {
    ee_mass = -1;
    ee_pt = -1;
    ee_deltaR = -1;
    for (int ii = 0; ii < 2; ii++) {
      e_pt[ii] = -1;
      e_eta[ii] = -1;
      e_phi[ii] = -1;
      e_pfIso[ii] = -1;
      e_dxy[ii] = -1;
      e_dz[ii] = -1;
      e_trkChi2[ii] = -1;
      e_trkNdof[ii] = -1;
    }
    /*vtx_prob = -1;
    vtx_pos[0] = -1; vtx_pos[1] = -1; vtx_pos[2] = -1;
    vtx_posError[0] = -1; vtx_posError[1] = -1; vtx_posError[2] = -1;
    vtx_chi2 = -1;
    npv = -1;
    pv_pos[0] = -1; pv_pos[1] = -1; pv_pos[2] = -1;*/

    e_ID[0].clear();
    e_ID[1].clear();
    l1Result_.clear();
    hltResult_.clear();

    slimmedPhoton_pt.clear();
    slimmedPhoton_eta.clear();
    slimmedPhoton_phi.clear();
    slimmedPhoton_mass.clear();
    slimmedPhoton_sigmaIetaIeta.clear();
    slimmedPhoton_hOverE.clear();
    slimmedPhoton_ecalIso.clear();
    slimmedPhoton_hcalIso.clear();
    slimmedPhoton_trkIso.clear();
    slimmedPhoton_r9.clear();

    pfCandPhoton_deltaR.clear();
    pfCandPhoton_iso.clear();
    pfCandPhoton_pt.clear();
    pfCandPhoton_eta.clear();
    pfCandPhoton_phi.clear();
    pfCandPhoton_energy.clear();
    pfCandPhoton_et.clear();
    pfCandPhoton_et2.clear();

    tree->Fill();
  }
}

// ------------ method called once each job just before starting event loop  ------------
void EEGammaTreeMaker::beginJob() {
    edm::Service<TFileService> fs;
    tree = fs->make<TTree>("tree_0"      , "tree_0");

    tree->Branch("eventNum"            , &eventNum                    , "eventNum/I");
    tree->Branch("lumiSec"             , &lumiSec                     , "lumiSec/I");
    tree->Branch("runNum"              , &runNum                      , "runNum/I");

    tree->Branch("ee_mass"                , &ee_mass                        , "ee_mass/F"    );
    tree->Branch("ee_pt"                  , &ee_pt                          , "ee_pt/F"      );
    tree->Branch("ee_deltaR"              , &ee_deltaR                      , "ee_deltaR/F"      );

    tree->Branch("e1_pt"                 , &e_pt[0]                         , "e1_pt/F"     );
    tree->Branch("e1_eta"                , &e_eta[0]                        , "e1_eta/F"    );
    tree->Branch("e1_phi"                , &e_phi[0]                        , "e1_phi/F"    );
    tree->Branch("e1_pfIso"              , &e_pfIso[0]                      , "e1_pfIso/F"  );

    tree->Branch("e2_pt"                 , &e_pt[1]                         , "e2_pt/F"     );
    tree->Branch("e2_eta"                , &e_eta[1]                        , "e2_eta/F"    );
    tree->Branch("e2_phi"                , &e_phi[1]                        , "e2_phi/F"    );
    tree->Branch("e2_pfIso"              , &e_pfIso[1]                      , "e2_pfIso/F"  );

    tree->Branch("e1_dxy"              , &e_dxy[0]                      , "e1_dxy/F"  );
    tree->Branch("e1_dz"               , &e_dz[0]                       , "e1_dz/F"  );
    tree->Branch("e1_trkChi2"          , &e_trkChi2[0]                  , "e1_trkChi2/F"  );
    tree->Branch("e1_trkNdof"          , &e_trkNdof[0]                  , "e1_trkNdof/F"  );

    tree->Branch("e2_dxy"              , &e_dxy[1]                      , "e2_dxy/F"  );
    tree->Branch("e2_dz"               , &e_dz[1]                       , "e2_dz/F"  );
    tree->Branch("e2_trkChi2"          , &e_trkChi2[1]                  , "e2_trkChi2/F"  );
    tree->Branch("e2_trkNdof"          , &e_trkNdof[1]                  , "e2_trkNdof/F"  );

    //tree->Branch("rho"                 , &rho                         , "rho/F"     );
    /*
    tree->Branch("probVtx"            , &vtx_prob                      , "probVtx/F"  );
    tree->Branch("vtxX"               , &vtx_pos[0]                         , "vtxX/F"  );
    tree->Branch("vtxY"               , &vtx_pos[1]                          , "vtxY/F"  );
    tree->Branch("vtxZ"               , &vtx_pos[2]                          , "vtxZ/F"  );
    tree->Branch("vtxXError"          , &vtx_posError[0]                     , "vtxXError/F"  );
    tree->Branch("vtxYError"          , &vtx_posError[1]                    , "vtxYError/F"  );
    tree->Branch("vtxZError"          , &vtx_posError[2]                    , "vtxZError/F"  );
    tree->Branch("vtx_chi2"           , &vtx_chi2                     , "vtx_chi2/F"  );

    tree->Branch("npv"                , &npv                          , "npv/I"  );
    tree->Branch("pvX"                , &pv_pos[0]                          , "pvX/F"  );
    tree->Branch("pvY"                , &pv_pos[1]                          , "pvY/F"  );
    tree->Branch("pvZ"                , &pv_pos[2]                          , "pvZ/F"  );
    */
    tree->Branch("e1_ID", "std::vector<bool>", &e_ID[0], 32000, 0);
    tree->Branch("e2_ID", "std::vector<bool>", &e_ID[1], 32000, 0);

    tree->Branch("l1Result", "std::vector<bool>"             ,&l1Result_, 32000, 0  );
    tree->Branch("hltResult", "std::vector<bool>"             ,&hltResult_, 32000, 0  );

    tree->Branch("slimmedPhoton_pt", "std::vector<float>", &slimmedPhoton_pt, 32000, 0);
    tree->Branch("slimmedPhoton_eta", "std::vector<float>", &slimmedPhoton_eta, 32000, 0);
    tree->Branch("slimmedPhoton_phi", "std::vector<float>", &slimmedPhoton_phi, 32000, 0);
    tree->Branch("slimmedPhoton_mass", "std::vector<float>", &slimmedPhoton_mass, 32000, 0);
    tree->Branch("slimmedPhoton_sigmaIetaIeta", "std::vector<float>", &slimmedPhoton_sigmaIetaIeta, 32000, 0);
    tree->Branch("slimmedPhoton_hOverE", "std::vector<float>", &slimmedPhoton_hOverE, 32000, 0);
    tree->Branch("slimmedPhoton_ecalIso", "std::vector<float>", &slimmedPhoton_ecalIso, 32000, 0);
    tree->Branch("slimmedPhoton_hcalIso", "std::vector<float>", &slimmedPhoton_hcalIso, 32000, 0);
    tree->Branch("slimmedPhoton_trkIso", "std::vector<float>", &slimmedPhoton_trkIso, 32000, 0);
    tree->Branch("slimmedPhoton_r9", "std::vector<float>", &slimmedPhoton_r9, 32000, 0);

    tree->Branch("pfCandPhoton_deltaR", "std::vector<float>", &pfCandPhoton_deltaR, 32000, 0);
    tree->Branch("pfCandPhoton_iso", "std::vector<float>", &pfCandPhoton_iso, 32000, 0);
    tree->Branch("pfCandPhoton_pt", "std::vector<float>", &pfCandPhoton_pt, 32000, 0);
    tree->Branch("pfCandPhoton_eta", "std::vector<float>", &pfCandPhoton_eta, 32000, 0);
    tree->Branch("pfCandPhoton_phi", "std::vector<float>", &pfCandPhoton_phi, 32000, 0);
    tree->Branch("pfCandPhoton_energy", "std::vector<float>", &pfCandPhoton_energy, 32000, 0);
    tree->Branch("pfCandPhoton_et", "std::vector<float>", &pfCandPhoton_et, 32000, 0);
    tree->Branch("pfCandPhoton_et2", "std::vector<float>", &pfCandPhoton_et2, 32000, 0);


    tree->Branch("matchedGenIdx"              , &matchedGenIdx                    , "matchedGenIdx/I"  );

    tree->Branch("gen_motherID", "std::vector<int>", &gen_motherID, 32000, 0);
    //tree->Branch("gen_grandMotherID", "std::vector<int>", &gen_grandMotherID, 32000, 0);
    tree->Branch("gen_ee_mass", "std::vector<float>", &gen_ee_mass, 32000, 0);

    tree->Branch("gen_e1_pt", "std::vector<float>", &gen_e_pt[0], 32000, 0);
    tree->Branch("gen_e1_eta", "std::vector<float>", &gen_e_eta[0], 32000, 0);
    tree->Branch("gen_e1_phi", "std::vector<float>", &gen_e_phi[0], 32000, 0);
    tree->Branch("gen_e1_recoMatch", "std::vector<bool>", &gen_e_recoMatch[0], 32000, 0);

    tree->Branch("gen_e2_pt", "std::vector<float>", &gen_e_pt[1], 32000, 0);
    tree->Branch("gen_e2_eta", "std::vector<float>", &gen_e_eta[1], 32000, 0);
    tree->Branch("gen_e2_phi", "std::vector<float>", &gen_e_phi[1], 32000, 0);
    tree->Branch("gen_e2_recoMatch", "std::vector<bool>", &gen_e_recoMatch[1], 32000, 0);

    tree->Branch("gen_matchedDaughtersIDs", "std::vector<std::vector<int> >", &gen_matchedDaughtersIDs, 32000, 0);
    tree->Branch("gen_matchedPhotonPt", "std::vector<std::vector<float> >", &gen_matchedPhotonPt, 32000, 0);
    tree->Branch("gen_matchedPhotonEta", "std::vector<std::vector<float> >", &gen_matchedPhotonEta, 32000, 0);
    tree->Branch("gen_matchedPhotonPhi", "std::vector<std::vector<float> >", &gen_matchedPhotonPhi, 32000, 0);

    tree->Branch("gen_isEta2EE", "std::vector<bool>", &gen_isEta2EE, 32000, 0);
    tree->Branch("gen_isEta2EEGamma", "std::vector<bool>", &gen_isEta2EEGamma, 32000, 0);
    tree->Branch("gen_isEtaPrime2EE", "std::vector<bool>", &gen_isEtaPrime2EE, 32000, 0);
    tree->Branch("gen_isEtaPrime2EEGamma", "std::vector<bool>", &gen_isEtaPrime2EEGamma, 32000, 0);
    tree->Branch("gen_isOmega2EE", "std::vector<bool>", &gen_isOmega2EE, 32000, 0);
    tree->Branch("gen_isOmega2Pi0EE", "std::vector<bool>", &gen_isOmega2Pi0EE, 32000, 0);
    tree->Branch("gen_isRho2EE", "std::vector<bool>", &gen_isRho2EE, 32000, 0);
    tree->Branch("gen_isPhi2EE", "std::vector<bool>", &gen_isPhi2EE, 32000, 0);
    tree->Branch("gen_isPhi2KK", "std::vector<bool>", &gen_isPhi2KK, 32000, 0);
}

// ------------ method called once each job just after ending the event loop  ------------
void EEGammaTreeMaker::endJob() {
  // please remove this method if not needed
}

void EEGammaTreeMaker::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
    // HLT paths
    triggerPathsVector.push_back("HLT_DoubleEle*_eta1p22_mMax6_v*");
    triggerPathsVector.push_back("HLT_DoubleEle*_eta1p22_mMax6_trkHits10_v*");
    triggerPathsVector.push_back("HLT_DoubleMu4_3_LowMass_v*");
    triggerPathsVector.push_back("HLT_DoubleMu4_LowMass_Displaced_v*");
    HLTConfigProvider hltConfig;
    bool changedConfig = false;
    hltConfig.init(iRun, iSetup, triggerResultsTag.process(), changedConfig);

    for (size_t i = 0; i < triggerPathsVector.size(); i++) {
        triggerPathsMap[triggerPathsVector[i]] = -1;
    }

    for(size_t i = 0; i < triggerPathsVector.size(); i++){
        TPRegexp pattern(triggerPathsVector[i]);
        for(size_t j = 0; j < hltConfig.triggerNames().size(); j++){
            std::string pathName = hltConfig.triggerNames()[j];
            if(TString(pathName).Contains(pattern)){
                triggerPathsMap[triggerPathsVector[i]] = j;
            }
        }
    }
}

void EEGammaTreeMaker::endRun(edm::Run const&, edm::EventSetup const&) {
}

void EEGammaTreeMaker::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) {
}

void EEGammaTreeMaker::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void EEGammaTreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EEGammaTreeMaker);
