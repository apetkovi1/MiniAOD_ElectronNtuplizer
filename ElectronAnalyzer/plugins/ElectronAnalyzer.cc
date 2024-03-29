#include <memory>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include <vector>
#include <algorithm>
#include <map>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include <FWCore/ParameterSet/interface/FileInPath.h>

#include "/interface/ElectronEffArea.h"

float Ele_Eff_Area(float SCeta_stuff); //  Added the stuff needed for RelIso calculations

float zero = 0.0;

class ElectronAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit ElectronAnalyzer(const edm::ParameterSet&);
      ~ElectronAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<pat::ElectronCollection> elecCollToken;
      edm::InputTag elecSrc_;
      edm::EDGetTokenT<double> theRhoToken;
      //edm::EDGetTokenT<float> theRhoToken;
      edm::InputTag rhoSrc_;
      typedef std::vector<PileupSummaryInfo> PileupSummaryInfoCollection;
      edm::EDGetTokenT<PileupSummaryInfoCollection> pileupToken;
      edm::InputTag pileupSrc_;

      TTree *electron_tree;
      std::vector<float> ele_pt,ele_eta,ele_phi,scl_eta,ele_oldsigmaietaieta,ele_oldsigmaiphiiphi,ele_oldcircularity,ele_oldr9,ele_scletawidth,ele_sclphiwidth,ele_he,ele_oldhe,
      ele_kfchi2,ele_gsfchi2,ele_fbrem,ele_ep,ele_eelepout,ele_IoEmIop,ele_deltaetain,ele_deltaphiin,ele_deltaetaseed,
      ele_psEoverEraw,ele_pfPhotonIso,ele_pfChargedHadIso,ele_pfNeutralHadIso,ele_PFPUIso,ElectronMVAEstimatorRun2Fall17IsoV2Values,ElectronMVAEstimatorRun2Fall17IsoV1Values,
      ElectronMVAEstimatorRun2Fall17NoIsoV1Values,ElectronMVAEstimatorRun2Fall17NoIsoV2Values, gen_pt, gen_eta, gen_phi, mother_pt, mother_eta, mother_phi, relISO_a, relISO_a_corr;
      std::vector<int> ele_kfhits, ele_chi2_hits,ele_gsfhits, ele_expected_inner_hits,ele_charge,ele_mother;
      float Diele_mass, rho, dR, relISO_0, relISO_1,relISO_0_corr, relISO_1_corr;
      int numele, PFnumele;
      TLorentzVector P,P0,P1;
      std::vector<bool> ele_isPF, cutBasedElectronID_Fall17_94X_V2_veto, cutBasedElectronID_Fall17_94X_V2_loose, cutBasedElectronID_Fall17_94X_V2_medium, cutBasedElectronID_Fall17_94X_V2_tight,
      mvaEleID_Fall17_iso_V2_wp90, mvaEleID_Fall17_iso_V2_wp80, mvaEleID_Fall17_noIso_V2_wp90, mvaEleID_Fall17_noIso_V2_wp80,mvaEleID_Fall17_noIso_V2_wpLoose_unsopported,
      mvaEleID_Fall17_iso_V2_wpHZZ_unsopported;

};

ElectronAnalyzer::ElectronAnalyzer(const edm::ParameterSet& iConfig):
elecSrc_(iConfig.getUntrackedParameter<edm::InputTag>("elecSrc")),
rhoSrc_(iConfig.getUntrackedParameter<edm::InputTag>("rhoSrc")),
pileupSrc_(iConfig.getUntrackedParameter<edm::InputTag>("pileupSrc"))

{
   elecCollToken = consumes<pat::ElectronCollection>(elecSrc_);
   theRhoToken = consumes<double>(rhoSrc_);
   pileupToken = consumes<PileupSummaryInfoCollection>(pileupSrc_);

   edm::Service<TFileService> fs;
   electron_tree = fs->make<TTree>("Events", "Events");

   electron_tree->Branch("numele",&numele);
   electron_tree->Branch("PFnumele",&PFnumele);
   electron_tree->Branch("ele_eta",&ele_eta);
   electron_tree->Branch("ele_phi",&ele_phi);
   electron_tree->Branch("ele_isPF",&ele_isPF);
   electron_tree->Branch("Diele_mass",&Diele_mass);

   //BDT ID branches
   electron_tree->Branch("scl_eta",&scl_eta);
   electron_tree->Branch("ele_pt",&ele_pt);
   electron_tree->Branch("ele_oldsigmaietaieta",&ele_oldsigmaietaieta);
   electron_tree->Branch("ele_oldsigmaiphiiphi",&ele_oldsigmaiphiiphi);
   electron_tree->Branch("ele_oldcircularity",&ele_oldcircularity);
   electron_tree->Branch("ele_oldr9",&ele_oldr9 );
   electron_tree->Branch("ele_scletawidth",&ele_scletawidth);
   electron_tree->Branch("ele_sclphiwidth",&ele_sclphiwidth);
   electron_tree->Branch("ele_he",&ele_he);
   electron_tree->Branch("ele_oldhe",&ele_oldhe);
   electron_tree->Branch("ele_kfhits",&ele_kfhits);
   electron_tree->Branch("ele_kfchi2",&ele_kfchi2 );
   electron_tree->Branch("ele_gsfchi2",&ele_gsfchi2);
   electron_tree->Branch("ele_chi2_hits",&ele_chi2_hits);
   electron_tree->Branch("ele_fbrem",&ele_fbrem);
   electron_tree->Branch("ele_gsfhits",&ele_gsfhits);
   electron_tree->Branch("ele_expected_inner_hits",&ele_expected_inner_hits);
   electron_tree->Branch("ele_ep",&ele_ep);
   electron_tree->Branch("ele_eelepout",&ele_eelepout);
   electron_tree->Branch("ele_IoEmIop",&ele_IoEmIop);
   electron_tree->Branch("ele_deltaetain",&ele_deltaetain);
   electron_tree->Branch("ele_deltaphiin",&ele_deltaphiin);
   electron_tree->Branch("ele_deltaetaseed",&ele_deltaetaseed);
   electron_tree->Branch("ele_psEoverEraw",&ele_psEoverEraw);
   electron_tree->Branch("rho",&rho);
   electron_tree->Branch("ele_charge",&ele_charge);

   //isolation variables
   electron_tree->Branch("ele_pfPhotonIso",&ele_pfPhotonIso);
   electron_tree->Branch("ele_pfChargedHadIso",&ele_pfChargedHadIso);
   electron_tree->Branch("ele_pfNeutralHadIso",&ele_pfNeutralHadIso);
   electron_tree->Branch("ele_PFPUIso",&ele_PFPUIso);
   electron_tree->Branch("relISO_0",&relISO_0);
   electron_tree->Branch("relISO_1",&relISO_1);
   electron_tree->Branch("relISO_a",&relISO_a);
   electron_tree->Branch("relISO_0_corr",&relISO_0_corr);
   electron_tree->Branch("relISO_1_corr",&relISO_1_corr);
   electron_tree->Branch("relISO_a_corr",&relISO_a_corr);

   //Already implemented ID
   electron_tree->Branch("cutBasedElectronID_Fall17_94X_V2_veto",&cutBasedElectronID_Fall17_94X_V2_veto);
   electron_tree->Branch("cutBasedElectronID_Fall17_94X_V2_loose",&cutBasedElectronID_Fall17_94X_V2_loose);
   electron_tree->Branch("cutBasedElectronID_Fall17_94X_V2_medium",&cutBasedElectronID_Fall17_94X_V2_medium);
   electron_tree->Branch("cutBasedElectronID_Fall17_94X_V2_tight",&cutBasedElectronID_Fall17_94X_V2_tight);
   electron_tree->Branch("mvaEleID_Fall17_iso_V2_wp90",&mvaEleID_Fall17_iso_V2_wp90);
   electron_tree->Branch("mvaEleID_Fall17_iso_V2_wp80",&mvaEleID_Fall17_iso_V2_wp80);
   electron_tree->Branch("mvaEleID_Fall17_noIso_V2_wp90",&mvaEleID_Fall17_noIso_V2_wp90);
   electron_tree->Branch("mvaEleID_Fall17_noIso_V2_wp80",&mvaEleID_Fall17_noIso_V2_wp80);

   //V2 and V1 MVA scores
   electron_tree->Branch("ElectronMVAEstimatorRun2Fall17IsoV2Values",&ElectronMVAEstimatorRun2Fall17IsoV2Values);
   electron_tree->Branch("ElectronMVAEstimatorRun2Fall17IsoV1Values",&ElectronMVAEstimatorRun2Fall17IsoV1Values);
   electron_tree->Branch("ElectronMVAEstimatorRun2Fall17NoIsoV2Values",&ElectronMVAEstimatorRun2Fall17NoIsoV2Values);
   electron_tree->Branch("ElectronMVAEstimatorRun2Fall17NoIsoV1Values",&ElectronMVAEstimatorRun2Fall17NoIsoV1Values);

   //unsopported IDs
   electron_tree->Branch("mvaEleID_Fall17_noIso_V2_wpLoose_unsopported",&mvaEleID_Fall17_noIso_V2_wpLoose_unsopported);
   electron_tree->Branch("mvaEleID_Fall17_iso_V2_wpHZZ_unsopported",&mvaEleID_Fall17_iso_V2_wpHZZ_unsopported);

   //Gen Info
   electron_tree->Branch("gen_pt",&gen_pt);
   electron_tree->Branch("gen_eta",&gen_eta);
   electron_tree->Branch("gen_phi",&gen_phi);
   electron_tree->Branch("mother_pt",&mother_pt);
   electron_tree->Branch("mother_eta",&mother_eta);
   electron_tree->Branch("mother_phi",&mother_phi);
   electron_tree->Branch("ele_mother",&ele_mother);

   //deltaR
   electron_tree->Branch("dR",&dR);
}


ElectronAnalyzer::~ElectronAnalyzer()
{


}

void
ElectronAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   edm::Handle< std::vector<pat::Electron>> electrons;
   iEvent.getByToken(elecCollToken, electrons);

   edm::Handle<double> rhoHandle;
   iEvent.getByToken(theRhoToken, rhoHandle);

   edm::Handle<PileupSummaryInfoCollection> pileupSummaryInfos;
   iEvent.getByToken(pileupToken, pileupSummaryInfos);

   numele=0;
   PFnumele=0;
   Diele_mass=0;
   rho=0;
   dR=0;
   relISO_0=0;
   relISO_1=0;
   relISO_0_corr=0;
   relISO_1_corr=0;
   relISO_a.clear();
   relISO_a_corr.clear();
   ele_eta.clear();
   ele_phi.clear();
   ele_isPF.clear();
   scl_eta.clear();
   ele_pt.clear();
   ele_mother.clear();
   ele_oldsigmaietaieta.clear();
   ele_oldsigmaiphiiphi.clear();
   ele_oldcircularity.clear();
   ele_oldr9.clear();
   ele_scletawidth.clear();
   ele_sclphiwidth.clear();
   ele_he.clear();
   ele_oldhe.clear();
   ele_kfhits.clear();
   ele_kfchi2.clear();
   ele_gsfchi2.clear();
   ele_chi2_hits.clear();
   ele_fbrem.clear();
   ele_gsfhits.clear();
   ele_expected_inner_hits.clear();
   ele_ep.clear();
   ele_eelepout.clear();
   ele_IoEmIop.clear();
   ele_deltaetain.clear();
   ele_deltaphiin.clear();
   ele_deltaetaseed.clear();
   ele_psEoverEraw.clear();
   ele_pfPhotonIso.clear();
   ele_pfChargedHadIso.clear();
   ele_pfNeutralHadIso.clear();
   ele_PFPUIso.clear();
   ele_charge.clear();
   cutBasedElectronID_Fall17_94X_V2_veto.clear();
   cutBasedElectronID_Fall17_94X_V2_loose.clear();
   cutBasedElectronID_Fall17_94X_V2_medium.clear();
   cutBasedElectronID_Fall17_94X_V2_tight.clear();
   mvaEleID_Fall17_iso_V2_wp80.clear();
   mvaEleID_Fall17_iso_V2_wp90.clear();
   mvaEleID_Fall17_noIso_V2_wp80.clear();
   mvaEleID_Fall17_noIso_V2_wp90.clear();
   ElectronMVAEstimatorRun2Fall17IsoV1Values.clear();
   ElectronMVAEstimatorRun2Fall17IsoV2Values.clear();
   ElectronMVAEstimatorRun2Fall17NoIsoV1Values.clear();
   ElectronMVAEstimatorRun2Fall17NoIsoV2Values.clear();
   mvaEleID_Fall17_noIso_V2_wpLoose_unsopported.clear();
   mvaEleID_Fall17_iso_V2_wpHZZ_unsopported.clear();
   gen_eta.clear();
   gen_phi.clear();
   gen_pt.clear();
   mother_eta.clear();
   mother_phi.clear();
   mother_pt.clear();

   /*for (PileupSummaryInfoCollection::const_iterator pileupSummaryInfo = pileupSummaryInfos->begin(); pileupSummaryInfo != pileupSummaryInfos->end(); ++pileupSummaryInfo)
    {
      if ( pileupSummaryInfo->getBunchCrossing() == 0 )
      {
	      std::cout<<pileupSummaryInfo->getTrueNumInteractions()<<std::endl;  //getTrueNumInteractions, getPU_NumInteractions
      }
    }*/

   for (auto it = electrons->cbegin(); it != electrons->cend(); ++it)
   {
     numele++;
     if(it->isPF())
     {
     PFnumele++;
     ele_eta.push_back(it->eta());
     ele_phi.push_back(it->phi());
     ele_isPF.push_back(it->isPF());
     scl_eta.push_back(it->superCluster()->eta());
     ele_pt.push_back(it->pt());
     ele_oldsigmaietaieta.push_back(it->full5x5_sigmaIetaIeta());
     ele_oldsigmaiphiiphi.push_back(it->full5x5_sigmaIphiIphi());
     ele_oldcircularity.push_back(1.0-(it->full5x5_e1x5())/(it->full5x5_e5x5()));
     ele_oldr9.push_back(it->full5x5_r9());
     ele_scletawidth.push_back(it->superCluster()->etaWidth());
     ele_sclphiwidth.push_back(it->superCluster()->phiWidth());
     ele_he.push_back(it->hadronicOverEm());
     ele_oldhe.push_back(it->full5x5_hcalOverEcal());
     ele_kfhits.push_back((it->closestCtfTrackRef().isAvailable() && it->closestCtfTrackRef().isNonnull()) ? it->closestCtfTrackRef()->hitPattern().trackerLayersWithMeasurement() : -1);
     ele_kfchi2.push_back((it->closestCtfTrackRef().isAvailable() && it->closestCtfTrackRef().isNonnull()) ? it->closestCtfTrackRef()->normalizedChi2() : 0);
     ele_gsfchi2.push_back(it->gsfTrack()->normalizedChi2());
     ele_chi2_hits.push_back(it->gsfTrack()->normalizedChi2());
     ele_fbrem.push_back(it->fbrem());
     ele_gsfhits.push_back(it->gsfTrack()->hitPattern().trackerLayersWithMeasurement());
     ele_expected_inner_hits.push_back(it->gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS));
     //ele_conversionVertexFitProbability.push_back(it->convVtxFitProb()) ;
     ele_ep.push_back(it->eSuperClusterOverP());
     ele_eelepout.push_back(it->eEleClusterOverPout());
     ele_IoEmIop.push_back(1.0/(it->ecalEnergy())-1.0/(it->trackMomentumAtVtx().R()));
     ele_deltaetain.push_back(it->deltaEtaSuperClusterTrackAtVtx());
     ele_deltaphiin.push_back(it->deltaPhiSuperClusterTrackAtVtx());
     ele_deltaetaseed.push_back(it->deltaEtaSeedClusterTrackAtCalo());
     ele_psEoverEraw.push_back((it->superCluster()->preshowerEnergy())/(it->superCluster()->rawEnergy()));
     ele_pfPhotonIso.push_back(it->pfIsolationVariables().sumPhotonEt);
     ele_pfChargedHadIso.push_back(it->pfIsolationVariables().sumChargedHadronPt);
     ele_pfNeutralHadIso.push_back(it->pfIsolationVariables().sumNeutralHadronEt);
     ele_PFPUIso.push_back(it->pfIsolationVariables().sumPUPt);
     ele_charge.push_back(it->charge());
     cutBasedElectronID_Fall17_94X_V2_veto.push_back(it->electronID("cutBasedElectronID-Fall17-94X-V2-veto"));
     cutBasedElectronID_Fall17_94X_V2_loose.push_back(it->electronID("cutBasedElectronID-Fall17-94X-V2-loose"));
     cutBasedElectronID_Fall17_94X_V2_medium.push_back(it->electronID("cutBasedElectronID-Fall17-94X-V2-medium"));
     cutBasedElectronID_Fall17_94X_V2_tight.push_back(it->electronID("cutBasedElectronID-Fall17-94X-V2-tight"));
     mvaEleID_Fall17_iso_V2_wp90.push_back(it->electronID("mvaEleID-Fall17-iso-V2-wp90"));
     mvaEleID_Fall17_iso_V2_wp80.push_back(it->electronID("mvaEleID-Fall17-iso-V2-wp80"));
     mvaEleID_Fall17_noIso_V2_wp90.push_back(it->electronID("mvaEleID-Fall17-noIso-V2-wp90"));
     mvaEleID_Fall17_noIso_V2_wp80.push_back(it->electronID("mvaEleID-Fall17-noIso-V2-wp80"));
     ElectronMVAEstimatorRun2Fall17IsoV1Values.push_back(it->userFloat("ElectronMVAEstimatorRun2Fall17IsoV1Values"));
     ElectronMVAEstimatorRun2Fall17IsoV2Values.push_back(it->userFloat("ElectronMVAEstimatorRun2Fall17IsoV2Values"));
     ElectronMVAEstimatorRun2Fall17NoIsoV1Values.push_back(it->userFloat("ElectronMVAEstimatorRun2Fall17NoIsoV1Values"));
     ElectronMVAEstimatorRun2Fall17NoIsoV2Values.push_back(it->userFloat("ElectronMVAEstimatorRun2Fall17NoIsoV2Values"));
     mvaEleID_Fall17_noIso_V2_wpLoose_unsopported.push_back(it->electronID("mvaEleID-Fall17-noIso-V2-wpLoose"));
     mvaEleID_Fall17_iso_V2_wpHZZ_unsopported.push_back(it->electronID("mvaEleID-Fall17-iso-V2-wpHZZ"));

     if ((it->genParticleRef ()).isNonnull ())
     {
       ele_mother.push_back(it->genParticle()->mother(0)->pdgId());
       gen_pt.push_back(it->genParticle()->pt());
       gen_eta.push_back(it->genParticle()->eta());
       gen_phi.push_back(it->genParticle()->phi());
       mother_pt.push_back(it->genParticle()->mother(0)->pt());
       mother_eta.push_back(it->genParticle()->mother(0)->eta());
       mother_phi.push_back(it->genParticle()->mother(0)->phi());
     }
     }
   }

   if(PFnumele==2 && ele_charge[0]==(-1.0)*ele_charge[1])
   {
     P0.SetPtEtaPhiM(ele_pt[0],ele_eta[0],ele_phi[0],0.00051);
     P1.SetPtEtaPhiM(ele_pt[1],ele_eta[1],ele_phi[1],0.00051);
     P=P0+P1;
     Diele_mass=P.M();
     rho=*(rhoHandle.product());
     dR=deltaR(ele_eta[0],ele_phi[0],ele_eta[1],ele_phi[1]);
     relISO_0 = (ele_pfPhotonIso[0]+ele_pfChargedHadIso[0]+ele_pfNeutralHadIso[0])/ele_pt[0];
     relISO_1 = (ele_pfPhotonIso[1]+ele_pfChargedHadIso[1]+ele_pfNeutralHadIso[1])/ele_pt[1];
     relISO_a = {relISO_0,relISO_1};
     
     relISO_0_corr = (ele_pfChargedHadIso[0]+std::max(zero,ele_pfPhotonIso[0]+ele_pfNeutralHadIso[0]-rho*Ele_Eff_Area(scl_eta[0])))/ele_pt[0];
     relISO_1_corr = (ele_pfChargedHadIso[1]+std::max(zero,ele_pfPhotonIso[1]+ele_pfNeutralHadIso[1]-rho*Ele_Eff_Area(scl_eta[1])))/ele_pt[1];
     relISO_a_corr = {relISO_0_corr,relISO_1_corr};

     
     
     /*
     if(dR < 0.3){

      relISO_0_woEle = (ele_pfPhotonIso[0]+ele_pfChargedHadIso[0]+ele_pfNeutralHadIso[0]-ele_pt[1])/ele_pt[0];
      relISO_1_woEle = (ele_pfPhotonIso[1]+ele_pfChargedHadIso[1]+ele_pfNeutralHadIso[1]-ele_pt[0])/ele_pt[1];
      relISO_a_woEle = {relISO_0_woEle,relISO_1_woEle};

     }
     else{

      relISO_0_woEle = (ele_pfPhotonIso[0]+ele_pfChargedHadIso[0]+ele_pfNeutralHadIso[0])/ele_pt[0];
      relISO_1_woEle = (ele_pfPhotonIso[1]+ele_pfChargedHadIso[1]+ele_pfNeutralHadIso[1])/ele_pt[1];
      relISO_a_woEle = {relISO_0_woEle,relISO_1_woEle};

     }*/


     if(Diele_mass>1 && Diele_mass<5)
     electron_tree->Fill();
   }

}

void
ElectronAnalyzer::beginJob()
{
}

void
ElectronAnalyzer::endJob()
{
}

void
ElectronAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

DEFINE_FWK_MODULE(ElectronAnalyzer);
