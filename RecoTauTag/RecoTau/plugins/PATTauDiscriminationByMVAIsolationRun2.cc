
/*
 * \class PATTauDiscriminationByMVAIsolationRun2
 * 
 * MVA based discriminator against jet -> tau fakes
 * 
 * Adopted from RecoTauTag/RecoTau/plugins/PFRecoTauDiscriminationByMVAIsolationRun2.cc
 * to enable computation of MVA isolation on MiniAOD
 * 
 * \author Alexander Nehrkorn, RWTH Aachen
 */

// todo 1: remove leadingTrackChi2 as input variable from:
//           - here
//           - TauPFEssential
//           - PFRecoTauDiscriminationByMVAIsolationRun2
//           - Training of BDT

#include "RecoTauTag/RecoTau/interface/TauDiscriminationProducerBase.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/PATTauDiscriminator.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "RecoTauTag/RecoTau/interface/PFRecoTauClusterVariables.h"

#include "CondFormats/EgammaObjects/interface/GBRForest.h"
#include "CondFormats/DataRecord/interface/GBRWrapperRcd.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include <TFile.h>

#include <iostream>

using namespace pat;

namespace
{
  const GBRForest* loadMVAfromFile(const edm::FileInPath& inputFileName, const std::string& mvaName, std::vector<TFile*>& inputFilesToDelete)
  {
    if ( inputFileName.location() == edm::FileInPath::Unknown ) throw cms::Exception("PATTauDiscriminationByIsolationMVARun2::loadMVA")
      << " Failed to find File = " << inputFileName << " !!\n"; 
    TFile* inputFile = new TFile(inputFileName.fullPath().data());
	
    //const GBRForest* mva = dynamic_cast<GBRForest*>(inputFile->Get(mvaName.data())); // CV: dynamic_cast<GBRForest*> fails for some reason ?!
    const GBRForest* mva = (GBRForest*)inputFile->Get(mvaName.data());
    if ( !mva )
      throw cms::Exception("PATTauDiscriminationByIsolationMVARun2::loadMVA")
        << " Failed to load MVA = " << mvaName.data() << " from file = " << inputFileName.fullPath().data() << " !!\n";
	
     inputFilesToDelete.push_back(inputFile);
	  
     return mva;
  }
	
  const GBRForest* loadMVAfromDB(const edm::EventSetup& es, const std::string& mvaName)
  {
    edm::ESHandle<GBRForest> mva;
    es.get<GBRWrapperRcd>().get(mvaName, mva);
    return mva.product();
  }
}

class PATTauDiscriminationByMVAIsolationRun2 : public PATTauDiscriminationProducerBase
{
  public:
    explicit PATTauDiscriminationByMVAIsolationRun2(const edm::ParameterSet& cfg)
      : PATTauDiscriminationProducerBase(cfg),
        moduleLabel_(cfg.getParameter<std::string>("@module_label")),
	mvaReader_(nullptr),
	mvaInput_(nullptr),
	category_output_()
    {
       mvaName_ = cfg.getParameter<std::string>("mvaName");
       loadMVAfromDB_ = cfg.getParameter<bool>("loadMVAfromDB");
       if ( !loadMVAfromDB_ ) {
	   inputFileName_ = cfg.getParameter<edm::FileInPath>("inputFileName");
       }
       std::string mvaOpt_string = cfg.getParameter<std::string>("mvaOpt");
       if      ( mvaOpt_string == "oldDMwoLT" ) mvaOpt_ = kOldDMwoLT;
       else if ( mvaOpt_string == "oldDMwLT"  ) mvaOpt_ = kOldDMwLT;
       else if ( mvaOpt_string == "newDMwoLT" ) mvaOpt_ = kNewDMwoLT;
       else if ( mvaOpt_string == "newDMwLT"  ) mvaOpt_ = kNewDMwLT;
       else if ( mvaOpt_string == "DBoldDMwLT"  ) mvaOpt_ = kDBoldDMwLT;
       else if ( mvaOpt_string == "DBnewDMwLT"  ) mvaOpt_ = kDBnewDMwLT;
       else if ( mvaOpt_string == "PWoldDMwLT"  ) mvaOpt_ = kPWoldDMwLT;
       else if ( mvaOpt_string == "PWnewDMwLT"  ) mvaOpt_ = kPWnewDMwLT;
       else if ( mvaOpt_string == "DBoldDMwLTwGJ" ) mvaOpt_ = kDBoldDMwLTwGJ;
       else if ( mvaOpt_string == "DBnewDMwLTwGJ" ) mvaOpt_ = kDBnewDMwLTwGJ;
       else throw cms::Exception("PATTauDiscriminationByMVAIsolationRun2")
         << " Invalid Configuration Parameter 'mvaOpt' = " << mvaOpt_string << " !!\n";
		    
       if      ( mvaOpt_ == kOldDMwoLT || mvaOpt_ == kNewDMwoLT ) mvaInput_ = new float[6];
       else if ( mvaOpt_ == kOldDMwLT  || mvaOpt_ == kNewDMwLT  ) mvaInput_ = new float[12];
       else if ( mvaOpt_ == kDBoldDMwLT || mvaOpt_ == kDBnewDMwLT ||
                 mvaOpt_ == kPWoldDMwLT || mvaOpt_ == kPWnewDMwLT ||
                 mvaOpt_ == kDBoldDMwLTwGJ || mvaOpt_ == kDBnewDMwLTwGJ) mvaInput_ = new float[23];
       else assert(0);

       chargedIsoPtSums_ = cfg.getParameter<std::string>("srcChargedIsoPtSum");
       neutralIsoPtSums_ = cfg.getParameter<std::string>("srcNeutralIsoPtSum");
       puCorrPtSums_ = cfg.getParameter<std::string>("srcPUcorrPtSum");
       photonPtSumOutsideSignalCone_ = cfg.getParameter<std::string>("srcPhotonPtSumOutsideSignalCone");
       footprintCorrection_ = cfg.getParameter<std::string>("srcFootprintCorrection");
		  
       verbosity_ = cfg.getParameter<int>("verbosity");

       produces<pat::PATTauDiscriminator>("category");
    }  
		
    void beginEvent(const edm::Event&, const edm::EventSetup&) override;
		
    double discriminate(const TauRef&) const override;
		
    void endEvent(edm::Event&) override;
		
    ~PATTauDiscriminationByMVAIsolationRun2() override
    {
      if(!loadMVAfromDB_) delete mvaReader_;
      delete[] mvaInput_;
      for ( std::vector<TFile*>::iterator it = inputFilesToDelete_.begin();
	    it != inputFilesToDelete_.end(); ++it ) {
        delete (*it);
      }
    }
    	
  private:
		
    std::string moduleLabel_;
		
    std::string mvaName_;
    bool loadMVAfromDB_;
    edm::FileInPath inputFileName_;
    const GBRForest* mvaReader_;
    enum { kOldDMwoLT, kOldDMwLT, kNewDMwoLT, kNewDMwLT, kDBoldDMwLT, kDBnewDMwLT, kPWoldDMwLT, kPWnewDMwLT, kDBoldDMwLTwGJ, kDBnewDMwLTwGJ };
    int mvaOpt_;
    float* mvaInput_;

    std::string chargedIsoPtSums_;
    std::string neutralIsoPtSums_;
    std::string puCorrPtSums_;
    std::string photonPtSumOutsideSignalCone_;
    std::string footprintCorrection_;
		
    edm::Handle<TauCollection> taus_;
    std::unique_ptr<pat::PATTauDiscriminator> category_output_;		
    std::vector<TFile*> inputFilesToDelete_;
    TauIdMVAAuxiliaries clusterVariables_;
    	
    int verbosity_;
};

void PATTauDiscriminationByMVAIsolationRun2::beginEvent(const edm::Event& evt, const edm::EventSetup& es)
{
  if( !mvaReader_ ) {
    if ( loadMVAfromDB_ ) {
      mvaReader_ = loadMVAfromDB(es, mvaName_);
    } else {
      mvaReader_ = loadMVAfromFile(inputFileName_, mvaName_, inputFilesToDelete_);
    }
  }
	
  evt.getByToken(Tau_token, taus_);
  category_output_.reset(new pat::PATTauDiscriminator(TauRefProd(taus_)));
}

double PATTauDiscriminationByMVAIsolationRun2::discriminate(const TauRef& tau) const
{
  // CV: define dummy category index in order to use RecoTauDiscriminantCutMultiplexer module to appy WP cuts
  double category = 0.; 
  category_output_->setValue(tauIndex_, category);
	
  // CV: computation of MVA value requires presence of leading charged hadron
  if ( tau->leadChargedHadrCand().isNull() ) return 0.;
	
  int tauDecayMode = tau->decayMode();
	
  if ( ((mvaOpt_ == kOldDMwoLT || mvaOpt_ == kOldDMwLT || mvaOpt_ == kDBoldDMwLT || mvaOpt_ == kPWoldDMwLT || mvaOpt_ == kDBoldDMwLTwGJ)
        && (tauDecayMode == 0 || tauDecayMode == 1 || tauDecayMode == 2 || tauDecayMode == 10))
       ||
       ((mvaOpt_ == kNewDMwoLT || mvaOpt_ == kNewDMwLT || mvaOpt_ == kDBnewDMwLT || mvaOpt_ == kPWnewDMwLT || mvaOpt_ == kDBnewDMwLTwGJ)
        && (tauDecayMode == 0 || tauDecayMode == 1 || tauDecayMode == 2 || tauDecayMode == 5 || tauDecayMode == 6 || tauDecayMode == 10 || tauDecayMode == 11))
  ) {
	
    float chargedIsoPtSum = tau->tauID(chargedIsoPtSums_);
    float neutralIsoPtSum = tau->tauID(neutralIsoPtSums_);
    float puCorrPtSum     = tau->tauID(puCorrPtSums_);
    float photonPtSumOutsideSignalCone = tau->tauID(photonPtSumOutsideSignalCone_);
    float footprintCorrection = tau->tauID(footprintCorrection_);
		
    float decayDistX = tau->flightLength().x();
    float decayDistY = tau->flightLength().y();
    float decayDistZ = tau->flightLength().z();
    float decayDistMag = std::sqrt(decayDistX*decayDistX + decayDistY*decayDistY + decayDistZ*decayDistZ);
		
    // --- The following 5 variables differ slightly between AOD & MiniAOD
    //     because they are recomputed using packedCandidates saved in the tau
    float nPhoton = (float)clusterVariables_.tau_n_photons_total(*tau);
    float ptWeightedDetaStrip = clusterVariables_.tau_pt_weighted_deta_strip(*tau, tauDecayMode);
    float ptWeightedDphiStrip = clusterVariables_.tau_pt_weighted_dphi_strip(*tau, tauDecayMode);
    float ptWeightedDrSignal = clusterVariables_.tau_pt_weighted_dr_signal(*tau, tauDecayMode);
    float ptWeightedDrIsolation = clusterVariables_.tau_pt_weighted_dr_iso(*tau, tauDecayMode);
    // ---
    float leadingTrackChi2 = tau->leadingTrackNormChi2();
    float eRatio = clusterVariables_.tau_Eratio(*tau);

    // Difference between measured and maximally allowed Gottfried-Jackson angle
    float gjAngleDiff = -999;
    if ( tauDecayMode == 10 ) {
        double mTau = 1.77682;
        double mAOne = tau->p4().M();
        double pAOneMag = tau->p();
        double argumentThetaGJmax = (std::pow(mTau,2) - std::pow(mAOne,2) ) / ( 2 * mTau * pAOneMag );
        double argumentThetaGJmeasured = ( tau->p4().px() * decayDistX + tau->p4().py() * decayDistY + tau->p4().pz() * decayDistZ ) / ( pAOneMag * decayDistMag );
        if ( std::abs(argumentThetaGJmax) <= 1. && std::abs(argumentThetaGJmeasured) <= 1. ) {
            double thetaGJmax = std::asin( argumentThetaGJmax );
            double thetaGJmeasured = std::acos( argumentThetaGJmeasured );
            gjAngleDiff = thetaGJmeasured - thetaGJmax;
        }
    }
		
    if ( mvaOpt_ == kOldDMwoLT || mvaOpt_ == kNewDMwoLT ) {
      mvaInput_[0]  = std::log(std::max(1.f, (float)tau->pt()));
      mvaInput_[1]  = std::abs((float)tau->eta());
      mvaInput_[2]  = std::log(std::max(1.e-2f, chargedIsoPtSum));
      mvaInput_[3]  = std::log(std::max(1.e-2f, neutralIsoPtSum - 0.125f*puCorrPtSum));
      mvaInput_[4]  = std::log(std::max(1.e-2f, puCorrPtSum));
      mvaInput_[5]  = tauDecayMode;
    } else if ( mvaOpt_ == kOldDMwLT || mvaOpt_ == kNewDMwLT  ) {
      mvaInput_[0]  = std::log(std::max(1.f, (float)tau->pt()));
      mvaInput_[1]  = std::abs((float)tau->eta());
      mvaInput_[2]  = std::log(std::max(1.e-2f, chargedIsoPtSum));
      mvaInput_[3]  = std::log(std::max(1.e-2f, neutralIsoPtSum - 0.125f*puCorrPtSum));
      mvaInput_[4]  = std::log(std::max(1.e-2f, puCorrPtSum));
      mvaInput_[5]  = tauDecayMode;
      mvaInput_[6]  = std::copysign(+1.f, tau->dxy());
      mvaInput_[7]  = std::sqrt(std::min(1.f, std::abs(tau->dxy())));
      mvaInput_[8]  = std::min(10.f, std::abs(tau->dxy_Sig()));
      mvaInput_[9]  = ( tau->hasSecondaryVertex() ) ? 1. : 0.;
      mvaInput_[10] = std::sqrt(decayDistMag);
      mvaInput_[11] = std::min(10.f, tau->flightLengthSig());
    } else if ( mvaOpt_ == kDBoldDMwLT || mvaOpt_ == kDBnewDMwLT ) {
      mvaInput_[0]  = std::log(std::max(1.f, (float)tau->pt()));
      mvaInput_[1]  = std::abs((float)tau->eta());
      mvaInput_[2]  = std::log(std::max(1.e-2f, chargedIsoPtSum));
      mvaInput_[3]  = std::log(std::max(1.e-2f, neutralIsoPtSum));
      mvaInput_[4]  = std::log(std::max(1.e-2f, puCorrPtSum));
      mvaInput_[5]  = std::log(std::max(1.e-2f, photonPtSumOutsideSignalCone));
      mvaInput_[6]  = tauDecayMode;
      mvaInput_[7]  = std::min(30.f, nPhoton);
      mvaInput_[8]  = std::min(0.5f, ptWeightedDetaStrip);
      mvaInput_[9]  = std::min(0.5f, ptWeightedDphiStrip);
      mvaInput_[10] = std::min(0.5f, ptWeightedDrSignal);
      mvaInput_[11] = std::min(0.5f, ptWeightedDrIsolation);
      mvaInput_[12] = std::min(100.f, leadingTrackChi2);
      mvaInput_[13] = std::min(1.f, eRatio);
      mvaInput_[14]  = std::copysign(+1.f, tau->dxy());
      mvaInput_[15]  = std::sqrt(std::min(1.f, std::abs(tau->dxy())));
      mvaInput_[16]  = std::min(10.f, std::abs(tau->dxy_Sig()));
      mvaInput_[17]  = std::copysign(+1.f, tau->ip3d());
      mvaInput_[18]  = std::sqrt(std::min(1.f, std::abs(tau->ip3d())));
      mvaInput_[19]  = std::min(10.f, std::abs(tau->ip3d_Sig()));
      mvaInput_[20]  = ( tau->hasSecondaryVertex() ) ? 1. : 0.;
      mvaInput_[21] = std::sqrt(decayDistMag);
      mvaInput_[22] = std::min(10.f, tau->flightLengthSig());
    } else if ( mvaOpt_ == kPWoldDMwLT || mvaOpt_ == kPWnewDMwLT ) {
      mvaInput_[0]  = std::log(std::max(1.f, (float)tau->pt()));
      mvaInput_[1]  = std::abs((float)tau->eta());
      mvaInput_[2]  = std::log(std::max(1.e-2f, chargedIsoPtSum));
      mvaInput_[3]  = std::log(std::max(1.e-2f, neutralIsoPtSum));
      mvaInput_[4]  = std::log(std::max(1.e-2f, footprintCorrection));
      mvaInput_[5]  = std::log(std::max(1.e-2f, photonPtSumOutsideSignalCone));
      mvaInput_[6]  = tauDecayMode;
      mvaInput_[7]  = std::min(30.f, nPhoton);
      mvaInput_[8]  = std::min(0.5f, ptWeightedDetaStrip);
      mvaInput_[9]  = std::min(0.5f, ptWeightedDphiStrip);
      mvaInput_[10] = std::min(0.5f, ptWeightedDrSignal);
      mvaInput_[11] = std::min(0.5f, ptWeightedDrIsolation);
      mvaInput_[12] = std::min(100.f, leadingTrackChi2);
      mvaInput_[13] = std::min(1.f, eRatio);
      mvaInput_[14]  = std::copysign(+1.f, tau->dxy());
      mvaInput_[15]  = std::sqrt(std::min(1.f, std::abs(tau->dxy())));
      mvaInput_[16]  = std::min(10.f, std::abs(tau->dxy_Sig()));
      mvaInput_[17]  = std::copysign(+1.f, tau->ip3d());
      mvaInput_[18]  = std::sqrt(std::min(1.f, std::abs(tau->ip3d())));
      mvaInput_[19]  = std::min(10.f, std::abs(tau->ip3d_Sig()));
      mvaInput_[20]  = ( tau->hasSecondaryVertex() ) ? 1. : 0.;
      mvaInput_[21] = std::sqrt(decayDistMag);
      mvaInput_[22] = std::min(10.f, tau->flightLengthSig());
    } else if ( mvaOpt_ == kDBoldDMwLTwGJ || mvaOpt_ == kDBnewDMwLTwGJ ) {
      mvaInput_[0]  = std::log(std::max(1.f, (float)tau->pt()));
      mvaInput_[1]  = std::abs((float)tau->eta());
      mvaInput_[2]  = std::log(std::max(1.e-2f, chargedIsoPtSum));
      mvaInput_[3]  = std::log(std::max(1.e-2f, neutralIsoPtSum));
      mvaInput_[4]  = std::log(std::max(1.e-2f, puCorrPtSum));
      mvaInput_[5]  = std::log(std::max(1.e-2f, photonPtSumOutsideSignalCone));
      mvaInput_[6]  = tauDecayMode;
      mvaInput_[7]  = std::min(30.f, nPhoton);
      mvaInput_[8]  = std::min(0.5f, ptWeightedDetaStrip);
      mvaInput_[9]  = std::min(0.5f, ptWeightedDphiStrip);
      mvaInput_[10] = std::min(0.5f, ptWeightedDrSignal);
      mvaInput_[11] = std::min(0.5f, ptWeightedDrIsolation);
      mvaInput_[12] = std::min(1.f, eRatio);
      mvaInput_[13]  = std::copysign(+1.f, tau->dxy());
      mvaInput_[14]  = std::sqrt(std::min(1.f, std::abs(tau->dxy())));
      mvaInput_[15]  = std::min(10.f, std::abs(tau->dxy_Sig()));
      mvaInput_[16]  = std::copysign(+1.f, tau->ip3d());
      mvaInput_[17]  = std::sqrt(std::min(1.f, std::abs(tau->ip3d())));
      mvaInput_[18]  = std::min(10.f, std::abs(tau->ip3d_Sig()));
      mvaInput_[19]  = ( tau->hasSecondaryVertex() ) ? 1. : 0.;
      mvaInput_[20] = std::sqrt(decayDistMag);
      mvaInput_[21] = std::min(10.f, tau->flightLengthSig());
      mvaInput_[22] = std::max(-1.f, gjAngleDiff);
    }
		
    double mvaValue = mvaReader_->GetClassifier(mvaInput_);
    if ( verbosity_ ) {
      edm::LogPrint("PATTauDiscByMVAIsolRun2") << "<PATTauDiscriminationByMVAIsolationRun2::discriminate>:";
      edm::LogPrint("PATTauDiscByMVAIsolRun2") << " tau: Pt = " << tau->pt() << ", eta = " << tau->eta();
      edm::LogPrint("PATTauDiscByMVAIsolRun2") << " isolation: charged = " << chargedIsoPtSum << ", neutral = " << neutralIsoPtSum << ", PUcorr = " << puCorrPtSum;
      edm::LogPrint("PATTauDiscByMVAIsolRun2") << " decay mode = " << tauDecayMode;
      edm::LogPrint("PATTauDiscByMVAIsolRun2") << " impact parameter: distance = " << tau->dxy() << ", significance = " << tau->dxy_Sig();
      edm::LogPrint("PATTauDiscByMVAIsolRun2") << " has decay vertex = " << tau->hasSecondaryVertex() << ":"
                                               << " distance = " << decayDistMag << ", significance = " << tau->flightLengthSig();
      edm::LogPrint("PATTauDiscByMVAIsolRun2") << "--> mvaValue = " << mvaValue;
    }
    return mvaValue;
  } else {
    return -1.;
  }
}

void PATTauDiscriminationByMVAIsolationRun2::endEvent(edm::Event& evt)
{
  // add all category indices to event
  evt.put(std::move(category_output_), "category");
}

DEFINE_FWK_MODULE(PATTauDiscriminationByMVAIsolationRun2);
