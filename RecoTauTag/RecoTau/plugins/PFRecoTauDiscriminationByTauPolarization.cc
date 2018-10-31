#include <FWCore/ParameterSet/interface/ConfigurationDescriptions.h>
#include <FWCore/ParameterSet/interface/ParameterSetDescription.h>
#include "RecoTauTag/RecoTau/interface/TauDiscriminationProducerBase.h"

/* class PFRecoTauDiscriminationByTauPolarization
 * created : May 26 2010,
 * contributors : Sami Lehti (sami.lehti@cern.ch ; HIP, Helsinki)
 */

using namespace reco;
using namespace std;
using namespace edm;

class PFRecoTauDiscriminationByTauPolarization :
  public PFTauDiscriminationProducerBase  {
  public:
    explicit PFRecoTauDiscriminationByTauPolarization(
        const ParameterSet& iConfig)
      :PFTauDiscriminationProducerBase(iConfig) {  // retrieve quality cuts
        rTauMin = iConfig.getParameter<double>("rtau");
        booleanOutput = iConfig.getParameter<bool>("BooleanOutput");
      }

    ~PFRecoTauDiscriminationByTauPolarization() override{}

    void beginEvent(const Event&, const EventSetup&) override;
    double discriminate(const PFTauRef&) const override;

    static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);
  private:
    bool booleanOutput;
    double rTauMin;
};

void PFRecoTauDiscriminationByTauPolarization::beginEvent(
    const Event& event, const EventSetup& eventSetup){}

double
PFRecoTauDiscriminationByTauPolarization::discriminate(const PFTauRef& tau) const{

  double rTau = 0;
  // rtau for PFTau has to be calculated for leading PF charged hadronic candidate
  // calculating it from leadingTrack can (and will) give rtau > 1!
  if(tau.isNonnull() && tau->p() > 0
      && tau->leadPFChargedHadrCand().isNonnull()) {
    rTau = tau->leadPFChargedHadrCand()->p()/tau->p();
  }

  if(booleanOutput) return ( rTau > rTauMin ? 1. : 0. );
  return rTau;
}

void
PFRecoTauDiscriminationByTauPolarization::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // pfRecoTauDiscriminationByTauPolarization
  edm::ParameterSetDescription desc;
  desc.add<double>("rtau", 0.8);
  desc.add<edm::InputTag>("PVProducer", edm::InputTag("offlinePrimaryVertices"));
  desc.add<bool>("BooleanOutput", true);
  desc.add<edm::InputTag>("PFTauProducer", edm::InputTag("pfRecoTauProducer"));
  {
    edm::ParameterSetDescription psd0;
    {
      edm::ParameterSetDescription psd1;
      psd1.add<double>("maxDeltaZ", 0.4);
      psd1.add<double>("minTrackPt", 0.5);
      psd1.add<double>("minTrackVertexWeight", -1.0);
      psd1.add<double>("maxTrackChi2", 100.0);
      psd1.add<unsigned int>("minTrackPixelHits", 0);
      psd1.add<double>("minGammaEt", 1.0);
      psd1.add<unsigned int>("minTrackHits", 3);
      psd1.add<double>("minNeutralHadronEt", 30.0);
      psd1.add<double>("maxTransverseImpactParameter", 0.1);
      psd0.add<edm::ParameterSetDescription>("signalQualityCuts", psd1);
    }
    {
      edm::ParameterSetDescription psd1;
      psd1.add<double>("minTrackPt", 0.5);
      psd1.add<double>("minTrackVertexWeight", -1.0);
      psd1.add<double>("maxTrackChi2", 100.0);
      psd1.add<unsigned int>("minTrackPixelHits", 0);
      psd1.add<double>("minGammaEt", 1.0);
      psd1.add<unsigned int>("minTrackHits", 3);
      psd1.add<double>("maxTransverseImpactParameter", 0.1);
      psd0.add<edm::ParameterSetDescription>("vxAssocQualityCuts", psd1);
    }
    psd0.add<std::string>("leadingTrkOrPFCandOption", "leadPFCand");
    {
      edm::ParameterSetDescription psd1;
      psd1.add<double>("maxDeltaZ", 0.2);
      psd1.add<double>("minTrackPt", 1.0);
      psd1.add<double>("minTrackVertexWeight", -1.0);
      psd1.add<double>("maxTrackChi2", 100.0);
      psd1.add<unsigned int>("minTrackPixelHits", 0);
      psd1.add<double>("minGammaEt", 1.5);
      psd1.add<unsigned int>("minTrackHits", 8);
      psd1.add<double>("maxTransverseImpactParameter", 0.03);
      psd0.add<edm::ParameterSetDescription>("isolationQualityCuts", psd1);
    }
    psd0.add<std::string>("pvFindingAlgo", "closestInDeltaZ");
    psd0.add<edm::InputTag>("primaryVertexSrc", edm::InputTag("offlinePrimaryVertices"));
    psd0.add<bool>("vertexTrackFiltering", false);
    psd0.add<bool>("recoverLeadingTrk", false);
    desc.add<edm::ParameterSetDescription>("qualityCuts", psd0);
  }
  {
    edm::ParameterSetDescription psd0;
    psd0.add<std::string>("BooleanOperator", "and");
    {
      edm::ParameterSetDescription psd1;
      psd1.add<double>("cut", 0.5);
      psd1.add<edm::InputTag>("Producer", edm::InputTag("pfRecoTauDiscriminationByLeadingTrackFinding"));
      psd0.add<edm::ParameterSetDescription>("leadTrack", psd1);
    }
    desc.add<edm::ParameterSetDescription>("Prediscriminants", psd0);
  }
  descriptions.add("pfRecoTauDiscriminationByTauPolarization", desc);
}

DEFINE_FWK_MODULE(PFRecoTauDiscriminationByTauPolarization);
