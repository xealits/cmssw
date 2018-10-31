#include "RecoTauTag/RecoTau/interface/TauDiscriminationProducerBase.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include <FWCore/ParameterSet/interface/ConfigurationDescriptions.h>
#include <FWCore/ParameterSet/interface/ParameterSetDescription.h>
#include "RecoTauTag/RecoTau/interface/RecoTauQualityCuts.h"
#include "RecoTauTag/RecoTau/interface/RecoTauVertexAssociator.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

/* class PFRecoTauDiscriminationByNProngs
 * created : August 30 2010,
 * contributors : Sami Lehti (sami.lehti@cern.ch ; HIP, Helsinki)
 * based on H+ tau ID by Lauri Wendland
 * Modified April 16 2014 by S.Lehti
 */

using namespace reco;
using namespace std;
using namespace edm;

class PFRecoTauDiscriminationByNProngs : public PFTauDiscriminationProducerBase  {
    public:
	explicit PFRecoTauDiscriminationByNProngs(const ParameterSet&);
      	~PFRecoTauDiscriminationByNProngs() override{}

	void beginEvent(const edm::Event&, const edm::EventSetup&) override;
	double discriminate(const reco::PFTauRef&) const override;

        static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);
    private:
	std::unique_ptr<tau::RecoTauQualityCuts> qcuts_;
	std::unique_ptr<tau::RecoTauVertexAssociator> vertexAssociator_;

	uint32_t minN,maxN;
	bool booleanOutput;
	edm::ParameterSet qualityCuts;
};

PFRecoTauDiscriminationByNProngs::PFRecoTauDiscriminationByNProngs(const ParameterSet& iConfig):
  PFTauDiscriminationProducerBase(iConfig),
  qualityCuts(iConfig.getParameterSet("qualityCuts"))
{
        minN          = iConfig.getParameter<uint32_t>("MinN");
	maxN          = iConfig.getParameter<uint32_t>("MaxN");
        booleanOutput = iConfig.getParameter<bool>("BooleanOutput");

	qcuts_.reset(new tau::RecoTauQualityCuts(qualityCuts.getParameterSet("signalQualityCuts")));
	vertexAssociator_.reset(new tau::RecoTauVertexAssociator(qualityCuts,consumesCollector()));
}

void PFRecoTauDiscriminationByNProngs::beginEvent(const Event& iEvent, const EventSetup& iSetup){
	vertexAssociator_->setEvent(iEvent);
}

double PFRecoTauDiscriminationByNProngs::discriminate(const PFTauRef& tau) const{

	reco::VertexRef pv = vertexAssociator_->associatedVertex(*tau);
	const PFCandidatePtr leadingTrack = tau->leadPFChargedHadrCand();

	uint np = 0;
	if(leadingTrack.isNonnull() && pv.isNonnull()){
	    qcuts_->setPV(pv);
	    qcuts_->setLeadTrack(tau->leadPFChargedHadrCand());

	    for(auto const& cand : tau->signalPFChargedHadrCands() ) {
	        if ( qcuts_->filterCandRef(cand) ) np++;
 	    }
	}

	bool accepted = false;
	if(maxN == 0){
	    if(np == 1 || np == 3) accepted = true;
	}else{
	    if(np >= minN && np <= maxN) accepted = true;
	}

	if(!accepted) np = 0;
	if(booleanOutput) return accepted;
	return np;
}

void
PFRecoTauDiscriminationByNProngs::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // pfRecoTauDiscriminationByNProngs
  edm::ParameterSetDescription desc;
  desc.add<unsigned int>("MinN", 1);
  desc.add<bool>("BooleanOutput", true);
  desc.add<edm::InputTag>("PFTauProducer", edm::InputTag("combinatoricRecoTaus"));
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
    desc.add<edm::ParameterSetDescription>("Prediscriminants", psd0);
  }
  desc.add<unsigned int>("MaxN", 0);
  descriptions.add("pfRecoTauDiscriminationByNProngs", desc);
}

DEFINE_FWK_MODULE(PFRecoTauDiscriminationByNProngs);
