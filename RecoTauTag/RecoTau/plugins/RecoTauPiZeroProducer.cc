/*
 * RecoTauPiZeroProducer
 *
 * Author: Evan K. Friis, UC Davis
 *
 * Associates reconstructed PiZeros to PFJets.  The PiZeros are built using one
 * or more RecoTauBuilder plugins.  Any overlaps (PiZeros sharing constituents)
 * are removed, with the best PiZero candidates taken.  The 'best' are defined
 * via the input list of RecoTauPiZeroQualityPlugins, which form a
 * lexicograpical ranking.
 *
 */

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/ptr_container/ptr_list.hpp>
#include <algorithm>
#include <functional>

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include <FWCore/ParameterSet/interface/ConfigurationDescriptions.h>
#include <FWCore/ParameterSet/interface/ParameterSetDescription.h>

#include "RecoTauTag/RecoTau/interface/RecoTauPiZeroPlugins.h"
#include "RecoTauTag/RecoTau/interface/RecoTauCleaningTools.h"
#include "RecoTauTag/RecoTau/interface/RecoTauCommonUtilities.h"

#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/TauReco/interface/JetPiZeroAssociation.h"
#include "DataFormats/TauReco/interface/RecoTauPiZero.h"
#include "DataFormats/Common/interface/Association.h"

#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

class RecoTauPiZeroProducer : public edm::stream::EDProducer<> {
  public:
    typedef reco::tau::RecoTauPiZeroBuilderPlugin Builder;
    typedef reco::tau::RecoTauPiZeroQualityPlugin Ranker;

    explicit RecoTauPiZeroProducer(const edm::ParameterSet& pset);
    ~RecoTauPiZeroProducer() override {}
    void produce(edm::Event& evt, const edm::EventSetup& es) override;
    void print(const std::vector<reco::RecoTauPiZero>& piZeros,
               std::ostream& out);

    static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);

  private:
    typedef boost::ptr_vector<Builder> builderList;
    typedef boost::ptr_vector<Ranker> rankerList;
    typedef boost::ptr_vector<reco::RecoTauPiZero> PiZeroVector;
    typedef boost::ptr_list<reco::RecoTauPiZero> PiZeroList;

    typedef reco::tau::RecoTauLexicographicalRanking<rankerList,
            reco::RecoTauPiZero> PiZeroPredicate;

    builderList builders_;
    rankerList rankers_;
    std::auto_ptr<PiZeroPredicate> predicate_;
    double piZeroMass_;

    // Output selector
    std::auto_ptr<StringCutObjectSelector<reco::RecoTauPiZero> >
      outputSelector_;

    //consumes interface
    edm::EDGetTokenT<reco::CandidateView> cand_token;

    double minJetPt_;
    double maxJetAbsEta_;

    int verbosity_;
};

RecoTauPiZeroProducer::RecoTauPiZeroProducer(const edm::ParameterSet& pset) 
{
  cand_token = consumes<reco::CandidateView>( pset.getParameter<edm::InputTag>("jetSrc"));
  minJetPt_ = pset.getParameter<double>("minJetPt");
  maxJetAbsEta_ = pset.getParameter<double>("maxJetAbsEta");

  typedef std::vector<edm::ParameterSet> VPSet;
  // Get the mass hypothesis for the pizeros
  piZeroMass_ = pset.getParameter<double>("massHypothesis");

  // Get each of our PiZero builders
  const VPSet& builders = pset.getParameter<VPSet>("builders");

  for (VPSet::const_iterator builderPSet = builders.begin();
      builderPSet != builders.end(); ++builderPSet) {
    // Get plugin name
    const std::string& pluginType =
      builderPSet->getParameter<std::string>("plugin");
    // Build the plugin
    builders_.push_back(RecoTauPiZeroBuilderPluginFactory::get()->create(
          pluginType, *builderPSet, consumesCollector()));
  }

  // Get each of our quality rankers
  const VPSet& rankers = pset.getParameter<VPSet>("ranking");
  for (VPSet::const_iterator rankerPSet = rankers.begin();
      rankerPSet != rankers.end(); ++rankerPSet) {
    const std::string& pluginType =
      rankerPSet->getParameter<std::string>("plugin");
    rankers_.push_back(RecoTauPiZeroQualityPluginFactory::get()->create(
          pluginType, *rankerPSet));
  }

  // Build the sorting predicate
  predicate_ = std::auto_ptr<PiZeroPredicate>(new PiZeroPredicate(rankers_));

  // now all producers apply a final output selection
  std::string selection = pset.getParameter<std::string>("outputSelection");
  if (!selection.empty()) {
    outputSelector_.reset(
        new StringCutObjectSelector<reco::RecoTauPiZero>(selection));
  }

  verbosity_ = pset.getParameter<int>("verbosity");

  produces<reco::JetPiZeroAssociation>();
}

void RecoTauPiZeroProducer::produce(edm::Event& evt, const edm::EventSetup& es) 
{
  // Get a view of our jets via the base candidates
  edm::Handle<reco::CandidateView> jetView;
  evt.getByToken(cand_token, jetView);

  // Give each of our plugins a chance at doing something with the edm::Event
  for(auto& builder : builders_) {
    builder.setup(evt, es);
  }

  // Convert the view to a RefVector of actual PFJets
  reco::PFJetRefVector jetRefs =
      reco::tau::castView<reco::PFJetRefVector>(jetView);
  // Make our association
  std::unique_ptr<reco::JetPiZeroAssociation> association;

  if (!jetRefs.empty()) {
    edm::Handle<reco::PFJetCollection> pfJetCollectionHandle;
    evt.get(jetRefs.id(), pfJetCollectionHandle);
    association = std::make_unique<reco::JetPiZeroAssociation>(reco::PFJetRefProd(pfJetCollectionHandle));
  } else {
    association = std::make_unique<reco::JetPiZeroAssociation>();
  }

  // Loop over our jets
  for(auto const& jet : jetRefs) {

    if(jet->pt() - minJetPt_ < 1e-5) continue;
    if(std::abs(jet->eta()) - maxJetAbsEta_ > -1e-5) continue;
    // Build our global list of RecoTauPiZero
    PiZeroList dirtyPiZeros;

    // Compute the pi zeros from this jet for all the desired algorithms
    for(auto const& builder : builders_) {
      try {
        PiZeroVector result(builder(*jet));
        dirtyPiZeros.transfer(dirtyPiZeros.end(), result);
      } catch ( cms::Exception &exception) {
        edm::LogError("BuilderPluginException")
            << "Exception caught in builder plugin " << builder.name()
            << ", rethrowing" << std::endl;
        throw exception;
      }
    }
    // Rank the candidates according to our quality plugins
    dirtyPiZeros.sort(*predicate_);

    // Keep track of the photons in the clean collection
    std::vector<reco::RecoTauPiZero> cleanPiZeros;
    std::set<reco::CandidatePtr> photonsInCleanCollection;
    while (!dirtyPiZeros.empty()) {
      // Pull our candidate pi zero from the front of the list
      std::auto_ptr<reco::RecoTauPiZero> toAdd(
          dirtyPiZeros.pop_front().release());
      // If this doesn't pass our basic selection, discard it.
      if (!(*outputSelector_)(*toAdd)) {
        continue;
      }
      // Find the sub-gammas that are not already in the cleaned collection
      std::vector<reco::CandidatePtr> uniqueGammas;
      std::set_difference(toAdd->daughterPtrVector().begin(),
                          toAdd->daughterPtrVector().end(),
                          photonsInCleanCollection.begin(),
                          photonsInCleanCollection.end(),
                          std::back_inserter(uniqueGammas));
      // If the pi zero has no unique gammas, discard it.  Note toAdd is deleted
      // when it goes out of scope.
      if (uniqueGammas.empty()) {
        continue;
      } else if (uniqueGammas.size() == toAdd->daughterPtrVector().size()) {
        // Check if it is composed entirely of unique gammas.  In this case
        // immediately add it to the clean collection.
        photonsInCleanCollection.insert(toAdd->daughterPtrVector().begin(),
                                        toAdd->daughterPtrVector().end());
        cleanPiZeros.push_back(*toAdd);
      } else {
        // Otherwise update the pizero that contains only the unique gammas and
        // add it back into the sorted list of dirty PiZeros
        toAdd->clearDaughters();
        // Add each of the unique daughters back to the pizero
        for(auto const& gamma : uniqueGammas) {
          toAdd->addDaughter(gamma);
        }
        // Update the four vector
        AddFourMomenta p4Builder_;
        p4Builder_.set(*toAdd);
        // Put this pi zero back into the collection of sorted dirty pizeros
        PiZeroList::iterator insertionPoint = std::lower_bound(
            dirtyPiZeros.begin(), dirtyPiZeros.end(), *toAdd, *predicate_);
        dirtyPiZeros.insert(insertionPoint, toAdd);
      }
    }
    // Apply the mass hypothesis if desired
    if (piZeroMass_ >= 0) {
      for( auto& cleanPiZero: cleanPiZeros )
         { cleanPiZero.setMass(this->piZeroMass_);};
    }
    // Add to association
    if ( verbosity_ >= 2 ) {
      print(cleanPiZeros, std::cout);
    }
    association->setValue(jet.key(), cleanPiZeros);
  }
  evt.put(std::move(association));
}

// Print some helpful information
void RecoTauPiZeroProducer::print(
    const std::vector<reco::RecoTauPiZero>& piZeros, std::ostream& out) {
  const unsigned int width = 25;
  for(auto const& piZero : piZeros) {
    out << piZero;
    out << "* Rankers:" << std::endl;
    for (rankerList::const_iterator ranker = rankers_.begin();
        ranker != rankers_.end(); ++ranker) {
      out << "* " << std::setiosflags(std::ios::left)
        << std::setw(width) << ranker->name()
        << " " << std::resetiosflags(std::ios::left)
        << std::setprecision(3) << (*ranker)(piZero);
      out << std::endl;
    }
  }
}

void
RecoTauPiZeroProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // common parameter descriptions
  edm::ParameterSetDescription vpsd_ranking;
  vpsd_ranking.add<std::string>("selectionPassFunction");
  vpsd_ranking.add<double>("selectionFailValue");
  vpsd_ranking.add<std::string>("selection");
  vpsd_ranking.add<std::string>("name");
  vpsd_ranking.add<std::string>("plugin");

  {
    // ak4PFJetsLegacyTaNCPiZeros
    edm::ParameterSetDescription desc;
    desc.add<double>("massHypothesis", 0.136);
    {
      std::vector<edm::ParameterSet> ranking_vactor;
      ranking_vactor.reserve(1);
      {
        edm::ParameterSet temp2;
        temp2.addParameter<std::string>("selectionPassFunction", "abs(mass() - 0.13579)");
        temp2.addParameter<double>("selectionFailValue", 1000);
        temp2.addParameter<std::string>("selection", "mass() < 0.2");
        temp2.addParameter<std::string>("name", "PFTDM");
        temp2.addParameter<std::string>("plugin", "RecoTauPiZeroStringQuality");
        ranking_vactor.push_back(temp2);
      }
      desc.addVPSet("ranking", vpsd_ranking, ranking_vactor);
    }
    desc.add<int>("verbosity", 0);
    desc.add<double>("maxJetAbsEta", 2.5);
    desc.add<std::string>("outputSelection", "pt > 1.5");
    desc.add<double>("minJetPt", 14.0);
    desc.add<edm::InputTag>("jetSrc", edm::InputTag("ak4PFJets"));
    {
      edm::ParameterSetDescription vpsd_builders;
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
        vpsd_builders.add<edm::ParameterSetDescription>("qualityCuts", psd0);
      }
      vpsd_builders.add<std::string>("name", "1");
      vpsd_builders.add<std::string>("plugin", "RecoTauPiZeroTrivialPlugin");

      vpsd_builders.addOptional<bool>("makeCombinatoricStrips");
      vpsd_builders.addOptional<int>("maxStripBuildIterations");
      vpsd_builders.addOptional<double>("minGammaEtStripAdd");
      vpsd_builders.addOptional<double>("minGammaEtStripSeed");
      vpsd_builders.addOptional<double>("minStripEt");
      vpsd_builders.addOptional<std::vector<int>>("stripCandidatesParticleIds");
      vpsd_builders.addOptional<bool>("updateStripAfterEachDaughter");
      vpsd_builders.addOptional<int>("verbosity");
      vpsd_builders.addOptional<bool>("applyElecTrackQcuts");

      std::vector<edm::ParameterSet> builders_vector;
      builders_vector.reserve(2);
      {
        edm::ParameterSet temp2;
        {
          edm::ParameterSet ps2;
          {
            edm::ParameterSet ps3;
            ps3.addParameter<double>("maxDeltaZ", 0.4);
            ps3.addParameter<double>("minTrackPt", 0.5);
            ps3.addParameter<double>("minTrackVertexWeight", -1.0);
            ps3.addParameter<double>("maxTrackChi2", 100.0);
            ps3.addParameter<unsigned int>("minTrackPixelHits", 0);
            ps3.addParameter<double>("minGammaEt", 1.0);
            ps3.addParameter<unsigned int>("minTrackHits", 3);
            ps3.addParameter<double>("minNeutralHadronEt", 30.0);
            ps3.addParameter<double>("maxTransverseImpactParameter", 0.1);
            ps2.addParameter<edm::ParameterSet>("signalQualityCuts", ps3);
          }
          {
            edm::ParameterSet ps3;
            ps3.addParameter<double>("minTrackPt", 0.5);
            ps3.addParameter<double>("minTrackVertexWeight", -1.0);
            ps3.addParameter<double>("maxTrackChi2", 100.0);
            ps3.addParameter<unsigned int>("minTrackPixelHits", 0);
            ps3.addParameter<double>("minGammaEt", 1.0);
            ps3.addParameter<unsigned int>("minTrackHits", 3);
            ps3.addParameter<double>("maxTransverseImpactParameter", 0.1);
            ps2.addParameter<edm::ParameterSet>("vxAssocQualityCuts", ps3);
          }
          ps2.addParameter<std::string>("leadingTrkOrPFCandOption", "leadPFCand");
          {
            edm::ParameterSet ps3;
            ps3.addParameter<double>("maxDeltaZ", 0.2);
            ps3.addParameter<double>("minTrackPt", 1.0);
            ps3.addParameter<double>("minTrackVertexWeight", -1.0);
            ps3.addParameter<double>("maxTrackChi2", 100.0);
            ps3.addParameter<unsigned int>("minTrackPixelHits", 0);
            ps3.addParameter<double>("minGammaEt", 1.5);
            ps3.addParameter<unsigned int>("minTrackHits", 8);
            ps3.addParameter<double>("maxTransverseImpactParameter", 0.03);
            ps2.addParameter<edm::ParameterSet>("isolationQualityCuts", ps3);
          }
          ps2.addParameter<std::string>("pvFindingAlgo", "closestInDeltaZ");
          ps2.addParameter<edm::InputTag>("primaryVertexSrc", edm::InputTag("offlinePrimaryVertices"));
          ps2.addParameter<bool>("vertexTrackFiltering", false);
          ps2.addParameter<bool>("recoverLeadingTrk", false);
          temp2.addParameter<edm::ParameterSet>("qualityCuts", ps2);
        }
        temp2.addParameter<std::string>("name", "1");
        temp2.addParameter<std::string>("plugin", "RecoTauPiZeroTrivialPlugin");
        builders_vector.push_back(temp2);
      }
      {
        edm::ParameterSet temp2;
        temp2.addParameter<double>("maxMass", -1.0);
        temp2.addParameter<std::string>("plugin", "RecoTauPiZeroCombinatoricPlugin");
        temp2.addParameter<double>("minMass", 0.0);
        {
          edm::ParameterSet ps2;
          {
            edm::ParameterSet ps3;
            ps3.addParameter<double>("maxDeltaZ", 0.4);
            ps3.addParameter<double>("minTrackPt", 0.5);
            ps3.addParameter<double>("minTrackVertexWeight", -1.0);
            ps3.addParameter<double>("maxTrackChi2", 100.0);
            ps3.addParameter<unsigned int>("minTrackPixelHits", 0);
            ps3.addParameter<double>("minGammaEt", 1.0);
            ps3.addParameter<unsigned int>("minTrackHits", 3);
            ps3.addParameter<double>("minNeutralHadronEt", 30.0);
            ps3.addParameter<double>("maxTransverseImpactParameter", 0.1);
            ps2.addParameter<edm::ParameterSet>("signalQualityCuts", ps3);
          }
          {
            edm::ParameterSet ps3;
            ps3.addParameter<double>("minTrackPt", 0.5);
            ps3.addParameter<double>("minTrackVertexWeight", -1.0);
            ps3.addParameter<double>("maxTrackChi2", 100.0);
            ps3.addParameter<unsigned int>("minTrackPixelHits", 0);
            ps3.addParameter<double>("minGammaEt", 1.0);
            ps3.addParameter<unsigned int>("minTrackHits", 3);
            ps3.addParameter<double>("maxTransverseImpactParameter", 0.1);
            ps2.addParameter<edm::ParameterSet>("vxAssocQualityCuts", ps3);
          }
          ps2.addParameter<std::string>("leadingTrkOrPFCandOption", "leadPFCand");
          {
            edm::ParameterSet ps3;
            ps3.addParameter<double>("maxDeltaZ", 0.2);
            ps3.addParameter<double>("minTrackPt", 1.0);
            ps3.addParameter<double>("minTrackVertexWeight", -1.0);
            ps3.addParameter<double>("maxTrackChi2", 100.0);
            ps3.addParameter<unsigned int>("minTrackPixelHits", 0);
            ps3.addParameter<double>("minGammaEt", 1.5);
            ps3.addParameter<unsigned int>("minTrackHits", 8);
            ps3.addParameter<double>("maxTransverseImpactParameter", 0.03);
            ps2.addParameter<edm::ParameterSet>("isolationQualityCuts", ps3);
          }
          ps2.addParameter<std::string>("pvFindingAlgo", "closestInDeltaZ");
          ps2.addParameter<edm::InputTag>("primaryVertexSrc", edm::InputTag("offlinePrimaryVertices"));
          ps2.addParameter<bool>("vertexTrackFiltering", false);
          ps2.addParameter<bool>("recoverLeadingTrk", false);
          temp2.addParameter<edm::ParameterSet>("qualityCuts", ps2);
        }
        temp2.addParameter<unsigned int>("choose", 2);
        temp2.addParameter<unsigned int>("maxInputGammas", 10);
        temp2.addParameter<std::string>("name", "2");
        builders_vector.push_back(temp2);
      }
      desc.addVPSet("builders", vpsd_builders, builders_vector);
    }
    descriptions.add("ak4PFJetsLegacyTaNCPiZeros", desc);
  }

  {
    // ak4PFJetsRecoTauGreedyPiZeros
    edm::ParameterSetDescription desc;
    desc.add<double>("massHypothesis", 0.136);
    {
      std::vector<edm::ParameterSet> ranking_vector;
      ranking_vector.reserve(1);
      {
        edm::ParameterSet temp2;
        temp2.addParameter<std::string>("selectionPassFunction", "-1*numberOfDaughters()");
        temp2.addParameter<double>("selectionFailValue", 1000);
        temp2.addParameter<std::string>("selection", "algoIs(\"kStrips\") || (abs(eta()) < 1.5 & abs(mass() - 0.13579) < 0.05) || (abs(eta()) > 1.5 & mass < 0.2)");
        temp2.addParameter<std::string>("name", "Greedy");
        temp2.addParameter<std::string>("plugin", "RecoTauPiZeroStringQuality");
        ranking_vector.push_back(temp2);
      }
      desc.addVPSet("ranking", vpsd_ranking, ranking_vector);
    }
    desc.add<int>("verbosity", 0);
    desc.add<double>("maxJetAbsEta", 2.5);
    desc.add<std::string>("outputSelection", "pt > 1.5");
    desc.add<double>("minJetPt", 14.0);
    desc.add<edm::InputTag>("jetSrc", edm::InputTag("ak4PFJets"));
    {
      edm::ParameterSetDescription vpsd_builders;
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
        vpsd_builders.add<edm::ParameterSetDescription>("qualityCuts", psd0);
      }
      vpsd_builders.add<int>("maxInputStrips", 5);
      vpsd_builders.add<std::string>("name", "cs");
      vpsd_builders.add<std::string>("plugin", "RecoTauPiZeroStripPlugin");
      vpsd_builders.add<double>("stripMassWhenCombining", 0.0);
      vpsd_builders.add<double>("stripPhiAssociationDistance", 0.2);
      vpsd_builders.add<double>("stripEtaAssociationDistance", 0.05);

      vpsd_builders.addOptional<bool>("makeCombinatoricStrips");
      vpsd_builders.addOptional<int>("maxStripBuildIterations");
      vpsd_builders.addOptional<double>("minGammaEtStripAdd");
      vpsd_builders.addOptional<double>("minGammaEtStripSeed");
      vpsd_builders.addOptional<double>("minStripEt");
      vpsd_builders.addOptional<std::vector<int>>("stripCandidatesParticleIds");
      vpsd_builders.addOptional<bool>("updateStripAfterEachDaughter");
      vpsd_builders.addOptional<int>("verbosity");
      vpsd_builders.addOptional<bool>("applyElecTrackQcuts");
      std::vector<edm::ParameterSet> builders_vector;
      builders_vector.reserve(1);
      {
        edm::ParameterSet temp2;
        temp2.addParameter<int>("maxInputStrips", 5);
        temp2.addParameter<std::string>("name", "cs");
        temp2.addParameter<double>("stripPhiAssociationDistance", 0.2);
        temp2.addParameter<std::string>("plugin", "RecoTauPiZeroStripPlugin");
        {
          edm::ParameterSet ps2;
          {
            edm::ParameterSet ps3;
            ps3.addParameter<double>("maxDeltaZ", 0.4);
            ps3.addParameter<double>("minTrackPt", 0.5);
            ps3.addParameter<double>("minTrackVertexWeight", -1.0);
            ps3.addParameter<double>("maxTrackChi2", 100.0);
            ps3.addParameter<unsigned int>("minTrackPixelHits", 0);
            ps3.addParameter<double>("minGammaEt", 1.0);
            ps3.addParameter<unsigned int>("minTrackHits", 3);
            ps3.addParameter<double>("minNeutralHadronEt", 30.0);
            ps3.addParameter<double>("maxTransverseImpactParameter", 0.1);
            ps2.addParameter<edm::ParameterSet>("signalQualityCuts", ps3);
          }
          {
            edm::ParameterSet ps3;
            ps3.addParameter<double>("minTrackPt", 0.5);
            ps3.addParameter<double>("minTrackVertexWeight", -1.0);
            ps3.addParameter<double>("maxTrackChi2", 100.0);
            ps3.addParameter<unsigned int>("minTrackPixelHits", 0);
            ps3.addParameter<double>("minGammaEt", 1.0);
            ps3.addParameter<unsigned int>("minTrackHits", 3);
            ps3.addParameter<double>("maxTransverseImpactParameter", 0.1);
            ps2.addParameter<edm::ParameterSet>("vxAssocQualityCuts", ps3);
          }
          ps2.addParameter<std::string>("leadingTrkOrPFCandOption", "leadPFCand");
          {
            edm::ParameterSet ps3;
            ps3.addParameter<double>("maxDeltaZ", 0.2);
            ps3.addParameter<double>("minTrackPt", 1.0);
            ps3.addParameter<double>("minTrackVertexWeight", -1.0);
            ps3.addParameter<double>("maxTrackChi2", 100.0);
            ps3.addParameter<unsigned int>("minTrackPixelHits", 0);
            ps3.addParameter<double>("minGammaEt", 1.5);
            ps3.addParameter<unsigned int>("minTrackHits", 8);
            ps3.addParameter<double>("maxTransverseImpactParameter", 0.03);
            ps2.addParameter<edm::ParameterSet>("isolationQualityCuts", ps3);
          }
          ps2.addParameter<std::string>("pvFindingAlgo", "closestInDeltaZ");
          ps2.addParameter<edm::InputTag>("primaryVertexSrc", edm::InputTag("offlinePrimaryVertices"));
          ps2.addParameter<bool>("vertexTrackFiltering", false);
          ps2.addParameter<bool>("recoverLeadingTrk", false);
          temp2.addParameter<edm::ParameterSet>("qualityCuts", ps2);
        }
        temp2.addParameter<double>("stripMassWhenCombining", 0.0);
        temp2.addParameter<bool>("makeCombinatoricStrips", true);
        temp2.addParameter<std::vector<int>>("stripCandidatesParticleIds", {
          2,
          4,
        });
        temp2.addParameter<double>("stripEtaAssociationDistance", 0.05);
        builders_vector.push_back(temp2);
      }
      desc.addVPSet("builders", vpsd_builders, builders_vector);
    }
    descriptions.add("ak4PFJetsRecoTauGreedyPiZeros", desc);
  }

  {
    // ak4PFJetsRecoTauPiZeros
    edm::ParameterSetDescription desc;
    desc.add<double>("massHypothesis", 0.136);
    {
      std::vector<edm::ParameterSet> ranking_vector;
      ranking_vector.reserve(3);
      {
        edm::ParameterSet temp2;
        temp2.addParameter<std::string>("selectionPassFunction", "abs(mass() - 0.13579)");
        temp2.addParameter<double>("selectionFailValue", 1000);
        temp2.addParameter<std::string>("selection", "abs(eta()) < 1.5 & abs(mass() - 0.13579) < 0.05");
        temp2.addParameter<std::string>("name", "nearPiZeroMass");
        temp2.addParameter<std::string>("plugin", "RecoTauPiZeroStringQuality");
        ranking_vector.push_back(temp2);
      }
      {
        edm::ParameterSet temp2;
        temp2.addParameter<std::string>("selectionPassFunction", "abs(mass() - 0.13579)");
        temp2.addParameter<double>("selectionFailValue", 1000);
        temp2.addParameter<std::string>("selection", "abs(eta()) > 1.5 & mass() < 0.2");
        temp2.addParameter<std::string>("name", "nearPiZeroMass");
        temp2.addParameter<std::string>("plugin", "RecoTauPiZeroStringQuality");
        ranking_vector.push_back(temp2);
      }
      {
        edm::ParameterSet temp2;
        temp2.addParameter<std::string>("selectionPassFunction", "abs(mass() - 0.13579)");
        temp2.addParameter<double>("selectionFailValue", 1000);
        temp2.addParameter<std::string>("selection", "algoIs(\"kStrips\")");
        temp2.addParameter<std::string>("name", "InStrip");
        temp2.addParameter<std::string>("plugin", "RecoTauPiZeroStringQuality");
        ranking_vector.push_back(temp2);
      }
      desc.addVPSet("ranking", vpsd_ranking, ranking_vector);
    }
    desc.add<int>("verbosity", 0);
    desc.add<double>("maxJetAbsEta", 2.5);
    desc.add<std::string>("outputSelection", "pt > 1.5");
    desc.add<double>("minJetPt", 14.0);
    desc.add<edm::InputTag>("jetSrc", edm::InputTag("ak4PFJets"));
    {
      edm::ParameterSetDescription vpsd_builders;
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
        vpsd_builders.add<edm::ParameterSetDescription>("qualityCuts", psd0);
      }
      vpsd_builders.add<std::string>("name", "2");
      vpsd_builders.add<std::string>("plugin", "RecoTauPiZeroCombinatoricPlugin");
      vpsd_builders.add<double>("maxMass", -1.0);
      vpsd_builders.add<double>("minMass", 0.0);
      vpsd_builders.add<unsigned int>("choose", 2);
      vpsd_builders.addOptional<unsigned int>("maxInputGammas");

      vpsd_builders.addOptional<bool>("makeCombinatoricStrips");
      vpsd_builders.addOptional<int>("maxStripBuildIterations");
      vpsd_builders.addOptional<double>("minGammaEtStripAdd");
      vpsd_builders.addOptional<double>("minGammaEtStripSeed");
      vpsd_builders.addOptional<double>("minStripEt");
      vpsd_builders.addOptional<std::vector<int>>("stripCandidatesParticleIds");
      vpsd_builders.addOptional<bool>("updateStripAfterEachDaughter");
      vpsd_builders.addOptional<int>("verbosity");
      vpsd_builders.addOptional<bool>("applyElecTrackQcuts");
      std::vector<edm::ParameterSet> builders_vector;
      builders_vector.reserve(2);
      {
        edm::ParameterSet temp2;
        temp2.addParameter<double>("maxMass", -1.0);
        temp2.addParameter<std::string>("plugin", "RecoTauPiZeroCombinatoricPlugin");
        temp2.addParameter<double>("minMass", 0.0);
        {
          edm::ParameterSet ps2;
          {
            edm::ParameterSet ps3;
            ps3.addParameter<double>("maxDeltaZ", 0.4);
            ps3.addParameter<double>("minTrackPt", 0.5);
            ps3.addParameter<double>("minTrackVertexWeight", -1.0);
            ps3.addParameter<double>("maxTrackChi2", 100.0);
            ps3.addParameter<unsigned int>("minTrackPixelHits", 0);
            ps3.addParameter<double>("minGammaEt", 1.0);
            ps3.addParameter<unsigned int>("minTrackHits", 3);
            ps3.addParameter<double>("minNeutralHadronEt", 30.0);
            ps3.addParameter<double>("maxTransverseImpactParameter", 0.1);
            ps2.addParameter<edm::ParameterSet>("signalQualityCuts", ps3);
          }
          {
            edm::ParameterSet ps3;
            ps3.addParameter<double>("minTrackPt", 0.5);
            ps3.addParameter<double>("minTrackVertexWeight", -1.0);
            ps3.addParameter<double>("maxTrackChi2", 100.0);
            ps3.addParameter<unsigned int>("minTrackPixelHits", 0);
            ps3.addParameter<double>("minGammaEt", 1.0);
            ps3.addParameter<unsigned int>("minTrackHits", 3);
            ps3.addParameter<double>("maxTransverseImpactParameter", 0.1);
            ps2.addParameter<edm::ParameterSet>("vxAssocQualityCuts", ps3);
          }
          ps2.addParameter<std::string>("leadingTrkOrPFCandOption", "leadPFCand");
          {
            edm::ParameterSet ps3;
            ps3.addParameter<double>("maxDeltaZ", 0.2);
            ps3.addParameter<double>("minTrackPt", 1.0);
            ps3.addParameter<double>("minTrackVertexWeight", -1.0);
            ps3.addParameter<double>("maxTrackChi2", 100.0);
            ps3.addParameter<unsigned int>("minTrackPixelHits", 0);
            ps3.addParameter<double>("minGammaEt", 1.5);
            ps3.addParameter<unsigned int>("minTrackHits", 8);
            ps3.addParameter<double>("maxTransverseImpactParameter", 0.03);
            ps2.addParameter<edm::ParameterSet>("isolationQualityCuts", ps3);
          }
          ps2.addParameter<std::string>("pvFindingAlgo", "closestInDeltaZ");
          ps2.addParameter<edm::InputTag>("primaryVertexSrc", edm::InputTag("offlinePrimaryVertices"));
          ps2.addParameter<bool>("vertexTrackFiltering", false);
          ps2.addParameter<bool>("recoverLeadingTrk", false);
          temp2.addParameter<edm::ParameterSet>("qualityCuts", ps2);
        }
        temp2.addParameter<unsigned int>("choose", 2);
        temp2.addParameter<unsigned int>("maxInputGammas", 10);
        temp2.addParameter<std::string>("name", "2");
        builders_vector.push_back(temp2);
      }
      {
        edm::ParameterSet temp2;
        temp2.addParameter<std::string>("name", "s");
        {
          edm::ParameterSet ps2;
          ps2.addParameter<std::string>("function", "TMath::Min(0.3, TMath::Max(0.05, [0]*TMath::Power(pT, -[1])))");
          ps2.addParameter<double>("par1", 0.707716);
          ps2.addParameter<double>("par0", 0.352476);
          temp2.addParameter<edm::ParameterSet>("stripPhiAssociationDistance", ps2);
        }
        temp2.addParameter<std::string>("plugin", "RecoTauPiZeroStripPlugin3");
        temp2.addParameter<double>("minGammaEtStripAdd", 1.0);
        temp2.addParameter<double>("minGammaEtStripSeed", 1.0);
        {
          edm::ParameterSet ps2;
          {
            edm::ParameterSet ps3;
            ps3.addParameter<double>("maxDeltaZ", 0.4);
            ps3.addParameter<double>("minTrackPt", 0.5);
            ps3.addParameter<double>("minTrackVertexWeight", -1.0);
            ps3.addParameter<double>("maxTrackChi2", 100.0);
            ps3.addParameter<unsigned int>("minTrackPixelHits", 0);
            ps3.addParameter<double>("minGammaEt", 1.0);
            ps3.addParameter<unsigned int>("minTrackHits", 3);
            ps3.addParameter<double>("minNeutralHadronEt", 30.0);
            ps3.addParameter<double>("maxTransverseImpactParameter", 0.1);
            ps2.addParameter<edm::ParameterSet>("signalQualityCuts", ps3);
          }
          {
            edm::ParameterSet ps3;
            ps3.addParameter<double>("minTrackPt", 0.5);
            ps3.addParameter<double>("minTrackVertexWeight", -1.0);
            ps3.addParameter<double>("maxTrackChi2", 100.0);
            ps3.addParameter<unsigned int>("minTrackPixelHits", 0);
            ps3.addParameter<double>("minGammaEt", 1.0);
            ps3.addParameter<unsigned int>("minTrackHits", 3);
            ps3.addParameter<double>("maxTransverseImpactParameter", 0.1);
            ps2.addParameter<edm::ParameterSet>("vxAssocQualityCuts", ps3);
          }
          ps2.addParameter<std::string>("leadingTrkOrPFCandOption", "leadPFCand");
          {
            edm::ParameterSet ps3;
            ps3.addParameter<double>("maxDeltaZ", 0.2);
            ps3.addParameter<double>("minTrackPt", 1.0);
            ps3.addParameter<double>("minTrackVertexWeight", -1.0);
            ps3.addParameter<double>("maxTrackChi2", 100.0);
            ps3.addParameter<unsigned int>("minTrackPixelHits", 0);
            ps3.addParameter<double>("minGammaEt", 1.5);
            ps3.addParameter<unsigned int>("minTrackHits", 8);
            ps3.addParameter<double>("maxTransverseImpactParameter", 0.03);
            ps2.addParameter<edm::ParameterSet>("isolationQualityCuts", ps3);
          }
          ps2.addParameter<std::string>("pvFindingAlgo", "closestInDeltaZ");
          ps2.addParameter<edm::InputTag>("primaryVertexSrc", edm::InputTag("offlinePrimaryVertices"));
          ps2.addParameter<bool>("vertexTrackFiltering", false);
          ps2.addParameter<bool>("recoverLeadingTrk", false);
          temp2.addParameter<edm::ParameterSet>("qualityCuts", ps2);
        }
        {
          edm::ParameterSet ps2;
          ps2.addParameter<std::string>("function", "TMath::Min(0.15, TMath::Max(0.05, [0]*TMath::Power(pT, -[1])))");
          ps2.addParameter<double>("par1", 0.658701);
          ps2.addParameter<double>("par0", 0.197077);
          temp2.addParameter<edm::ParameterSet>("stripEtaAssociationDistance", ps2);
        }
        temp2.addParameter<int>("maxStripBuildIterations", -1);
        temp2.addParameter<bool>("updateStripAfterEachDaughter", false);
        temp2.addParameter<bool>("makeCombinatoricStrips", false);
        temp2.addParameter<bool>("applyElecTrackQcuts", false);
        temp2.addParameter<std::vector<int>>("stripCandidatesParticleIds", {
          2,
          4,
        });
        temp2.addParameter<double>("minStripEt", 1.0);
        temp2.addParameter<int>("verbosity", 0);
        builders_vector.push_back(temp2);
      }
      desc.addVPSet("builders", vpsd_builders, builders_vector);
    }
    descriptions.add("ak4PFJetsRecoTauPiZeros", desc);
  }

  {
    // ak4PFJetsLegacyHPSPiZeros
    edm::ParameterSetDescription desc;
    desc.add<double>("massHypothesis", 0.136);
    {
      std::vector<edm::ParameterSet> ranking_vector;
      ranking_vector.reserve(1);
      {
        edm::ParameterSet temp2;
        temp2.addParameter<std::string>("selectionPassFunction", "abs(mass() - 0.13579)");
        temp2.addParameter<double>("selectionFailValue", 1000);
        temp2.addParameter<std::string>("selection", "algoIs(\"kStrips\")");
        temp2.addParameter<std::string>("name", "InStrip");
        temp2.addParameter<std::string>("plugin", "RecoTauPiZeroStringQuality");
        ranking_vector.push_back(temp2);
      }
      desc.addVPSet("ranking", vpsd_ranking, ranking_vector);
    }
    desc.add<int>("verbosity", 0);
    desc.add<double>("maxJetAbsEta", 2.5);
    desc.add<std::string>("outputSelection", "pt > 0");
    desc.add<double>("minJetPt", 14.0);
    desc.add<edm::InputTag>("jetSrc", edm::InputTag("ak4PFJets"));
    {
      edm::ParameterSetDescription vpsd_builders;
      {
        edm::ParameterSetDescription psd0;
        psd0.add<std::string>("function", "TMath::Min(0.3, TMath::Max(0.05, [0]*TMath::Power(pT, -[1])))");
        psd0.add<double>("par1", 0.707716);
        psd0.add<double>("par0", 0.352476);
        vpsd_builders.add<edm::ParameterSetDescription>("stripPhiAssociationDistance", psd0);
      }
      {
        edm::ParameterSetDescription psd0;
        psd0.add<std::string>("function", "TMath::Min(0.15, TMath::Max(0.05, [0]*TMath::Power(pT, -[1])))");
        psd0.add<double>("par1", 0.658701);
        psd0.add<double>("par0", 0.197077);
        vpsd_builders.add<edm::ParameterSetDescription>("stripEtaAssociationDistance", psd0);
      }
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
        vpsd_builders.add<edm::ParameterSetDescription>("qualityCuts", psd0);
      }
      vpsd_builders.add<std::string>("name", "s");
      vpsd_builders.add<std::string>("plugin", "RecoTauPiZeroStripPlugin3");

      vpsd_builders.addOptional<bool>("makeCombinatoricStrips");
      vpsd_builders.addOptional<int>("maxStripBuildIterations");
      vpsd_builders.addOptional<double>("minGammaEtStripAdd");
      vpsd_builders.addOptional<double>("minGammaEtStripSeed");
      vpsd_builders.addOptional<double>("minStripEt");
      vpsd_builders.addOptional<std::vector<int>>("stripCandidatesParticleIds");
      vpsd_builders.addOptional<bool>("updateStripAfterEachDaughter");
      vpsd_builders.addOptional<int>("verbosity");
      vpsd_builders.addOptional<bool>("applyElecTrackQcuts");

      std::vector<edm::ParameterSet> builders_vector;
      builders_vector.reserve(1);
      {
        edm::ParameterSet temp2;
        {
          edm::ParameterSet ps2;
          {
            edm::ParameterSet ps3;
            ps3.addParameter<double>("maxDeltaZ", 0.4);
            ps3.addParameter<double>("minTrackPt", 0.5);
            ps3.addParameter<double>("minTrackVertexWeight", -1.0);
            ps3.addParameter<double>("maxTrackChi2", 100.0);
            ps3.addParameter<unsigned int>("minTrackPixelHits", 0);
            ps3.addParameter<double>("minGammaEt", 1.0);
            ps3.addParameter<unsigned int>("minTrackHits", 3);
            ps3.addParameter<double>("minNeutralHadronEt", 30.0);
            ps3.addParameter<double>("maxTransverseImpactParameter", 0.1);
            ps2.addParameter<edm::ParameterSet>("signalQualityCuts", ps3);
          }
          {
            edm::ParameterSet ps3;
            ps3.addParameter<double>("minTrackPt", 0.5);
            ps3.addParameter<double>("minTrackVertexWeight", -1.0);
            ps3.addParameter<double>("maxTrackChi2", 100.0);
            ps3.addParameter<unsigned int>("minTrackPixelHits", 0);
            ps3.addParameter<double>("minGammaEt", 1.0);
            ps3.addParameter<unsigned int>("minTrackHits", 3);
            ps3.addParameter<double>("maxTransverseImpactParameter", 0.1);
            ps2.addParameter<edm::ParameterSet>("vxAssocQualityCuts", ps3);
          }
          ps2.addParameter<std::string>("leadingTrkOrPFCandOption", "leadPFCand");
          {
            edm::ParameterSet ps3;
            ps3.addParameter<double>("maxDeltaZ", 0.2);
            ps3.addParameter<double>("minTrackPt", 1.0);
            ps3.addParameter<double>("minTrackVertexWeight", -1.0);
            ps3.addParameter<double>("maxTrackChi2", 100.0);
            ps3.addParameter<unsigned int>("minTrackPixelHits", 0);
            ps3.addParameter<double>("minGammaEt", 1.5);
            ps3.addParameter<unsigned int>("minTrackHits", 8);
            ps3.addParameter<double>("maxTransverseImpactParameter", 0.03);
            ps2.addParameter<edm::ParameterSet>("isolationQualityCuts", ps3);
          }
          ps2.addParameter<std::string>("pvFindingAlgo", "closestInDeltaZ");
          ps2.addParameter<edm::InputTag>("primaryVertexSrc", edm::InputTag("offlinePrimaryVertices"));
          ps2.addParameter<bool>("vertexTrackFiltering", false);
          ps2.addParameter<bool>("recoverLeadingTrk", false);
          temp2.addParameter<edm::ParameterSet>("qualityCuts", ps2);
        }
        {
          edm::ParameterSet ps2;
          ps2.addParameter<std::string>("function", "TMath::Min(0.3, TMath::Max(0.05, [0]*TMath::Power(pT, -[1])))");
          ps2.addParameter<double>("par1", 0.707716);
          ps2.addParameter<double>("par0", 0.352476);
          temp2.addParameter<edm::ParameterSet>("stripPhiAssociationDistance", ps2);
        }
        {
          edm::ParameterSet ps2;
          ps2.addParameter<std::string>("function", "TMath::Min(0.15, TMath::Max(0.05, [0]*TMath::Power(pT, -[1])))");
          ps2.addParameter<double>("par1", 0.658701);
          ps2.addParameter<double>("par0", 0.197077);
          temp2.addParameter<edm::ParameterSet>("stripEtaAssociationDistance", ps2);
        }
        temp2.addParameter<std::string>("name", "s");
        temp2.addParameter<std::string>("plugin", "RecoTauPiZeroStripPlugin3");
        temp2.addParameter<double>("minGammaEtStripAdd", 1.0);
        temp2.addParameter<double>("minGammaEtStripSeed", 1.0);
        temp2.addParameter<int>("maxStripBuildIterations", -1);
        temp2.addParameter<bool>("updateStripAfterEachDaughter", false);
        temp2.addParameter<bool>("makeCombinatoricStrips", false);
        temp2.addParameter<bool>("applyElecTrackQcuts", false);
        temp2.addParameter<std::vector<int>>("stripCandidatesParticleIds", {
          2,
          4,
        });
        temp2.addParameter<double>("minStripEt", 1.0);
        temp2.addParameter<int>("verbosity", 0);
        builders_vector.push_back(temp2);
      }
      desc.addVPSet("builders", vpsd_builders, builders_vector);
    }
    descriptions.add("ak4PFJetsLegacyHPSPiZeros", desc);
    descriptions.add("ak4PFJetsLegacyHPSPiZerosBoosted", desc); // this one is generated in configs with a strange procedure
    // RecoTauTag/Configuration/python/boostedHPSPFTaus_cfi.py
    //    process.PATTauSequenceBoosted = cloneProcessingSnippet(process,process.PATTauSequence, "Boosted", addToTask = True)
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RecoTauPiZeroProducer);
