/*
 * =============================================================================
 *       Filename:  RecoTauGenericImpactParameterSignificancePlugin.cc
 *
 *    Description:  Add the IP significance of the lead track w.r.t to the PV.
 *                  to a PFTau.
 *        Created:  10/31/2010 13:32:14
 *
 *         Authors:  Evan K. Friis (UC Davis), evan.klose.friis@cern.ch,
 *                   Simone Gennai, Ludovic Houchu
 *
 * =============================================================================
 */

#include "RecoTauTag/RecoTau/interface/RecoTauBuilderPlugins.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

namespace reco { namespace tau {

template<class TauType>
class RecoTauGenericImpactParameterSignificancePlugin : public RecoTauModifierPlugin<TauType> {
  public:
    explicit RecoTauGenericImpactParameterSignificancePlugin(
						      const edm::ParameterSet& pset,edm::ConsumesCollector &&iC);
    virtual ~RecoTauGenericImpactParameterSignificancePlugin() {}
    void operator()(TauType& tau) const override;
    virtual void beginEvent() override;
  private:
    RecoTauVertexAssociator vertexAssociator_;
    const TransientTrackBuilder *builder_;
};

template<class TauType>
RecoTauGenericImpactParameterSignificancePlugin<TauType>
::RecoTauGenericImpactParameterSignificancePlugin(const edm::ParameterSet& pset,edm::ConsumesCollector &&iC)
  :RecoTauModifierPlugin<TauType>(pset,std::move(iC)),
   vertexAssociator_(pset.getParameter<edm::ParameterSet>("qualityCuts"),std::move(iC)){}

template<class TauType>
void RecoTauGenericImpactParameterSignificancePlugin<TauType>::beginEvent() {
  vertexAssociator_.setEvent(*this->evt());
  // Get tranisent track builder.
  edm::ESHandle<TransientTrackBuilder> myTransientTrackBuilder;
  this->evtSetup()->get<TransientTrackRecord>().get("TransientTrackBuilder",
                                              myTransientTrackBuilder);
  builder_= myTransientTrackBuilder.product();
}

namespace 
{
  // JAN - FIXME - this method is needed multiple times

  inline const reco::TrackBaseRef getTrack(const Candidate& cand)
  {
    const PFCandidate* pfCandPtr = dynamic_cast<const PFCandidate*>(&cand);
    if (pfCandPtr) {
      if      ( pfCandPtr->trackRef().isNonnull()    ) return reco::TrackBaseRef(pfCandPtr->trackRef());
      // else if ( pfCandPtr->gsfTrackRef().isNonnull() ) return reco::TrackBaseRef(pfCandPtr->gsfTrackRef());
      else return reco::TrackBaseRef();
    }
    // JAN - FIXME: Add method for miniAOD PackedCandidate
    return reco::TrackBaseRef();
  }
}

template<class TauType>
void RecoTauGenericImpactParameterSignificancePlugin<TauType>::operator()(TauType& tau) const {
  // Get the transient lead track
  if (tau.leadPFChargedHadrCand().isNonnull()) {
    TrackBaseRef leadTrack = getTrack(*tau.leadPFChargedHadrCand());
    if (leadTrack.isNonnull()) {
      const TransientTrack track = builder_->build(*leadTrack);
      GlobalVector direction(tau.jetRef()->px(), tau.jetRef()->py(),
                             tau.jetRef()->pz());
      VertexRef pv = vertexAssociator_.associatedVertex(tau);
      // Compute the significance
      std::pair<bool,Measurement1D> ipsig =
          IPTools::signedImpactParameter3D(track, direction, *pv);
      if (ipsig.first)
        tau.setleadPFChargedHadrCandsignedSipt(ipsig.second.significance());
    }
  }
}

template class RecoTauGenericImpactParameterSignificancePlugin<reco::PFTau>;
typedef RecoTauGenericImpactParameterSignificancePlugin<reco::PFTau> RecoTauImpactParameterSignificancePlugin;

template class RecoTauGenericImpactParameterSignificancePlugin<reco::PFBaseTau>;
typedef RecoTauGenericImpactParameterSignificancePlugin<reco::PFBaseTau> RecoBaseTauImpactParameterSignificancePlugin;

}} // end namespace reco::tau
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_EDM_PLUGIN(RecoTauModifierPluginFactory,
    reco::tau::RecoTauImpactParameterSignificancePlugin,
    "RecoTauImpactParameterSignificancePlugin");
DEFINE_EDM_PLUGIN(RecoBaseTauModifierPluginFactory,
    reco::tau::RecoBaseTauImpactParameterSignificancePlugin,
    "RecoBaseTauImpactParameterSignificancePlugin");
