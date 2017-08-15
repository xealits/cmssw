#include "RecoTauTag/RecoTau/interface/RecoTauBuilderPlugins.h"
#include "RecoTauTag/RecoTau/interface/RecoTauCommonUtilities.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

namespace reco { namespace tau {

// Update our reference to the PFCandidates & PVs
void RecoTauBuilderPlugin::beginEvent() {
  vertexAssociator_.setEvent(*evt());
  evt()->getByToken(pfCand_token, pfCands_);
}

}}  // end namespace reco::tau

#include "FWCore/Framework/interface/MakerMacros.h"
EDM_REGISTER_PLUGINFACTORY(RecoTauBuilderPluginFactory,
                           "RecoTauBuilderPluginFactory");
EDM_REGISTER_PLUGINFACTORY(RecoTauModifierPluginFactory,
                           "RecoTauModifierPluginFactory");
EDM_REGISTER_PLUGINFACTORY(RecoTauCleanerPluginFactory,
                           "RecoTauCleanerPluginFactory");
EDM_REGISTER_PLUGINFACTORY(RecoBaseTauBuilderPluginFactory,
                           "RecoBaseTauBuilderPluginFactory");
EDM_REGISTER_PLUGINFACTORY(RecoBaseTauModifierPluginFactory,
                           "RecoBaseTauModifierPluginFactory");
EDM_REGISTER_PLUGINFACTORY(RecoBaseTauCleanerPluginFactory,
                           "RecoBaseTauCleanerPluginFactory");
