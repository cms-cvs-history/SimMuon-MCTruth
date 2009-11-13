#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "SimMuon/MCTruth/test/testMyAssociator.h"
#include "SimMuon/MCTruth/interface/MuonAssociatorByHits.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

testMyAssociator::testMyAssociator(const edm::ParameterSet& parset) :
  tracksTag(parset.getParameter< edm::InputTag >("tracksTag")),
  tpTag(parset.getParameter< edm::InputTag >("tpTag")),
  parset_(parset)
{
  LogTrace("testMyAssociator") << "constructing  testMyAssociator" << parset_.dump();
}

testMyAssociator::~testMyAssociator() {
}

void testMyAssociator::beginJob() {

  associatorByHits = new MuonAssociatorByHits::MuonAssociatorByHits(parset_);
}

void testMyAssociator::analyze(const edm::Event& event, const edm::EventSetup& setup)
{  
  edm::Handle<edm::View<reco::Track> > trackCollectionH;
  LogTrace("testMyAssociator") << "getting reco::Track collection "<<tracksTag;
  event.getByLabel(tracksTag,trackCollectionH);

  const edm::View<reco::Track>  trackCollection = *(trackCollectionH.product()); 
  LogTrace("testMyAssociator") << "...size = "<<trackCollection.size();

  edm::Handle<TrackingParticleCollection>  TPCollectionH ;
  LogTrace("testMyAssociator") << "getting TrackingParticle collection "<<tpTag;
  event.getByLabel(tpTag,TPCollectionH);

  const TrackingParticleCollection tPC   = *(TPCollectionH.product());
  LogTrace("testMyAssociator") << "...size = "<<tPC.size();

  edm::LogVerbatim("testMyAssociator") <<"\n === Event ID = "<< event.id()<<" ===";
    
  //RECOTOSIM 
  edm::LogVerbatim("testMyAssociator") 
    << "\n                      ****************** Reco To Sim ****************** ";

  reco::RecoToSimCollection recSimColl = 
    associatorByHits->associateRecoToSim(trackCollectionH,TPCollectionH,&event,&setup);

  edm::LogVerbatim("testMyAssociator")
    << "\n There are " << trackCollection.size() << " reco::Tracks "<< "("<<recSimColl.size()<<" matched) \n";

  for(edm::View<reco::Track>::size_type i=0; i<trackCollection.size(); ++i) {
    edm::RefToBase<reco::Track> track(trackCollectionH, i);
    
    if(recSimColl.find(track) != recSimColl.end()) {
      std::vector<std::pair<TrackingParticleRef, double> > recSimAsso = recSimColl[track];
      
      for (std::vector<std::pair<TrackingParticleRef, double> >::const_iterator IT = recSimAsso.begin(); 
	   IT != recSimAsso.end(); ++IT) {
	TrackingParticleRef trpart = IT->first;
	double purity = IT->second;
	edm::LogVerbatim("testMyAssociator") <<"reco::Track #" << int(i) << " with pt = " << track->pt()
					     << " associated to TrackingParticle #" <<trpart.key()
					     << " (pt = " << trpart->pt() << ") with Quality = " << purity;
      }
    } else {
      edm::LogVerbatim("testMyAssociator") << "reco::Track #" << int(i) << " with pt=" << track->pt()
					   << " NOT associated to any TrackingParticle" << "\n";		   
    } 
  }
  
  //SIMTORECO
  edm::LogVerbatim("testMyAssociator") 
    << "\n                      ****************** Sim To Reco ****************** ";

  reco::SimToRecoCollection simRecColl =
    associatorByHits->associateSimToReco(trackCollectionH,TPCollectionH,&event,&setup);
  
  edm::LogVerbatim("testMyAssociator")
    << "\n There are " << tPC.size() << " TrackingParticles "<<"("<<simRecColl.size()<<" matched) \n";

  bool any_trackingParticle_matched = false;

  for (TrackingParticleCollection::size_type i=0; i<tPC.size(); i++) {
    TrackingParticleRef trpart(TPCollectionH, i);
    
    std::vector<std::pair<edm::RefToBase<reco::Track>, double> > simRecAsso;
    if(simRecColl.find(trpart) != simRecColl.end()) { 
      simRecAsso = (std::vector<std::pair<edm::RefToBase<reco::Track>, double> >) simRecColl[trpart];
      
      for (std::vector<std::pair<edm::RefToBase<reco::Track>, double> >::const_iterator IT = simRecAsso.begin(); 
	   IT != simRecAsso.end(); ++IT) {
	edm::RefToBase<reco::Track> track = IT->first;
	double quality = IT->second;
	any_trackingParticle_matched = true;

	// find the purity from RecoToSim association (set purity = -1 for unmatched recoToSim)
	double purity = -1.;
	if(recSimColl.find(track) != recSimColl.end()) {
	  std::vector<std::pair<TrackingParticleRef, double> > recSimAsso = recSimColl[track];
	  for (std::vector<std::pair<TrackingParticleRef, double> >::const_iterator ITS = recSimAsso.begin(); 
	       ITS != recSimAsso.end(); ++ITS) {
	    TrackingParticleRef tp = ITS->first;
	    if (tp == trpart) purity = ITS->second;
	  }
	}

	edm::LogVerbatim("testMyAssociator") <<"TrackingParticle #" << int(i)<< " with pt = " << trpart->pt()
					     << " associated to reco::Track #" <<track.key()
					     << " (pt = " << track->pt() << ") with Quality = " << quality
	                                     << " and Purity = "<< purity;
      }
    } 
    
    //    else 
    //      {
    //	LogTrace("testReader") << "TrackingParticle #" << int(i)<< " with pt = " << trpart->pt() 
    //			       << " NOT associated to any reco::Track" ;
    //      }    
  }
  
  if (!any_trackingParticle_matched) {
    edm::LogVerbatim("testMyAssociator") << "NO TrackingParticle associated to ANY input reco::Track !" << "\n";
  }
  
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(testMyAssociator);


