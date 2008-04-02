#ifndef testMyAssociator_h
#define testMyAssociator_h

#include <memory>
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include <iostream>
#include <string>
#include <map>

class MuonAssociatorByHits;

class testMyAssociator : public edm::EDAnalyzer {
  
 public:
  testMyAssociator(const edm::ParameterSet&);
  virtual ~testMyAssociator();
  virtual void beginJob( const edm::EventSetup& );  
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  
 private:
  edm::InputTag tracksTag, tpTag, simtracksTag;
  const edm::ParameterSet  parset_;

  MuonAssociatorByHits * associatorByHits;
  
};

#endif
