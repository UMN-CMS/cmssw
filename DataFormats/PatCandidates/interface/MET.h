//
// $Id: MET.h,v 1.10 2008/04/04 18:22:08 srappocc Exp $
//

#ifndef DataFormats_PatCandidates_MET_h
#define DataFormats_PatCandidates_MET_h

/**
  \class    pat::MET MET.h "DataFormats/PatCandidates/interface/MET.h"
  \brief    Analysis-level MET class

   MET implements an analysis-level missing energy class as a 4-vector
   within the 'pat' namespace.

  \author   Steven Lowette
  \version  $Id: MET.h,v 1.10 2008/04/04 18:22:08 srappocc Exp $
*/


#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"

namespace pat {


  typedef reco::CaloMET METType;


  class MET : public PATObject<METType> {

    public:

      MET();
      MET(const METType & aMET);
      MET(const edm::RefToBase<METType> & aMETRef);
      MET(const edm::Ptr<METType> & aMETRef);
      virtual ~MET();

      virtual MET * clone() const { return new MET(*this); }
      const reco::GenMET * genMET() const;

      void setGenMET(const reco::GenMET & gm);

    protected:

      std::vector<reco::GenMET> genMET_;

  };


}

#endif
