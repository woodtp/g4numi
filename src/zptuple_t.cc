//
// zpdata_t.cc
//
//  Jasmine, July 2005
//  This is a class that defines the hadmmtuple_t object that is used to store the g4numi output
//   in a root tree.
//------------------

#include "zptuple_t.hh"

ClassImp(zptuple_t)

  zptuple_t::zptuple_t():
  run(-1),
  evtno(-10),
  xposatz(-1.e4),
  yposatz(-1.e4),
  zposatz(-1.e4),
  xmomatz(-1.e4),
  ymomatz(-1.e4),
  zmomatz(-1.e4),
  matilen(-1.e4),
  field(-1.e4),
  pathlength(-1.e4),
  ptypeatz(-10),
  pidtype(-10),
  zpoint(-10000)
{


}
  
    

  //-----------------------------------------------------------------------------------
zptuple_t::~zptuple_t()
{
  // zptuple_t destructor
}

