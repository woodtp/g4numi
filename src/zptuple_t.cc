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
  matilen(-1.e4),
  field(-1.e4),
  pidtype(-10),
  xposatz(-10000),
  yposatz(-10000),
  zposatz(-10000),
  xmomatz(-10000),
  ymomatz(-10000),
  zmomatz(-10000),
  ptypeatz(-10000)
{
  
}
  
    

  //-----------------------------------------------------------------------------------
zptuple_t::~zptuple_t()
{
  // zptuple_t destructor
}

