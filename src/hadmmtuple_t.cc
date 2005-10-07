//
// data_t.cc
//
//  ADM, July 2005
//  This is a class that defines the hadmmtuple_t object that is used to store the g4numi output
//   in a root tree.
//------------------

#include "hadmmtuple_t.hh"

ClassImp(hadmmtuple_t)

hadmmtuple_t::hadmmtuple_t():
  run(-1),       
  evtno(-1), 
  mtgthpos(-1.e4),
  mtgtvpos(-1.e4),
  mtgthsig(-1.e4),
  mtgtvsig(-1.e4),
  ptype(-1),
  hmmenergy(-1.e4),
  hmmxpos(-1.e4),
  hmmypos(-1.e4),
  hmmzpos(-1.e4),
  hmmpx(-1.e4),
  hmmpy(-1.e4),
  hmmpz(-1.e4)
{
  //hadmmtuple_t constructor
}

//-----------------------------------------------------------------------------------
hadmmtuple_t::~hadmmtuple_t()
{
  // hadmmtuple_t destructor
}

