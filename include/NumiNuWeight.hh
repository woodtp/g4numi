// $Id: NumiNuWeight.hh,v 1.1.4.1 2010/08/19 19:50:54 minervacvs Exp $

#ifndef NUMINUWEIGHT_H
#define NUMINUWEIGHT_H 1

#include "Rtypes.h"
#include "data_t.hh"
#include <vector>

class NumiNuWeight
{
  public:
    NumiNuWeight();
    ~NumiNuWeight();
 
    double GetWeight(const data_t* nudata, 
		     const std::vector<double> xdet,
		     double& nu_wght, double& nu_energy);
 
};
 
#endif
