#include "G4numiNtp.hh"

G4numiNtp::G4numiNtp()
{}

G4numiNtp::~G4numiNtp()
{}

//void G4numiNtp::SetTree(TTree &t, Int_t entry)
void G4numiNtp::SetTree(TTree &t)
{
  entry=0;

  l_run=t.GetLeaf("run");
  run=Int_t(l_run->GetValue(entry));
  l_evtno = t.GetLeaf("evtno");
  evtno= Int_t(l_evtno->GetValue(entry));
  l_Ndxdz = t.GetLeaf("Ndxdz");
  Ndxdz= l_Ndxdz->GetValue(entry);
  l_Ndydz = t.GetLeaf("Ndydz");
  Ndydz= l_Ndydz->GetValue(entry);
  l_Npz = t.GetLeaf("Npz");
  Npz= l_Npz->GetValue(entry);
  l_Nenergy = t.GetLeaf("Nenergy");
  Nenergy= l_Nenergy->GetValue(entry);
  l_Norig = t.GetLeaf("Norig");
  Norig= Int_t(l_Norig->GetValue(entry));
  l_Ndecay = t.GetLeaf("Ndecay");
  Ndecay= Int_t(l_Ndecay->GetValue(entry));
  l_Ntype = t.GetLeaf("Ntype");
  Ntype= Int_t(l_Ntype->GetValue(entry));
  l_Vx = t.GetLeaf("Vx");
  Vx= l_Vx->GetValue(entry);
  l_Vy = t.GetLeaf("Vy");
  Vy= l_Vy->GetValue(entry);
  l_Vz = t.GetLeaf("Vz");
  Vz= l_Vz->GetValue(entry);
  l_pdpx = t.GetLeaf("pdPx");
  pdpx= l_pdpx->GetValue(entry);
  l_pdpy = t.GetLeaf("pdPy");
  pdpy= l_pdpy->GetValue(entry);
  l_pdpz = t.GetLeaf("pdPz");
  pdpz= l_pdpz->GetValue(entry);
  l_ppdxdz = t.GetLeaf("ppdxdz");
  ppdxdz= l_ppdxdz->GetValue(entry);
  l_ppdydz = t.GetLeaf("ppdydz");
  ppdydz= l_ppdydz->GetValue(entry);
  l_pppz = t.GetLeaf("pppz");
  pppz= l_pppz->GetValue(entry);
  l_ppenergy = t.GetLeaf("ppenergy");
  ppenergy= l_ppenergy->GetValue(entry);
  l_ppmedium = t.GetLeaf("ppmedium");
  ppmedium= Int_t(l_ppmedium->GetValue(entry));
  l_ptype = t.GetLeaf("ptype");
  ptype= Int_t(l_ptype->GetValue(entry));
  l_ppvx = t.GetLeaf("ppvx");
  ppvx= l_ppvx->GetValue(entry);
  l_ppvy = t.GetLeaf("ppvy");
  ppvy= l_ppvy->GetValue(entry);
  l_ppvz = t.GetLeaf("ppvz");
  ppvz= l_ppvz->GetValue(entry);
  l_muparpx = t.GetLeaf("muparpx");
  muparpx= l_muparpx->GetValue(entry);
  l_muparpy = t.GetLeaf("muparpy");
  muparpy= l_muparpy->GetValue(entry);
  l_muparpz = t.GetLeaf("muparpz");
  muparpz= l_muparpz->GetValue(entry);
  l_Necm = t.GetLeaf("Necm");
  Necm= l_Necm->GetValue(entry);
  l_Nimpwt = t.GetLeaf("Nimpwt");
  Nimpwt= l_Nimpwt->GetValue(entry);
  l_xpoint = t.GetLeaf("xpoint");
  xpoint= l_xpoint->GetValue(entry);
  l_ypoint = t.GetLeaf("ypoint");
  ypoint= l_ypoint->GetValue(entry);
  l_zpoint = t.GetLeaf("zpoint");
  zpoint= l_zpoint->GetValue(entry);
  l_tvx = t.GetLeaf("tvx");
  tvx= l_tvx->GetValue(entry);
  l_tvy = t.GetLeaf("tvy");
  tvy= l_tvy->GetValue(entry);
  l_tvz = t.GetLeaf("tvz");
  tvz= l_tvz->GetValue(entry);
  l_tpx = t.GetLeaf("tpx");
  tpx= l_tpx->GetValue(entry);
  l_tpy = t.GetLeaf("tpy");
  tpy= l_tpy->GetValue(entry);
  l_tpz = t.GetLeaf("tpz");
  tpz= l_tpz->GetValue(entry);
  l_tptype = t.GetLeaf("tptype");
  tptype= Int_t(l_tptype->GetValue(entry));
  l_tgen = t.GetLeaf("tgen");
  tgen= Int_t(l_tgen->GetValue(entry));
  /*l_tgptype = t.GetLeaf("tgptype");
  tgptype= l_tgptype->GetValue(entry);
  l_tgppx = t.GetLeaf("tgppx");
  tgppx= l_tgppx->GetValue(entry);
  l_tgppy = t.GetLeaf("tgppy");
  tgppy= l_tgppy->GetValue(entry);
  l_tgppz = t.GetLeaf("tgppz");
  tgppz= l_tgppz->GetValue(entry);
  l_tprivx = t.GetLeaf("tprivx");
  tprivx= l_tprivx->GetValue(entry);
  l_tprivy = t.GetLeaf("tprivy");
  tprivy= l_tprivy->GetValue(entry);
  l_tprivz = t.GetLeaf("tprivz");
  tprivz= l_tprivz->GetValue(entry);
  l_beamx = t.GetLeaf("beamx");
  beamx= l_beamx->GetValue(entry);
  l_beamy = t.GetLeaf("beamy");
  beamy= l_beamy->GetValue(entry);
  l_beamz = t.GetLeaf("beamz");
  beamz= l_beamz->GetValue(entry);
  l_beampx = t.GetLeaf("beampx");
  beampx= l_beampx->GetValue(entry);
  l_beampy = t.GetLeaf("beampy");
  beampy= l_beampy->GetValue(entry);
  l_beampz = t.GetLeaf("beampz");
  beampz= l_beampz->GetValue(entry);
  */

  //index of array with 11 near values and 2 far values
  entry=0;
  
  l_NdxdzNear=t.GetLeaf("NdxdzNear");
  NdxdzNear=l_NdxdzNear->GetValue(entry);
  l_NdydzNear=t.GetLeaf("NdydzNear");
  NdydzNear=l_NdydzNear->GetValue(entry);
  l_NenergyN=t.GetLeaf("NenergyN");
  NenergyN=l_NenergyN->GetValue(entry);  
  l_NWtNear = t.GetLeaf("NWtNear");
  NWtNear= l_NWtNear->GetValue(entry);
  NWtNear_one=1.0;
  l_NdxdzFar = t.GetLeaf("NdxdzFar");
  NdxdzFar= l_NdxdzFar->GetValue(entry);
  l_NdydzFar = t.GetLeaf("NdydzFar");
  NdydzFar= l_NdydzFar->GetValue(entry);
  l_NenergyF = t.GetLeaf("NenergyF");
  NenergyF= l_NenergyF->GetValue(entry);
  l_NWtFar = t.GetLeaf("NWtFar");
  NWtFar= l_NWtFar->GetValue(entry);
  

  //l_ = t.GetLeaf("");
  //= l_->GetValue(entry);

}
