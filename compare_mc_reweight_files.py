#! python
from ROOT import *
import glob
import cmath

# pull a histogram from a file
def get_histo(histname,filename):
    f = TFile(filename)
    h = f.Get(histname)
    h.SetDirectory(0)
    return h


momentums=[120,110,100,90,80,70,60,50,40,31,20,12]

hchi2=TH1F("hchi2",";#chi^{2}/ndf",100,0,10)
hchi2prob=TH1F("hchi2prob",";#chi^{2} probability",50,0,1)
hratio=TH1F("hratio",";#chi^{2} probability",40,0.98,1.02)

particles=["pip","pim","kap","kam","prt"]
dirB="/minerva/app/users/kordosky/cmtuser/Minerva_v10r6p9/MParamFiles/data/Reweight/new_RW_files/"
for part in particles:
    fA="%s_FTFP_BERT.root"%(part)
    fB=dirB+fA
    for mom in momentums:
        hname="xF_pT_%iGeV"%(mom)
        hA=get_histo(hname,fA)
        hB=get_histo(hname,fB)
        
#        chi2=hA.Chi2Test(hB,"WW")
#        ndf=(hA.GetNbinsX()*hA.GetNbinsY())
#        normchi2=chi2/ndf
#        hchi2.Fill(normchi2)
#        hchi2prob.Fill(TMath.Prob(chi2,ndf))
#        print "%s %s has chi2
        ratio=hB.Integral()/hA.Integral()
        hratio.Fill(ratio)
        ohno=""
        if ratio>1.01 or ratio<0.99:
            ohno="OH NO!!!"
        print "%s %s  %s  has ratio = %f"%(ohno,fA, hname,ratio)
#        hr = hB.Clone(hA.GetName()+"_rat")
#        hr.Divide(hA)
        
#        hr.SetMaximum(1.05)
#        hr.SetMinimum(0.95)
#        hr.Draw("colz")
#        raw_input("Press Enter to continue...")
        
c1 = TCanvas()
hratio.Draw()
TPython.Prompt()
