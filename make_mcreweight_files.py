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

# add 2d xfpt histograms from each file
# needs: list of files, histo_name, returns xf,pt histo
def sum_histograms(file_list, histo_name):
    hsum=TH2D()
    hsum.SetDirectory(0)
    cntr=0
    for f in file_list:
        h=get_histo(histo_name,f)
#        print "working on file %i",cntr
        if cntr==0:
            hsum=h
        else:
            hsum.Add(h)
        cntr=cntr+1
    return hsum



def yield_to_xsec(yield_histo,mom,mproj,mproduced,nincident,rho,A,dx):
    """ copy histo and convert from yield to cross section, returns xf,pt histo
    needs to know details: projectile momentum, projectile mass, produced particle mass,number of incident projectiles, target density, A and thickness.  The returned histo has a temporary name which should be reassigned."""

    NA=6.022e23 # atoms/mol
    cm2_per_mb=1e-27 # conversion factor
    sigma_factor=1.0/(rho/A * NA * dx * cm2_per_mb) # mb

    h=yield_histo.Clone("xsec_temp")
    h.SetDirectory(0)
    h.Sumw2()
    for ix in range(1,h.GetNbinsX()+1):
        for ipt in range(1,h.GetNbinsY()+1):
            ptval=h.GetYaxis().GetBinCenter(ipt)
            dpt=h.GetYaxis().GetBinWidth(ipt)
            ptup=ptval+dpt/2.0
            ptdown=ptval-dpt/2.0
            dpt2=abs(ptup**2 - ptdown**2)
            xfval=h.GetXaxis().GetBinCenter(ix)
            dxf=h.GetXaxis().GetBinWidth(ix)
            E = get_E_lab(mom,mproj,xfval,ptval,mproduced)
            pzup=get_pz_lab(mom,mproj,xfval+dxf/2.0,ptval,mproduced)
            pzdown=get_pz_lab(mom,mproj,xfval-dxf/2.0,ptval,mproduced)
            dpz=abs(pzup-pzdown)
            dp3=cmath.pi*dpz*dpt2
            y=h.GetBinContent(ix,ipt)# the yield for this bin
            err_y=h.GetBinError(ix,ipt)# the error on the yield
            
            invxs=E/dp3 * (y/nincident)*sigma_factor
            err_invxs=E/dp3 * (err_y/nincident)*sigma_factor
            h.SetBinContent(ix,ipt,invxs)
            h.SetBinError(ix,ipt,err_invxs)
    return h


def get_E_lab(mom,mproj,xf,pt,m):
    """ get the produced particle energy in the lab frame 
    based on projectile momentum and mass, 
    and xf, pt, and m of produced particle """
    mtarg=(0.938272+0.939565)/2.0 # average of proton and neutron
    beam_ener = (mom**2+mproj**2)**0.5
    beta= mom/(beam_ener+mtarg)
    gamma = 1.0/(1-beta**2)**0.5
    Ecm = (mproj**2 + mtarg**2 + 2*beam_ener*mtarg)**0.5
    plcm=xf*Ecm/2.0
    Ecm=(plcm**2 + pt**2 + m**2)**0.5
    E=gamma*(Ecm+beta*plcm)
    return E

def get_pz_lab(mom,mproj,xf,pt,m):
    """ get the produced particle pz in the lab frame
    based on projectile momentum and mass, 
    and xf, pt, and m of produced particle """

    mtarg=(0.938272+0.939565)/2.0 # average of proton and neutron
    beam_ener = (mom**2+mproj**2)**0.5
    beta= mom/(beam_ener+mtarg)
    gamma = 1.0/(1-beta**2)**0.5
    Ecm = (mproj**2 + mtarg**2 + 2*beam_ener*mtarg)**0.5
    plcm=xf*Ecm/2.0
    Ecm=(plcm**2 + pt**2 + m**2)**0.5
    pz=gamma*(beta*Ecm+plcm)
    return pz

# write xsec histo to a file, needs histo, histo name, file name
# must update the file, not recreate it (makes the program more modular)

def sum_qeinfo(filelist):
    tot=0
    el_like=0
    qe_like=0
    frag_like=0
    meson_prod=0
    nincident=len(filelist)*26.2e6
    pip=0.0
    pim=0.0
    for name in filelist:
        f=open(name)
        line1=f.readline()
        line2=f.readline()
        line3=f.readline()
        line4=f.readline()
        line5=f.readline()
    
        v=line3.split()
        tot=tot+int(v[0])
        el_like=el_like+int(v[2])
        qe_like=qe_like+int(v[3])
        frag_like=frag_like+int(v[4])
        meson_prod=meson_prod+int(v[5])

        w=line4.split()
        pip=pip+float(w[6])
        x=line5.split()
        pim=pim+float(x[6])
    rho=1.78
    dx=0.7
    A=12.01
    NA=6.022e23
    cm2_per_mb=1e-27
    sigma_factor=1.0/(rho/A * NA * dx * cm2_per_mb) # mb
    sigma_tot = tot/nincident*sigma_factor
    print "Total reaction cross-section: %3.3f mb"%(sigma_tot)
    print "  EL    QEL    FRAG   INEL"
    print "%3.3f%% %3.3f%% %3.3f%% %3.3f%%"%(el_like/float(tot)*100, qe_like/float(tot)*100, frag_like/float(tot)*100, meson_prod/float(tot)*100)
    
    print "Inelastic (meson production) cross-section: %3.3f mb"%(sigma_tot*meson_prod/float(tot))
    print "EL + QEL + FRAG cross-section: %3.3f mb"%(sigma_tot*(el_like+qe_like+frag_like)/float(tot))
    
    print "average pi+ multiplicity per production event: %3.3f"%(pip/float(len(filelist)))
    print "average pi- multiplicity per production event: %3.3f"%(pim/float(len(filelist)))
    return meson_prod

#h=get_h1("pip/xF197_pip","/minerva/data/users/kordosky/hp/FTFP_BERT/histos/Yields_FTFP_BERT_pC0158GeV_mc07.root")
#h.Print()

#particles=["neu","kshort","klong"]
#masses=[0.939565,0.497614,0.497614]
particles=["pip","pim","kap","kam","prt","neu","kshort","klong"]
masses=[0.13957,0.13957,0.493667,0.493667,0.938272,0.939565,0.497614,0.497614]

momentums=[120,110,100,90,80,70,60,50,40,31,20,12]
#momentums=[120,110]
for part,mproduced in zip(particles,masses):
    outname="%s_FTFP_BERT.root"%(part)
    fout=TFile(outname,"RECREATE")
    for mom in momentums:
        print "Working on ", part, " production  with beam momentum ", mom, "GeV/c"        
        indir="/minerva/data/users/kordosky/hp/FTFP_BERT/histos/"
        file_pat="Yields_FTFP_BERT_pC%04iGeV_mc*.root"%(mom)
        print indir+file_pat
        filelist=glob.glob(indir+file_pat)
        pot_per_file=26.2e6
        rho=1.78 #g/cc
        A=12.01 #g/mol
        dx=0.7 #cm
        mproj=0.938272
        nfiles=len(filelist)
        nincident=nfiles*pot_per_file
        hin_name="xFpT_%s"%(part)
        hout_name="xF_pT_%iGeV"%(mom)
        print hin_name,hout_name
        hyield=sum_histograms(filelist, hin_name)
        hx=yield_to_xsec(hyield,mom,mproj,mproduced,nincident,rho,A,dx)
        hx.SetName(hout_name)
        hx.SetDirectory(fout)
        qe_pat="QEInfo_FTFP_BERT_pC%04iGeV_mc*.txt"%(mom)
        if part=="neu":
            qelist=glob.glob(indir+qe_pat)
            nprod=sum_qeinfo(qelist)
            neut_in="dndxf_neu_prod_cut"
            neut_out="dndxf_%iGeV"%(mom)
            hneut=sum_histograms(filelist,neut_in)
            hneut.Sumw2()
            hneut.Scale(1.0/(nprod*hneut.GetBinWidth(1)))
            hneut.SetName(neut_out)
            hneut.SetDirectory(fout)
            neut_in2="dndxf_neu_prod"
            neut_out2="dndxf_nocut_%iGeV"%(mom)
            hneut2=sum_histograms(filelist,neut_in2)
            hneut2.Sumw2()
            hneut2.Scale(1.0/(nprod*hneut2.GetBinWidth(1)))
            hneut2.SetName(neut_out2)
            hneut2.SetDirectory(fout)


        fout.ls()
        fout.Write()

    fout.Close()

#h.Draw("colz")
#hx.Draw("colz")
#TPython.Prompt()
