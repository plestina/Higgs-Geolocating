#! /usr/bin/env python
import sys
import os
import re
import math
from ROOT import *
import ROOT
from array import array
import lib.util.MiscTools as misc
import optparse

#import errno


def parseOptions():
    pass
    #usage = ('usage: %prog [options] \n'
             #+ '%prog -h for help')
    #parser = optparse.OptionParser(usage)
    ##parser.add_option('-p', action='store_true', dest='plot', default=False ,help='do plot')
    ##parser.add_option('-p', action='store_true', dest='plot', default=False ,help='do plot')
    #parser.add_option('-d', '--dir', dest='run_on_dir', type='string', default="",    help='Path with 2D templates')
    ### store options and arguments as global variables
    #global opt, args
    #(opt, args) = parser.parse_args()
def refurbishedHist(hist,histName,norm=-1):

    #nBinsX = hist.GetXaxis().GetNbins()
    #nBinsY = hist.GetYaxis().GetNbins()
    
    #hist2D = ROOT.TH2F()
    hist2D = hist
    hist2D.SetTitle(histName)
    return hist2D

def unfoldedHist(hist,histName,norm=-1):

    nBinsX = hist.GetXaxis().GetNbins()
    nBinsY = hist.GetYaxis().GetNbins()
    totalNbins = nBinsX*nBinsY
    print ">>>>>> Bins in 1D hist "+str(histName)+": ",totalNbins
    
    hist1D = ROOT.TH1F(histName,histName,totalNbins,0,totalNbins)

    i=0
    for j in xrange(nBinsX):
	for k in xrange(nBinsY):
	    
	    binContent = hist.GetBinContent(j+1,k+1)
	    hist1D.SetBinContent(i+1,binContent)
	    #print hist1D.GetBinContent(i+1)
	    i+=1
    
    if norm > 0: hist1D.Scale(norm/hist1D.Integral("width"))

    #figs = self.outputDir+"/figs"
    figs = ""
    canv = ROOT.TCanvas("c","c",700,700)
    canv.cd()
    hist1D.SetFillColor(kBlack)
    if histName.startswith("ggH"):
	hist1D.SetFillColor(kOrange+10)
    if histName.startswith("bkg2d_qqzz"):
	hist1D.SetFillColor(kAzure-9)
    if histName.startswith("bkg2d_ggzz"):
	hist1D.SetFillColor(kAzure-8)
    if histName.startswith("bkg2d_zjets"):
	hist1D.SetFillColor(kGreen-5)
    hist1D.Draw("HIST")
    #canv.SaveAs(figs+histName+".eps")
    #canv.SaveAs(figs+histName+".png")

		
    return hist1D

    
def unfoldedHist_Flat(hist,histName,norm=-1):

    nBinsX = hist.GetXaxis().GetNbins()
    nBinsY = hist.GetYaxis().GetNbins()
    totalNbins = nBinsX*nBinsY
    print ">>>>>> Bins in 1D hist "+str(histName)+": ",totalNbins
    
    hist1D = ROOT.TH1F(histName,histName,totalNbins,0,totalNbins)

    i=0
    binContent = 1
    for j in xrange(nBinsX):
	for k in xrange(nBinsY):
	    
	    hist1D.SetBinContent(i+1,binContent)
	    #print hist1D.GetBinContent(i+1)
	    i+=1
    
    if norm > 0: hist1D.Scale(norm/hist1D.Integral("width"))

    #figs = self.outputDir+"/figs"
    figs = ""
    canv = ROOT.TCanvas("c","c",700,700)
    canv.cd()
    hist1D.SetFillColor(kBlack)
    if histName.startswith("ggH"):
	hist1D.SetFillColor(kOrange+10)
    if histName.startswith("bkg2d_qqzz"):
	hist1D.SetFillColor(kAzure-9)
    if histName.startswith("bkg2d_ggzz"):
	hist1D.SetFillColor(kAzure-8)
    if histName.startswith("bkg2d_zjets"):
	hist1D.SetFillColor(kGreen-5)
    hist1D.Draw("HIST")
    canv.SaveAs(figs+histName+".eps")
    canv.SaveAs(figs+histName+".png")

		
    return hist1D

        
    
    
def make_projected_hist(template_file, dir_2D ):
    
    print "------------- Projecting  templates to X"
    #dir_2D = '/afs/cern.ch/work/r/roko/Stat/CMSSW_611_JCP/src/HZZ4L_Combination/CombinationPy/CreateDatacards/Sandbox/templates_compare2Pedja/Templates2D_D0M_Dint13/'
    dir_2D_projected = dir_2D+"/projected/"
    misc.make_sure_path_exists(dir_2D_projected)
    #f=ROOT.TFile('/afs/cern.ch/work/r/roko/Stat/CMSSW_611_JCP/src/HZZ4L_Combination/CombinationPy/CreateDatacards/Sandbox/templates2D_D0M_Dint13/Dsignal_2e2mu.root')
    #f=ROOT.TFile('/afs/cern.ch/work/r/roko/Stat/CMSSW_611_JCP/src/HZZ4L_Combination/CombinationPy/CreateDatacards/Sandbox/Templates2D_D0Ph_Dint12/DSignal_4l.root')
    f=ROOT.TFile(dir_2D +"/"+template_file,"READ")
    print f.GetName()
    #f.cd()
    #ggH_proj = f.Get('ggH_shape').ProjectionX("ggH_shape_D0M",0,-1,"e")
    allowedTermNames = ['ggH',    'gg0Ph','gg0M' ,'ggInt_12P','ggInt_12N','ggInt_13P','ggInt_13N','ggInt_23P','ggInt_23N',
                        'qqZZ','ggZZ','zjets']
    #allowedTermNames = ['ggH',   'gg0M', 'ggInt_13P','ggInt_13N']
    new_names = {
    'gg0M' : 'gg0Ph',
    'ggInt_13P':'ggInt_12P',
    'ggInt_13N':'ggInt_12N',
    }
    fnew=ROOT.TFile(dir_2D_projected +"/"+template_file,"RECREATE")
    print fnew.GetName()
    projected = None
    #th2 = None
    found_th2 = False
    th2 = None
    for term in allowedTermNames:
        print "Searching for :", term
        try :
              th2_name = f.Get('{0}_shape'.format(term)).GetTitle()
        except ReferenceError:
              print "Doesn't exist : ", term
        else:
              print "----", th2_name 
              found_th2=True
              print "Found term for axis definiton: ", term
              th2 = f.Get('{0}_shape'.format(term))
              break
        #if found_th2 : break
              
          
    print "Axis x = {0} y = {1}".format(th2.GetXaxis().GetTitle(), th2.GetYaxis().GetTitle())
    for term in allowedTermNames:
        try:
            #th2 = f.Get('{0}_shape'.format(term))
            #projected = f.Get('{0}_shape'.format(term)).GetName()
            projected = f.Get('{0}_shape'.format(term)).ProjectionX("{0}_shape_D0M".format(term))

        except:
            print "@@@@ Shape {0}_shape doesn't exist.".format(term)
        else:
            print 'Working with {0}_shape'.format(term)
            fnew.cd()
            
            if 'Dbackground_ggZZ' in template_file: 
                    print "---------------> Changing name to : ggZZ_shape"
                    projected.SetNameTitle("ggZZ_shape","ggZZ_shape_projected")
            elif 'Dbackground_ZJets' in template_file: 
                    print "---------------> Changing name to : zjets_shape"
                    projected.SetNameTitle("zjets_shape","zjets_shape_projected")
                    
            else: 
                term_name = term
                if copy_13_templates_to_12 and term in new_names.keys():
                        term_name = new_names[term]
                projected.SetNameTitle('{0}_shape'.format(term_name),"{0}_shape_projected".format(term))
            print "New template name/title: {0}/{1}".format(projected.GetName(),projected.GetTitle())
            c1 =ROOT.TCanvas("Templates","Templates",800,800)
            c1.cd()
            ROOT.gPad.SetRightMargin(0.085)
            projected.Draw("hist")
            misc.make_sure_path_exists(dir_2D_projected+"/fig")
            for ext in [".png", ".pdf", ".root",".C"]: 
                c1.SaveAs(dir_2D_projected +"/fig/"+template_file+"_"+term+ext)            
            projected.Write()
            
    if not 'Dbackground' in template_file: 
                copy_tree(f, fnew,"factors")
    f.Close()
    fnew.Close()   
    
    
def make_unfolded_hist(template_file, dir_2D):
    
    print "------------- Unfolding  templates"
    #dir_2D = '/afs/cern.ch/work/r/roko/Stat/CMSSW_611_JCP/src/HZZ4L_Combination/CombinationPy/CreateDatacards/Sandbox/templates_compare2Pedja/Templates2D_D0M_Dint13/'
    dir_2D_unfolded = dir_2D+"/unfolded/"
    misc.make_sure_path_exists(dir_2D_unfolded)
    #f=ROOT.TFile('/afs/cern.ch/work/r/roko/Stat/CMSSW_611_JCP/src/HZZ4L_Combination/CombinationPy/CreateDatacards/Sandbox/templates2D_D0M_Dint13/Dsignal_2e2mu.root')
    #f=ROOT.TFile('/afs/cern.ch/work/r/roko/Stat/CMSSW_611_JCP/src/HZZ4L_Combination/CombinationPy/CreateDatacards/Sandbox/Templates2D_D0Ph_Dint12/DSignal_4l.root')
    f=ROOT.TFile(dir_2D +"/"+template_file,"READ")
    print f.GetName()
    #f.cd()
    #ggH_proj = f.Get('ggH_shape').ProjectionX("ggH_shape_D0M",0,-1,"e")
    allowedTermNames = ['ggH',	  'gg0Ph','gg0M' ,'ggInt_12P','ggInt_12N','ggInt_13P','ggInt_13N','ggInt_23P','ggInt_23N',
			'qqZZ','ggZZ','zjets']
    new_names = {
        'gg0M' : 'gg0Ph',
        'ggInt_13P':'ggInt_12P',
        'ggInt_13N':'ggInt_12N',
        }
			
			
    #allowedTermNames = ['ggH',	  'gg0M', 'ggInt_13P','ggInt_13N']
    fnew=ROOT.TFile(dir_2D_unfolded +"/"+template_file,"RECREATE")
    print fnew.GetName()
    unfolded = None
    #th2 = None
    found_th2 = False
    th2 = None
    for term in allowedTermNames:
	print "Searching for :", term
	try :
	      th2_name = f.Get('{0}_shape'.format(term)).GetTitle()
	except ReferenceError:
	      print "Doesn't exist : ", term
	else:
	      print "----", th2_name 
	      found_th2=True
	      print "Found term for axis definiton: ", term
	      th2 = f.Get('{0}_shape'.format(term))
	      break
	#if found_th2 : break
	      
	  
    print "Axis x = {0} y = {1}".format(th2.GetXaxis().GetTitle(), th2.GetYaxis().GetTitle())
    for term in allowedTermNames:
	try:
	    #th2 = f.Get('{0}_shape'.format(term))
	    unfolded = f.Get('{0}_shape'.format(term)).GetName()
	except:
	    print "@@@@ Shape {0}_shape doesn't exist.".format(term)
	else:
	    print 'Working with {0}_shape'.format(term)
	    unfolded = unfoldedHist(f.Get('{0}_shape'.format(term)),'{0}_shape_unfolded'.format(term),norm=-1)
	    fnew.cd()
	    
	    if 'Dbackground_ggZZ' in template_file: 
                print "---------------> Changing name to : ggZZ_shape"
                unfolded.SetNameTitle("ggZZ_shape","ggZZ_shape_unfolded")
	    elif 'Dbackground_ZJets' in template_file: 
                print "---------------> Changing name to : zjets_shape"
                unfolded.SetNameTitle("zjets_shape","zjets_shape_unfolded")
	    else: 
                term_name = term
                if copy_13_templates_to_12 and term in new_names.keys():
                        term_name = new_names[term]
                unfolded.SetNameTitle('{0}_shape'.format(term_name),"{0}_shape_unfolded".format(term))
                
            print "New template name/title: {0}/{1}".format(unfolded.GetName(),unfolded.GetTitle())
            c1 =ROOT.TCanvas("Templates","Templates",800,800)
            c1.cd()
            ROOT.gPad.SetRightMargin(0.085)
            unfolded.Draw("hist")
            misc.make_sure_path_exists(dir_2D_unfolded+"/fig")
            for ext in [".png", ".pdf", ".root",".C"]: c1.SaveAs(dir_2D_unfolded +"/fig/"+template_file+"_"+term+ext)            
            unfolded.Write()
	    
            if not 'Dbackground' in template_file: 
                copy_tree(f, fnew,"factors")
    
    
    f.Close()
    fnew.Close()
   
    
def make_refurbished_hist(template_file, dir_2D):
    
    print "------------- Unfolding  templates"
    #dir_2D = '/afs/cern.ch/work/r/roko/Stat/CMSSW_611_JCP/src/HZZ4L_Combination/CombinationPy/CreateDatacards/Sandbox/templates_compare2Pedja/Templates2D_D0M_Dint13/'
    dir_2D_refurbished = dir_2D+"/refurbished/"
    misc.make_sure_path_exists(dir_2D_refurbished)
    #f=ROOT.TFile('/afs/cern.ch/work/r/roko/Stat/CMSSW_611_JCP/src/HZZ4L_Combination/CombinationPy/CreateDatacards/Sandbox/templates2D_D0M_Dint13/Dsignal_2e2mu.root')
    #f=ROOT.TFile('/afs/cern.ch/work/r/roko/Stat/CMSSW_611_JCP/src/HZZ4L_Combination/CombinationPy/CreateDatacards/Sandbox/Templates2D_D0Ph_Dint12/DSignal_4l.root')
    f=ROOT.TFile(dir_2D +"/"+template_file,"READ")
    print f.GetName()
    #f.cd()
    #ggH_proj = f.Get('ggH_shape').ProjectionX("ggH_shape_D0M",0,-1,"e")
    allowedTermNames = ['ggH',    'gg0Ph','gg0M' ,'ggInt_12P','ggInt_12N','ggInt_13P','ggInt_13N','ggInt_23P','ggInt_23N',
                        'qqZZ','ggZZ','zjets']
    new_names = {
        'gg0M' : 'gg0Ph',
        'ggInt_13P':'ggInt_12P',
        'ggInt_13N':'ggInt_12N',
        }
                        
                        
    #allowedTermNames = ['ggH',   'gg0M', 'ggInt_13P','ggInt_13N']
    fnew=ROOT.TFile(dir_2D_refurbished +"/"+template_file,"RECREATE")
    print fnew.GetName()
    refurbished = None
    #th2 = None
    found_th2 = False
    th2 = None
    for term in allowedTermNames:
        print "Searching for :", term
        try :
              th2_name = f.Get('{0}_shape'.format(term)).GetTitle()
        except ReferenceError:
              print "Doesn't exist : ", term
        else:
              print "----", th2_name 
              found_th2=True
              print "Found term for axis definiton: ", term
              th2 = f.Get('{0}_shape'.format(term))
              break
        #if found_th2 : break
              
          
    print "Axis x = {0} y = {1}".format(th2.GetXaxis().GetTitle(), th2.GetYaxis().GetTitle())
    for term in allowedTermNames:
        try:
            #th2 = f.Get('{0}_shape'.format(term))
            refurbished = f.Get('{0}_shape'.format(term)).GetName()
        except:
            print "@@@@ Shape {0}_shape doesn't exist.".format(term)
        else:
            print 'Working with {0}_shape'.format(term)
            refurbished = refurbishedHist(f.Get('{0}_shape'.format(term)),'{0}_shape_refurbished'.format(term),norm=-1)
            #refurbished = f.Get('{0}_shape'.format(term))
            #refurbished.SetTitle('{0}_shape_refurbished'.format(term))
            fnew.cd()
            
            if 'Dbackground_ggZZ' in template_file: 
                print "---------------> Changing name to : ggZZ_shape"
                refurbished.SetNameTitle("ggZZ_shape","ggZZ_shape_refurbished")
            elif 'Dbackground_ZJets' in template_file: 
                print "---------------> Changing name to : zjets_shape"
                refurbished.SetNameTitle("zjets_shape","zjets_shape_refurbished")
            else: 
                term_name = term
                if copy_13_templates_to_12 and term in new_names.keys():
                        term_name = new_names[term]
                refurbished.SetNameTitle('{0}_shape'.format(term_name),"{0}_shape_refurbished".format(term))
                
            print "New template name/title: {0}/{1}".format(refurbished.GetName(),refurbished.GetTitle())
            c1 =ROOT.TCanvas("Templates","Templates",800,800)
            c1.cd()
            ROOT.gPad.SetRightMargin(0.085)
            refurbished.Draw("colz")
            misc.make_sure_path_exists(dir_2D_refurbished+"/fig")
            for ext in [".png", ".pdf", ".root",".C"]: c1.SaveAs(dir_2D_refurbished +"/fig/"+template_file+"_"+term+ext)            
            refurbished.Write()
            
            if not 'Dbackground' in template_file: 
                copy_tree(f, fnew,"factors")
    
    
    f.Close()
    fnew.Close()
    
def copy_tree(f_src, f_dest,tree_name):
    #Get old file, old tree and set top branch address
    f_src.cd()
    oldtree = f_src.Get(tree_name)
    oldtree.SetBranchStatus("*",1)
    #clone of old tree in new file
    f_dest.cd()
    newtree = oldtree.CloneTree();
    newtree.Print()
    newtree.Write()

    
if __name__ == "__main__":
     
    # parse the arguments and options
    global opt, args
    parseOptions()
    print "Starting projection and unfolding ..."
    
    copy_13_templates_to_12 = False
    #make_projection()
    
    copy_map = {
            ###Discriminant templates for k2k1_ratio
            'DSignal_D0M_Dint13_2e2mu.root':['Dsignal_2e2mu.root','Dsignal_superMELA_2e2mu.root'],
            'DSignal_D0M_Dint13_4e.root':['Dsignal_4e.root','Dsignal_superMELA_4e.root'],
            'DSignal_D0M_Dint13_4mu.root':['Dsignal_4mu.root','Dsignal_superMELA_4mu.root'],

            'DBackground_D0M_Dint13_qqZZ_2e2mu.root':['Dbackground_qqZZ_2e2mu.root','Dbackground_qqZZ_superMELA_2e2mu.root',],
            'DBackground_D0M_Dint13_qqZZ_4e.root':['Dbackground_qqZZ_4e.root',      'Dbackground_qqZZ_superMELA_4e.root'],
            'DBackground_D0M_Dint13_qqZZ_4mu.root':['Dbackground_qqZZ_4mu.root','Dbackground_qqZZ_superMELA_4mu.root',
                                                    'Dbackground_ggZZ_2e2mu.root','Dbackground_ggZZ_4e.root','Dbackground_ggZZ_4mu.root',
                                                    'Dbackground_ZJetsCR_AllChans.root'],
                                                    
            ###Discriminant templates for k2k1_ratio
            'DSignal_D0Ph_Dint12_2e2mu.root':['Dsignal_2e2mu.root','Dsignal_superMELA_2e2mu.root'],
            'DSignal_D0Ph_Dint12_4e.root':['Dsignal_4e.root','Dsignal_superMELA_4e.root'],
            'DSignal_D0Ph_Dint12_4mu.root':['Dsignal_4mu.root','Dsignal_superMELA_4mu.root'],

            'DBackground_D0Ph_Dint12_qqZZ_2e2mu.root':['Dbackground_qqZZ_2e2mu.root','Dbackground_qqZZ_superMELA_2e2mu.root',],
            'DBackground_D0Ph_Dint12_qqZZ_4e.root':['Dbackground_qqZZ_4e.root',      'Dbackground_qqZZ_superMELA_4e.root'],
            'DBackground_D0Ph_Dint12_qqZZ_4mu.root':['Dbackground_qqZZ_4mu.root','Dbackground_qqZZ_superMELA_4mu.root',
                                                    'Dbackground_ggZZ_2e2mu.root','Dbackground_ggZZ_4e.root','Dbackground_ggZZ_4mu.root',
                                                    'Dbackground_ZJetsCR_AllChans.root'],
                                                    
        }
    
    
    files_for_DC = [
      'Dsignal_2e2mu.root',
      'Dsignal_4e.root',
      'Dsignal_4mu.root',
      'Dsignal_superMELA_2e2mu.root',
      'Dsignal_superMELA_4e.root',
      'Dsignal_superMELA_4mu.root',
      
      'Dbackground_qqZZ_2e2mu.root',
      'Dbackground_qqZZ_4e.root',
      'Dbackground_qqZZ_4mu.root',
      'Dbackground_qqZZ_superMELA_2e2mu.root',
      'Dbackground_qqZZ_superMELA_4e.root',
      'Dbackground_qqZZ_superMELA_4mu.root',
      
      'Dbackground_ggZZ_2e2mu.root',
      'Dbackground_ggZZ_4e.root',
      'Dbackground_ggZZ_4mu.root',
      'Dbackground_ZJetsCR_AllChans.root'
      ]
    dir_2D = [
        '/afs/cern.ch/work/r/roko/Stat/CMSSW_611_JCP/src/HZZ4L_Combination/CombinationPy/CreateDatacards/Sandbox/templates_compare2Pedja/LF_RECO_8TeV/Templates2D_D0M_Dint13/',
        #'/afs/cern.ch/work/r/roko/Stat/CMSSW_611_JCP/src/HZZ4L_Combination/CombinationPy/CreateDatacards/Sandbox/templates_compare2Pedja/LF_GEN_8TeV/Templates2D_D0Ph_Dint12/',
        #'/afs/cern.ch/work/r/roko/Stat/CMSSW_611_JCP/src/HZZ4L_Combination/CombinationPy/CreateDatacards/Sandbox/templates_compare2Pedja/LF_RECO_8TeV/Templates2D_D0Ph_Dint12/',
#'/afs/cern.ch/work/r/roko/Stat/CMSSW_611_JCP/src/HZZ4L_Combination/CombinationPy/CreateDatacards/Sandbox/templates_compare2Pedja/LF_RECO_8TeV/Templates2D_D0M_Dint13_for_k2k1/',
#'/afs/cern.ch/work/r/roko/Stat/CMSSW_611_JCP/src/HZZ4L_Combination/CombinationPy/CreateDatacards/Sandbox/templates_compare2Pedja/LF_GEN_8TeV/Templates2D_D0M_Dint13_test_copy_tree/',
        #'/afs/cern.ch/work/r/roko/Stat/CMSSW_611_JCP/src/HZZ4L_Combination/CombinationPy/CreateDatacards/Sandbox/templates_compare2Pedja/LF_GEN_8TeV/Templates2D_D0M_Dint13/'
        ]  
    #copy to proper name
    import shutil
    for origTmplFile in copy_map.keys():
        for newFile in copy_map[origTmplFile]: 
            if os.path.exists(origTmplFile):
                shutil.copy(origTmplFile, newFile)
    
    #make projection and unfolded
    #if opt.run_on_dir!="" and os.path.isdir(opt.run_on_dir):
        #dir_2D=[]
        #dir_2D.append(opt.run_on_dir)
        
    for theDir in dir_2D:
        for file_name in files_for_DC:
            print "------------------------------ INPUT FILE:", file_name
            #make_unfolded_hist(file_name, theDir)
            #make_projected_hist(file_name, theDir)
            make_refurbished_hist(file_name, theDir)
	
    #make_unfolded_hist_bkg()
    #sys.exit()
    
    
    