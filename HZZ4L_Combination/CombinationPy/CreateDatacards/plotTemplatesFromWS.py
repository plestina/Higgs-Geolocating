#! /usr/bin/env python
import sys
import os
import re
import math
from ROOT import *
import ROOT
from array import array

def parseOptions():
    pass
    #usage = ('usage: %prog [options] \n'
             #+ '%prog -h for help')
    #parser = optparse.OptionParser(usage)
    
    ##parser.add_option('-d', '--is2D',   dest='is2D',       type='int',    default=1,     help='is2D (default:1)')
    ##parser.add_option('-a', '--append', dest='appendName', type='string', default="",    help='append name for cards dir')
    #parser.add_option('-t', action='store_true', dest='t2w', default=False ,help='do text2workspace')
    #parser.add_option('-g', action='store_true', dest='gen', default=False ,help='do generate asimov')
    #parser.add_option('-a', action='store_true', dest='addasimov', default=False ,help='do add asimov dataset')
    #parser.add_option('-f', action='store_true', dest='fit', default=False ,help='do MultiDimFit')
    #parser.add_option('-p', action='store_true', dest='plot', default=False ,help='do plot')
    
    
    ## store options and arguments as global variables
    #global opt, args
    #(opt, args) = parser.parse_args()

def make_plots(ws_name):
    
    ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
    print "Ploting templates"
    f=ROOT.TFile(ws_name)
    f.cd()
    w = f.Get('w')
    var_vals = {
        'MH':126,
        'k2k1_ratio':0,
        'k3k1_ratio':0,
        }
    for v in var_vals.keys():
        w.var(v).setVal(var_vals[v])
    
    #rrv_Djcp = w.var("CMS_zz4l_Djcp")
    
    plot = w.var("CMS_zz4l_Djcp").frame(RooFit.Name("Djcp_frame"),RooFit.Title("Templates"))
    
    
    
    tot_pdf = w.pdf('pdf_binch3_nuis')
    tot_pdf.plotOn(plot)
    #tot_pdf.plotOn(plot, RooFit.Components("shapeSig_ggH_ch3"), RooFit.LineColor(kRed))
    #tot_pdf.plotOn(plot, RooFit.Components("shapeSig_gg0Ph_ch3"), RooFit.LineColor(209))
    #tot_pdf.plotOn(plot, RooFit.Components("shapeSig_ggInt_12N_ch3"), RooFit.LineColor(222))
    
    
    plot.Draw()
    th2 = tot_pdf.createHistogram("CMS_zz4l_Djcp,CMS_zz4l_Djcp_int")
    c1 =ROOT.TCanvas("cc1","Templates",1000,800)
    c1.cd()
    th2.Draw('colz')

    
    
    
    
    
    
    #plot = rrv_Djcp.frame()
    #w.pdf("ggH" ).plotOn(plot, RooFit.LineColor(kRed),RooFit.LineStyle(kSolid) )
    #w.pdf("gg0Ph").plotOn(plot, RooFit.LineColor(kBlue),RooFit.LineStyle(kDashed) )
    #w.pdf("ggInt_12N").plotOn(plot, RooFit.LineColor(kGreen),RooFit.LineStyle(kDotted))
    #plot.Draw()
    
    

 #// Change line color to red
  #gauss.plotOn(frame1,LineColor(kRed)) 

  #// Change line style to dashed
  #gauss.plotOn(frame2,LineStyle(kDashed)) 

  #// Filled shapes in green color
  #gauss.plotOn(frame3,DrawOption("F"),FillColor(kOrange),MoveToBack()) 

  #//
  #gauss.plotOn(frame4,Range(-8,3),LineColor(kMagenta)) 

    
    
    
    
    
    plots = ROOT.TFile("plots.root","RECREATE")
    canv = ROOT.TCanvas("template","template",600,600)
    canv.cd()
    plot.Draw()
    canv.SaveAs("plot.png")
    plots.cd()
    plot.Write()
    plots.Close()
	   
    #print 'FINALSTATE----',finalstate+"_"+period    
    #print 'GGH:',w.function("ggH_norm").getVal()
    #print 'QQH',w.function("qqH_norm").getVal()
    #print 'WH',w.function("WH_norm").getVal()
    #print 'ZH',w.function("ZH_norm").getVal()
    #print 'ttH',w.function("ttH_norm").getVal()

    f.Close()



if __name__ == "__main__":
     
    # parse the arguments and options
    global opt, args
    parseOptions()
    print "Starting plots..."
    #ws_name = '/afs/cern.ch/work/r/roko/Stat/CMSSW_611_JCP/src/HZZ4L_Combination/CombinationPy/CreateDatacards/cards_2D.k2k1.0.PedjaIn.IntOn.D0Ph.Dint12.14TeV.M4l131to141.RECO/HCG/126/workspaceWithAsimov_k2k1_ratio_0_lumi_25.root'
    ws_name = '/afs/cern.ch/work/r/roko/Stat/CMSSW_611_JCP/src/HZZ4L_Combination/CombinationPy/CreateDatacards/cards_test_2Dtemplates/HCG/126/workspaceWithAsimov_k3k1_ratio_0_lumi_25.root'
    make_plots(ws_name)
    
    #sys.exit()