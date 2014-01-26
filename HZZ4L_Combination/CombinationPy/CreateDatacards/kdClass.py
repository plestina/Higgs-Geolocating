#! /usr/bin/env python
import sys
import os
import re
import math
from scipy.special import erf
from ROOT import *
import ROOT
from array import array
from datacardClass import *


class kdClass(datacardClass):


    def setKD(self):

        self.discVarName = ""

        if self.useMEKD:
            self.discVarName = "mekdLLD"
        else:
            self.discVarName = "melaLD"
        print "discVarName ",self.discVarName

        templateSigName = "{0}/Dsignal_{1}.root".format(self.templateDir ,self.appendName)
        self.sigTempFile = ROOT.TFile(templateSigName)
        
        self.sigTemplate = self.sigTempFile.Get("h_mzzD")
        self.sigTemplate_Up = self.sigTempFile.Get("h_mzzD_up")
        self.sigTemplate_Down = self.sigTempFile.Get("h_mzzD_dn")
            

        
        #Set Bins
        dBins = self.sigTemplate.GetYaxis().GetNbins()
        dLow = self.sigTemplate.GetYaxis().GetXmin()
        dHigh = self.sigTemplate.GetYaxis().GetXmax()
        self.D = ROOT.RooRealVar(self.discVarName,self.discVarName,dLow,dHigh)
        self.D.setBins(dBins)
        print "discVarName ", self.discVarName, " dLow " , dLow, " dHigh ", dHigh, " dBins ", dBins


    ## --------------------------- DATASET --------------------------- ##
    def fetchDatasetKD(self):

        dataFileDir = "CMSdata"
        dataTreeName = "data_obs" 
        dataFileName = "{0}/hzz{1}_{2}.root".format(dataFileDir,self.appendName,self.lumi)
        if (self.DEBUG): print dataFileName," ",dataTreeName 
        data_obs_file = ROOT.TFile(dataFileName)
        
        print data_obs_file.Get(dataTreeName)
        
        if not (data_obs_file.Get(dataTreeName)):
            print "File, \"",dataFileName,"\", or tree, \"",dataTreeName,"\", not found" 
            print "Exiting..."
            sys.exit()
            
        self.data_obs_tree = data_obs_file.Get(dataTreeName)
        self.data_obs = ROOT.RooDataSet()
        self.datasetName = "data_obs_{0}".format(self.appendName)
            
        if (self.is2D == 0):
            if(self.bIncludingError):
                self.data_obs = ROOT.RooDataSet(self.datasetName,self.datasetName,self.data_obs_tree,ROOT.RooArgSet(self.CMS_zz4l_mass, self.RelErr))
            else:
                self.data_obs = ROOT.RooDataSet(self.datasetName,self.datasetName,self.data_obs_tree,ROOT.RooArgSet(self.CMS_zz4l_mass))
            
        if (self.is2D == 1):
            if(self.bIncludingError):
                self.data_obs = ROOT.RooDataSet(self.datasetName,self.datasetName,self.data_obs_tree,ROOT.RooArgSet(self.CMS_zz4l_mass,self.D,self.RelErr))
            else:
                self.data_obs = ROOT.RooDataSet(self.datasetName,self.datasetName,self.data_obs_tree,ROOT.RooArgSet(self.CMS_zz4l_mass,self.D))
                
        if (self.is2D == 2):
            self.data_obs = ROOT.RooDataSet(self.datasetName,self.datasetName,self.data_obs_tree,ROOT.RooArgSet(self.SD))
            



    def makeKDAnalysis(self):

        TemplateName = "sigTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.sigTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(self.CMS_zz4l_mass,self.D),self.sigTemplate)
        TemplateName = "sigTempDataHist_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.sigTempDataHist_Up = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(self.CMS_zz4l_mass,self.D),self.sigTemplate_Up)
        TemplateName = "sigTempDataHist_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.sigTempDataHist_Down = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(self.CMS_zz4l_mass,self.D),self.sigTemplate_Down)

        
        TemplateName = "sigTemplatePdf_ggH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.sigTemplatePdf_ggH = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.CMS_zz4l_mass,self.D),self.sigTempDataHist)
        TemplateName = "sigTemplatePdf_ggH_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.sigTemplatePdf_ggH_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.CMS_zz4l_mass,self.D),self.sigTempDataHist_Up)
        TemplateName = "sigTemplatePdf_ggH_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.sigTemplatePdf_ggH_Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.CMS_zz4l_mass,self.D),self.sigTempDataHist_Down)
        
        TemplateName = "sigTemplatePdf_VBF_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.sigTemplatePdf_VBF = ROOT.RooHistPdf(TemplateName,TemplateName,RooArgSet(self.CMS_zz4l_mass,self.D),self.sigTempDataHist)
        TemplateName = "sigTemplatePdf_VBF_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.sigTemplatePdf_VBF_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.CMS_zz4l_mass,self.D),self.sigTempDataHist_Up)
        TemplateName = "sigTemplatePdf_VBF_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.sigTemplatePdf_VBF_Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.CMS_zz4l_mass,self.D),self.sigTempDataHist_Down)
        
        TemplateName = "sigTemplatePdf_WH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.sigTemplatePdf_WH = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.CMS_zz4l_mass,self.D),self.sigTempDataHist)
        TemplateName = "sigTemplatePdf_WH_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.sigTemplatePdf_WH_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.CMS_zz4l_mass,self.D),self.sigTempDataHist_Up)
        TemplateName = "sigTemplatePdf_WH_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.sigTemplatePdf_WH_Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.CMS_zz4l_mass,self.D),self.sigTempDataHist_Down)
        
        TemplateName = "sigTemplatePdf_ZH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.sigTemplatePdf_ZH = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.CMS_zz4l_mass,self.D),self.sigTempDataHist)
        TemplateName = "sigTemplatePdf_ZH_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.sigTemplatePdf_ZH_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.CMS_zz4l_mass,self.D),self.sigTempDataHist_Up)
        TemplateName = "sigTemplatePdf_ZH_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.sigTemplatePdf_ZH_Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.CMS_zz4l_mass,self.D),self.sigTempDataHist_Down)
        
        TemplateName = "sigTemplatePdf_ZH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.sigTemplatePdf_ttH = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.CMS_zz4l_mass,self.D),self.sigTempDataHist)
        TemplateName = "sigTemplatePdf_ZH_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.sigTemplatePdf_ttH_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.CMS_zz4l_mass,self.D),self.sigTempDataHist_Up)
        TemplateName = "sigTemplatePdf_ZH_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.sigTemplatePdf_ttH_Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.CMS_zz4l_mass,self.D),self.sigTempDataHist_Down)
        
        self.funcList_ggH = ROOT.RooArgList()  
        self.funcList_VBF = ROOT.RooArgList()
        self.funcList_WH  = ROOT.RooArgList()
        self.funcList_ZH  = ROOT.RooArgList()
        self.funcList_ttH = ROOT.RooArgList()

        if(self.sigMorph):
            
            self.funcList_ggH.add(self.sigTemplatePdf_ggH)
            self.funcList_ggH.add(self.sigTemplatePdf_ggH_Up)
            self.funcList_ggH.add(self.sigTemplatePdf_ggH_Down)  
            
            self.funcList_VBF.add(self.sigTemplatePdf_VBF)
            self.funcList_VBF.add(self.sigTemplatePdf_VBF_Up)
            self.funcList_VBF.add(self.sigTemplatePdf_VBF_Down)  
            
            self.funcList_WH.add(self.sigTemplatePdf_WH)
            self.funcList_WH.add(self.sigTemplatePdf_WH_Up)
            self.funcList_WH.add(self.sigTemplatePdf_WH_Down)  
            
            self.funcList_ZH.add(self.sigTemplatePdf_ZH)
            self.funcList_ZH.add(self.sigTemplatePdf_ZH_Up)
            self.funcList_ZH.add(self.sigTemplatePdf_ZH_Down)  
            
            self.funcList_ttH.add(self.sigTemplatePdf_ttH)
            self.funcList_ttH.add(self.sigTemplatePdf_ttH_Up)
            self.funcList_ttH.add(self.sigTemplatePdf_ttH_Down)  

        else:
            
            self.funcList_ggH.add(self.sigTemplatePdf_ggH)
            self.funcList_VBF.add(self.sigTemplatePdf_VBF)
            self.funcList_WH.add(self.sigTemplatePdf_WH)
            self.funcList_ZH.add(self.sigTemplatePdf_ZH)
            self.funcList_ttH.add(self.sigTemplatePdf_ttH)
   

        morphSigVarName = "CMS_zz4l_sigMELA_{0:.0f}".format(self.channel)
        self.alphaMorphSig = ROOT.RooRealVar(morphSigVarName,morphSigVarName,0,-20,20)
        if(self.sigMorph):
            self.alphaMorphSig.setConstant(False)
        else:
            self.alphaMorphSig.setConstant(True)
            
        self.morphVarListSig = ROOT.RooArgList()
        
        if(self.sigMorph):
            self.morphVarListSig.add(self.alphaMorphSig)  ## just one morphing for all signal processes (fully correlated systematics)


        
        TemplateName = "sigTemplateMorphPdf_ggH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.sigTemplateMorphPdf_ggH = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,self.CMS_zz4l_mass,self.D,True,self.funcList_ggH,self.morphVarListSig,1.0,1)
        
        TemplateName = "sigTemplateMorphPdf_VBF_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.sigTemplateMorphPdf_VBF = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,self.CMS_zz4l_mass,self.D,True,self.funcList_VBF,self.morphVarListSig,1.0,1)
        
        TemplateName = "sigTemplateMorphPdf_WH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.sigTemplateMorphPdf_WH = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,self.CMS_zz4l_mass,self.D,True,self.funcList_WH,self.morphVarListSig,1.0,1)
        
        TemplateName = "sigTemplateMorphPdf_ZH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.sigTemplateMorphPdf_ZH = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,self.CMS_zz4l_mass,self.D,True,self.funcList_ZH,self.morphVarListSig,1.0,1)
        
        TemplateName = "sigTemplateMorphPdf_ttH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.sigTemplateMorphPdf_ttH = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,self.CMS_zz4l_mass,self.D,True,self.funcList_ttH,self.morphVarListSig,1.0,1)


        self.sig2d_ggH = ROOT.RooProdPdf("sig2d_ggH","sig2d_ggH",ROOT.RooArgSet(self.getVariable(self.sig_ggH_HM,self.sig_ggH,self.isHighMass)),ROOT.RooFit.Conditional(ROOT.RooArgSet(self.sigTemplateMorphPdf_ggH),ROOT.RooArgSet(self.D)))
        self.sig2d_VBF = ROOT.RooProdPdf("sig2d_VBF","sig2d_VBF",ROOT.RooArgSet(self.getVariable(self.sig_VBF_HM,self.sig_VBF,self.isHighMass)),ROOT.RooFit.Conditional(ROOT.RooArgSet(self.sigTemplateMorphPdf_VBF),ROOT.RooArgSet(self.D)))
        self.sig2d_WH = ROOT.RooProdPdf("sig2d_WH","sig2d_WH",ROOT.RooArgSet(self.getVariable(self.sig_WH_HM,self.sig_WH,self.isHighMass)),ROOT.RooFit.Conditional(ROOT.RooArgSet(self.sigTemplateMorphPdf_WH),ROOT.RooArgSet(self.D)))
        self.sig2d_ZH = ROOT.RooProdPdf("sig2d_ZH","sig2d_ZH",ROOT.RooArgSet(self.getVariable(self.sig_ZH_HM,self.sig_ZH,self.isHighMass)),ROOT.RooFit.Conditional(ROOT.RooArgSet(self.sigTemplateMorphPdf_ZH),ROOT.RooArgSet(self.D)))
        self.sig2d_ttH = ROOT.RooProdPdf("sig2d_ttH","sig2d_ttH",ROOT.RooArgSet(self.getVariable(self.sig_ttH_HM,self.sig_ttH,self.isHighMass)),ROOT.RooFit.Conditional(ROOT.RooArgSet(self.sigTemplateMorphPdf_ttH),ROOT.RooArgSet(self.D)))
                
        self.sigCB2d_ggH = ROOT.RooProdPdf("sigCB2d_ggH","sigCB2d_ggH",ROOT.RooArgSet(self.getVariable(self.sig_ggHErr,self.signalCB_ggH,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(self.sigTemplateMorphPdf_ggH),ROOT.RooArgSet(self.D)))
        self.sigCB2d_VBF = ROOT.RooProdPdf("sigCB2d_VBF","sigCB2d_VBF",ROOT.RooArgSet(self.getVariable(self.sig_VBFErr,self.signalCB_VBF,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(self.sigTemplateMorphPdf_VBF),ROOT.RooArgSet(self.D)))
        self.sigCB2d_WH = ROOT.RooProdPdf("sigCB2d_WH","sigCB2d_WH",ROOT.RooArgSet(self.getVariable(self.sig_WHErr,self.signalCB_WH,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(self.sigTemplateMorphPdf_WH),ROOT.RooArgSet(self.D)))
        self.sigCB2d_ZH = ROOT.RooProdPdf("sigCB2d_ZH","sigCB2d_ZH",ROOT.RooArgSet(self.getVariable(self.sig_ZHErr,self.signalCB_ZH,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(self.sigTemplateMorphPdf_ZH),ROOT.RooArgSet(self.D)))
        self.sigCB2d_ttH = ROOT.RooProdPdf("sigCB2d_ttH","sigCB2d_ttH",ROOT.RooArgSet(self.getVariable(self.sig_ttHErr,self.signalCB_ttH,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(self.sigTemplateMorphPdf_ttH),ROOT.RooArgSet(self.D)))



      ## ----------------- 2D BACKGROUND SHAPES --------------- ##
        templateBkgName = "{0}/Dbackground_qqZZ_{1}.root".format(self.templateDir ,self.appendName)
        bkgTempFile = ROOT.TFile(templateBkgName)
        self.bkgTemplate = bkgTempFile.Get("h_mzzD")
        self.bkgTemplate_Up = bkgTempFile.Get("h_mzzD_up")
        self.bkgTemplate_Down = bkgTempFile.Get("h_mzzD_dn")
        
        TemplateName = "bkgTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        self.bkgTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(self.CMS_zz4l_mass,self.D),self.bkgTemplate)
        TemplateName = "bkgTempDataHist_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.bkgTempDataHist_Up = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(self.CMS_zz4l_mass,self.D),self.bkgTemplate_Up)
        TemplateName = "bkgTempDataHist_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.bkgTempDataHist_Down = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(self.CMS_zz4l_mass,self.D),self.bkgTemplate_Down)
        
        templateggBkgName = "{0}/Dbackground_ggZZ_{1}.root".format(self.templateDir ,self.appendName)
        ggbkgTempFile = ROOT.TFile(templateggBkgName)
        self.ggbkgTemplate = ggbkgTempFile.Get("h_mzzD")
        TemplateName = "ggbkgTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        self.ggbkgTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(self.CMS_zz4l_mass,self.D),self.ggbkgTemplate)
        
        TemplateName = "bkgTemplatePdf_qqzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        self.bkgTemplatePdf_qqzz = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.CMS_zz4l_mass,self.D),self.bkgTempDataHist)
        TemplateName = "bkgTemplatePdf_ggzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        self.bkgTemplatePdf_ggzz = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.CMS_zz4l_mass,self.D),self.ggbkgTempDataHist)


        ### NEW MEKD
        if self.useMEKD :
            templateBkgName = "{0}/Dbackground_ZX_{1}.root".format(self.templateDir ,self.appendName)
        else:
            templateBkgName = "{0}/Dbackground_qqZZ_{1}.root".format(self.templateDir ,self.appendName)

        print "Using {0} as ZX KD template".format(templateBkgName)

        bkgTempFileZX = ROOT.TFile(templateBkgName)
        self.bkgTemplateZX = bkgTempFileZX.Get("h_mzzD")
        self.bkgTemplateZX_Up = bkgTempFileZX.Get("h_mzzD_up")
        self.bkgTemplateZX_Down = bkgTempFileZX.Get("h_mzzD_dn")

        TemplateName = "bkgTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.bkgTempDataHistZX = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(self.CMS_zz4l_mass,self.D),self.bkgTemplateZX)
        TemplateName = "bkgTempDataHist_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.bkgTempDataHistZX_Up = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(self.CMS_zz4l_mass,self.D),self.bkgTemplateZX_Up)
        TemplateName = "bkgTempDataHist_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.bkgTempDataHistZX_Down = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(self.CMS_zz4l_mass,self.D),self.bkgTemplateZX_Down)
        
        TemplateName = "bkgTemplatePdf_zjets_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        self.bkgTemplatePdf_zjets = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.CMS_zz4l_mass,self.D),self.bkgTempDataHistZX)
        TemplateName = "bkgTemplatePdf_zjets_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        self.bkgTemplatePdf_zjets_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.CMS_zz4l_mass,self.D),self.bkgTempDataHistZX_Up)
        TemplateName = "bkgTemplatePdf_zjets_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        self.bkgTemplatePdf_zjets_Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.CMS_zz4l_mass,self.D),self.bkgTempDataHistZX_Down)
        
        self.funcList_zjets = ROOT.RooArgList()  
        morphBkgVarName = "CMS_zz4l_bkgMELA"    
        self.alphaMorphBkg = ROOT.RooRealVar(morphBkgVarName,morphBkgVarName,0,-20,20)
        self.morphVarListBkg = ROOT.RooArgList()
        
        if(self.bkgMorph):
            self.funcList_zjets.add(self.bkgTemplatePdf_zjets)
            self.funcList_zjets.add(self.bkgTemplatePdf_zjets_Up)
            self.funcList_zjets.add(self.bkgTemplatePdf_zjets_Down)  
            self.alphaMorphBkg.setConstant(False)
            self.morphVarListBkg.add(self.alphaMorphBkg)  
        else:
            self.funcList_zjets.add(self.bkgTemplatePdf_zjets)
            self.alphaMorphBkg.setConstant(True)


        TemplateName = "bkgTemplateMorphPdf_qqzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        self.bkgTemplateMorphPdf_qqzz = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,self.CMS_zz4l_mass,self.D,True,ROOT.RooArgList(self.bkgTemplatePdf_qqzz),ROOT.RooArgList(),1.0,1)
        TemplateName = "bkgTemplateMorphPdf_ggzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        self.bkgTemplateMorphPdf_ggzz = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,self.CMS_zz4l_mass,self.D,True,ROOT.RooArgList(self.bkgTemplatePdf_ggzz),ROOT.RooArgList(),1.0,1)
        TemplateName = "bkgTemplateMorphPdf_zjets_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        self.bkgTemplateMorphPdf_zjets = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,self.CMS_zz4l_mass,self.D,True,self.funcList_zjets,self.morphVarListBkg,1.0,1)


        self.bkg2d_qqzz = ROOT.RooProdPdf("bkg2d_qqzz","bkg2d_qqzz",ROOT.RooArgSet(self.getVariable(self.bkg_qqzzErr,self.bkg_qqzz,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(self.bkgTemplateMorphPdf_qqzz),ROOT.RooArgSet(self.D)))
        self.bkg2d_ggzz = ROOT.RooProdPdf("bkg2d_ggzz","bkg2d_ggzz",ROOT.RooArgSet(self.getVariable(self.bkg_ggzzErr,self.bkg_ggzz,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(self.bkgTemplateMorphPdf_ggzz),ROOT.RooArgSet(self.D)))
        self.bkg2d_zjets = ROOT.RooProdPdf("bkg2d_zjets","bkg2d_zjets",ROOT.RooArgSet(self.getVariable(self.bkg_zjetsErr,self.bkg_zjets,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(self.bkgTemplateMorphPdf_zjets),ROOT.RooArgSet(self.D)))

            
    ## --------------------------- WORKSPACE -------------------------- ##
    def writeWorkspaceKD(self):
            
        endsInP5 = False
        tmpMH = self.mH
        if (math.fabs(math.floor(tmpMH)-self.mH) > 0.001): endsInP5 = True
        if (self.DEBUG): print "ENDS IN P5  ",endsInP5
        
        self.name_Shape = ""
        self.name_ShapeWS = ""
        self.name_ShapeWS2 = ""
        self.name_ShapeWSXSBR = ""
        
        if (endsInP5):
            self.name_Shape = "{0}/HCG/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts)
        else:
            self.name_Shape = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts)
            
        if (endsInP5):
            self.name_ShapeWS = "{0}/HCG/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts)
        else:
            self.name_ShapeWS = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts)
            
        if (endsInP5):
            self.name_ShapeWSXSBR = "{0}/HCG_XSxBR/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts)
        else:
            self.name_ShapeWSXSBR = "{0}/HCG_XSxBR/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts)

        self.name_ShapeWS2 = "hzz4l_{0}S_{1:.0f}TeV.input.root".format(self.appendName,self.sqrts)
        
        if(self.DEBUG): print self.name_Shape,"  ",self.name_ShapeWS2
        
        self.w = ROOT.RooWorkspace("w","w")
        
        self.w.importClassCode(RooqqZZPdf_v2.Class(),True)
        self.w.importClassCode(RooggZZPdf_v2.Class(),True)
        self.w.importClassCode(RooRelBWUFParam.Class(),True)
        self.w.importClassCode(RooDoubleCB.Class(),True)
        self.w.importClassCode(RooFormulaVar.Class(),True)
        if self.isHighMass :
            self.w.importClassCode(RooRelBWHighMass.Class(),True)
                
        getattr(self.w,'import')(self.data_obs,ROOT.RooFit.Rename("data_obs")) 
                
        if(self.bUseCBnoConvolution) :
            if self.is2D == 1:

                self.sigCB2d_ggH.SetNameTitle("ggH","ggH")
                self.sigCB2d_VBF.SetNameTitle("qqH","qqH")
                self.sigCB2d_WH.SetNameTitle("WH","WH")
                self.sigCB2d_ZH.SetNameTitle("ZH","ZH")
                self.sigCB2d_ttH.SetNameTitle("ttH","ttH")
                
                getattr(self.w,'import')(self.sigCB2d_ggH, ROOT.RooFit.RecycleConflictNodes())
                getattr(self.w,'import')(self.sigCB2d_VBF, ROOT.RooFit.RecycleConflictNodes())
                getattr(self.w,'import')(self.sigCB2d_WH, ROOT.RooFit.RecycleConflictNodes())
                getattr(self.w,'import')(self.sigCB2d_ZH, ROOT.RooFit.RecycleConflictNodes())
                getattr(self.w,'import')(self.sigCB2d_ttH, ROOT.RooFit.RecycleConflictNodes())

        else:

            if self.is2D == 1:
                self.sig2d_ggH.SetNameTitle("ggH","ggH")
                self.sig2d_VBF.SetNameTitle("qqH","qqH")
                self.sig2d_WH.SetNameTitle("WH","WH")
                self.sig2d_ZH.SetNameTitle("ZH","ZH")
                self.sig2d_ttH.SetNameTitle("ttH","ttH")     
                
                getattr(self.w,'import')(self.sig2d_ggH,ROOT.RooFit.RecycleConflictNodes())
                getattr(self.w,'import')(self.sig2d_VBF,ROOT.RooFit.RecycleConflictNodes())
                getattr(self.w,'import')(self.sig2d_WH,ROOT.RooFit.RecycleConflictNodes())
                getattr(self.w,'import')(self.sig2d_ZH,ROOT.RooFit.RecycleConflictNodes())
                getattr(self.w,'import')(self.sig2d_ttH,ROOT.RooFit.RecycleConflictNodes())            

 
        getattr(self.w,'import')(self.bkg2d_qqzz,ROOT.RooFit.RecycleConflictNodes())
        getattr(self.w,'import')(self.bkg2d_ggzz,ROOT.RooFit.RecycleConflictNodes())
        getattr(self.w,'import')(self.bkg2d_zjets,ROOT.RooFit.RecycleConflictNodes())
                
        getattr(self.w,'import')(self.rfvSigRate_ggH, ROOT.RooFit.RecycleConflictNodes())
        getattr(self.w,'import')(self.rfvSigRate_VBF, ROOT.RooFit.RecycleConflictNodes())
        getattr(self.w,'import')(self.rfvSigRate_WH, ROOT.RooFit.RecycleConflictNodes())
        getattr(self.w,'import')(self.rfvSigRate_ZH, ROOT.RooFit.RecycleConflictNodes())
        getattr(self.w,'import')(self.rfvSigRate_ttH, ROOT.RooFit.RecycleConflictNodes())

            
        self.w.writeToFile(self.name_ShapeWS)
        self.w.writeToFile(self.name_ShapeWSXSBR)
        

    def prepareDatacardKD(self):

        name_Shape = ""
        name_ShapeWS = ""

        ## Write Datacards
        if (self.endsInP5(self.mH)):
            name_Shape = "{0}/HCG/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts)
        else:
            name_Shape = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts)
        fo = open( name_Shape, "wb")
        self.WriteDatacardKD(fo, self.name_ShapeWS2, self.rates, self.data_obs.numEntries())
        self.systematics.WriteSystematics(fo, self.inputs)
        self.systematics.WriteShapeSystematics(fo,self.inputs)
        fo.close()

        if (self.endsInP5(self.mH)):
            name_Shape = "{0}/HCG_XSxBR/{2:.1f}/hzz4l_{1}S_{3:.0f}TeV.txt".format(self.outputDir,self.appendName,self.mH,self.sqrts)	
        else:
            name_Shape = "{0}/HCG_XSxBR/{2:.0f}/hzz4l_{1}S_{3:.0f}TeV.txt".format(self.outputDir,self.appendName,self.mH,self.sqrts)
        fo = open( name_Shape, "wb" )
        self.WriteDatacardKD(fo,self.name_ShapeWS2, self.rates, self.data_obs.numEntries())
        self.systematics_forXSxBR.WriteSystematics(fo, self.inputs)
        self.systematics_forXSxBR.WriteShapeSystematics(fo,self.inputs)
        fo.close()




        
    def WriteDatacardKD(self,file,nameWS,theRates,obsEvents):

        numberSig = self.numberOfSigChan(self.inputs)
        numberBg  = self.numberOfBgChan(self.inputs)

        file.write("imax 1\n")
        file.write("jmax {0}\n".format(numberSig+numberBg-1))
        file.write("kmax *\n")
        
        file.write("------------\n")
        file.write("shapes * * {0} w:$PROCESS \n".format(nameWS))
        file.write("------------\n")
        
        file.write("bin a{0} \n".format(self.channel))
        file.write("observation {0} \n".format(obsEvents))
        
        file.write("------------\n")
        file.write("## mass window [{0},{1}] \n".format(self.inputs['low_M'],self.inputs['high_M']))
        file.write("bin ")        

        channelList=['ggH','qqH','WH','ZH','ttH','qqZZ','ggZZ','zjets','ttbar','zbb']
        channelName2D=['ggH','qqH','WH','ZH','ttH','bkg2d_qqzz','bkg2d_ggzz','bkg2d_zjets','bkg2d_ttbar','bkg2d_zbb']

        if self.inputs["all"]:
            print "ALL CHANNELS --- KDCLASS --- WriteDatacardKD\n"
            channelList=['ggH','qqZZ','ggZZ','zjets','ttbar','zbb']
            channelName2D=['ggH','bkg2d_qqzz','bkg2d_ggzz','bkg2d_zjets','bkg2d_ttbar','bkg2d_zbb']
          
        for chan in channelList:
            if self.inputs[chan]:
                file.write("a{0} ".format(self.channel))
            else:
                if chan.startswith("ggH") and self.inputs["all"] :
                    file.write("a{0} ".format(self.channel))


        file.write("\n")
                                        
        file.write("process ")

        i=0
        for chan in channelList:
            print chan
            if self.inputs[chan]:
                print "passed first if"
                file.write("{0} ".format(channelName2D[i]))
                print 'writing in card index=',i,'  chan=',chan
            else:
                if self.inputs["all"]:
                    print 'writing in card index=',i,'  chan=',chan
                    file.write("{0} ".format(channelName2D[i]))
                    print channelName2D[i]
            i+=1


        
        file.write("\n")
            
        processLine = "process "

        for x in range(-numberSig+1,1):
            processLine += "{0} ".format(x)

        for y in range(1,numberBg+1):
            processLine += "{0} ".format(y)

        file.write(processLine)
        file.write("\n")
            
        file.write("rate ")
        for chan in channelList:
            if self.inputs[chan] or self.inputs["all"]:
                file.write("{0:.4f} ".format(theRates[chan]))
        file.write("\n")
        file.write("------------\n")


