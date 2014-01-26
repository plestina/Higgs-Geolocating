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



class superkdClass(datacardClass):


    def setSuperKD(self):
      
	self.sigMorph = False      
	self.killBackground = True
      
        self.altHypothesis = self.inputs['altHypothesis']

        self.discVarName = "CMS_zz4l_pseudoKD"
        if self.altHypothesis == 'gg0-':
            self.discVarName = "CMS_zz4l_pseudoKD"
        elif self.altHypothesis == 'gg2m+':
            self.discVarName = "CMS_zz4l_graviKD"
        elif self.altHypothesis == 'gg0h+':
            self.discVarName = "CMS_zz4l_p0hplusKD"
        elif self.altHypothesis == 'qq1+':
            self.discVarName = "CMS_zz4l_p1plusKD"
        elif self.altHypothesis == 'qq1-':
            self.discVarName = "CMS_zz4l_p1minusKD"
        elif self.altHypothesis == 'qq2m+':
            self.discVarName = "CMS_zz4l_qqgraviKD"
        elif self.altHypothesis == '2h+':
            self.discVarName = "CMS_zz4l_p2hplusKD"
        elif self.altHypothesis == '2h-':
            self.discVarName = "CMS_zz4l_p2hminusKD"
        elif self.altHypothesis == '2b+':
            self.discVarName = "CMS_zz4l_p2bplusKD"
        elif self.altHypothesis == 'PI2+':
            self.discVarName = "CMS_zz4l_gravi_ProdIndepKD"
        elif self.altHypothesis == 'PI1-':
            self.discVarName = "CMS_zz4l_p1minus_ProdIndepKD"
        elif self.altHypothesis == 'PI1+':
            self.discVarName = "CMS_zz4l_p1plus_ProdIndepKD"
        else :
            self.discVarName = "CMS_zz4l_pseudoKD"
        print '>>>>>> SuperKD 2D PDFS using discriminat named :',self.discVarName
        
        templateSigName = "{0}/Dsignal_{1}.root".format(self.templateDir ,self.appendName)
        self.sigTempFile = ROOT.TFile(templateSigName)
    
        self.sigTemplate = self.sigTempFile.Get("h_superDpsD")
        self.sigTemplate_syst1Up = self.sigTempFile.Get("h_superDpsD_LeptScaleUp")
        self.sigTemplate_syst1Down = self.sigTempFile.Get("h_superDpsD_LeptScaleDown")
        self.sigTemplate_syst2Up = self.sigTempFile.Get("h_superDpsD_LeptSmearUp")
        self.sigTemplate_syst2Down = self.sigTempFile.Get("h_superDpsD_LeptSmearDown")
        
        templateSigName_ALT = "{0}/Dsignal{2}_{1}.root".format(self.templateDir,self.appendName, self.appendHypType)
        print '>>>>>> Taking 2D template for ALT signal from ',templateSigName_ALT
        self.sigTempFile_ALT = ROOT.TFile(templateSigName_ALT)
        self.sigTemplate_ALT = self.sigTempFile_ALT.Get("h_superDpsD")
        self.sigTemplate_ALT_syst1Up = self.sigTempFile_ALT.Get("h_superDpsD_LeptScaleUp")
        self.sigTemplate_ALT_syst1Down = self.sigTempFile_ALT.Get("h_superDpsD_LeptScaleDown")
        self.sigTemplate_ALT_syst2Up = self.sigTempFile_ALT.Get("h_superDpsD_LeptSmearUp")
        self.sigTemplate_ALT_syst2Down = self.sigTempFile_ALT.Get("h_superDpsD_LeptSmearDown")
            
        print '>>>>>> Signal Templates Files: ',templateSigName, ', ',templateSigName_ALT
        
        #### Set Proper Binning
        self.xbins_sig  = ROOT.RooBinning(self.sigTemplate.GetXaxis().GetNbins(),self.sigTemplate.GetXaxis().GetXbins().GetArray())
        self.ybins_sig  = ROOT.RooBinning(self.sigTemplate.GetYaxis().GetNbins(),self.sigTemplate.GetYaxis().GetXbins().GetArray()) 

        ##Set Bins
        dBins = self.sigTemplate.GetYaxis().GetNbins()
        dLow = self.sigTemplate.GetXaxis().GetXmin()
        dHigh = self.sigTemplate.GetXaxis().GetXmax()
        self.D = ROOT.RooRealVar(self.discVarName,self.discVarName,dLow,dHigh)
        self.D.setBins(dBins)
        self.D.setBinning(self.ybins_sig)
        print '>>>>>> discVarName: ', self.discVarName
        print '>>>>>> bins [low,high]: ',dBins,'['+str(dLow)+','+str(dHigh)+']'
        self.ybins_sig.Print()

        #################################
        self.superDiscVarName = "CMS_zz4l_smd"
        dBins = self.sigTemplate_ALT.GetXaxis().GetNbins()
        dLow = self.sigTemplate_ALT.GetXaxis().GetXmin()
        dHigh = self.sigTemplate_ALT.GetXaxis().GetXmax()
        self.SD = ROOT.RooRealVar(self.superDiscVarName,self.superDiscVarName,dLow,dHigh)
        self.SD.setBins(dBins)
        self.SD.setBinning(self.xbins_sig)
        
        print '>>>>>> superDiscVarName: ', self.superDiscVarName
        print '>>>>>> bins [low,high]: ',dBins,'['+str(dLow)+','+str(dHigh)+']'
        self.xbins_sig.Print()

    ## --------------------------- DATASET --------------------------- ##
    def fetchDatasetSuperKD(self):

        self.dataFileDir = "CMSdata"
        self.dataTreeName = "data_obs" 
        self.dataFileName = "{0}/hzz{1}_{2}.root".format(self.dataFileDir,self.appendName,self.lumi)
        if (self.DEBUG): print self.dataFileName," ",self.dataTreeName 
        self.data_obs_file = ROOT.TFile(self.dataFileName)
        
        print self.data_obs_file.Get(self.dataTreeName)
        
        if not (self.data_obs_file.Get(self.dataTreeName)):
            print "File, \"",self.dataFileName,"\", or tree, \"",self.dataTreeName,"\", not found" 
            print "Exiting..."
            sys.exit()
            
        self.data_obs_tree = self.data_obs_file.Get(self.dataTreeName)
        self.data_obs = ROOT.RooDataSet()
        self.datasetName = "data_obs_{0}".format(self.appendName)
            
        self.data_obs = ROOT.RooDataSet(self.datasetName,self.datasetName,self.data_obs_tree,ROOT.RooArgSet(self.CMS_zz4l_mass,self.SD,self.D),'CMS_zz4l_mass>106.0&&CMS_zz4l_mass<141.0').reduce(ROOT.RooArgSet(self.SD,self.D))
                


    def makeSuperKDAnalysis(self):

        TemplateName = "sigTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.sigTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(self.SD,self.D),ROOT.RooFit.Import(self.sigTemplate,kFALSE))
        TemplateName = "sigTempDataHist_syst1Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.sigTempDataHist_syst1Up = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(self.SD,self.D),ROOT.RooFit.Import(self.sigTemplate_syst1Up,kFALSE))
        TemplateName = "sigTempDataHist_syst1Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.sigTempDataHist_syst1Down =ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(self.SD,self.D),ROOT.RooFit.Import(self.sigTemplate_syst1Down,kFALSE))
        TemplateName = "sigTempDataHist_syst2Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.sigTempDataHist_syst2Up = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(self.SD,self.D),ROOT.RooFit.Import(self.sigTemplate_syst1Up,kFALSE))
        TemplateName = "sigTempDataHist_syst2Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.sigTempDataHist_syst2Down = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(self.SD,self.D),ROOT.RooFit.Import(self.sigTemplate_syst1Down,kFALSE))
        
        TemplateName = "sigTemplatePdf_ggH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.sigTemplatePdf_ggH = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.SD,self.D),self.sigTempDataHist)
        TemplateName = "sigTemplatePdf_ggH_syst1Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.sigTemplatePdf_ggH_syst1Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.SD,self.D),self.sigTempDataHist_syst1Up)
        TemplateName = "sigTemplatePdf_ggH_syst1Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.sigTemplatePdf_ggH_syst1Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.SD,self.D),self.sigTempDataHist_syst1Down)
        TemplateName = "sigTemplatePdf_ggH_syst2Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.sigTemplatePdf_ggH_syst2Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.SD,self.D),self.sigTempDataHist_syst2Up)
        TemplateName = "sigTemplatePdf_ggH_syst2Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.sigTemplatePdf_ggH_syst2Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.SD,self.D),self.sigTempDataHist_syst2Down)
        
        #########################################################
        
        TemplateName = "sigTempDataHist_{0:.0f}_{1:.0f}{2}".format(self.channel,self.sqrts, self.appendHypType)
        self.sigTempDataHist_ALT = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(self.SD,self.D),ROOT.RooFit.Import(self.sigTemplate_ALT,kFALSE))
        TemplateName = "sigTempDataHist_syst1Up_{0:.0f}_{1:.0f}{2}".format(self.channel,self.sqrts, self.appendHypType)
        self.sigTempDataHist_ALT_syst1Up = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(self.SD,self.D),ROOT.RooFit.Import(self.sigTemplate_ALT_syst1Up,kFALSE))
        TemplateName = "sigTempDataHist_syst1Down_{0:.0f}_{1:.0f}{2}".format(self.channel,self.sqrts, self.appendHypType)
        self.sigTempDataHist_ALT_syst1Down = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(self.SD,self.D),ROOT.RooFit.Import(self.sigTemplate_ALT_syst1Down,kFALSE))
        TemplateName = "sigTempDataHist_syst2Up_{0:.0f}_{1:.0f}{2}".format(self.channel,self.sqrts, self.appendHypType)
        self.sigTempDataHist_ALT_syst2Up = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(self.SD,self.D),ROOT.RooFit.Import(self.sigTemplate_ALT_syst2Up,kFALSE))
        TemplateName = "sigTempDataHist_syst2Down_{0:.0f}_{1:.0f}{2}".format(self.channel,self.sqrts, self.appendHypType)
        self.sigTempDataHist_ALT_syst2Down = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(self.SD,self.D),ROOT.RooFit.Import(self.sigTemplate_ALT_syst2Down,kFALSE))
        
        TemplateName = "sigTemplatePdf_ggH{2}_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts, self.appendHypType)
        self.sigTemplatePdf_ggH_ALT = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.SD,self.D),self.sigTempDataHist_ALT)
        TemplateName = "sigTemplatePdf_ggH{2}_syst1Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts, self.appendHypType)
        self.sigTemplatePdf_ggH_ALT_syst1Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.SD,self.D),self.sigTempDataHist_ALT_syst1Up)
        TemplateName = "sigTemplatePdf_ggH{2}_syst1Down_{0:.0f}_{1:.0f}{2}".format(self.channel,self.sqrts, self.appendHypType)
        self.sigTemplatePdf_ggH_ALT_syst1Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.SD,self.D),self.sigTempDataHist_ALT_syst1Down)
        TemplateName = "sigTemplatePdf_ggH{2}_syst2Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts, self.appendHypType)
        self.sigTemplatePdf_ggH_ALT_syst2Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.SD,self.D),self.sigTempDataHist_ALT_syst2Up)
        TemplateName = "sigTemplatePdf_ggH{2}_syst2Down_{0:.0f}_{1:.0f}{2}".format(self.channel,self.sqrts, self.appendHypType)
        self.sigTemplatePdf_ggH_ALT_syst2Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.SD,self.D),self.sigTempDataHist_ALT_syst2Down)


        
        ###Shape systematics for signal
        self.funcList_ggH = ROOT.RooArgList()
        self.funcList_ggH_ALT = ROOT.RooArgList()

        if self.sigMorph:
            self.funcList_ggH.add(self.sigTemplatePdf_ggH)
            self.funcList_ggH.add(self.sigTemplatePdf_ggH_syst1Up)
            self.funcList_ggH.add(self.sigTemplatePdf_ggH_syst1Down)
            self.funcList_ggH.add(self.sigTemplatePdf_ggH_syst2Up)
            self.funcList_ggH.add(self.sigTemplatePdf_ggH_syst2Down)  
            
            self.funcList_ggH_ALT.add(self.sigTemplatePdf_ggH_ALT)
            self.funcList_ggH_ALT.add(self.sigTemplatePdf_ggH_ALT_syst1Up)
            self.funcList_ggH_ALT.add(self.sigTemplatePdf_ggH_ALT_syst1Down)
            self.funcList_ggH_ALT.add(self.sigTemplatePdf_ggH_ALT_syst2Up)
            self.funcList_ggH_ALT.add(self.sigTemplatePdf_ggH_ALT_syst2Down)
        else:
            self.funcList_ggH.add(self.sigTemplatePdf_ggH)
            self.funcList_ggH_ALT.add(self.sigTemplatePdf_ggH_ALT)
        


        #############################################################################
        self.morphVarListSig1 = ROOT.RooArgList()
        morphSigVarName = "CMS_zz4l_smd_leptScale_sig_{0:.0f}".format(self.channel)
        self.syst1MorphSig = ROOT.RooRealVar(morphSigVarName,morphSigVarName,0,-20,20)
        morphSigVarName = "CMS_zz4l_smd_leptResol_sig_{0:.0f}".format(self.channel)
        self.syst2MorphSig = ROOT.RooRealVar(morphSigVarName,morphSigVarName,0,-20,20)
        if self.sigMorph:
            self.syst1MorphSig.setConstant(False)
            self.syst2MorphSig.setConstant(False)
            self.morphVarListSig1.add(self.syst1MorphSig)
            ### just one morphing for all signal processes (fully correlated systs)
            ### self.morphVarListSig2.add(self.syst2MorphSig)
            self.morphVarListSig1.add(self.syst2MorphSig)

        else:
            self.syst1MorphSig.setConstant(True)
            self.syst2MorphSig.setConstant(True)
                
            
        TemplateName = "sigTemplateMorphPdf_ggH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.sigTemplateMorphPdf_ggH = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,self.SD,self.D,False,self.funcList_ggH,self.morphVarListSig1,1.0,1)
        
        TemplateName = "sigTemplateMorphPdf_ggH{2}_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts,self.appendHypType)
        self.sigTemplateMorphPdf_ggH_ALT = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,self.SD,self.D,False,self.funcList_ggH_ALT,self.morphVarListSig1,1.0,1)
        
        self.sigCB2d_ggH = self.signalCB_ggH.Clone("sigCB2d_ggH")
        self.sigCB2d_ggH_ALT = self.signalCB_ggH.Clone("sigCB2d_ggH{0}".format(self.appendHypType))

        ## ----------------- SuperKD 2D BACKGROUND SHAPES --------------- ##
        templateBkgName = "{0}/Dbackground_qqZZ_{1}.root".format(self.templateDir ,self.appendName)
        self.bkgTempFile = ROOT.TFile(templateBkgName)
        self.bkgTemplate = self.bkgTempFile.Get("h_superDpsD")
        self.bkgTemplateqqZZ = self.bkgTempFile.Get("h_superDpsD")
        

        TemplateName = "bkgTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.bkgTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(self.SD,self.D),ROOT.RooFit.Import(self.bkgTemplate,kFALSE))
        TemplateName = "bkgTemplatePdf_qqzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.bkgTemplatePdf_qqzz = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.SD,self.D),self.bkgTempDataHist)


        if self.reflectZX:
            
            TemplateName = "zjetsbkgTempDataHist_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.bkgTempDataHist_zjetsUp = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(self.SD,self.D),ROOT.RooFit.Import(self.bkgTemplate,kFALSE))
            TemplateName = "bkgTemplatePdf_zjets_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.bkgTemplatePdf_zjets_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.SD,self.D),self.bkgTempDataHist_zjetsUp)
        

            templateBkgName = "{0}/Dbackground_ZJetsCR_AllChans.root".format(self.templateDir)
            self.zjetsTempFile = ROOT.TFile(templateBkgName)
            self.zjetsTemplate = self.zjetsTempFile.Get("h_superDpsD")
            TemplateName = "zjetsbkgTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.bkgTempDataHist_zjets = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(self.SD,self.D),ROOT.RooFit.Import(self.zjetsTemplate,kFALSE))
            TemplateName = "bkgTemplatePdf_zjets_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.bkgTemplatePdf_zjets = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.SD,self.D),self.bkgTempDataHist_zjets)

            self.bkgTempDataHist_zjetsDown1 = self.reflectSystematics(self.zjetsTemplate,self.bkgTemplateqqZZ)
            TemplateName = "zjetsbkgTempDataHist_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.bkgTempDataHist_zjetsDown = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(self.SD,self.D),ROOT.RooFit.Import(self.bkgTempDataHist_zjetsDown1,kFALSE))
            TemplateName = "bkgTemplatePdf_zjets_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.bkgTemplatePdf_zjets_Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.SD,self.D),self.bkgTempDataHist_zjetsDown)
        

        else:

            #Nominal Z+X
            TemplateName = "zjetsbkgTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.bkgTempDataHist_zjets = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(self.SD,self.D),ROOT.RooFit.Import(self.bkgTemplate,kFALSE))
            TemplateName = "bkgTemplatePdf_zjets_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.bkgTemplatePdf_zjets = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.SD,self.D),self.bkgTempDataHist_zjets)
            

            #Up Z+X
            TemplateName = "zjetsbkgTempDataHist_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.bkgTempDataHist_zjetsUp = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(self.SD,self.D),ROOT.RooFit.Import(self.bkgTemplate,kFALSE))
            TemplateName = "bkgTemplatePdf_zjets_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.bkgTemplatePdf_zjets_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.SD,self.D),self.bkgTempDataHist_zjetsUp)
        

            #Down Z+X
            templateBkgName = "{0}/Dbackground_ZJetsCR_AllChans.root".format(self.templateDir)
            self.zjetsTempFile = ROOT.TFile(templateBkgName)
            self.zjetsTemplateDown = self.zjetsTempFile.Get("h_superDpsD")
            TemplateName = "zjetsbkgTempDataHist_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.bkgTempDataHist_zjetsDown = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(self.SD,self.D),ROOT.RooFit.Import(self.zjetsTemplateDown,kFALSE))
            TemplateName = "bkgTemplatePdf_zjets_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.bkgTemplatePdf_zjets_Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.SD,self.D),self.bkgTempDataHist_zjetsDown)




        templateggBkgName = "{0}/Dbackground_ggZZ_{1}.root".format(self.templateDir ,self.appendName)
        self.ggbkgTempFile = ROOT.TFile(templateggBkgName)
        self.ggbkgTemplate = self.ggbkgTempFile.Get("h_superDpsD")
        TemplateName = "ggbkgTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.ggbkgTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(self.SD,self.D),ROOT.RooFit.Import(self.ggbkgTemplate,kFALSE))
        TemplateName = "bkgTemplatePdf_ggzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.bkgTemplatePdf_ggzz = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(self.SD,self.D),self.ggbkgTempDataHist)


        #################################
        
        self.funcList_zjets = ROOT.RooArgList()
        morphBkgVarName =  "CMS_zz4l_smd_zjets_bkg_{0:.0f}".format(self.channel)
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
        self.bkgTemplateMorphPdf_qqzz = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,self.SD,self.D,False,ROOT.RooArgList(self.bkgTemplatePdf_qqzz),ROOT.RooArgList(),1.0,1)
        TemplateName = "bkgTemplateMorphPdf_ggzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.bkgTemplateMorphPdf_ggzz = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,self.SD,self.D,False,ROOT.RooArgList(self.bkgTemplatePdf_ggzz),ROOT.RooArgList(),1.0,1)
        
        TemplateName = "bkgTemplateMorphPdf_zjets_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.bkgTemplateMorphPdf_zjets = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,self.SD,self.D,False,self.funcList_zjets,self.morphVarListBkg,1.0,1)
        
        self.bkg1d_qqzz = self.bkg_qqzz.Clone("bkg1d_qqzz")
        self.bkg1d_ggzz = self.bkg_ggzz.Clone("bkg1d_ggzz")
        self.bkg1d_zjets = self.bkg_zjets.Clone("bkg1d_zjets")

            
    ## --------------------------- WORKSPACE -------------------------- ##
    def writeWorkspaceSuperKD(self):
        
        endsInP5 = False
        tmpMH = self.mH
        if (math.fabs(math.floor(tmpMH)-self.mH) > 0.001): endsInP5 = True
        if (self.DEBUG): print "ENDS IN P5  ",endsInP5
        
        self.name_Shape = ""
        self.name_ShapeWS = ""
        self.name_ShapeWS2 = ""

        if (endsInP5):
            self.name_Shape = "{0}/HCG/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts)
        else:
            self.name_Shape = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts)
            
        if (endsInP5):
            self.name_ShapeWS = "{0}/HCG/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts)
        else:
            self.name_ShapeWS = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts)
            
        self.name_ShapeWS2 = "hzz4l_{0}S_{1:.0f}TeV.input.root".format(self.appendName,self.sqrts)
            
        if(self.DEBUG): print self.name_Shape,"  ",self.name_ShapeWS2
        
        #Workspace
        self.w = ROOT.RooWorkspace("w","w")
        #Class code
        self.w.importClassCode(RooqqZZPdf_v2.Class(),True)
        self.w.importClassCode(RooggZZPdf_v2.Class(),True)
        self.w.importClassCode(RooRelBWUFParam.Class(),True)
        self.w.importClassCode(RooDoubleCB.Class(),True)
        self.w.importClassCode(RooFormulaVar.Class(),True)
        #Data
        getattr(self.w,'import')(self.data_obs,ROOT.RooFit.Rename("data_obs"))
        
        self.sigCB2d_ggH.SetNameTitle("ORIGggH","ORIGggH")
        
        getattr(self.w,'import')(self.sigCB2d_ggH, ROOT.RooFit.RecycleConflictNodes())
        
        self.sigTemplateMorphPdf_ggH.SetNameTitle("ggH","ggH")
        getattr(self.w,'import')(self.sigTemplateMorphPdf_ggH, ROOT.RooFit.RecycleConflictNodes())
        #save syst templates individually
        systTempName=("ggHCMS_zz4l_leptScale_sig_{0}_{1:.0f}_Up").format(self.channel,self.sqrts)
        self.sigTemplatePdf_ggH_syst1Up.SetNameTitle(systTempName,systTempName)
        systTempName=("ggHCMS_zz4l_leptScale_sig_{0}_{1:.0f}_Down").format(self.channel,self.sqrts)
        self.sigTemplatePdf_ggH_syst1Down.SetNameTitle(systTempName,systTempName)
        
        self.sigCB2d_ggH_ALT.SetNameTitle("ORIGggH{0}".format(self.appendHypType),"ggH{0}".format(self.appendHypType))
        self.sigTemplateMorphPdf_ggH_ALT.SetNameTitle("ggH{0}".format(self.appendHypType),"ggH{0}".format(self.appendHypType))
        getattr(self.w,'import')(self.sigCB2d_ggH_ALT, ROOT.RooFit.RecycleConflictNodes())
        getattr(self.w,'import')(self.sigTemplatePdf_ggH_ALT, ROOT.RooFit.RecycleConflictNodes())
        getattr(self.w,'import')(self.sigTemplateMorphPdf_ggH_ALT, ROOT.RooFit.RecycleConflictNodes())
        
        
        
        self.bkg1d_qqzz.SetNameTitle("ORIGbkg_qqzz","ORIGbkg_qqzz")
        self.bkg1d_ggzz.SetNameTitle("ORIGbkg_ggzz","ORIGbkg_ggzz")
        self.bkg1d_zjets.SetNameTitle("ORIGbkg_zjets","ORIGbkg_zjets")
        
        self.bkgTemplateMorphPdf_qqzz.SetNameTitle("bkg2d_qqzz","bkg2d_qqzz")
        self.bkgTemplateMorphPdf_ggzz.SetNameTitle("bkg2d_ggzz","bkg2d_ggzz")
        
        self.bkgTemplateMorphPdf_zjets.SetNameTitle("bkg2d_zjets","bkg2d_zjets")
        
        getattr(self.w,'import')(self.bkgTemplateMorphPdf_qqzz,ROOT.RooFit.RecycleConflictNodes())
        getattr(self.w,'import')(self.bkgTemplateMorphPdf_ggzz,ROOT.RooFit.RecycleConflictNodes())
        getattr(self.w,'import')(self.bkgTemplateMorphPdf_zjets,ROOT.RooFit.RecycleConflictNodes())
        
        ### save signal rates
        getattr(self.w,'import')(self.rfvSigRate_ggH, ROOT.RooFit.RecycleConflictNodes())
        self.rfvSigRate_ggH_ALT = ROOT.RooFormulaVar(self.rfvSigRate_ggH,"ggH{0}_norm".format(self.appendHypType))
#        print '>>>>>> Compare Signal Rates: SM=',self.rfvSigRate_ggH.getVal()," ALT=",self.rfvSigRate_ggH_ALT.getVal()

        self.rrvJHUgen_ggH_ALT = ROOT.RooRealVar("jhuGen_ALT","jhuGen_ALT",1)
        if self.altHypothesis == 'gg0-':
            self.rrvJHUgen_ggH_ALT.setVal(float(self.inputs['jhuGen_0minus_yield']))
            self.sigRate_ggH_ALT_Shape = self.rrvJHUgen_ggH_ALT.getVal()*self.rrv_SMggH_ratio.getVal()
        elif self.altHypothesis == 'gg0h+':
            self.rrvJHUgen_ggH_ALT.setVal(float(self.inputs['jhuGen_0hplus_yield']))
            self.sigRate_ggH_ALT_Shape = self.rrvJHUgen_ggH_ALT.getVal()*self.rrv_SMggH_ratio.getVal()
#        elif self.altHypothesis == 'gg2m+': ### 0+ vs gg2m+ with or without production
#            self.rrvJHUgen_ggH_ALT.setVal(float(self.inputs['jhuGen_gg2mplus_yield']))
#            self.sigRate_ggH_ALT_Shape = self.rrvJHUgen_ggH_ALT.getVal()*self.rrv_SMggH_ratio.getVal()
#        elif self.altHypothesis == '2h+':
#            self.rrvJHUgen_ggH_ALT.setVal(float(self.inputs['jhuGen_2hplus_yield']))
#            self.sigRate_ggH_ALT_Shape = self.rrvJHUgen_ggH_ALT.getVal()*self.rrv_SMggH_ratio.getVal()
#        elif self.altHypothesis == '2h-': ### 0+ vs 2h-
#            self.rrvJHUgen_ggH_ALT.setVal(float(self.inputs['jhuGen_2hminus_yield']))
#            self.sigRate_ggH_ALT_Shape = self.rrvJHUgen_ggH_ALT.getVal()*self.rrv_SMggH_ratio.getVal()
#        elif self.altHypothesis == '2b+':
#            self.rrvJHUgen_ggH_ALT.setVal(float(self.inputs['jhuGen_2bplus_yield']))
#            self.sigRate_ggH_ALT_Shape = self.rrvJHUgen_ggH_ALT.getVal()*self.rrv_SMggH_ratio.getVal()
        else:
            self.rrvJHUgen_ggH_ALT.setVal(self.rates['ggH']*self.calcTotalYieldCorr(self.channel,self.sqrts,self.altHypothesis))
            self.sigRate_ggH_ALT_Shape = self.rrvJHUgen_ggH_ALT.getVal()
                    

        self.rfvSigRate_ggH_ALT = ROOT.RooFormulaVar("ggH{0}_norm".format(self.appendHypType),"@0",ROOT.RooArgList(self.one))
        print '>>>>>> Compare signal rates: STD=',self.rfvSigRate_ggH.getVal(),"   ALT=",self.rfvSigRate_ggH_ALT.getVal()
        print '>>>>>> Compare signal rates: STD=',self.rates['ggH'],"   ALT=",self.sigRate_ggH_ALT_Shape
        self.rates['ggH{0}'.format(self.appendHypType)] = self.sigRate_ggH_ALT_Shape
        getattr(self.w,'import')(self.rfvSigRate_ggH_ALT, ROOT.RooFit.RecycleConflictNodes())
                    
        self.w.writeToFile(self.name_ShapeWS)
        
                
        
        
    def writeWorkspaceUnfoldedSuperKD(self):
            
        endsInP5 = False
        tmpMH = self.mH
        if (math.fabs(math.floor(tmpMH)-self.mH) > 0.001): endsInP5 = True
        if (self.DEBUG): print "ENDS IN P5  ",endsInP5
        
        self.name_Shape = ""
        self.name_ShapeWS = ""
        self.name_ShapeWS2 = ""
        
        if (endsInP5):
            self.name_Shape = "{0}/HCG/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts)
        else:
            self.name_Shape = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts)
            
        if (endsInP5):
            self.name_ShapeWS = "{0}/HCG/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts)
        else:
            self.name_ShapeWS = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts)
            
        self.name_ShapeWS2 = "hzz4l_{0}S_{1:.0f}TeV.input.root".format(self.appendName,self.sqrts)
        
        if(self.DEBUG): print self.name_Shape,"  ",self.name_ShapeWS2


        workspaceFile = ROOT.TFile(self.name_ShapeWS,"RECREATE")
        workspaceFile.cd()
        
        #SM
        theYield = self.rates['ggH']
        hist = self.unfoldedHist(self.sigTemplate,'ggH',theYield)
        hist.Write()
        hist = self.unfoldedHist(self.sigTemplate_syst1Up,'ggH_CMS_zz4l_smd_leptScale_sig_{0}Up'.format(self.channel),theYield)
        hist.Write()
        hist = self.unfoldedHist(self.sigTemplate_syst1Down,'ggH_CMS_zz4l_smd_leptScale_sig_{0}Down'.format(self.channel),theYield)
        hist.Write()
        hist = self.unfoldedHist(self.sigTemplate_syst2Up,'ggH_CMS_zz4l_smd_leptResol_sig_{0}Up'.format(self.channel),theYield)
        hist.Write()
        hist = self.unfoldedHist(self.sigTemplate_syst2Down,'ggH_CMS_zz4l_smd_leptResol_sig_{0}Down'.format(self.channel),theYield)
        hist.Write()

        #ALT
        self.rrvJHUgen_ggH_ALT = ROOT.RooRealVar("jhuGen_ALT","jhuGen_ALT",1)
        if self.altHypothesis == 'gg0-':
            self.rrvJHUgen_ggH_ALT.setVal(float(self.inputs['jhuGen_0minus_yield']))
            self.sigRate_ggH_ALT_Shape = self.rrvJHUgen_ggH_ALT.getVal()*self.rrv_SMggH_ratio.getVal()
        elif self.altHypothesis == 'gg0h+':
            self.rrvJHUgen_ggH_ALT.setVal(float(self.inputs['jhuGen_0hplus_yield']))
            self.sigRate_ggH_ALT_Shape = self.rrvJHUgen_ggH_ALT.getVal()*self.rrv_SMggH_ratio.getVal()
#        elif self.altHypothesis == 'gg2m+': ### 0+ vs gg2m+ with or without production
#            self.rrvJHUgen_ggH_ALT.setVal(float(self.inputs['jhuGen_gg2mplus_yield']))
#            self.sigRate_ggH_ALT_Shape = self.rrvJHUgen_ggH_ALT.getVal()*self.rrv_SMggH_ratio.getVal()
#        elif self.altHypothesis == '2h+':
#            self.rrvJHUgen_ggH_ALT.setVal(float(self.inputs['jhuGen_2hplus_yield']))
#            self.sigRate_ggH_ALT_Shape = self.rrvJHUgen_ggH_ALT.getVal()*self.rrv_SMggH_ratio.getVal()
#        elif self.altHypothesis == '2h-': ### 0+ vs 2h-
#            self.rrvJHUgen_ggH_ALT.setVal(float(self.inputs['jhuGen_2hminus_yield']))
#            self.sigRate_ggH_ALT_Shape = self.rrvJHUgen_ggH_ALT.getVal()*self.rrv_SMggH_ratio.getVal()
#        elif self.altHypothesis == '2b+':
#            self.rrvJHUgen_ggH_ALT.setVal(float(self.inputs['jhuGen_2bplus_yield']))
#            self.sigRate_ggH_ALT_Shape = self.rrvJHUgen_ggH_ALT.getVal()*self.rrv_SMggH_ratio.getVal()
        else:
            self.rrvJHUgen_ggH_ALT.setVal(self.rates['ggH']*self.calcTotalYieldCorr(self.channel,self.sqrts,self.altHypothesis))
            self.sigRate_ggH_ALT_Shape = self.rrvJHUgen_ggH_ALT.getVal()

        self.rfvSigRate_ggH_ALT = ROOT.RooFormulaVar("ggH{0}_norm".format(self.appendHypType),"@0",ROOT.RooArgList(self.one))
        print '>>>>>> Compare signal rates: STD=',self.rfvSigRate_ggH.getVal(),"   ALT=",self.rfvSigRate_ggH_ALT.getVal()
        print '>>>>>> Compare signal rates: STD=',self.rates['ggH'],"   ALT=",self.sigRate_ggH_ALT_Shape
        self.rates['ggH{0}'.format(self.appendHypType)] = self.sigRate_ggH_ALT_Shape

        theYield = self.rates['ggH{0}'.format(self.appendHypType)]
        hist = self.unfoldedHist(self.sigTemplate_ALT,'ggH{0}'.format(self.appendHypType),theYield)
        hist.Write()
        hist = self.unfoldedHist(self.sigTemplate_ALT_syst1Up,'ggH{0}_CMS_zz4l_smd_leptScale_sig_{1}Up'.format(self.appendHypType,self.channel),theYield)
        hist.Write()
        hist = self.unfoldedHist(self.sigTemplate_ALT_syst1Down,'ggH{0}_CMS_zz4l_smd_leptScale_sig_{1}Down'.format(self.appendHypType,self.channel),theYield)
        hist.Write()
        hist = self.unfoldedHist(self.sigTemplate_ALT_syst2Up,'ggH{0}_CMS_zz4l_smd_leptResol_sig_{1}Up'.format(self.appendHypType,self.channel),theYield)
        hist.Write()
        hist = self.unfoldedHist(self.sigTemplate_ALT_syst2Down,'ggH{0}_CMS_zz4l_smd_leptResol_sig_{1}Down'.format(self.appendHypType,self.channel),theYield)
        hist.Write()
        
        #qqZZ
        theYield = self.rates['qqZZ']
        hist = self.unfoldedHist(self.bkgTemplate,'bkg2d_qqzz',theYield)
        hist.Write()
        
        #ggZZ
        theYield = self.rates['ggZZ']
        hist = self.unfoldedHist(self.ggbkgTemplate,'bkg2d_ggzz',theYield)
        hist.Write()
        
        #zjets
        theYield = self.rates['zjets']
        if self.reflectZX:
            hist = self.unfoldedHist(self.zjetsTemplate,'bkg2d_zjets',theYield)
            hist.Write()
            hist = self.unfoldedHist(self.bkgTemplate,'bkg2d_zjets_CMS_zz4l_smd_zjets_bkg_{0:.0f}Down'.format(self.channel),theYield)
            hist.Write()
            hist = self.unfoldedHist(self.bkgTemplate,'bkg2d_zjets_CMS_zz4l_smd_zjets_bkg_{0:.0f}Up'.format(self.channel),theYield)#for now
            hist.Write()
        else:
            hist = self.unfoldedHist(self.bkgTemplate,'bkg2d_zjets',theYield)
            hist.Write()
            hist = self.unfoldedHist(self.zjetsTemplateDown,'bkg2d_zjets_CMS_zz4l_smd_zjets_bkg_{0:.0f}Down'.format(self.channel),theYield)
            hist.Write()
            hist = self.unfoldedHist(self.bkgTemplate,'bkg2d_zjets_CMS_zz4l_smd_zjets_bkg_{0:.0f}Up'.format(self.channel),theYield)#for now
            hist.Write()

        #data
        histData2D = self.sigTemplate.Clone()
        for x in xrange(1,histData2D.GetNbinsX()+1):
            for y in xrange(1,histData2D.GetNbinsY()+1):
                histData2D.SetBinContent(x,y,0)

        myfile = ROOT.TFile(self.dataFileName,"READ")
        myTree = myfile.Get(self.dataTreeName)
        
        for i in range( myTree.GetEntries() ):
            myTree.GetEntry(i);

            if myTree.CMS_zz4l_mass < self.inputs['low_M'] or myTree.CMS_zz4l_mass > self.inputs['high_M']: continue

            sd = getattr(myTree,self.superDiscVarName)
            pd = getattr(myTree,self.discVarName)
            print sd, pd
            histData2D.Fill(sd,pd)

        workspaceFile.cd()
        hist = self.unfoldedHist(histData2D,'data_obs',-1)
        hist.Write()
            
        workspaceFile.Close()


    def unfoldedHist(self,hist,histName,norm):

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

        figs = self.outputDir+"/figs"
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
        canv.SaveAs(figs+"/"+histName+"_"+self.appendName+".eps")
        canv.SaveAs(figs+"/"+histName+"_"+self.appendName+".png")

                    
        return hist1D




    def prepareDatacardSuperKD(self):

        name_Shape = ""
        name_ShapeWS = ""

        ## Write Datacards
        if (self.endsInP5(self.mH)):
            name_Shape = "{0}/HCG/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV{4}.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts,self.appendHypType)
        else:
            name_Shape = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV{4}.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts,self.appendHypType)
        fo = open( name_Shape, "wb")
        self.WriteDatacardSuperKD(fo, self.name_ShapeWS2, self.rates, self.data_obs.numEntries())
        self.systematics.WriteSystematics(fo, self.inputs)
        self.systematics.WriteShapeSystematics(fo,self.inputs)
        self.systematics.WriteSuperKDShapeSystematics(fo,self.inputs)
        
        fo.close()


        
    def WriteDatacardSuperKD(self,file,nameWS,theRates,obsEvents):

        numberSig = self.numberOfSigChan(self.inputs)
        numberBg  = self.numberOfBgChan(self.inputs)
        
        print "DEBUG: Number of signals={0} and backgrounds={1}. IsAltHyp={2} ".format(numberSig,numberBg,self.isAltSig)

        file.write("imax 1\n")
        file.write("jmax {0}\n".format(numberSig+numberBg-1))
        file.write("kmax *\n")
        file.write("------------\n")

        if self.inputs['unfold']:
            file.write("shapes * * {0} $PROCESS $PROCESS_$SYSTEMATIC \n".format(nameWS))
            file.write("shapes data_obs a{0} {1} data_obs \n".format(self.channel,nameWS))
        else:
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
            channelList=['ggH','ggH','qqZZ','ggZZ','zjets']
            channelName2D=['ggH','ggH{0}'.format(self.appendHypType),'bkg2d_qqzz','bkg2d_ggzz','bkg2d_zjets']

          
        for chan in channelList:
            if self.inputs[chan] or self.inputs['all']:
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
        j=0
        for chan in channelList:
            if self.inputs[chan] or self.inputs["all"]:
                if channelName2D[j].startswith('ggH{0}'.format(self.appendHypType)):
                    theYield = theRates['ggH{0}'.format(self.appendHypType)]
                    file.write("{0:.4f} ".format(theYield))
		elif self.killBackground and not channelName2D[j].startswith('ggH'):
		    file.write("1e-5 ")
                else:
                    file.write("{0:.4f} ".format(theRates[chan]))
            j+=1
            
        file.write("\n")
        file.write("------------\n")



    def calcTotalYieldCorr(self, channel, sqrts, spinHypCode):
        
        ### the parameters are the fraction of 2e2mu with respect to the total
        ### we modify the total assuming that the 2e2mu yield is constant
        ### among different spin hypotheses
                
        r = 1.0 

        #[4mu, 4e, 2e2mu]
        #corrFactor_0minus = [0.912885,0.858691,1.12188]
        #corrFactor_gg2plusm = [0.900717,0.865112,1.12826]
        #corrFactor_0hplus = [1.12826,0.947087,1.06477]
        #corrFactor_1plus  = [0.968514,0.882397,1.07114]
        #corrFactor_1minus = [0.939779,0.854791,1.1036]
        #corrFactor_qq2plusm = [0.884544,0.861437,1.1417]

        corrFactor_0minus_7 = [0.94178,0.847481,1.13082]
        corrFactor_0hplus_7 = [0.926585,0.934054,1.09296]
        corrFactor_1plus_7  = [0.973026,0.907252,1.09667]
        corrFactor_1minus_7 = [0.931127,0.89238,1.14357]
        corrFactor_qq2plusm_7 = [0.916151,0.854769,1.15253]
        corrFactor_gg2plusm_7 = [0.911745,0.866069,1.12853]
        corrFactor_2hplus_7 = [0.927213,0.918012,1.09863]
        corrFactor_2hminus_7 = [0.945386,0.903211,1.07699]
        corrFactor_2bplus_7 = [0.891317,0.869832,1.14833]
        
        corrFactor_0minus_8 = [0.914991,0.854913,1.1167]
        corrFactor_0hplus_8 = [0.938505,0.945478,1.06486]
        corrFactor_1plus_8  = [0.964351,0.874151,1.06928]
        corrFactor_1minus_8 = [0.937464,0.846397,1.09798]
        corrFactor_qq2plusm_8 = [0.889664,0.867378,1.12889]
        corrFactor_gg2plusm_8 = [0.911387,0.864769,1.11862]
        corrFactor_2hplus_8 = [0.920798,0.928968,1.08486]
        corrFactor_2hminus_8 = [0.946448,0.927209,1.06996]
        corrFactor_2bplus_8 = [0.876706,0.886499,1.13396]


        if sqrts == 7:

            if spinHypCode == 'SM':
                r = 1.0
            elif spinHypCode == 'gg0-':
                r = 1#corrFactor_0minus[channel-1]
            elif spinHypCode == 'gg0h+':
                r = 1#corrFactor_0hplus[channel-1]
            elif spinHypCode == 'qq1-' or spinHypCode == 'PI1-':
                r = corrFactor_1minus_7[channel-1]
            elif spinHypCode == 'qq1+' or spinHypCode == 'PI1+':
                r = corrFactor_1plus_7[channel-1]
            elif spinHypCode == 'gg2m+' or spinHypCode == 'PI2+':
                r = corrFactor_gg2plusm_7[channel-1]
            elif spinHypCode == 'qq2m+':
                r = corrFactor_qq2plusm_7[channel-1]
            elif spinHypCode == '2h+':
                r = corrFactor_2hplus_7[channel-1]
            elif spinHypCode == '2h-':
                r = corrFactor_2hminus_7[channel-1]
            elif spinHypCode == '2b+':
                r = corrFactor_2bplus_7[channel-1]
            else:
                raise RuntimeError,'spinHypCode '+str(spinHypCode)+' is unknown! Please choose: SM, gg0-, gg0h+, qq1-, qq1+, gg2m+, qq2m+, 2h+, 2h-, 2b+'

        elif sqrts == 8:            

            if spinHypCode == 'SM':
                r = 1.0
            elif spinHypCode == 'gg0-':
                r = 1#corrFactor_0minus[channel-1]
            elif spinHypCode == 'gg0h+':
                r = 1#corrFactor_0hplus[channel-1]
            elif spinHypCode == 'qq1-' or spinHypCode == 'PI1-':
                r = corrFactor_1minus_8[channel-1]
            elif spinHypCode == 'qq1+' or spinHypCode == 'PI1+':
                r = corrFactor_1plus_8[channel-1]
            elif spinHypCode == 'gg2m+' or spinHypCode == 'PI2+':
                r = corrFactor_gg2plusm_8[channel-1]
            elif spinHypCode == 'qq2m+':
                r = corrFactor_qq2plusm_8[channel-1]
            elif spinHypCode == '2h+':
                r = corrFactor_2hplus_8[channel-1]
            elif spinHypCode == '2h-':
                r = corrFactor_2hminus_8[channel-1]
            elif spinHypCode == '2b+':
                r = corrFactor_2bplus_8[channel-1]
            else:
                raise RuntimeError,'spinHypCode '+str(spinHypCode)+' is unknown! Please choose: SM, gg0-, gg0h+, qq1-, qq1+, gg2m+, qq2m+, 2h+, 2h-, 2b+'
        else:
            raise RuntimeError,'spinHypCorrFactor: Unknown sqrts! '+str(sqrts)

        
        print '>>>>>> Acceptance*BR correction for ALT hypothesis: ',r
        return r



                                                                    
    def reflectSystematics(self,nomShape,altShape):

        if(nomShape.GetNbinsX()!=altShape.GetNbinsX() or nomShape.GetNbinsY()!=altShape.GetNbinsY()):
            print nomShape.GetNbinsX(),altShape.GetNbinsX(),nomShape.GetNbinsY(),altShape.GetNbinsY()
            print "AHHHHHHHHHHH, templates don't have the same binning!!!!"
            return 0
        
        newAltShape = ROOT.TH2D(altShape)
        
        for x in range(1,nomShape.GetNbinsX()):
            for y in range(1,nomShape.GetNbinsY()):
                delta=altShape.GetBinContent(x,y)-nomShape.GetBinContent(x,y)
                newAltShape.SetBinContent(x,y,nomShape.GetBinContent(x,y)-delta)
                if(newAltShape.GetBinContent(x,y)<0.0):
                    newAltShape.SetBinContent(x,y,0.0)
                    # done with loop over y bins
                    #done with loop over x bins
                    
        newAltShape.Scale(1.0/newAltShape.Integral())
                    
        #check that no bins are zero
                    
        for x in range(1,newAltShape.GetNbinsX()):
            for y in range(1,newAltShape.GetNbinsY()):
                if(newAltShape.GetBinContent(x,y)<0.000001):
                    newAltShape.SetBinContent(x,y,0.000001)
                    # done with loop over y bins
                    #done with loop over x bins
                    
        newAltShape.Scale(1.0/newAltShape.Integral())
                    
        return newAltShape
        if self.inputs["all"]:
            print "ALL CHANNELS --- KDCLASS --- WriteDatacardKD\n"
            channelList=['ggH','ggH','qqZZ','ggZZ','zjets']
            channelName2D=['ggH','ggH{0}'.format(self.appendHypType),'bkg2d_qqzz','bkg2d_ggzz','bkg2d_zjets']

          
        for chan in channelList:
            if self.inputs[chan] or self.inputs['all']:
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
        j=0
        for chan in channelList:
            if self.inputs[chan] or self.inputs["all"]:
                if channelName2D[j].startswith('ggH{0}'.format(self.appendHypType)):
                    theYield = theRates['ggH{0}'.format(self.appendHypType)]
                    file.write("{0:.4f} ".format(theYield))
                else:
                    file.write("{0:.4f} ".format(theRates[chan]))
            j+=1
            
        file.write("\n")
        file.write("------------\n")



    def calcTotalYieldCorr(self, channel, sqrts, spinHypCode):
        
        ### the parameters are the fraction of 2e2mu with respect to the total
        ### we modify the total assuming that the 2e2mu yield is constant
        ### among different spin hypotheses
                
        r = 1.0 

        #[4mu, 4e, 2e2mu]
        #corrFactor_0minus = [0.912885,0.858691,1.12188]
        #corrFactor_gg2plusm = [0.900717,0.865112,1.12826]
        #corrFactor_0hplus = [1.12826,0.947087,1.06477]
        #corrFactor_1plus  = [0.968514,0.882397,1.07114]
        #corrFactor_1minus = [0.939779,0.854791,1.1036]
        #corrFactor_qq2plusm = [0.884544,0.861437,1.1417]

        corrFactor_0minus_7 = [0.94178,0.847481,1.13082]
        corrFactor_0hplus_7 = [0.926585,0.934054,1.09296]
        corrFactor_1plus_7  = [0.973026,0.907252,1.09667]
        corrFactor_1minus_7 = [0.931127,0.89238,1.14357]
        corrFactor_qq2plusm_7 = [0.916151,0.854769,1.15253]
        corrFactor_gg2plusm_7 = [0.911745,0.866069,1.12853]
        corrFactor_2hplus_7 = [0.927213,0.918012,1.09863]
        corrFactor_2hminus_7 = [0.945386,0.903211,1.07699]
        corrFactor_2bplus_7 = [0.891317,0.869832,1.14833]
        
        corrFactor_0minus_8 = [0.914991,0.854913,1.1167]
        corrFactor_0hplus_8 = [0.938505,0.945478,1.06486]
        corrFactor_1plus_8  = [0.964351,0.874151,1.06928]
        corrFactor_1minus_8 = [0.937464,0.846397,1.09798]
        corrFactor_qq2plusm_8 = [0.889664,0.867378,1.12889]
        corrFactor_gg2plusm_8 = [0.911387,0.864769,1.11862]
        corrFactor_2hplus_8 = [0.920798,0.928968,1.08486]
        corrFactor_2hminus_8 = [0.946448,0.927209,1.06996]
        corrFactor_2bplus_8 = [0.876706,0.886499,1.13396]


        if sqrts == 7:

            if spinHypCode == 'SM':
                r = 1.0
            elif spinHypCode == 'gg0-':
                r = 1#corrFactor_0minus[channel-1]
            elif spinHypCode == 'gg0h+':
                r = 1#corrFactor_0hplus[channel-1]
            elif spinHypCode == 'qq1-' or spinHypCode == 'PI1-':
                r = corrFactor_1minus_7[channel-1]
            elif spinHypCode == 'qq1+' or spinHypCode == 'PI1+':
                r = corrFactor_1plus_7[channel-1]
            elif spinHypCode == 'gg2m+' or spinHypCode == 'PI2+':
                r = corrFactor_gg2plusm_7[channel-1]
            elif spinHypCode == 'qq2m+':
                r = corrFactor_qq2plusm_7[channel-1]
            elif spinHypCode == '2h+':
                r = corrFactor_2hplus_7[channel-1]
            elif spinHypCode == '2h-':
                r = corrFactor_2hminus_7[channel-1]
            elif spinHypCode == '2b+':
                r = corrFactor_2bplus_7[channel-1]
            else:
                raise RuntimeError,'spinHypCode '+str(spinHypCode)+' is unknown! Please choose: SM, gg0-, gg0h+, qq1-, qq1+, gg2m+, qq2m+, 2h+, 2h-, 2b+'

        elif sqrts == 8:            

            if spinHypCode == 'SM':
                r = 1.0
            elif spinHypCode == 'gg0-':
                r = 1#corrFactor_0minus[channel-1]
            elif spinHypCode == 'gg0h+':
                r = 1#corrFactor_0hplus[channel-1]
            elif spinHypCode == 'qq1-' or spinHypCode == 'PI1-':
                r = corrFactor_1minus_8[channel-1]
            elif spinHypCode == 'qq1+' or spinHypCode == 'PI1+':
                r = corrFactor_1plus_8[channel-1]
            elif spinHypCode == 'gg2m+' or spinHypCode == 'PI2+':
                r = corrFactor_gg2plusm_8[channel-1]
            elif spinHypCode == 'qq2m+':
                r = corrFactor_qq2plusm_8[channel-1]
            elif spinHypCode == '2h+':
                r = corrFactor_2hplus_8[channel-1]
            elif spinHypCode == '2h-':
                r = corrFactor_2hminus_8[channel-1]
            elif spinHypCode == '2b+':
                r = corrFactor_2bplus_8[channel-1]
            else:
                raise RuntimeError,'spinHypCode '+str(spinHypCode)+' is unknown! Please choose: SM, gg0-, gg0h+, qq1-, qq1+, gg2m+, qq2m+, 2h+, 2h-, 2b+'
        else:
            raise RuntimeError,'spinHypCorrFactor: Unknown sqrts! '+str(sqrts)

        
        print '>>>>>> Acceptance*BR correction for ALT hypothesis: ',r
        return r



                                                                    
    def reflectSystematics(self,nomShape,altShape):

        if(nomShape.GetNbinsX()!=altShape.GetNbinsX() or nomShape.GetNbinsY()!=altShape.GetNbinsY()):
            print nomShape.GetNbinsX(),altShape.GetNbinsX(),nomShape.GetNbinsY(),altShape.GetNbinsY()
            print "AHHHHHHHHHHH, templates don't have the same binning!!!!"
            return 0
        
        newAltShape = ROOT.TH2D(altShape)
        
        for x in range(1,nomShape.GetNbinsX()):
            for y in range(1,nomShape.GetNbinsY()):
                delta=altShape.GetBinContent(x,y)-nomShape.GetBinContent(x,y)
                newAltShape.SetBinContent(x,y,nomShape.GetBinContent(x,y)-delta)
                if(newAltShape.GetBinContent(x,y)<0.0):
                    newAltShape.SetBinContent(x,y,0.0)
                    # done with loop over y bins
                    #done with loop over x bins
                    
        newAltShape.Scale(1.0/newAltShape.Integral())
                    
        #check that no bins are zero
                    
        for x in range(1,newAltShape.GetNbinsX()):
            for y in range(1,newAltShape.GetNbinsY()):
                if(newAltShape.GetBinContent(x,y)<0.000001):
                    newAltShape.SetBinContent(x,y,0.000001)
                    # done with loop over y bins
                    #done with loop over x bins
                    
        newAltShape.Scale(1.0/newAltShape.Integral())
                    
        return newAltShape


