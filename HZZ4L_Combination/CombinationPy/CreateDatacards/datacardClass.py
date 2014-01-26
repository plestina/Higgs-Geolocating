#! /usr/bin/env python
import sys
import os
import re
import math
from scipy.special import erf
from ROOT import *
import ROOT
from array import array
from systematicsClass import *
from lib.util.Logger import *
## ------------------------------------
##  card and workspace class
## ------------------------------------

class datacardClass(object):


    #variables
    isFSR = True
    DEBUG = False
    ID_4mu = 1
    ID_4e = 2
    ID_2e2mu = 3

    def __init__(self):
	if self.DEBUG: level = 10
	else: level = 20
	self.log = Logger().getLogger(self.__class__.__name__,level)
	
    def setVariables(self, theMH,theOutputDir,theInputs,theTemplateDir,theMassError,theis2D,theUseMEKD):

        ## --------------- SETTINGS AND DECLARATIONS --------------- ##
        self.mH = theMH
        self.lumi = theInputs['lumi']
        self.sqrts = theInputs['sqrts']
        self.channel = theInputs['decayChannel']
        self.sigMorph = theInputs['useCMS_zz4l_sigMELA']
        self.bkgMorph = theInputs['useCMS_zz4l_bkgMELA']
        self.templateDir = theTemplateDir
        self.outputDir = theOutputDir        
        self.bUseCBnoConvolution = False
                
        self.bIncludingError = theMassError
        self.is2D = theis2D
        self.useMEKD = theUseMEKD
        
        self.all_chan = theInputs['all']
        self.ggH_chan = theInputs['ggH']
        self.qqH_chan = theInputs['qqH']
        self.WH_chan = theInputs['WH']
        self.ZH_chan = theInputs['ZH']
        self.ttH_chan = theInputs['ttH']
        self.qqZZ_chan = theInputs['qqZZ']
        self.ggZZ_chan = theInputs['ggZZ']
        self.zjets_chan = theInputs['zjets']
        self.ttbar_chan = theInputs['ttbar']
        self.zbb_chan = theInputs['zbb']
        self.isAltSig = theInputs['doHypTest']
        self.appendHypType = theInputs['altHypLabel']
        self.inputs = theInputs

        self.reflectZX = False
        #self.reflectZX = True #mod-roko
        
        if (self.channel == self.ID_4mu): self.appendName = '4mu'
        elif (self.channel == self.ID_4e): self.appendName = '4e'
        elif (self.channel == self.ID_2e2mu): self.appendName = '2e2mu'
        else:
            raise RuntimeError, "Input Error: Unknown channel! (4mu = 1, 4e = 2, 2e2mu = 3)"

            
    # return trueVar if testStatement else return falseVar
    def getVariable(self,trueVar,falseVar,testStatement):

        if (testStatement):
            return trueVar
        else:
            return falseVar
        
                                

    # cross section filter for 7 TeV efficiency
    def csFilter(self,hmass):

        a = 80.85
        b = 50.42
        
        f = 0.5 + 0.5*erf( (hmass - a)/b )
        
        return f

    # cs x br function 
    def makeXsBrFunction(self,signalProc,rrvMH):
            
        procName = "ggH"
        if(signalProc == 0): procName = "ggH" #dummy, when you sum up all the 5 chans
        if(signalProc == 1): procName = "ggH"
        if(signalProc == 2): procName = "qqH"
        if(signalProc == 3): procName = "WH"
        if(signalProc == 4): procName = "ZH"
        if(signalProc == 5): procName = "ttH"

        
        channelName = ""
        if (self.channel == self.ID_4mu): channelName = "4mu"
        elif (self.channel == self.ID_4e): channelName = "4e"
        elif (self.channel == self.ID_2e2mu): channelName = "2e2mu"
        else:
            raise RuntimeError, "Input Error: Unknown channel! (4mu = 1, 4e = 2, 2e2mu = 3)" 

     
        
        self.myCSWrhf = HiggsCSandWidth()
        
        histXsBr = ROOT.TH1F("hsmxsbr_{0}_{1}".format(procName,channelName),"", 8905, 109.55, 1000.05)
                
        for i in range(1,8906):
            
            mHVal = histXsBr.GetBinCenter(i)
            BR = 0.0 
            if (self.channel == self.ID_2e2mu):
                BR = self.myCSWrhf.HiggsBR(13,mHVal)
            else:
                BR = self.myCSWrhf.HiggsBR(12,mHVal)

            if (signalProc == 3 or signalProc == 4 or signalProc == 5):
                #overwrite BR if VH,ttH sample
                #these samples have inclusive Z decay
                BR = self.myCSWrhf.HiggsBR(11,mHVal)
                
            if (signalProc==0):
                totXs=0
                for ch in range(1,6):
                    totXs+=self.myCSWrhf.HiggsCS(ch, mHVal, self.sqrts)
                histXsBr.SetBinContent(i, totXs * BR)
            else:
                histXsBr.SetBinContent(i, self.myCSWrhf.HiggsCS(signalProc, mHVal, self.sqrts) * BR)

            #print '\nmakeXsBrFunction : procName=',procName,'   signalProc=',signalProc,'  mH (input)=',rrvMH.getVal(),
            #print '   CS=',self.myCSWrhf.HiggsCS(signalProc, mHVal, self.sqrts),'   BR=',BR
            
        rdhname = "rdhXsBr_{0}_{1}_{2}".format(procName,self.channel,self.sqrts)
        rdhXsBr = RooDataHist(rdhname,rdhname, ROOT.RooArgList(rrvMH), histXsBr)  
        
        return rdhXsBr
    
    
    # main datacard and workspace function
    def makeMassShapesYields(self,theMH,theOutputDir,theInputs,theTemplateDir,theMassError,theis2D,theUseMEKD):
        
        self.setVariables(theMH,theOutputDir,theInputs,theTemplateDir, theMassError,theis2D,theUseMEKD)

        self.myCSW = HiggsCSandWidth()
        self.widthHVal =  self.myCSW.HiggsWidth(0,self.mH)
        self.isHighMass = False
        if(self.widthHVal < 0.12):
            self.bUseCBnoConvolution = True
        if self.mH >= 390:
            if self.inputs['useHighMassReweightedShapes']:
                self.isHighMass = True
            else:
                print ">>>>>> useHighMassReweightedShapes set to FALSE, using non-reweighted shapes!"

        ## ---------------- SET PLOTTING STYLE ---------------- ## 
        ROOT.setTDRStyle(True)
        ROOT.gStyle.SetPalette(1)
        ROOT.gStyle.SetPadLeftMargin(0.16)        
        
        ## ------------------------- SYSTEMATICS CLASSES ----------------------------- ##
        self.systematics = systematicsClass( self.mH, False, self.isFSR, self.inputs)
        self.systematics_forXSxBR = systematicsClass( self.mH, True, self.isFSR, self.inputs)

        ## -------------------------- SIGNAL SHAPE ----------------------------------- ##
        bins = 1000
        if(self.bUseCBnoConvolution): bins = 200
        if(self.bIncludingError): bins = 200

        self.CMS_zz4l_mass = ROOT.RooRealVar("CMS_zz4l_mass","CMS_zz4l_mass",self.inputs['low_M'],self.inputs['high_M'])
        self.CMS_zz4l_mass.setBins(bins,"fft") 

        #self.LUMI = ROOT.RooRealVar("LUMI_{0:.0f}".format(self.sqrts),"LUMI_{0:.0f}".format(self.sqrts),self.lumi)
        #self.LUMI.setConstant(True)
    
        self.MH = ROOT.RooRealVar("MH","MH",self.mH)
        self.MH.setConstant(True)
    
	# n2, alpha2 are right side parameters of DoubleCB
	# n, alpha are left side parameters of DoubleCB

        n_CB_d = 0.0
        alpha_CB_d = 0.0
        n2_CB_d = 0.0
        alpha2_CB_d = 0.0
        mean_CB_d = 0.0
        sigma_CB_d = 0.0
        mean_BW_d = self.mH
        gamma_BW_d = 0.0
        
#        if(self.all_chan):
#            self.rdhXsBrFuncV_1 = self.makeXsBrFunction(0,self.MH)
#        else:
        self.rdhXsBrFuncV_1 = self.makeXsBrFunction(1,self.MH)
        rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}".format("ggH",self.channel,self.sqrts)
        self.rhfXsBrFuncV_1 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), self.rdhXsBrFuncV_1, 1)
        
        self.rdhXsBrFuncV_2 = self.makeXsBrFunction(2,self.MH)
        rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}".format("VBF",self.channel,self.sqrts)
        self.rhfXsBrFuncV_2 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), self.rdhXsBrFuncV_2, 1)
        
        self.rdhXsBrFuncV_3 = self.makeXsBrFunction(3,self.MH)
        rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}".format("WH",self.channel,self.sqrts)
        self.rhfXsBrFuncV_3 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), self.rdhXsBrFuncV_3, 1)
        
        self.rdhXsBrFuncV_4 = self.makeXsBrFunction(4,self.MH)
        rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}".format("ZH",self.channel,self.sqrts)
        self.rhfXsBrFuncV_4 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), self.rdhXsBrFuncV_4, 1)
        
        self.rdhXsBrFuncV_5 = self.makeXsBrFunction(5,self.MH)
        rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}".format("ttH",self.channel,self.sqrts)
        self.rhfXsBrFuncV_5 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), self.rdhXsBrFuncV_5, 1)

    
        ## -------- Variable Definitions -------- ##
        name = "CMS_zz4l_mean_e_sig"
        self.CMS_zz4l_mean_e_sig = ROOT.RooRealVar(name,"CMS_zz4l_mean_e_sig",0.0,-10.0,10.0)
        name = "CMS_zz4l_sigma_e_sig"
        self.CMS_zz4l_sigma_e_sig = ROOT.RooRealVar(name,"CMS_zz4l_sigma_sig",3.0,0.0,30.0)
        
        name = "CMS_zz4l_mean_m_sig"
        self.CMS_zz4l_mean_m_sig = ROOT.RooRealVar(name,"CMS_zz4l_mean_sig",0.0,-10.0,10.0)
        name = "CMS_zz4l_sigma_m_sig"
        self.CMS_zz4l_sigma_m_sig = ROOT.RooRealVar(name,"CMS_zz4l_sigma_sig",3.0,0.0,30.0)
        
        name = "CMS_zz4l_alpha2_{0}_{1:.0f}".format(self.channel,self.sqrts)
        self.CMS_zz4l_alpha2 = ROOT.RooRealVar(name,"CMS_zz4l_alpha2",1.,-10.,10.)
        name = "CMS_zz4l_n2_sig_{0}_{1:.0f}".format(self.channel,self.sqrts)
        self.CMS_zz4l_n2 = ROOT.RooRealVar(name,"CMS_zz4l_n2",2.,-10.,10.)
        name = "CMS_zz4l_alpha_{0}_{1:.0f}".format(self.channel,self.sqrts)
        self.CMS_zz4l_alpha = ROOT.RooRealVar(name,"CMS_zz4l_alpha",1.,-10.,10.)
        name = "CMS_zz4l_n_sig_{0}_{1:.0f}".format(self.channel,self.sqrts)
        self.CMS_zz4l_n = ROOT.RooRealVar(name,"CMS_zz4l_n",2.,-10.,10.)
        name = "CMS_zz4l_mean_BW_{0}_{1:.0f}".format(self.channel,self.sqrts)
        self.CMS_zz4l_mean_BW = ROOT.RooRealVar(name,"CMS_zz4l_mean_BW",self.mH,self.inputs['low_M'],self.inputs['high_M'])
        name = "interf_ggH"
        #name = "CMS_zz4l_gamma_sig_{0}_{1:.0f}".format(self.channel,self.sqrts)
        self.CMS_zz4l_gamma = ROOT.RooRealVar(name,"CMS_zz4l_gamma",10.,0.001,1000.)
        name = "CMS_zz4l_widthScale_{0}_{1:.0f}".format(self.channel,self.sqrts)
        self.CMS_zz4l_widthScale = ROOT.RooRealVar(name,"CMS_zz4l_widthScale",1.0)
        self.one = ROOT.RooRealVar("one","one",1.0)
        self.one.setConstant(True)
    
        self.CMS_zz4l_mean_BW.setVal( mean_BW_d )
        self.CMS_zz4l_gamma.setVal(0)
        self.CMS_zz4l_mean_e_sig.setVal(0)
        self.CMS_zz4l_sigma_e_sig.setVal(0)
        self.CMS_zz4l_mean_m_sig.setVal(0)
        self.CMS_zz4l_sigma_m_sig.setVal(0)
        self.CMS_zz4l_alpha.setVal(0)
        self.CMS_zz4l_n.setVal(0)
        self.CMS_zz4l_alpha2.setVal(0)
        self.CMS_zz4l_n2.setVal(0)
    
        self.CMS_zz4l_widthScale.setConstant(True)
        #self.CMS_zz4l_alpha.setConstant(True)  # also read from input file
        self.CMS_zz4l_mean_BW.setConstant(True)
        #CMS_zz4l_gamma_BW.setConstant(True)

        print ">>>>>> mean_BW ", self.CMS_zz4l_mean_BW.getVal()
        print ">>>>>> gamma_BW ", self.CMS_zz4l_gamma.getVal()
        print ">>>>>> mean_e_sig ", self.CMS_zz4l_mean_e_sig.getVal()
        print ">>>>>> sigma_e ", self.CMS_zz4l_sigma_e_sig.getVal()
        print ">>>>>> mean_m_sig ", self.CMS_zz4l_mean_m_sig.getVal()
        print ">>>>>> sigma_m ", self.CMS_zz4l_sigma_m_sig.getVal()
        print ">>>>>> alpha ", self.CMS_zz4l_alpha.getVal()
        print ">>>>>> n ", self.CMS_zz4l_n.getVal()
        print ">>>>>> alpha2 ", self.CMS_zz4l_alpha2.getVal()
        print ">>>>>> n2 ", self.CMS_zz4l_n2.getVal()

                                                                
        ## -------------------- RooFormulaVar's -------------------- ##
        self.rfv_n_CB = ROOT.RooFormulaVar()
        self.rfv_alpha_CB = ROOT.RooFormulaVar()
        self.rfv_n2_CB = ROOT.RooFormulaVar()
        self.rfv_alpha2_CB = ROOT.RooFormulaVar()
        self.rfv_mean_CB = ROOT.RooFormulaVar()
        self.rfv_sigma_CB = ROOT.RooFormulaVar()
        
        name = "CMS_zz4l_n_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
        if self.isHighMass :
            self.rfv_n_CB = ROOT.RooFormulaVar(name,"("+self.inputs['n_CB_shape_HM']+")"+"*(1+@1)",ROOT.RooArgList(self.MH,self.CMS_zz4l_n))
        else :
            self.rfv_n_CB = ROOT.RooFormulaVar(name,"("+self.inputs['n_CB_shape']+")"+"*(1+@1)",ROOT.RooArgList(self.MH,self.CMS_zz4l_n))

        name = "CMS_zz4l_alpha_{0:.0f}_centralValue".format(self.channel)
        if self.isHighMass :
            self.rfv_alpha_CB = ROOT.RooFormulaVar(name,self.inputs['alpha_CB_shape_HM'], ROOT.RooArgList(self.MH))
        else :
            self.rfv_alpha_CB = ROOT.RooFormulaVar(name,self.inputs['alpha_CB_shape'], ROOT.RooArgList(self.MH))

        name = "CMS_zz4l_n2_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
        #if self.isHighMass : rfv_n2_CB = ROOT.RooFormulaVar(name,"("+self.inputs['n2_CB_shape_HM']+")"+"*(1+@1)",ROOT.RooArgList(self.MH,CMS_zz4l_n2))
        #else : rfv_n2_CB = ROOT.RooFormulaVar(name,"("+self.inputs['n2_CB_shape']+")"+"*(1+@1)",ROOT.RooArgList(self.MH,CMS_zz4l_n2))
        if self.isHighMass :
            self.rfv_n2_CB = ROOT.RooFormulaVar(name,"("+self.inputs['n2_CB_shape_HM']+")",ROOT.RooArgList(self.MH))
        else :
            self.rfv_n2_CB = ROOT.RooFormulaVar(name,"("+self.inputs['n2_CB_shape']+")",ROOT.RooArgList(self.MH))

        name = "CMS_zz4l_alpha2_{0:.0f}_centralValue".format(self.channel)
        if self.isHighMass :
            self.rfv_alpha2_CB = ROOT.RooFormulaVar(name,self.inputs['alpha2_CB_shape_HM'], ROOT.RooArgList(self.MH))
        else :
            self.rfv_alpha2_CB = ROOT.RooFormulaVar(name,self.inputs['alpha2_CB_shape'], ROOT.RooArgList(self.MH))

        name = "CMS_zz4l_mean_sig_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
        if (self.channel == self.ID_4mu) :

            if self.isHighMass :
                self.rfv_mean_CB = ROOT.RooFormulaVar(name,"("+self.inputs['mean_CB_shape_HM']+")"+"+@0*@1", ROOT.RooArgList(self.MH, self.CMS_zz4l_mean_m_sig))
            else :
                self.rfv_mean_CB = ROOT.RooFormulaVar(name,"("+self.inputs['mean_CB_shape']+")"+"+@0*@1", ROOT.RooArgList(self.MH, self.CMS_zz4l_mean_m_sig))

        elif (self.channel == self.ID_4e) :

            if self.isHighMass :
                self.rfv_mean_CB = ROOT.RooFormulaVar(name,"("+self.inputs['mean_CB_shape_HM']+")"+"+@0*@1", ROOT.RooArgList(self.MH, self.CMS_zz4l_mean_e_sig))
            else :
                self.rfv_mean_CB = ROOT.RooFormulaVar(name,"("+self.inputs['mean_CB_shape']+")"+"+@0*@1", ROOT.RooArgList(self.MH, self.CMS_zz4l_mean_e_sig))

        elif (self.channel == self.ID_2e2mu) :

            if self.isHighMass :
                self.rfv_mean_CB = ROOT.RooFormulaVar(name,"("+self.inputs['mean_CB_shape_HM']+")"+"+ (@0*@1 + @0*@2)/2", ROOT.RooArgList(self.MH, self.CMS_zz4l_mean_m_sig,self.CMS_zz4l_mean_e_sig))
            else :
                self.rfv_mean_CB = ROOT.RooFormulaVar(name,"("+self.inputs['mean_CB_shape']+")"+"+ (@0*@1 + @0*@2)/2", ROOT.RooArgList(self.MH, self.CMS_zz4l_mean_m_sig,self.CMS_zz4l_mean_e_sig))
        

        name = "CMS_zz4l_sigma_sig_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
        if (self.channel == self.ID_4mu) :
            if self.isHighMass :
                self.rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+self.inputs['sigma_CB_shape_HM']+")"+"*(1+@1)", ROOT.RooArgList(self.MH, self.CMS_zz4l_sigma_m_sig))
            else :
                self.rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+self.inputs['sigma_CB_shape']+")"+"*(1+@1)", ROOT.RooArgList(self.MH, self.CMS_zz4l_sigma_m_sig))
        elif (self.channel == self.ID_4e) :

            if self.isHighMass :
                self.rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+self.inputs['sigma_CB_shape_HM']+")"+"*(1+@1)", ROOT.RooArgList(self.MH, self.CMS_zz4l_sigma_e_sig))
            else :
                self.rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+self.inputs['sigma_CB_shape']+")"+"*(1+@1)", ROOT.RooArgList(self.MH, self.CMS_zz4l_sigma_e_sig))
        elif (self.channel == self.ID_2e2mu) :

            if self.isHighMass :
                self.rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+self.inputs['sigma_CB_shape_HM']+")"+"*TMath::Sqrt((1+@1)*(1+@2))", ROOT.RooArgList(self.MH, self.CMS_zz4l_sigma_m_sig,self.CMS_zz4l_sigma_e_sig))
            else :
                self.rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+self.inputs['sigma_CB_shape']+")"+"*TMath::Sqrt((1+@1)*(1+@2))", ROOT.RooArgList(self.MH, self.CMS_zz4l_sigma_m_sig,self.CMS_zz4l_sigma_e_sig))


        name = "CMS_zz4l_gamma_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
        self.rfv_gamma_BW = ROOT.RooFormulaVar(name,"("+self.inputs['gamma_BW_shape_HM']+")"+"*(1+@1*0.05)",ROOT.RooArgList(self.MH,self.CMS_zz4l_gamma))

        print ">>>>>> n_CB ", self.rfv_n_CB.getVal()
        print ">>>>>> alpha_CB ", self.rfv_alpha_CB.getVal()
        print ">>>>>> n2_CB ", self.rfv_n2_CB.getVal()
        print ">>>>>> alpha2_CB ", self.rfv_alpha2_CB.getVal()
        print ">>>>>> mean_CB ", self.rfv_mean_CB.getVal()
        print ">>>>>> sigma_CB ", self.rfv_sigma_CB.getVal()
        print ">>>>>> gamma_BW ", self.rfv_gamma_BW.getVal()    

        
        self.CMS_zz4l_mean_sig_NoConv = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_NoConv_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts),"@0+@1", ROOT.RooArgList(self.rfv_mean_CB, self.MH))


        ## ----------------- NEED MASS ERRORS FOR CB ------------------ ##
        self.defineMassErr()

        ## --------------------- SHAPE FUNCTIONS ---------------------- ##
        ##ggH
        self.signalCB_ggH = ROOT.RooDoubleCB("signalCB_ggH","signalCB_ggH",self.CMS_zz4l_mass, self.getVariable(self.CMS_zz4l_mean_sig_NoConv,self.rfv_mean_CB, self.bUseCBnoConvolution), self.getVariable(self.CMS_zz4l_massErr,self.rfv_sigma_CB,self.bIncludingError),self.rfv_alpha_CB,self.rfv_n_CB, self.rfv_alpha2_CB, self.rfv_n2_CB)
        #Low mass pdf
        self.signalBW_ggH = ROOT.RooRelBWUFParam("signalBW_ggH", "signalBW_ggH",self.CMS_zz4l_mass,self.CMS_zz4l_mean_BW,self.CMS_zz4l_widthScale)
        self.sig_ggH =  ROOT.RooFFTConvPdf("sig_ggH","BW (X) CB",self.CMS_zz4l_mass,self.signalBW_ggH,self.signalCB_ggH,2)
        #High mass pdf
        self.signalBW_ggH_HM = ROOT.RooRelBWHighMass("signalBW_ggH", "signalBW_ggH",self.CMS_zz4l_mass,self.CMS_zz4l_mean_BW,self.rfv_gamma_BW)
        self.sig_ggH_HM =  ROOT.RooFFTConvPdf("sig_ggH","BW (X) CB",self.CMS_zz4l_mass,self.signalBW_ggH_HM,self.signalCB_ggH, 2)
  
        ##VBF
        self.signalCB_VBF = ROOT.RooDoubleCB("signalCB_VBF","signalCB_VBF",self.CMS_zz4l_mass,self.getVariable(self.CMS_zz4l_mean_sig_NoConv,self.rfv_mean_CB,self.bUseCBnoConvolution),self.getVariable(self.CMS_zz4l_massErr,self.rfv_sigma_CB, self.bIncludingError),self.rfv_alpha_CB,self.rfv_n_CB, self.rfv_alpha2_CB,self.rfv_n2_CB)
        #Low mass pdf
        self.signalBW_VBF = ROOT.RooRelBWUFParam("signalBW_VBF", "signalBW_VBF",self.CMS_zz4l_mass,self.CMS_zz4l_mean_BW,self.CMS_zz4l_widthScale)
        self.sig_VBF = ROOT.RooFFTConvPdf("sig_VBF","BW (X) CB",self.CMS_zz4l_mass,self.signalBW_VBF,self.signalCB_VBF, 2)
        #High mass pdf
        self.signalBW_VBF_HM = ROOT.RooRelBWHighMass("signalBW_VBF", "signalBW_VBF",self.CMS_zz4l_mass,self.CMS_zz4l_mean_BW,self.rfv_gamma_BW)
        self.sig_VBF_HM = ROOT.RooFFTConvPdf("sig_VBF","BW (X) CB",self.CMS_zz4l_mass,self.signalBW_VBF_HM,self.signalCB_VBF, 2)
                       
        ##WH
        self.signalCB_WH = ROOT.RooDoubleCB("signalCB_WH","signalCB_WH",self.CMS_zz4l_mass,self.getVariable(self.CMS_zz4l_mean_sig_NoConv,self.rfv_mean_CB,self.bUseCBnoConvolution),self.getVariable(self.CMS_zz4l_massErr,self.rfv_sigma_CB, self.bIncludingError),self.rfv_alpha_CB,self.rfv_n_CB,self.rfv_alpha2_CB,self.rfv_n2_CB)
        #Low mass pdf
        self.signalBW_WH = ROOT.RooRelBWUFParam("signalBW_WH", "signalBW_WH",self.CMS_zz4l_mass,self.CMS_zz4l_mean_BW,self.CMS_zz4l_widthScale)
        self.sig_WH = ROOT.RooFFTConvPdf("sig_WH","BW (X) CB",self.CMS_zz4l_mass,self.signalBW_WH,self.signalCB_WH, 2)
        #High mass pdf
        self.signalBW_WH_HM = ROOT.RooRelBWHighMass("signalBW_WH", "signalBW_WH",self.CMS_zz4l_mass,self.CMS_zz4l_mean_BW,self.rfv_gamma_BW)
        self.sig_WH_HM = ROOT.RooFFTConvPdf("sig_WH","BW (X) CB",self.CMS_zz4l_mass,self.signalBW_WH_HM,self.signalCB_WH, 2)

        ##ZH
        self.signalCB_ZH = ROOT.RooDoubleCB("signalCB_ZH","signalCB_ZH",self.CMS_zz4l_mass,self.getVariable(self.CMS_zz4l_mean_sig_NoConv,self.rfv_mean_CB,self.bUseCBnoConvolution),self.getVariable(self.CMS_zz4l_massErr,self.rfv_sigma_CB, self.bIncludingError),self.rfv_alpha_CB,self.rfv_n_CB,self.rfv_alpha2_CB,self.rfv_n2_CB)
        #Low mass pdf
        self.signalBW_ZH = ROOT.RooRelBWUFParam("signalBW_ZH", "signalBW_ZH",self.CMS_zz4l_mass,self.CMS_zz4l_mean_BW,self.CMS_zz4l_widthScale)
        self.sig_ZH = ROOT.RooFFTConvPdf("sig_ZH","BW (X) CB",self.CMS_zz4l_mass,self.signalBW_ZH,self.signalCB_ZH, 2)
        #High mass pdf
        self.signalBW_ZH_HM = ROOT.RooRelBWHighMass("signalBW_ZH", "signalBW_ZH",self.CMS_zz4l_mass,self.CMS_zz4l_mean_BW,self.rfv_gamma_BW)
        self.sig_ZH_HM = ROOT.RooFFTConvPdf("sig_ZH","BW (X) CB",self.CMS_zz4l_mass,self.signalBW_ZH_HM,self.signalCB_ZH, 2)

        ##ttH
        self.signalCB_ttH = ROOT.RooDoubleCB("signalCB_ttH","signalCB_ttH",self.CMS_zz4l_mass,self.getVariable(self.CMS_zz4l_mean_sig_NoConv,self.rfv_mean_CB,self.bUseCBnoConvolution),self.getVariable(self.CMS_zz4l_massErr,self.rfv_sigma_CB,self.bIncludingError),self.rfv_alpha_CB,self.rfv_n_CB,self.rfv_alpha2_CB,self.rfv_n2_CB)
        #Low mass pdf
        self.signalBW_ttH = ROOT.RooRelBWUFParam("signalBW_ttH", "signalBW_ttH",self.CMS_zz4l_mass,self.CMS_zz4l_mean_BW,self.CMS_zz4l_widthScale)
        self.sig_ttH = ROOT.RooFFTConvPdf("sig_ttH","BW (X) CB",self.CMS_zz4l_mass,self.signalBW_ttH,self.signalCB_ttH, 2) 
        #High mass pdf
        self.signalBW_ttH_HM = ROOT.RooRelBWHighMass("signalBW_ttH", "signalBW_ttH",self.CMS_zz4l_mass,self.CMS_zz4l_mean_BW,self.rfv_gamma_BW)
        self.sig_ttH_HM = ROOT.RooFFTConvPdf("sig_ttH","BW (X) CB",self.CMS_zz4l_mass,self.signalBW_ttH_HM,self.signalCB_ttH, 2)
        
        
        ## Buffer fraction for cyclical behavior
        self.sig_ggH.setBufferFraction(0.2)
        self.sig_VBF.setBufferFraction(0.2)
        self.sig_WH.setBufferFraction(0.2)
        self.sig_ZH.setBufferFraction(0.2)
        self.sig_ttH.setBufferFraction(0.2)
        
        self.sig_ggH_HM.setBufferFraction(0.2)
        self.sig_VBF_HM.setBufferFraction(0.2)
        self.sig_WH_HM.setBufferFraction(0.2)
        self.sig_ZH_HM.setBufferFraction(0.2)
        self.sig_ttH_HM.setBufferFraction(0.2)


        ## -------------------------- BACKGROUND SHAPES ---------------------------------- ##
    
        ## qqZZ contribution
        name = "CMS_qqzzbkg_a0_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        self.CMS_qqzzbkg_a0 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a0",115.3,0.,200.)
        name = "CMS_qqzzbkg_a1_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        self.CMS_qqzzbkg_a1 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a1",21.96,0.,200.)
        name = "CMS_qqzzbkg_a2_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        self.CMS_qqzzbkg_a2 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a2",122.8,0.,200.)
        name = "CMS_qqzzbkg_a3_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        self.CMS_qqzzbkg_a3 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a3",0.03479,0.,1.)
        name = "CMS_qqzzbkg_a4_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        self.CMS_qqzzbkg_a4 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a4",185.5,0.,200.)
        name = "CMS_qqzzbkg_a5_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        self.CMS_qqzzbkg_a5 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a5",12.67,0.,200.)
        name = "CMS_qqzzbkg_a6_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        self.CMS_qqzzbkg_a6 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a6",34.81,0.,100.)
        name = "CMS_qqzzbkg_a7_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        self.CMS_qqzzbkg_a7 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a7",0.1393,0.,1.)
        name = "CMS_qqzzbkg_a8_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        self.CMS_qqzzbkg_a8 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a8",66.,0.,200.)
        name = "CMS_qqzzbkg_a9_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        self.CMS_qqzzbkg_a9 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a9",0.07191,0.,1.)
        name = "CMS_qqzzbkg_a10_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts )
        self.CMS_qqzzbkg_a10 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a10",94.11,0.,200.)
        name = "CMS_qqzzbkg_a11_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts )
        self.CMS_qqzzbkg_a11 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a11",-5.111,-100.,100.)
        name = "CMS_qqzzbkg_a12_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts )
        self.CMS_qqzzbkg_a12 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a12",4834,0.,10000.)
        name = "CMS_qqzzbkg_a13_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts )
        self.CMS_qqzzbkg_a13 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a13",0.2543,0.,1.)
        

        if (self.DEBUG) :
            print "qqZZshape_a0 = ",self.inputs['qqZZshape_a0']
            print "qqZZshape_a1 = ",self.inputs['qqZZshape_a1']
            print "qqZZshape_a2 = ",self.inputs['qqZZshape_a2']
            print "qqZZshape_a3 = ",self.inputs['qqZZshape_a3']
            print "qqZZshape_a4 = ",self.inputs['qqZZshape_a4']
            print "qqZZshape_a5 = ",self.inputs['qqZZshape_a5']
            print "qqZZshape_a6 = ",self.inputs['qqZZshape_a6']
            print "qqZZshape_a7 = ",self.inputs['qqZZshape_a7']
            print "qqZZshape_a8 = ",self.inputs['qqZZshape_a8']
            print "qqZZshape_a9 = ",self.inputs['qqZZshape_a9']
            print "qqZZshape_a10 = ",self.inputs['qqZZshape_a10']
            print "qqZZshape_a11 = ",self.inputs['qqZZshape_a11']
            print "qqZZshape_a12 = ",self.inputs['qqZZshape_a12']
            print "qqZZshape_a13 = ",self.inputs['qqZZshape_a13']

        
        self.CMS_qqzzbkg_a0.setVal(self.inputs['qqZZshape_a0'])
        self.CMS_qqzzbkg_a1.setVal(self.inputs['qqZZshape_a1'])
        self.CMS_qqzzbkg_a2.setVal(self.inputs['qqZZshape_a2'])
        self.CMS_qqzzbkg_a3.setVal(self.inputs['qqZZshape_a3'])
        self.CMS_qqzzbkg_a4.setVal(self.inputs['qqZZshape_a4'])
        self.CMS_qqzzbkg_a5.setVal(self.inputs['qqZZshape_a5'])
        self.CMS_qqzzbkg_a6.setVal(self.inputs['qqZZshape_a6'])
        self.CMS_qqzzbkg_a7.setVal(self.inputs['qqZZshape_a7'])
        self.CMS_qqzzbkg_a8.setVal(self.inputs['qqZZshape_a8'])
        self.CMS_qqzzbkg_a9.setVal(self.inputs['qqZZshape_a9'])
        self.CMS_qqzzbkg_a10.setVal(self.inputs['qqZZshape_a10'])
        self.CMS_qqzzbkg_a11.setVal(self.inputs['qqZZshape_a11'])
        self.CMS_qqzzbkg_a12.setVal(self.inputs['qqZZshape_a12'])
        self.CMS_qqzzbkg_a13.setVal(self.inputs['qqZZshape_a13'])
        
        self.CMS_qqzzbkg_a0.setConstant(True)
        self.CMS_qqzzbkg_a1.setConstant(True)
        self.CMS_qqzzbkg_a2.setConstant(True)
        self.CMS_qqzzbkg_a3.setConstant(True)
        self.CMS_qqzzbkg_a4.setConstant(True)
        self.CMS_qqzzbkg_a5.setConstant(True)
        self.CMS_qqzzbkg_a6.setConstant(True)
        self.CMS_qqzzbkg_a7.setConstant(True)
        self.CMS_qqzzbkg_a8.setConstant(True)
        self.CMS_qqzzbkg_a9.setConstant(True)
        self.CMS_qqzzbkg_a10.setConstant(True)
        self.CMS_qqzzbkg_a11.setConstant(True)
        self.CMS_qqzzbkg_a12.setConstant(True)
        self.CMS_qqzzbkg_a13.setConstant(True)
        
        self.bkg_qqzz = ROOT.RooqqZZPdf_v2("bkg_qqzzTmp","bkg_qqzzTmp",self.CMS_zz4l_mass,self.CMS_qqzzbkg_a0,self.CMS_qqzzbkg_a1,self.CMS_qqzzbkg_a2,self.CMS_qqzzbkg_a3,self.CMS_qqzzbkg_a4,self.CMS_qqzzbkg_a5,self.CMS_qqzzbkg_a6,self.CMS_qqzzbkg_a7,self.CMS_qqzzbkg_a8,self.CMS_qqzzbkg_a9,self.CMS_qqzzbkg_a10,self.CMS_qqzzbkg_a11,self.CMS_qqzzbkg_a12,self.CMS_qqzzbkg_a13)
        
        ## ggZZ contribution
        name = "CMS_ggzzbkg_a0_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        self.CMS_ggzzbkg_a0 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a0",115.3,0.,200.)
        name = "CMS_ggzzbkg_a1_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        self.CMS_ggzzbkg_a1 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a1",21.96,0.,200.)
        name = "CMS_ggzzbkg_a2_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        self.CMS_ggzzbkg_a2 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a2",122.8,0.,200.)
        name = "CMS_ggzzbkg_a3_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        self.CMS_ggzzbkg_a3 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a3",0.03479,0.,1.)
        name = "CMS_ggzzbkg_a4_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts )
        self.CMS_ggzzbkg_a4 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a4",185.5,0.,200.)
        name = "CMS_ggzzbkg_a5_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        self.CMS_ggzzbkg_a5 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a5",12.67,0.,200.)
        name = "CMS_ggzzbkg_a6_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        self.CMS_ggzzbkg_a6 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a6",34.81,0.,100.)
        name = "CMS_ggzzbkg_a7_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        self.CMS_ggzzbkg_a7 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a7",0.1393,0.,1.)
        name = "CMS_ggzzbkg_a8_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        self.CMS_ggzzbkg_a8 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a8",66.,0.,200.)
        name = "CMS_ggzzbkg_a9_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts )
        self.CMS_ggzzbkg_a9 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a9",0.07191,0.,1.)
        
        self.CMS_ggzzbkg_a0.setVal(self.inputs['ggZZshape_a0'])
        self.CMS_ggzzbkg_a1.setVal(self.inputs['ggZZshape_a1'])
        self.CMS_ggzzbkg_a2.setVal(self.inputs['ggZZshape_a2'])
        self.CMS_ggzzbkg_a3.setVal(self.inputs['ggZZshape_a3'])
        self.CMS_ggzzbkg_a4.setVal(self.inputs['ggZZshape_a4'])
        self.CMS_ggzzbkg_a5.setVal(self.inputs['ggZZshape_a5'])
        self.CMS_ggzzbkg_a6.setVal(self.inputs['ggZZshape_a6'])
        self.CMS_ggzzbkg_a7.setVal(self.inputs['ggZZshape_a7'])
        self.CMS_ggzzbkg_a8.setVal(self.inputs['ggZZshape_a8'])
        self.CMS_ggzzbkg_a9.setVal(self.inputs['ggZZshape_a9'])
        
        self.CMS_ggzzbkg_a0.setConstant(True)
        self.CMS_ggzzbkg_a1.setConstant(True)
        self.CMS_ggzzbkg_a2.setConstant(True)
        self.CMS_ggzzbkg_a3.setConstant(True)
        self.CMS_ggzzbkg_a4.setConstant(True)
        self.CMS_ggzzbkg_a5.setConstant(True)
        self.CMS_ggzzbkg_a6.setConstant(True)
        self.CMS_ggzzbkg_a7.setConstant(True)
        self.CMS_ggzzbkg_a8.setConstant(True)
        self.CMS_ggzzbkg_a9.setConstant(True)

        if (self.DEBUG) :
            print "ggZZshape_a0 = ",self.inputs['ggZZshape_a0']
            print "ggZZshape_a1 = ",self.inputs['ggZZshape_a1']
            print "ggZZshape_a2 = ",self.inputs['ggZZshape_a2']
            print "ggZZshape_a3 = ",self.inputs['ggZZshape_a3']
            print "ggZZshape_a4 = ",self.inputs['ggZZshape_a4']
            print "ggZZshape_a5 = ",self.inputs['ggZZshape_a5']
            print "ggZZshape_a6 = ",self.inputs['ggZZshape_a6']
            print "ggZZshape_a7 = ",self.inputs['ggZZshape_a7']
            print "ggZZshape_a8 = ",self.inputs['ggZZshape_a8']
            print "ggZZshape_a9 = ",self.inputs['ggZZshape_a9']
                   
        
        self.bkg_ggzz = ROOT.RooggZZPdf_v2("bkg_ggzzTmp","bkg_ggzzTmp",self.CMS_zz4l_mass,self.CMS_ggzzbkg_a0,self.CMS_ggzzbkg_a1,self.CMS_ggzzbkg_a2,self.CMS_ggzzbkg_a3,self.CMS_ggzzbkg_a4,self.CMS_ggzzbkg_a5,self.CMS_ggzzbkg_a6,self.CMS_ggzzbkg_a7,self.CMS_ggzzbkg_a8,self.CMS_ggzzbkg_a9)

    ######################
        ## Reducible backgrounds
        self.val_meanL_3P1F = float(theInputs['zjetsShape_mean_3P1F'])
        self.val_sigmaL_3P1F = float(theInputs['zjetsShape_sigma_3P1F'])
        self.val_normL_3P1F = float(theInputs['zjetsShape_norm_3P1F'])
        
        self.val_meanL_2P2F = float(theInputs['zjetsShape_mean_2P2F'])
        self.val_sigmaL_2P2F = float(theInputs['zjetsShape_sigma_2P2F'])
        self.val_normL_2P2F = float(theInputs['zjetsShape_norm_2P2F'])
        self.val_pol0_2P2F = float(theInputs['zjetsShape_pol0_2P2F'])
        self.val_pol1_2P2F = float(theInputs['zjetsShape_pol1_2P2F'])
        
        self.val_meanL_2P2F_2 = float(theInputs['zjetsShape_mean_2P2F_2e2mu'])
        self.val_sigmaL_2P2F_2 = float(theInputs['zjetsShape_sigma_2P2F_2e2mu'])
        self.val_normL_2P2F_2 = float(theInputs['zjetsShape_norm_2P2F_2e2mu'])

        if (self.channel == self.ID_4mu):
            name = "mlZjet_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.mlZjet = ROOT.RooRealVar(name,"mean landau Zjet",self.val_meanL_2P2F)
            name = "slZjet_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.slZjet = ROOT.RooRealVar(name,"sigma landau Zjet",self.val_sigmaL_2P2F)
            print "mean 4mu: ",self.mlZjet.getVal()
            print "sigma 4mu: ",self.slZjet.getVal()
            self.bkg_zjets = ROOT.RooLandau("bkg_zjetsTmp","bkg_zjetsTmp",self.CMS_zz4l_mass,self.mlZjet,self.slZjet)
        elif (self.channel == self.ID_4e):
            name = "mlZjet_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.mlZjet_2p2f = ROOT.RooRealVar(name,"mean landau Zjet 2p2f",self.val_meanL_2P2F)
            name = "slZjet_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.slZjet_2p2f = ROOT.RooRealVar(name,"sigma landau Zjet 2p2f",self.val_sigmaL_2P2F)
            name = "nlZjet_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.nlZjet_2p2f = ROOT.RooRealVar(name,"norm landau Zjet 2p2f",self.val_normL_2P2F)
            name = "p0Zjet_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.p0Zjet_2p2f = ROOT.RooRealVar(name,"p0 Zjet 2p2f",self.val_pol0_2P2F)
            name = "p1Zjet_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.p1Zjet_2p2f = ROOT.RooRealVar(name,"p1 Zjet 2p2f",self.val_pol1_2P2F)
            print "mean 2p2f 4e: ",self.mlZjet_2p2f.getVal()
            print "sigma 2p2f 4e: ",self.slZjet_2p2f.getVal()
            print "norm 2p2f 4e: ",self.nlZjet_2p2f.getVal()
            print "pol0 2p2f 4e: ",self.p0Zjet_2p2f.getVal()
            print "pol1 2p2f 4e: ",self.p1Zjet_2p2f.getVal()
            self.bkg_zjets_2p2f = ROOT.RooGenericPdf("bkg_zjetsTmp_2p2f","bkg_zjetsTmp_2p2f","(TMath::Landau(@0,@1,@2))*@3*(1.+ TMath::Exp(@4+@5*@0))",RooArgList(self.CMS_zz4l_mass,self.mlZjet_2p2f,self.slZjet_2p2f,self.nlZjet_2p2f,self.p0Zjet_2p2f,self.p1Zjet_2p2f))
            
            name = "mlZjet_3p1f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.mlZjet_3p1f = ROOT.RooRealVar(name,"mean landau Zjet 3p1f",self.val_meanL_3P1F)
            name = "slZjet_3p1f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.slZjet_3p1f = ROOT.RooRealVar(name,"sigma landau Zjet 3p1f",self.val_sigmaL_3P1F)
            name = "nlZjet_3p1f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.nlZjet_3p1f = ROOT.RooRealVar(name,"norm landau Zjet 3p1f",self.val_normL_3P1F)
            print "mean 3p1f 4e: ",self.mlZjet_3p1f.getVal()
            print "sigma 3p1f 4e: ",self.slZjet_3p1f.getVal()
            print "norm 3p1f 4e: ",self.nlZjet_3p1f.getVal()
            self.bkg_zjets_3p1f = ROOT.RooLandau("bkg_zjetsTmp_3p1f","bkg_zjetsTmp_3p1f",self.CMS_zz4l_mass,self.mlZjet_3p1f,self.slZjet_3p1f)
            
            self.bkg_zjets = ROOT.RooAddPdf("bkg_zjetsTmp","bkg_zjetsTmp",ROOT.RooArgList(self.bkg_zjets_2p2f,self.bkg_zjets_3p1f),ROOT.RooArgList(self.nlZjet_2p2f,self.nlZjet_3p1f))
            
        elif (self.channel == self.ID_2e2mu):
            name = "mlZjet_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.mlZjet_2p2f = ROOT.RooRealVar(name,"mean landau Zjet 2p2f",self.val_meanL_2P2F)
            name = "slZjet_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.slZjet_2p2f = ROOT.RooRealVar(name,"sigma landau Zjet 2p2f",self.val_sigmaL_2P2F)
            name = "nlZjet_2p2f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.nlZjet_2p2f = ROOT.RooRealVar(name,"norm landau Zjet 2p2f",self.val_normL_2P2F)
            print "mean 2p2f 2mu2e: ",self.mlZjet_2p2f.getVal()
            print "sigma 2p2f 2mu2e: ",self.slZjet_2p2f.getVal()
            print "norm 2p2f 2mu2e: ",self.nlZjet_2p2f.getVal()
            self.bkg_zjets_2p2f = ROOT.RooLandau("bkg_zjetsTmp_2p2f","bkg_zjetsTmp_2p2f",self.CMS_zz4l_mass,self.mlZjet_2p2f,self.slZjet_2p2f)
            
            name = "mlZjet_2p2f_2_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.mlZjet_2p2f_2 = ROOT.RooRealVar(name,"mean landau Zjet 2p2f 2e2mu",self.val_meanL_2P2F_2)
            name = "slZjet_2p2f_2_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.slZjet_2p2f_2 = ROOT.RooRealVar(name,"sigma landau Zjet 2p2f 2e2mu",self.val_sigmaL_2P2F_2)
            name = "nlZjet_2p2f_2_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.nlZjet_2p2f_2 = ROOT.RooRealVar(name,"norm landau Zjet 2p2f 2e2mu",self.val_normL_2P2F_2)
            print "mean 2p2f 2e2mu: ",self.mlZjet_2p2f_2.getVal()
            print "sigma 2p2f 2e2mu: ",self.slZjet_2p2f_2.getVal()
            print "norm 2p2f 2e2mu: ",self.nlZjet_2p2f_2.getVal()
            self.bkg_zjets_2p2f_2 = ROOT.RooLandau("bkg_zjetsTmp_2p2f_2","bkg_zjetsTmp_2p2f_2",self.CMS_zz4l_mass,self.mlZjet_2p2f_2,self.slZjet_2p2f_2)
            
            name = "mlZjet_3p1f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.mlZjet_3p1f = ROOT.RooRealVar(name,"mean landau Zjet 3p1f",self.val_meanL_3P1F)
            name = "slZjet_3p1f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.slZjet_3p1f = ROOT.RooRealVar(name,"sigma landau Zjet 3p1f",self.val_sigmaL_3P1F)
            name = "nlZjet_3p1f_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.nlZjet_3p1f = ROOT.RooRealVar(name,"norm landau Zjet 3p1f",self.val_normL_3P1F)
            print "mean 3p1f 2mu2e: ",self.mlZjet_3p1f.getVal()
            print "sigma 3p1f 2mu2e: ",self.slZjet_3p1f.getVal()
            print "norm 3p1f 2mu2e: ",self.nlZjet_3p1f.getVal()
            self.bkg_zjets_3p1f = ROOT.RooLandau("bkg_zjetsTmp_3p1f","bkg_zjetsTmp_3p1f",self.CMS_zz4l_mass,self.mlZjet_3p1f,self.slZjet_3p1f)
            
            self.bkg_zjets = ROOT.RooAddPdf("bkg_zjetsTmp","bkg_zjetsTmp",ROOT.RooArgList(self.bkg_zjets_2p2f,self.bkg_zjets_3p1f,self.bkg_zjets_2p2f_2),ROOT.RooArgList(self.nlZjet_2p2f,self.nlZjet_3p1f,self.nlZjet_2p2f_2))
            
    ######################

#         ## Reducible backgrounds
#         self.val_meanL = float(self.inputs['zjetsShape_mean'])
#         self.val_sigmaL = float(self.inputs['zjetsShape_sigma'])
#         self.val_poly0L = float(self.inputs['zjetsShape_p0'])
#         self.val_poly1L = float(self.inputs['zjetsShape_p1'])
#         self.val_poly2L = float(self.inputs['zjetsShape_p2'])
#         self.val_poly3L = float(self.inputs['zjetsShape_p3'])
#         self.val_poly4L = float(self.inputs['zjetsShape_p4'])


#         name = "mlZjet_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
#         self.mlZjet = ROOT.RooRealVar(name,"mean landau Zjet",self.val_meanL)
#         name = "slZjet_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
#         self.slZjet = ROOT.RooRealVar(name,"sigma landau Zjet",self.val_sigmaL)
#         name = "bkg_zjetsLandau_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
#         self.bkg_zjets = ROOT.RooLandau(name,"bkg_zjetsTmpLandau",self.CMS_zz4l_mass,self.mlZjet,self.slZjet)

#         name = "p0Zjet_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
#         self.p0Zjet = ROOT.RooRealVar(name,"p0 Zjet",self.val_poly0L)
#         name = "p1Zjet_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
#         self.p1Zjet = ROOT.RooRealVar(name,"p1 Zjet",self.val_poly1L)
#         name = "p2Zjet_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
#         self.p2Zjet = ROOT.RooRealVar(name,"p2 Zjet",self.val_poly2L)
#         name = "p3Zjet_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
#         self.p3Zjet = ROOT.RooRealVar(name,"p3 Zjet",self.val_poly3L)
#         name = "p4Zjet_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
#         self.p4Zjet = ROOT.RooRealVar(name,"p4 Zjet",self.val_poly4L)
#         name = "bkg_zjetsPoly_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
#         self.bkg_zjetsPoly = ROOT.RooPolynomial(name,"bkg_zjetsTmpPoly",self.CMS_zz4l_mass,RooArgList(self.p0Zjet,self.p1Zjet,self.p2Zjet,self.p3Zjet,self.p4Zjet))
#         name = "bkg_zjets_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
#         #self.bkg_zjets = ROOT.RooProdPdf(name,"bkg_zjetsTmp",self.bkg_zjetsLandau,self.bkg_zjetsPoly)
        

        ## ----------------- for massErr ------------------ ##
        self.makeMassErr()


	####  ----------------------- mekd  parametrized double gaussian stuffs  -------------------------
	#if self.bMEKD: 
	#	name = "mekd_qqZZ_a0_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
	#	print "mekd_qqZZ_a0_shape=",self.inputs['mekd_qqZZ_a0_shape'] 
	#	print "mekd_sig_a0_shape=",self.inputs['mekd_sig_a0_shape'] 
	#	mekd_qqZZ_a0 = ROOT.RooFormulaVar(name,"("+self.inputs['mekd_qqZZ_a0_shape']+")", ROOT.RooArgList(CMS_zz4l_mass))
	#	name = "mekd_qqZZ_a1_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
	#	mekd_qqZZ_a1 = ROOT.RooFormulaVar(name,"("+self.inputs['mekd_qqZZ_a1_shape']+")", ROOT.RooArgList(CMS_zz4l_mass))
	#	name = "mekd_qqZZ_a2_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
	#	mekd_qqZZ_a2 = ROOT.RooFormulaVar(name,"("+self.inputs['mekd_qqZZ_a2_shape']+")", ROOT.RooArgList(CMS_zz4l_mass))
	#	name = "mekd_qqZZ_a3_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
	#	mekd_qqZZ_a3 = ROOT.RooFormulaVar(name,"("+self.inputs['mekd_qqZZ_a3_shape']+")", ROOT.RooArgList(CMS_zz4l_mass))
	#	name = "mekd_qqZZ_a4_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
	#	mekd_qqZZ_a4 = ROOT.RooFormulaVar(name,"("+self.inputs['mekd_qqZZ_a4_shape']+")", ROOT.RooArgList(CMS_zz4l_mass))
	#	bkgTemplateMorphPdf_qqzz = ROOT.RooGenericPdf("mekd_qqZZ", "mekd_qqZZ", "@3*exp((-(@0-@1)^2)/(2*@2^2))/@2+(1-@3)*exp((-(@0-@4)^2)/(2*@5^2))/@5", ROOT.RooArgList(MEKD,mekd_qqZZ_a0, mekd_qqZZ_a1, mekd_qqZZ_a2, mekd_qqZZ_a3, mekd_qqZZ_a4))
	#	bkgTemplateMorphPdf_ggzz = ROOT.RooGenericPdf("mekd_ggZZ", "mekd_ggZZ", "@3*exp((-(@0-@1)^2)/(2*@2^2))/@2+(1-@3)*exp((-(@0-@4)^2)/(2*@5^2))/@5", ROOT.RooArgList(MEKD,mekd_qqZZ_a0, mekd_qqZZ_a1, mekd_qqZZ_a2, mekd_qqZZ_a3, mekd_qqZZ_a4))
	#	bkgTemplateMorphPdf_zjets= ROOT.RooGenericPdf("mekd_zjets", "mekd_zjets", "@3*exp((-(@0-@1)^2)/(2*@2^2))/@2+(1-@3)*exp((-(@0-@4)^2)/(2*@5^2))/@5", ROOT.RooArgList(MEKD,mekd_qqZZ_a0, mekd_qqZZ_a1, mekd_qqZZ_a2, mekd_qqZZ_a3, mekd_qqZZ_a4))
	#	m = 100
	#	while m >= 100 and m < 150:
	#		CMS_zz4l_mass.setVal(m)
	#		m = m + 0.1
	#		if mekd_qqZZ_a2.getVal() < 0 : print m, mekd_qqZZ_a2.getVal() 
	#		if mekd_qqZZ_a2.getVal() > 1 : print m, mekd_qqZZ_a2.getVal() 
	#	print "\n \n mekd_qqZZ_a1 channel ",self.channel
	#	m = 100
	#	while m >= 100 and m < 150:
	#		CMS_zz4l_mass.setVal(m)
	#		m = m + 0.1
	#		if mekd_qqZZ_a1.getVal() <= 0 : print m, mekd_qqZZ_a1.getVal() 
	#	print "\n \n mekd_qqZZ_a4 channel ",self.channel
	#	m = 100
	#	while m >= 100 and m < 150:
	#		CMS_zz4l_mass.setVal(m)
	#		m = m + 0.1
	#		if mekd_qqZZ_a4.getVal() <= 0 : print m, mekd_qqZZ_a4.getVal()
        #        print "self.DEBUG Mingshui ", mekd_qqZZ_a0.getVal(), mekd_qqZZ_a1.getVal(), mekd_qqZZ_a2.getVal(), mekd_qqZZ_a3.getVal(), mekd_qqZZ_a4.getVal()
	####  ----------------------- end mekd -----------------------------------------------------------

 
        ## ----------------------- PLOTS FOR SANITY CHECKS -------------------------- ##
        if(self.DEBUG):
            czz = ROOT.TCanvas( "czz", "czz", 750, 700 )
            czz.cd()
            zzframe_s = self.CMS_zz4l_mass.frame(45)
            if self.bUseCBnoConvolution: super(RooDoubleCB,self.signalCB_ggH).plotOn(zzframe_s, ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineColor(1) )
            elif self.isHighMass : super(ROOT.RooFFTConvPdf,self.sig_ggH_HM).plotOn(zzframe_s, ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineColor(1) )
            else : super(ROOT.RooFFTConvPdf,self.sig_ggH).plotOn(zzframe_s, ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineColor(1) )
            super(ROOT.RooqqZZPdf_v2,self.bkg_qqzz).plotOn(zzframe_s, ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineColor(4) )
            super(ROOT.RooggZZPdf_v2,self.bkg_ggzz).plotOn(zzframe_s, ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineColor(6) )
            if (self.channel == self.ID_4mu):
                super(ROOT.RooLandau,self.bkg_zjets).plotOn(zzframe_s, ROOT.RooFit.LineStyle(2), ROOT.RooFit.LineColor(6) )
            else:
                super(ROOT.RooAddPdf,self.bkg_zjets).plotOn(zzframe_s, ROOT.RooFit.LineStyle(2), ROOT.RooFit.LineColor(6) )
            zzframe_s.Draw()
            figName = "{0}/figs/mzz_{1}_{2}.png".format(self.outputDir, self.mH, self.appendName)
            czz.SaveAs(figName)
            del czz
        
        ## ------------------- LUMI -------------------- ##
        self.rrvLumi = ROOT.RooRealVar("cmshzz4l_lumi","cmshzz4l_lumi",self.lumi)  

        
        ## ----------------------- SIGNAL RATES ----------------------- ##
        #self.CMS_zz4l_mass.setRange("shape",self.inputs['low_M'],self.inputs['high_M'])
        self.CMS_zz4l_mass.setRange("shape",121,131)
        
        fr_low_M = self.inputs['low_M']
        fr_high_M = self.inputs['high_M']        
        if (self.mH >= 450): 
            fr_low_M = 100
            fr_high_M = 1000
        if (self.mH >= 750):
            fr_low_M = 100
            fr_high_M = 1400
            


        self.CMS_zz4l_mass.setRange("fullrangesignal",fr_low_M,fr_high_M)
        self.CMS_zz4l_mass.setRange("fullrange",100,1400)
        
        self.rfvCsFilter = RooFormulaVar()
        filterName = "cmshzz4l_csFilter_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        #if(self.sqrts == 7): 
        #    self.rfvCsFilter = ROOT.RooFormulaVar(filterName,"0.5+0.5*TMath::Erf((@0 - 80.85)/50.42)", ROOT.RooArgList(self.MH) )
        #else:
        self.rfvCsFilter = ROOT.RooFormulaVar(filterName,"@0",ROOT.RooArgList(self.one))

        if(self.DEBUG):
            print ">>>>>>  rfvCsFilter = ",self.rfvCsFilter.getVal()

#         sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_a1".format(self.channel,self.sqrts)
#         self.rrva1 = ROOT.RooRealVar(sigEffName,sigEffName, self.inputs['sigEff_a1'])
#         sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_a2".format(self.channel,self.sqrts)
#         self.rrva2 = ROOT.RooRealVar(sigEffName,sigEffName, self.inputs['sigEff_a2'])
#         sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_a3".format(self.channel,self.sqrts)
#         self.rrva3 = ROOT.RooRealVar(sigEffName,sigEffName, self.inputs['sigEff_a3'])
#         sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_a4".format(self.channel,self.sqrts)
#         self.rrva4 = ROOT.RooRealVar(sigEffName,sigEffName, self.inputs['sigEff_a4'])
#         sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_b1".format(self.channel,self.sqrts)
#         self.rrvb1 = ROOT.RooRealVar(sigEffName,sigEffName, self.inputs['sigEff_b1'])
#         sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_b2".format(self.channel,self.sqrts)
#         self.rrvb2 = ROOT.RooRealVar(sigEffName,sigEffName, self.inputs['sigEff_b2'])
#         sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_b3".format(self.channel,self.sqrts)
#         self.rrvb3 = ROOT.RooRealVar(sigEffName,sigEffName, self.inputs['sigEff_b3'])
#         sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_g1".format(self.channel,self.sqrts)
#         self.rrvg1 = ROOT.RooRealVar(sigEffName,sigEffName, self.inputs['sigEff_g1'])
#         sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_g2".format(self.channel,self.sqrts)
#         self.rrvg2 = ROOT.RooRealVar(sigEffName,sigEffName, self.inputs['sigEff_g2'])
#         sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_g3".format(self.channel,self.sqrts)
#         self.rrvg3 = ROOT.RooRealVar(sigEffName,sigEffName, self.inputs['sigEff_g3'])

#         if(self.DEBUG):
#             print "sigEff_a1 = ",self.inputs['sigEff_a1']
#             print "sigEff_a2 = ",self.inputs['sigEff_a2']
#             print "sigEff_a3 = ",self.inputs['sigEff_a3']
#             print "sigEff_a4 = ",self.inputs['sigEff_a4']
#             print "sigEff_b1 = ",self.inputs['sigEff_b1']
#             print "sigEff_b2 = ",self.inputs['sigEff_b2']
#             print "sigEff_b3 = ",self.inputs['sigEff_b3']

    
#         sigEffName = "hzz4lsignaleff_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)

#         self.listEff = ROOT.RooArgList(self.rrva1,self.rrva2,self.rrva3,self.rrva4,self.rrvb1,self.rrvb2,self.rrvb3,self.MH)
#         self.listEff.add(self.rrvg1)
#         self.listEff.add(self.rrvg2)
#         self.listEff.add(self.rrvg3)
#         self.rfvSigEff = ROOT.RooFormulaVar(sigEffName,"(@0+@1*TMath::Erf((@7-@2)/@3))*(@4+@5*@7+@6*@7*@7)+@8*TMath::Gaus(@7,@9,@10)",self.listEff) #ROOT.RooArgList(rrva1,rrva2,rrva3,rrva4,rrvb1,rrvb2,rrvb3,self.MH,rrvg1,rrvg2,rrvg3))
#         #from TF1 *polyFunc= new TF1("polyFunc","([0]+[1]*TMath::Erf( (x-[2])/[3] ))*([4]+[5]*x+[6]*x*x)+[7]*TMath::Gaus(x,[8],[9])", 110., xMax);
        
#         ## following printout is needed ,  dont remove it
#         print ">>>>>>  sigeff ",self.rfvSigEff.getVal()




#########################
        sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_a1".format(self.channel,self.sqrts)
        self.rrva1 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_a1'])
        sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_a2".format(self.channel,self.sqrts)
        self.rrva2 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_a2'])
        sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_a3".format(self.channel,self.sqrts)
        self.rrva3 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_a3'])
        sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_a4".format(self.channel,self.sqrts)
        self.rrva4 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_a4'])
        sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_b1".format(self.channel,self.sqrts)
        self.rrvb1 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_b1'])
        sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_b2".format(self.channel,self.sqrts)
        self.rrvb2 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_b2'])
        sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_b3".format(self.channel,self.sqrts)
        self.rrvb3 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_b3'])
        sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_g1".format(self.channel,self.sqrts)
        self.rrvg1 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_g1'])
        sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_g2".format(self.channel,self.sqrts)
        self.rrvg2 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_g2'])
        sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_g3".format(self.channel,self.sqrts)
        self.rrvg3 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_g3'])
        
        sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_a1".format(self.channel,self.sqrts)
        self.rrva1_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHa1'])
        sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_a2".format(self.channel,self.sqrts)
        self.rrva2_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHa2'])
        sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_a3".format(self.channel,self.sqrts)
        self.rrva3_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHa3'])
        sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_a4".format(self.channel,self.sqrts)
        self.rrva4_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHa4'])
        sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_b1".format(self.channel,self.sqrts)
        self.rrvb1_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHb1'])
        sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_b2".format(self.channel,self.sqrts)
        self.rrvb2_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHb2'])
        sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_b3".format(self.channel,self.sqrts)
        self.rrvb3_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHb3'])
        sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_g1".format(self.channel,self.sqrts)
        self.rrvg1_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHg1'])
        sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_g2".format(self.channel,self.sqrts)
        self.rrvg2_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHg2'])
        sigEffName = "hzz4lqqHeff_{0:.0f}_{1:.0f}_g3".format(self.channel,self.sqrts)
        self.rrvg3_qqh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_qqHg3'])
        
        sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_a1".format(self.channel,self.sqrts)
        self.rrva1_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHa1'])
        sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_a2".format(self.channel,self.sqrts)
        self.rrva2_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHa2'])
        sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_a3".format(self.channel,self.sqrts)
        self.rrva3_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHa3'])
        sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_a4".format(self.channel,self.sqrts)
        self.rrva4_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHa4'])
        sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_b1".format(self.channel,self.sqrts)
        self.rrvb1_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHb1'])
        sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_b2".format(self.channel,self.sqrts)
        self.rrvb2_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHb2'])
        sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_b3".format(self.channel,self.sqrts)
        self.rrvb3_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHb3'])
        sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_g1".format(self.channel,self.sqrts)
        self.rrvg1_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHg1'])
        sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_g2".format(self.channel,self.sqrts)
        self.rrvg2_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHg2'])
        sigEffName = "hzz4lZHeff_{0:.0f}_{1:.0f}_g3".format(self.channel,self.sqrts)
        self.rrvg3_zh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ZHg3'])
        
        sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_a1".format(self.channel,self.sqrts)
        self.rrva1_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHa1'])
        sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_a2".format(self.channel,self.sqrts)
        self.rrva2_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHa2'])
        sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_a3".format(self.channel,self.sqrts)
        self.rrva3_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHa3'])
        sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_a4".format(self.channel,self.sqrts)
        self.rrva4_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHa4'])
        sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_b1".format(self.channel,self.sqrts)
        self.rrvb1_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHb1'])
        sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_b2".format(self.channel,self.sqrts)
        self.rrvb2_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHb2'])
        sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_b3".format(self.channel,self.sqrts)
        self.rrvb3_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHb3'])
        sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_g1".format(self.channel,self.sqrts)
        self.rrvg1_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHg1'])
        sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_g2".format(self.channel,self.sqrts)
        self.rrvg2_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHg2'])
        sigEffName = "hzz4lWHeff_{0:.0f}_{1:.0f}_g3".format(self.channel,self.sqrts)
        self.rrvg3_wh = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_WHg3'])
        
        sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_a1".format(self.channel,self.sqrts)
        self.rrva1_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHa1'])
        sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_a2".format(self.channel,self.sqrts)
        self.rrva2_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHa2'])
        sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_a3".format(self.channel,self.sqrts)
        self.rrva3_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHa3'])
        sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_a4".format(self.channel,self.sqrts)
        self.rrva4_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHa4'])
        sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_b1".format(self.channel,self.sqrts)
        self.rrvb1_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHb1'])
        sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_b2".format(self.channel,self.sqrts)
        self.rrvb2_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHb2'])
        sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_b3".format(self.channel,self.sqrts)
        self.rrvb3_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHb3'])
        sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_g1".format(self.channel,self.sqrts)
        self.rrvg1_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHg1'])
        sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_g2".format(self.channel,self.sqrts)
        self.rrvg2_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHg2'])
        sigEffName = "hzz4lttHeff_{0:.0f}_{1:.0f}_g3".format(self.channel,self.sqrts)
        self.rrvg3_tth = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_ttHg3'])
        
        if self.DEBUG:
            print "sigEff_a1 = ",theInputs['sigEff_a1']
            print "sigEff_a2 = ",theInputs['sigEff_a2']
            print "sigEff_a3 = ",theInputs['sigEff_a3']
            print "sigEff_a4 = ",theInputs['sigEff_a4']
            print "sigEff_b1 = ",theInputs['sigEff_b1']
            print "sigEff_b2 = ",theInputs['sigEff_b2']
            print "sigEff_b3 = ",theInputs['sigEff_b3']
            print "sigEff_g1 = ",theInputs['sigEff_g1']
            print "sigEff_g2 = ",theInputs['sigEff_g2']
            print "sigEff_g3 = ",theInputs['sigEff_g3']
            
            print "sigEff_qqHa1 = ",theInputs['sigEff_qqHa1']
            print "sigEff_qqHa2 = ",theInputs['sigEff_qqHa2']
            print "sigEff_qqHa3 = ",theInputs['sigEff_qqHa3']
            print "sigEff_qqHa4 = ",theInputs['sigEff_qqHa4']
            print "sigEff_qqHb1 = ",theInputs['sigEff_qqHb1']
            print "sigEff_qqHb2 = ",theInputs['sigEff_qqHb2']
            print "sigEff_qqHb3 = ",theInputs['sigEff_qqHb3']
            print "sigEff_qqHg1 = ",theInputs['sigEff_qqHg1']
            print "sigEff_qqHg2 = ",theInputs['sigEff_qqHg2']
            print "sigEff_qqHg3 = ",theInputs['sigEff_qqHg3']
            
            print "sigEff_ZHa1 = ",theInputs['sigEff_ZHa1']
            print "sigEff_ZHa2 = ",theInputs['sigEff_ZHa2']
            print "sigEff_ZHa3 = ",theInputs['sigEff_ZHa3']
            print "sigEff_ZHa4 = ",theInputs['sigEff_ZHa4']
            print "sigEff_ZHb1 = ",theInputs['sigEff_ZHb1']
            print "sigEff_ZHb2 = ",theInputs['sigEff_ZHb2']
            print "sigEff_ZHb3 = ",theInputs['sigEff_ZHb3']
            print "sigEff_ZHg1 = ",theInputs['sigEff_ZHg1']
            print "sigEff_ZHg2 = ",theInputs['sigEff_ZHg2']
            print "sigEff_ZHg3 = ",theInputs['sigEff_ZHg3']
            
            print "sigEff_WHa1 = ",theInputs['sigEff_WHa1']
            print "sigEff_WHa2 = ",theInputs['sigEff_WHa2']
            print "sigEff_WHa3 = ",theInputs['sigEff_WHa3']
            print "sigEff_WHa4 = ",theInputs['sigEff_WHa4']
            print "sigEff_WHb1 = ",theInputs['sigEff_WHb1']
            print "sigEff_WHb2 = ",theInputs['sigEff_WHb2']
            print "sigEff_WHb3 = ",theInputs['sigEff_WHb3']
            print "sigEff_WHg1 = ",theInputs['sigEff_WHg1']
            print "sigEff_WHg2 = ",theInputs['sigEff_WHg2']
            print "sigEff_WHg3 = ",theInputs['sigEff_WHg3']
            
            print "sigEff_ttHa1 = ",theInputs['sigEff_ttHa1']
            print "sigEff_ttHa2 = ",theInputs['sigEff_ttHa2']
            print "sigEff_ttHa3 = ",theInputs['sigEff_ttHa3']
            print "sigEff_ttHa4 = ",theInputs['sigEff_ttHa4']
            print "sigEff_ttHb1 = ",theInputs['sigEff_ttHb1']
            print "sigEff_ttHb2 = ",theInputs['sigEff_ttHb2']
            print "sigEff_ttHb3 = ",theInputs['sigEff_ttHb3']
            print "sigEff_ttHg1 = ",theInputs['sigEff_ttHg1']
            print "sigEff_ttHg2 = ",theInputs['sigEff_ttHg2']
            print "sigEff_ttHg3 = ",theInputs['sigEff_ttHg3']
            
            
        sigEffName_ggH = "hzz4lggHeff_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigEffName_qqH = "hzz4lqqHeff_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigEffName_WH = "hzz4lWHeff_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigEffName_ZH = "hzz4lZHeff_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigEffName_ttH = "hzz4lttHeff_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        
        self.listEff = ROOT.RooArgList(self.rrva1,self.rrva2,self.rrva3,self.rrva4,self.rrvb1,self.rrvb2,self.rrvb3,self.MH)
        self.listEff.add(self.rrvg1)
        self.listEff.add(self.rrvg2)
        self.listEff.add(self.rrvg3)

        self.listEff_qqh = ROOT.RooArgList(self.rrva1_qqh,self.rrva2_qqh,self.rrva3_qqh,self.rrva4_qqh,self.rrvb1_qqh,self.rrvb2_qqh,self.rrvb3_qqh,self.MH)
        self.listEff_qqh.add(self.rrvg1_qqh)
        self.listEff_qqh.add(self.rrvg2_qqh)
        self.listEff_qqh.add(self.rrvg3_qqh)

        self.listEff_wh = ROOT.RooArgList(self.rrva1_wh,self.rrva2_wh,self.rrva3_wh,self.rrva4_wh,self.rrvb1_wh,self.rrvb2_wh,self.rrvb3_wh,self.MH)
        self.listEff_wh.add(self.rrvg1_wh)
        self.listEff_wh.add(self.rrvg2_wh)
        self.listEff_wh.add(self.rrvg3_wh)

        self.listEff_zh = ROOT.RooArgList(self.rrva1_zh,self.rrva2_zh,self.rrva3_zh,self.rrva4_zh,self.rrvb1_zh,self.rrvb2_zh,self.rrvb3_zh,self.MH)
        self.listEff_zh.add(self.rrvg1_zh)
        self.listEff_zh.add(self.rrvg2_zh)
        self.listEff_zh.add(self.rrvg3_zh)

        self.listEff_tth = ROOT.RooArgList(self.rrva1_tth,self.rrva2_tth,self.rrva3_tth,self.rrva4_tth,self.rrvb1_tth,self.rrvb2_tth,self.rrvb3_tth,self.MH)
        self.listEff_tth.add(self.rrvg1_tth)
        self.listEff_tth.add(self.rrvg2_tth)
        self.listEff_tth.add(self.rrvg3_tth)
      
	
	
	#self.rrva1 = rrva1
	#self.rrva2 = rrva2
	#self.rrva3 = rrva3
	#self.rrva4 = rrva4
	#self.rrvb1 = rrvb1
	#self.rrvb2 = rrvb2
	#self.rrvb3 = rrvb3
	#self.rrvg1 = rrvg1
	#self.rrvg2 = rrvg2
	#self.rrvg3 = rrvg3

	#self.listEffRokoTest = ROOT.RooArgList(self.rrva1,self.rrva2,self.rrva3,self.rrva4,self.rrvb1,self.rrvb2,self.rrvb3,self.MH)
        #self.listEffRokoTest.add(self.rrvg1)
        #self.listEffRokoTest.add(self.rrvg2)
        #self.listEffRokoTest.add(self.rrvg3)
	
	#self.listEffRokoTest2 = listEff.clone("listEffRokoTest2")
      
        #self.rfvSigEff_ggH = ROOT.RooFormulaVar(sigEffName_ggH,"(@0+@1*TMath::Erf((@7-@2)/@3))*(@4+@5*@7+@6*@7*@7)+@8*TMath::Gaus(@7,@9,@10)",self.self.listEffRokoTest2) #ROOT.RooArgList(rrva1,rrva2,rrva3,rrva4,rrvb1,rrvb2,rrvb3,self.MH,rrvg1,rrvg2,rrvg3))
        self.rfvSigEff_ggH = ROOT.RooFormulaVar(sigEffName_ggH,"(@0+@1*TMath::Erf((@7-@2)/@3))*(@4+@5*@7+@6*@7*@7)+@8*TMath::Gaus(@7,@9,@10)",self.listEff) #ROOT.RooArgList(rrva1,rrva2,rrva3,rrva4,rrvb1,rrvb2,rrvb3,self.MH,rrvg1,rrvg2,rrvg3))
        self.rfvSigEff_qqH = ROOT.RooFormulaVar(sigEffName_qqH,"(@0+@1*TMath::Erf((@7-@2)/@3))*(@4+@5*@7+@6*@7*@7)+@8*TMath::Gaus(@7,@9,@10)",self.listEff_qqh)
        self.rfvSigEff_ZH = ROOT.RooFormulaVar(sigEffName_ZH,"(@0+@1*TMath::Erf((@7-@2)/@3))*(@4+@5*@7+@6*@7*@7)+@8*TMath::Gaus(@7,@9,@10)",self.listEff_zh)
        self.rfvSigEff_WH = ROOT.RooFormulaVar(sigEffName_WH,"(@0+@1*TMath::Erf((@7-@2)/@3))*(@4+@5*@7+@6*@7*@7)+@8*TMath::Gaus(@7,@9,@10)",self.listEff_wh)
        self.rfvSigEff_ttH = ROOT.RooFormulaVar(sigEffName_ttH,"(@0+@1*TMath::Erf((@7-@2)/@3))*(@4+@5*@7+@6*@7*@7)+@8*TMath::Gaus(@7,@9,@10)",self.listEff_tth)
        #from TF1 *polyFunc= new TF1("polyFunc","([0]+[1]*TMath::Erf( (x-[2])/[3] ))*([4]+[5]*x+[6]*x*x)+[7]*TMath::Gaus(x,[8],[9])", 110., xMax);
        
        ## following printout is needed ,  dont remove it
        print " @@@@@@@@ ggHeff ",self.rfvSigEff_ggH.getVal()
        print " @@@@@@@@ qqHeff ",self.rfvSigEff_qqH.getVal()
        print " @@@@@@@@ ZHeff ",self.rfvSigEff_ZH.getVal()
        print " @@@@@@@@ WHeff ",self.rfvSigEff_WH.getVal()
        print " @@@@@@@@ ttHeff ",self.rfvSigEff_ttH.getVal()


#######################

        #if self.isAltSig:
        #    CS_ggH = self.myCSW.HiggsCS(0,self.mH,self.sqrts)
        #else:
        CS_ggH = self.myCSW.HiggsCS(1,self.mH,self.sqrts)
        CS_VBF = self.myCSW.HiggsCS(2,self.mH,self.sqrts)
        CS_WH = self.myCSW.HiggsCS(3,self.mH,self.sqrts)
        CS_ZH = self.myCSW.HiggsCS(4,self.mH,self.sqrts)
        CS_ttH = self.myCSW.HiggsCS(5,self.mH,self.sqrts)
    
        BRH2e2mu = self.myCSW.HiggsBR(13,self.mH)
        BRH4mu = self.myCSW.HiggsBR(12,self.mH)
        BRH4e = self.myCSW.HiggsBR(12,self.mH)
        self.BR = 0.0
        if( self.channel == self.ID_4mu ): BR = BRH4mu
        if( self.channel == self.ID_4e ): BR = BRH4e
        if( self.channel == self.ID_2e2mu ): BR = BRH2e2mu
    
#        sigEfficiency = self.rfvSigEff.getVal()
        BRZZ = self.myCSW.HiggsBR(11,self.mH)
    
        sigEfficiency_ggH = self.rfvSigEff_ggH.getVal()
        sigEfficiency_qqH = self.rfvSigEff_qqH.getVal()
        sigEfficiency_ZH = self.rfvSigEff_ZH.getVal()
        sigEfficiency_WH = self.rfvSigEff_WH.getVal()
        sigEfficiency_ttH = self.rfvSigEff_ttH.getVal()

        if(self.DEBUG):
            print "CS_ggH: ",CS_ggH,", CS_VBF: ",CS_VBF,", CS_WH: ",CS_WH,", CS_ZH: ",CS_ZH
            print ", CS_ttH: ",CS_ttH,", BRH2e2mu: ",BRH2e2mu,", BRH4mu: ",BRH4mu,", BRH4e: ",BRH4e, ", BRZZ: ",BRZZ

        csCorrection = 1.0
        #if(self.sqrts == 7): csCorrection = self.csFilter(self.mH)

        ## SIG YIELDS
        sigRate_ggH = csCorrection*CS_ggH*BR*sigEfficiency_ggH*1000.*self.lumi
        sigRate_VBF = csCorrection*CS_VBF*BR*sigEfficiency_qqH*1000.*self.lumi
        sigRate_WH = csCorrection*CS_WH*BRZZ*sigEfficiency_WH*1000.*self.lumi
        sigRate_ZH = csCorrection*CS_ZH*BRZZ*sigEfficiency_ZH*1000.*self.lumi
        sigRate_ttH = csCorrection*CS_ttH*BRZZ*sigEfficiency_ttH*1000.*self.lumi
       
        normalizationSignal = self.signalCB_ggH.createIntegral( ROOT.RooArgSet(self.CMS_zz4l_mass), ROOT.RooFit.Range("fullrangesignal") ).getVal()
      
            
        print ">>>>>>  Norm Signal",normalizationSignal
        
        sclFactorSig_ggH = sigRate_ggH/normalizationSignal
        sclFactorSig_VBF = sigRate_VBF/normalizationSignal
        sclFactorSig_WH = sigRate_WH/normalizationSignal
        sclFactorSig_ZH = sigRate_ZH/normalizationSignal
        sclFactorSig_ttH = sigRate_ttH/normalizationSignal

        integral_ggH = 0.0
        integral_VBF = 0.0
        integral_WH  = 0.0
        integral_ZH  = 0.0
        integral_ttH = 0.0

        if self.isHighMass :
            integral_ggH = self.sig_ggH_HM.createIntegral( ROOT.RooArgSet(self.CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        else :
            integral_ggH = self.getVariable(self.signalCB_ggH.createIntegral( ROOT.RooArgSet(self.CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),self.sig_ggH.createIntegral( ROOT.RooArgSet(self.CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),self.bUseCBnoConvolution)

        if self.isHighMass :
            integral_VBF = self.sig_VBF_HM.createIntegral( ROOT.RooArgSet(self.CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        else :
            integral_VBF = self.getVariable(self.signalCB_VBF.createIntegral( ROOT.RooArgSet(self.CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),self.sig_VBF.createIntegral( ROOT.RooArgSet(self.CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),self.bUseCBnoConvolution)

        if self.isHighMass :
            integral_WH = self.sig_WH_HM.createIntegral( ROOT.RooArgSet(self.CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        else :
            integral_WH = self.getVariable(self.signalCB_WH.createIntegral( ROOT.RooArgSet(self.CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),self.sig_WH.createIntegral( ROOT.RooArgSet(self.CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),self.bUseCBnoConvolution)

        if self.isHighMass :
            integral_ZH = self.sig_ZH_HM.createIntegral( ROOT.RooArgSet(self.CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        else :
            integral_ZH = self.getVariable(self.signalCB_ZH.createIntegral( ROOT.RooArgSet(self.CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),self.sig_ZH.createIntegral( ROOT.RooArgSet(self.CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),self.bUseCBnoConvolution)

        if self.isHighMass :
            integral_ttH = self.sig_ttH_HM.createIntegral( ROOT.RooArgSet(self.CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        else :
            integral_ttH = self.getVariable(self.signalCB_ttH.createIntegral( ROOT.RooArgSet(self.CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),self.sig_ttH.createIntegral( ROOT.RooArgSet(self.CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),self.bUseCBnoConvolution)
        
        sigRate_ggH_Shape = sclFactorSig_ggH*integral_ggH
        sigRate_VBF_Shape = sclFactorSig_VBF*integral_VBF
        sigRate_WH_Shape = sclFactorSig_WH*integral_WH
        sigRate_ZH_Shape = sclFactorSig_ZH*integral_ZH
        sigRate_ttH_Shape = sclFactorSig_ttH*integral_ttH
        
        
        normSigName = "cmshzz4l_normalizationSignal_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.rrvNormSig = ROOT.RooRealVar()

        
        if self.isHighMass :
            self.rrvNormSig = ROOT.RooRealVar(normSigName,normSigName, self.sig_ggH_HM.createIntegral(ROOT.RooArgSet(self.CMS_zz4l_mass)).getVal())
        else :
            self.rrvNormSig = ROOT.RooRealVar(normSigName,normSigName, self.getVariable(self.signalCB_ggH.createIntegral(ROOT.RooArgSet(self.CMS_zz4l_mass)).getVal(),self.sig_ggH.createIntegral(ROOT.RooArgSet(self.CMS_zz4l_mass)).getVal(),self.bUseCBnoConvolution))
        self.rrvNormSig.setConstant(True)

      ##rfvSigRate_ggH = ROOT.RooFormulaVar("ggH_norm","@0*@1*@2*1000*{0}*{2}/{1}".format(self.lumi,rrvNormSig.getVal(),self.getVariable(signalCB_ggH.createIntegral(RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),sig_ggH.createIntegral(RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),self.bUseCBnoConvolution)),ROOT.RooArgList(rfvCsFilter,rfvSigEff, rhfXsBrFuncV_1))

        #self.rfvSigRate_ggH = ROOT.RooFormulaVar("ggH_norm","@0*@1*@2*1000*{0}*{2}/{1}".format(self.lumi,self.rrvNormSig.getVal(),integral_ggH),ROOT.RooArgList(self.rfvCsFilter,self.rfvSigEff, self.rhfXsBrFuncV_1))

        print ">>>>>> Compare Integrals: integral_ggH=",integral_ggH,"  ; calculated=",self.getVariable(self.signalCB_ggH.createIntegral(RooArgSet(self.CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),self.sig_ggH.createIntegral(RooArgSet(self.CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),self.bUseCBnoConvolution)
        
        #self.rfvSigRate_VBF = ROOT.RooFormulaVar("qqH_norm","@0*@1*@2*1000*{0}*{2}/{1}".format(self.lumi,self.rrvNormSig.getVal(),integral_VBF),
						   #ROOT.RooArgList(self.rfvCsFilter,self.rfvSigEff_qqH, self.rhfXsBrFuncV_2))
                         

        #self.rfvSigRate_WH = ROOT.RooFormulaVar("WH_norm","@0*@1*@2*1000*{0}*{2}/{1}".format(self.lumi,self.rrvNormSig.getVal(),integral_WH),ROOT.RooArgList(self.rfvCsFilter,self.rfvSigEff_WH, self.rhfXsBrFuncV_3))
        
        
        
        #self.rfvSigRate_ZH = ROOT.RooFormulaVar("ZH_norm","@0*@1*@2*1000*{0}*{2}/{1}".format(self.lumi,self.rrvNormSig.getVal(),integral_ZH),ROOT.RooArgList(self.rfvCsFilter,self.rfvSigEff_ZH, self.rhfXsBrFuncV_4))
                         

        #self.rfvSigRate_ttH = ROOT.RooFormulaVar("ttH_norm","@0*@1*@2*1000*{0}*{2}/{1}".format(self.lumi,self.rrvNormSig.getVal(),integral_ttH),ROOT.RooArgList(self.rfvCsFilter,self.rfvSigEff_ttH, self.rhfXsBrFuncV_5))
                         
        #if self.DEBUG:
            #print self.signalCB_ggH.createIntegral(ROOT.RooArgSet(self.CMS_zz4l_mass)).getVal(),"   ",self.sig_ggH.createIntegral(ROOT.RooArgSet(self.CMS_zz4l_mass)).getVal()
            #print self.signalCB_ggH.createIntegral(ROOT.RooArgSet(self.CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),"   ",self.sig_ggH.createIntegral(ROOT.RooArgSet(self.CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal()


        #self.rfvSigRate_ggH_temp = ROOT.RooFormulaVar("ggH_norm_temp","@0*@1*@2*1000*{0}*{2}/{1}".format(self.lumi,self.rrvNormSig.getVal(),integral_ggH),ROOT.RooArgList(self.rfvCsFilter,self.rfvSigEff_ggH, self.rhfXsBrFuncV_1))
        
        
        
        
        #self.rfvSigRate_VBF = ROOT.RooFormulaVar("qqH_norm","@0*@1*@2*1000*@3*{1}/{0}".format(self.rrvNormSig.getVal(),integral_VBF), ROOT.RooArgList(self.rfvCsFilter,self.rfvSigEff_qqH, self.rhfXsBrFuncV_2,self.rrvLumi))
        #self.rfvSigRate_WH = ROOT.RooFormulaVar("WH_norm","@0*@1*@2*1000*@3*{1}/{0}".format(self.rrvNormSig.getVal(),integral_WH),ROOT.RooArgList(self.rfvCsFilter,self.rfvSigEff_WH, self.rhfXsBrFuncV_3,self.rrvLumi))
        #self.rfvSigRate_ZH = ROOT.RooFormulaVar("ZH_norm","@0*@1*@2*1000*@3*{1}/{0}".format(self.rrvNormSig.getVal(),integral_ZH),ROOT.RooArgList(self.rfvCsFilter,self.rfvSigEff_ZH, self.rhfXsBrFuncV_4,self.rrvLumi))
        #self.rfvSigRate_ttH = ROOT.RooFormulaVar("ttH_norm","@0*@1*@2*1000*@3*{1}/{0}".format(self.rrvNormSig.getVal(),integral_ttH),ROOT.RooArgList(self.rfvCsFilter,self.rfvSigEff_ttH, self.rhfXsBrFuncV_5,self.rrvLumi))
        #self.rfvSigRate_ggH_temp = ROOT.RooFormulaVar("ggH_norm_temp","@0*@1*@2*1000*@3*{1}/{0}".format(self.rrvNormSig.getVal(),integral_ggH),ROOT.RooArgList(self.rfvCsFilter,self.rfvSigEff_ggH, self.rhfXsBrFuncV_1,self.rrvLumi))
        postfix_norm = "{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        
        self.rfvSigRate_VBF = ROOT.RooFormulaVar("qqH_norm_{0}".format(postfix_norm),"@0*@1*@2*1000*@3*{1}/{0}".format(self.rrvNormSig.getVal(),integral_VBF), ROOT.RooArgList(self.rfvCsFilter,self.rfvSigEff_qqH, self.rhfXsBrFuncV_2,self.rrvLumi))
        self.rfvSigRate_WH = ROOT.RooFormulaVar("WH_norm_{0}".format(postfix_norm),"@0*@1*@2*1000*@3*{1}/{0}".format(self.rrvNormSig.getVal(),integral_WH),ROOT.RooArgList(self.rfvCsFilter,self.rfvSigEff_WH, self.rhfXsBrFuncV_3,self.rrvLumi))
        self.rfvSigRate_ZH = ROOT.RooFormulaVar("ZH_norm_{0}".format(postfix_norm),"@0*@1*@2*1000*@3*{1}/{0}".format(self.rrvNormSig.getVal(),integral_ZH),ROOT.RooArgList(self.rfvCsFilter,self.rfvSigEff_ZH, self.rhfXsBrFuncV_4,self.rrvLumi))
        self.rfvSigRate_ttH = ROOT.RooFormulaVar("ttH_norm_{0}".format(postfix_norm),"@0*@1*@2*1000*@3*{1}/{0}".format(self.rrvNormSig.getVal(),integral_ttH),ROOT.RooArgList(self.rfvCsFilter,self.rfvSigEff_ttH, self.rhfXsBrFuncV_5,self.rrvLumi))
        self.rfvSigRate_ggH_temp = ROOT.RooFormulaVar("ggH_norm_temp_{0}".format(postfix_norm),"@0*@1*@2*1000*@3*{1}/{0}".format(self.rrvNormSig.getVal(),integral_ggH),ROOT.RooArgList(self.rfvCsFilter,self.rfvSigEff_ggH, self.rhfXsBrFuncV_1,self.rrvLumi))
        
        

        ##set ggH yield to output of jhuGen and correct for vbf+zh+wh+tth
        self.rrvJHUgen_SMggH = ROOT.RooRealVar("jhuGen_SM","jhuGen_SM",float(theInputs['jhuGen_SM_yield']))
        self.rrv_SMggH_ratio = ROOT.RooRealVar("ggH_ratio","ggH_ratio",(self.rfvSigRate_ggH_temp.getVal()+self.rfvSigRate_VBF.getVal()+self.rfvSigRate_WH.getVal()+self.rfvSigRate_ZH.getVal()+self.rfvSigRate_ttH.getVal())/self.rfvSigRate_ggH_temp.getVal())
        self.rfvSigRate_ggH = ROOT.RooFormulaVar("ggH_norm","@0",ROOT.RooArgList(self.one))
                
        if self.all_chan:
            print " >>>>>> Requested to sum up over the 5 chans: the norm in rfvSigRate_ggH should be the sum of the values of sigRate_XYZ_Shape variables:"
        print ">>>>>> rfvSigRate_ggH_temp = ",self.rfvSigRate_ggH_temp.getVal()
        print ">>>>>> Norm Sig = ",self.rrvNormSig.getVal()
        print ">>>>>> rfvSigRate_ggH = ",self.rfvSigRate_ggH.getVal()
        print ">>>>>> sigRate_ggH_Shape = ",sigRate_ggH_Shape
        print ">>>>>> rfvSigRate_VBF = ",self.rfvSigRate_VBF.getVal()
        print ">>>>>> sigRate_VBF_Shape = ",sigRate_VBF_Shape
        print ">>>>>> rfvSigRate_WH = ",self.rfvSigRate_WH.getVal()
        print ">>>>>> sigRate_WH_Shape = ",sigRate_WH_Shape
        print ">>>>>> rfvSigRate_ZH = ",self.rfvSigRate_ZH.getVal()
        print ">>>>>> sigRate_ZH_Shape = ",sigRate_ZH_Shape
        print ">>>>>> rfvSigRate_ttH = ",self.rfvSigRate_ttH.getVal()
        print ">>>>>> sigRate_ttH_Shape = ",sigRate_ttH_Shape
        print ">>>>>> Sum of sigRate_XYZ_Shape =" ,sigRate_ggH_Shape+sigRate_VBF_Shape+sigRate_WH_Shape+sigRate_ZH_Shape+sigRate_ttH_Shape
        print ">>>>>> rrv_SMggH_ratio = ",self.rrv_SMggH_ratio.getVal()
        print ">>>>>> rrvJHUgen_SMggH = ",self.rrvJHUgen_SMggH.getVal()
             
        ## ----------------------- BACKGROUND RATES ----------------------- ##
        ## rates per lumi for scaling
        bkgRate_qqzz = self.inputs['qqZZ_rate']/self.inputs['qqZZ_lumi']
        bkgRate_ggzz = self.inputs['ggZZ_rate']/self.inputs['qqZZ_lumi']
        bkgRate_zjets = self.inputs['zjets_rate']/self.inputs['zjets_lumi']
        
        ## Get Normalizations
        normalizationBackground_qqzz = self.bkg_qqzz.createIntegral( ROOT.RooArgSet(self.CMS_zz4l_mass), ROOT.RooFit.Range("fullrange") ).getVal()
        normalizationBackground_ggzz = self.bkg_ggzz.createIntegral( ROOT.RooArgSet(self.CMS_zz4l_mass), ROOT.RooFit.Range("fullrange") ).getVal()
        normalizationBackground_zjets = self.bkg_zjets.createIntegral( ROOT.RooArgSet(self.CMS_zz4l_mass), ROOT.RooFit.Range("fullrange") ).getVal()
        
        sclFactorBkg_qqzz = self.lumi*bkgRate_qqzz/normalizationBackground_qqzz
        sclFactorBkg_ggzz = self.lumi*bkgRate_ggzz/normalizationBackground_ggzz
        sclFactorBkg_zjets = self.lumi*bkgRate_zjets/normalizationBackground_zjets
               
        bkgRate_qqzz_Shape = sclFactorBkg_qqzz * self.bkg_qqzz.createIntegral( ROOT.RooArgSet(self.CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        bkgRate_ggzz_Shape = sclFactorBkg_ggzz * self.bkg_ggzz.createIntegral( ROOT.RooArgSet(self.CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        bkgRate_zjets_Shape = sclFactorBkg_zjets * self.bkg_zjets.createIntegral( ROOT.RooArgSet(self.CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        
        if(self.DEBUG):
            print "Shape signal rate: ",sigRate_ggH_Shape,", background rate: ",bkgRate_qqzz_Shape,", ",bkgRate_zjets_Shape," in ",self.inputs['low_M']," - ",self.inputs['high_M']
            self.CMS_zz4l_mass.setRange("lowmassregion",100.,160.)
            bkgRate_qqzz_lowmassregion = sclFactorBkg_qqzz * self.bkg_qqzz.createIntegral( ROOT.RooArgSet(self.CMS_zz4l_mass), ROOT.RooFit.Range("lowmassregion") ).getVal()
            bkgRate_ggzz_lowmassregion = sclFactorBkg_ggzz * self.bkg_ggzz.createIntegral( ROOT.RooArgSet(self.CMS_zz4l_mass), ROOT.RooFit.Range("lowmassregion") ).getVal()
            bkgRate_zjets_lowmassregion = sclFactorBkg_zjets * self.bkg_zjets.createIntegral( ROOT.RooArgSet(self.CMS_zz4l_mass), ROOT.RooFit.Range("lowmassregion") ).getVal()
            lowmassyield = bkgRate_qqzz_lowmassregion + bkgRate_ggzz_lowmassregion + bkgRate_zjets_lowmassregion
            print "low mass yield: ",lowmassyield
        

        self.systematics.setSystematics(bkgRate_qqzz_Shape,bkgRate_ggzz_Shape,bkgRate_zjets_Shape)
        self.systematics_forXSxBR.setSystematics(bkgRate_qqzz_Shape,bkgRate_ggzz_Shape,bkgRate_zjets_Shape)

        if self.isAltSig:
            #            sigRate_ggH_Shape = self.rfvSigRate_ggH.getVal()
            #            print ">>>>>> Forcing sigRate_ggH_Shape = ",self.rfvSigRate_ggH.getVal()
            if self.inputs['altHypothesis'] == 'gg0-' or self.inputs['altHypothesis'] == 'gg0h+':
                sigRate_ggH_Shape = self.rrvJHUgen_SMggH.getVal()*self.rrv_SMggH_ratio.getVal()
            else:
                sigRate_ggH_Shape = self.rfvSigRate_ggH_temp.getVal()+self.rfvSigRate_VBF.getVal()+self.rfvSigRate_WH.getVal()+self.rfvSigRate_ZH.getVal()+self.rfvSigRate_ttH.getVal()


        ## SET RATES TO 1 
        ## DC RATES WILL BE MULTIPLIED
        ## BY RATES IMPORTED TO WS
      
	#if not self.inputs['unfold']:#mod-roko
	    #sigRate_ggH_Shape = 1
	    #sigRate_VBF_Shape = 1
	    #sigRate_WH_Shape = 1
	    #sigRate_ZH_Shape = 1
	    #sigRate_ttH_Shape = 1
	    ##end mod-roko 
        
        ## If the channel is not declared in inputs, set rate = 0
        if not self.ggH_chan and not self.all_chan:  sigRate_ggH_Shape = 0
        if not self.qqH_chan:  sigRate_VBF_Shape = 0
        if not self.WH_chan:   sigRate_WH_Shape = 0
        if not self.ZH_chan:   sigRate_ZH_Shape = 0
        if not self.ttH_chan:  sigRate_ttH_Shape = 0
        
        if not self.qqZZ_chan and not self.all_chan:  bkgRate_qqzz_Shape = 0
        if not self.ggZZ_chan and not self.all_chan:  bkgRate_ggzz_Shape = 0
        if not self.zjets_chan and not self.all_chan: bkgRate_zjets_Shape = 0
        
        self.rates = {}
        self.rates['ggH'] = sigRate_ggH_Shape
        self.rates['qqH'] = sigRate_VBF_Shape
        self.rates['WH']  = sigRate_WH_Shape
        self.rates['ZH']  = sigRate_ZH_Shape
        self.rates['ttH'] = sigRate_ttH_Shape
        
        self.rates['qqZZ']  = bkgRate_qqzz_Shape
        self.rates['ggZZ']  = bkgRate_ggzz_Shape
        self.rates['zjets'] = bkgRate_zjets_Shape
        self.rates['ttbar'] = 0
        self.rates['zbb']   = 0
        
        

    ## --------------------------- DATASET --------------------------- ##
    def fetchDataset(self):

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
        datasetName = "data_obs_{0}".format(self.appendName)
            
        if(self.bIncludingError):
            self.data_obs = ROOT.RooDataSet(datasetName,datasetName,self.data_obs_tree,ROOT.RooArgSet(self.CMS_zz4l_mass, self.RelErr))
        else:
            self.data_obs = ROOT.RooDataSet(datasetName,datasetName,self.data_obs_tree,ROOT.RooArgSet(self.CMS_zz4l_mass))
		


    def endsInP5(self,massH):
        EndsInP5 = False
        tmpMH = massH
        if (math.fabs(math.floor(tmpMH)-massH) > 0.001): EndsInP5 = True
        if (self.DEBUG): print "ENDS IN P5  ",EndsInP5
        return EndsInP5



    ## --------------------------- WORKSPACE -------------------------- ##
    def writeWorkspace(self):
            
        self.name_ShapeWS2 = ""
        self.name_ShapeWSXSBR = ""
        
        if (self.endsInP5(self.mH)):
            self.name_ShapeWS = "{0}/HCG/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts)
        else:
            self.name_ShapeWS = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts)
            
        if (self.endsInP5(self.mH)):
            self.name_ShapeWSXSBR = "{0}/HCG_XSxBR/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts)
        else:
            self.name_ShapeWSXSBR = "{0}/HCG_XSxBR/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts)

        self.name_ShapeWS2 = "hzz4l_{0}S_{1:.0f}TeV.input.root".format(self.appendName,self.sqrts)
        
        self.w = ROOT.RooWorkspace("w","w")
        
        self.w.importClassCode(RooqqZZPdf_v2.Class(),True)
        self.w.importClassCode(RooggZZPdf_v2.Class(),True)
        self.w.importClassCode(RooRelBWUFParam.Class(),True)
        self.w.importClassCode(RooDoubleCB.Class(),True)
        self.w.importClassCode(RooFormulaVar.Class(),True)
        if self.isHighMass :
            self.w.importClassCode(RooRelBWHighMass.Class(),True)

        getattr(self.w,'import')(self.data_obs,ROOT.RooFit.Rename("data_obs")) ### Should this be renamed?


        getattr(self.w,'import')(self.rfvSigRate_ggH, ROOT.RooFit.RecycleConflictNodes())  #mod-roko
        getattr(self.w,'import')(self.rfvSigRate_VBF, ROOT.RooFit.RecycleConflictNodes())
        getattr(self.w,'import')(self.rfvSigRate_WH, ROOT.RooFit.RecycleConflictNodes())
        getattr(self.w,'import')(self.rfvSigRate_ZH, ROOT.RooFit.RecycleConflictNodes())
        getattr(self.w,'import')(self.rfvSigRate_ttH, ROOT.RooFit.RecycleConflictNodes())
        
                                                

    
        if(self.bUseCBnoConvolution) :
            if not self.bIncludingError:
                self.signalCB_ggH.SetNameTitle("ggH","ggH")
                self.signalCB_VBF.SetNameTitle("qqH","qqH")
                self.signalCB_WH.SetNameTitle("WH","WH")
                self.signalCB_ZH.SetNameTitle("ZH","ZH")
                self.signalCB_ttH.SetNameTitle("ttH","ttH")
                
                getattr(self.w,'import')(self.signalCB_ggH, ROOT.RooFit.RecycleConflictNodes())
                getattr(self.w,'import')(self.signalCB_VBF, ROOT.RooFit.RecycleConflictNodes())
                getattr(self.w,'import')(self.signalCB_WH, ROOT.RooFit.RecycleConflictNodes())
                getattr(self.w,'import')(self.signalCB_ZH, ROOT.RooFit.RecycleConflictNodes())
                getattr(self.w,'import')(self.signalCB_ttH, ROOT.RooFit.RecycleConflictNodes())
            else:
                self.sig_ggHErr.SetNameTitle("ggH","ggH")
                self.sig_VBFErr.SetNameTitle("qqH","qqH")
                self.sig_WHErr.SetNameTitle("WH","WH")
                self.sig_ZHErr.SetNameTitle("ZH","ZH")
                self.sig_ttHErr.SetNameTitle("ttH","ttH")
                
                getattr(self.w,'import')(self.sig_ggHErr, ROOT.RooFit.RecycleConflictNodes())
                getattr(self.w,'import')(self.sig_VBFErr, ROOT.RooFit.RecycleConflictNodes())
                getattr(self.w,'import')(self.sig_WHErr, ROOT.RooFit.RecycleConflictNodes())
                getattr(self.w,'import')(self.sig_ZHErr, ROOT.RooFit.RecycleConflictNodes())
                getattr(self.w,'import')(self.sig_ttHErr, ROOT.RooFit.RecycleConflictNodes())
                        
        else:
            
            if self.isHighMass:
                self.sig_ggH_HM.SetNameTitle("ggH","ggH")
                self.sig_VBF_HM.SetNameTitle("qqH","qqH")
                self.sig_WH_HM.SetNameTitle("WH","WH")
                self.sig_ZH_HM.SetNameTitle("ZH","ZH")
                self.sig_ttH_HM.SetNameTitle("ttH","ttH")
                
                getattr(self.w,'import')(self.sig_ggH_HM, ROOT.RooFit.RecycleConflictNodes())
                getattr(self.w,'import')(self.sig_VBF_HM, ROOT.RooFit.RecycleConflictNodes())
                getattr(self.w,'import')(self.sig_WH_HM, ROOT.RooFit.RecycleConflictNodes())
                getattr(self.w,'import')(self.sig_ZH_HM, ROOT.RooFit.RecycleConflictNodes())
                getattr(self.w,'import')(self.sig_ttH_HM, ROOT.RooFit.RecycleConflictNodes())
                
            else:
                self.sig_ggH.SetNameTitle("ggH","ggH")
                self.sig_VBF.SetNameTitle("qqH","qqH")
                self.sig_WH.SetNameTitle("WH","WH")
                self.sig_ZH.SetNameTitle("ZH","ZH")
                self.sig_ttH.SetNameTitle("ttH","ttH")
                
                getattr(self.w,'import')(self.sig_ggH, ROOT.RooFit.RecycleConflictNodes())
                getattr(self.w,'import')(self.sig_VBF, ROOT.RooFit.RecycleConflictNodes())
                getattr(self.w,'import')(self.sig_WH, ROOT.RooFit.RecycleConflictNodes())
                getattr(self.w,'import')(self.sig_ZH, ROOT.RooFit.RecycleConflictNodes())
                getattr(self.w,'import')(self.sig_ttH, ROOT.RooFit.RecycleConflictNodes())
                    
                    
        if not self.bIncludingError:
            self.bkg_qqzz.SetNameTitle("bkg_qqzz","bkg_qqzz")
            self.bkg_ggzz.SetNameTitle("bkg_ggzz","bkg_ggzz")
            self.bkg_zjets.SetNameTitle("bkg_zjets","bkg_zjets")
            getattr(self.w,'import')(self.bkg_qqzz, ROOT.RooFit.RecycleConflictNodes())
            getattr(self.w,'import')(self.bkg_ggzz, ROOT.RooFit.RecycleConflictNodes())
            getattr(self.w,'import')(self.bkg_zjets, ROOT.RooFit.RecycleConflictNodes())
        else:
            self.bkg_qqzzErr.SetNameTitle("bkg_qqzz","bkg_qqzz")
            self.bkg_ggzzErr.SetNameTitle("bkg_ggzz","bkg_ggzz")
            self.bkg_zjetsErr.SetNameTitle("bkg_zjets","bkg_zjets")
            getattr(self.w,'import')(self.bkg_qqzzErr, ROOT.RooFit.RecycleConflictNodes())
            getattr(self.w,'import')(self.bkg_ggzzErr, ROOT.RooFit.RecycleConflictNodes())
            getattr(self.w,'import')(self.bkg_zjetsErr, ROOT.RooFit.RecycleConflictNodes())
            

	if self.DEBUG : self.w.Print()
        self.w.writeToFile(self.name_ShapeWS)
        self.w.writeToFile(self.name_ShapeWSXSBR)
        
        

    ## --------------------------- DATACARDS -------------------------- ##
    def prepareDatacard(self):

        name_Shape = ""
        name_ShapeWS = ""

        if (self.endsInP5(self.mH)):
            name_Shape = "{0}/HCG/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts)
        else:
            name_Shape = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts)


        ## Write Datacards
        fo = open( name_Shape, "wb")
        self.WriteDatacard(fo, self.name_ShapeWS2, self.rates, self.data_obs.numEntries(), self.is2D )
        self.systematics.WriteSystematics(fo, self.inputs)
        self.systematics.WriteShapeSystematics(fo,self.inputs)
        fo.close()
        
        ## forXSxBR
        if (self.endsInP5(self.mH)):
            name_Shape = "{0}/HCG_XSxBR/{2:.1f}/hzz4l_{1}S_{3:.0f}TeV.txt".format(self.outputDir,self.appendName,self.mH,self.sqrts)	
        else:
            name_Shape = "{0}/HCG_XSxBR/{2:.0f}/hzz4l_{1}S_{3:.0f}TeV.txt".format(self.outputDir,self.appendName,self.mH,self.sqrts)
            
        fo = open( name_Shape, "wb" )
        self.WriteDatacard(fo,self.name_ShapeWS2, self.rates, self.data_obs.numEntries(), self.is2D )
        self.systematics_forXSxBR.WriteSystematics(fo, self.inputs)
        self.systematics_forXSxBR.WriteShapeSystematics(fo,self.inputs)
        fo.close()

            
    def WriteDatacard(self,file,nameWS,theRates,obsEvents,is2D,isAltCard=False,AltLabel=""):

        numberSig = self.numberOfSigChan(self.inputs)
        numberBg  = self.numberOfBgChan(self.inputs)

        print "DEBUG: Number of signals={0} and backgrounds={1}. IsAltHyp={3} ".format(numberSig,numberOfBgChan,self.isAltSig)
        
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
        channelName1D=['ggH','qqH','WH','ZH','ttH','bkg_qqzz','bkg_ggzz','bkg_zjets','bkg_ttbar','bkg_zbb']

        if self.inputs["all"]:
            channelList=['ggH','qqZZ','ggZZ','zjets','ttbar','zbb']
            if isAltCard :
		
		print 'THIS IS AN ALTERNATIVE CARD !!!!'
		channelList=['ggH','ggH','qqZZ','ggZZ','zjets','ttbar','zbb']
		channelName1D=['ggH','ggH{0}'.format(AltLabel),'bkg_qqzz','bkg_ggzz','bkg_zjets','bkg_ttbar','bkg_zbb']
		channelName2D=['ggH','ggH{0}'.format(AltLabel),'bkg2d_qqzz','bkg2d_ggzz','bkg2d_zjets','bkg2d_ttbar','bkg2d_zbb']
		channelName2D=[      'ggH{0}'.format(AltLabel),'bkg2d_qqzz','bkg2d_ggzz','bkg2d_zjets','bkg2d_ttbar','bkg2d_zbb']
            else:
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
            if self.inputs[chan]:
                file.write("{0} ".format(channelName1D[i]))
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
            if self.inputs[chan] or (chan.startswith("ggH") and self.inputs["all"]):
                file.write("{0:.4f} ".format(theRates[chan]))
        file.write("\n")
        file.write("------------\n")


        
    def numberOfSigChan(self,myinputs):

        counter=0

        if myinputs['ggH']: counter+=1
        if myinputs['qqH']: counter+=1
        if myinputs['WH']:  counter+=1
        if myinputs['ZH']:  counter+=1
        if myinputs['ttH']: counter+=1
        if myinputs['all']: counter+=1
        if self.isAltSig:
            if myinputs['all']: counter+=1  #mod-roko
            #begin mod-roko  -- added to count different contributions to X width
            #if myinputs['ggH']: counter+=1
            #if myinputs['gg0Ph']: counter+=1
            #if myinputs['ggInt_13P']: counter+=1
            #if myinputs['ggInt_13N']: counter+=1
            #end mod-roko
            
        else:
            if myinputs['all']: counter+=2
               
        return counter

    def numberOfBgChan(self,myinputs):

        counter=0

        if myinputs['qqZZ']:  counter+=1
        if myinputs['ggZZ']:  counter+=1
        if myinputs['zjets']: counter+=1
        if myinputs['ttbar']: counter+=1
        if myinputs['zbb']:   counter+=1

        if self.isAltSig:
            if myinputs['all']: counter+=3
        else:
            if myinputs['all']: counter+=5
                    
        
        return counter



    def defineMassErr(self):

        ## -------------------- DEFINE VARIABLES ------------------- ##
        mrelerrVarName = "CMS_zz4l_massRelErr"
        self.RelErr = ROOT.RooRealVar(mrelerrVarName,mrelerrVarName,0.002,0.2)
        self.RelErr.setBins(100)
        self.CMS_zz4l_massErr = ROOT.RooFormulaVar()
        if (self.channel == self.ID_4mu) :
            self.CMS_zz4l_massErr = ROOT.RooFormulaVar("CMS_zz4l_massErr","@0*@1*(1+@2)", ROOT.RooArgList(self.CMS_zz4l_mass, self.RelErr, self.CMS_zz4l_sigma_m_sig))
        elif (self.channel == self.ID_4e) :
            self.CMS_zz4l_massErr = ROOT.RooFormulaVar("CMS_zz4l_massErr","@0*@1*(1+@2)", ROOT.RooArgList(self.CMS_zz4l_mass, self.RelErr, self.CMS_zz4l_sigma_e_sig))
        elif (self.channel == self.ID_2e2mu) :
            self.CMS_zz4l_massErr = ROOT.RooFormulaVar("CMS_zz4l_massErr","@0*@1*TMath::Sqrt((1+@2)*(1+@3))", ROOT.RooArgList(self.CMS_zz4l_mass, self.RelErr, self.CMS_zz4l_sigma_m_sig, self.CMS_zz4l_sigma_e_sig))


    def makeMassErr(self):
        
        name = "CMS_zz4l_massErrS_ln_kappa_{0:.0f}".format(self.channel)
	self.rfv_EBE_sig_ln_kappa = ROOT.RooFormulaVar(name, "("+self.inputs['relerr_ggH_gs_sigma']+")", ROOT.RooArgList(self.MH))
	name = "CMS_zz4l_massErrS_ln_mean_{0:.0f}".format(self.channel)
	self.rfv_EBE_sig_ln_mean = ROOT.RooFormulaVar(name, "("+self.inputs['relerr_ggH_gs_mean']+")", ROOT.RooArgList(self.MH))
	self.EBE_sig_ln = ROOT.RooLognormal("errLN_ggH","errLN_ggH", self.RelErr, self.rfv_EBE_sig_ln_mean, self.rfv_EBE_sig_ln_kappa)
	if self.channel!=1:
            self.EBE_sig_ln = ROOT.RooGaussian("errGaus_ggH","errGaus_ggH", self.RelErr, self.rfv_EBE_sig_ln_mean, self.rfv_EBE_sig_ln_kappa)
	name = "CMS_zz4l_massErrS_ld_sigma_{0:.0f}".format(self.channel)
	self.rfv_EBE_sig_ld_sigma = ROOT.RooFormulaVar(name, "("+self.inputs['relerr_ggH_ld_mean']+")", ROOT.RooArgList(self.MH))
	name = "CMS_zz4l_massErrS_ld_mean_{0:.0f}".format(self.channel)
	self.rfv_EBE_sig_ld_mean = ROOT.RooFormulaVar(name, "("+self.inputs['relerr_ggH_ld_sigma']+")", ROOT.RooArgList(self.MH)) 
	self.EBE_sig_ld = ROOT.RooLandau("errLD_ggH","errLD_ggH", self.RelErr, self.rfv_EBE_sig_ld_mean, self.rfv_EBE_sig_ld_sigma)
	name = "CMS_zz4l_massErrS_ld_frac_{0:.0f}".format(self.channel)
	self.rfv_EBE_sig_frac = ROOT.RooFormulaVar(name, "("+self.inputs['relerr_ggH_ld_frac']+")", ROOT.RooArgList(self.MH)) 
	self.pdfErrS = ROOT.RooAddPdf("pdfErrS","pdfErrS", self.EBE_sig_ld, self.EBE_sig_ln, self.rfv_EBE_sig_frac)

	name = "CMS_zz4l_massErrZZ_ln_kappa_{0:.0f}".format(self.channel)
	self.rfv_EBE_zz_ln_kappa = ROOT.RooFormulaVar(name, "("+self.inputs['relerr_qqzz_gs_sigma']+")", ROOT.RooArgList(self.CMS_zz4l_mass)) 
	name = "CMS_zz4l_massErrZZ_ln_mean_{0:.0f}".format(self.channel)
	self.rfv_EBE_zz_ln_mean = ROOT.RooFormulaVar(name, "("+self.inputs['relerr_qqzz_gs_mean']+")", ROOT.RooArgList(self.CMS_zz4l_mass))
	self.EBE_zz_ln = ROOT.RooLognormal("errLN_qqzz","errLN_qqzz", self.RelErr, self.rfv_EBE_zz_ln_mean, self.rfv_EBE_zz_ln_kappa)
	if self.channel!=1:
            self.EBE_zz_ln = ROOT.RooGaussian("errGaus_qqzz","errGaus_qqzz", self.RelErr, self.rfv_EBE_zz_ln_mean, self.rfv_EBE_zz_ln_kappa)	
	name = "CMS_zz4l_massErrZZ_ld_sigma_{0:.0f}".format(self.channel)
	self.rfv_EBE_zz_ld_sigma = ROOT.RooFormulaVar(name, "("+self.inputs['relerr_qqzz_ld_mean']+")", ROOT.RooArgList(self.CMS_zz4l_mass))
	name = "CMS_zz4l_massErrZZ_ld_mean_{0:.0f}".format(self.channel)
	self.rfv_EBE_zz_ld_mean = ROOT.RooFormulaVar(name, "("+self.inputs['relerr_qqzz_ld_sigma']+")", ROOT.RooArgList(self.CMS_zz4l_mass)) 
	self.EBE_zz_ld = ROOT.RooLandau("errLD_qqzz","errLD_qqzz", self.RelErr, self.rfv_EBE_zz_ld_mean, self.rfv_EBE_zz_ld_sigma)
	name = "CMS_zz4l_massErrZZ_ld_frac_{0:.0f}".format(self.channel)
	self.rfv_EBE_zz_frac = ROOT.RooFormulaVar(name, "("+self.inputs['relerr_qqzz_ld_frac']+")", ROOT.RooArgList(self.MH))
	self.pdfErrZZ = ROOT.RooAddPdf("pdfErrZZ","pdfErrZZ", self.EBE_zz_ld, self.EBE_zz_ln, self.rfv_EBE_zz_frac)

	name = "CMS_zz4l_massErrZX_ln_kappa_{0:.0f}".format(self.channel);
	self.rfv_EBE_zx_ln_kappa = ROOT.RooFormulaVar(name, "("+self.inputs['relerr_zx_gs_sigma']+")", ROOT.RooArgList(self.CMS_zz4l_mass)); 
	name = "CMS_zz4l_massErrZX_ln_mean_{0:.0f}".format(self.channel);
	self.rfv_EBE_zx_ln_mean = ROOT.RooFormulaVar(name, "("+self.inputs['relerr_zx_gs_mean']+")", ROOT.RooArgList(self.CMS_zz4l_mass)); 
	self.EBE_zx_ln = ROOT.RooLognormal("errLN_zx","errLN_zx", self.RelErr, self.rfv_EBE_zx_ln_mean, self.rfv_EBE_zx_ln_kappa);	
	if self.channel!=1:
            self.EBE_zx_ln = ROOT.RooGaussian("errGaus_zx","errGaus_zx", self.RelErr, self.rfv_EBE_zx_ln_mean, self.rfv_EBE_zx_ln_kappa);
	name = "CMS_zz4l_massErrZX_ld_sigma_{0:.0f}".format(self.channel);
	self.rfv_EBE_zx_ld_sigma = ROOT.RooFormulaVar(name, "("+self.inputs['relerr_zx_ld_mean']+")", ROOT.RooArgList(self.CMS_zz4l_mass)); 
	name = "CMS_zz4l_massErrZX_ld_mean_{0:.0f}".format(self.channel);
	self.rfv_EBE_zx_ld_mean = ROOT.RooFormulaVar(name, "("+self.inputs['relerr_zx_ld_sigma']+")", ROOT.RooArgList(self.CMS_zz4l_mass)); 
	self.EBE_zx_ld = ROOT.RooLandau("errLD_zx","errLD_zx", self.RelErr, self.rfv_EBE_zx_ld_mean, self.rfv_EBE_zx_ld_sigma);	
	name = "CMS_zz4l_massErrZX_ld_frac_{0:.0f}".format(self.channel);
	self.rfv_EBE_zx_frac = ROOT.RooFormulaVar(name, "("+self.inputs['relerr_zx_ld_frac']+")", ROOT.RooArgList(self.MH)); 
	self.pdfErrZX = ROOT.RooAddPdf("pdfErrZX","pdfErrZX", self.EBE_zx_ld, self.EBE_zx_ln, self.rfv_EBE_zx_frac);


	self.sig_ggHErr = ROOT.RooProdPdf("sig_ggHErr","BW (X) CB * pdfErr", ROOT.RooArgSet(self.signalCB_ggH), ROOT.RooFit.Conditional(ROOT.RooArgSet(self.pdfErrS), ROOT.RooArgSet(self.RelErr)));
	self.sig_VBFErr = ROOT.RooProdPdf("sig_VBFErr","BW (X) CB * pdfErr", ROOT.RooArgSet(self.signalCB_VBF), ROOT.RooFit.Conditional(ROOT.RooArgSet(self.pdfErrS), ROOT.RooArgSet(self.RelErr)));
	self.sig_WHErr = ROOT.RooProdPdf("sig_WHErr","BW (X) CB * pdfErr", ROOT.RooArgSet(self.signalCB_WH), ROOT.RooFit.Conditional(ROOT.RooArgSet(self.pdfErrS), ROOT.RooArgSet(self.RelErr)));
	self.sig_ZHErr = ROOT.RooProdPdf("sig_ZHErr","BW (X) CB * pdfErr", ROOT.RooArgSet(self.signalCB_ZH), ROOT.RooFit.Conditional(ROOT.RooArgSet(self.pdfErrS), ROOT.RooArgSet(self.RelErr)));
	self.sig_ttHErr = ROOT.RooProdPdf("sig_ttHErr","BW (X) CB * pdfErr", ROOT.RooArgSet(self.signalCB_ttH), ROOT.RooFit.Conditional(ROOT.RooArgSet(self.pdfErrS), ROOT.RooArgSet(self.RelErr)));
	
        self.bkg_qqzzErr = ROOT.RooProdPdf("bkg_qqzzErr","bkg_qqzzErr", ROOT.RooArgSet(self.bkg_qqzz), ROOT.RooFit.Conditional(ROOT.RooArgSet(self.pdfErrZZ), ROOT.RooArgSet(self.RelErr)))
        self.bkg_ggzzErr = ROOT.RooProdPdf("bkg_ggzzErr","bkg_ggzzErr", ROOT.RooArgSet(self.bkg_ggzz), ROOT.RooFit.Conditional(ROOT.RooArgSet(self.pdfErrZZ), ROOT.RooArgSet(self.RelErr)))
        self.bkg_zjetsErr = ROOT.RooProdPdf("bkg_zjetsErr","bkg_zjetsErr", ROOT.RooArgSet(self.bkg_zjets), ROOT.RooFit.Conditional(ROOT.RooArgSet(self.pdfErrZX), ROOT.RooArgSet(self.RelErr)));
        





