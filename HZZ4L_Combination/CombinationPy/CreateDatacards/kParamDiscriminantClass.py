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
#from lib.util.Logger import *



class kParamDiscriminantClass(datacardClass):
  
    
    def __init__(self):
        #self.DEBUG=True
	#if self.DEBUG: level = 10
	#else: level = 20
	level = 10
	self.log = Logger().getLogger(self.__class__.__name__,level)
	self.DEBUG = (self.log.getEffectiveLevel()<20)
	print "++++++++++++++++++++++++++++++++++++++++++++++++++++ DEBUG = ", self.DEBUG
	
    def add_postfix_to_keys(self,mydict,postfix):
        return dict((k+postfix,f(v) if hasattr(v,'keys') else v) for k,v in mydict.items())
        
    def setTermNames(self, termNames=[]):
	self.allowedTermNames = ['ggH','gg0Ph','gg0M' ,'ggInt_12P','ggInt_12N','ggInt_13P','ggInt_13N','ggInt_23P','ggInt_23N']
	if len(termNames)>0 : 
	  self.termNames=[]
	  for term in termNames:
	      if term in self.allowedTermNames: self.termNames.append(term)
	      else :
		  raise RuntimeError,"Term {0} is not allowed for these datacards. Please check for typos but also in the allowedTermNames".format(term)
	  self.log.info("Setting termNames collection to user provided = {0}".format(self.termNames))
	  
	else :
	  self.termNames = ['ggH','gg0M','ggInt_13P','ggInt_13N']
	  self.log.info("Setting termNames collection to default = {0}".format(self.termNames))

    def TFile_safe_open(self, file_name, access='READ'):
        rootfile = ROOT.TFile.Open(file_name,access)
        if not rootfile:
            raise IOError, 'The file {0} either doesn\'t exist or cannot be open'.format(file_name)
        return rootfile
        
        
        
    def setSuperKD(self):
	self.isTemplate2D=True
	self.killBackground=True
	self.isRokoTest=True
	#self.DEBUG=False
	self.bkgMorph = False
	self.sigMorph = False
	if self.inputs['unfolded']: self.isTemplate2D=False
	
        self.pureStateDiscVarName = "CMS_zz4l_Djcp"
        print '>>>>>> Discriminant 1D PDFS using discriminat named :',self.pureStateDiscVarName
        
        templateSigName = "{0}/Dsignal_{1}.root".format(self.templateDir ,self.appendName)
        self.sigTempFile = self.TFile_safe_open(templateSigName)
        
        
        try:
	    self.termNames
	except:
	    self.setTermNames()
    	#self.termNames = ['ggH','gg0M']
	self.sigTemplate = {}
	self.sigTemplate_syst1Up = {}
	self.sigTemplate_syst1Down = {}
	self.sigTemplate_syst2Up = {}
	self.sigTemplate_syst2Down = {}
	for term in self.termNames:
	    self.sigTemplate[term] 		= self.sigTempFile.Get("{0}_shape".format(term))
	    #now we don't have systematic, so we provide the same templates
	    if self.isRokoTest:
		self.sigTemplate_syst1Up[term] 	= self.sigTempFile.Get("{0}_shape".format(term))
		self.sigTemplate_syst1Down[term] 	= self.sigTempFile.Get("{0}_shape".format(term))
		self.sigTemplate_syst2Up[term] 	= self.sigTempFile.Get("{0}_shape".format(term))
		self.sigTemplate_syst2Down[term] 	= self.sigTempFile.Get("{0}_shape".format(term))
	    else:
		self.sigTemplate_syst1Up[term] 	= self.sigTempFile.Get("{0}_shape_LeptScaleUp".format(term))
		self.sigTemplate_syst1Down[term] 	= self.sigTempFile.Get("{0}_shape_LeptScaleDown".format(term))
		self.sigTemplate_syst2Up[term] 	= self.sigTempFile.Get("{0}_shape_LeptSmearUp".format(term))
		self.sigTemplate_syst2Down[term] 	= self.sigTempFile.Get("{0}_shape_LeptSmearDown".format(term))
	    
	    
		
	self.log.info('Signal Templates File: {0}'.format(templateSigName))
        
        #### Set Proper Binning
        if self.DEBUG:
	    print "DEBUG: N bins in X = "   + str(self.sigTemplate['ggH'].GetXaxis().GetNbins())
	    print "DEBUG: Bin eges in X = " + str(self.sigTemplate['ggH'].GetXaxis().GetXbins().GetArray())
	   
        #self.xbins_sig  = ROOT.RooBinning(self.sigTemplate['ggH'].GetXaxis().GetNbins(),self.sigTemplate['ggH'].GetXaxis().GetXbins().GetArray())
	self.xbins_sig  = ROOT.RooBinning(self.sigTemplate['ggH'].GetXaxis().GetNbins(), self.sigTemplate['ggH'].GetXaxis().GetXmin(),self.sigTemplate['ggH'].GetXaxis().GetXmax())
					   
        if self.isTemplate2D: 
	    #self.ybins_sig  = ROOT.RooBinning(self.sigTemplate['ggH'].GetYaxis().GetNbins(), self.sigTemplate['ggH'].GetYaxis().GetXbins().GetArray())
	    self.ybins_sig  = ROOT.RooBinning(self.sigTemplate['ggH'].GetYaxis().GetNbins(), self.sigTemplate['ggH'].GetYaxis().GetXmin(),self.sigTemplate['ggH'].GetYaxis().GetXmax())
        ##Set Bins for discriminant D
        dBins = self.sigTemplate['ggH'].GetXaxis().GetNbins()
        dLow = self.sigTemplate['ggH'].GetXaxis().GetXmin()
        dHigh = self.sigTemplate['ggH'].GetXaxis().GetXmax()
        self.D0 = ROOT.RooRealVar(self.pureStateDiscVarName,self.pureStateDiscVarName,dLow,dHigh)
        self.D0.setBins(dBins)
        self.D0.setBinning(self.xbins_sig)
        self.log.info('pureStateDiscVarName: {0} and bins [low,high]: {1} [{2},{3}]'.format(self.pureStateDiscVarName, dBins,dLow,dHigh))
        if self.DEBUG:
            self.xbins_sig.Print()

        

        #################################
        if self.isTemplate2D : 
	    #todo: this has to be changed dependant on 2nd dimension in templates
	    self.interferenceDiscVarName = "CMS_zz4l_Djcp_int"
	    dBins = self.sigTemplate['ggH'].GetYaxis().GetNbins()
	    dLow = self.sigTemplate['ggH'].GetYaxis().GetXmin()
	    dHigh = self.sigTemplate['ggH'].GetYaxis().GetXmax()
	    self.D1 = ROOT.RooRealVar(self.interferenceDiscVarName,self.interferenceDiscVarName,dLow,dHigh)
	    self.D1.setBins(dBins)
	    self.D1.setBinning(self.ybins_sig)
	    
	    self.log.info('interferenceDiscVarName: {0} and bins [low,high]: {1} [{2},{3}]'.format(self.interferenceDiscVarName, dBins,dLow,dHigh))
	    if self.DEBUG:
                self.ybins_sig.Print()

        
    ## --------------------------- DATASET --------------------------- ##
    def fetchDatasetSuperKD(self):

        self.dataFileDir = "CMSdata"
        self.dataTreeName = "data_obs" 
        self.dataFileName = "{0}/hzz{1}_{2}.root".format(self.dataFileDir,self.appendName,self.lumi)
        
        #if self.isRokoTest : self.dataFileName = "{0}/hzz{1}_{2:.2f}.root".format("Sandbox/CMSdata/DSignal_Dmix90_4l_0M_QuasiData",self.appendName,self.lumi)
        #if self.isRokoTest : self.dataFileName = "{0}/hzz{1}_19.79.root".format("Sandbox/CMSdata/DSignal_Dmix90_4l_0M_QuasiData",self.appendName)
        if self.isRokoTest : 
            self.dataFileName = "Sandbox/data_obs.random.2D.root"
        
        if (self.DEBUG): 
            print self.dataFileName," ",self.dataTreeName 
        self.data_obs_file = self.TFile_safe_open(self.dataFileName)
        
        print self.data_obs_file.Get(self.dataTreeName)
        
        if not (self.data_obs_file.Get(self.dataTreeName)):
            self.log.error("File, \"{0}\", or tree, \"{1}\", not found".format(self.dataFileName,self.dataTreeName))
            print "Exiting..."
            sys.exit()
            
        self.data_obs_tree = self.data_obs_file.Get(self.dataTreeName)
        self.data_obs = ROOT.RooDataSet()
        self.datasetName = "data_obs_{0}".format(self.appendName)
                       
	self.data_obs = ROOT.RooDataSet(self.datasetName,self.datasetName,self.data_obs_tree,
					ROOT.RooArgSet(self.CMS_zz4l_mass,self.D0),'CMS_zz4l_mass>121.0&&CMS_zz4l_mass<131.0').reduce(ROOT.RooArgSet(self.D0))
	
	if self.isRokoTest:
	    self.data_obs = ROOT.RooDataSet(self.datasetName,self.datasetName,self.data_obs_tree,ROOT.RooArgSet(self.D0)).reduce(ROOT.RooArgSet(self.D0))
	  
	
        self.log.debug("Tree entries after reweighting= {0}: {1} : {2}".format(self.data_obs.isWeighted(), self.data_obs.sumEntries(), self.data_obs.numEntries()))
	
        if self.isTemplate2D: 
	    self.data_obs = ROOT.RooDataSet(self.datasetName,self.datasetName,self.data_obs_tree,ROOT.RooArgSet(self.CMS_zz4l_mass,self.D1,self.D0),'CMS_zz4l_mass>121.0&&CMS_zz4l_mass<131.0').reduce(ROOT.RooArgSet(self.D1,self.D0))
                

    def makeSuperKDAnalysis(self):

	self.ralTemplDimensions= ROOT.RooArgList(self.D0)
	self.rasTemplDimensions = ROOT.RooArgSet(self.D0)
	if self.isTemplate2D : 
	    self.ralTemplDimensions= ROOT.RooArgList(self.D0,self.D1)
            self.rasTemplDimensions = ROOT.RooArgSet(self.D0,self.D1)
            
	    
	#termNames = ['ggH','gg0M','ggInt_13P','ggInt_13N']
	
	self.sigTempDataHist = {}
	self.sigTempDataHist_syst1Up = {}
	self.sigTempDataHist_syst1Down = {}
	self.sigTempDataHist_syst2Up = {}
	self.sigTempDataHist_syst2Down = {}
	
	self.sigTemplatePdf = {}
	self.sigTemplatePdf_syst1Up = {}
	self.sigTemplatePdf_syst1Down = {}
	self.sigTemplatePdf_syst2Up = {}
	self.sigTemplatePdf_syst2Down = {}
	
	self.funcList = {}
	self.sigTemplateMorphPdf = {}
	self.sigCB2d = {}
	#create postfix for RooDataHist and RooHistPdf
	postfix={}
	for term in self.termNames:
	    postfix[term] = "{2}_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts,term)
	    
        self.log.debug("Postfixes for RooDataHist and RooHistPdf = {0}".format(postfix))
	
	for term in self.termNames:
	    self.log.debug("Making signal template for {0}".format(term))
	    #create RooDataHist
	    self.sigTempDataHist[term] 	   = ROOT.RooDataHist("sigTempDataHist_{0}".format(postfix[term]),"sigTempDataHist_{0}".format(postfix[term]),self.ralTemplDimensions,ROOT.RooFit.Import(self.sigTemplate[term],kFALSE))
	    self.sigTempDataHist_syst1Up[term]   = ROOT.RooDataHist("sigTempDataHist_syst1Up_{0}".format(postfix[term]),"sigTempDataHist_syst1Up_{0}".format(postfix[term]),self.ralTemplDimensions,ROOT.RooFit.Import(self.sigTemplate_syst1Up[term],kFALSE))
	    self.sigTempDataHist_syst1Down[term] = ROOT.RooDataHist("sigTempDataHist_syst1Down_{0}".format(postfix[term]),"sigTempDataHist_syst1Down_{0}".format(postfix[term]),self.ralTemplDimensions,ROOT.RooFit.Import(self.sigTemplate_syst1Down[term],kFALSE))
	    self.sigTempDataHist_syst2Up[term]   = ROOT.RooDataHist("sigTempDataHist_syst2Up_{0}".format(postfix[term]),"sigTempDataHist_syst2Up_{0}".format(postfix[term]),self.ralTemplDimensions,ROOT.RooFit.Import(self.sigTemplate_syst1Up[term],kFALSE))
	    self.sigTempDataHist_syst2Down[term] = ROOT.RooDataHist("sigTempDataHist_syst2Down_{0}".format(postfix[term]),"sigTempDataHist_syst2Down_{0}".format(postfix[term]),self.ralTemplDimensions,ROOT.RooFit.Import(self.sigTemplate_syst1Down[term],kFALSE))

            #self.sigTempDataHist[term]     = ROOT.RooDataHist("sigTempDataHist_{0}".format(postfix[term]),"sigTempDataHist_{0}".format(postfix[term]),self.ralTemplDimensions,ROOT.RooFit.Import(self.sigTemplate[term],kTRUE))
            #self.sigTempDataHist_syst1Up[term]   = ROOT.RooDataHist("sigTempDataHist_syst1Up_{0}".format(postfix[term]),"sigTempDataHist_syst1Up_{0}".format(postfix[term]),self.ralTemplDimensions,ROOT.RooFit.Import(self.sigTemplate_syst1Up[term],kTRUE))
            #self.sigTempDataHist_syst1Down[term] = ROOT.RooDataHist("sigTempDataHist_syst1Down_{0}".format(postfix[term]),"sigTempDataHist_syst1Down_{0}".format(postfix[term]),self.ralTemplDimensions,ROOT.RooFit.Import(self.sigTemplate_syst1Down[term],kTRUE))
            #self.sigTempDataHist_syst2Up[term]   = ROOT.RooDataHist("sigTempDataHist_syst2Up_{0}".format(postfix[term]),"sigTempDataHist_syst2Up_{0}".format(postfix[term]),self.ralTemplDimensions,ROOT.RooFit.Import(self.sigTemplate_syst1Up[term],kTRUE))
            #self.sigTempDataHist_syst2Down[term] = ROOT.RooDataHist("sigTempDataHist_syst2Down_{0}".format(postfix[term]),"sigTempDataHist_syst2Down_{0}".format(postfix[term]),self.ralTemplDimensions,ROOT.RooFit.Import(self.sigTemplate_syst1Down[term],kTRUE))
            
	    
	    
	    
	    #create RooDataHistPdf
	    self.sigTemplatePdf[term] 		 = ROOT.RooHistPdf("sigTemplatePdf_{0}".format(postfix[term]),"sigTemplatePdf_{0}".format(postfix[term]),self.rasTemplDimensions,self.sigTempDataHist[term])
	    self.sigTemplatePdf_syst1Up[term]   = ROOT.RooHistPdf("sigTemplatePdf_syst1Up_{0}".format(postfix[term]),"sigTemplatePdf_syst1Up_{0}".format(postfix[term]),self.rasTemplDimensions,self.sigTempDataHist_syst1Up[term])
	    self.sigTemplatePdf_syst1Down[term] = ROOT.RooHistPdf("sigTemplatePdf_syst1Down_{0}".format(postfix[term]),"sigTemplatePdf_syst1Down_{0}".format(postfix[term]),self.rasTemplDimensions,self.sigTempDataHist_syst1Down[term])
	    self.sigTemplatePdf_syst2Up[term]   = ROOT.RooHistPdf("sigTemplatePdf_syst2Up_{0}".format(postfix[term]),"sigTemplatePdf_syst2Up_{0}".format(postfix[term]),self.rasTemplDimensions,self.sigTempDataHist_syst1Up[term])
	    self.sigTemplatePdf_syst2Down[term] = ROOT.RooHistPdf("sigTemplatePdf_syst2Down_{0}".format(postfix[term]),"sigTemplatePdf_syst2Down_{0}".format(postfix[term]),self.rasTemplDimensions,self.sigTempDataHist_syst1Down[term])
        
	    ###Shape systematics for signal
	    
	    self.funcList[term] = ROOT.RooArgList()

	    if self.sigMorph:
		self.funcList[term].add(self.sigTemplatePdf[term])
		self.funcList[term].add(self.sigTemplatePdf_syst1Up[term])
		self.funcList[term].add(self.sigTemplatePdf_syst1Down[term])
		self.funcList[term].add(self.sigTemplatePdf_syst2Up[term])
		self.funcList[term].add(self.sigTemplatePdf_syst2Down[term])  
	    else:
		self.funcList[term].add(self.sigTemplatePdf[term])
	    


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

	    if self.isTemplate2D: 
		#self.sigTemplateMorphPdf[term] = ROOT.FastVerticalInterpHistPdf2D("sigTemplateMorphPdf_{0}".format(postfix[term]),"sigTemplateMorphPdf_{0}".format(postfix[term]),self.D1,self.D0,False,self.funcList[term],self.morphVarListSig1,1.0,1)
		self.sigTemplateMorphPdf[term] = ROOT.FastVerticalInterpHistPdf2D("sigTemplateMorphPdf_{0}".format(postfix[term]),"sigTemplateMorphPdf_{0}".format(postfix[term]),self.D0,self.D1,False,self.funcList[term],self.morphVarListSig1,1.0,1)
	    else :
		self.sigTemplateMorphPdf[term] = ROOT.FastVerticalInterpHistPdf("sigTemplateMorphPdf_{0}".format(postfix[term]),"sigTemplateMorphPdf_{0}".format(postfix[term]),self.D0,self.funcList[term],self.morphVarListSig1,1.0,1)
		
	    self.sigCB2d[term] = self.signalCB_ggH.Clone("sigCB2d_{0}".format(term))

	    #end "for term in self.termNames:"
	
	    
        ## ----------------- SuperKD 2D BACKGROUND SHAPES --------------- ##
        templateBkgName = "{0}/Dbackground_qqZZ_{1}.root".format(self.templateDir ,self.appendName)
        self.bkgTempFile = self.TFile_safe_open(templateBkgName)
        self.bkgTemplate = self.bkgTempFile.Get("qqZZ_shape")  #mod-roko
        self.bkgTemplateqqZZ = self.bkgTempFile.Get("qqZZ_shape")
        
        #if self.isRokoTest:
	    #self.bkgTemplate = self.bkgTempFile.Get("ggH_shape")  #mod-roko
	    #self.bkgTemplateqqZZ = self.bkgTempFile.Get("ggH_shape")
        

        TemplateName = "bkgTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.bkgTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,self.ralTemplDimensions,ROOT.RooFit.Import(self.bkgTemplate,kFALSE))
        TemplateName = "bkgTemplatePdf_qqzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.bkgTemplatePdf_qqzz = ROOT.RooHistPdf(TemplateName,TemplateName,self.rasTemplDimensions,self.bkgTempDataHist)


        if self.reflectZX:
            
            TemplateName = "zjetsbkgTempDataHist_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.bkgTempDataHist_zjetsUp = ROOT.RooDataHist(TemplateName,TemplateName,self.ralTemplDimensions,ROOT.RooFit.Import(self.bkgTemplate,kFALSE))
            TemplateName = "bkgTemplatePdf_zjets_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.bkgTemplatePdf_zjets_Up = ROOT.RooHistPdf(TemplateName,TemplateName,self.rasTemplDimensions,self.bkgTempDataHist_zjetsUp)
        

            templateBkgName = "{0}/Dbackground_ZJetsCR_AllChans.root".format(self.templateDir)
            self.zjetsTempFile = self.TFile_safe_open(templateBkgName)
            
            self.zjetsTemplate = self.zjetsTempFile.Get("zjets_shape")
            #if not self.isRokoTest : self.zjetsTemplate = self.zjetsTempFile.Get("zjets_shape")
            #else : self.zjetsTemplate = self.zjetsTempFile.Get("ggH_shape")
            
            TemplateName = "zjetsbkgTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.bkgTempDataHist_zjets = ROOT.RooDataHist(TemplateName,TemplateName,self.ralTemplDimensions,ROOT.RooFit.Import(self.zjetsTemplate,kFALSE))
            TemplateName = "bkgTemplatePdf_zjets_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.bkgTemplatePdf_zjets = ROOT.RooHistPdf(TemplateName,TemplateName,self.rasTemplDimensions,self.bkgTempDataHist_zjets)

            self.bkgTempDataHist_zjetsDown1 = self.reflectSystematics(self.zjetsTemplate,self.bkgTemplateqqZZ)
            TemplateName = "zjetsbkgTempDataHist_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.bkgTempDataHist_zjetsDown = ROOT.RooDataHist(TemplateName,TemplateName,self.ralTemplDimensions,ROOT.RooFit.Import(self.bkgTempDataHist_zjetsDown1,kFALSE))
            TemplateName = "bkgTemplatePdf_zjets_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.bkgTemplatePdf_zjets_Down = ROOT.RooHistPdf(TemplateName,TemplateName,self.rasTemplDimensions,self.bkgTempDataHist_zjetsDown)
        

        else:

            #Nominal Z+X
            TemplateName = "zjetsbkgTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.bkgTempDataHist_zjets = ROOT.RooDataHist(TemplateName,TemplateName,self.ralTemplDimensions,ROOT.RooFit.Import(self.bkgTemplate,kFALSE))
            TemplateName = "bkgTemplatePdf_zjets_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.bkgTemplatePdf_zjets = ROOT.RooHistPdf(TemplateName,TemplateName,self.rasTemplDimensions,self.bkgTempDataHist_zjets)
            

            #Up Z+X
            TemplateName = "zjetsbkgTempDataHist_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.bkgTempDataHist_zjetsUp = ROOT.RooDataHist(TemplateName,TemplateName,self.ralTemplDimensions,ROOT.RooFit.Import(self.bkgTemplate,kFALSE))
            TemplateName = "bkgTemplatePdf_zjets_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.bkgTemplatePdf_zjets_Up = ROOT.RooHistPdf(TemplateName,TemplateName,self.rasTemplDimensions,self.bkgTempDataHist_zjetsUp)
        

            #Down Z+X
            templateBkgName = "{0}/Dbackground_ZJetsCR_AllChans.root".format(self.templateDir)
            self.zjetsTempFile = self.TFile_safe_open(templateBkgName)
            
            #self.zjetsTemplateDown = self.zjetsTempFile.Get("h_superDpsD")
            self.zjetsTemplateDown = self.zjetsTempFile.Get("zjets_shape")
            #if not self.isRokoTest : self.zjetsTemplateDown = self.zjetsTempFile.Get("zjets_shape")
            #else : self.zjetsTemplateDown = self.zjetsTempFile.Get("ggH_shape")
            
            TemplateName = "zjetsbkgTempDataHist_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.bkgTempDataHist_zjetsDown = ROOT.RooDataHist(TemplateName,TemplateName,self.ralTemplDimensions,ROOT.RooFit.Import(self.zjetsTemplateDown,kFALSE))
            TemplateName = "bkgTemplatePdf_zjets_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
            self.bkgTemplatePdf_zjets_Down = ROOT.RooHistPdf(TemplateName,TemplateName,self.rasTemplDimensions,self.bkgTempDataHist_zjetsDown)




        templateggBkgName = "{0}/Dbackground_ggZZ_{1}.root".format(self.templateDir ,self.appendName)
        self.ggbkgTempFile = self.TFile_safe_open(templateggBkgName)
        #self.ggbkgTemplate = self.ggbkgTempFile.Get("h_superDpsD")
        self.ggbkgTemplate = self.ggbkgTempFile.Get("ggZZ_shape")
	#if not self.isRokoTest : self.ggbkgTemplate = self.ggbkgTempFile.Get("ggZZ_shape")
	#else : self.ggbkgTemplate = self.ggbkgTempFile.Get("ggH_shape")
    
        
        TemplateName = "ggbkgTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.ggbkgTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,self.ralTemplDimensions,ROOT.RooFit.Import(self.ggbkgTemplate,kFALSE))
        TemplateName = "bkgTemplatePdf_ggzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        self.bkgTemplatePdf_ggzz = ROOT.RooHistPdf(TemplateName,TemplateName,self.rasTemplDimensions,self.ggbkgTempDataHist)


        #################################
        
        self.funcList_zjets = ROOT.RooArgList()
        morphBkgVarName =  "CMS_zz4l_smd_zjets_bkg_{0:.0f}".format(self.channel)
        self.alphaMorphBkg = ROOT.RooRealVar(morphBkgVarName,morphBkgVarName,0,-20,20)
        self.morphVarListBkg = ROOT.RooArgList()
        
        if self.bkgMorph:
            self.funcList_zjets.add(self.bkgTemplatePdf_zjets)
            self.funcList_zjets.add(self.bkgTemplatePdf_zjets_Up)
            self.funcList_zjets.add(self.bkgTemplatePdf_zjets_Down)
            self.alphaMorphBkg.setConstant(False) 
            self.morphVarListBkg.add(self.alphaMorphBkg)
        else:
            self.funcList_zjets.add(self.bkgTemplatePdf_zjets)
            self.alphaMorphBkg.setConstant(True) 
            
	if self.isTemplate2D: 
	    TemplateName = "bkgTemplateMorphPdf_qqzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
	    self.bkgTemplateMorphPdf_qqzz = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,self.D1,self.D0,False,ROOT.RooArgList(self.bkgTemplatePdf_qqzz),ROOT.RooArgList(),1.0,1)
	    TemplateName = "bkgTemplateMorphPdf_ggzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
	    self.bkgTemplateMorphPdf_ggzz = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,self.D1,self.D0,False,ROOT.RooArgList(self.bkgTemplatePdf_ggzz),ROOT.RooArgList(),1.0,1)
	    
	    TemplateName = "bkgTemplateMorphPdf_zjets_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
	    self.bkgTemplateMorphPdf_zjets = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,self.D1,self.D0,False,self.funcList_zjets,self.morphVarListBkg,1.0,1)
        else :
	    TemplateName = "bkgTemplateMorphPdf_qqzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
	    self.bkgTemplateMorphPdf_qqzz = ROOT.FastVerticalInterpHistPdf(TemplateName,TemplateName,self.D0,ROOT.RooArgList(self.bkgTemplatePdf_qqzz),ROOT.RooArgList(),1.0,1)
	    TemplateName = "bkgTemplateMorphPdf_ggzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
	    self.bkgTemplateMorphPdf_ggzz = ROOT.FastVerticalInterpHistPdf(TemplateName,TemplateName,self.D0,ROOT.RooArgList(self.bkgTemplatePdf_ggzz),ROOT.RooArgList(),1.0,1)
	    
	    TemplateName = "bkgTemplateMorphPdf_zjets_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
	    self.bkgTemplateMorphPdf_zjets = ROOT.FastVerticalInterpHistPdf(TemplateName,TemplateName,self.D0,self.funcList_zjets,self.morphVarListBkg,1.0,1)
	    
        self.bkg1d_qqzz  = self.bkg_qqzz.Clone("bkg1d_qqzz")
        self.bkg1d_ggzz  = self.bkg_ggzz.Clone("bkg1d_ggzz")
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
        try:
	    self.data_obs
	except:
	    self.log.warn("No dataset \'data_obs\' will be imported to workspace.")
	else:    
	    getattr(self.w,'import')(self.data_obs,ROOT.RooFit.Rename("data_obs"))
        
        
        #termNames = ['ggH','gg0M','ggInt_13P','ggInt_13N']
           
        for term in self.termNames:
	    #self.sigCB2d[term].SetNameTitle("ORIG{0}".format(term),"ORIG{0}".format(term))
	    #getattr(self.w,'import')(self.sigCB2d[term], ROOT.RooFit.RecycleConflictNodes())
	    
	    self.sigTemplateMorphPdf[term].SetNameTitle(term,term)  #maybe they don't have SetNameTitle Mmethod
	    getattr(self.w,'import')(ROOT.RooArgSet(self.sigTemplateMorphPdf[term]), ROOT.RooFit.RecycleConflictNodes())
	    
	    #getattr(self.w,'import')(self.sigTemplatePdf[term], ROOT.RooFit.RecycleConflictNodes())

	self.log.debug("Ended importing signal pdfs for terms: {0}".format(self.termNames))
	
        #save syst templates individually
        #todo: 
        systTempName=("ggHCMS_zz4l_leptScale_sig_{0}_{1:.0f}_Up").format(self.channel,self.sqrts)
        self.sigTemplatePdf_syst1Up['ggH'].SetNameTitle(systTempName,systTempName)
        systTempName=("ggHCMS_zz4l_leptScale_sig_{0}_{1:.0f}_Down").format(self.channel,self.sqrts)
        self.sigTemplatePdf_syst1Down['ggH'].SetNameTitle(systTempName,systTempName)
        
        self.bkg1d_qqzz.SetNameTitle("ORIGbkg_qqzz","ORIGbkg_qqzz")
        self.bkg1d_ggzz.SetNameTitle("ORIGbkg_ggzz","ORIGbkg_ggzz")
        self.bkg1d_zjets.SetNameTitle("ORIGbkg_zjets","ORIGbkg_zjets")
        
        self.bkgTemplateMorphPdf_qqzz.SetNameTitle("bkg2d_qqzz","bkg2d_qqzz")
        self.bkgTemplateMorphPdf_ggzz.SetNameTitle("bkg2d_ggzz","bkg2d_ggzz")
        self.bkgTemplateMorphPdf_zjets.SetNameTitle("bkg2d_zjets","bkg2d_zjets")
        
        getattr(self.w,'import')(self.bkgTemplateMorphPdf_qqzz ,ROOT.RooFit.RecycleConflictNodes())
        getattr(self.w,'import')(self.bkgTemplateMorphPdf_ggzz ,ROOT.RooFit.RecycleConflictNodes())
        getattr(self.w,'import')(self.bkgTemplateMorphPdf_zjets,ROOT.RooFit.RecycleConflictNodes())
        
        ### save signal rates
        
        #mod-roko
	self.k3k1_ratio= ROOT.RooRealVar("k3k1_ratio","k3k1_ratio",0,-100,100)
	self.k3k1_ratio.setConstant(True)

	self.k2k1_ratio= ROOT.RooRealVar("k2k1_ratio","k2k1_ratio",0,-100,100)
	self.k2k1_ratio.setConstant(True)
	
	
	
	self.rfvSigRate = {}
	
	postfix = "{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
	self.rfvSigRate_all = ROOT.RooFormulaVar("all_norm_{0}".format(postfix),"@0+@1+@2+@3+@4",ROOT.RooArgList(self.rfvSigRate_ggH_temp,self.rfvSigRate_VBF,self.rfvSigRate_WH,self.rfvSigRate_ZH,self.rfvSigRate_ttH))
	#set to 1 since the rate will be in the datacard explicitly
	#self.rfvSigRate_all = ROOT.RooFormulaVar("all_norm","@0",ROOT.RooArgList(self.one))
	

	
	if self.DEBUG : self.rfvSigRate_all.Print()
	rate_input = {
		      #'gg0M'  : 'jhuGen_0minus_yield',
		      #'gg0Ph' :  'jhuGen_0hplus_yield'
		      }	   
		      
	#self.rfvGeoLocNorms = self._setupGeolocationNormalization(self.termNames)
	self._setupGeolocationNormalization(self.termNames)
        
	for term in self.termNames:
            self.log.debug('Building and importing the final norms into worksppace for term: {0}'.format(term))
	    #self.rfvSigRate[term] = ROOT.RooFormulaVar("{0}_norm".format(term),"@0*@1",ROOT.RooArgList(self.rfvSigRate_all,self.rfvGeoLocNorms[term] ))  #mod-roko -- multiply by gama lambda dependant factor
	    #to avoid putting all the SM norms per production channel into WS, we just put the final value as number. 
	    self.rfvGeoLocNorm_tmp = self.w.function("geoloc_{0}_norm_{1}".format(term, postfix))
	    self.rfvSigRate[term] = ROOT.RooFormulaVar("{0}_norm".format(term),"{0}/{1}*@0*(@1)".format(self.rfvSigRate_all.getVal(), self.rrvLumi.getVal()),
                                                                                                      ROOT.RooArgList(self.rrvLumi,self.rfvGeoLocNorm_tmp ))  
	    if self.DEBUG : 
                    self.rfvSigRate[term].Print()
	    getattr(self.w,'import')(self.rfvSigRate[term], ROOT.RooFit.RecycleConflictNodes())
	    
	    #fix for JHUGen rates
	    self.rates[term] = 1
	    	    ####################################
 
	    #if term in rate_input.keys():
		#self.rates[term] = float(self.inputs[rate_input[term]]) * self.rrv_SMggH_ratio.getVal()
	    #else:
		#self.rates[term] = self.rrvJHUgen_SMggH.getVal()*self.rrv_SMggH_ratio.getVal() #the same as for ggH 
		
	self.bkg2d_ggzz_norm  =  ROOT.RooFormulaVar("bkg2d_ggzz_norm","@0/{0}".format(self.rrvLumi.getVal()),ROOT.RooArgList(self.rrvLumi))  
	self.bkg2d_qqzz_norm  =  ROOT.RooFormulaVar("bkg2d_qqzz_norm","@0/{0}".format(self.rrvLumi.getVal()),ROOT.RooArgList(self.rrvLumi))  
	self.bkg2d_zjets_norm =  ROOT.RooFormulaVar("bkg2d_zjets_norm","@0/{0}".format(self.rrvLumi.getVal()),ROOT.RooArgList(self.rrvLumi)) 
	self.log.debug('Importing the background scaling with luminosity...')
	getattr(self.w,'import')(self.bkg2d_ggzz_norm, ROOT.RooFit.RecycleConflictNodes())
	getattr(self.w,'import')(self.bkg2d_qqzz_norm, ROOT.RooFit.RecycleConflictNodes())
	getattr(self.w,'import')(self.bkg2d_zjets_norm, ROOT.RooFit.RecycleConflictNodes())
	
  
	if self.DEBUG: 
	    print 20*"----" 
	    self.w.Print()
	    print 20*"----" 
        
        self.w.writeToFile(self.name_ShapeWS)
        
                
    
    def prepareDatacardSuperKD(self):

        name_Shape = ""
        name_ShapeWS = ""

        ## Write Datacards
        if (self.endsInP5(self.mH)):
            name_Shape = "{0}/HCG/{1:.1f}/hzz4l_{2}S_{3:.0f}TeV{4}.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts,self.appendHypType)
        else:
            name_Shape = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV{4}.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts,self.appendHypType)
        with open( name_Shape, "wb") as fo:
            self.log.info('Writing textual datacards...')
            self.WriteDatacardSuperKD(fo, self.name_ShapeWS2, self.rates, self.data_obs.numEntries())
            #self.systematics.WriteSystematics(fo, self.inputs)
            #self.systematics.WriteShapeSystematics(fo,self.inputs)
            #self.systematics.WriteSuperKDShapeSystematics(fo,self.inputs)
        if self.DEBUG:    
            #read the datacard again and dump its contents 
            self.log.debug('Dumping the contents of datacard {0}'.format(name_Shape))
            with open( name_Shape, "r") as fo:
                print fo.read()


        
    def WriteDatacardSuperKD(self,file,nameWS,theRates,obsEvents):

        numberSig = self.numberOfSigChan(self.inputs)
        numberSig = len(self.termNames)  #mod-roko --> this is temporary
        numberBg  = self.numberOfBgChan(self.inputs)
        numberBg  = 3 #mod-roko --> this is temporary before changing the sumation function
        numberBg  = 1 #mod-roko --> this is temporary before changing the sumation function
        

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

	#termNames = ['ggH','gg0M','ggInt_13P','ggInt_13N']

        channelList=['ggH','qqH','WH','ZH','ttH','qqZZ','ggZZ','zjets','ttbar','zbb']
        channelName2D=['ggH','qqH','WH','ZH','ttH','bkg2d_qqzz','bkg2d_ggzz','bkg2d_zjets','bkg2d_ttbar','bkg2d_zbb']

        if self.inputs["all"]:
            print "ALL CHANNELS --- KDCLASS --- WriteDatacardKD\n"
            #channelList=['ggH','ggH','qqZZ','ggZZ','zjets']
            #channelName2D=['ggH','ggH{0}'.format(self.appendHypType),'bkg2d_qqzz','bkg2d_ggzz','bkg2d_zjets']
	    #channelList=['ggH','ggH','ggH','ggH','qqZZ','ggZZ','zjets','ttbar','zbb']
	    #channelName2D=['ggH','gg0M','ggInt_13P','ggInt_13N','bkg2d_qqzz','bkg2d_ggzz','bkg2d_zjets','bkg2d_ttbar','bkg2d_zbb']
            #channelList=['ggH','ggH','ggH','ggH','qqZZ','ggZZ','zjets']
	    #channelName2D=['ggH','gg0M','ggInt_13P','ggInt_13N','bkg2d_qqzz','bkg2d_ggzz','bkg2d_zjets']
	    channelList=[]
	    channelName2D=[]
	    for term in self.termNames:
		channelList.append('ggH')
		channelName2D.append(term)
	    
	    #channelList+=['qqZZ','ggZZ','zjets']
	    #channelName2D+=['bkg2d_qqzz','bkg2d_ggzz','bkg2d_zjets']
            channelList+=['qqZZ']
	    channelName2D+=['bkg2d_qqzz']
	     
          
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
        scaleFactor = {}
        #scaleFactor['8TeVTo14TevXsection_sig'] = 2.6
        scaleFactor['8TeVTo14TevXsection_sig'] = 2.54
        scaleFactor['8TeVTo14TevXsection_bkg'] = 1.9
        #scaleFactor['ExpectedToObserved_sig'] = 18.7/22.50768 #pedjas expected / my expected @
        scaleFactor['ExpectedToObserved_sig'] = 18.7/20.694727031574804
        #scaleFactor['ExpectedToObserved_sig'] = 1 #my expected @
        scaleFactor['ExpectedToObserved_bkg'] = 1 #my expected @
        
        
        scale_sig=1
        scale_bkg=1
        try:
	    self.inputs['user_option']
	except:
	    print "No user option defined..."
	    scale_sig=1
	    scale_bkg=1
	else: 
	    if "14TeV_rescale" in self.inputs['user_option'] :
		scale_sig  = scaleFactor['8TeVTo14TevXsection_sig']*scaleFactor['ExpectedToObserved_sig']
		scale_bkg  = scaleFactor['8TeVTo14TevXsection_bkg']*scaleFactor['ExpectedToObserved_bkg']
	    else:
		scale_sig=1
		scale_bkg=1
	
	    if "UniformBkg" in self.inputs['user_option'] : 
                self.killBackground = False
	    if "NoDenominator"in self.inputs['user_option'] : 
                self.removeGeolocDenominator = True
        
        for chan in channelList:
            if self.inputs[chan] or self.inputs["all"]:
                #if channelName2D[j].startswith('ggH{0}'.format(self.appendHypType)):
                    #theYield = theRates['ggH{0}'.format(self.appendHypType)]
                    #file.write("{0:.4f} ".format(theYield))
                if channelName2D[j] in self.termNames:
		    file.write("{0:.4f} ".format(theRates[channelName2D[j]]*scale_sig))
		    
		elif self.killBackground:
		    file.write("1e-5 ")
                else:
		    #file.write("{0:.4f} ".format(theRates[chan]*scale_bkg))
		    file.write("{0:.4f} ".format((theRates['qqZZ']+theRates['ggZZ']+theRates['zjets'])*scale_bkg)) #the rate of qqZZ+ggZZ+zjets goes under qqZZ
            j+=1
            
        file.write("\n")
        file.write("------------\n")
        
        
        #test systematics mod-roko
        file.write("\n")
            
        processLine = "test_sys lnN "

        for x in range(-numberSig+1,1):
            processLine += "1.00001 "

        for y in range(1,numberBg+1):
            processLine += "1.00001 "

        file.write(processLine)
        file.write("\n")
        
        
        #Channel list will be asked in systematicsClass.py
        #currently we are putting all the signal terms to 
        #have the same systematic as ggH
        self.inputs['terms_for_systematics'] = channelList  
        
        
            
  
    def _setupGeolocationNormalization(self, termNames=None):
	"""
	Normalization in front of terms in parametrization a la Geolocating paper.
	It depends on k1,k3,lambda, gamma and channel(thru acceptance)
	"""
	if termNames==None:
            termNames=[]
	
	def _readFactorsFromFile(file_name, termNames, channel = self.channel):
		    
	    if len(termNames)==0 or not isinstance(termNames, list): 
		raise RuntimeError, "_setupGeolocationNormalization:: You have to provide a list of terms for which I can calculate the normalization."
	    
	    sigTempFile = self.TFile_safe_open(file_name)
	    factorNames= ['gamma11','gamma22','gamma33','gamma12','gamma13','gamma23',
		      'lambda13_cosP','lambda13_cosN','lambda13_sinP','lambda13_sinN',
		      'lambda12_cosP','lambda12_cosN','lambda12_sinP','lambda12_sinN',
		      'lambda23_cosP','lambda23_cosN','lambda23_sinP','lambda23_sinN']
		      
	    needed_factors = {
	      'ggH'	  :   'gamma11,gamma12',
	      'gg0Ph'    :   'gamma22,gamma12',     # @0 = k2k1_ratio
	      'gg0M'     :   'gamma33,gamma12',     # @1 = k3k1_ratio
	      'ggInt_12P':   'lambda12_cosP',  
	      'ggInt_12N':   'lambda12_cosN',  
	      
	      'ggInt_13P':   'lambda13_cosP',  
	      'ggInt_13N':   'lambda13_cosN',  
	      'ggInt_23P':   'lambda23_cosP',  
	      'ggInt_23N':   'lambda23_cosN'
	    }		    
	    factors_to_read = []
	    for term in termNames: 
		if term in needed_factors.keys(): 
                    #factors_to_read.append(needed_factors[term])
                    this_term_factors_to_read = [factor.strip() for factor in needed_factors[term].split(',')]
                    factors_to_read+=this_term_factors_to_read
		else:
		    raise KeyError, "_setupGeolocationNormalization:: The term \"{0}\" is unknown. Please check spelling. Allowed terms are : {1}".format(term, str(needed_factors.keys()))
            factors_to_read = list(set(factors_to_read))  #make list with unique entries
            self.log.debug('Factors to read from file {0}: {1}'.format(file_name, str(factors_to_read)))
	    
	    #initialize factor dictionary	
	    factors={}
	    sigTempFile.cd()
	    tree=sigTempFile.Get('factors')
            
            tree.GetEntry(0)
            for fn in factorNames:
                if fn in factors_to_read:
                    
                    try:
                        #factors[fn] = eval('tree.{0}'.format(fn))
                        factor_value = eval('tree.{0}'.format(fn))
                        factors[fn] = float('{0:.4f}'.format(factor_value))  #jsut to have a factor value up to 3 decimals precirsion
                    except:
                        self.log.warn('The factors tree does not contain the {0} branch. Please check the tree or the spelling of the variable! Currently setting value to 0.'.format(fn))
                        factors[fn] = 0
                else:
                    factors[fn]=0
	   
	   
	    #if channel==3:
	      ##gamma 33 = 0.040 #2e2mu
	      #factors['gamma33']=0.040
	    #else:
	      #factors['gamma33']=0.034
	      
	    #providional until 
	    #self.log.warn('Current values used for interference 12 are provisional. Change them when you get templates from Pedja.')
	    #factors['gamma11']=1  
            #factors['gamma22']=0.090
            #if channel==3:
            #factors['lambda12_cosN']*=0.9
            #factors['gamma22']*=1.1
            #factors['gamma12'] = -factors['lambda12_cosN']
            #factors['lambda12_cosN'] = -factors['gamma12']
            #factors['lambda12_cosP'] = factors['lambda12_cosN'] = 0.0280214
	
	
            #change name and make RooRealVar
            appendName = "_{0}_{1:.0f}".format(channel, self.sqrts)
            try:
                mean_lamda_13 = (factors['lambda13_cosN']+factors['lambda13_cosP'])/2.0
            except KeyError:
                pass
            else:
                factors['lambda13_cosP'] = mean_lamda_13
                factors['lambda13_cosN'] = mean_lamda_13
            for fn in factors.keys():
                factor_name = '{0}{1}'.format(fn,appendName)
                factors[fn] = ROOT.RooRealVar(factor_name, factor_name, factors[fn])
                #factors[fn].setConstant()
           
            factors = self.add_postfix_to_keys(factors,appendName)
            #self.log.debug('Factors after adding postfix: {0}'.format(factors))
            return factors
	    
		    
	def _getFinalStateYieldRatios():
	    #pass
	    #f.s.==self.appendName
	    #r = N(f.s.)/N_TOT
	    #normFile = self.TFile_safe_open("ws_norm.root")
	    #all_norm = {}
	    #tot_norm = 0
	    #self.log.debug('getFinalStateYieldRatio -> getting total normalization for the 3 final states')
	    #for fs in [1,2,3]:
		#appendName = "{0}_{1:.0f}".format(fs, self.sqrts)
		#print "NAME WS = ws_norm_{0}".format(appendName)
		#w = normFile.Get("ws_norm_{0}".format(appendName))
		#w.var('cmshzz4l_lumi').setVal(self.lumi)
		##if self.DEBUG: w.Print()
		#all_norm[fs] = w.function("all_norm_{0}".format(appendName)).getVal()
		#tot_norm += all_norm[fs]
		##self.log.debug('getFinalStateYieldRatio {0} = {1}'.format(appendName, str(all_norm[fs])) )
		#self.log.debug('getFinalStateYieldRatio -> all_norm_{0} = {1}'.format(appendName, str(all_norm[fs])) )
		
	    #self.log.debug('getFinalStateYieldRatio -> tot_norm_{0} = {1}'.format(appendName, tot_norm) )	
	    #self.log.debug('getFinalStateYieldRatio -> channel : tot_norm_{0} = {1}'.format(appendName, tot_norm) )	
	    #assert tot_norm > 0, "Total norm has to be grater than 0 to do the ratio."
	    #yield_ratios = {channel : all_norm[channel]/tot_norm  for channel in [1,2,3]}  
	    #self.log.debug('getFinalStateYieldRatio -> yield ratios = {0}'.format(yield_ratios) )	
	    yield_ratios = {1: 0.3543010974318054, 2: 0.18414932432211456, 3: 0.46154957824608}
	    #yield_ratios = {1: 0.354, 2: 0.184, 3: 0.462}
	    return yield_ratios
	    
	    
	def _getDenominatorSegment(channel,den_base):
	    import os
	    appendName = "{0}_{1:.0f}".format(channel, self.sqrts)
	    #yield_ratio = _getFinalStateYieldRatios(channel)
	    yield_ratio = yield_ratios[channel]
	    channelName = ['4mu', '4e','2e2mu']
	    file_name = "{1}/Dsignal_{0}.root".format(channelName[channel-1], os.path.dirname(templateSigName))
	    self.factors.update(_readFactorsFromFile(file_name, termNames,channel))
	    
            #ROOT.RooRealVar('{0}_{1}'.format(factor,appendName),factors[factor])

	    #den_segment_geoloc = den_base.format('gamma11_{0}'.format(appendName),'gamma22_{0}'.format(appendName),'gamma33_{0}'.format(appendName),'gamma12_{0}'.format(appendName),'gamma13_{0}'.format(appendName),'gamma23_{0}'.format(appendName))
	    #den_segment_geoloc = den_base.format('1','gamma22_{0}'.format(appendName),'gamma33_{0}'.format(appendName),'gamma12_{0}'.format(appendName),'0'.format(appendName),'gamma23_{0}'.format(appendName))
	    den_segment_geoloc = den_base.format('gamma22_{0}'.format(appendName),'gamma33_{0}'.format(appendName),'gamma12_{0}'.format(appendName),'gamma23_{0}'.format(appendName))
	    den_segment = "{0:.3f}*({1})".format(yield_ratio, den_segment_geoloc)
	    self.log.debug('Denominator segment for ch={0} : {1}'.format(channel, den_segment))
	    return den_segment
	    
	    
	
	
	#read file
	templateSigName = "{0}/Dsignal_{1}.root".format(self.templateDir ,self.appendName)
	
	#factors = self._readFactorsFromFile(templateSigName, termNames)
	appendName = "{0}_{1:.0f}".format(self.channel, self.sqrts)
	self.factors={}
	self.factors.update(_readFactorsFromFile(templateSigName, termNames, self.channel))
        #mean_lamda_13 = (self.factors['lambda13_cosN_{0}'.format(appendName)].getVal()+self.factors['lambda13_cosP_{0}'.format(appendName)].getVal())/2.0
        #self.factors['lambda13_cosP_{0}'.format(appendName)].setVal(mean_lamda_13)
        #self.factors['lambda13_cosN_{0}'.format(appendName)].setVal(mean_lamda_13)
        
	self.log.debug('Factors from file {0}: {1}'.format(templateSigName, str(factors)))
	
	if self.DEBUG: 
	    for fn in self.factors.keys(): 
                #print "{0} = {1}".format(fn, factors[fn])
                self.factors[fn].setConstant()
                print "{0} = {1:.4f}".format(fn, self.factors[fn].getVal())
                #getattr(self.w,'import')(self.factors[fn],ROOT.RooFit.RecycleConflictNodes())
	
	yield_ratios = _getFinalStateYieldRatios()
	#todo: make denominator a function of self.termNames, so to excludde parts which are =0. Use dictionary...
	#todo: check if we have to add factor of 2 in front of gamma12 (= gamma21)... Ask Pedja!!!
	#denominator_geoloc = "{0}+{1}*@0*@0+{2}*@1*@1+({3}*@0)+({4}*@1)+({5}*@0*@1)"
	#denominator_geoloc = "{0}+{1}*@0*@0+{2}*@1*@1+2*(({3}*@0)+({4}*@1)+({5}*@0*@1))" ### added factor 2 for mixed terms
	#denominator_geoloc = "{0}+{1}*@0*@0+{2}*@1*@1+(({3}*@0)+({4}*@1)+({5}*@0*@1))" ### added factor 2 for mixed terms
	#denominator_geoloc = "{0:.3f}+{1:.3f}*@0*@0+{2:.3f}*@1*@1+2*(({3:.3f}*@0)+({4:.3f}*@1)+({5:.3f}*@0*@1))" ### added factor 2 for mixed terms
	#denominator_geoloc = "{0:.3f}+{1:.3f}*@0*@0+{2:.3f}*@1*@1+1*(({3:.3f}*@0)+({4:.3f}*@1)+({5:.3f}*@0*@1))" ### added factor 2 for mixed terms
	#denominator = denominator_geoloc.format(self.factors['gamma11'],self.factors['gamma22'],self.factors['gamma33'],
                                                #self.factors['gamma12'],self.factors['gamma13'],self.factors['gamma23'])
                                                
        #denominator_geoloc = "1+{0}*@0*@0+{1}*@1*@1+(({2}*@0)+({3}*@0*@1))" ### added factor 2 for mixed terms
                                                
                                                
        denominator_geoloc = "1+{0}*@0*@0+{1}*@1*@1+2*(({2}*@0)+({3}*@0*@1))" ### removed insertion of numbers that are not needed, i.e gamma11 and gamma13
                                                
                                                
	denominator = "{0}+{1}+{2}".format(_getDenominatorSegment(1,denominator_geoloc),_getDenominatorSegment(2,denominator_geoloc),_getDenominatorSegment(3,denominator_geoloc))	
	self.log.debug('denominator = {0}'.format(denominator))
	#effectively denominator is = "g11+g33*@1*@1"
	
	
	
	removeGeolocDenominator = False								
	try: 
            self.inputs['user_option']
	except:
	    print "No user option defined... Will not remove Geolocating denominator."
	    removeGeolocDenominator = False    
	else: 
	    if "NoDenominator" in self.inputs['user_option'] : 
                removeGeolocDenominator = True
	
	if removeGeolocDenominator: 
            denominator = "1"
	
	self.log.debug("Geolocating Denominator = {0}".format(denominator))
	
	#r_fs  = self._getFinalStateYieldRatios(self.rfvSigRate_all)
	r_fs  = yield_ratios[self.channel]
	#r_fs = 1
	self.log.debug('FinalState yield ratio: FS = {0} value = {1}'.format(self.channel, r_fs))
	
        #mean_lamda_13 = (self.factors['lambda13_cosN_{0}'.format(appendName)].getVal()+self.factors['lambda13_cosP_{0}'.format(appendName)].getVal())/2
        #self.factors['lambda13_cosP_{0}'.format(appendName)].setVal(mean_lamda_13)
        #self.factors['lambda13_cosN_{0}'.format(appendName)].setVal(mean_lamda_13)
        
        
	nominator={
          'ggH'      :           '{0}'.format('gamma11_{0}'.format(appendName)),
          'gg0Ph'    :     '{0}*@0*@0'.format('gamma22_{0}'.format(appendName)),     # @0 = k2k1_ratio
          'gg0M'     :     '{0}*@1*@1'.format('gamma33_{0}'.format(appendName)),     # @1 = k3k1_ratio
          'ggInt_12P':        '{0}*@0'.format('lambda12_cosP_{0}'.format(appendName)),  
          'ggInt_12N':   '{0}*@0*(-1)'.format('lambda12_cosN_{0}'.format(appendName)),
          'ggInt_13P':        '{0}*@1'.format('lambda13_cosP_{0}'.format(appendName)),  
          'ggInt_13N':   '{0}*@1*(-1)'.format('lambda13_cosN_{0}'.format(appendName)),  
          'ggInt_23P':     '{0}*@0*@1'.format('lambda23_cosP_{0}'.format(appendName)),  
          'ggInt_23N':'{0}*@0*@1*(-1)'.format('lambda23_cosN_{0}'.format(appendName))  
        }
	
	self.rfvNorms={}									
	self.rfvNorms_denom={}
	self.rfvNorms_nom={}
	
	
	postfix = "{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        for one_ch in [1,2,3]:
            denominator_formula = _getDenominatorSegment(one_ch,denominator_geoloc)
            self.ral_denominator = ROOT.RooArgList(self.k2k1_ratio,self.k3k1_ratio)
            for fn in self.factors.keys():
                if fn in denominator_formula:
                    self.ral_denominator.add(self.factors[fn])
            if self.DEBUG:
                self.log.debug('RooArgList for denominator:')
                self.ral_denominator.Print()    
            
            #denominator segment depends upon channel and sqrts
            self.rfvNorms_denom['ch{0}'.format(one_ch)] = ROOT.RooFormulaVar("geoloc_denom_segmentCh{0}_{1:.0f}".format(one_ch,self.sqrts),denominator_formula,self.ral_denominator)
            #getattr(self.w,'import')(self.rfvNorms_denom['ch{0}'.format(one_ch)],ROOT.RooFit.RecycleConflictNodes())
                
	
	
	
	for term in termNames:
	    
            self.ral_nominator = ROOT.RooArgList(self.k2k1_ratio,self.k3k1_ratio)
            for fn in self.factors.keys():
                if fn in nominator[term]:
                    self.ral_nominator.add(self.factors[fn])
            if self.DEBUG:
                self.log.debug('RooArgList for nominator:')
                self.ral_nominator.Print()    
            #nominator segment depends upon term, channel,sqrts
            self.rfvNorms_nom[term] = ROOT.RooFormulaVar("geoloc_nomin_{0}_{1}".format(term, postfix),"{0}".format(nominator[term]),self.ral_nominator)    
            #getattr(self.w,'import')(self.rfvNorms_nom[term],ROOT.RooFit.RecycleConflictNodes())

            self.ral_geolocation_norm = ROOT.RooArgList(self.rfvNorms_nom[term], 
                                                        self.rfvNorms_denom['ch1'],
                                                        self.rfvNorms_denom['ch2'],
                                                        self.rfvNorms_denom['ch3'])
            self.rfvNorms[term] = ROOT.RooFormulaVar("geoloc_{0}_norm_{1}".format(term, postfix),"(@0)/(@1+@2+@3)", self.ral_geolocation_norm)
            #self.rfvNorms[term] = ROOT.RooFormulaVar("geolocation_{0}_norm_{1}".format(term, postfix),"(@0)", self.ral_geolocation_norm)
            getattr(self.w,'import')(self.rfvNorms[term],ROOT.RooFit.RecycleConflictNodes())
	    #rfvNorms[term] = ROOT.RooFormulaVar("geolocation_{0}_norm_{1}".format(term, postfix),"({0})/({1})".format(nominator[term], denominator),ROOT.RooArgList(self.k2k1_ratio, self.k3k1_ratio))
	    #rfvNorms[term] = ROOT.RooFormulaVar("geolocation_{0}_norm_{1}".format(term, postfix),"({0}*{2})/({1})".format(nominator[term], denominator, r_fs),ROOT.RooArgList(self.k2k1_ratio, self.k3k1_ratio))
	    if self.DEBUG:
		self.rfvNorms[term].Print()
		
            #if self.DEBUG : 
                #for i_k2k1_ratio in range(-3,4):
                    #self.k2k1_ratio.setVal(i_k2k1_ratio)
                    #print "k2k1_ratio = ", self.k2k1_ratio.getVal()
                    #self.rfvNorms[term].Print('V')
		


                                                                    
    def reflectSystematics(self,nomShape,altShape):
        assert (nomShape.GetNbinsX()==altShape.GetNbinsX()) and  (nomShape.GetNbinsY()==altShape.GetNbinsY()),\
                "AHHHHHHHHHHH, templates don't have the same binning!!!! X: {0} vs {1} and Y: {2} vs {3}".format(nomShape.GetNbinsX(),altShape.GetNbinsX(),nomShape.GetNbinsY(),altShape.GetNbinsY())
        #if(nomShape.GetNbinsX()!=altShape.GetNbinsX() or nomShape.GetNbinsY()!=altShape.GetNbinsY()):
            #print nomShape.GetNbinsX(),altShape.GetNbinsX(),nomShape.GetNbinsY(),altShape.GetNbinsY()
            #print "AHHHHHHHHHHH, templates don't have the same binning!!!!"
            #return 0
        
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


