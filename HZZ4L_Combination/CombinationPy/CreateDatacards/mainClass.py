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
from kdClass import *
from superkdClass import *
from kParamDiscriminantClass import *
class mainClass():
    
    ID_4mu = 1
    ID_4e  = 2
    ID_2e2mu = 3
    
    

    def __init__(self):

        ROOT.gSystem.AddIncludePath("-I$ROOFITSYS/include/")
        ROOT.gSystem.AddIncludePath("-Iinclude/")
        ROOT.gROOT.ProcessLine(".L include/tdrstyle.cc")
        ROOT.gSystem.Load("libRooFit")
        ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
        ROOT.gSystem.Load("include/HiggsCSandWidth_cc.so")
        ROOT.gSystem.Load("include/HiggsCSandWidthSM4_cc.so")
        

    def makeCardsWorkspaces(self, theMH, theOutputDir, theInputs,theTemplateDir,theMassError,theis2D,theUseMEKD,altHypo):


        ## ------------------ CHECK CHANNEL ------------------- ##
        if (theInputs['decayChannel'] == self.ID_4mu): appendName = '4mu'
        elif (theInputs['decayChannel'] == self.ID_4e): appendName = '4e'
        elif (theInputs['decayChannel'] == self.ID_2e2mu): appendName = '2e2mu'
        else: print "Input Error: Unknown channel! (4mu = 1, 4e = 2, 2e2mu = 3)"
        
        print '>>>>>> Decay Channel: ',theInputs['decayChannel'],' = ',appendName

        ## ----------------- WIDTH AND RANGES ----------------- ##
        #myCSW = HiggsCSandWidth()
        #widthHVal =  myCSW.HiggsWidth(0,theMH)

        #self.windowVal = max( widthHVal, 1.0)
        #self.windowVal = max( widthHVal, 1.0)
        #lowside = 100.0
        #highside = 1000.0
        #if (theMH >= 275):
            #lowside = 180.0
            #highside = 650.0
        #if (theMH >= 350):
            #lowside = 200.0
            #highside = 900.0
        #if (theMH >= 500):
            #lowside = 250.0
            #highside = 1000.0
        #if (theMH >= 700):
            #lowside = 350.0
            #highside = 1400.0
                        
        #self.low_M = max( (theMH - 20.*self.windowVal), lowside)
        #self.high_M = min( (theMH + 15.*self.windowVal), highside)
               
        #print '>>>>>> Higgs Width: ',widthHVal,' Window --- Low: ', self.low_M,' High: ', self.high_M

        ### add to the inputs
        #theInputs['low_M'] = self.low_M
        #theInputs['high_M'] = self.high_M
        
        

        ## ------------- Signal Separation ----------- ##
        theInputs['altHypothesis'] = altHypo
        if theInputs['doHypTest'] and theInputs['altHypLabel']=="" :
            theInputs['altHypLabel'] = "_ALT"

        if theInputs['unfold'] and not theInputs['doHypTest']:
            raise RuntimeError, "Cannot unfold 2D into 1D if not doing Hypothesis test!"

        if theInputs['doHypTest']:
            print '----------------- Running Signal Hypothesis Test Cards -----------------'
            print '>>>>>> Alt Hypothesis: ',theInputs['altHypothesis']
            print '>>>>>> Alt Label: ',theInputs['altHypLabel']

        if theInputs['unfold']:
            print '>>>>>> Unfolding 2D into 1D hists'
                        
        if theInputs['doHypTest'] and not theInputs['all']:
            raise RuntimeError, "You asked to prepare DC and WS for Hyp Test but you did not want to sum over all signal channels. This is forbidden. Check inputs ! (it should have already send you this error message, strange that  you are here...)"
        
        if (theInputs['doHypTest'] and not (theis2D==1)):
            raise RuntimeError, "Cannot perform hypothesis testing without a 2D analysis, feature not supported yet. Exiting."
        
        
        myDatacardClass = datacardClass()

        if( theis2D == 0 ):

            myDatacardClass.makeMassShapesYields(theMH,theOutputDir,theInputs,theTemplateDir,theMassError,theis2D,theUseMEKD)
            myDatacardClass.fetchDataset()
            myDatacardClass.writeWorkspace()
            myDatacardClass.prepareDatacard()

        if( theis2D == 1 ):

	    useSuperKD = False
            if theInputs['doHypTest'] and useSuperKD:

                myDatacardClass2D = superkdClass()
                myDatacardClass2D.makeMassShapesYields(theMH,theOutputDir,theInputs,theTemplateDir,theMassError,theis2D,theUseMEKD)
                
                myDatacardClass2D.setSuperKD()
                myDatacardClass2D.fetchDatasetSuperKD()            
                myDatacardClass2D.makeSuperKDAnalysis()
                if theInputs['unfold']:
                    myDatacardClass2D.writeWorkspaceUnfoldedSuperKD() 
                else:
                    myDatacardClass2D.writeWorkspaceSuperKD()
		
                myDatacardClass2D.prepareDatacardSuperKD()  
                
            elif theInputs['doHypTest']:
                ### TODO remove the swich for hypothess testing because we are doing MLF now
	      
		print "--------> Using kParamDiscriminantClass to prepare datacards. <----------"
                myDatacardClass2D = kParamDiscriminantClass()
                myDatacardClass2D.makeMassShapesYields(theMH,theOutputDir,theInputs,theTemplateDir,theMassError,theis2D,theUseMEKD)
		try:
		    theInputs['termNames']
		except NameError:	
		    print "@@@@ Setting termNames collection to default!"
		    myDatacardClass2D.setTermNames()
		else:
		    print "@@@@ Setting termNames collection to user values!"
		    myDatacardClass2D.setTermNames(theInputs['termNames'])
		    
                myDatacardClass2D.setSuperKD()
                myDatacardClass2D.fetchDatasetSuperKD()            
                myDatacardClass2D.makeSuperKDAnalysis()
                if theInputs['unfold']:
                    myDatacardClass2D.writeWorkspaceUnfoldedSuperKD() 
                else:
		    #myDatacardClass2D.writeWorkspace()

                    myDatacardClass2D.writeWorkspaceSuperKD()
		
                myDatacardClass2D.prepareDatacardSuperKD()    
                

            else:
            
                myDatacardClass2D = kdClass()
                myDatacardClass2D.makeMassShapesYields(theMH,theOutputDir,theInputs,theTemplateDir,theMassError,theis2D,theUseMEKD)
                
                myDatacardClass2D.setKD()
                myDatacardClass2D.fetchDatasetKD()            
                myDatacardClass2D.makeKDAnalysis()
                myDatacardClass2D.writeWorkspaceKD()
                myDatacardClass2D.prepareDatacardKD()
                
            

            
        

