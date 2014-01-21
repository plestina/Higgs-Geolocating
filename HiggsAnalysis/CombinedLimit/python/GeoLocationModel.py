from HiggsAnalysis.CombinedLimit.PhysicsModel import *

# mainly based off of FloatingXSHiggs
class GeoLocationModel(PhysicsModel):
    """combine module to setup geolocationg parametrization for fitting"""
    def __init__(self,mode):
        PhysicsModel.__init__(self)        
        self.muAsPOI    = False
        self.muFloating = False
        self.pois      =  ['k3k1_ratio','k2k1_ratio']
        self.mode      = mode
	if ( self.mode == 'K1andK3' ): 
            self.pois      =  ['k3k1_ratio']
	if ( self.mode == 'K1andK2' ): 
            self.pois      =  ['k2k1_ratio']
        print "GeoLocationModel:  ", self.mode
        self.k_ratio_ranges = {'k3k1_ratio':['0','100'],
				'k2k1_ratio':['0','100']}
	self.mHRange=[]
        
        self.verbose = False
        
        
    #things coming in from the command line
    def setPhysicsOptions(self,physOptions):
        """make the POI"""
        print "PhysicsOptions: "+ str(physOptions)
        for po in physOptions:
            if po.startswith("higgsMassRange="):
                self.mHRange = po.replace("higgsMassRange=","").split(",")
                if len(self.mHRange) != 2:
                    raise RuntimeError, "Higgs mass range definition requires two extrema"
                elif float(self.mHRange[0]) >= float(self.mHRange[1]):
                    raise RuntimeError, "Higgs mass range: Extrema for Higgs mass range defined with inverterd order. Second must be larger the first"                
            if po.startswith("poi="):
                self.pois = po.replace("poi=","").split(",")
                
            if po == "muAsPOI": 
                print "Will consider the signal strength as a parameter of interest"
                self.muAsPOI = True
                self.muFloating = True
            if po == "muFloating": 
                print "Will consider the signal strength as a floating parameter (as a parameter of interest if --PO muAsPOI is specified, as a nuisance otherwise)"
                self.muFloating = True    
            #process the relevant POIs
            for poi in self.pois:
                if po.startswith("range_%s"%poi):
                    self.k_ratio_ranges[poi] = po.replace\
                                                      ("range_%s="%poi,"").\
                                                      split(",")
		    print "Range: {0},{1}".format(self.k_ratio_ranges[poi][0],self.k_ratio_ranges[poi][1])
		    
		    if len(self.k_ratio_ranges[poi]) != 2:
			raise RuntimeError, "K-parameter ratios range definition requires two extrema"
		    if float(self.k_ratio_ranges[poi][0]) >= float(self.k_ratio_ranges[poi][1]):
			raise RuntimeError, "K-parameter ratios range: Extrema for K-parameter ratios range defined with inverterd order. Second must be larger the first"

      
    def doParametersOfInterest(self):
      
      

        #if ( self.mode == 'K1andK3' ):
            #poi='k3k1_ratio'
        #elif ( self.mode == 'K1andK2' ):
            #poi='k2k1_ratio'
        #elif ( self.mode == 'K1andK2andK3' ):
            #poi='k2k1_ratio,k3k1_ratio'
        #else:
            #raise RuntimeError('InvalidKParameterChoice',
                               #'We can only use k3k1_ratio, k2k1_ratio'\
                               #' as POIs right now!')      
      
	poi = ",".join(self.pois)
	print "POI 1: "+str(poi)
	
        #if self.pois: poi = self.pois
         # --- Higgs Mass as other parameter ----
        if self.modelBuilder.out.var("MH"):
            if len(self.mHRange):
                print 'MH will be left floating within', self.mHRange[0], 'and', self.mHRange[1]
                self.modelBuilder.out.var("MH").setRange(float(self.mHRange[0]),float(self.mHRange[1]))
                self.modelBuilder.out.var("MH").setConstant(False)
                poi+=',MH'
            else:
                print 'MH will be assumed to be', self.options.mass
                self.modelBuilder.out.var("MH").removeRange()
                self.modelBuilder.out.var("MH").setVal(self.options.mass)
        else:
            if len(self.mHRange):
                print 'MH will be left floating within', self.mHRange[0], 'and', self.mHRange[1]
                self.modelBuilder.doVar("MH[%s,%s]" % (self.mHRange[0],self.mHRange[1]))
                poi+=',MH'
            else:
                print 'MH (not there before) will be assumed to be', self.options.mass
                self.modelBuilder.doVar("MH[%g]" % self.options.mass)

        print "POI 2: "+str(poi)       
        self.pois = poi.split(",")        
	for thepoi in self.pois:
	    if thepoi=="k3k1_ratio" or thepoi=="k2k1_ratio" :
		lower = self.k_ratio_ranges[thepoi][0]
		upper = self.k_ratio_ranges[thepoi][1]
		if self.modelBuilder.out.var(thepoi):
		    self.modelBuilder.out.var(thepoi).setRange(float(lower),float(upper))
		    self.modelBuilder.out.var(thepoi).setConstant(False)
		else:  
		    self.modelBuilder.doVar('%s[%s,%s]'%(thepoi,lower,upper))
        
	#set other ratio to const=0
        if ( self.mode == 'K1andK3' ):
            self.modelBuilder.out.var('k2k1_ratio').setVal(0)
            self.modelBuilder.out.var('k2k1_ratio').setConstant(True)
        elif ( self.mode == 'K1andK2' ):
            self.modelBuilder.out.var('k3k1_ratio').setVal(0)
            self.modelBuilder.out.var('k3k1_ratio').setConstant(True)  
        elif ( self.mode == 'K1andK2andK3' ):
            print "Model with nonzero k1, k2 and k3."
           
        else:
            raise RuntimeError('InvalidKParameterChoice',
                               'We can only use k3k1_ratio, k2k1_ratio'\
                               ' as POIs right now!')
        
	self.modelBuilder.doVar("r[1,0,4]");
	if self.muAsPOI:
	    print 'Treating r as a POI'
	    poi += ",r"

        print "Model POIs : "+str(poi)
        print poi
        self.modelBuilder.doSet('POI',poi)
        
        # display the glory of our work
        self.modelBuilder.out.Print()

        
        
    def getYieldScale(self,bin,process):
	if not self.DC.isSignal[process]: return 1
	return "r"
        
        
        
K1andK3Model 		= GeoLocationModel('K1andK3')
K1andK2Model 		= GeoLocationModel('K1andK2')
K1andK2andK3Model 	= GeoLocationModel('K1andK2andK3')
        
        
        
    
