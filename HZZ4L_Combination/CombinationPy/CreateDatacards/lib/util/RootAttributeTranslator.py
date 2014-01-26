#__author__ = "Roko Plestina"
#__version__ = "0.0"
#__maintainer__ = "Roko Plestina"
#__email__ = "roko.plestina@cern.ch"
#__status__ = "Quick and dirty"	
import ROOT
from ROOT import TFormula


class RootAttributeTranslator(object):
    """
    Class used to translate typical ROOT constant names provided from external configuration files.
    The translator will return the integer value for the attribute which can then be passed to another
    ROOT object...
    """

    def __init__(self):
	pass	 
      
    def init_colors(self):
        self.Colors = { 'kWhite':0,   'kBlack':1,   'kGray':920,
			  'kRed':632, 'kGreen':416, 'kBlue':600, 'kYellow':400, 'kMagenta':616, 'kCyan':432,
			  'kOrange':800, 'kSpring':820, 'kTeal':840, 'kAzure':860, 'kViolet':880, 'kPink':900 }      

    def init_markerstyles(self):
	self.MarkerStyles =  {'kDot':1, 'kPlus':2, 'kStar':3, 'kCircle':4, 'kMultiply':5,
				'kFullDotSmall':6, 'kFullDotMedium':7, 'kFullDotLarge':8,
				'kFullCircle':20, 'kFullSquare':21, 'kFullTriangleUp':22,
				'kFullTriangleDown':23, 'kOpenCircle':24, 'kOpenSquare':25,
				'kOpenTriangleUp':26, 'kOpenDiamond':27, 'kOpenCross':28,
				'kFullStar':29, 'kOpenStar':30, 'kOpenTriangleDown':32,
				'kFullDiamond':33, 'kFullCross':34}
      
      
    def init_linestyles(self):
	self.LineStyles = { 'kSolid': 1, 'kDashed':2, 'kDotted':3, 'kDashDotted':4 }
	
    def init_fillstyles(self):
	"""
	Fill styles as defined at 
	http://root.cern.ch/root/html/TAttFill.html
	Include more from there...
	"""
      
        self.FillStyles = {'kDotted' : 3001, 
			    'kHatched' : 3004, 'kHatchedRaise':3004, 'kHatchedFall':3005, 'kHatchedVertical':3006, 'kHatchedHorizontal':3007}
			    
			    
	
    def translate_all(self,user_dict):
	"""
	Translate a full dictionary with user provided values to 
	integer values that can be used for Root plots. 
	"""
	#{color : 1,  line : 2}
	for att in user_dict.keys():
	    #update dictionary with value translated
	    if 'color' in att.lower():
		user_dict[att] = self.translate('color',user_dict[att])
	    if 'linestyle'==att.lower(): 
		user_dict[att] = self.translate('linestyle',user_dict[att])
	    if 'markerstyle'==att.lower(): 
		user_dict[att] = self.translate('markerstyle',user_dict[att])
	    if 'fillstyle'==att.lower():
		user_dict[att] = self.translate('fillstyle',user_dict[att])
		
	return user_dict    
	    
	
      
      
    def translate(self,att_type, user_att):
      
	if not isinstance(user_att,str): return user_att
	
	attribute_dict = {}
	
	if 'color' in att_type.lower(): 		
	    try:
		self.Colors
	    except:
		self.init_colors()
			
	    attribute_dict = self.Colors
	    
	elif 'markerstyle' in att_type.lower(): 	
	    try:
		self.MarkerStyles
	    except:
		self.init_markerstyles()
		
	    attribute_dict = self.MarkerStyles
	    
	elif 'linestyle' in att_type.lower(): 	
	    try:
		self.LineStyles
	    except:
		self.init_linestyles()
	    attribute_dict = self.LineStyles
	    
	elif 'fillstyle' in att_type.lower(): 	
	    try:
		self.FillStyles
	    except:
		self.init_fillstyles()
	    attribute_dict = self.FillStyles
	    
	    
	for attribute in attribute_dict.keys():
	    if attribute in user_att: user_att = user_att.replace(attribute,str(attribute_dict[attribute]))
	print "User attribute translated = ", user_att
	
	#return int(ROOT.TFormula("f", user_att).Eval(0))
	return int(eval(user_att))
	