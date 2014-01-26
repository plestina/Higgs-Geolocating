#import sys
#import os
#import re
#import math
from ROOT import *
import ROOT
from array import array
#import shutil
#import yaml
import pprint

import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir)))
from lib.util.RootAttributeTranslator import *
from lib.util.Logger import *
from lib.util.UniversalConfigParser import *




class PlotPolisher(object):
  DEBUG = False
  def __init__(self):
      self.pp = pprint.PrettyPrinter(indent=4)
      
  def setStyle(self,root_object, style={}):
      """
      Set color, marker, fill, style from dictionary
      """
      #check keys that are presentand set boolfor them
      #set attribute toggles to be offinitially
      att_toggle = {'linecolor':0, 'linewidth':0,'linestyle':0, 'markercolor':0, 'markersize':0, 'markerstyle':0, 'fillcolor':0, 'fillstyle':0}
      if self.DEBUG: 
	  print "@@@@ Style in PlotPolisher: ", 
	  self.pp.pprint(style)
      
      
      self.n_colors = 0
      for att in style.keys(): 
	  #try:
	    #att_toggle[att]
	  #except ValueError:
	    #print "The atribute %(att)s cannot be attached to your plot.Not implemented. Check in setStyle()." %locals()
	  #else:  
	    #att_toggle[att] = True
	    
	  if att in att_toggle.keys() : att_toggle[att] = True
	    
	  if 'color' in att: 
		self.n_colors+=1
		theColor = style[att]
	    
      
      if att_toggle['linecolor'] :   root_object.SetLineColor(style['linecolor'])
      if att_toggle['linewidth'] :   root_object.SetLineWidth(style['linewidth'])
      if att_toggle['linestyle'] :   root_object.SetLineStyle(style['linestyle']) 
      if att_toggle['markercolor'] : root_object.SetMarkerColor(style['markercolor'])
      if att_toggle['markersize'] :  root_object.SetMarkerSize(style['markersize'])
      if att_toggle['markerstyle'] : root_object.SetMarkerStyle(style['markerstyle'])
      if att_toggle['fillcolor'] :   root_object.SetFillColor(style['fillcolor']) 
      if att_toggle['fillstyle'] :   root_object.SetFillStyle(style['fillstyle']) 
      
      #set color of only one is provided
      if self.n_colors==1: 
	  root_object.SetLineColor(theColor)
	  root_object.SetMarkerColor(theColor)
	  root_object.SetFillColor(theColor)
	  
  def toggle_options(self,setup_keys,att_toggle):
      """
      Turns on options that are provided by user.
      """
      if self.DEBUG: print setup_keys, att_toggle
      for att in setup_keys: 
	  if self.DEBUG: print "Att = ", att
	  if att in att_toggle.keys() : att_toggle[att] = True
	    
      if self.DEBUG : print "Options used: ", att_toggle
      return att_toggle	
      
  def arrangeAxis(self, root_object,setup={}):
      att_toggle = {'x_axis':0, 'y_axis':0}
      att_toggle = self.toggle_options(setup.keys(),att_toggle)
      if self.DEBUG: print "Options used (after toggle function): ", att_toggle
      
      
      axis_att_toggle = {'title':0,'range':0 }
      if att_toggle['x_axis'] : 
	 axis = root_object.GetXaxis()
	 this_setup = setup['x_axis']
	 x_axis_att_toggle = self.toggle_options(this_setup.keys(),axis_att_toggle )
	 if x_axis_att_toggle['title']: axis.SetTitle(this_setup['title'])
	 if x_axis_att_toggle['range']: axis.SetRangeUser(this_setup['range'][0],this_setup['range'][1])
      
      if att_toggle['y_axis'] : 
	 axis = root_object.GetYaxis()
	 this_setup = setup['y_axis']
	 this_axis_att_toggle = self.toggle_options(this_setup.keys(),axis_att_toggle )
	 if this_axis_att_toggle['title']: axis.SetTitle(this_setup['title'])
	 if this_axis_att_toggle['range']: axis.SetRangeUser(this_setup['range'][0],this_setup['range'][1])
	 
  #def add_text(self, canvas, text_list):
  def add_text(self, text_list):
	 t = ROOT.TLatex()
	 
	 rat = RootAttributeTranslator()
	 for text_entry in text_list:
	    t.SetTextFont( 62 )
	    t.SetTextSize( 0.025 )
	    t.SetTextAlign( 10 )
	    t.SetTextColor( kBlack )
	    text_entry = rat.translate_all(text_entry)
	    this_att_toggle = self.toggle_options(text_entry.keys(),{'text':0, 'position':0, 'size':0, 'font':0, 'color':0, 'angle':0, 'align':0})
	    if this_att_toggle['size'] :  t.SetTextSize(float(text_entry['size']))
	    if this_att_toggle['font'] :  t.SetTextFont(int(text_entry['font']))
	    if this_att_toggle['color'] : t.SetTextColor(int(text_entry['color']))
	    if this_att_toggle['align'] : t.SetTextAlign(int(text_entry['align']))
	    if this_att_toggle['position'] :
		x,y = text_entry['position'][:2]
		#x = gPad.GetFrame().GetX1()
		#y = gPad.GetFrame().GetY2()
		if self.DEBUG: print "@@@@ Text coordinates: ",x,y
		
		text_entry['text'] = text_entry['text'].replace('\\','#')
		if self.DEBUG: print "@@@@ Painting text : ", text_entry['position'][:2], text_entry['text']
		
		if 'NDC' or 'ndc' in text_entry['position']: 
		    t.SetNDC(1)
		    t.DrawLatex(x,y,text_entry['text'])
		    t.SetNDC(0)
		else: t.DrawLatex(x,y,text_entry['text'])
		
  def arrangeCanvas(self, canvas, setup={}):
      """
      Setup canvas, axis, legend title and position...
      """
      #check keys that are presentand set boolfor them
      #set attribute toggles to be offinitially
      att_toggle = {'grid':0, 'legend':0,'x_axis':0, 'y_axis':0,'title':0,'size':0, 'save_ext':0}
      att_toggle = self.toggle_options(setup.keys(),att_toggle)
      if self.DEBUG: print "Options used (after toggle function): ", att_toggle
      
      if self.DEBUG: 
	  print "@@@@ Setup in PlotPolisher: ", 
	  self.pp.pprint(setup)
      
      
      if att_toggle['grid'] :   
	 if setup['grid']==1 : canvas.SetGrid(bool(setup['grid']))
	 if 'x' in str(setup['grid']).lower(): canvas.SetGridx(1)
	 if 'y' in str(setup['grid']).lower(): canvas.SetGridy(1)
      
      axis_att_toggle = {'log':0}
      if att_toggle['x_axis'] : 
	 this_setup = setup['x_axis']
	 x_axis_att_toggle = self.toggle_options(this_setup.keys(),axis_att_toggle )
	 if x_axis_att_toggle['log']:   canvas.SetLogx(bool(this_setup['log']))
      
      if att_toggle['y_axis'] : 
	 this_setup = setup['y_axis']
	 this_axis_att_toggle = self.toggle_options(this_setup.keys(),axis_att_toggle )
	 if this_axis_att_toggle['log']:   canvas.SetLogy(bool(this_setup['log']))
	 
	
      if att_toggle['legend'] : 
	 this_setup = setup['legend']
	 this_att_toggle = self.toggle_options(this_setup.keys(),{'title':0, 'position':0})
	 self.leg = ROOT.TLegend(0.18,0.15,0.8,0.4);
	 self.leg.SetBorderSize(0)
	 self.leg.SetFillColor(-1)
	 if this_att_toggle['position'] :
	    assert len(this_setup['position'])==4, "The legend corners should have x_min,y_min,x_max,y_max values. There shouldbe 4 values."
	    x1,y1,x2,y2 = this_setup['position']
	    self.leg.SetX1NDC(x1)
	    self.leg.SetY1NDC(y1)
	    self.leg.SetX2NDC(x2)
	    self.leg.SetY2NDC(y2)
	 if this_att_toggle['title'] : self.leg.SetHeader(str(this_setup['title']))
	 
	 
class SimplePlotter(PlotPolisher):
  DEBUG = True
  def __init__(self,name = "My Graph"):
      self.style = {'linecolor' : kBlack, 'linestyle':1, 'linewidth':2, 'markersize':0.5, 'markerstyle':20}
      self.xtitle = "x"
      self.ytitle = "y"
      self.pp = pprint.PrettyPrinter(indent=4)
      
  def getGraph(self, X, Y, user_style = {}):
      #make default style but recieve an updated dict for the style
      self.style.update(user_style)
      print "@@@@ Graph Style = "+str(self.style)
      npoints = len(X)
      #self.npoints = 15
      self.gr = ROOT.TGraph(npoints,X,Y)
      self.gr.GetXaxis().SetTitle(self.xtitle)
      self.gr.GetYaxis().SetTitle(self.ytitle)
      
      #self.gr.GetXaxis().SetMoreLogLabels(1)
      #self.gr.GetYaxis().SetMoreLogLabels(1)
      
      
      if self.DEBUG: print "@@@@ Style in SimplePlotter: ", self.style
      self.setStyle(self.gr, self.style)
      return self.gr
      
  def getTH1(self,source, user_style={}) :   
      pass
      return hist

  def getTH2(self,filename, histname, user_style={}) :   
      #WARNING This function is not working yet.
      pass
      print "getTH2: Running on file:{0}  template:{1}".format(filename, histname)
      f = TFile(filename,"READ")
      import copy
      self.th2 = copy.deepcopy(f.Get(histname))
      print "Xaxis bins: ", self.th2.GetXaxis().GetNbins()
      print "Yaxis bins: ", self.th2.GetYaxis().GetNbins()
      canv = TCanvas("test","test", 300,200)
      canv.cd()
      self.th2.Draw("colz") 
      
      import os
      print "The filename path exists = ",os.path.exists(filename)
      print "Getting the histogram {0} from file:{1}".format(histname,filename)
      
      return self.th2
      
      
