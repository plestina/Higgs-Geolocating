#! /usr/bin/env python
#import sys
#import os
#import re
#import math
from ROOT import *
import ROOT
from array import array
import shutil
#import yaml
import pprint
from lib.util.RootAttributeTranslator import RootAttributeTranslator 
from lib.util.Logger import Logger
from lib.util.UniversalConfigParser import UniversalConfigParser
from lib.plotting.RootPlotters import PlotPolisher,SimplePlotter
from lib.plotting.RootPlottersBase import RootPlottersBase


class LimitVsLumi(PlotPolisher,RootPlottersBase):  
  def __init__(self,name = "limitVsLumi_plot" ):
      self.log = Logger().getLogger(self.__class__.__name__, 10)
      self.name = name
      #ROOT.gSystem.AddIncludePath("-I$ROOFITSYS/include/");
      ROOT.gROOT.ProcessLine(".L tdrstyle.cc")
      from ROOT import setTDRStyle
      ROOT.setTDRStyle(True)
      
      ROOT.gStyle.SetPalette(1)
      ROOT.gStyle.SetOptStat(0)
      #self.copy_to_web_dir = False
      #self.webdir = ""
      self.pp = pprint.PrettyPrinter(indent=4)
      
      
      
  def setName(self, newname): self.name = newname

   
  def makePlot(self, data):
      #self.ytitle = "U.L. @95% for |k_{3}/k_{1}|" 
      #self.xtitle = "L (fb^{-1})"

      self.c1 =ROOT.TCanvas("cc1","95% CL limits",800,800)
      self.c1.cd()
      ROOT.gPad.SetRightMargin(0.0085)
      
      try:
	  data['content']
      except:
	  raise ValueError, "Canvas \'content\' dictionary is not provided in config file."
	  
      try:
	  data['setup']
      except:
	  print "@@@@ Canvas \'setup\' dictionary is not provided. "
	  self.setup_exist=False
      else:
	  self.arrangeCanvas(self.c1, data['setup'])
	  self.setup_exist=True
	  
      
      try:
	  self.leg
      except AttributeError:
	  self.leg_exist=False
      else:
	  self.leg_exist=True
      
      try:
	  data['setup']['add_text']
      except AttributeError:
	  self.add_text_exist=False
      else:
	  self.add_text_exist=True
      
      i=0
      rat = RootAttributeTranslator()
      for theGraph in data['content']:
	    #theGraph = data[key]
	    style=rat.translate_all(theGraph['style'])
	    #style=theGraph['style']
	    
	    #x_values = array('d', self.getColumnFromTable(theGraph['table'],'lumi'))
	    #y_values = array('d', self.getColumnFromTable(theGraph['table'],'ul95'))
	    x_values = array('d', self._getColumnFromTable(theGraph['table'],'lumi'))
            y_values = array('d', self._getColumnFromTable(theGraph['table'],'ul95'))
            
	    self.log.debug("(x,z) = "+str(zip(x_values,y_values)))
	    
	    #theGraph['graph'] = self.getGraph(x_values,y_values, style)
	    gr = SimplePlotter()
	    theGraph['graph'] = gr.getGraph(x_values,y_values, style)
	    if self.setup_exist: self.arrangeAxis(theGraph['graph'],data['setup'])
	    if self.leg_exist:
		self.leg.AddEntry(theGraph['graph'],theGraph['legend']['text'],theGraph['legend']['opt']);
	    try:
		theGraph['draw_opt']
	    except ValueError:
		draw_opt = "LP"
	    else:
		draw_opt = str(theGraph['draw_opt'])
		
	    if i==0: draw_opt += "A"
	    else :   draw_opt += "same"
	    i+=1
	    theGraph['graph'].Draw(draw_opt)
	 
      if self.leg_exist: self.leg.Draw()
      if self.add_text_exist: self.add_text(data['setup']['add_text'])
      
      
      #plot_name = "limitVsLumi_2D_k3k1_0_logy"
      plot_name = self.name
      self.save_extensions = ['png','pdf','eps']
      try:
	  data['setup']['save_ext']
      except ValueError:
	  self.log.info("No extensions are provided in setup. Using default: ", self.save_extensions)
      else:
	  self.save_extensions = list(data['setup']['save_ext'])
      self.save(self.c1, plot_name, self.save_extensions)
  
   

      
  def read_table(self, file_name):
      # read it
      import csv
      table=[]
      with open(file_name, 'r') as csvfile:
	  reader = csv.reader(csvfile,skipinitialspace=True,delimiter=" ")
	  #table = [[float(e) for e in r] for r in reader]
	  table = [[str(e) for e in r] for r in reader]
	  #table = [[e for e in r] for r in reader]
      return table
      
      
  def getColumnFromTable(self, table_file,col_name):
      col_names_list = {'lumi':0,'bf':1,'ul68':2,'ul95':3,'wf':4}
      filename = table_file.replace("_",".")+".lumi.bf.ul68.ul95.wf.tab"
      self.log.debug("Reading column {0}={1} from file = {2} ...".format(col_name,col_names_list[col_name],filename))
      table = self.read_table(filename)
      assert len(table[0]) == len(col_names_list), "Wrong number of columns in table: %(filename)s " %locals()
      #print table
      column = [row[col_names_list[col_name]] for row in table]
      self.log.debug("Values {0} = {1}".format(col_name,column))
      return column
      
     
  def _getColumnFromTable(self, table_file,col_name):
      filename = table_file.replace("_",".")+".tab"
      
      
      table = self.read_table(filename)
      #print table
      #col_names_list = {'lumi':0,'bf':1,'ul68':2,'ul95':3,'wf':4}
      col_names_list = dict(zip(table[0],[i for i in range(len(table[0]))])) #e.g. {'lumi':0,'bf':1,'ul68':2,'ul95':3,'wf':4}
      self.log.debug("Reading column {0}={1} from file = {2} ...".format(col_name,col_names_list[col_name],filename))
      for irow in range(len(table))[1:] : 
            assert len(table[0]) == len(table[irow]), "In table {0}, wrong number of columns in row {1}. Should be {2} ".format(filename, len(table[irow]), len(table[0]))
      #print table
      column = [float(row[col_names_list[col_name.lower()]]) for row in table[1:]]
      self.log.debug("Values {0} = {1}".format(col_name,column))
      return column
      
      
      
    
    
    
    
if __name__ == "__main__":
     
    # parse the arguments and options
    #global opt, args
    #parseOptions()
    print "Starting plots..."
    
    pp = pprint.PrettyPrinter(indent=4)
    
    DEBUG = True
    cfg_reader = UniversalConfigParser(cfg_type="YAML",file_list = "plotLimitVsLumi.yaml")
    cfg_reader.setLogLevel(10) 
    
    plots_data = cfg_reader.get_dict()
    
    if DEBUG: 
	  print "plots_data = ", 
	  pp.pprint(plots_data)
    
    plotter = LimitVsLumi()
    plotter.setCopyToWebDir(True,"/afs/cern.ch/user/r/roko/www/html/Geoloc/")
    ##1str plot
    #plotter.setName("plot_4l_unfolding_withandwitout_bkg_14Tev")
    #plotter.makePlot(plots_data['plot_4l_unfolding_withandwitout_bkg_14Tev'])
    #plotter.setName("plot_4l_unfolding_with_bkg_14Tev_int_effect")
    #plotter.makePlot(plots_data['plot_4l_unfolding_with_bkg_14Tev_int_effect'])
    #plotter.setName("plot_4l_unfolding_withBkg_newDenominator_14Tev")
    #plotter.makePlot(plots_data['plot_4l_unfolding_withBkg_newDenominator_14Tev'])
    
    #plotter.setName("plot_4l_unfolding_withBkg_newDenominator_14Tev_fixWS")
    #plotter.makePlot(plots_data['plot_4l_unfolding_withBkg_newDenominator_14Tev_fixWS'])
    
    for plot_name in plots_data['list_of_plots_to_do']:
	print 'Make plot: ',plot_name
	plotter.setName(plot_name)
	plotter.makePlot(plots_data[plot_name])
   
    
    
    