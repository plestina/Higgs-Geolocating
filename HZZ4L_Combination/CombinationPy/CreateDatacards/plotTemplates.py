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
from lib.util.RootAttributeTranslator import *
from lib.util.Logger import *
from lib.util.UniversalConfigParser import *
from lib.plotting.RootPlotters import *
import lib.util.MiscTools as misctools

class TemplatesPlotter(PlotPolisher):  
  def __init__(self,name = "templates" ):
      self.log = Logger().getLogger(self.__class__.__name__, 10)
      self.name = name
      #ROOT.gSystem.AddIncludePath("-I$ROOFITSYS/include/");
      ROOT.gROOT.ProcessLine(".L tdrstyle.cc")
      from ROOT import setTDRStyle
      ROOT.setTDRStyle(True)
      
      ROOT.gStyle.SetPalette(1)
      ROOT.gStyle.SetOptStat(0)
      self.copy_to_web_dir = False
      self.webdir = ""
      self.pp = pprint.PrettyPrinter(indent=4)
      self.k2k1_ratio=0
      self.k3k1_ratio=0
      self.total_pdf=0
      
  def set_k2k1_and_k3k1(self, k2k1, k3k1):
      self.k2k1_ratio = k2k1
      self.k3k1_ratio = k3k1
      
  def setName(self, newname): self.name = newname
  def get_scale_factors(self):
            
        factors = {'lambda12_cosP' : 0.0,
            'lambda13_cosP' : 0.0,
            'lambda23_cosN' : 0.0,
            'lambda23_cosP' : 0.0,
            'gamma33' : 0.034,
            'gamma13' : 0.0,
            'gamma12' : -0.269303399267,
            'gamma11' : 1.0,
            'lambda13_cosN' : 0.0,
            'lambda12_sinN' : 0.0,
            'lambda12_cosN' : 0.538606798534,
            'gamma22' : 0.09,
            'gamma23' : 0.0,
            }
            
        nominator={
            'ggH_shape'      :             '{0}'.format(factors['gamma11']),
            'gg0Ph_shape'    :     '{0}*{1}*{1}'.format(factors['gamma22'],self.k2k1_ratio),     # @0 = k2k1_ratio
            'gg0M_shape'     :     '{0}*{1}*{1}'.format(factors['gamma33'],self.k3k1_ratio),     # @1 = k3k1_ratio
            'ggInt_12P_shape':         '{0}*{1}'.format(factors['lambda12_cosP'],self.k2k1_ratio),  
            'ggInt_12N_shape':    '{0}*{1}*(-1)'.format(factors['lambda12_cosN'],self.k2k1_ratio),
            'ggInt_13P_shape':         '{0}*{1}'.format(factors['lambda13_cosP'],self.k3k1_ratio),  
            'ggInt_13N_shape':    '{0}*{1}*(-1)'.format(factors['lambda13_cosN'],self.k3k1_ratio),  
            'ggInt_23P_shape':     '{0}*{1}*{2}'.format(factors['lambda23_cosP'],self.k2k1_ratio,self.k3k1_ratio),  
            'ggInt_23N_shape':'{0}*{1}*{2}*(-1)'.format(factors['lambda23_cosN'],self.k2k1_ratio,self.k3k1_ratio)  
            }
        return nominator


  def set_total_pdf(self, plot, root_file):
        f = TFile(root_file,"READ")
        nominator= self.get_scale_factors()    
        for i_pdf in range(len(plot['histos'])):
            if i_pdf==0:
                total_pdf = f.Get(plot['histos'][i_pdf])
                total_pdf.Scale(eval(nominator[plot['histos'][i_pdf]]))
                self.log.debug('TOTAL_PDF -> Picking up the first term {0} and the scale {1}. The histo is TH2: {2}'.format(plot['histos'][i_pdf], eval(nominator[plot['histos'][i_pdf]]), isinstance(total_pdf,TH2)))

            if plot['histos'][i_pdf]=="TOTAL_PDF":
                break
            another_template = f.Get(plot['histos'][i_pdf])
            another_template.Scale(float(eval(nominator[plot['histos'][1]]))) 
            self.log.debug('TOTAL_PDF -> Adding up term {0} and the scale {1} to histo TH2:{2}'.format(plot['histos'][i_pdf], eval(nominator[plot['histos'][i_pdf]]),isinstance(total_pdf,TH2)))
            total_pdf.Add(another_template)
        self.total_pdf = total_pdf    
        self.log.debug('TOTAL_PDF -> Added all terms and now the histo is TH2:{0} ... returning the value.'.format(isinstance(self.total_pdf,TH2)))
        
  def makePlot(self, data):

      self.c1 =ROOT.TCanvas("cc1","Templates",1000,800)
      self.c1.cd()
      ROOT.gPad.SetRightMargin(0.2)
      
      try:
	  data['content']
      except KeyError:
	  raise KeyError, "Canvas \'content\' dictionary is not provided in config file."
	  
      try:
	  data['setup']
      except KeyError:
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
      print data['setup']['add_text']
      try:
	  data['setup']['add_text']
      except KeyError:
	  self.add_text_exist=False
      else:
	  self.add_text_exist=True
	  
      
      
      i=0
      #rat = RootAttributeTranslator()
      for plot in data['content']:
        try:
            plot['POI']
        except KeyError:
            self.plot_total_pdf=False
        else:
            self.plot_total_pdf=True  
            plot['histos'].append("TOTAL_PDF")
            
            
        
  
        for fs in data['setup']['final_states']:
            print 100*"-"
            theTemplFile = plot['file']+"_"+fs+".root"
            interference_tag = "13"
            
            if ('ggInt_13P_shape' in plot['histos']) and ('ggInt_13N_shape' in plot['histos']): 
                plot['histos'].append("ggInt_13_shape")
                interference_tag = "13"
            elif ('ggInt_12P_shape' in plot['histos']) and ('ggInt_12N_shape' in plot['histos']): 
                plot['histos'].append("ggInt_12_shape")
                interference_tag = "12"
            if self.plot_total_pdf:
                self.set_total_pdf(plot, theTemplFile)
                
                self.log.debug('TOTAL_PDF -> Returned histo is TH2:{0}'.format(isinstance(self.total_pdf,TH2)))
            for theTemplate in plot['histos']:
               
                doSum = False
                if theTemplate == "ggInt_{0}_shape".format(interference_tag) : 
                    doSum = True
                    
                self.log.debug("Running on file:{0}  template:{1}".format(theTemplFile, theTemplate))
                if DEBUG:
                    pp.pprint(plot)
                th2 = SimplePlotter()
                
                #plot['th2'] = th2.getTH2(theTemplFile, theTemplate)
                if theTemplate=="TOTAL_PDF":
                    plot['th2'] = self.total_pdf
                    self.log.debug('TOTAL_PDF -> The addition is over. Going to plot now.' )
                elif not doSum:
                    f = TFile(theTemplFile,"READ")
                    plot['th2'] = f.Get(theTemplate)
                else:
                    self.log.debug('Summing up the interference template, which has been divided into positive and negative part.')
                    f = TFile(theTemplFile,"READ")
                    th2_P = f.Get("ggInt_{0}P_shape".format(interference_tag))
                    th2_N = f.Get("ggInt_{0}N_shape".format(interference_tag))
                        
                    #get lamda factors
                    tree=f.Get('factors')
                    factors={}
                    for fn in ['lambda{0}_cosP'.format(interference_tag),'lambda{0}_cosN'.format(interference_tag)]:
                        factors[fn] = array('d',[0])
                        tree.SetBranchAddress(fn,factors[fn])
                    tree.GetEntry(0)
                    for fn in factors.keys(): factors[fn] = factors[fn][0]
                    print "Lambdas:",factors['lambda{0}_cosP'.format(interference_tag)], factors['lambda{0}_cosN'.format(interference_tag)]
                    th2_P.Scale(factors['lambda{0}_cosP'.format(interference_tag)])
                    th2_N.Scale(factors['lambda{0}_cosN'.format(interference_tag)])
                    th2_P.Add(th2_N, -1)
                    plot['th2'] = th2_P
                 
                 
                 
                    
                if self.setup_exist: self.arrangeAxis(plot['th2'],data['setup'])
                if self.leg_exist:
                    self.leg.AddEntry(plot['th2'],plot['legend']['text'],plot['legend']['opt']);
                try:
                    plot['draw_opt']
                except KeyError:
                    draw_opt = "COLZ"
                else:
                    draw_opt = str(plot['draw_opt'])
                if draw_opt.lower()=="surf":
                    self.c1.SetTheta(16.56)
                    self.c1.SetPhi(57.83133)

                plot['th2'].Draw(draw_opt)
                
                if self.leg_exist: self.leg.Draw()
                
                
                template_name={
                    'ggH_shape':"T_{11}",
                    'gg0M_shape':'T_{33}',
                    'ggInt_13P_shape':'T_{13}',
                    'ggInt_13N_shape':'T_{13}',
                    'ggInt_13_shape':'T_{13}',
                    'gg0Ph_shape':'T_{22}',
                    'ggInt_12P_shape':'T_{12}',
                    'ggInt_12N_shape':'T_{12}',
                    'ggInt_12_shape':'T_{12}',
                    'qqZZ_shape':'T_{qqZZ}',
                    'TOTAL_PDF' : 'Total PDF'
                    }
                    
                if self.add_text_exist: 
                        data['setup']['add_text'][0]['text'] = "{0}(X,Y), {1}, GEN level, 8TeV".format(template_name[theTemplate], fs)
                        if theTemplate=="TOTAL_PDF":
                            data['setup']['add_text'][0]['text']+=", Parameters: k2/k1={0}, k3/k1={1}".format(self.k2k1_ratio,self.k3k1_ratio)
                        self.add_text(data['setup']['add_text'])
                
                
                plot_name = "{0}_{1}_{2}".format(theTemplFile.replace("/","_"), theTemplate, draw_opt)
                #plot_name = self.name
                self.save_extensions = ['png','pdf','eps']
                try:
                    data['setup']['save_ext']
                except KeyError:
                    self.log.info("No extensions are provided in setup. Using default: ", self.save_extensions)
                else:
                    self.save_extensions = list(data['setup']['save_ext'])
                self.save(self.c1, plot_name, self.save_extensions)
    
 
            
        
      
  def setCopyToWebDir(self,doCopy=False,webdir=""):
      if doCopy:
	  self.copy_to_web_dir = True
	  if webdir!="":
	      self.webdir = webdir
	  else:
	      raise ValueError, "You have to provide a webdir path if you want to copy the files."
      else:
	self.copy_to_web_dir = False
	self.webdir = ""
      return 0
	

  def doCopyToWebDir(self,plot_name, newname=""):
      if newname=="":
	    newname = plot_name
      if self.webdir!="" :
          misctools.make_sure_path_exists(self.webdir)
          if not os.path.exists("{0}/index.php".format(self.webdir)) : shutil.copy("/afs/cern.ch/user/r/roko/www/html/index.php",self.webdir)
	  shutil.copy(plot_name,self.webdir+"/"+newname)
	  self.log.info("Copied {0} to webdir {1}".format(plot_name,self.webdir+"/"+newname))
      else :
	  raise ValueError, "You have to provide a webdir path if you want to copy the files."      
      return 0

      
  def save(self, canv, plot_name, extensions=[]):     
      #extensions = ['.png','.pdf','.eps','.root']
      if len(extensions)==0: extensions=['']
      for ext in extensions:
	postfix = "."+ext
	if ext=='': postfix=''
	canv.SaveAs(plot_name+postfix)
	self.log.debug("Saving to: {0}.*".format(plot_name))
	if self.copy_to_web_dir : 
	  self.doCopyToWebDir(plot_name+postfix)
	  

      
      
    
    
    
    
if __name__ == "__main__":
     
    # parse the arguments and options
    #global opt, args
    #parseOptions()
    print "Starting plots..."
    
    pp = pprint.PrettyPrinter(indent=4)
    
    DEBUG = True
    cfg_reader = UniversalConfigParser(cfg_type="YAML",file_list = "plotTemplates.yaml")
    cfg_reader.setLogLevel(10) 
    
    plots_data = cfg_reader.get_dict()
    
    if DEBUG: 
	  print "plots_data = ", 
	  pp.pprint(plots_data)
    
    plotter = TemplatesPlotter()
    plotter.setCopyToWebDir(True,"/afs/cern.ch/user/r/roko/www/html/Geoloc/Templates/")
    plotter.set_k2k1_and_k3k1(3,0)
    for plot_name in plots_data['list_of_plots_to_do']:
	print 'Make plot: ',plot_name
	plotter.setName(plot_name)
	plotter.makePlot(plots_data[plot_name])
   
    
    
    