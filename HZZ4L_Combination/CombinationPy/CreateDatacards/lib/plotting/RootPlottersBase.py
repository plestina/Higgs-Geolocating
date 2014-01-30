#! /usr/bin/env python
from ROOT import *
import ROOT
from array import array
import shutil, os
import pprint
from lib.util.RootAttributeTranslator import *
from lib.util.Logger import *
import lib.util.MiscTools as misctools

class RootPlottersBase(object):  
  """Class as a base class for many plotters containing the structure, common functions...
  """
  
  def __init__(self,name = "plotters_base_functionality" ):
        self.log = Logger().getLogger(self.__class__.__name__, 10)
        self.name = name
        #ROOT.gSystem.AddIncludePath("-I$ROOFITSYS/include/");
        #ROOT.gROOT.ProcessLine(".L tdrstyle.cc")
        #from ROOT import setTDRStyle
        #ROOT.setTDRStyle(True)

        #ROOT.gStyle.SetPalette(1)
        #ROOT.gStyle.SetOptStat(0)
        self.copy_to_web_dir = False
        self.webdir = ""
        self.save_extensions = ['png','pdf','eps']
        self.pp = pprint.PrettyPrinter(indent=4)
      
  def setName(self, newname): self.name = newname
        
  def makePlot(self, data):
        print "This is a default method for plotters. It has to be implemented in derived classes"
        pass
    
  def   get_TTree(tree_name, file_name):
        rootfile = ROOT.TFile.Open(root_file_name,'READ')
        t = rootfile.Get('limit')
        return t
        
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
        if len(extensions)==0: 
            extensions=['']
        for ext in extensions:
            postfix = "."+ext
            if ext=='': postfix=''
            canv.SaveAs(plot_name+postfix)
            self.log.debug("Saving to: {0}.*".format(plot_name))
            if self.copy_to_web_dir : 
                self.doCopyToWebDir(plot_name+postfix)
            

      
      
    
    
    
    
    