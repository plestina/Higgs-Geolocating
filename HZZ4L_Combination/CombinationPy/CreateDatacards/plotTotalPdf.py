#! /usr/bin/env python
from ROOT import *
import ROOT
import shutil
#import yaml
import pprint
from array import array
import lib.util.MiscTools as misctools
from lib.util.RootAttributeTranslator import RootAttributeTranslator 
from lib.util.Logger import *
from lib.util.UniversalConfigParser import UniversalConfigParser 

from lib.plotting.RootPlotters import PlotPolisher
from lib.plotting.RootPlottersBase import RootPlottersBase
import copy
import collections as cols

                        
class TotalPdfPlotter(PlotPolisher, RootPlottersBase):
  """Plots PDF from input templates without using nusance. 
     This is because, sometimes it happens that the templates
     are negative for some values of parameter. In that case
     the RooFit gives errors. 
  """
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
      
  def get_scale_factors(self):
            
        factors = {'lambda12_cosP' : 0.0,
            'lambda13_cosP' : 0.0,
            'lambda23_cosN' : 0.0,
            'lambda23_cosP' : 0.0,
            'gamma33' : 0.034,
            'gamma13' : 0.0,
            'gamma12' : -0.290992,
            'gamma11' : 1.0,
            'lambda13_cosN' : 0.0,
            'lambda12_sinN' : 0.0,
            'lambda12_cosN' : 0.581984,
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
            
        denominator  =   "1"
        #denominator_geoloc = "({0}+{1}*@0*@0+{2}*@1*@1+2*(({3}*@0)+({4}*@1)+({5}*@0*@1)))" ### added factor 2 for mixed terms
        ##denominator = denominator_geoloc.format(factors['gamma11'],factors['gamma22'],factors['gamma33'],
                                                ##factors['gamma12'],factors['gamma13'],factors['gamma23'])
        #denominator = "{0}+{1}+{2}".format(_getDenominatorSegment(1,denominator_geoloc),_getDenominatorSegment(2,denominator_geoloc),_getDenominatorSegment(3,denominator_geoloc))      
        #denominator = "0.354301097432*(1.0+0.0*@0*@0+0.034*@1*@1+(0.0*@0)+(0.0*@1)+(0.0*@0*@1))
                      #+0.184149324322*(1.0+0.0*@0*@0+0.034*@1*@1+(0.0*@0)+(0.0*@1)+(0.0*@0*@1))
                      #+0.461549578246*(1.0+0.0*@0*@0+0.040*@1*@1+(0.0*@0)+(0.0*@1)+(0.0*@0*@1))"
        norm={}
        for nom in nominator.keys():    
            norm[nom] = eval("({0})/({1})".format(nominator[nom], denominator))
            self.log.debug('Norm for {0}  = {1}'.format(nom, norm[nom]))    
            
        return norm

  def make_animated_gif(self, list_for_animated_gif, name_animated_gif, delay=100, loop=True):
        #import lib.util.MiscTools as misctools
        try:
            delay = int(delay)
        except ValueError:
            raise ValueError, "Delay must be an integer."
        
        loop_string = ""
        if loop:
            loop_string = " -loop 0"
            
        self.log.info('Creating animated gif: {0}'.format(name_animated_gif))    
        misctools.processCmd("convert -delay {2} {0}{3} {1} ".format(" ".join("%s"%(one_gif) for one_gif in list_for_animated_gif), name_animated_gif, delay, loop_string))    
        self.doCopyToWebDir(name_animated_gif)
        
  def makePlot(self, data):

        #self.c_big = ROOT.TCanvas("c_big","Template Components",1200,1200)
        #self.c_big.Divide(len())
        self.c1 =ROOT.TCanvas("cc1","Templates",1000,800)
        #pdf_components = {}
        self.c1.cd()
        ROOT.gPad.SetRightMargin(0.2)
        
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
        print data['setup']['add_text']
        try:
            data['setup']['add_text']
        except AttributeError:
            self.add_text_exist=False
        else:
            self.add_text_exist=True
        
        i=0
        for plot in data['content']:
            
                list_for_animated_gif=[]
                list_for_animated_gif_components=[]
                
                
                try:
                    plot['POI']
                except KeyError:
                    raise KeyError, "You chould provide POI dictionary with values in your confguration."
                else:
                    self.plot_total_pdf=True  
                    if isinstance(plot['POI'], dict):
                        if "k2k1_ratio" in plot['POI'].keys():
                            k2k1_list = list(plot['POI']['k2k1_ratio'])
                            self.log.debug('k2k1_list = {0}'.format(", ".join([str(k) for k in k2k1_list])))
                        if "k3k1_ratio" in plot['POI'].keys():
                            k3k1_list = list(plot['POI']['k3k1_ratio'])
                            self.log.debug('k3k1_list = {0}'.format(", ".join([str(k) for k in k3k1_list])))
                
                for k2k1 in k2k1_list:
                    for k3k1 in k3k1_list:
                        self.set_k2k1_and_k3k1(k2k1, k3k1)
                        self.log.debug('Seting parameters to k2k1 = {0} and k3k1 = {1}'.format(self.k2k1_ratio,self.k3k1_ratio) )
                        
                        n_components = len(plot['histos'])+1 #we include the total pdf into computation
                        self.c_big = ROOT.TCanvas("c_big","Template Components",int(500*n_components),400)
                        #self.c_big = ROOT.TCanvas("c_big","Template Components",1600,800)
                        self.c_big.Divide(n_components,1)
                        
                        
                        self.log.debug("Canvas is divided to N = {0}".format(n_components))
                        pdf_components = cols.OrderedDict()
                        pdf_components_minmax = cols.OrderedDict()
                        #pdf_components_order = []
                        
                        for fs in data['setup']['final_states']:
                            print 200*"_"
                            theTemplFile = plot['file']+"_"+fs+".root"
                            
                            if DEBUG:
                                pp.pprint(plot)
                        
                            f = TFile(theTemplFile,"READ")
                            self.log.debug("Running on file: {0}".format(theTemplFile))
                            norm= self.get_scale_factors()
                            
                            for i_pdf in range(len(plot['histos'])):
                                if i_pdf==0:
                                    total_pdf = f.Get(plot['histos'][i_pdf])
                                    total_pdf.Scale(norm[plot['histos'][i_pdf]])
                                    pdf_components[plot['histos'][i_pdf]] = copy.deepcopy(total_pdf)
                                    self.log.debug('TOTAL_PDF -> Picking up the first term {0} and the scale {1}. The histo is TH2: {2}'.format(plot['histos'][i_pdf], norm[plot['histos'][i_pdf]], isinstance(total_pdf,TH2)))
                                    #pdf_components_order.append(plot['histos'][i_pdf])

                                another_template = f.Get(plot['histos'][i_pdf])
                                another_template.Scale(float(norm[plot['histos'][i_pdf]])) 
                                pdf_components[plot['histos'][i_pdf]] = copy.deepcopy(another_template)
                                #pdf_components_order.append(plot['histos'][i_pdf])
                                
                                self.log.debug('TOTAL_PDF -> Adding up term {0} and the scale {1} to histo TH2:{2}'.format(plot['histos'][i_pdf], norm[plot['histos'][i_pdf]],isinstance(total_pdf,TH2)))
                                total_pdf.Add(another_template)
                                
                            self.total_pdf = total_pdf    
                            self.log.debug('TOTAL_PDF -> Added all terms and now the histo is TH2:{0}'.format(isinstance(self.total_pdf,TH2)))
                            bins_x, bins_y = self.total_pdf.GetXaxis().GetNbins(), self.total_pdf.GetYaxis().GetNbins()
                            self.log.debug('Size of template (x,y) = ({0},{1})'.format(bins_x, bins_y))
                            
                            pdf_components["TotalPDF"] = copy.deepcopy(self.total_pdf)

                            
                            for x_bin in range(1,bins_x+1):
                                for y_bin in range(1,bins_y+1):
                                    bin_content = self.total_pdf.GetBinContent(x_bin, y_bin)
                                    #self.log.debug("Bin (i,j) = {0},{1} value = {2}".format(x_bin,y_bin, bin_content))
                                    if bin_content <= 0:
                                        raise ValueError, "!!! PDF CANNOT BE 0 or NEGATIVE !!!"
                                    assert bin_content > 0, "!!! PDF CANNOT BE 0 or NEGATIVE !!! Check your tamplates. Bin (i,j) = {0},{1} value = {2}".format(x_bin,y_bin, bin_content)   
                            
                            ### settup for the plot: text, legend and save    
                            self.c1.cd()
                            if self.setup_exist: 
                                self.arrangeAxis(self.total_pdf,data['setup'])
                            if self.leg_exist:
                                self.leg.AddEntry(self.total_pdf,plot['legend']['text'],plot['legend']['opt']);
                            try:
                                plot['draw_opt']
                            except ValueError:
                                draw_opt = "COLZ"
                            else:
                                draw_opt = str(plot['draw_opt'])
                            if draw_opt.lower()=="surf":
                                self.c1.SetTheta(16.56)
                                self.c1.SetPhi(57.83133)
                            self.total_pdf.SetMaximum(0.033)
                            self.total_pdf.SetMinimum(-0.005)
                            self.total_pdf.Draw(draw_opt)
                            
                            if self.leg_exist: self.leg.Draw()
                            
                                            
                            if self.add_text_exist: 
                                    data['setup']['add_text'][0]['text'] = "Total PDF(X,Y), {0}, GEN level, 8TeV".format(fs)
                                    data['setup']['add_text'][0]['text']+=", k2/k1 ={0: >6.2f}, k3/k1 ={1: >6.2f}".format(self.k2k1_ratio,self.k3k1_ratio)
                                    self.add_text(data['setup']['add_text'])
                            
                            ### saving the plot and copyng to webdir(if configured)
                            plot_name = "TotalPdf_{0}_{1}_k2k1_{2}_k3k1_{3}".format(theTemplFile.replace("/","_"), draw_opt,self.k2k1_ratio,self.k3k1_ratio)
                            #plot_name = self.name
                            self.save_extensions = ['png','pdf','eps', 'gif']
                            try:
                                data['setup']['save_ext']
                            except KeyError:
                                self.log.info("No extensions are provided in setup. Using default: ", self.save_extensions)
                            else:
                                self.save_extensions = list(data['setup']['save_ext'])
                            self.save(self.c1, plot_name, self.save_extensions)
                            list_for_animated_gif.append(plot_name+".gif")
                            
                            #### Plot in the same canvas
                            i_div = 1
                            print "COMPONENTS :", pdf_components.keys()
                            
                            for comp in pdf_components.keys():
                                    self.c_big.cd(i_div)
                                    ROOT.gPad.SetRightMargin(0.2)
                                    ROOT.gPad.SetLeftMargin(0.1)
                                    
                                    i_div+=1
                                    ### settup for the plot: text, legend and save    
                                    if self.setup_exist: 
                                        self.arrangeAxis(pdf_components[comp],data['setup'])
                                    if self.leg_exist:
                                        self.leg.AddEntry(pdf_components[comp],plot['legend']['text'],plot['legend']['opt']);
                                    try:
                                        plot['draw_opt']
                                    except ValueError:
                                        draw_opt = "COLZ"
                                    else:
                                        draw_opt = str(plot['draw_opt'])
                                    if draw_opt.lower()=="surf":
                                        self.c_big.SetTheta(16.56)
                                        self.c_big.SetPhi(57.83133)
                                    pdf_components[comp].SetMaximum(0.033)
                                    pdf_components[comp].SetMinimum(-0.033)
                                    
                                    pdf_components[comp].GetYaxis().SetTitleOffset(0.7)
                                    if (i_div-1)>1:
                                        pdf_components[comp].GetYaxis().SetLabelOffset(999)
                                        pdf_components[comp].GetYaxis().SetLabelSize(0)
                                        pdf_components[comp].GetYaxis().SetTitle("")
                                    
                                    pdf_components[comp].Draw(draw_opt)
                                    
                                    if self.leg_exist: self.leg.Draw()
                                    
                                    if self.add_text_exist: 
                                            data['setup']['add_text'][0]['text'] = "{1}(X,Y), {0}, GEN level, 8TeV".format(fs,comp )
                                            data['setup']['add_text'][0]['text']+=", k2/k1 ={0: >6.2f}, k3/k1 ={1: >6.2f}".format(self.k2k1_ratio,self.k3k1_ratio)
                                            self.add_text(data['setup']['add_text'])
                                    
                            ### saving the plot and copyng to webdir(if configured)
                            plot_name = "TotalPdfComponents_{0}_{1}_k2k1_{2}_k3k1_{3}".format(theTemplFile.replace("/","_"), draw_opt,self.k2k1_ratio,self.k3k1_ratio)
                            #plot_name = self.name
                            self.save_extensions = ['png','pdf','eps', 'gif']
                            try:
                                data['setup']['save_ext']
                            except KeyError:
                                self.log.info("No extensions are provided in setup. Using default: ", self.save_extensions)
                            else:
                                self.save_extensions = list(data['setup']['save_ext'])
                            self.save(self.c_big, plot_name, self.save_extensions)
                            list_for_animated_gif_components.append(plot_name+".gif")
                            
                            
                            
                name_animated_gif = "TotalPdf_{0}_{1}_ANIMATED_k2k1_{2}_k3k1_{3}.gif".format(theTemplFile.replace("/","_"),draw_opt,"_".join("{0}".format(one_k2k1) for one_k2k1 in k2k1_list), "_".join("{0}".format(one_k3k1) for one_k3k1 in k3k1_list))                            
                self.make_animated_gif(list_for_animated_gif, name_animated_gif)      
                name_animated_gif = "TotalPdfComponents_{0}_{1}_ANIMATED_k2k1_{2}_k3k1_{3}.gif".format(theTemplFile.replace("/","_"),draw_opt,"_".join("{0}".format(one_k2k1) for one_k2k1 in k2k1_list), "_".join("{0}".format(one_k3k1) for one_k3k1 in k3k1_list))                            
                self.make_animated_gif(list_for_animated_gif_components, name_animated_gif)      
                            
                            
                ##make animated gif            
                #with open("animated_gif.list","w") as f_list_for_animated_gif:
                    ##list_for_animated_gif.write()
                    #f_list_for_animated_gif.write("\n".join("%s"%(one_gif) for one_gif in list_for_animated_gif))
                    
                #name_animated_gif = "TotalPdf_{0}_{1}_ANIMATED_k2k1_{2}_k3k1_{3}.gif".format(theTemplFile.replace("/","_"),draw_opt,"_".join("{0}".format(one_k2k1) for one_k2k1 in k2k1_list), "_".join("{0}".format(one_k3k1) for one_k3k1 in k3k1_list))
                
                ##import lib.util.MiscTools as misc
                ##misc.processCmd('echo PROCESSING COMMAND')                                 
                #misctools.processCmd("convert -delay 100 {0} -loop 0 {1} ".format(" ".join("%s"%(one_gif) for one_gif in list_for_animated_gif), name_animated_gif))    
                #self.doCopyToWebDir(name_animated_gif)

                
                
                
if __name__ == "__main__":
     
    # parse the arguments and options
    #global opt, args
    #parseOptions()
    print "Starting plots..."
    
    pp = pprint.PrettyPrinter(indent=4)
    
    DEBUG = True
    cfg_reader = UniversalConfigParser(cfg_type="YAML",file_list = "plotTotalPdf.yaml")
    cfg_reader.setLogLevel(10) 
    
    plots_data = cfg_reader.get_dict()
    
    if DEBUG: 
	  print "plots_data = ", 
	  pp.pprint(plots_data)
    
    plotter = TotalPdfPlotter()
    plotter.setCopyToWebDir(True,"/afs/cern.ch/user/r/roko/www/html/Geoloc/TotalPdfs/")
    
    for plot_name in plots_data['list_of_plots_to_do']:
	print 'Make plot: ',plot_name
	plotter.setName(plot_name)
	plotter.makePlot(plots_data[plot_name])
   
    
    
    