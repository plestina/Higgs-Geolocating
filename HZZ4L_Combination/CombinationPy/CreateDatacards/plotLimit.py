from ROOT import *
import ROOT
from array import array
import shutil
#import yaml
import pprint
#from lib.util.RootAttributeTranslator import RootAttributeTranslator 
from lib.util.Logger import Logger
from lib.util.UniversalConfigParser import UniversalConfigParser
from lib.plotting.RootPlotters import PlotPolisher,SimplePlotter
from lib.plotting.RootPlottersBase import RootPlottersBase
import yaml
import os


class LikelihoodFitContourPlot(PlotPolisher,RootPlottersBase):  

    def __init__(self,name = "limit_plot" ):
        
        self.log = Logger().getLogger(self.__class__.__name__, 10)
        self.name = name
        #ROOT.gROOT.ProcessLine(".L /afs/cern.ch/work/r/roko/Stat/CMSSW_611_JCP/src/HZZ4L_Combination/CombinationPy/CreateDatacards/tdrstyle.cc")
        #from ROOT import setTDRStyle
        #ROOT.setTDRStyle(True)
        ROOT.gROOT.SetBatch(1)
        ROOT.gStyle.SetPalette(1)
        ROOT.gStyle.SetOptStat(0)
        #self.copy_to_web_dir = False
        #self.webdir = ""
        self.pp = pprint.PrettyPrinter(indent=4)
        self.title = ''
        self.canv = {}
        self.axis_nice_names = {
                                'k2k1_ratio': ('k_{2}/k_{1}', 'AxisK2K1'),
                                'k3k1_ratio': ('k_{3}/k_{1}', 'AxisK3K1'),
                                }
        
        
    def set_title(self,title):
        self.title = title
        self.log.debug('Title of plot: {0}'.format(self.title))
        
    def get_limits_dict(self):
        try:
            self.limits_dict
        except:
            raise RuntimeError, 'Please call the make_plot function before...'
        return self.limits_dict
    def _get_limit_intervals(self,ll_list, ul_list, x_min=-99999, x_max=99999):
        self.log.debug('Getting limit intervals...')
        all_tuple_list_sorted = [(element, 'll') for element in ll_list] + [(element, 'ul') for element in ul_list]
        all_tuple_list_sorted.sort(key=lambda tup: tup[0])
        #print all_tuple_list_sorted 
        if len(all_tuple_list_sorted ) > 0:
            if all_tuple_list_sorted[0][1]=='ul':
                all_tuple_list_sorted.insert(0, (x_min,'ll'))
                #print 'adding in front ', all_tuple_list_sorted 
            if all_tuple_list_sorted[-1][1]=='ll':    
                all_tuple_list_sorted.append((x_max,'ul'))
                #print 'adding to back ', all_tuple_list_sorted
        assert len(all_tuple_list_sorted)%2==0, 'Should be an even number of limits : {0} '.format(all_tuple_list_sorted)
        ind = 0
        llul_pairs = []
            
        while ind < len(all_tuple_list_sorted):
            a=all_tuple_list_sorted[ind][1]
            if (ind+1) == len(all_tuple_list_sorted):
                self.log.warn('The last limit provided ({0}) will be ignored.'.format(all_tuple_list_sorted[ind][0]))
                break
            b=all_tuple_list_sorted[ind+1][1]
            if (a,b) == ('ll','ul'): 
                #assert all_tuple_list_sorted[ind][1]=='ll'  and all_tuple_list_sorted[ind+1][1]=='ul', 'Is this physical? You should have paires in (ll, ul) ordering.'
                llul_pairs.append((all_tuple_list_sorted[ind][0], all_tuple_list_sorted[ind+1][0]))
                ind+=2
            else:
                self.log.warn('It seems strange that you have 2 limits of same type consequtively. Ignoring limit at {0}.'.format(all_tuple_list_sorted[ind][0]))
                ind+=1    
        self.log.info('Made list of limit pairs: {0}'.format(llul_pairs))
        return llul_pairs

    def putLimitLinesOnPlot(self, fit_result_dict, c):
        #   #1 sigma and 2 sigma lines
        
        self.log.debug('Putting limit lines on the plot.')
        self.pp.pprint(fit_result_dict)
        c.cd()
        bestFit = fit_result_dict['BF']
        bf_arr = TArrow()
        bf_arr.SetLineColor(kRed)
        bf_arr.SetLineWidth(2)
        bf_arr.DrawArrow(bestFit,-0.5,bestFit,-0.25,0.02,"|->")

        gra = c.GetPrimitive('Graph').GetXaxis()
        x_min = gra.GetXmin()
        x_max = gra.GetXmax()

        line = TLine()
        line.SetLineStyle(kDashed)
        #line.SetLineStyle(kSolid)
        line.SetLineWidth(2)
        #  #68% C.L        
        deltaNLL = 1
        line.SetLineColor(kRed)
        llul_pairs = self._get_limit_intervals(fit_result_dict['LL68'], fit_result_dict['UL68'],x_min, x_max )
        for pair in llul_pairs:
            line.DrawLine(pair[0],deltaNLL,pair[1],deltaNLL)
            line.DrawLine(pair[0],0,pair[0],deltaNLL)
            line.DrawLine(pair[1],0,pair[1],deltaNLL)
        
        #95 %CL
        line.SetLineColor(kBlue)
        deltaNLL = 3.84
        llul_pairs = self._get_limit_intervals(fit_result_dict['LL95'], fit_result_dict['UL95'],x_min, x_max )
        for pair in llul_pairs:
            line.DrawLine(pair[0],deltaNLL,pair[1],deltaNLL)
            line.DrawLine(pair[0],0,pair[0],deltaNLL)
            line.DrawLine(pair[1],0,pair[1],deltaNLL)

          
    def putBrasilianFlagsOnPlot(self, fit_result_dict, c):
        #   #1 sigma and 2 sigma lines
        
        self.log.debug('Putting limit lines on the plot.')
        self.pp.pprint(fit_result_dict)
        c.cd()
        bestFit = fit_result_dict['BF']
        bf_arr = TArrow()
        bf_arr.SetLineColor(kRed)
        bf_arr.SetLineWidth(2)
        bf_arr.DrawArrow(bestFit,-0.5,bestFit,-0.25,0.02,"|->")

        gra_x = c.GetPrimitive('Graph').GetXaxis()
        x_min = gra_x.GetXmin()
        x_max = gra_x.GetXmax()
        gra_y = c.GetPrimitive('Graph').GetYaxis()
        #y_max = gra_y.GetXmax()
        y_max = 5
        
        
        box = TBox()
        
        #95 %CL
        box.SetFillStyle(3002)
        box.SetLineColor(kBlue)
        box.SetFillColor(kYellow)
        llul_pairs = self._get_limit_intervals(fit_result_dict['LL95'], fit_result_dict['UL95'],x_min, x_max )
        for pair in llul_pairs:
            box.DrawBox(pair[0],0,pair[1],y_max)
            
        box.SetFillColor(kGreen)
        llul_pairs = self._get_limit_intervals(fit_result_dict['LL68'], fit_result_dict['UL68'],x_min, x_max )
        for pair in llul_pairs:
            box.DrawBox(pair[0],0,pair[1],y_max)
            
        gr = c.GetPrimitive('Graph')
        gr.Draw('Psame')
        


    def putInformationOnPlot(self, fit_result_dict, c, x, y, POI="|k_3/k_1|"):

        latex = TLatex()
        latex.SetTextSize(0.025)
        latex.SetTextAlign(10)  #align at special bottom

        max_y = 5.0
        max_y_2 = max_y - 6*(max_y)/100.0
        max_y_3 = max_y - 10*(max_y)/100.0
        max_y_4 = max_y - 14*(max_y)/100.0
        self.log.debug('Positions max_y: {0} {1} {2} {3}'.format(max_y, max_y_2, max_y_3, max_y_4))


        min_x = x.GetXmin() + 2*(x.GetXmax() - x.GetXmin())/100
        self.log.debug("min_x = {0} max_y = {1}".format(min_x, max_y))
        fit_info = '{0}'.format(self.title)

        import string
        ll68 = string.join(['{0:.2f}'.format(limit) for limit in fit_result_dict['LL68']],',')
        ll95 = string.join(['{0:.2f}'.format(limit) for limit in fit_result_dict['LL95']],',')
        ul68 = string.join(['{0:.2f}'.format(limit) for limit in fit_result_dict['UL68']],',')
        ul95 = string.join(['{0:.2f}'.format(limit) for limit in fit_result_dict['UL95']],',')
        
        
        self.log.debug('ll68: {0}'.format(ll68))
        self.log.debug('ll95: {0}'.format(ll95))
        self.log.debug('ul68: {0}'.format(ul68))
        self.log.debug('ul95: {0}'.format(ul95))
                      
        if len(fit_result_dict['LL68'])==0 :  ll68 = 'N.A.'
        if len(fit_result_dict['LL95'])==0 :  ll95 = 'N.A.'
        if len(fit_result_dict['UL68'])==0 :  ul68 = 'N.A.'
        if len(fit_result_dict['UL95'])==0 :  ul95 = 'N.A.'
        
        fit_info_2 = "Best fit {0} = {1: .3f}".format(POI, fit_result_dict['BF'])
        fit_info_3 = "U.L. 68(95)% = {0}({1})".format(ul68, ul95)
        fit_info_4 = "L.L. 68(95)% = {0}({1})".format(ll68, ll95)

        print  fit_info  

        c.cd()
        latex.DrawLatex(min_x, max_y,fit_info)
        latex.DrawLatex(min_x, max_y_2,fit_info_2)
        latex.DrawLatex(min_x, max_y_3,fit_info_3)
        latex.DrawLatex(min_x, max_y_4,fit_info_4)

  
   
   
    def _setup_graph(self, graph, POI, dims=1):
        assert dims==1, 'Currently we have only one POI axis supported.Sorry :('
        
        graph.SetMarkerSize(2.0)
        #   TString ytitle = "|k_3/k_1|"
        graph.SetNameTitle('Graph','')
        graph.GetXaxis().SetTitle(self.axis_nice_names[POI][0])
        graph.GetYaxis().SetTitle("-2 #Delta ln L")
        graph.GetXaxis().SetLabelSize(0.04)
        graph.GetYaxis().SetRangeUser(0,5)
        graph.SetMarkerStyle(20)
        graph.SetMarkerSize(.5)

        graph.GetYaxis().SetTitleOffset(1.3)
        graph.GetXaxis().SetTitleOffset(1.3)
        


    def make_plot(self, file_name, POI,title=''):
        if title:
            self.set_title(title)
        
        self.fit_results = FitResultReader(POI, [file_name], combine_method='MultiDimFit')
        self.fit_results.set_files([file_name])
        self.fit_results.set_POI(POI)

        gr_xnll = self.fit_results.get_graph('{0}:2*deltaNLL'.format(POI))
        self._setup_graph(gr_xnll, POI, dims=1)
        #ROOT.gStyle.SetPalette(1)
        #   gr_xnll.SetMarkerStyle(34) 
        c1 = TCanvas("cc1","CANVAS1",600,600)
        c1.cd()
        #c1.SetLogy()
        gPad.SetRightMargin(0.0085)

        self.limits_dict = self.fit_results.get_results_dict(POI, option='standard')
        
        with open('{0}.limits'.format(file_name), 'w') as limits_file:
             #limits_file.write(yaml.dump(self.limits_dict, default_flow_style=True))
             for value_id in self.limits_dict.keys():
                 limits_file.write('{0} : {1}\n'.format(value_id, self.limits_dict[value_id]))
                 
             self.log.info('Writing fit results into {0}'.format(limits_file.name))

        #   if (POI != "k3k1_ratio") 
        #     gr_xnll.Draw("AP")
        #     c1.SaveAs(Form("%s.png",file_name.Data()))
        #   

        #
        gr_xnll.Draw("AP")
        self.putBrasilianFlagsOnPlot(self.limits_dict, c1)
        self.putLimitLinesOnPlot(self.limits_dict, c1)
        self.putInformationOnPlot(self.limits_dict,c1, gr_xnll.GetXaxis(), gr_xnll.GetYaxis(), self.axis_nice_names[POI][0])
        #gr_xnll.Draw('same')
        
         #plot_name = "limitVsLumi_2D_k3k1_0_logy"
        plot_name = self.name
        self.save_extensions = ['png','gif', 'root']
        self.save(c1, '{0}.{1}'.format(file_name, self.axis_nice_names[POI][1]), self.save_extensions)

        
        #TCanvas *c2=new TCanvas("cc2","CANVAS2",600,600)
        #c2.cd()
        #gPad.SetRightMargin(0.0085)

        #TGraph * gr_xnll_arctan = getRescaledGraph(gr_xnll,0.03676929746947648)
        #gr_xnll_arctan.GetXaxis().SetTitle("1/#pi arctan(#sqrt#gamma_33 k_3/k_1)")
        #gr_xnll_arctan.GetYaxis().SetTitle("-2 #Delta ln L")
        #gr_xnll_arctan.GetYaxis().SetTitleOffset(1.3)
        #gr_xnll_arctan.GetXaxis().SetTitleOffset(1.3)
        #gr_xnll_arctan.SetTitle("")
        #gr_xnll_arctan.GetXaxis().SetLabelSize(0.04)
        #gr_xnll_arctan.SetMarkerStyle(20)
        #gr_xnll_arctan.SetMarkerSize(.5)
        #gr_xnll_arctan.GetYaxis().SetRangeUser(0,5)
        #gr_xnll_arctan.Draw("AP")
        ##     max_y = gr_xnll_arctan.GetYaxis().GetXmax()
        ##     min_x = gr_xnll_arctan.GetXaxis().GetXmin()
        ##     latex.DrawLatex(min_x, max_y,fit_info)
        ##     putInformationOnPlot(t, POI, c2, gr_xnll_arctan.GetXaxis(), gr_xnll_arctan.GetYaxis())
        #putInformationOnPlot(fitres, c2, gr_xnll_arctan.GetXaxis(), gr_xnll_arctan.GetYaxis())
        #putLimitLinesOnPlot(&getRescaledFitResult(fitres,0.03676929746947648), c2)

        #c2.SaveAs(Form("%s.AxisArctanGamK3K1.png",file_name.Data()))




        #TCanvas *c3=new TCanvas("cc3","CANVAS3",600,600)
        #c3.cd()
        #gPad.SetRightMargin(0.0085)

        #TGraph * gr_xnll_arctan_2 = getRescaledGraph(gr_xnll,1)
        #gr_xnll_arctan_2.GetXaxis().SetTitle("1/#pi arctan(k_3/k_1)")
        #gr_xnll_arctan_2.GetYaxis().SetTitle("-2 #Delta ln L")
        #gr_xnll_arctan_2.GetYaxis().SetTitleOffset(1.3)
        #gr_xnll_arctan_2.GetXaxis().SetTitleOffset(1.3)
        #gr_xnll_arctan_2.SetTitle("")
        #gr_xnll_arctan_2.GetXaxis().SetLabelSize(0.04)
        #gr_xnll_arctan_2.SetMarkerStyle(20)
        #gr_xnll_arctan_2.SetMarkerSize(.5)
        #gr_xnll_arctan_2.GetYaxis().SetRangeUser(0,5)
        #gr_xnll_arctan_2.Draw("AP")
        #putInformationOnPlot(fitres, c3, gr_xnll_arctan_2.GetXaxis(), gr_xnll_arctan_2.GetYaxis())
        #putLimitLinesOnPlot(&getRescaledFitResult(fitres,1), c3)
        #c3.SaveAs(Form("%s.AxisArctanK3K1.png",file_name.Data()))


        ##fa3 result
        #TCanvas *c4=new TCanvas("cc4","CANVAS4",600,600)
        #c4.cd()
        #gPad.SetRightMargin(0.0085)

        #TGraph * gr_xnll_fa3 = getfa3RescaledGraph(gr_xnll)
        #gr_xnll_fa3.GetXaxis().SetTitle("f_a3")
        #gr_xnll_fa3.GetYaxis().SetTitle("-2 #Delta ln L")
        #gr_xnll_fa3.GetYaxis().SetTitleOffset(1.3)
        #gr_xnll_fa3.GetXaxis().SetTitleOffset(1.3)
        #gr_xnll_fa3.SetTitle("")
        #gr_xnll_fa3.GetXaxis().SetLabelSize(0.04)
        #gr_xnll_fa3.SetMarkerStyle(20)
        #gr_xnll_fa3.SetMarkerSize(.5)
        #gr_xnll_fa3.GetYaxis().SetRangeUser(0,5)
        #gr_xnll_fa3.Draw("AP")
        ##     max_y = gr_xnll_fa3.GetYaxis().GetXmax()
        ## #     min_x = gr_xnll_fa3.GetXaxis().GetXmin()
        ##     latex.DrawLatex(min_x, max_y,fit_info)
        ##     putInformationOnPlot(t, POI, c4, gr_xnll_fa3.GetXaxis(), gr_xnll_fa3.GetYaxis())
        #putInformationOnPlot(getfa3RescaledFitResult(fitres), c4, gr_xnll_fa3.GetXaxis(), gr_xnll_fa3.GetYaxis(),"f_a3")
        #putLimitLinesOnPlot(&getfa3RescaledFitResult(fitres), c4)
        #c4.SaveAs(Form("%s.AxisFa3.png",file_name.Data()))



        #if (POI == "k2k1_ratio") 

        #gr_xnll.Draw("AP")
        #putInformationOnPlot(fitres,c1, gr_xnll.GetXaxis(), gr_xnll.GetYaxis())
        #putLimitLinesOnPlot(&fitres, c1)
        #c1.SaveAs(Form("%s.AxisK2K1.png",file_name.Data()))
        #c1.SaveAs(Form("%s.AxisK2K1.root",file_name.Data()))


  
  






def return_filenames(rootdir,searchstring) :
  fileslist = []
  for folders,subs,files in os.walk(rootdir):
      for fi in files:
          fullpath = os.path.join(folders,fi)
          if ".root" in fi: fileslist.append(fullpath)
  if searchstring: fileslist=filter(lambda x: searchstring in x,fileslist)
  return fileslist

def AreSame(a, b, tolerance=None):
    if tolerance==None:
        import sys
        tolerance = sys.float_info.epsilon
    return abs(a - b) < tolerance
    
def belongsTo(value, rangeStart, rangeEnd):
    if value >= rangeStart and value <= rangeEnd:
        return True
    return False
    
    
    
  
class FitResultReader(object):
    """Reads the tree with fit information and gives back the 
       information relevant for plotng the limits or measurements. 
    """
    def __init__(self, POIs=None, file_names=[], combine_method=None):
        self.log = Logger().getLogger(self.__class__.__name__, 10)
       
        self.combine_method = combine_method
        self.combine_result = None
        self.global_best_fit = {}  #best fit from all the files. 
        self.global_best_fit_dict= {}  #contains one best fit per input file
        self.ll_values_dict = {}
        self.ul_values_dict = {}
        self.contours = None
        self.contour_graph= {}
        self._has_parsed_combine_result_already = False
        self.set_POI(POIs)
        self.set_files(file_names)
        
    def set_POI(self, POIs):
        """Provide a list of POIs in terms of python list 
           of strings or string with POIs separated by ";:*, "
        """
        self.POI = []
        assert isinstance(POIs, list) or isinstance(POIs, str), "POIs should be provided either as list of strings or as string with \";,: \" as delimiters. "
        if isinstance(POIs, list):
            self.POI = POIs
        elif isinstance(POIs, str):
            import re
            POIs = re.sub('[;,: ]+',':',POIs) #pois can be split by ";:*, " - we don't care
            self.POI = POIs.split(":")
        
        for poi in self.POI:
            self.global_best_fit.update({poi : None})
            self.log.debug('Initializing the best fit dictionary {0}'.format(self.global_best_fit))
            
        self.log.debug('Setting POI list to {0}'.format(self.POI))
        
    def set_files(self, file_names, start_dir = "."):
        """Set the list of output files from combine that will be used. 
           One can use the start dir and filename pattern to run on all
           files that are found recursively on the path.
        """
        self.file_list = []
        assert isinstance(file_names, list) or isinstance(file_names, str), "File names should be provided either as list of strings or as string with \";,: \" as delimiters. "
        assert isinstance(start_dir, str), "The start directory should be provided as string."
        if isinstance(file_names, list):
            self.file_list = file_names
        elif isinstance(file_names, str):
            raise ValueError, 'Please provide a list of strings for the file names. The option with strings only doesn\'t work for the moment. :('
            self.file_list = return_filenames(start_dir, self.file_list)
            
        print 'File list = ', self.file_list    
        self.log.debug('Loaded {0} files.'.format(len(self.file_list)))
        self._has_parsed_combine_result_already = False  #has to be set to False so that the limits and best fits are recalculated when new file is set. 

    def _get_crossings_for_limits(self, list_of_segments, cl=0.68):
        """Internal function for getting the Y values from given
           TGraph objets. It is used to get POI limits for particular 
           confidence level.
        """
        assert belongsTo(float(cl),0,1), "Confidence level has to be given in interval [0,1]"
        quantileExpected = 1 - cl
        values=[]
        for seg in list_of_segments:
            #ll_seg is a TGraph
            xmin = TMath.MinElement(seg.GetN(),seg.GetX())
            xmax = TMath.MaxElement(seg.GetN(),seg.GetX())
            if belongsTo(quantileExpected, xmin, xmax):
                values.append(seg.Eval(quantileExpected))
        return values
        
    def get_graph(self, contour_axis="x:y:z", dims = 1):
        """Returns the full likelihood scan graph.
           Specify contour_axis= "2*deltaNLL" or "1-quantileExpected"
           The last axis provided should be 
        """
        import re
        contour_axis = re.sub('[;,:]+',':',contour_axis) #can be split by ";:*, " - we don't care
        try:
            self.contour_graph[contour_axis]
        except KeyError:
            contour_axis_list = contour_axis.split(":")
            dims = len(contour_axis_list)
            assert 1<dims<=3, "We can accept 2 to 3 axis for the graph. You provided {0}.Please behave :)".format(dims)
            assert ('deltaNLL' in contour_axis_list[-1] or 'quantileExpected' in contour_axis_list[-1]), 'Your last axis has to contain either deltaNLL or quantileExpected.'
            import copy
            required_branches = copy.deepcopy(contour_axis_list)
            #Solve to accept even a formula as the axis.
            if 'deltaNLL' in contour_axis_list[-1]:
                #self.log.debug('deltaNLL in contour_axis_list: {0}'.format(contour_axis_list[-1]) )
                contour_axis_list[-1] = str(contour_axis_list[-1]).replace('deltaNLL','t.deltaNLL')
                required_branches[-1] = 'deltaNLL'
                #self.log.debug('deltaNLL is an estimator: {0}'.format(contour_axis_list[-1]) )
                
            elif 'quantileExpected' in contour_axis_list[-1]:
                contour_axis_list[-1] = contour_axis_list[-1].replace('quantileExpected', 't.quantileExpected')
                required_branches[-1] = 'quantileExpected'
                #self.log.debug('quantileExpected is an estimator: {0}'.format(contour_axis_list[-1]) )
                
            self.log.debug('Changing names of pois for evaluation of formula later: N_poi = {0} N_axis= {1}'.format(len(self.POI),len(contour_axis_list[:-1])))
            for poi_id in range(len(self.POI)):
                for axis_id in range(len(contour_axis_list[:-1])):
                    self.log.debug('Changing names of pois for evaluation of formula later: poi_id = {0} axis_id = {1}'.format(poi_id, axis_id))
                    contour_axis_list[axis_id] = contour_axis_list[axis_id].replace(self.POI[poi_id],'t.{0}'.format(self.POI[poi_id]))
                    required_branches[axis_id] = self.POI[poi_id]
                    
            self.log.debug('Contour axis list changed for evaluation of formula to {0}'.format(contour_axis_list))
                    
            
            if dims==2:
                self.contour_graph[contour_axis] = TGraph()
            elif dims==3:
                self.contour_graph[contour_axis] = TGraph2D()
            self.contour_graph[contour_axis].SetNameTitle(contour_axis,contour_axis.replace(':',';') )
            
                       
            self.log.debug('Graph {0} is being created from the tree in {1}'.format(contour_axis, self.file_list[0]))
            
            rootfile = TFile.Open(self.file_list[0],'READ')
            if not rootfile:
                raise IOError, 'The file {0} either doesn\'t exist or cannot be open'.format(self.file_list[0])
            t = rootfile.Get('limit')
            
            self.log.debug('Required branches are : {0}'.format(required_branches))
            for axis in range(dims):
                assert t.GetListOfBranches().FindObject(required_branches[axis]), "The branch \"{0}\" doesn't exist.".format(required_branches[axis])
                
            t.SetBranchStatus("*", False)    
            for axis in range(dims):
                t.SetBranchStatus(required_branches[axis], True)
            
            
            for en in range(1,t.GetEntriesFast()):
                t.GetEntry(en)
                if dims==2:
                    X = eval(contour_axis_list[0])
                    Y = eval(contour_axis_list[1])
                    if en%100 == 0:
                        #self.log.debug('Inputs X={0} Y={1}'.format(eval('t.{0}'.format(required_branches[0])),eval('t.{0}'.format(required_branches[1]))))
                        self.log.debug('Entry={2} X={0} Y={1}'.format(X,Y, en))
                    self.contour_graph[contour_axis].SetPoint(en-1,X,Y)
                    
                elif dims==3:
                    X = eval(contour_axis_list[0])
                    Y = eval(contour_axis_list[1])
                    Z = eval(contour_axis_list[2])
                    if en%100 == 0:
                        #t.Show()
                        self.log.debug('Entry={3} X={0} Y={1} Z={2}'.format(X,Y,Z, en))
                    self.contour_graph[contour_axis] = TGraph2D(en-1,X,Y,Z)
                
        return self.contour_graph[contour_axis]
        
        
    def ll_values(self, POI, cl = 0.68):
        """returns a list of lower limits for a given level for a given POI
        """
        if not self._has_parsed_combine_result_already:
            self._parse_combine_result()
        assert belongsTo(float(cl),0,1), "Confidence level has to be given in interval [0,1]"
        cl_name = "{1}_CL@{0:.2f}".format(cl, POI)
        try:
            self.ll_values_dict[cl_name]
        except KeyError:
            self.ll_values_dict[cl_name] = self._get_crossings_for_limits(self.raising_segments[POI], float(cl))
            self.log.debug('Creating limit for C.L.@{0}'.format(cl))
        else:    
            self.log.debug('Returning existing limit for C.L.@{0}'.format(cl))
        return self.ll_values_dict[cl_name]
    
    def ul_values(self, POI, cl = 0.68):
        """returns a list of lower limits for a given level for a given POI
        """
        if not self._has_parsed_combine_result_already:
            self._parse_combine_result()
        assert belongsTo(float(cl),0,1), "Confidence level has to be given in interval [0,1]"
        cl_name = "{1}_CL@{0:.2f}".format(cl, POI)
        try:
            self.ul_values_dict[cl_name]
        except KeyError:
            self.ul_values_dict[cl_name] = self._get_crossings_for_limits(self.falling_segments[POI], float(cl))
            self.log.debug('Creating limit for C.L.@{0}'.format(cl))
        else:    
            self.log.debug('Returning existing limit for C.L.@{0}'.format(cl))
            
        return self.ul_values_dict[cl_name]
        
    def best_fit(self, POI):
        """Get the best fit value fora particular POI
        """
        if not self._has_parsed_combine_result_already:
            self._parse_combine_result()
        try:
            self.global_best_fit[POI]
        except KeyError:
            raise KeyError, 'The POI name \"{0}\" is invalid.'.format(POI)
        else:    
            return self.global_best_fit[POI]
        
        
    def is_set_ll(self,POI, cl=0.68):
        assert belongsTo(float(cl),0,1), "Confidence level has to be given in interval [0,1]"
        cl_name = "{1}_CL@{0:.2f}".format(cl, POI)
        try:
            self.ll_values_dict[cl_name]
        except KeyError:
            raise KeyError, 'The POI name \"{0}\" is invalid.'.format(POI)
        else:    
            return (len(self.ll_values_dict[cl_name])>0)

    def is_set_ul(self,POI, cl=0.68):
        assert belongsTo(float(cl),0,1), "Confidence level has to be given in interval [0,1]"
        cl_name = "{1}_CL@{0:.2f}".format(cl, POI)
        try:
            self.ul_values_dict[cl_name]
        except KeyError:
            raise KeyError, 'The POI name \"{0}\" is invalid.'.format(POI)
        else:    
            return (len(self.ul_values_dict[cl_name])>0)
        
    def is_set_best_fit(self,POI):
        try:
            self.global_best_fit[POI]
        except KeyError:
            raise KeyError, 'The POI name \"{0}\" is invalid.'.format(POI)
        else:    
            return (self.global_best_fit[POI]!=None)
            
    def get_results_dict(self, POI, option='standard'):
        """Returns a dict with best fit values and limits at 68(95)%
        """
        self.log.info('Compiling the fit results dictionary...')
        if option.lower()=='standard':
            import collections
            self.limits_dict = collections.OrderedDict()
            self.limits_dict['BF']  = self.best_fit(POI)
            self.limits_dict['LL68']= self.ll_values(POI, 0.68)
            self.limits_dict['LL95']= self.ll_values(POI, 0.95)
            self.limits_dict['UL68']= self.ul_values(POI, 0.68)
            self.limits_dict['UL95']= self.ul_values(POI, 0.95)
            
            self.log.debug('Limits are: {0}'.format(self.limits_dict))    
            return self.limits_dict
        else:
            print 'The option {0} is still not implemented. Do you want to volonteer? :)'.format(option)
            
            
        
    def contours(self, cl = 0.68):
        """Return list of TGraph contours with a given confidence level.
        """
        pass
        return self.contours
        
    def set_combine_method(self,combine_method):
        """Set method in order to know how the limit trees look like. 
        """
        self.combine_method = combine_method
        
    def _get_TGraph_from_segment(self, segment):
        """Create TGraph from list of tuples(qe, poi, dNLL)
        """
        qe_vals  = array('d', [t[0] for t in segment])
        poi_vals = array('d', [t[1] for t in segment])
        return TGraph(len(segment), qe_vals, poi_vals)
    
    def _parse_combine_result(self, combine_method="MultiDimFit"):
        """Parsing the combine result and filling the contour information.
           Should be run on first demand of any information.
        """
        self._has_parsed_combine_result_already = True
        self.log.debug("Parsing the combine output files for method = {0}".format(combine_method))
        assert len(self.file_list)>0, "There is no files to read."
        #assert len(self.file_list)==1, "Currently, we allow only one file from combine output. The implementation is still lacking this feature..."
        
        #containers for TGraph raising and falling segments
        self.falling_segments =  {poi:[] for poi in self.POI}
        self.raising_segments = {poi:[] for poi in self.POI}
        
        for poi in self.POI:
            for root_file_name in self.file_list:
                root_file_name =  self.file_list[0]
                #get the TTree and enable relevant branches (POI, quantileExpected, deltaNLL)
                self.log.debug("Parsing the combine output file = {0}".format(root_file_name))
                rootfile = ROOT.TFile.Open(root_file_name,'READ')
                if not rootfile:
                    raise IOError, 'The file {0} either doesn\'t exist or cannot be open'.format(root_file_name)
                t = rootfile.Get('limit')
                assert t.GetListOfBranches().FindObject(poi), "The branch \"{0}\" doesn't exist.".format()
                
                #don't read uninteresting branches
                t.SetBranchStatus("*", False)
                t.SetBranchStatus("quantileExpected", True)
                t.SetBranchStatus(poi, True)
                t.SetBranchStatus("deltaNLL", True)
                
                
                #get the best fit
                if self.global_best_fit[poi] == None:
                    t.GetEntry(0)
                    if AreSame(t.quantileExpected,1) and AreSame(t.deltaNLL, 0):  #This is true if combine has found a good minimum.
                        self.global_best_fit[poi] = eval('t.{0}'.format(poi))
                        self.global_best_fit_dict.update({poi : {'best_fit' : eval('t.{0}'.format(poi)), 'quantileExpected' : t.quantileExpected, '2*deltaNLL' : 2*t.deltaNLL}} )
                        self.log.debug("Global best fit in file {0} = {1}".format(root_file_name, self.global_best_fit_dict))
                
                
                is_raising = False
                #we want to find the first trend of the function
                ien = 1
                tmp_list_qe_poi_dNLL = []
                while True:
                    t.GetEntry(ien)
                    if ien==1:
                        qe_prev = t.quantileExpected
                    if t.quantileExpected > qe_prev:
                        is_raising=True
                        break
                    elif t.quantileExpected < qe_prev:
                        is_raising=False
                        break
                    qe_prev = t.quantileExpected    
                    ien+=1
                self.log.debug('Detected trend of first interval: Raising = {0}, Falling = {1}'.format(is_raising, (not is_raising)))
                        
                self.log.debug('Searching for intervals at 68(95)% C.L.')
                tmp_list_qe_poi_dNLL = []
                is_trend_changed = False
                n_entries = t.GetEntriesFast()
                for en in range(1,n_entries+1):  #add +1 to range so that to save all the entries into the segments
                    if en < n_entries:  #if not passed the last entry
                        t.GetEntry(en)
                        # set x=quantileExpected and y=POI and create segments
                        qe, poi_value, dNLL  = t.quantileExpected, eval('t.{0}'.format(poi)), t.deltaNLL
                        
                        #check if trend is change, and then change the bool is_raising which will show which vector should be filled
                        if (en > 2):
                            if ((qe < qe_prev) and is_raising) or ((qe > qe_prev) and not is_raising):
                                is_trend_changed = True
                                self.log.debug('Trend of the segment has changed at entry {0}.'.format(en))
                                self.log.debug('**********************************************')
                        
                    #fill tmp_list_qe_poi_dNLL until we see that the trend is changed
                    #then put that list into self.raising_segments or self.falling_segments and clear the tmp_list_qe_poi_dNLL
                    if is_trend_changed or en==n_entries:
                        if is_raising:
                            
                            self.raising_segments[poi].append(self._get_TGraph_from_segment(tmp_list_qe_poi_dNLL))
                            #self.log.debug('Appending segment to self.raising_segments[{1}]: {0}'.format(tmp_list_qe_poi_dNLL, poi))
                        else:
                            self.falling_segments[poi].append(self._get_TGraph_from_segment(tmp_list_qe_poi_dNLL))
                            #self.log.debug('Appending segment to self.falling_segments[{1}]: {0}'.format(tmp_list_qe_poi_dNLL, poi))
                        
                        self.log.debug('Raising_segments state: {0}'.format(self.raising_segments))    
                        self.log.debug('Falling_segments state: {0}'.format(self.falling_segments))    
                        #delete all elements from the tmp list    
                        del tmp_list_qe_poi_dNLL[:]
                        #change the state of raising/falling interval    
                        is_raising = (not is_raising)
                        is_trend_changed = False
                        
                    if en < n_entries:    
                        #fill tmp_list_qe_poi_dNLL
                        tmp_list_qe_poi_dNLL.append((qe, poi_value, dNLL))
                        #self.log.debug('Entry = {2} Raising = {0}, Falling = {1}'.format(is_raising, (not is_raising), en))
                        #self.log.debug('Entry = {4} Filling tmp_list_qe_poi_dNLL (size={0}) with : qe = {1}, poi_value = {2}, dNLL = {3}'.format(len(tmp_list_qe_poi_dNLL),qe, poi_value, dNLL, en ))
                        qe_prev = qe
                        #is_trend_changed = False




#def getArcTanParam(self, k3k1,  gamma=1):
    #return (1/TMath.Pi())*TMath.ATan(TMath.Sqrt(gamma)*k3k1)

#FitResult_t getRescaledFitResult(FitResult_t fitres,  gamma)
##     rescale all values in fit result
    #FitResult_t nfr
    #nfr.bf_global = getArcTanParam(fitres.bf_global,gamma)
    #nfr.bf = getArcTanParam(fitres.bf,gamma)
    #nfr.wf = getArcTanParam(fitres.wf,gamma)
    
    # ulimit68 = fitres.ul68
    # ulimit95 = fitres.ul95
    # llimit68 = fitres.ll68
    # llimit95 = fitres.ll95
    #nfr.ll68 = (AreSame(llimit68,default_fr.ll68)) ? default_fr.ll68 :  getArcTanParam(fitres.ll68,gamma)
    #nfr.ll95 = (AreSame(llimit95,default_fr.ll95)) ? default_fr.ll95 :  getArcTanParam(fitres.ll95,gamma)
    #nfr.ul68 = (AreSame(ulimit68,default_fr.ul68)) ? default_fr.ul68 :  getArcTanParam(fitres.ul68,gamma)
    #nfr.ul95 = (AreSame(ulimit95,default_fr.ul95)) ? default_fr.ul95 :  getArcTanParam(fitres.ul95,gamma)
#return nfr    


#FitResult_t getfa3RescaledFitResult(FitResult_t fitres)
##     rescale all values in fit result
    #FitResult_t nfr= -99999,-99999,-99999,99999,99999,-99999,-99999 
    #nfr.bf_global = fa3(fitres.bf_global,0)
    #nfr.bf = fa3(fitres.bf,0)
    #nfr.wf = fa3(fitres.wf,0)
    
    # ulimit68 = fitres.ul68
    # ulimit95 = fitres.ul95
    # llimit68 = fitres.ll68
    # llimit95 = fitres.ll95
    #nfr.ll68 = (AreSame(llimit68,default_fr.ll68)) ? default_fr.ll68 :  fa3(fitres.ll68,0)
    #nfr.ll95 = (AreSame(llimit95,default_fr.ll95)) ? default_fr.ll95 :  fa3(fitres.ll95,0)
    #nfr.ul68 = (AreSame(ulimit68,default_fr.ul68)) ? default_fr.ul68 :  fa3(fitres.ul68,0)
    #nfr.ul95 = (AreSame(ulimit95,default_fr.ul95)) ? default_fr.ul95 :  fa3(fitres.ul95,0)
#return nfr    



#TGraph * getfa3RescaledGraph(TGraph *gr)
  #Int_t N = gr.GetN()
  #print "N = " N  
  # *x
  #vector<> new_x
  # *y
  
  #x = gr.GetX()
  #y = gr.GetY()
##   print "N = " N  
  #for (Int_t i=0 i<N i++)
    #new_x.push_back(fa3(x[i],0))
##     print "x, new_x = "  x[i]  " " new_x[i] 
  
  
  #TGraph * gr_ret = new TGraph(N,&(new_x.front()),y)
##   TGraph * gr_ret = new TGraph(N,x,y)
  #return gr_ret



#def getRescaledGraph(self, gr, gamma):
    #N = gr.GetN()
    #print "N = ",N  
    #x = gr.GetX()
    #y = gr.GetY()
    #new_x = []
    #for i in range(N):
        #new_x[i] = self.getArcTanParam(x[i],gamma)
        ##print "x= {0} new_x = {1}".format(x[i],new_x[i]) 
    #gr_ret = TGraph(N,array('d', new_x),y)
    #return gr_ret



#def getGamma33(self, channel=0):
    #g=0.038
    #if channel==0:  
        #g = 0.040
    ##         g = (0.46154957824608*0.040+0.3543010974318054*0.034+0.18414932432211456*0.034 ) #weighted mean in case of calculation for all channels
    #elif channel==3: 
        #g = 0.040
    #else:
        #g=0.034
    ##print "@@@@@ gamma33 (ch={0})={1}".format(channel,g)
    #return g

#def fa3(self, k3k1, channel):
    #sgn = 1
    #if k3k1<0: 
        #sgn = -1
        ##print "@@@@ fa3Calculator: " 
    #g = self.getGamma33(channel)
    #return sgn * g * k3k1 * k3k1/(1 + g * k3k1 * k3k1)


#if __name__ == "__main__":
     
    #print "Starting plots..."
    
    #pp = pprint.PrettyPrinter(indent=4)
    
    #DEBUG = True
    ##cfg_reader = UniversalConfigParser(cfg_type="YAML",file_list = "plotLimitVsLumi.yaml")
    ##cfg_reader.setLogLevel(10) 
    
    ##plots_data = cfg_reader.get_dict()
    
    ##if DEBUG: 
          ##print "plots_data = ", 
          ##pp.pprint(plots_data)
    
    #plotter = LikelihoodFitContourPlot()
    #plotter.setCopyToWebDir(False,"/afs/cern.ch/user/r/roko/www/html/Geoloc/")
    #plotter.make_plot()
    
