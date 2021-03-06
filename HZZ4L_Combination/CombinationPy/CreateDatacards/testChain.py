#!/usr/bin/python
#-----------------------------------------------
# Roko Plestina, 2013-2014
#-----------------------------------------------
import sys, os, pwd, commands
import optparse, shlex, re
import math
import errno
import shutil
import pprint
from array import array
from lib.util.Logger import *

def parseOptions():

    usage = ('usage: %prog [options] \n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)
    parser.add_option('-c', '--cfg', dest='config_filename', type='string', default="testChain.yaml",    help='Name of the file with full configuration')
    parser.add_option('-D','--datacard', dest='datacard', type='int', default=1, help='Toggle the datacard creation (0 or 1). Deafault [1]')
    parser.add_option('-r', '--run', dest='run_data_name', type='string', default="",    help='Name of the section to be run')
    parser.add_option('-p', '--print', action='store_true', dest='print_config_and_exit', default=False ,help='Just print the configuration and exit')
    parser.add_option('-w', '--www', action='store_true', dest='do_copy_to_webdir', default=False, help='Copy plots and relevant stuff to a webdir')
    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

      
class ChainProcessor(object):
    """
    Processes the configuration set for running the creation 
    of Datacards, combine fit and then plot unit
    """
    
    

    def __init__(self, run_data_name, run_dict):
        """Set the values from configuration"""
        
        self.log = Logger().getLogger(self.__class__.__name__, 10)
        self.run_dict = run_dict
        self.create_cards_dir="/afs/cern.ch/work/r/roko/Stat/CMSSW_611_JCP/src/HZZ4L_Combination/CombinationPy/CreateDatacards/"
        self.version= "v1"
        self.sqrts_dc="8TeV"
        self.sqrts_plot="14"
        self.fs="2e2mu"
        self.lumi="19.79"
        self.lumi_zfill="19.79"
        self.discriminant="D_{0-}(90^{o})"
        self.MY_CUR_WWW_SUBDIR=""
        self.append_name_base = run_data_name
        self.MY_CUR_WWW_SUBDIR = self.append_name_base
        self.POIs = list(self.run_dict['POI'])
        self.poi = self.run_dict['POI_setup']
        self.termNames =  self.run_dict['termNames']
        self.templatedir = self.run_dict['templatedir']
        self.discriminant = self.run_dict['discriminant_name']
        self.additional_options = self.run_dict['additional_options']
        self.fs  = self.run_dict['final_state']
        self.www_dir = "/afs/cern.ch/user/r/roko/www/html/Geoloc/{0}".format(self.MY_CUR_WWW_SUBDIR)
        self.do_copy_to_webdir = True
        print "Setting web directory to: {0}".format(self.www_dir)
        self.do_cmd = ""
        #self.poi= {"k2k1_ratio":{"switch":0, "range":"-100,100", "name":"k2k1", "nice_name":"k_{2}/k_{1}", "value":"0"}, 
              #"k3k1_ratio":{"switch":0, "range":"-100,100", "name":"k3k1", "nice_name":"k_{3}/k_{1}", "value":"0"}}
        self.setup_POI_related_stuff()
            
        
        
    def setup_POI_related_stuff(self):
        """Sets up the information which is related 
           to POIs and is needed to run the commands. 
        """
        def _set_POI_range(self,this_poi):
            """We use this function to set ranges to parameters because we set them lumi-dependant.
            """
            if this_poi=="k3k1_ratio":
                theRange=""
                if float(self.lumi) >  50: 
                    theRange = "-10,10"
                elif float(self.lumi) >  20: 
                    theRange = "-20,20"
                else : 
                    theRange = "-30,30"
                self.poi["k3k1_ratio"]["range"] = theRange
            elif this_poi=="k2k1_ratio":
                pass        
        
        #make list out of POI
        if isinstance(self.run_dict['POI'],list):
            pass
        else:
            tmp_list=[]
            tmp_list.append(self.run_dict['POI'])
            self.run_dict['POI'] = tmp_list
        
        self.POIs = self.run_dict['POI']
        self.log.debug('POI: {0}'.format(self.POIs))
        assert len(self.POIs)>0 and len(self.POIs)<3 and ("k2k1_ratio" or "k3k1_ratio" in self.POIs), \
            "You should provide at least one POI (k2k1_ratio, k3k1_ratio), while you provided {0}".format(self.POIs)
        
        if "k3k1_ratio" in self.POIs and "k2k1_ratio" in self.POIs:
            self.poi_physics_model = "HiggsAnalysis.CombinedLimit.GeoLocationModel:K1andK2andK3Model"
            self.poi_ranges_string_t2w = "--PO range_k3k1_ratio={0} --PO range_k2k1_ratio={1}".format(self.poi["k3k1_ratio"]['range'],self.poi["k2k1_ratio"]['range'] )
            _set_POI_range(self,"k3k1_ratio")
            _set_POI_range(self,"k2k1_ratio")
            self.poi_ranges_string_fit = "k3k1_ratio={0}:k2k1_ratio={1}".format(self.poi["k3k1_ratio"]['range'],self.poi["k2k1_ratio"]['range'])
            self.poi_name_value_plot = "{0}={2} {1}={3}".format(self.poi["k3k1_ratio"]['nice_name'],self.poi["k2k1_ratio"]['nice_name'],self.poi["k3k1_ratio"]['value'], self.poi["k2k1_ratio"]['value'])
            self.poi_name_value = "k3k1_ratio={0},k2k1_ratio={1}".format(self.poi["k3k1_ratio"]['value'],self.poi["k2k1_ratio"]['value'])
            self.poi_name_value_filename = "k3k1_ratio_{0}_k2k1_ratio_{1}".format(self.poi["k3k1_ratio"]['value'],self.poi["k2k1_ratio"]['value'])
            self.pois = "k2k1_ratio,k3k1_ratio"
            self.poi_n_points = str(self.poi["k2k1_ratio"]['n_scan_points']*self.poi["k3k1_ratio"]['n_scan_points'])
            self.log.debug('Set up the {0}'.format(self.poi_physics_model))

            
        elif "k2k1_ratio" in self.POIs:
            self.poi_physics_model = "HiggsAnalysis.CombinedLimit.GeoLocationModel:K1andK2Model"
            self.poi_ranges_string_t2w = "--PO range_k2k1_ratio={0}".format(self.poi["k2k1_ratio"]['range'] )
            _set_POI_range(self,"k2k1_ratio")
            self.poi_ranges_string_fit = "k2k1_ratio={0}".format(self.poi["k2k1_ratio"]['range'])
            self.poi_name_value_plot = "{0}={1}".format(self.poi["k2k1_ratio"]['nice_name'],self.poi["k2k1_ratio"]['value'])
            self.poi_name_value = "k2k1_ratio={0}".format(self.poi["k2k1_ratio"]['value'])
            self.poi_name_value_filename = "k2k1_ratio_{0}".format(self.poi["k2k1_ratio"]['value'])
            self.pois = "k2k1_ratio"
            self.poi_n_points = str(self.poi["k2k1_ratio"]['n_scan_points'])
            self.log.debug('Set up the {0}'.format(self.poi_physics_model))

            
        elif "k3k1_ratio" in self.POIs:
            self.poi_physics_model = "HiggsAnalysis.CombinedLimit.GeoLocationModel:K1andK3Model"
            self.poi_ranges_string_t2w = "--PO range_k3k1_ratio={0}".format(self.poi["k3k1_ratio"]['range'] )
            _set_POI_range(self,"k3k1_ratio")
            self.poi_ranges_string_fit = "k3k1_ratio={0}".format(self.poi["k3k1_ratio"]['range'])
            self.poi_name_value_plot = "{0}={1}".format(self.poi["k3k1_ratio"]['nice_name'],self.poi["k3k1_ratio"]['value'])
            self.poi_name_value = "k3k1_ratio={0}".format(self.poi["k3k1_ratio"]['value'])
            self.poi_name_value_filename = "k3k1_ratio_{0}".format(self.poi["k3k1_ratio"]['value'])
            self.pois = "k3k1_ratio"
            self.poi_n_points = str(self.poi["k3k1_ratio"]['n_scan_points'])
            self.log.debug('Set up the {0}'.format(self.poi_physics_model))
        
        self.log.debug("poi_physics_model ={0}".format(self.poi_physics_model))
        self.log.debug("poi_ranges_string_t2w ={0}".format(self.poi_ranges_string_t2w))
        self.log.debug("poi_ranges_string_fit ={0}".format(self.poi_ranges_string_fit))
        self.log.debug("poi_name_value_plot ={0}".format(self.poi_name_value_plot))
        self.log.debug("poi_name_value ={0}".format(self.poi_name_value))
        self.log.debug("poi_name_value_filename ={0}".format(self.poi_name_value_filename))
        self.log.debug("poi_n_points ={0}".format(self.poi_n_points))
        self.log.debug("pois ={0}".format(self.pois))

        
    def process(self, do_cmd=""):
        self.do_cmd = do_cmd
        self.setup_POI_related_stuff()

        import lib.util.MiscTools as misc
        #print self.__dict__
        self.fs_expanded = self.fs
        if self.fs.lower() == "4l":
            self.fs_expanded = "4e,4#mu,2e2#mu"
        os.environ['PLOT_TAG'] = "Asimov data %(poi_name_value_plot)s | Discrim. %(discriminant)s | L=%(lumi)s fb^{-1} @ %(sqrts_plot)s TeV | Fin. state = %(fs_expanded)s" %self.__dict__

        cmd = {}
        cmd['createCards']    = "rm -r cards_%(append_name_base)s; python makeDCsandWSs.py -b -i SM_inputs_8TeV_CJLST -a %(append_name_base)s -t %(templatedir)s --terms %(termNames)s %(additional_options)s" %self.__dict__
        cmd['combCards']      = "rm -rf hzz4l_4lS_%(sqrts_dc)s_ALT.txt; combineCards.py hzz4l_4muS_%(sqrts_dc)s_ALT.txt hzz4l_4eS_%(sqrts_dc)s_ALT.txt hzz4l_2e2muS_%(sqrts_dc)s_ALT.txt> hzz4l_4lS_%(sqrts_dc)s_ALT.txt" %self.__dict__
        cmd['t2w'] 		= "text2workspace.py hzz4l_%(fs)sS_%(sqrts_dc)s_ALT.txt -m 126 -P %(poi_physics_model)s %(poi_ranges_string_t2w)s  --PO muFloating -o combine.ws.%(fs)s.%(version)s.root" %self.__dict__
        cmd['gen'] 		= "combine -M GenerateOnly combine.ws.%(fs)s.%(version)s.root -m 126 -t -1 --expectSignal=1 --saveToys --setPhysicsModelParameters %(poi_name_value)s,cmshzz4l_lumi=%(lumi)s -S 0" %self.__dict__
        cmd['addasimov'] 	= "root -b -l -q %(create_cards_dir)s/addToyDataset.C\(\\\"combine.ws.%(fs)s.%(version)s.root\\\",\\\"higgsCombineTest.GenerateOnly.mH126.123456.root\\\",\\\"toy_asimov\\\",\\\"workspaceWithAsimov_%(poi_name_value_filename)s_lumi_%(lumi)s.root\\\"\)" %self.__dict__
        cmd['fit'] 		= "combine -M MultiDimFit workspaceWithAsimov_%(poi_name_value_filename)s_lumi_%(lumi)s.root --algo=grid --points %(poi_n_points)s  -m 126 -n .asimov.%(fs)s.%(poi_name_value_filename)s.lumi_%(lumi_zfill)s -D toys/toy_asimov -S 0 --setPhysicsModelParameters cmshzz4l_lumi=%(lumi)s --setPhysicsModelParameterRanges %(poi_ranges_string_fit)s" %self.__dict__
        cmd['plot'] 	       = "root -l -b -q %(create_cards_dir)s/plotLimit.C\(\\\"higgsCombine.asimov.%(fs)s.%(poi_name_value_filename)s.lumi_%(lumi_zfill)s.MultiDimFit.mH126.root\\\",\\\"%(pois)s\\\",\\\"environ PLOT_TAG\\\" \)" %self.__dict__
        
        if "createCards" in do_cmd:
            print "--------------------------------------------------------"
            print "Cards for Lumi=%(lumi)s and tag=%(append_name_base)s terms=%(termNames)s template=%(templatedir)s discriminant=%(discriminant)s additional_options=%(additional_options)s" %self.__dict__
            misc.processCmd(cmd['createCards'])
            print "--------------------------------------------------------"
        if "combCards" in do_cmd: misc.processCmd(cmd['combCards']) 
        if "t2w" in do_cmd: misc.processCmd(cmd['t2w']) 
        if "gen" in do_cmd: misc.processCmd(cmd['gen']) 
        if "addasimov" in do_cmd: misc.processCmd(cmd['addasimov']) 
        if "fit" in do_cmd: misc.processCmd(cmd['fit']) 
        if "plot" in do_cmd: 
            misc.processCmd(cmd['plot']) 
            import lib.util.MiscTools as misc
            misc.make_sure_path_exists(self.www_dir)
            shutil.copy("/afs/cern.ch/user/r/roko/www/html/index.php",self.www_dir)
            for self.x in ['AxisK3K1', 'AxisArctanGamK3K1', 'AxisArctanK3K1', 'AxisFa3', 'AxisK2K1']:
                #print "cp higgsCombine.asimov.{fs}.{poi_name_value}.MultiDimFit.mH126.root*png {www_dir}/exp.{fs}.{poi_name_value_filename}.mH126.{x}.png".format(**self.__dict__)
                #shutil.copy("higgsCombine.asimov.{fs}.{poi_name_value_filename}.lumi_{lumi_zfill}.MultiDimFit.mH126.root.{x}.png".format(**self.__dict__), \
                            #"{www_dir}/exp.{fs}.{poi_name_value_filename}.mH126.{x}.{lumi_zfill}ifb.png".format(**self.__dict__))
                file_2_copy = "higgsCombine.asimov.{fs}.{poi_name_value_filename}.lumi_{lumi_zfill}.MultiDimFit.mH126.root.{x}.png".format(**self.__dict__)            
                if os.path.exists(file_2_copy):
                    shutil.copy(file_2_copy,"{www_dir}/exp.{fs}.{poi_name_value_filename}.mH126.{x}.{lumi_zfill}ifb.png".format(**self.__dict__))
            print "{www_dir}/?match=exp.{fs}.{poi_name_value_filename}.mH126.*.png".format(**self.__dict__).replace("/afs/cern.ch/user/r/roko/www/html/","https://roko.web.cern.ch/roko/")
            
    def get_table_row(self, col_names=False) :
        if col_names: 
            return "LUMI BF UL68 UL95 WF\n".lower()
        import lib.util.MiscTools as misc
        limit_file = "higgsCombine.asimov.%(fs)s.%(poi_name_value_filename)s.lumi_%(lumi_zfill)s.MultiDimFit.mH126.root.limits" %self.__dict__
        self.log.debug("Searching for limits info in file: {0}".format(limit_file))
        BF   = misc.grep("BF", limit_file)[0].split(":")[1].strip()
        UL68 = misc.grep("UL68", limit_file)[0].split(":")[1].strip()
        UL95 = misc.grep("UL95", limit_file)[0].split(":")[1].strip()
        WF   = misc.grep("WF", limit_file)[0].split(":")[1].strip()
        table_raw = "{0} {1} {2} {3} {4}\n".format(self.lumi, BF, UL68,UL95, WF)
        print "Table raw: ", table_raw
        return table_raw
        
        
    def set_lumi(self, lumi, n_digits=0) : 
        self.lumi = str(lumi)
        self.lumi_zfill = self.lumi.zfill(n_digits)
        
    def set_sqrts_dc(self, sqrts) : self.sqrts_dc = str(sqrts)
    def set_sqrts_plot(self, sqrts) : self.sqrts_plot = str(sqrts)
    def set_poi_value(self, poi_name, value) : self.poi[poi_name]['value']= str(value)
    def get_webdir(self): return self.www_dir
    def set_do_copy_to_webdir(self, do_copy=True): self.do_copy_to_webdir = do_copy
    def get_sqrts_dc(self): return self.sqrts_dc
    def get_fs(self): return self.fs
    def get_cmd(self): return self.do_cmd
    def get_poi_name_value(self) : return self.poi_name_value
    def get_poi_name_value_filename(self) : return self.poi_name_value_filename
    def get_run_data(self): return self.run_dict
    
    

if __name__ == "__main__":
    """
    This script can create datacards, make toy datasets, fit and plt results.
    Finally it is coping the relevalnt output to the web directory...
    
    """
    # parse the arguments and options
    global opt, args
    parseOptions()
    
    log = Logger().getLogger(__name__, 10)

    #read configuration
    import lib.util.UniversalConfigParser as ucp
    cfg_reader = ucp.UniversalConfigParser(cfg_type="YAML",file_list = opt.config_filename)
    cfg_reader.setLogLevel(10)
    full_config = cfg_reader.get_dict()
    run_data = full_config["COMMON"]
    run_data.update(full_config[opt.run_data_name])
    pp = pprint.PrettyPrinter(indent=4)
    #pp.pprint(full_config)
        
    print "@@@@ THE RUN DATA: YAML -> DICT"
    pp.pprint(run_data)
    if opt.print_config_and_exit:
        quit()

   ############################
    chain_process = ChainProcessor(opt.run_data_name, run_data)
    #common_data = full_config["COMMON"]
    #print "@@@@ THE COMMON DATA:"
    #pp.pprint(common_data)
    chain_process.set_sqrts_dc(str(run_data['sqrts_dc']))
    chain_process.set_sqrts_plot(str(run_data['sqrts_plot']))
    chain_process.set_do_copy_to_webdir(opt.do_copy_to_webdir)
        
    for poi in list(run_data['POI']):
        log.debug('Setting POI value for {0}'.format(poi))
        chain_process.set_poi_value(str(poi), str(run_data[poi]))
        
    
    lumi_list = run_data['lumi_list']
  
    this_dir=os.getcwd()
    cards_dir = "cards_{0}/HCG/126".format(opt.run_data_name)
    #out_file_name_base = "{0}.{1}.asimov.k3k1.{2}".format(opt.run_data_name, run_data['final_state'],run_data['k3k1_ratio'])
    out_file_name_base = "{0}.{1}.asimov.{2}".format(opt.run_data_name, run_data['final_state'],chain_process.get_poi_name_value_filename())
    if 'info' not in run_data.keys(): 
        run_data['info']='Not available...'
    log.info(run_data['info'])
    info_name = out_file_name_base+".info"
    tab_name = out_file_name_base+".tab"
    log.debug("tab_name={0}".format(tab_name))
    
    import yaml
    with open('data_config.yaml', 'w') as fdata:
        fdata.write(yaml.dump(run_data, default_flow_style=False))
        log.info('Writing used configuration to data_config.yaml')
    
    
    if opt.datacard:
        chain_process.process("createCards")
        try:
            os.remove(info_name)
        except OSError , exception:
            if exception.errno != errno.EEXIST:
                pass
            
        try:
            os.remove(tab_name)
        except OSError , exception:
            if exception.errno != errno.EEXIST:
                pass
    os.chdir(cards_dir)        
    log.debug("Current dir = {0}".format(os.getcwd()))
    chain_process.process("combCards,t2w")
    os.chdir(this_dir)
    log.debug("Current dir = {0}".format(os.getcwd()))        
    print "--------------------------------------------------------"
    
    finfo = open(info_name, "w")
    finfo.write(run_data['info'])
    finfo.close()
    
    f = open(tab_name, "w")
    f.write(chain_process.get_table_row(col_names=1))
    
    for lumi in lumi_list:
        chain_process.set_lumi(lumi, 4)
        
        os.chdir(cards_dir)
        log.debug("Current dir = {0}".format(os.getcwd()))
        chain_process.process("gen,addasimov,fit,plot")
        #chain_process.process("plot")
        table_raw = chain_process.get_table_row()
        os.chdir(this_dir)
        log.debug("Current dir = {0}".format(os.getcwd()))        
        f.write(table_raw)
        #print "%({lumi} %({bf} %({ul68} %({ul95} %({wf}" >> %(tab_name
        
    f.close()
    #copy info, table and DC to webdir
    shutil.copy("{0}/{1}".format(this_dir,tab_name),"{0}/{1}".format(cards_dir, tab_name))
    shutil.copy("{0}/{1}".format(this_dir,info_name),"{0}/{1}".format(cards_dir, info_name))
    shutil.copy("{0}/data_config.yaml".format(this_dir,info_name),"{0}/data_config.yaml".format(cards_dir))
    
    
    if opt.do_copy_to_webdir:
        import lib.util.MiscTools as misc
        webdir = chain_process.get_webdir()
        misc.make_sure_path_exists(webdir)
        shutil.copy("/afs/cern.ch/user/r/roko/www/html/index.php",webdir)
        
        file_list_for_moving = {
            "{0}/{1}".format(this_dir,tab_name) : "{0}/{1}".format(webdir, "results.tab"),
            "{0}/{1}".format(this_dir,info_name) : "{0}/{1}".format(webdir, "basic.info"),
            "{0}/data_config.yaml".format(this_dir) : "{0}/data_config.yaml".format(webdir)
            }
            
        for file in file_list_for_moving.values():
            if os.path.exists(file) and os.path.isfile(file): 
                os.remove(file)
        for file in file_list_for_moving.keys():
            shutil.move(file,file_list_for_moving[file])
        
        if opt.datacard:
            shutil.copy("{0}/hzz4l_4lS_{1}_ALT.txt".format(cards_dir,chain_process.get_sqrts_dc()),"{0}/{1}".format(webdir, "combine.card"))
            shutil.copy("{0}/combine.ws.{1}.v1.root".format(cards_dir,chain_process.get_fs()),"{0}/{1}".format(webdir, "combine.ws.root"))
            
        ############################

        
        
        
        
