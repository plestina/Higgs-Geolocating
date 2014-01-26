#! /usr/bin/env python
#import pprint
from Logger import *



class UniversalConfigParser(object):
    """
    Wraper for XML, YAML, JSON, SimpleConfigFile which can return dictionary.
    
    """
    #from lib.util.Logger import *
    supported_cfg_types = ['yaml','xml','json','ini']
    def __init__(self, cfg_type=None, file_list = None):
	self.log = Logger().getLogger(self.__class__.__name__, level=10) #set to DEBUG level
	
	self.cfg_type = self.set_cfg_type(cfg_type)
	self.file_list = self.set_files(file_list)
	self.cfg_dict = {}
	pass
      
      
    def _get_cfg_type(self,file_name):
      
	if self.cfg_type == None : 
	    self.log.debug('The type of CFG file is being be figured out from file extension.')
	    import os.path
	    import re
	    cfg_type = re.sub(r".",'',os.path.splitext()[1]).lower()
	    
	assert cfg_type in self.supported_cfg_types,"Cannont figure out the configuration file type. Supported CFG types are only {0}".format(self.supported_cfg_types) 
	return cfg_type
	
    def set_cfg_type(self,cfg_type):
	self.log.debug('Checking the validity of CFG type.')
	
	if cfg_type == None: return None    
	try: cfg_type  = cfg_type.lower()
	except AttributeError: 
	    raise AttributeError, "You have to provide a string as CFG type."
	assert cfg_type in self.supported_cfg_types,"CFG type not supported. Supported CFG types are only {0}".format(self.supported_cfg_types) 
	return cfg_type
	 
    def setLogLevel(self,level): self.log.setLevel(level)
    
    def set_files(self, file_list):
	"""
	Check if list is list or a CSV string. In case of string, create a list. 
	"""
	if file_list==None: return []
	import types
	isString = isinstance(file_list, types.StringTypes) 
	isList = isinstance(file_list, list) 
	assert isString or isList, "You should provide a list of files as list or as CVS string!"
	if isList: return file_list
	if isString :
	  import re
	  file_list_converted = re.sub(r'\s', '', file_list).split(',') #remove all whitespaces
	  return file_list_converted
      
    def item(self, item_name):
	"""
	Get one item from cfg dictionary. If the item is nested, then use "." to set the path to the item.
	"""
	self.log.info('Not implemented yet... Sorry!')
	pass
      
    def get_dict(self):
	"""
	Returns the full configuration red 
	"""
	self.log.debug('Getting dictionary from config files: %s', str(self.file_list))
	for cfg_file in self.file_list:
	    """
	    We want to append dictionaries from all the config files.
	    """
	    if self.cfg_type == None: self.cfg_type = self._get_cfg_type(cfg_file)
	    self.log.debug('Updating dictionary from config file in the order provided: %s',str(cfg_file) )
	    if self.cfg_type.lower() in ['yaml', "yml"]: self._get_dict_yaml(cfg_file)
	    elif self.cfg_type.lower() == 'xml': self._get_dict_xml(cfg_file)
	    elif self.cfg_type.lower() == 'json': self._get_dict_json(cfg_file)
	    elif self.cfg_type.lower() == 'ini': self._get_dict_ini(cfg_file)
	    
	return self.cfg_dict
	
    def _get_dict_yaml(self,file_name):
	self.log.debug('Reading yaml configuration and updating dictionary.')  
	import yaml
	fd = open(file_name, 'r')
	self.cfg_dict.update(yaml.load(fd))
	fd.close()
    
    def _get_dict_xml(self,file_name):
        print "Not implemented yet."
        pass

    def _get_dict_json(self,file_name):
        print "Not implemented yet."
        pass

	
    def _get_dict_ini(self,file_name):
        print "Not implemented yet."
        pass

 