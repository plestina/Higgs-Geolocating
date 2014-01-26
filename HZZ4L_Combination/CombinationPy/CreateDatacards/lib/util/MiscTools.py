import sys, os, pwd, commands
import optparse, shlex, re
import math
import errno
import shutil
from array import array


##define function for processing of os command
def processCmd(cmd):
    print cmd
    status, output = commands.getstatusoutput(cmd)
    print output
    if status !=0:
        print 'Error in processing command:\n   ['+cmd+'] \nExiting...'
        #print output
        sys.exit()


def make_sure_path_exists(path) :
    try:
        os.makedirs(path)
    except OSError, exception:
        if exception.errno != errno.EEXIST:
            raise
        
def grep(what, where_list):       
    import re, sys, glob
    #for arg in sys.argv[2:]:
    where_list = [f.strip() for f in where_list.split(",")]
    #print "----- Grep for what={0} and where={1}".format(what, where_list)
    grep_list=[]
    for arg in where_list:
        #print glob.glob(arg.strip())
        for file in glob.glob(arg.strip()):
            for line in open(file, 'r'):
                    if re.search(what, line):
                        grep_list.append(line.strip())        
                        #print "line"line,
    return grep_list