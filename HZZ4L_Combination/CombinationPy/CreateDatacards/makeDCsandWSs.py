#!/usr/bin/python
#-----------------------------------------------
# Latest update: 2012.08.30
# by Matt Snowball
#-----------------------------------------------
import sys, os, pwd, commands
import optparse, shlex, re
import math
from ROOT import *
import ROOT
from array import array
from mainClass import *
from inputReader import *

#define function for parsing options
def parseOptions():

    usage = ('usage: %prog [options] datasetList\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)
    
    parser.add_option('-i', '--input', dest='inputDir', type='string', default="",    help='inputs directory')
    parser.add_option('-d', '--is2D',   dest='is2D',       type='int',    default=1,     help='is2D (default:1)')
    parser.add_option('-a', '--append', dest='appendName', type='string', default="",    help='append name for cards dir')
    parser.add_option('-b', action='store_true', dest='noX', default=True ,help='no X11 windows')
    parser.add_option('-t', '--templateDir', type='string', dest='templateDir', default="templates2D" ,help='directory with 2D template histos')
    parser.add_option('-e', '--massError',   dest='massError',       type='int',    default=0,     help='massError (default:0)')
    parser.add_option('-u', '--mekd',   dest='mekd',       type='int',    default=0,     help='mekd double gaussian inputs (default:0)')
    parser.add_option('-s', '--altHypo',   dest='hypothesis', type='string',  default='gg0-',  help='Alt Hypo:gg0-, gg0+, qq1-, qq1+, gg2+m, qq2+m')
    parser.add_option('--unfold', action='store_true', dest='unfold', default=False ,help='unfold 2D to 1D histograms')
    parser.add_option('--unfolded', action='store_true', dest='unfolded', default=False ,help='unfold 2D to 1D histograms but use external unfolding')
    parser.add_option('-L', '--lumi', dest='lumi', type='string', default="UseFromInputFile",    help='lumi for datacards')
    parser.add_option('--terms', dest='termNames', type='string', default="UseDeafaultTerms",    help='pass collection of term names separated by colon ":"')
    parser.add_option('--user_option', dest='user_option', type='string', default="UseFromInputFile",    help='Arbitrary user options for datacards')

    
    
    # store options and arguments as global variables
    global opt, args
    (opt, args) = parser.parse_args()

    if (opt.is2D != 0 and opt.is2D != 1):
        print 'The input '+opt.is2D+' is unkown for is2D.  Please choose 0 or 1. Exiting...'
        sys.exit()

    if (opt.appendName == ''):
        print 'Please pass an append name for the cards directory! Exiting...'
        sys.exit()

    if (opt.inputDir == ''):
        print 'Please pass an input directory! Exiting...'
        sys.exit()


    
# define make directory function
def makeDirectory(subDirName):
    if (not os.path.exists(subDirName)):
        cmd = 'mkdir -p '+subDirName
        status, output = commands.getstatusoutput(cmd)
        if status !=0:
            print 'Error in creating submission dir '+subDirName+'. Exiting...'
            sys.exit()
    else:
        print 'Directory '+subDirName+' already exists. Exiting...'
        sys.exit()


#define function for processing of os command
def processCmd(cmd):
#    print cmd
    status, output = commands.getstatusoutput(cmd)
    if status !=0:
        #print 'Error in processing command:\n   ['+cmd+'] \nExiting...'
        print output
        sys.exit()



def creationLoop(directory):
    global opt, args
    
#    startMass=[ 110.0, 140.0, 160.0, 290.0, 350.0, 400.0, 600.0 ]
#    stepSizes=[ 0.5, 0.5, 2.0, 5.0, 10.0, 20.0, 50.0 ]
#    endVal=[ 60, 40, 65, 12, 5, 10, 9 ]

#    startMass=[ 110.0, 124.5, 126.5, 130.0, 160.0, 290.0, 350.0, 400.0, 600.0 ]
#    stepSizes=[ 0.5, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0 ]
#    endVal=[ 29, 20, 7,  30, 65, 12, 5, 10, 9 ]

#    startMass=[ 400.0, 600.0 ]
#    stepSizes=[ 20.0, 50.0 ]
#    endVal=[      10,   9 ]


    startMass=[ 126.0 ]
    stepSizes=[ 0.5 ]
    endVal=[ 1 ]

    myClass = mainClass()
    
    myReader4e = inputReader(opt.inputDir+"/inputs_4e.txt")
    myReader4e.readInputs()
    theInputs4e = myReader4e.getInputs()
    theInputs4e['unfold'] = opt.unfold
    theInputs4e['unfolded'] = opt.unfolded

    myReader4mu = inputReader(opt.inputDir+"/inputs_4mu.txt")
    myReader4mu.readInputs()
    theInputs4mu = myReader4mu.getInputs()
    theInputs4mu['unfold'] = opt.unfold
    theInputs4mu['unfolded'] = opt.unfolded
    
    myReader2e2mu = inputReader(opt.inputDir+"/inputs_2e2mu.txt")
    myReader2e2mu.readInputs()
    theInputs2e2mu = myReader2e2mu.getInputs()
    theInputs2e2mu['unfold'] = opt.unfold
    theInputs2e2mu['unfolded'] = opt.unfolded
    
      
    if opt.user_option != "UseFromInputFile": 
	theInputs4e['user_option'] = opt.user_option
	theInputs4mu['user_option'] = opt.user_option
	theInputs2e2mu['user_option'] = opt.user_option
	

    
    if opt.lumi != "UseFromInputFile": 
	theInputs4e['lumi'] = float(opt.lumi)
	theInputs4mu['lumi'] = float(opt.lumi)
	theInputs2e2mu['lumi'] = float(opt.lumi)
	

	
    if opt.termNames != "UseDeafaultTerms":
	termNames = opt.termNames.split(",")
	print "@@@@ Using termNames = {0}".format(str(termNames))
	
	theInputs4e['termNames'] = termNames
	theInputs4mu['termNames'] = termNames 
	theInputs2e2mu['termNames'] = termNames
		
    
    
    a=0
    while (a < len(startMass) ):
	
	c = 0
        while (c < endVal[a] ): 
            
            mStart = startMass[a]
            step = stepSizes[a]
            mh = mStart + ( step * c ) 
            mhs = str(mh).replace('.0','')
            makeDirectory(directory+'/HCG/'+mhs)
            makeDirectory(directory+'/HCG_XSxBR/'+mhs)

            print mh
            myClass.makeCardsWorkspaces(mh,directory,theInputs4e,opt.templateDir,opt.massError,opt.is2D,opt.mekd,opt.hypothesis)
            myClass.makeCardsWorkspaces(mh,directory,theInputs4mu,opt.templateDir,opt.massError,opt.is2D,opt.mekd,opt.hypothesis)
            myClass.makeCardsWorkspaces(mh,directory,theInputs2e2mu,opt.templateDir,opt.massError,opt.is2D,opt.mekd,opt.hypothesis)
                          
            c += 1
            

	a += 1






# the main procedure
def makeDCsandWSs():
    
    # parse the arguments and options
    global opt, args
    parseOptions()

    if (opt.appendName != ''):
        dirName = 'cards_'+opt.appendName
    

    subdir = ['HCG','HCG_XSxBR','figs']

    for d in subdir:
        makeDirectory(dirName+'/'+d)
        

    creationLoop(dirName)
    

    sys.exit()



# run the create_RM_cfg() as main()
if __name__ == "__main__":
    makeDCsandWSs()


