#!/usr/bin/env python

import WWAnalysis.AnalysisStep.scaleFactors_cff
from WWAnalysis.AnalysisStep.stepOneDatasets_cff import stepOneDatasets
from WWAnalysis.AnalysisStep.stepTwoDatasets_cff import stepTwoDatasets
from optparse import OptionParser
from glob import glob
import re, os, os.path
import pdb

parser = OptionParser(usage="%prog [options] block [base path (to switch to pattern)]")
parser.add_option("-p", "--pattern",     dest="pattern",  default='%(id)s.%(name)s/*.root', help="Pattern for file names")
parser.add_option("-P", "--base-path",   dest="basepath", help="base path for searching for files manually (mutually exclusive to -c)")
parser.add_option("-q", "--query",       dest="query",    default='find file where dataset = %(dataset)s', help="Pattern for file names")
parser.add_option("-d", "--db",          dest="db",       default='usePhys', help="using local physics db, usePhys (=2) of usePhys1");
parser.add_option("-j", "--json",        dest="json",     default='certifiedLatinos.52X',  help="JSON file, under WWAnalysis/Misc/Jsons")
parser.add_option("-n", "--num",         dest="num",      type='int',  default=1, help="number of file per job")
parser.add_option("-c", "--crab",        dest="crab",     action='store_true',  help="setup a folder for each dataset ready for crab submission (mutually exclusive to -P)")
parser.add_option("-2", "--two",         dest="two",      action='store_true',  help="set it up to also run step 2")
parser.add_option("-s", "--sched",       dest="sched",    default='glidein', help="change scheduler (default = glidein)")
parser.add_option("-S", "--serv",        dest="serv",     action='store_false',  default=True, help="turn server off")
parser.add_option("-e", "--events-per",  dest="eventsper",  type='int',  default=100000, help="number of events per job")
parser.add_option("-l", "--lumis-per",   dest="lumisper",   type='int',  default=1000,   help="number of lumis per job")
parser.add_option("-a", "--additional",  dest="additional",  default='',  help="additional instructions to step3")
parser.add_option('-t', '--tag',         dest='tag',      default='', help='additional tag for the crab folders')

(options, args) = parser.parse_args()

# need at least a sample to act on
if len(args) < 1:
    raise RuntimeError, "hey, I need a sample to generate for you ..." 

#mutually exclusive
if options.basepath != None and options.crab: 
    raise RuntimeError, "base-path and crab are mutually exclusive!!!" 

# check if the argument is in the scale factor file, exit otherwise
if not hasattr(WWAnalysis.AnalysisStep.scaleFactors_cff, args[0]):
    print "Sample block %s not found in scaleFactors_cff" % args[0]
    for k in dir(WWAnalysis.AnalysisStep.scaleFactors_cff) :
        if "__" not in k: print "\t",k
    raise RuntimeError, "Sample block %s not found in scaleFactors_cff" % args[0]

# check if user wants to search for files manually instead of looking them up 
# with DBS, need a valid path
if options.basepath and not os.path.exists(options.basepath): 
    raise RuntimeError, "Path %s does not exist" % options.basepath

dictToUse = stepTwoDatasets if not options.two else stepOneDatasets
for id,list in getattr(WWAnalysis.AnalysisStep.scaleFactors_cff, args[0]).items(): 
    idn = re.sub('[^0-9]','',id)
    nick = list[0]+('_'+options.tag if options.tag else '')
    print "ID %s: name %s" % (id, nick); 
    arg3  = "scale="+str(list[1]) if len(list) == 2 else "json="+options.json
    additional = options.additional
#     print additional
    db = options.db
    numDB = 1
    if db == "usePhys": 
        numDB = 2

    if options.basepath: 
        pattern = options.basepath+"/"+(options.pattern % {'name':nick, 'id':id})
        os.system("cmsSplit.pl step3.py label=%(dataset)s id=%(id)s two=%(two)s %(additional)s %(arg3)s -a --bash --files=%(pattern)s --label=%(id2)s_%(dataset)s --fj %(num)d" % {
                    'dataset':nick, 'id':idn, 'additional':additional, 'arg3':arg3, 'pattern':pattern, 'num':options.num, 'id2':id, 'two':'True' if options.two else 'False'
                  })
    elif options.crab: 
        query = options.query % {'dataset':dictToUse[id]}
        os.system("mkdir -p %s.%s" % (id,nick) )
        #generate proper config file
        os.system('cp -a ' + os.getenv("CMSSW_BASE")+'/src/WWAnalysis/AnalysisStep/test/step3/step3.py %s.%s/%s.py' % (id,nick,nick))
        #generate the crab.cfg file
        crabFile = open('%s.%s/crab.cfg' % (id,nick), 'w')
        print >> crabFile, """
[CRAB]
jobtype                    = cmssw
scheduler                  = %(sched)s
use_server                 = %(server)d
                          
[CMSSW]                   
datasetpath                = %(dataset)s
pset                       = %(name)s.py
%(isData)s_per_job         = %(eventsper)d
total_number_of_%(isData)s = -1
output_file                = tree_%(id)s_%(name)s.root
dbs_url                    = http://cmsdbsprod.cern.ch/cms_dbs_ph_analysis_0%(numDB)s/servlet/DBSServlet
%(lumimask)s

pycfg_params               = print id=%(id)s label=%(name)s two=True %(arg3)s outputFile=tree_%(id)s_%(name)s.root %(additional)s 
                                                    
[USER]                   

ui_working_dir             = ./crab_0_S3_ID%(id)s_%(name)s
return_data                = 1
# copy_data                  = 0
# storage_element            = T2_US_UCSD
# local_stage_out            = 0
# publish_data               = 0
# publish_data_name          = R42X_%(hwwtag)s_ID%(id)s_%(name)s
# dbs_url_for_publication    = https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_0%(numDB)s_writer/servlet/DBSServlet
        """ % { 
            "sched": options.sched,
            "server": 1 if options.serv else 0,
            "dataset": dictToUse[id],
            "name": nick,
            "numDB": numDB,
            "id": id,
            "arg3": arg3,
            'additional':additional,
            "eventsper": options.eventsper if len(list) ==2 else options.lumisper,
            "isData": 'events' if len(list) == 2 else 'lumis',
            "hwwtag": os.popen("showtags | grep -e \<WWAnalysis\> | head -n 1 | awk '{print $1}'","r").read().strip(),
            "lumimask": "lumi_mask = %s" % os.getenv("CMSSW_BASE")+"/src/WWAnalysis/Misc/Jsons/"+options.json+".json" if len(list)!=2 else ""
        }
        crabFile.close()
    else:
        query = options.query % {'dataset':dictToUse[id]}
        os.system('cmsSplit.pl step3.py label=%(dataset)s %(id)s two=%(two)s %(arg3)s -a --bash --usePhys --dbsql="%(query)s" --label=%(id2)s_%(dataset)s --fj %(num)d' % {
                'dataset':nick, 'id':idn, 'arg3':arg3, 'query':query, 'num':options.num, 'id2':id, 'two':'True' if options.two else 'False'
              })
