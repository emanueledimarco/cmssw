#!/usr/bin/env python

import WWAnalysis.AnalysisStep.scaleFactors_cff
from WWAnalysis.AnalysisStep.stepOneDatasets_cff import stepOneDatasets
from WWAnalysis.AnalysisStep.stepTwoDatasets_cff import stepTwoDatasets
from optparse import OptionParser
from glob import glob
import re, os, os.path

parser = OptionParser(usage="%prog [options] block [base path (to switch to pattern)]")
parser.add_option("-p", "--pattern",   dest="pattern",  default='%(id)s.%(name)s/*.root', help="Pattern for file names");
parser.add_option("-P", "--base-path", dest="basepath", help="base path for searching for files manually (mutually exclusive to -c)");
parser.add_option("-q", "--query",     dest="query",    default='find file where dataset = %(dataset)s', help="Pattern for file names");
parser.add_option("-j", "--json",      dest="json",     default='certifiedLatinos.42X',  help="JSON file, under WWAnalysis/Misc/Jsons")
parser.add_option("-n", "--num",       dest="num",      type='int',  default=1, help="number of file per job")
parser.add_option("-c", "--crab",      dest="crab",     action='store_true',  help="setup a folder for each dataset ready for crab submission (mutually exclusive to -P)")
parser.add_option("-2", "--two",       dest="two",      action='store_true',  help="set it up to also run step 2")
parser.add_option("-s", "--sched",     dest="sched",    default='glidein', help="change scheduler (default = glidein)")
parser.add_option("-S", "--serv",      dest="serv",     action='store_false',  default=True, help="turn server off")
parser.add_option("-e", "--events-per",      dest="eventsper",  type='int',  default=100000, help="number of events per job")
parser.add_option("-l", "--lumis-per",       dest="lumisper",   type='int',  default=1000,   help="number of lumis per job")
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
    #print "ID %s: name %s" % (id, list[0]); continue
    arg3  = list[1] if len(list) == 2 else options.json
    if options.basepath: 
        pattern = options.basepath+"/"+(options.pattern % {'name':list[0], 'id':id})
        os.system("cmsSplit.pl step3.py %(dataset)s %(id)s %(arg3)s %(two)s -a --bash --files=%(pattern)s --label=%(id2)s_%(dataset)s --fj %(num)d" % {
                    'dataset':list[0], 'id':idn, 'arg3':arg3, 'pattern':pattern, 'num':options.num, 'id2':id, 'two':'true' if options.two else 'false'
                  })
    elif options.crab: 
        query = options.query % {'dataset':dictToUse[id]}
        os.system("mkdir -p %s.%s" % (id,list[0]) )
        #generate proper config file
        os.system('cmsSplit.pl step3.py %(dataset)s %(id)s %(arg3)s %(two)s -a --bash --usePhys --dbsql="%(query)s" --label=%(id2)s_%(dataset)s --fj 10000' % {
                'dataset':list[0], 'id':idn, 'arg3':arg3, 'query':query, 'id2':id, 'two':'true' if options.two else 'false'
              })
        #generate the crab.cfg file
        crabFile = open('%s.%s/crab.cfg' % (id,list[0]), 'w')
        print >> crabFile, """
[CRAB]
jobtype                    = cmssw
scheduler                  = %(sched)s
use_server                 = %(server)d
                          
[CMSSW]                   
datasetpath                = %(dataset)s
pset                       = %(name)s.py
%(isData)s_per_job             = %(eventsper)d
total_number_of_%(isData)s     = -1
output_file                = tree_%(id)s_%(name)s_job1.root
dbs_url                    = http://cmsdbsprod.cern.ch/cms_dbs_ph_analysis_02/servlet/DBSServlet
%(lumimask)s
                          
[USER]                    
ui_working_dir             = ./crab_0_S3_ID%(id)s_%(name)s
return_data                = 1
copy_data                  = 0
storage_element            = T2_US_UCSD
local_stage_out            = 0
publish_data               = 0
publish_data_name          = R42X_%(hwwtag)s_ID%(id)s_%(name)s
dbs_url_for_publication    = https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_02_writer/servlet/DBSServlet
        """ % { 
            "sched": options.sched,
            "server": 1 if options.serv else 0,
            "dataset": dictToUse[id],
            "name": list[0],
            "id": id,
            "eventsper": options.eventsper if len(list) ==2 else options.lumisper,
            "isData": 'events' if len(list) == 2 else 'lumis',
            "hwwtag": os.popen("showtags | grep WWAnalysis | head -n 1 | awk '{print $1}'","r").read().strip(),
            "lumimask": "lumi_mask = %s" % os.getenv("CMSSW_BASE")+"/src/WWAnalysis/Misc/Jsons/"+arg3+".json" if len(list)!=2 else ""
        }
        crabFile.close()
        #put them both in the right dir
        os.system("mv step3_%(id)s_%(name)s_job1.py %(id)s.%(name)s/%(name)s.py" % {'id':id,'name':list[0]} )
        #remove the rest of the shit?
        os.system("rm step3_%(id)s_%(name)s*" % {'id':id,'name':list[0]} )
    else:
        query = options.query % {'dataset':dictToUse[id]}
        os.system('cmsSplit.pl step3.py %(dataset)s %(id)s %(arg3)s %(two)s -a --bash --usePhys --dbsql="%(query)s" --label=%(id2)s_%(dataset)s --fj %(num)d' % {
                'dataset':list[0], 'id':idn, 'arg3':arg3, 'query':query, 'num':options.num, 'id2':id, 'two':'true' if options.two else 'false'
              })
