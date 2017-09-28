#!/usr/bin/env python
import sys, os.path
inputdir=sys.argv[1]
cfgfile=sys.argv[2]

if len(sys.argv)<3:
    print "Need at least write_scalecorr_iovs.py <inputdir> <cfgfile>"
    exit()

sqlitefile = "testEPC" if len(sys.argv)==3 else sys.argv[3]
tag="Ele_Scales_PromptReco_Moriond17" if len(sys.argv)==4 else sys.argv[4]

sqlitefile=sqlitefile.split('.')[0]

import glob
files = (glob.glob(inputdir+"/scales_*.txt"))

if os.path.isfile(sqlitefile+'_EB.db') or os.path.isfile(sqlitefile+'_EE.db'):
    print "The output DB files ",sqlitefile," for EB and EE exist. This has to be removed in order not to append spurious IOVs.\n"
    print '\nDo you agree? [y/N]\n'
    if raw_input()!='y':
        print 'Aborting'
        exit()
    else: os.system('rm '+sqlitefile+'_EB.db '+sqlitefile+'_EE.db')

for f in files:
    print "Importing txt file: ",f," into the DB..."
    rmin,rmax=((os.path.basename(f)).split('.')[0]).split('_')[2:]
    subdet='EB' if 'EB' in f else 'EE'
    tagdet = tag+'_'+subdet
    dbfile=sqlitefile+"_"+subdet+".db"
    cmd="cmsRun "+cfgfile+" "+rmin+" "+f+" "+dbfile+" "+tagdet
    print cmd
    os.system(cmd)
