#!/bin/env python

import os

import ROOT
import sys
import optparse
import pdb

def probe2latino(chans,ifile,ofile,dryrun=False):

    ROOT.gROOT.cd()
#     ROOT.gDebug=10
    latino = ROOT.TChain('probe_tree')
    for chan in chans:
        
        tree2add = ifile+'/'+chan+'Tree/probe_tree'
        latino.Add(tree2add)
        print ' - ',latino.GetEntries(),tree2add

    if dryrun:
        return

    latino.SetName('latino')
    latino.SetTitle('latino')
    latino.Merge(ofile)
    latino.SetDirectory(0x0)


from optparse import OptionParser
usage='''
 %prog

 Converts step3 trees into standard latino files

'''

parser = OptionParser(usage=usage)
parser.add_option('-d', dest='odir', help='output directory')
parser.add_option('-a', '--all', dest='best', action='store_false', help='Store all candidates (only the best by default)', default=True)
parser.add_option('-n', dest='dry',  action='store_true',  help='output directory',          default=False)

(opts, args) = parser.parse_args()


filenames = args
odir=''


chans=['ellell'] if opts.best else ['elel','elmu','muel','mumu']

if opts.odir:
    odir = opts.odir
    os.system('mkdir -p '+odir)

for f in filenames:
    ofile = os.path.join(odir,'latino'+os.path.basename(f).replace('tree_',''))

    print ofile

    probe2latino(chans,f,ofile,opts.dry)

print '...and DONE'




