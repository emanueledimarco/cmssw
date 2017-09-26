#!/usr/bin/env python
# usage: python convertEcalElfTablesToECP.py Moriond17_23Jan_ele_scales.dat test
import sys, os.path

infile = sys.argv[1]
processed=0
for line in open(infile,'r').readlines():
    fields = line.split()
    if len(fields)==0 or (fields[0])[0]=="#": continue
    bins = fields[0].split('-')
    ranges = {}
    for b in bins:
        binvar = b.split('_')[0]
        if any(var in binvar for var in ['absEta','Et']):
            ranges[binvar] = b.split('_')[1:]
        elif any(var in binvar for var in ['bad','lowR9']):
            ranges['r9'] = ['-1','0.94']
        elif any(var in binvar for var in ['gold','highR9']):
            ranges['r9'] = ['0.94','2']
        elif 'gainEle' in binvar:
            ranges['gain'] = [float(b.split('_')[1])-0.5,float(b.split('_')[1])+0.5]
        else: print "ERROR! binvar ",binvar," unknown. Not adding it."
    startParsIdx = 4 if 'runNumber' in fields else 1
    pars = fields[startParsIdx:]
    variables=ranges.keys()
    variables.sort()
    subdet = 'EB' if float(ranges['absEta'][0])<1.5 else 'EE'
    if 'runNumber' in fields:
        outfilename = '_'.join(['scales',subdet,fields[2],fields[3]])+".txt"
    else:
        outfilename = '_'.join(['smearings',subdet])+".txt"
    newIOV = not os.path.isfile(outfilename) 
    if newIOV:
        outf = open(outfilename,'w')
        outf.write(('{%d' % len(bins))+''.join(['%20s' % v for v in variables])+
                   ('%10d' % len(pars))+'          '+'         '.join(['par_%d' % p for p in xrange(len(pars))])+'}\n')
        outf.close()
    outfapp = open(outfilename,'a')
    for v in variables:
        r = ranges[v]
        outfapp.write('%10s%10s' % (r[0],r[1])+'    ')
    outfapp.write(''.join(['%10s' % p for p in pars])+'\n')
    outfapp.close()
    processed += 1
    if processed % 100==0: print "Processed ",processed," lines."

print "DONE."
