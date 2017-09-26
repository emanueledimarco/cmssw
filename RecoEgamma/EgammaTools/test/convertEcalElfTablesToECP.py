#!/usr/bin/env python
# usage: python convertEcalElfTablesToECP.py Moriond17_23Jan_ele_scales.dat test
import sys, os.path

infile, outfile = sys.argv[1], sys.argv[2]
outf = {'EB':open(outfile+"EB.txt",'w'), 'EE':open(outfile+"EE.txt",'w')}
lines = {'EB':0,'EE':0}
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
    if lines[subdet]==0:
        outf[subdet].write(('{%d' % len(bins))+''.join(['%20s' % v for v in variables])+
                           ('%10d' % len(pars))+'          '+'         '.join(['par_%d' % p for p in xrange(len(pars))])+'}\n')
    for v in variables:
        r = ranges[v]
        outf[subdet].write('%10s%10s' % (r[0],r[1])+'    ')
    outf[subdet].write(''.join(['%10s' % p for p in pars])+'\n')
    lines[subdet] += 1

print "DONE. Output written into ",outfile
