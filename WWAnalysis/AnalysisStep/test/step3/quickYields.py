#!/usr/bin/env python
from WWAnalysis.AnalysisStep.mcAnalysis import *

def sorted(x):
    all = [i for i in x]
    all.sort()
    return all

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] mca.txt process cuts%%d.txt")
    addMCAnalysisOptions(parser)
    parser.add_option("--jets",         dest="jets",     type="int",    default=0,            help="Jet bin");
    parser.add_option("--masses",       dest="massfile", type="string", default="masses",     help="File containing masses to consider");
    parser.add_option("--scale",        dest="scale",    type="float",  default=(0.0, 0.0),   nargs=2, help="Scale factor and _absolute_ uncertainty on it");
    parser.add_option("--gamma",        dest="gamma",    type="int",    default=0,           help="Extract out from the scale uncertainty a gamma factor with the specified number of observed events");
    parser.add_option("-v",             dest="verbose",  action="store_true", help="Print out report")
    parser.remove_option("--out");
    parser.add_option("-o", "--out",    dest="out",     default="{process}Card_{ll}_{jets}j.txt",  help="Output file pattern; ll = (ee,mm,...); channel = (elel, mumu,...)");
    (options, args) = parser.parse_args()
    if not options.verbose: options.final = True
#     (mcafile, process, cutspattern) = args
    mcafile = args[0]
    process = args[1]
    cutspattern = args[2:]
    print cutspattern

    tty = MCAnalysis(mcafile,options)
    if options.scale[0] != 0: tty.scaleProcess(process, 1.0)

    if options.gamma:
        if options.scale[0] == 0: raise RuntimeError, "Can't use gamma with no scale"
        if options.scale[1] == 0: raise RuntimeError, "Can't use gamma with no scale uncertainty"
        gammaRelErr2 = 1.0/(options.gamma+1.)
        scaleRelErr2 = pow(options.scale[1]/options.scale[0], 2)
        if scaleRelErr2 < gammaRelErr2: 
            raise RuntimeError, "Relative uncertainty on scale %f doesn't cover the stat. uncertainty from the gamma %f" % (sqrt(scaleRelErr2), sqrt(gammaRelErr2))
        options.scale = ( options.scale[0], options.scale[0]*sqrt(scaleRelErr2-gammaRelErr2) )

    reports = {}; channels = None
    for line in open(options.massfile, "r"):
        mass = int(line)
        print "Running mass ",mass," from ", cutspattern[0] % mass if '%d' in cutspattern[0] else cutspattern[0]
        cf  = CutsFile(cutspattern[0] % mass if '%d' in cutspattern[0] else cutspattern[0], options)
        for cutFile in cutspattern[1:]:
            temp = CutsFile(cutFile%mass if '%d' in cutFile else cutFile,options)
            for cut in temp.cuts():
                cf.add(cut[0],cut[1])
        report = tty.getYields(cf, process=process)
        if options.verbose: tty.prettyPrint(report)
        reports[mass] = report[process][-1][1]
        channels = [ n for n,y,e in reports[mass] ]

    for ik,K in enumerate(channels):
        if K == "all": continue
        ll = K if K == "of" or K == "sf" else (K[0]+K[2])
        outname = options.out.format(process=process, channel=K, ll=ll, jets=options.jets)
        print "Writing to ",outname
        out = open(outname, 'w')
        if options.gamma: out.write("#m(H) N alpha err(alpha)\n")
        else:             out.write("#m(H) yield err(yield)\n")
        for m in sorted(reports.keys()):
            y0, e0 = reports[m][ik][1:]
            if options.scale[0] != 0:
                y = y0 * options.scale[0]
                e = hypot(e0*options.scale[0], y0*options.scale[1])
            else: 
                y,e = y0,e0
            if options.gamma:
                a = y/options.gamma; ea = e/options.gamma
                out.write("%3d   %3d  %6.4f  %6.4f\n" % (m,options.gamma,a,ea))
                print "    %3d   %3d  %6.4f  %6.4f  " % (m,options.gamma,a,ea)
            else:
                out.write("%3d   %6.3f  %6.3f\n" % (m,y,e))
                print "    %3d   %6.3f  %6.3f  " % (m,y,e)
        print
