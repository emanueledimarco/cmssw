from WWAnalysis.AnalysisStep.mcAnalysis import *
import sys

# class TableMaker:
#     def __init__(self,samples,options):
# #         super(TableMaker,self).__init__(samples,options)
#         self._options = options
#         self._allData = {}
#         self._ncol = options.ncol
#         self._fOut = open(options.tex.replace('%m',str(options.mass)).replace('%l',str(int(options.lumi*1000))),'w') if options.tex!= "" else sys.stdout 
#         self._order   = []
#         for line in open(samples,'r'):
#             field = line.split(':')
#             rootfile = options.path+"/"+"tree_%s.root" % field[1].strip()
#             if len(field) != 2 or field[0][0] == '#': continue
#             fieldName = field[0].replace('%m',str(options.mass))
#             if fieldName in self._order:                               self._allData[fieldName].append(TreeToYield(rootfile % options.mass if rootfile.find('%d') != -1 else rootfile,options))
#             else                      : self._order.append(fieldName); self._allData[fieldName] =     [TreeToYield(rootfile % options.mass if rootfile.find('%d') != -1 else rootfile,options)]
#     def getYields(self,cuts):
#         ret = {}
#         for key in self._allData:
#             ret[key] = self._getYields(self._allData[key],cuts)
#         return ret
#     def printTable(self,report):
#         self._printHeader()
#         if self._options.final: self._printSummary(report)
#         else:                   self._printBreakdown(report)
#         self._printFooter()
#     def _printHeader(self):
#         print  >> self._fOut, "\\begin{table}[h!]\\begin{center}"
#         print  >> self._fOut, "\\caption{{ {0} \\label{{ {1} }} }}".format(self._options.caption,self._options.tex.split('.')[0] if self._options.tex!= "" else "labelME" )
#         self._fOut.write("\\%s{\\begin{tabular}{|c|" %self._options.size)
#         for i in range(self._ncol) if self._ncol != 0 else self._order: self._fOut.write("c|")
#         print  >> self._fOut, "}"
#     def _printFooter(self):
#         print  >> self._fOut, "\\end{tabular}}\\end{center}\\end{table}"
#     def _printSummary(self,report):
#         S = len(self._order)                          #nSamples
#         C = self._ncol if self._ncol != 0 else S+1    #nCol
#         R = S / (C-1) + ((S%(C-1))>0)                 #n "row";  well nRow is really R*(nCuts+1)
#         for r in range(R):
#             print >> self._fOut, "\\hline"
#             for sample in self._order[r*(C-1):min((r+1)*(C-1),S)]: 
#                 print >> self._fOut, " & {0}".format(sample),
#             print >> self._fOut, " \\\\ \\hline"
#             for i in range(len(report[self._order[0]][0][1])):
#                 if report[self._order[0]][0][1][i][0] == "all" : print >> self._fOut, "\\hline"
#                 print >> self._fOut, report[self._order[0]][0][1][i][0].replace('el','e').replace('mu','$\\mu$'),
#                 for sample in self._order[r*(C-1):min((r+1)*(C-1),S)]: 
#                     fmt = " & ${0:.2f}\\pm{1:.2f}$"
#                     if report[sample][-1][1][i][1] >  100: fmt = " & ${0:.1f}\\pm{1:.1f}$"
#                     if report[sample][-1][1][i][1] > 1000: fmt = " & ${0:.0f}\\pm{1:.0f}$"
#                     if sample == 'data': fmt = " & ${0:.0f}$" 
#                     print >> self._fOut, fmt.format(report[sample][-1][1][i][1],report[sample][-1][1][i][2]),
#                 print >> self._fOut, "\\\\"
#             print >> self._fOut, "\\hline"
#     def _printBreakdown(self,report):
#         S = len(self._order)                          #nSamples
#         C = self._ncol if self._ncol != 0 else S+1    #nCol
#         R = S / (C-1) + ((S%(C-1))>0)                 #n "row";  well nRow is really R*(nCuts+1)
#         for r in range(R):
#             print >> self._fOut, "\\hline"
#             for sample in self._order[r*(C-1):min((r+1)*(C-1),S)]: 
#                 print >> self._fOut, " & {0}".format(sample),
#             print >> self._fOut, " \\\\ \\hline"
#             for i in range(len(report[self._order[0]])):
#                 if report[self._order[0]][i][0] == "entry point": continue
#                 print >> self._fOut, report[self._order[0]][i][0],
#                 for sample in self._order[r*(C-1):min((r+1)*(C-1),S)]: 
#                     for channel in report[sample][i][1]: 
#                         fmt = " & ${0:.1f}\\pm{1:.1f}$"
#                         if channel[1] > 1000: fmt = " & ${0:.0f}\\pm{1:.0f}$"
#                         if sample == 'data': fmt = " & ${0:.0f}$" 
#                         if channel[0] == self._options.hypo: print >> self._fOut, fmt.format(channel[1],channel[2]),
#                 print >> self._fOut, "\\\\"
#             print >> self._fOut, "\\hline"
#     def _getYields(self,ttylist,cuts):
#         return mergeReports([tty.getYields(cuts) for tty in ttylist])
        
class TableMaker(MCAnalysis):
    def __init__(self,samples,options):
        MCAnalysis.__init__(self,samples,options)
        self._ncol = options.ncol
        self._fOut = open(options.tex.replace('%m',str(options.mass)).replace('%l',str(int(options.lumi*1000))),'w') if options.tex!= "" else sys.stdout 
        self._twikiOut = open(options.twiki.replace('%m',str(options.mass)).replace('%l',str(int(options.lumi*1000))),'w') if options.twiki != "" else None
        self._order   = []
        for line in open(samples,'r'):
            if len(line.strip()) == 0 or line.strip()[0] == '#': continue
            field = [f.strip() for f in line.split(':')]
            if ("%d" in field[0]): field[0] = field[0] % options.mass
            if field[0] not in self._order: self._order.append(field[0]) 
    def printTable(self,report):
        self._printHeader()
        if self._options.final: self._printSummary(report)
        else:                   self._printBreakdown(report)
        self._printFooter()
    def _printHeader(self):
        print  >> self._fOut, "\\begin{table}[h!]\\begin{center}"
        print  >> self._fOut, "\\caption{{ {0} \\label{{ {1} }} }}".format(self._options.caption,self._options.tex.split('.')[0] if self._options.tex!= "" else "labelME" )
        self._fOut.write("\\%s{\\begin{tabular}{|c|" %self._options.size)
        for i in range(self._ncol) if self._ncol != 0 else self._order: self._fOut.write("c|")
        print  >> self._fOut, "}"
        if self._twikiOut != None:
            print  >> self._twikiOut, "<!-- \n{0} \n{1} \n-->".format(self._options.caption,self._options.twiki.split('.')[0])
    def _printFooter(self):
        print  >> self._fOut, "\\end{tabular}}\\end{center}\\end{table}"
    def _printSummary(self,report):
        S = len(self._order)                          #nSamples
        C = self._ncol if self._ncol != 0 else S+1    #nCol
        R = S / (C-1) + ((S%(C-1))>0)                 #n "row";  well nRow is really R*(nCuts+1)
        for r in range(R):
            print >> self._fOut, "\\hline"
            for sample in self._order[r*(C-1):min((r+1)*(C-1),S)]: 
                print >> self._fOut, " & {0}".format(sample),
            print >> self._fOut, " \\\\ \\hline"
            for i in range(len(report[self._order[0]][0][1])):
                if report[self._order[0]][0][1][i][0] == "all" : print >> self._fOut, "\\hline"
                print >> self._fOut, report[self._order[0]][0][1][i][0].replace('el','e').replace('mu','$\\mu$'),
                for sample in self._order[r*(C-1):min((r+1)*(C-1),S)]: 
                    fmt = " & ${0:.2f}\\pm{1:.2f}$"
                    if report[sample][-1][1][i][1] >  100: fmt = " & ${0:.1f}\\pm{1:.1f}$"
                    if report[sample][-1][1][i][1] > 1000: fmt = " & ${0:.0f}\\pm{1:.0f}$"
                    if 'data' in sample: fmt = " & ${0:.0f}$" 
                    print >> self._fOut, fmt.format(report[sample][-1][1][i][1],report[sample][-1][1][i][2]),
                print >> self._fOut, "\\\\"
            print >> self._fOut, "\\hline"
        if self._twikiOut != None:
            for r in range(R):
                print >> self._twikiOut, "|  *chan*  ",
                for sample in self._order[r*(C-1):min((r+1)*(C-1),S)]: 
                    print >> self._twikiOut, "|  *{0}*  ".format(sample),
                print >> self._twikiOut, "|"
                for i in range(len(report[self._order[0]][0][1])):
                    print >> self._twikiOut, "|  *{0}*  ".format(report[self._order[0]][0][1][i][0].replace('el','e').replace('mu','&mu;')),
                    for sample in self._order[r*(C-1):min((r+1)*(C-1),S)]: 
                        fmt = "|  {0:.2f}&plusmn;{1:.2f}  "
                        if report[sample][-1][1][i][1] >  100: fmt = "|  {0:.1f}&plusmn;{1:.1f}  "
                        if report[sample][-1][1][i][1] > 1000: fmt = "|  {0:.0f}&plusmn;{1:.0f}  "
                        if 'data' in sample: fmt = "|  {0:.0f}  " 
                        print >> self._twikiOut, fmt.format(report[sample][-1][1][i][1],report[sample][-1][1][i][2]),
                    print >> self._twikiOut, "|"
    def _printBreakdown(self,report):
        S = len(self._order)                          #nSamples
        C = self._ncol if self._ncol != 0 else S+1    #nCol
        R = S / (C-1) + ((S%(C-1))>0)                 #n "row";  well nRow is really R*(nCuts+1)
        for r in range(R):
            print >> self._fOut, "\\hline"
            for sample in self._order[r*(C-1):min((r+1)*(C-1),S)]: 
                print >> self._fOut, " & {0}".format(sample),
            print >> self._fOut, " \\\\ \\hline"
            for i in range(len(report[self._order[0]])):
                if report[self._order[0]][i][0] == "entry point": continue
                print >> self._fOut, report[self._order[0]][i][0],
                for sample in self._order[r*(C-1):min((r+1)*(C-1),S)]: 
                    for channel in report[sample][i][1]: 
                        fmt = " & ${0:.1f}\\pm{1:.1f}$"
                        if channel[1] > 1000: fmt = " & ${0:.0f}\\pm{1:.0f}$"
                        if sample == 'data': fmt = " & ${0:.0f}$" 
                        if channel[0] == self._options.hypo: print >> self._fOut, fmt.format(channel[1],channel[2]),
                print >> self._fOut, "\\\\"
            print >> self._fOut, "\\hline"
        
if __name__ == "__main__":
    from optparse import OptionParser
    
    parser = OptionParser(usage="%prog [options] mcList.txt cuts.txt")
    addMCAnalysisOptions(parser)

    parser.add_option("-c", "--ncol",         dest="ncol",    type="int",    default="0",           help="number of columns in the table, gonna try to be smart here (default == nSamples)") 
    parser.add_option("-C", "--caption",      dest="caption", type="string", default="Caption Me!", help="caption at the top of the table") 
    parser.add_option("-H", "--hypo",         dest="hypo",    type="string", default="all",         help="which mode to use (ignored if -s)") 
    parser.add_option("-L", "--latex",        dest="tex",     type="string", default="",            help="tex out name (can will sub %m with mass, %l with lumi)")
    parser.add_option("-T", "--twiki",        dest="twiki",   type="string", default="",            help="twiki out name (can will sub %m with mass, %l with lumi)")
    parser.add_option("--size",               dest="size",    type="string", default="normalsize",  help="text size")


    (options, args) = parser.parse_args()
    tm = TableMaker(args[0],options)
    cf  = CutsFile(args[1],options)
    for cutFile in args[2:]:
        temp = CutsFile(cutFile,options)
        for cut in temp.cuts():
            cf.add(cut[0],cut[1])
    report = tm.getYields(cf)
    tm.printTable(report)
