#!/usr/bin/python
import sys
import os

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-n", "--ntoys", dest="ntoys", default=10, help="number of toys to generate")
(options, args) = parser.parse_args()

from ROOT import gROOT
from ROOT import gSystem
gSystem.Load("libRooFit");
gSystem.Load("../../../MLFit/workdir/libMLFit.so");
gROOT.LoadMacro('toyMC_OoT.cc')
from ROOT import Generate
Generate(int(options.ntoys))
