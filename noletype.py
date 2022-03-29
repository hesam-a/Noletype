#################################################################################
#
#
# Title: Noletype
# Description: Ionic typer for the polarizable AMOEBA force field
#
# Copyright:            Copyright (c) Hesam Arabzadeh, Chengwen Liu
#           Orlando Acevedo, Pengyu Ren, Wei Yang, Thomas Albrecht-Schoenzart 2022
#
# Noletype is a free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License version 3
# as published by the Free Software Foundation.
#
# Noletype is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR 
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# "Nole" part of the name of the package is from Florida State Universty that is 
# well-known to "Seminoles" and where the package is written.
#
#
#################################################################################

import warnings
import traceback
import os
import sys
import subprocess
import openbabel
import time
import copy
from packaging import version
import numpy as np
import itertools
import smtplib
import textwrap
import matplotlib.pyplot as plt

class NoleType():
    def __init__(self,ionname=None,charge=None,inputmoleculefolderpaths=None,qmpesonly=False,freq=False,optmaxcycle=500,useqmcomputedpes=False,deleteallnonqmfiles=False,ionqmbasisset='DEF2-TZVP',moleculeqmbasisset='DEF2-TZVP',relativistic=False, dkh=False,zora=False,picturechange=False,qmoptmethod='MP2',qmpesmethod='CCSD(T)',waterdimer=True,imiddimer=True,acetmdimer=True, use_psi4=True,use_orca=False,templatekeyfile=None,prmlib=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/ParameterFiles/'+'amoeba09.prm',basissetpath=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+'/'+'BasisSets/',ionxyz =' 0.000000    0.000000    0.000000',poltypepath=os.path.abspath(os.path.split(__file__)[0]),printoutput=False,noltypeini=True,prmstartidx=601,numproc="4",maxmem="4GB",tinkerdir=None,scratchdir="/scratch",paramhead=os.path.abspath(os.path.join(os.path.split(__file__)[0] , os.pardir))+ "/ParameterFiles/amoebabio18_header.prm",,minimizeexe='minimize.x',analyzeexe='analyze.x',Hartree2kcal_mol=627.5095,bashrcpath=None,amoebabioprmpath=None,helpfile='README.md',versionfile='version.md',sleeptime=30,hftol=1e-7, mp2zsolver='pople',mp2zmaxiter='1000', mp2ztol=1e-5,usewf=True,mdcimaxiter=1000,mdcimaxdiis=7,mdcilshift=0.3,mdcistol=2.5e-4,tinkerpolareps=1e-5,polarbound=[1, 2.5],polarweight=1.0, tholebound=[0.0, 0.39],tholeweight=1.0, vdwradiusbound=[2.0, 4.5],vdwradiusweight=1.0,vdwwellweight=1.0,vdwwellbound=[0.6,0.9],minimizationmethod='SLSQP', minimizationmaxiter=100,minimizationtol=1e-6):
        self.ionname = ionname;
        self.charge  = charge;
        self.inputmoleculefolderpaths = inputmoleculefolderpaths;
        self.qmpesonly = qmpesonly;
        self.freq = freq;
        self.optmaxcycle = optmaxcycle;
        self.useqmcomputedpes = useqmcomputedpes;
        self.deleteallnonqmfiles = deleteallnonqmfiles;
        self.ionqmbasisset = ionqmbasisset;
        self.moleculeqmbasisset = moleculeqmbasisset;
        self.relativistic = relativistic;
        self.dkh = dkh;
        self.zora = zora;
        self.picturechange= picturechange;
        self.qmoptmethod = qmoptmethod;
        self.qmpesmethod = qmespmethod;
        self.waterdimer = waterdimer;
        self.imiddimer = imiddimer;
        self.acetmdimer = acetmdimer;
        self.use_psi4 = use_psi4;
        self.use_orca = use_orca
        self.templatekeyfile = templatekeyfile;
        self.prmlib = prmlib;
        self.basissetpath = basissetpath;
        self.ionxyz = ionxyz;
        self.poltypepath = poltypepath;
        self.printoutput = printoutput;
        self.noltypeini = noltypeini;
        self.prmstartidx = prmstartidx;
        self.numproc = numproc;
        self.maxmem = maxmem;
        self.tinkerdir = tinkerdir;
        self.scratchdir = scratchdir;
        self.paramhead = paramhead;
        self.bashrcpath = bashrcpath;
        self.amoebabioprmpath = amoebabioprmpath;
        self.helpfile = helpfile;
        self.versionfile = versionfile;
        self.sleeptime = sleeptime;
        self.hftol = hftol;
        self.mp2zsolver = mp2zsolver;
        self.mp2zmaxiter = mp2zmaxiter;
        self.mp2ztol = mp2ztol;
        self.usewf = usewf;
        self.mdcimaxiter = mdcimaxiter;
        self.mdcimaxdiis = mdcimaxdiis;
        self.mdcilshift = mdcilshift;
        self.mdcistol = mdcistol;
        self.tinkerpolareps = tinkerpolareps;
        self.polarbound = polarbound;
        self.polarweight = polarweight;
        self.tholebound = tholebound;
        self.tholeweight = tholeweight;
        self.vdwradiusbound = vdwradiusbound;
        self.vdwradiusweight = vdwradiusweight;
        self.vdwwellweight = vdwwellweight;
        self.vdwwellbound = vdwwellbound;
        self.minimizationmethod = minimizationmethod; 
        self.minimizationmaxiter = minimizationmaxiter;
        self.minimizationtol = minimizationtol;

        opts, xargs = getopt.getopt(sys.argv[1:],'h',["help"])

        for i, j in opts:
            if o in ("-h", "--help"):
                self.copyright()
                self.usage()
                sys.exit(2)


if __name__ == '__main__':
    def RunNoltype():
        noletype=NoleType()
        try:
            noletype.main()
        except:
            traceback.print_exc(file=sys.stdout)
            text = str(traceback.format_exc())
            #if os.path.exists(poltype.scrtmpdirgau):
            #    shutil.rmtree(poltype.scrtmpdirgau)
            #if os.path.exists(poltype.scrtmpdirpsi4):
            #    shutil.rmtree(poltype.scrtmpdirpsi4)
            #if poltype.email!=None:
            #    password='amoebaisbest'
            #    fromaddr = 'poltypecrashreportnoreply@gmail.com'
            #    toaddr = poltype.email
            #    filename=poltype.logfname
            #    poltype.WriteToLog(text)
            #    poltype.WriteToLog('Poltype has crashed!')
            #    try:
            #        poltype.SendCrashReportEmail(text,fromaddr,toaddr,password,filename)
            #    except:
            #        pass
            raise ValueError('Houston, we have a problem. Buy a developer some coffee!')
    RunNoletype()
