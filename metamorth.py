#!/usr/bin/env python                                                                                                                                                                                       
import argparse, os, sys, time,subprocess, signal
sourcedir=os.path.dirname(os.path.abspath(__file__))
cwdir=os.getcwd()
sys.path.append(sourcedir)

from pythonmods import runsubprocess


def default_sigpipe():
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)


parser = argparse.ArgumentParser(description='Run pipeline scripts')
parser.add_argument('-s','--sequences', help='Genbank file (full), containing sequences for all-vs-all pairwise comparison', required=False)
parser.add_argument('-b','--besthits', help='Text file containing best hits or reciprocal best hits', required=False)
parser.add_argument('-o','--out', help='Output directory (required)', required=True)
parser.add_argument('-e','--evalue', help='BLAST e-value cutoff (default: 1e-6)', default=1e-6, type=float)
parser.add_argument('-i','--pident', help='BLAST percent identity cutoff (default: 40)', default=40, type=int)
parser.add_argument('-c','--qcovhsp', help='BLAST hsp query coverage cutoff (default: 80)', default=80, type=int)
parser.add_argument('-t','--threads', help='Number of threads to use (default: 1)', default=1, type=int)
parser.add_argument('--breakpoint', action='store_true', help='Calculate breakpoint distance statistics (default: do not calculate unless --besthits file is provided)')
args = parser.parse_args()
outputpath=os.path.relpath(args.out, cwdir)


if args.sequences==None and args.besthits==None:
    parser.error('as input, you must either provide --sequences or --besthits')

if args.sequences!=None:
    runsubprocess(['python','%s/getproteins.py'%sourcedir,outputpath, str(args.sequences)])  
    runsubprocess(['bash','%s/makeblastdbs.sh'%sourcedir,outputpath, str(args.threads), sourcedir])
    runsubprocess(['bash','%s/runblast.sh'%sourcedir,outputpath, str(args.evalue),str(args.threads)])
    runsubprocess(['bash','%s/reformatblast.sh'%sourcedir,outputpath,str(args.pident),str(args.qcovhsp)])
    runsubprocess(['Rscript','%s/getreciprocalhits.R'%sourcedir,outputpath])

    if args.breakpoint==True:
        rbhinput='metamorth'
        runsubprocess(['Rscript','%s/getbreakpointdistance.R'%sourcedir,outputpath,str(args.threads),rbhinput])

else:
    rbhinput='userprovided'
    runsubprocess(['mkdir -p %s/blast'%outputpath],shell=True)
    runsubprocess(['mkdir -p %s/output'%outputpath],shell=True)
    runsubprocess(['Rscript','%s/getbreakpointdistance.R'%sourcedir,outputpath,str(args.threads),rbhinput,str(args.besthits)])

