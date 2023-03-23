#!/usr/bin/env python

import os
import sys
import subprocess
import pandas as pd
import yaml

from string import Template

def getRuns(config):
    """parse metasheet for Run groupings"""
    ret = {}

    metadata = pd.read_table(config['metasheet'], index_col=0, sep=',', comment='#', skipinitialspace=True)
    f = metadata.to_csv().split() #make it resemble an actual file with lines
    for l in f[1:]:
        tmp = l.strip().split(",")
        #print(tmp)
        ret[tmp[0]] = tmp[1:]

    #print(ret)
    config['runs'] = ret
    return config

def addPy2Paths_Config(config):
    """ADDS the python2 paths to config"""
    """There are a few dependencies that are incompatible wtih this environment (macs2 - a ChIP-seq peak calling tool, and python2). These are in a separate conda environment (recipe here):
envs/chips_py2_env.yaml"""
    conda_root = subprocess.check_output('conda info --root',shell=True).decode('utf-8').strip()
    conda_path = os.path.join(conda_root, 'pkgs')
    config["python2_pythonpath"] = os.path.join(conda_root, 'envs', 'chips_py2', 'lib', 'python2.7', 'site-packages')
    
    if not "python2" in config or not config["python2"]:
        config["python2"] = os.path.join(conda_root, 'envs', 'chips_py2', 'bin', 'python2.7')

    if not "mdseqpos_path" in config or not config["mdseqpos_path"]:
        config["mdseqpos_path"] = os.path.join(conda_root, 'envs', 'chips_py2', 'bin', 'MDSeqPos.py')

    if not "macs2_path" in config or not config["macs2_path"]:
        config["macs2_path"] = os.path.join(conda_root, 'envs', 'chips_py2', 'bin', 'macs2')

def loadRef(config):
    """Adds the static reference paths found in config['ref']
    NOTE: if the elm is already defined, then we DO NOT clobber the value
    """
    f = open(config['ref'])
    ref_info = yaml.safe_load(f)
    f.close()
    #print(ref_info[config['assembly']])
    for (k,v) in ref_info[config['assembly']].items():
        #NO CLOBBERING what is user-defined!
        if k not in config:
            config[k] = v

#---------  CONFIG set up  ---------------
configfile: "chips/config.yaml"   # This makes snakemake load up yaml into config 
config = getRuns(config)
addPy2Paths_Config(config)

#NOW load ref.yaml - SIDE-EFFECT: loadRef CHANGES config
loadRef(config)
#-----------------------------------------

#------------------------------------------------------------------------------
# Handle replicates
#------------------------------------------------------------------------------
#used to define the replicate structure
_reps = {}
for run in config['runs'].keys():
    r = config['runs'][run]
    tmp = []
    for (rep, i) in enumerate(range(0, len(r), 2)):
        if r[i]: tmp.append("rep%s" % str(rep+1))
    _reps[run] = tmp
#print(_reps)

#NOTE: Template class allows for _ in the variable names, we want to DISALLOW
#that for replicates
#ref: http://stackoverflow.com/questions/2326757/string-templates-in-python-what-are-legal-characters
class RepTemplate(Template):
    idpattern = r'[a-z][a-z0-9]*'

#THIS helper fn is used in several of the modules peaks, ceas, frips
#Instead of an expand, we need this fn to create the CORRECT input-list
def _getRepInput(temp, suffix=""):
    """generalized input fn to get the replicate files
    CALLER passes in temp: a python string template that has the var runRep
    e.g. analysis/ceas/$runRep/$runRep_DHS_stats.txt
    Return: list of the string template filled with the correct runRep names
    """
    #print(temp)
    s = RepTemplate(temp)
    ls = []
    for run in config['runs'].keys():
        for rep in _reps[run]:
            #GENERATE Run name: concat the run and rep name
            runRep = "%s.%s" % (run, rep)
            ls.append(s.substitute(runRep=runRep,))
    #print(ls)
    return ls

#------------------------------------------------------------------------------
# TARGETS
#------------------------------------------------------------------------------

def all_targets(wildcards):
    _qdnaseq = config["cnv_analysis"]
    ls = []
    #IMPORT all of the module targets
    ls.extend(align_targets(wildcards))
    ls.extend(peaks_targets(wildcards))
    ls.extend(fastqc_targets(wildcards))
    ls.extend(conservation_targets(wildcards))
    ls.extend(ceas_targets(wildcards))
    ls.extend(frips_targets(wildcards))

    ls.extend(report_targets(wildcards))
    return ls

rule target:
    input: 
        all_targets,

    message: "Compiling all output"
if config['aligner'] == 'bwt2':
    include: "./modules/align_bwt2.snakefile"     # rules specific to Bowtie2
else:
    include: "./modules/align_bwa.snakefile"      # rules specific to BWA

include: "./modules/align_common.snakefile"  # common align rules
include: "./modules/peaks.snakefile"         # peak calling rules
include: "./modules/fastqc.snakefile"        # fastqc (sequence qual) rules
include: "./modules/conservation.snakefile"  # generate conservation plot
include: "./modules/ceas.snakefile"          # annotate peak regions
include: "./modules/frips.snakefile"         # fraction of reads in peaks

include: "./modules/report.snakefile"        # report module
