#REPORT module - must come last in includes!
import sys
from string import Template
from snakemake.utils import report
from snakemake.report import data_uri

#DEPENDENCIES
from tabulate import tabulate

_ReportTemplate = Template(open("/n/scratch3/users/b/brf340/20221105_novogene_partial_Reduced_B_A/chips/static/chips_report.txt").read())
_logfile = "analysis/logs/report.log"

def report_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    ls.append('analysis/report/sequencingStatsSummary.csv')
    ls.append('analysis/report/peaksSummary.csv')
    ls.append('analysis/report/report.html')
    return ls

def csvToSimpleTable(csv_file):
    """function to translate a .csv file into a reStructuredText simple table
    ref: http://docutils.sourceforge.net/docs/ref/rst/restructuredtext.html#simple-tables
    """
    #read in file
    f = open(csv_file)
    stats = [l.strip().split(",") for l in f]
    f.close()

    hdr = stats[0]
    rest = stats[1:]
    #RELY on the tabulate pkg
    ret = tabulate(rest, hdr, tablefmt="rst")
    return ret

def processRunInfo(run_info_file):
    """extracts the macs version and the fdr used first and second line"""
    f = open(run_info_file)
    ver = f.readline().strip()
    fdr = f.readline().strip()
    return (ver,fdr)

def sampleGCandContam_table(fastqc_stats, fastqc_gc_plots, contam_table):
    """generates rst formatted table that will contain the fastqc GC dist. plot
    for all samples
    hdr = array of header/column elms
    NOTE: we use the thumb nail image for the GC content plots
    """
    #READ in fastqc_stats
    f_stats = {}
    f = open(fastqc_stats)
    hdr = f.readline().strip().split(",")  #largely ignored
    for l in f:
        tmp = l.strip().split(",")
        #store GC content, 3 col, using sample names as keys
        f_stats[tmp[0]] = tmp[2]
    f.close()
    
    #READ in contam_panel
    contam = {}
    f = open(contam_table)
    hdr = f.readline().strip().split(',') #We'll use this!
    species = hdr[1:]
    for l in f:
        tmp = l.strip().split(",")
        #build dict, use sample name as key
        contam[tmp[0]] = zip(species, tmp[1:]) #dict of tupes, (species, %)
    f.close()

    #PUT GC plots into a dictionary--keys are sample names
    plots = {}
    for p in fastqc_gc_plots:
        #NOTE the path structure is analysis/fastqc/{sample}/png_filename
        #we take second to last
        tmp = p.split("/")[-2]
        plots[tmp] = p
    
    #build output
    ret=[]
    samples = sorted(list(f_stats.keys()))
    hdr = ["Sample", "GC median", "GC distribution"]
    hdr.extend(species)
    rest=[]

    for sample in samples:
        #HANDLE null values
        if plots[sample] and (plots[sample] != 'NA'):
            gc_plot = ".. image:: %s" % data_uri(plots[sample])
        else:
            gc_plot = "NA"
        #get rest of values and compose row
        contam_values = [v for (s, v) in contam[sample]]
        row = [sample, f_stats[sample], gc_plot]
        row.extend(contam_values)

        rest.append(row)
    ret = tabulate(rest, hdr, tablefmt="rst")
    return ret

rule report_all:
    input:
        report_targets

rule report:
    input:
        cfce_logo="/n/scratch3/users/b/brf340/20221105_novogene_partial_Reduced_B_A/chips/static/CFCE_Logo_Final.jpg",
        run_info="analysis/peaks/run_info.txt",
        map_stat="analysis/report/mapping.png",
        pbc_stat="analysis/report/pbc.png",
        peakFoldChange_png="analysis/report/peakFoldChange.png",
        conservPlots = sorted(_getRepInput("analysis/conserv/$runRep/$runRep_conserv_thumb.png")),
	samples_summary="analysis/report/sequencingStatsSummary.csv",
	runs_summary="analysis/report/peaksSummary.csv",
        fastqc_stats="analysis/fastqc/fastqc.csv",
        fastqc_gc_plots = expand("analysis/fastqc/{sample}/{sample}_perSeqGC_thumb.png", sample=config["samples"]),
    params:
        #OBSOLETE, but keeping around
        samples = config['samples']
    output: html="analysis/report/report.html"
    run:
        (macsVersion, fdr) = processRunInfo(input.run_info)
        samplesSummaryTable = csvToSimpleTable(input.samples_summary)
        runsSummaryTable = csvToSimpleTable(input.runs_summary)
        tmp = _ReportTemplate.substitute(cfce_logo=data_uri(input.cfce_logo),map_stat=data_uri(input.map_stat),pbc_stat=data_uri(input.pbc_stat),peakSummitsTable=peakSummitsTable,peakFoldChange_png=data_uri(input.peakFoldChange_png))
        report(tmp, output.html, metadata="Len Taing", **input)


rule samples_summary_table:
    input:
        fastqc = "analysis/fastqc/fastqc.csv",
        mapping = "analysis/align/mapping.csv",
        pbc = "analysis/frips/pbc.csv",
    output:
        "analysis/report/sequencingStatsSummary.csv"
    log: _logfile
    shell:
        "chips/modules/scripts/get_sampleSummary.py -f {input.fastqc} -m {input.mapping} -p {input.pbc} > {output} 2>>{log}"

rule runs_summary_table:
    input:
        peaks = "analysis/peaks/peakStats.csv",
        frips = "analysis/frips/frips.csv",
        dhs = "analysis/ceas/dhs.csv",
        meta = "analysis/ceas/meta.csv",
    output:
        "analysis/report/peaksSummary.csv"
    log: _logfile
    shell:
        "chips/modules/scripts/get_runsSummary.py -p {input.peaks} -f {input.frips} -d {input.dhs} -m {input.meta} -o {output} 2>>{log}"

######## PLOTS ######
rule plot_map_stat:
    input:
        "analysis/align/mapping.csv"
    output:
        "analysis/report/mapping.png"
    log: _logfile
    shell:
        "Rscript chips/modules/scripts/map_stats.R {input} {output}"

rule plot_pbc_stat:
    input:
        #"analysis/align/pbc.csv"
        "analysis/frips/pbc.csv"
    output:
        "analysis/report/pbc.png"
    log: _logfile
    shell:
        "Rscript chips/modules/scripts/plot_pbc.R {input} {output}"

rule plot_peakFoldChange:
    input: 
        "analysis/peaks/peakStats.csv"
    output:
        "analysis/report/peakFoldChange.png"
    log: _logfile
    shell:
        "Rscript chips/modules/scripts/plot_foldChange.R {input} {output}"
