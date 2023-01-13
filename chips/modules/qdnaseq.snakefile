#MODULE: qdnaseq- perform qdnaseq CNV analysis on sample bams
_logfile="analysis/logs/qdnaseq.log"

def qdnaseq_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    ls.append("analysis/qdnaseq/qdnaseq.bed")
    ls.append("analysis/qdnaseq/qdnaseq.igv")
    ls.append("analysis/qdnaseq/qdnaseq.txt")
    ls.append("analysis/qdnaseq/qdnaseq.pdf")
    ls.append("analysis/qdnaseq/qdnaseq_segmented.igv")
    ls.append("analysis/qdnaseq/qdnaseq_calls.igv")
    ls.append("analysis/qdnaseq/qdnaseq_genes.txt")
    ls.append("analysis/qdnaseq/qdnaseq_genes.igv")
    return ls

#NOTE: can not initialize _run_controls to [], otherwise getControls won't run
_run_controls = None
def getControls(wildcards):
    """Tries to grab only the control files for each run.
    NOTE: one control can be used for multiple runs--we only use it once"""
    tmp = []
    for (runs, ls) in config['runs'].items():
        if ls[1]:
            tmp.append(ls[1])
        if ls[3]:
            tmp.append(ls[3])

    #get only unique by making it a set
    _run_controls = set(tmp)
    #print(_run_controls)

    #this is wrong...try the other
    #ret = ["analysis/align/%s/%s_unique.sorted.bam" % (s, s) for s in _run_controls]
    ret = [".tmp/%s_unique.sorted.bam" % s for s in _run_controls]
    #print(ret)
    return ret

rule qdnaseq_all:
    input:
        qdnaseq_targets

rule qdnaseq_linkFiles:
    """Link files into one single (temporary) directory"""
    input:
        "analysis/align/{sample}/{sample}_unique.sorted.bam"
    output:
        temp(".tmp/{sample}_unique.sorted.bam")
    message: "QDNAseq: linking files"
    log: _logfile
    shell:
        "ln -s ../{input} {output} 2>>{log}"

rule qdnaseq:
    """performs qdnaseq analysis on ALL of the samples (made by linkFiles)"""
    input:
        getControls
    output:
        "analysis/qdnaseq/qdnaseq.bed",
        "analysis/qdnaseq/qdnaseq.igv",
        "analysis/qdnaseq/qdnaseq.txt",
        "analysis/qdnaseq/qdnaseq.pdf",
        "analysis/qdnaseq/qdnaseq_segmented.igv",
        "analysis/qdnaseq/qdnaseq_calls.igv",
    message: "QDNAseq: performing cnv analysis"
    log: _logfile
    params:
        name="qdnaseq",
        qbin="chips/static/qdnaseq/qdnaseq_hg19_50.bin",
        out="analysis/qdnaseq/"
    shell:
        "R CMD BATCH --vanilla '--args {params.name} .tmp {params.qbin} {params.out}' chips/modules/scripts/qdnaseq.R {log}"

rule qdnaseq_annotate:
    """Processes the segmented.igv file, which is region-based, and creates
    genes.igv file which is gene-based."""
    input:
        "analysis/qdnaseq/qdnaseq_segmented.igv",
    output:
        "analysis/qdnaseq/qdnaseq_genes.txt",
        "analysis/qdnaseq/qdnaseq_genes.igv",
    message: "QDNASEQ: annotating cnv analysis output"
    log: _logfile
    params:
        geneTable=config['geneTable'],
        outPath="analysis/qdnaseq/"
    shell:
        "chips/modules/scripts/qdnaseq_annotate.py -g {params.geneTable} -i {input} -o {params.outPath}"
