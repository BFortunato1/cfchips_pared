#MODULE: contamination module to check for sample contamination
_logfile="analysis/logs/contamination.log"
_bwa_q="5"
_bwa_l="32"
_bwa_k="2"
_bwa_threads=4

###############################################################################
# HELPERS
###############################################################################
_contaminationPanel= config['contamination_panel'] if 'contamination_panel' in config and config['contamination_panel'] else []

def extractIndexName(path):
    """Given a contamination panel path, e.g. /some/path/to/BWA/hg19.fa,
    this fn returns the index filename, e.g. hg19"""
    tmp = path.split("/")[-1]
    #lop off the suffix
    return ".".join(tmp.split(".")[:-1])
_contaminationNames = list(map(extractIndexName, _contaminationPanel))
_contaminationDict = dict(zip(_contaminationNames, _contaminationPanel))
###############################################################################

def contamination_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for sample in config["samples"]:
        for panel in _contaminationNames:
            ls.append("analysis/contam/%s/%s.%s.txt" % (sample, sample, panel))
            ls.append("analysis/contam/%s/%s_contamination.txt" % (sample, sample))
    ls.append("analysis/contam/contamination.csv")
    return ls

def get100k_sample(wildcards):
    """NOTE: there are two different 100k files, either _100k.fastq or
    _100k.bam.fastq.  The former comes from runs that have fastqs as input,
    while latter started w/ bams as inputs.  

    This fn distinguishes which. (based on fastqc.snakefile getFastqcInput
    """
    s = wildcards.sample
    first_file = config["samples"][wildcards.sample][0]
    if first_file.endswith('.bam'):
        ret = ["analysis/align/%s/%s_100k.bam.fastq" % (s,s)]
    else:
        ret = ["analysis/align/%s/%s_100k.fastq" % (s,s)]
    return ret

rule contamination_all:
    input:
        contamination_targets

rule contamination:
    """For each sample, run an alignment for each entry in the 
    contaminationPanel and get the mapping rate"""
    input:
        get100k_sample
    output:
        temp("analysis/contam/{sample}/{sample}.{panel}.bam")
    params:
        index=lambda wildcards: _contaminationDict[wildcards.panel],
        bwa_q = _bwa_q,
        bwa_l = _bwa_l,
        bwa_k = _bwa_k,
        sample = lambda wildcards: wildcards.sample,
        panel = lambda wildcards: wildcards.panel
    threads: _bwa_threads
    message: "CONTAMINATION: checking {params.sample} against {params.panel}"
    log: _logfile
    shell:
        "bwa mem -t {threads} {params.index} {input} | samtools view -Sb - > {output} 2>> {log}"

rule contaminationStats:
    """Extract the mapping stats for each species"""
    input:
        "analysis/contam/{sample}/{sample}.{panel}.bam"
    output:
        #TEMP
        "analysis/contam/{sample}/{sample}.{panel}.txt"
    params:
        sample = lambda wildcards: wildcards.sample,
        panel = lambda wildcards: wildcards.panel,
        #READ out the 5th row, the first element (and divide by 100/100000)
        awk_cmd = "\'BEGIN {RS=\'\\t\'}{print $23 / $1 * 100}\'"
    message: "CONTAMINATION: get mapping stats {params.sample}:{params.panel}"
    log: _logfile
    shell:
        "samtools flagstat {input} | awk {params.awk_cmd} > {output} 2>>{log}"

rule contaminationCollectStats:
    """Collect the mapping stats across the entire panel"""
    input:
        expand("analysis/contam/{{sample}}/{{sample}}.{panel}.txt", panel=_contaminationNames)
    output:
        "analysis/contam/{sample}/{sample}_contamination.txt"
    run:
        for (n, f) in zip(_contaminationNames, input):
            shell("per=$(cat {f}) && echo {n} $per >> {output}")

rule collect_allContamination:
    """Aggregate all of the contamination stats into one large table/panel"""
    input:
        expand("analysis/contam/{sample}/{sample}_contamination.txt", sample=config['samples'])
    message: "Contamination: collecting contamination panel"
    log: _logfile
    output:
        "analysis/contam/contamination.csv"
    run:
        files = " -f ".join(input)
        shell("chips/modules/scripts/contam_getStats.py -f {files} -o {output} 2>>{log}")
