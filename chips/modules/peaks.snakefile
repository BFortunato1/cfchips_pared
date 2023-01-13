#MODULE: PEAK CALLING using macs2

#PARAMETERS
_logfile="analysis/logs/peaks.log"
_macs_fdr="0.01"
_macs_keepdup="1"
_macs_extsize="146"
_macs_species="hs"

def getTreats(wildcards):
    tmp=[]
    #print(wildcards)
    rep_s = wildcards.rep
    rep_n = int(rep_s[-1])
    #print(rep_n)

    #USE the formula: treatment = 2^(rep-1); control = treatment+1
    treat = 2**(rep_n-1) if rep_n > 1 else 0
    r = config['runs'][wildcards.run]
    #print(r)
    #print(treat)
    if treat < len(r) and r[treat]:
        tmp = ["analysis/align/%s/%s_unique.sorted.bam" % (r[treat],r[treat])]
    #print("TREAT: %s" % tmp)
    if not tmp:
        #NOTE: I can't figure out a proper kill command so I'll do this
        tmp=["ERROR! Missing treatment file for run: %s, rep: %s" % (wildcards.run, rep_s)]
    return tmp

def getConts(wildcards):
    tmp=[]
    #print(wildcards)
    rep_s = wildcards.rep
    rep_n = int(rep_s[-1])

    #USE the formula: treatment = 2^(rep-1); control = treatment+1
    cont = 2**(rep_n-1) + 1 if rep_n > 1 else 1
    r = config['runs'][wildcards.run]
    #print(r)
    #print(cont)
    if cont < len(r) and r[cont]:
        tmp = ["analysis/align/%s/%s_unique.sorted.bam" % (r[cont],r[cont])]
    #print("CONT: %s" % tmp)
    return tmp

def checkBAMPE(wildcards):
    """Fn returns '-f BAMPE' IF the run's FIRST treatment replicate (sample) is
    Paired-END.
    NOTE: this flag is required for macs2 callpeak, AUTO detect format does not
    work with PE bams!
    """
    r = config['runs'][wildcards.run]
    #GET first treatement sample
    first = config['samples'][r[0]]
    ret = "-f BAMPE" if len(first) == 2 else ""
    return ret

#NOTE: using the _refs from chips.snakefile
def peaks_targets(wildcards):
    """Generates the targets for this module"""
    #print(wildcards)
    ls = []
    for run in config["runs"].keys():
        #NOTE: using the fact that _reps has this info parsed already!
        for rep in _reps[run]:
            #GENERATE Run name: concat the run and rep name
            runRep = "%s.%s" % (run, rep)
            ls.append("analysis/peaks/%s/%s_sorted_peaks.narrowPeak" % (runRep,runRep))
            ls.append("analysis/peaks/%s/%s_sorted_peaks.narrowPeak.bed" % (runRep,runRep))
            ls.append("analysis/peaks/%s/%s_sorted_summits.bed" % (runRep,runRep))
            ls.append("analysis/peaks/%s/%s_treat_pileup.bw" % (runRep,runRep))
#            ls.append("analysis/peaks/%s/%s_control_lambda.bw" % (runRep,runRep))
#            ls.append("analysis/peaks/%s/%s_treat_pileup.sorted.bdg.gz" % (runRep,runRep))
#            ls.append("analysis/peaks/%s/%s_control_lambda.sorted.bdg.gz" % (runRep,runRep))

    ls.append("analysis/peaks/peakStats.csv")
    ls.append("analysis/peaks/run_info.txt")
    ls.append("analysis/peaks/all_treatments.igv.xml")
    return ls

rule peaks_all:
    input:
        peaks_targets

rule macs2_callpeaks:
    input:
        treat=getTreats,
        cont=getConts
    output:
        "analysis/peaks/{run}.{rep}/{run}.{rep}_peaks.narrowPeak",
        "analysis/peaks/{run}.{rep}/{run}.{rep}_summits.bed",
        temp("analysis/peaks/{run}.{rep}/{run}.{rep}_treat_pileup.bdg"),
        temp("analysis/peaks/{run}.{rep}/{run}.{rep}_control_lambda.bdg"),
    params:
        fdr=_macs_fdr,
        keepdup=_macs_keepdup,
        extsize=_macs_extsize,
        species=_macs_species,
        outdir="analysis/peaks/{run}.{rep}/",
        name="{run}.{rep}",
        #handle PE alignments--need to add -f BAMPE to macs2 callpeaks
        BAMPE = lambda wildcards: checkBAMPE(wildcards),
        pypath="PYTHONPATH=%s" % config["python2_pythonpath"],
    message: "PEAKS: calling peaks with macs2"
    log:_logfile
    run:
        treatment = "-t %s" % input.treat if input.treat else "",
        control = "-c %s" % input.cont if input.cont else "",        
        shell("{params.pypath} {config[macs2_path]} callpeak --SPMR -B -q {params.fdr} --keep-dup {params.keepdup} -g {params.species} {params.BAMPE} --extsize {params.extsize} --nomodel {treatment} {control} --outdir {params.outdir} -n {params.name} 2>>{log}")

rule peakToBed:
    """Convert MACS's narrowPeak format, which is BED12 to BED5"""
    input:
        "analysis/peaks/{run}.{rep}/{run}.{rep}_sorted_peaks.narrowPeak"
    output:
        "analysis/peaks/{run}.{rep}/{run}.{rep}_sorted_peaks.narrowPeak.bed"
    message: "PEAKS: Converting peak file to bed file"
    log:_logfile
    shell:
        "cut -f1,2,3,4,9 {input} > {output} 2>>{log}"

rule sortSummits:
    input:
        "analysis/peaks/{run}.{rep}/{run}.{rep}_summits.bed"
    output:
        "analysis/peaks/{run}.{rep}/{run}.{rep}_sorted_summits.bed"
    message: "PEAKS: sorting the summits bed by score"
    log:_logfile
    shell:
        "sort -r -n -k 5 {input} > {output} 2>>{log}"

rule sortNarrowPeaks:
    input:
        "analysis/peaks/{run}.{rep}/{run}.{rep}_peaks.narrowPeak"
    output:
        "analysis/peaks/{run}.{rep}/{run}.{rep}_sorted_peaks.narrowPeak"
    message: "PEAKS: sorting the narrowPeaks by -log10qval (col9)"
    log:_logfile
    shell:
        "sort -r -n -k 9 {input} > {output} 2>>{log}"

rule sortBedgraphs:
    """Sort bed graphs--typically useful for converting bdg to bw"""
    input:
        "analysis/peaks/{run}.{rep}/{run}.{rep}_{suffix}.bdg"
    output:
        temp("analysis/peaks/{run}.{rep}/{run}.{rep}_{suffix}.sorted.bdg")
    params:
        #msg just for message below
        msg= lambda wildcards: "%s.%s_%s" % (wildcards.run, wildcards.rep, wildcards.suffix)
    message: "PEAKS: sorting bdg pileups {params.msg}"
    log:_logfile
    shell:
        "bedSort {input} {output} 2>>{log}"

rule bdgToBw:
    """Convert bedGraphs to BigWig"""
    input:
        "analysis/peaks/{run}.{rep}/{run}.{rep}_{suffix}.sorted.bdg"
    output:
        "analysis/peaks/{run}.{rep}/{run}.{rep}_{suffix}.bw"
    params:
        chroms=config['chrom_lens'],
        #msg just for message below
        msg= lambda wildcards: "%s.%s_%s" % (wildcards.run, wildcards.rep, wildcards.suffix)
    message: "PEAKS: Convert bedGraphs to BigWig {params.msg}"
    log:_logfile
    shell:
        "bedGraphToBigWig {input} {params.chroms} {output} 2>>{log}"

rule gzip_bdg:
    """Space saving rule to compress the bdg output"""
    input:
        bdg="analysis/peaks/{run}.{rep}/{run}.{rep}_{suffix}.sorted.bdg",
        #NOTE: the .bw is NOT used, but it helps ensure rule bdgToBw runs first
        bw="analysis/peaks/{run}.{rep}/{run}.{rep}_{suffix}.bw"
    output:
        temp("analysis/peaks/{run}.{rep}/{run}.{rep}_{suffix}.sorted.bdg.gz")
    params:
        #msg just for message below
        msg= lambda wildcards: "%s.%s" % (wildcards.run, wildcards.rep)
    message: "PEAKS: compressing sorted.bdg {params.msg}"
    log:_logfile
    shell:
        "gzip {input.bdg} 2>> {log}"

rule getPeaksStats:
    """Counts  number of peaks, # of 10FC, # of 20FC peaks for each sample"""
    input:
        #Generalized INPUT fn defined in chips.snakefile
        _getRepInput("analysis/peaks/$runRep/$runRep_sorted_peaks.narrowPeak")
    output:
        "analysis/peaks/peakStats.csv"
    message: "PEAKS: collecting peaks stats for each run"
    log:_logfile
    run:
        files = " -f ".join(input)
        shell("chips/modules/scripts/peaks_getPeakStats.py -f {files} -o {output} 2>>{log}")

rule macsRunInfo:
    """Dump the current version of macs and the fdr used into a text file 
    for the report"""
    params:
        fdr = _macs_fdr,
        pypath="PYTHONPATH=%s" % config["python2_pythonpath"],
    output:
        #MAKE temp
        "analysis/peaks/run_info.txt"
    message: "PEAKS/REPORT - collection macs version and fdr info"
    shell:
        "{params.pypath} {config[macs2_path]} --version 2> {output} && echo fdr {params.fdr} >> {output}"
    
rule generate_IGV_session:
    """Generates analysis/peaks/all_treatments.igv.xml, a igv session of all
    of the treatment.bw files"""
    input:
        _getRepInput("analysis/peaks/$runRep/$runRep_treat_pileup.bw")
    params:
        genome=config['assembly']
    log:_logfile
    output:
        "analysis/peaks/all_treatments.igv.xml"
    message: "PEAKS: generate IGV session for all treatment.bw files"
    run:
        treats = " -t ".join(input)
        shell("chips/modules/scripts/peaks_generateIGVSession.py -g {params.genome} -t {treats} -o {output} 2>>{log}")
        
