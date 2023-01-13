#MODULE: CEAS- annotating where the peaks fall (promoter, exon, intron, interg)
_logfile="analysis/logs/ceas.log"

#NOTE: using the _refs from chips.snakefile
def ceas_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config["runs"].keys():
        for rep in _reps[run]:
            #GENERATE Run name: concat the run and rep name
            runRep = "%s.%s" % (run, rep)
            ls.append("analysis/ceas/%s/%s_summary.txt" % (runRep,runRep))
            ls.append("analysis/ceas/%s/%s_DHS_peaks.bed" % (runRep,runRep))
            ls.append("analysis/ceas/%s/%s_DHS_stats.txt" % (runRep,runRep))
            ls.append("analysis/ceas/%s/%s_velcro_peaks.bed" % (runRep,runRep))
            ls.append("analysis/ceas/%s/%s_velcro_stats.txt" % (runRep,runRep))

    #ADD bam_regionStats
    for sample in config["samples"]:
        if config['exons']:
            ls.append("analysis/ceas/samples/%s/%s.exons" % (sample,sample))
        if config['promoters']:
            ls.append("analysis/ceas/samples/%s/%s.promoters" % (sample,sample))
        if config['DHS']:
            ls.append("analysis/ceas/samples/%s/%s.DHS" % (sample,sample))
    ls.append("analysis/ceas/samples/bamRegionStats.csv")
    ls.append("analysis/ceas/dhs.csv")
    ls.append("analysis/ceas/meta.csv")
    return ls

rule ceas_all:
    input:
        ceas_targets

rule ceas:
    """Annotate peak regions"""
    input:
        "analysis/peaks/{run}.{rep}/{run}.{rep}_summits.bed"
    output:
        promoter="analysis/ceas/{run}.{rep}/{run}.{rep}_summits_promoter.bed",
        exon="analysis/ceas/{run}.{rep}/{run}.{rep}_summits_exon.bed",
        intron="analysis/ceas/{run}.{rep}/{run}.{rep}_summits_intron.bed",
        intergenic="analysis/ceas/{run}.{rep}/{run}.{rep}_summits_intergenic.bed",
        summary="analysis/ceas/{run}.{rep}/{run}.{rep}_summary.txt",
    message: "CEAS: annotating peak regions"
    log: _logfile
    params:
        db=config['geneTable'],
        path="analysis/ceas/{run}.{rep}/",
        name= "{run}.{rep}_summits",
    shell:
        #TWO ways to run bedAnnotate: w/ basename param (-n) or w/o
        #For now we keep the -n explictly defined
        "chips/modules/scripts/bedAnnotate.v2.py -g {params.db} -b {input} -o {params.path} -n {params.name} > {output.summary} 2>>{log}"
        #"chips/modules/scripts/bedAnnotate.v2.py -g {params.db} -b {input} -o {params.path} > {output.summary} 2>>{log}"


rule takeTop5k:
    """Take the top 5000 sites"""
    input:
        "analysis/peaks/{run}.{rep}/{run}.{rep}_sorted_peaks.narrowPeak.bed"
    params:
        n=5000
    message: "DHS: Take top sites"
    log: _logfile
    output:
        temp('analysis/ceas/{run}.{rep}/{run}.{rep}_sorted_5k_peaks.bed')
    shell:
        "head -n {params.n} {input} > {output} 2>>{log}"

#------------------------------------------------------------------------------
rule DHS_intersectDHS:
    """Intersect PEAKS with DHS regions"""
    input:
        'analysis/ceas/{run}.{rep}/{run}.{rep}_sorted_5k_peaks.bed'
    params:
        dhs=config['DHS']
    message: "DHS: intersect PEAKS with DHS regions"
    log: _logfile
    output:
        'analysis/ceas/{run}.{rep}/{run}.{rep}_DHS_peaks.bed'
    run:
        #HACK! IN the CASE where no DHS is defined/available for species/assemb
        #USE AN EMPTY FILE.  CONSEQUENCE- everything downstream will count 0 
        #DHS regions, NOT N/A
        if config['DHS']:
            shell("intersectBed -wa -u -a {input} -b {params.dhs} > {output} 2>>{log}")
        else:
            #make empty file
            shell("touch .snakemake/null.dhs.txt")
            shell("intersectBed -wa -u -a {input} -b .snakemake/null.dhs.txt > {output} 2>>{log}")

rule DHS_stat:
    """collect DHS stats"""
    input:
        n='analysis/ceas/{run}.{rep}/{run}.{rep}_sorted_5k_peaks.bed',
        dhs='analysis/ceas/{run}.{rep}/{run}.{rep}_DHS_peaks.bed'
    message: "DHS: collecting stats"
    log: _logfile
    output:
        'analysis/ceas/{run}.{rep}/{run}.{rep}_DHS_stats.txt'
    shell:
        "wc -l {input.n} {input.dhs} > {output} 2>>{log}"

#------------------------------------------------------------------------------
rule VELCRO_intersectVelcro:
    """Intersect PEAKS with velcro regions"""
    input:
        'analysis/ceas/{run}.{rep}/{run}.{rep}_sorted_5k_peaks.bed'
    params:
        velcro=config['velcro_regions']
    message: "VELCRO: intersect PEAKS with velcro regions"
    log: _logfile
    output:
        'analysis/ceas/{run}.{rep}/{run}.{rep}_velcro_peaks.bed'
    run:
        #CHECK for the existence of this file!
        if params.velcro:
            shell("intersectBed -wa -u -a {input} -b {params.velcro} > {output} 2>>{log}")
        else:
            #No velcro file defined --> empty output
            shell("touch {output} && echo 'WARNING: no velcro region defined' >>{log}")
    
rule VELCRO_stat:
    """collect VELCRO stats"""
    input:
        n='analysis/ceas/{run}.{rep}/{run}.{rep}_sorted_5k_peaks.bed',
        velcro='analysis/ceas/{run}.{rep}/{run}.{rep}_velcro_peaks.bed'
    message: "VELCRO: collecting stats"
    log: _logfile
    output:
        'analysis/ceas/{run}.{rep}/{run}.{rep}_velcro_stats.txt'
    shell:
        "wc -l {input.n} {input.velcro} > {output} 2>>{log}"

rule bam_regionStat:
    """count the number of reads in promoter, exon, dhs--these regions
    are defined in the config.yaml"""
    input:
        "analysis/align/{sample}/{sample}_4M_unique_nonChrM.bam"
    params:
        bed = lambda wildcards: config[wildcards.region],
        #for use in message only
        msg = lambda wildcards: "%s:%s" % (wildcards.sample, wildcards.region)
    message: "CEAS: bam stat region {params.msg}"
    log: _logfile
    output:
        #make temp
        'analysis/ceas/samples/{sample}/{sample}.{region}'
    shell:
        "chips/modules/scripts/meta_bamRegionCount.sh -i {input} -b {params.bed} -o {output} 2>> {log}"

rule collect_BamRegionStats:
    """collect the BAM region stats into a single file"""
    input:
        #INPUT the stats directories--
        #hack just add all of the file so we dont get a missing input exception
        #and collect the directories down below
        dhs = expand("analysis/ceas/samples/{sample}/{sample}.DHS", sample=config['samples']),
        prom = expand("analysis/ceas/samples/{sample}/{sample}.promoters", sample=config['samples']),
        exon = expand("analysis/ceas/samples/{sample}/{sample}.exons", sample=config['samples'])
    message: "CEAS: collect bam region stats"
    log: _logfile
    output:
        'analysis/ceas/samples/bamRegionStats.csv'
    run:
        #REMOVE the file to get directories - 
        #Aribitrarily USE only the exon list
        dirs = ["/".join(d.split("/")[:-1]) for d in input.exon]
        d = " -d ".join(dirs)
        shell("chips/modules/scripts/ceas_collectBamRegStats.py -d {d} > {output} 2>>{log}")

rule collect_DHSstats:
    """collect the DHS stats into a single file"""
    input:
        #Generalized INPUT fn defined in chips.snakefile
        _getRepInput("analysis/ceas/$runRep/$runRep_DHS_stats.txt")
    message: "CEAS: collect DHS stats"
    log: _logfile
    output:
        'analysis/ceas/dhs.csv'
    run:
        files = " -f ".join(input)
        shell("chips/modules/scripts/peaks_getDHSstats.py -f {files} -o {output} 2>>{log}")

rule collect_CEASstats:
    """collect the CEAS stats into a single file"""
    input:
        #Generalized INPUT fn defined in chips.snakefile
        _getRepInput("analysis/ceas/$runRep/$runRep_summary.txt")
    message: "CEAS: collect CEAS stats"
    log: _logfile
    output:
        'analysis/ceas/meta.csv'
    run:
        files = " -f ".join(input)
        shell("chips/modules/scripts/ceas_getMetaStats.py -f {files} -o {output} 2>>{log}")
