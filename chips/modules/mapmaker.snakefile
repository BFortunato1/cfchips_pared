#MODULE: mapmaker- generating the analysis/mapmaker folder
#The folder should be a fully runnable version of the latest mapmaker USING
#the chips run information
_logfile="analysis/logs/mapmaker.log"

def mapmaker_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    ls.append("analysis/mapmaker")
    ls.append("analysis/mapmaker/config.yaml")
    ls.append("analysis/mapmaker/metasheet.csv")
    return ls

def make_relative(path):
    """GIVEN a path that is relative to the PROJECT directory, e.g.
    analysis/peaks/{run}.rep1/.., we convert the path to be relative to the 
    analysis/mapmaker directory.  In other words, we replace analysis with ..
    """
    return path.replace("analysis","..")

def getTreat(run):
    """GIVEN a run, will return the sample name of the first treatment"""
    treat1 = config['runs'][run][0]
    return treat1

rule mapmaker_all:
    input:
        mapmaker_targets

rule mapmaker:
    """Clones the mapmaker project"""
    message: "MAPMAKER: cloning the mapmaker project from github"
    log: _logfile
    output:
        "analysis/mapmaker",
        "analysis/mapmaker/.config.yaml",
        "analysis/mapmaker/.metasheet.csv"
    shell:
        "git clone https://bitbucket.org/cfce/cobra.git analysis/mapmaker > {log} 2>&1 && "
        "sleep 10 && "
        #NOTE: the two files below are OTHER targets of this snakefile, so
        #we trick snakemake into thinking they don't exist by renaming them
        "mv analysis/mapmaker/config.yaml analysis/mapmaker/.config.yaml && "
        "mv analysis/mapmaker/metasheet.csv analysis/mapmaker/.metasheet.csv"

rule mapmaker_config:
    """Tries to configure the mapmaker run based on the chips run info"""
    input:
        bed_files= expand("analysis/peaks/{run}.rep1/{run}.rep1_sorted_peaks.narrowPeak.bed", run=config['runs']),
        bw_files= expand("analysis/peaks/{run}.rep1/{run}.rep1_treat_pileup.bw", run=config['runs']),
        bam_files = expand("analysis/align/{sample}/{sample}_unique.sorted.dedup.bam", sample=[getTreat(r) for r in config['runs']]),
        igv_files= expand("analysis/peaks/{run}.rep1/{run}.rep1_treatment.igv.xml", run=config['runs']),
        config="analysis/mapmaker/.config.yaml"
    params:
        run_names=config['runs']
    message: "MAPMAKER: configuring mapmaker"
    log: _logfile
    output:
        "analysis/mapmaker/config.yaml"
    run:
        #TRANSFORM to analysis/mapmaker relative
        #input.bed_files = [make_relative(b) for b in bed_files]
        bed_files = " -b ".join([make_relative(b) for b in input.bed_files])
        bw_files = " -w ".join([make_relative(bw) for bw in input.bw_files])
        bam_files= " -a ".join([make_relative(bam) for bam in input.bam_files])
        igv_files= " -i ".join([make_relative(igv) for igv in input.igv_files])
        names = " -n ".join(params.run_names)

        shell("chips/modules/scripts/mapmaker_config.py -n {names} -b {bed_files} -w {bw_files} -a {bam_files} -i {igv_files} -c {input.config} -o {output}")

rule mapmaker_meta:
    """Tries to configure the mapmaker meta based on the chips run info"""
    input:
        meta="analysis/mapmaker/.metasheet.csv"
    params:
        run_names=config['runs']
    message: "MAPMAKER: configuring mapmaker"
    log: _logfile
    output:
        "analysis/mapmaker/metasheet.csv"
    run:
        names = " -n ".join(params.run_names)
        shell("chips/modules/scripts/mapmaker_meta.py -n {names} -o {output}")
