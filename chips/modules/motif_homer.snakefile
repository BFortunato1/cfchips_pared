#MODULE: homer motif module--performs motif analysis and generates motif table
import subprocess
_logfile="analysis/logs/motif.log"
_threads=8
_minPeaks = 500

#NOTE: using the _refs from chips.snakefile
def motif_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config["runs"].keys():
        for rep in _reps[run]:
            runRep = "%s.%s" % (run, rep)
            ls.append("analysis/motif/%s/results/homerResults.html" % runRep)
            ls.append("analysis/peaks/%s/%s_annotatePeaks.txt" % (runRep,runRep))
            ls.append("analysis/peaks/%s/%s_annotatePeaks.tsv" % (runRep,runRep))
            ls.append("analysis/peaks/%s/%s_annotatePeaks.csv" % (runRep,runRep))
    #ls.append("analysis/motif/motifSummary.csv")
    return ls

rule motif_all:
    input:
        motif_targets

def _createEmptyMotif(motif_html):
    """When the _sorted_5k_summits.bed has too few peaks, or is empty,
    we still want to create an emtpy homerResult.html
    INPUT: output paths of these files
    """
    #CHECK for dir existence:
    _path = "/".join(motif_html.split("/")[:-1])
    if not os.path.exists(_path):
        os.makedirs(_path, exist_ok=True)
    #Create an empty mdseqpos_index.html
    subprocess.call(['touch', motif_html])

rule motif_homer:
    """call HOMER on top 5k summits"""
    input:
        bed = "analysis/peaks/{run}.{rep}/{run}.{rep}_sorted_5k_summits.bed"
    output:
        results="analysis/motif/{run}.{rep}/results",
        html="analysis/motif/{run}.{rep}/results/homerResults.html",
    params:
        genome=config['motif_path'],
        size=600,
    message: "MOTIF: calling HOMER on top 5k summits"
    threads:_threads
    log: _logfile
    run:
        #check to see if _sorted_5k_summits.bed is valid
        wc = str(subprocess.check_output(['wc', '-l', input.bed]))
        #this returns a byte string--we need to eliminate the b' ... '
        #and then convert the first elm to int
        wc = int(wc[2:-1].split()[0])

        if wc >= _minPeaks:
            #PASS- run motif scan
            shell("findMotifsGenome.pl {input} {params.genome} {output.results} -size {params.size} -p {threads} -mask -seqlogo -preparsedDir {output.results} >>{log} 2>&1")
        else:
            #FAIL - create empty outputs
            _createEmptyMotif(output.html)

rule getMotifSummary:
    """Summarize the top hits for each run into a file"""
    input:
        #Generalized INPUT fn defined in chips.snakefile
        _getRepInput("analysis/motif/$runRep/results/homerResults.html")
    output:
        "analysis/motif/motifSummary.csv"
    message: "MOTIF: summarizing motif runs"
    log: _logfile
    run:
        files = " -m ".join(input)
        shell("chips/modules/scripts/motif_homerSummary.py -m {files} -o {output} 2>> {log}")

rule homer_annotatePeaks:
    """Annotate peak files.
    NOTE: only for motif_homer modules
    """
    input:
        bed = "analysis/peaks/{run}.{rep}/{run}.{rep}_sorted_peaks.narrowPeak.bed"
    output:
        "analysis/peaks/{run}.{rep}/{run}.{rep}_annotatePeaks.txt"
    params:
        genome=config['motif_path'],
    message: "MOTIF: homer annotatePeaks"
    log: _logfile
    shell:
        "annotatePeaks.pl {input} {params.genome} > {output}"

rule homer_processAnnPeaks:
    """Process peaks/{run}/{run}_annotatePeaks.txt files"""
    input:
        "analysis/peaks/{run}.{rep}/{run}.{rep}_annotatePeaks.txt"
    output:
        tsv="analysis/peaks/{run}.{rep}/{run}.{rep}_annotatePeaks.tsv",
        csv="analysis/peaks/{run}.{rep}/{run}.{rep}_annotatePeaks.csv",
    message: "MOTIF: Post-process homer annotatePeaks.txt file"
    log: _logfile
    shell:
        "chips/modules/scripts/motif_annPeaksTsvCsv.sh {input} {output.tsv} {output.csv}"
