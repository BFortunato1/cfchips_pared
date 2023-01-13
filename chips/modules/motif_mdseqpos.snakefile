#MODULE: motif module--uses MDSeqPos.py to perform motif analysis and generate 
#motif table
import subprocess
_logfile="analysis/logs/motif.log"

_minPeaks = 500

def motif_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config["runs"].keys():
        #NOTE: using the fact that _reps has this info parsed already!
        for rep in _reps[run]:
            #GENERATE Run name: concat the run and rep name
            runRep = "%s.%s" % (run, rep)
            ls.append("analysis/motif/%s/results/mdseqpos_index.html" % runRep)
            ls.append("analysis/motif/%s/results/motif_list.json" % runRep)
    ls.append("analysis/motif/motifSummary.csv")
    return ls

rule motif_all:
    input:
        motif_targets

def _createEmptyMotif(motif_html, motif_json):
    """When the _sorted_5k_summits.bed has too few peaks, or is empty,
    we still want to create an emtpy mdseqpos_index.html and an 
    appropriate motif_list.json

    INPUT: output paths of these files
    """
    #CHECK for dir existence:
    _path = "/".join(motif_html.split("/")[:-1])
    if not os.path.exists(_path):
        os.makedirs(_path, exist_ok=True)
    #Create an empty mdseqpos_index.html
    subprocess.call(['touch', motif_html])

    #create appropriate json
    #CHECK for dir existence:
    _path = "/".join(motif_json.split("/")[:-1])
    if not os.path.exists(_path):
        os.makedirs(_path, exist_ok=True)
    f = open(motif_json, 'w')
    f.write("{}\n")
    f.close()

rule motif:
    """call MDSeqpos on top 5k summits"""
    input:
        #KEY: Since motif analysis is costly, we're only running it on rep1
        bed = "analysis/peaks/{run}.{rep}/{run}.{rep}_sorted_5k_summits.bed"
    output:
        results="analysis/motif/{run}.{rep}/results",
        html="analysis/motif/{run}.{rep}/results/mdseqpos_index.html",
        json="analysis/motif/{run}.{rep}/results/motif_list.json",
    params:
        genome=config['motif_path'],
        pypath="PYTHONPATH=%s" % config["python2_pythonpath"],
        runName = "{run}.{rep}"
    message: "MOTIF: calling MDSeqPos on top 5k summits"
    log: _logfile
    run:
        #check to see if _sorted_5k_summits.bed is valid
        wc = str(subprocess.check_output(['wc', '-l', input.bed]))
        #this returns a byte string--we need to eliminate the b' ... '
        #and then convert the first elm to int
        wc = int(wc[2:-1].split()[0])

        if wc >= _minPeaks:
            #PASS- run motif scan
            shell("{params.pypath} {config[mdseqpos_path]} {input} {params.genome} -m cistrome.xml -d -O analysis/motif/{params.runName}/results 1>>{log}")
        else:
            #FAIL - create empty outputs
            _createEmptyMotif(output.html, output.json)

rule getMotifSummary:
    """Summarize the top hits for each run into a file"""
    input:
        _getRepInput("analysis/motif/$runRep/results/motif_list.json")
    output:
        "analysis/motif/motifSummary.csv"
    message: "MOTIF: summarizing motif runs"
    log: _logfile
    run:
        files = " -m ".join(input)
        shell("chips/modules/scripts/motif_getSummary.py -m {files} -o {output} 2>> {log}")
