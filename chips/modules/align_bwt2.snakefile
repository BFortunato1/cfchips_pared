#MODULE: Align fastq files to genome - BOWTIE2 specific calls
#PARAMETERS:
_logfile="analysis/logs/align.log"
#_bwa_q="5"
#_bwa_l="32"
#_bwa_k="2"
_bwt2_threads=8
_samtools_threads=4

def getFastq(wildcards):
    return config["samples"][wildcards.sample]

#def getFastq2(wildcards):
#    return config["samples"][wildcards.sample]

# def getMates(wildcards):
#     s = wildcards.sample
#     files = config["samples"][s]
#     return ["analysis/align/%s/%s_%s.sai" % (s,s,m) for m in range(len(files))]

# def getRunType(wildcards):
#     s = wildcards.sample
#     return "sampe" if len(config["samples"][s]) > 1 else "samse"

rule bwt2_aln:
    input:
        getFastq
    output:
        temp("analysis/align/{sample}/{sample}.sam")
    params:
        index=config['bwt2_index'],
    threads: _bwt2_threads
    message: "ALIGN: Running Bowtie2 alignment"
    log: _logfile
    run:
        if len(input) > 1:
            #PE
            _inputs=("-1 %s -2 %s" % (input[0], input[1]))
        else:
            #SE
            _inputs="-U %s" % input
        shell("module load gcc/6.2.0 bowtie2/2.3.4.3 && bowtie2 -p {threads} -x {params.index} --minins 50 --maxins 500 --dovetail {_inputs} -S {output} 2>>{log}")

rule bwt2_convert:
    input:
        sam="analysis/align/{sample}/{sample}.sam"
    output:
        temp("analysis/align/{sample}/{sample}.bam")
    threads: _samtools_threads
    message: "ALIGN: convert sam to bam"
    log: _logfile
    shell:
        "samtools view -@ {threads} -Sb {input} > {output}"

