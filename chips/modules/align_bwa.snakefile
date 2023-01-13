#MODULE: Align fastq files to genome - BWA specific calls
#PARAMETERS:
_logfile="analysis/logs/align.log"
_bwa_threads=8

def getFastq(wildcards):
    return config["samples"][wildcards.sample]

rule bwa_mem:
    input:
        getFastq
    output:
        temp("analysis/align/{sample}/{sample}.bam")
    params:
        index=config['bwa_index'],
    threads: _bwa_threads
    message: "ALIGN: Running BWA mem for alignment"
    log: _logfile
    shell:
        "bwa mem -t {threads} {params.index} {input} | samtools view -Sb - > {output} 2>>{log}"


