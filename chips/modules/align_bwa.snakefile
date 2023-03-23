#MODULE: Align fastq files to genome - BWA specific calls
#PARAMETERS:
_logfile="analysis/logs/align.log"
_bwa_threads=8

#Obtain fastq files
def getFastq(wildcards):
    return config["samples"][wildcards.sample]

#takes fastq files, aligns fragments against a reference genome ('bwa_index'). This is to make sense of an otherwise jumble of fragments. 
#This outputs to a sam file and then is piped to an aligned bam file using 'samtools view'.
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


