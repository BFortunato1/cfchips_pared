#MODULE: Align fastq files to genome - BOWTIE2 specific calls
#PARAMETERS:
_logfile="analysis/logs/align.log"
_bwt2_threads=8
_samtools_threads=4

#grabs fastq files
def getFastq(wildcards):
    return config["samples"][wildcards.sample]

#using bowtie2, fastq files are aligned against a reference index into sam files
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

#sam files are converted to bam files
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

