#MODULE: ichorCNA - estimate CNAs and tumor fraction 

#configfile: "config/config.ichorCNA.yaml"
#configfile: "chips/modules/config/config.ichorCNA.yaml"

def ichorCNA_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for sample in config["samples"]:
        ls.append("analysis/ichorCNA/%s/%s.cna.seg" % (sample,sample))
        ls.append("analysis/ichorCNA/readDepth/%s/%s.wig" % (sample,sample))
    return(ls)

rule ichorCNA_all:
    input: 
        ichorCNA_targets

rule target_windows:
    input:
        "analysis/align/{samples}/{samples}_unique.sorted.dedup.bam"
    output:
        counts = "analysis/ichorCNA/targets/{samples}/{samples}.counts.bed",
        filt = "analysis/ichorCNA/targets/{samples}/{samples}.filt.bed"
    params:
        win = config["target_windows"] 
    shell:
        """ 
        samtools view -H {input} | grep ":chr" | sed 's/.*SN://' | sed 's/LN://' > {output.counts}.order.tmp
        grep -v chrM {params.win} | bedtools sort -i - -g {output.counts}.order.tmp > {output.counts}.tmp.sorted.bed
        bedtools intersect -sorted -g {output.counts}.order.tmp -c -a {output.counts}.tmp.sorted.bed -b {input} > {output.counts}
        nr=`wc -l {output.counts} | cut -f1 -d ' ' | awk '{{print $1 * 0.9}}'`
        sort -k4,4n {output.counts} > {output.counts}.tmp #for some reason, this fails if piped with the next line
	head -n $nr {output.counts}.tmp | sort -k1,1 -k2,2n | bedtools merge -i - > {output.filt} 
        rm {output.counts}.tmp
        """

rule subset_bam:
    input:
        bam = "analysis/align/{samples}/{samples}_unique.sorted.dedup.bam",
        filt = "analysis/ichorCNA/targets/{samples}/{samples}.filt.bed",
        peaks = "analysis/peaks/{samples}.rep1/{samples}.rep1_sorted_peaks.narrowPeak.bed"
    output:
        bam = "analysis/ichorCNA/targets/{samples}/{samples}.filt.bam",
        finalfilt = "analysis/ichorCNA/targets/{samples}/{samples}.final.filter.bed"
    params:
        exclude = config["exclude"],
        chrsize = config["chrsize"]
    shell:
        """
        cat {input.peaks} {params.exclude} | cut -f1-3 | sort -k1,1 -k2,2n | bedtools slop -b 5000 -g {params.chrsize} -i - | bedtools merge -i - > {output.finalfilt}.tmp
        grep -v MT {input.filt} | bedtools complement -i - -g {params.chrsize} | bedtools slop -b 5000 -g {params.chrsize} | bedtools complement -i - -g {params.chrsize} | bedtools subtract -a - -b {output.finalfilt}.tmp > {output.finalfilt}
        samtools view -b -L {output.finalfilt} {input.bam} > {output.bam} 
        samtools index {output.bam}
        """

rule read_counter:
    input:
        "analysis/ichorCNA/targets/{samples}/{samples}.filt.bam"
    output:
        "analysis/ichorCNA/readDepth/{samples}/{samples}.wig"  
    params:
        binSize=config["binSize"],
        qual="20",
        chrs=config["chrs"]
    resources:
        mem=4
    log:
        "logs/readDepth/{samples}.log"
    shell:
        "readCounter {input} -c {params.chrs} -w {params.binSize} -q {params.qual} > {output} 2> {log}"

rule ichorCNA:
    input:
        tum="analysis/ichorCNA/readDepth/{tumor}/{tumor}.wig",
        exons_bed="analysis/ichorCNA/targets/{tumor}/{tumor}.final.filter.bed",
        #norm=lambda wildcards: "analysis/ichorCNA/readDepth/" + config["pairings"][wildcards.tumor] + ".bin" + str(config["binSize"]) + ".wig"
    output:
        #corrDepth="results/ichorCNA/{tumor}/{tumor}.correctedDepth.txt",
        #param="results/ichorCNA/{tumor}/{tumor}.params.txt",
        cna="analysis/ichorCNA/{tumor}/{tumor}.cna.seg",
        #segTxt="results/ichorCNA/{tumor}/{tumor}.seg.txt",
        #seg="results/ichorCNA/{tumor}/{tumor}.seg",
        #rdata="results/ichorCNA/{tumor}/{tumor}.RData"
    params:
        outDir="analysis/ichorCNA/{tumor}/",
        rscript=config["ichorCNA_rscript"],
        id="{tumor}",
        ploidy=config["ichorCNA_ploidy"],
        normal=config["ichorCNA_normal"],
        gcwig=config["ichorCNA_gcWig"],
        mapwig=config["ichorCNA_mapWig"],
#        normalpanel=config["ichorCNA_normalPanel"],
        estimateNormal=config["ichorCNA_estimateNormal"],
        estimatePloidy=config["ichorCNA_estimatePloidy"],
        estimateClonality=config["ichorCNA_estimateClonality"],
        scStates=config["ichorCNA_scStates"],
        maxCN=config["ichorCNA_maxCN"],
        includeHOMD=config["ichorCNA_includeHOMD"],
        chrs=config["ichorCNA_chrs"],
        chrTrain=config["ichorCNA_chrTrain"],
        genomeBuild=config["ichorCNA_genomeBuild"],
        genomeStyle=config["ichorCNA_genomeStyle"],
        centromere=config["ichorCNA_centromere"],
        fracReadsChrYMale=config["ichorCNA_fracReadsInChrYForMale"],
        minMapScore=config["ichorCNA_minMapScore"],
        maxFracGenomeSubclone=config["ichorCNA_maxFracGenomeSubclone"],
        maxFracCNASubclone=config["ichorCNA_maxFracCNASubclone"],
#        exons=config["ichorCNA_exons"],
        txnE=config["ichorCNA_txnE"],
        txnStrength=config["ichorCNA_txnStrength"],
        plotFileType=config["ichorCNA_plotFileType"],
        plotYlim=config["ichorCNA_plotYlim"],
        libdir=config["ichorCNA_libdir"]
    resources:
        mem=4
    log:
        "logs/ichorCNA/{tumor}.log"    
#    conda:
#        config['ichorCNA_env']
    shell:
        "Rscript {params.rscript} --id {params.id} --libdir {params.libdir} --WIG {input.tum} --gcWig {params.gcwig} --mapWig {params.mapwig} --ploidy \"{params.ploidy}\" --normal \"{params.normal}\" --maxCN {params.maxCN} --includeHOMD {params.includeHOMD} --chrs \"{params.chrs}\" --chrTrain \"{params.chrTrain}\" --genomeStyle {params.genomeStyle} --genomeBuild {params.genomeBuild} --estimateNormal {params.estimateNormal} --estimatePloidy {params.estimatePloidy} --estimateScPrevalence {params.estimateClonality} --scStates \"{params.scStates}\" --centromere {params.centromere} --exons.bed {input.exons_bed} --txnE {params.txnE} --txnStrength {params.txnStrength} --minMapScore {params.minMapScore} --fracReadsInChrYForMale {params.fracReadsChrYMale} --maxFracGenomeSubclone {params.maxFracGenomeSubclone} --maxFracCNASubclone {params.maxFracCNASubclone} --plotFileType {params.plotFileType} --plotYLim \"{params.plotYlim}\" --outDir {params.outDir} > {log} 2> {log}"
#        "Rscript {params.rscript} --id {params.id} --libdir {params.libdir} --WIG {input.tum} --gcWig {params.gcwig} --mapWig {params.mapwig} --normalPanel {params.normalpanel} --ploidy \"{params.ploidy}\" --normal \"{params.normal}\" --maxCN {params.maxCN} --includeHOMD {params.includeHOMD} --chrs \"{params.chrs}\" --chrTrain \"{params.chrTrain}\" --genomeStyle {params.genomeStyle} --genomeBuild {params.genomeBuild} --estimateNormal {params.estimateNormal} --estimatePloidy {params.estimatePloidy} --estimateScPrevalence {params.estimateClonality} --scStates \"{params.scStates}\" --centromere {params.centromere} --exons.bed {params.exons} --txnE {params.txnE} --txnStrength {params.txnStrength} --minMapScore {params.minMapScore} --fracReadsInChrYForMale {params.fracReadsChrYMale} --maxFracGenomeSubclone {params.maxFracGenomeSubclone} --maxFracCNASubclone {params.maxFracCNASubclone} --plotFileType {params.plotFileType} --plotYLim \"{params.plotYlim}\" --outDir {params.outDir} > {log} 2> {log}"


