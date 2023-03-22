# cfChips for O2- ChIP-seq analysis pipeline in snakemake specialized for cfCHIP data

Chips is dependent on two conda environments, *chips* and *chips_py2*.
0. **clone the chips source code**:
    `[git clone git@bitbucket.org:cfce/chips](https://github.com/scbaca/chips_o2.git)`

1. **installing chips**:
    After cloning the git repository, create the chips environment by doing this:
    - `cd chips`
    - `conda env create -f environment.yml -n chips`
2. **installing chips_py2**:
    Similarly, create the chips_py2 environment by doing this:
    - `conda env create -f environment_chipsPy2.yml`
3. **Post installation steps: running setupChips.py**
    If you are in the 'chips' directory, setupChips is found in a sub-directory, 'static/scripts'.  To run it:
    1. load the chips conda environment:
    `source activate chips`
    2. `python setupChips.py`
4. **Post installation steps: configuring homer**:
    NOTE: Chips uses the [homer](http://homer.ucsd.edu/homer/motif/index.html) software for motif analysis.  It also has the capability of using the [MDSeqPos](https://bitbucket.org/cistrome/cistrome-applications-harvard/src/c477732c5c88/mdseqpos/) motif finder for a similar analysis.  If you are interested in using MDSeqPos for motif analysis, please see **Appendix D**.
    - To activate/initialize homer:
    1. run the configure script:
    `configureHomer.pl -install`
    2. install the required assemblies:
    For human samples: `configureHomer.pl -install hg19`
    	      	       `configureHomer.pl -install hg38`
    For mouse samples: `configureHomer.pl -install mm9`

### Chips reference files
Chips comes pre-packaged with static reference files (e.g. bwa index, refSeq tables, etc.) for hg19 and mm9.  These files are located at 
/n/scratch3/users/s/scb20/_RESTORE/scb20/cfchip/ref_files
It is suggested that you create a symbolic link to this location.

# Using Chips
### Anatomy of a Chips project
All work in Chips is done in a **PROJECT** directory, which is simply a directory to contain a single Chips analysis run.  **PROJECT** directories can be named anything (and they usually start with a simple mkdir command, e.g. mkdir chips_for_paper),  but what is CRITICAL about a **PROJECT** directory is that you fill them with the following core components:
(We first lay out the directory structure and explain each element below)
> PROJECT/
> chips/
> data/  - *optional*
> config.yaml
> metasheet.csv
> ref.yaml - * **only if you are using chips OTHER THAN hg19 and mm9** *

The 'chips' directory contains all of the chips source code.  We'll explain how to download that directory below.  The 'data' directory is an optional directory that contains all of your raw data.  It is optional because those paths __may__ be fully established in the config.yaml, __however__ it is best practice to gather your raw data within 'data' using [symbolic links](https://www.cyberciti.biz/faq/creating-soft-link-or-symbolic-link/).

The *config.yaml* and *metasheet.csv* are configurations for your VIPER run (also explained below).

The ref.yaml file is explained in **Appendix E**.

After a successful **Chips** run, another 'analysis' folder is generated which contains all of the resulting output files.

### Setting up a Chips project
0. **creating the project directory**
    As explained above, the **PROJECT** directory is simply a directory to contain an entire Chips run.  **It can be named anything, but for this section, we'll simply call it 'PROJECT'**
    `mkdir PROJECT`
    `cd PROJECT`
1. **link data files**
    As explained above, creating a data directory is a place to gather all of your **raw data files (.fastq, .fastq.gz, .bam)**.  It is optional, but **highly recommended**.
    `mkdir data`
    And in 'data', copy over or make symbolic links to your raw data files
2. **clone chips**
    In your PROJECT directory:
    `git clone git@bitbucket.org:cfce/chips`
3. **creating config.yaml and metasheet.csv**
    a. **copy chips/config.yaml and chips/metasheet.csv into the PROJECT dir:**
    In the PROJECT directory:
    `cp chips/config.yaml .`
    `cp chips/metasheet.csv .`

    b. **setup config.yaml**:
        The config.yaml is where you define Chips run parameters and the ChIP-seq samples for analysis.
        
    1. **Set the assembly**: typically hg19 or mm9 (default: hg19)
            
    2. **Choose the aligner**: either bwa or bowtie2 (default: bwa)
    3. **Choose the motif software**: either homer or mdseqpos (default: homer)
    4. **Contamination Panel**:
        The contamination panel is a panel that Chips will check for "cross-species" contamination.  Out of the box, the config.yaml has hg19 and mm9 as assemblies to check.  **IF you would like to add other species/assemblies, simply add as many BWA indices as you would like** (bowtie assemblies will not work, even if your aligner is set to bowtie)
    5. **cnv_analysis**: Set to 'true' to enable copy number variation analysis
    6. **samples**:
        __The most important part of the config file is to define the samples for Chips analysis.__
        Each sample is given an arbitrary name, e.g. MCF7_ER, MCF7_input, etc.  **Sample names, however, can not start with a number, and cannot contain '.', '-' (dash--use underscores instead)** (POSSIBLY others).  For each sample, define the path to the raw data file (.fastq, .fastq.gz, .bam).  For paired-end samples, simply add another line to define the path to the second pair.
    
    c. **setup metasheet.csv**:
    The metasheet.csv is where you group the **samples** (defined in config.yaml) into Treatment, Control (and if applicable, replicates).  For Chips, each of these groupings is called a **run**.
    
    Open metasheet.csv in Excel or in a text-editor.You will notice the first (uncommented) line is:
    `RunName,Treat1,Cont1,Treat2,Cont2`
    
    **RunName**- arbitraty name for the run, e.g. *MCF7_ER_run*
    **Treat1**- The sample name that corresponds to treatment sample.  **It must exactly match the sample name used in config.yaml**
    **Cont1**- (optional) The input-control sample that should be paired with Treat1.
    **Treat2**- (optional) if you have replicates, define the treatment of the replicate here.
    **Cont2**- (optional) the input-control, if available, for Treat2
    
3. setting up refs
    - change config.yaml ref: "chips/ref.yaml"
    - linking to static refs
    - copying ref.yaml
### Running Chips
1. source activate chips
2. dry run (snakemake -s chips/chips.Snakefile -npr)
3. full run (snakemake -s chips/chips.Snakefile -pr --cores 6)

### Appendix A: System requirements
### Appendix B: Recommended requirements
### Appendix C: Installing Chips system-wide 
###### for system administrator, those who wish to share their Chips installation
### Appendix D: Installing the mdseqpos motif finder for chips
### Appendix E: Generating static reference files for Chips
- all of the required files
- using your own files
- supporting something new
- adding to ref.yaml
