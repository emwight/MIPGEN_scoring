MIPGEN_scoring
==============================

**Only compatible with -score_method "svr"**

An iterative algorithmic approach to optimize uniform coverage using MIPgen scoring. In this section, any references to "MIPgen" refer to the original code/documentation; references to "MIPgen_scoring" refer to modifications

The algorithm minimizes loss of uniform coverage, as determined by the probability of 0 reads in each base pair along the target. We consider each MIP as a Poisson process, the rate of which is determined by the given MIP's relative expected number of reads (as determined by MIPgen's SVR model) to all MIPs in the panel. We then consider each base pair along the target as conditional sampling from the Poisson process (i.e. transformed into a Multinomial distribution). The probability of each base pair having 0 reads is then calculated from a binomial pdf with parameters (p = the sum of rates for MIPs covering the given base pair; N = the user-provided number of expected reads we want to achieve uniformly with the panel). Finally via conditional independence, we assume little to no dependence between these probabilities, and so we sum them up across all base pairs in the target to get an overall loss of uniform coverage for the given panel.

(once algorithm is finalized, discuss here)

The below documentation from MIPgen is edited in **bold** with differences/additions from MIPgen_scoring.

MIPGEN
======

One stop MIP design and analysis

Use MIPgen to design custom mip panels for target enrichment of moderate to high complexity DNA targets ranging from 120 to 250bp in size

To compile MIPgen, you will need a C++ compiler, such as GCC (http://gcc.gnu.org/). That's it! For running design and analysis, see below for dependencies

**MIPgen dependencies: BWA, tabix, samtools**

**MIPgen_scoring additional dependency: GNU Scientific Library (GSL)**

**(For those newer to C++ like me, I installed it this way in Ubuntu:)**

    sudo apt-get install libgsl-dev

-----
BUILD MIPGEN
-----

To build MIPgen from source, simply enter the MIPGEN directory and type 'make'

The 'mipgen' executable can be used to perform designs

-----
MIPGEN PARAMETERS
-----

MIPgen accepts many options intended to maximize the flexibility of your designs

The only required parameters are the minimum and maximum sizes for the captured sequence, a prefix for output files and a BED file consisting of the coordinates of the desired targets, and the bwa- and samtools-indexed reference genome from which to pull the sequences

**NOTE: index your genome before running MIPgen or MIPgen_scoring**

    bwa index [reference_genome].fasta

Other options provide control over tiling, scoring and oligo design behavior

MIPgen offers three scoring methods based on two methods of scoring: SVR and logistic regression. Both methods offer similar performance, with slight improvements seen with the SVR scoring. The 'mixed' scoring option will perform designs primarily with the faster logistic regression score and finish with the more accurate SVR scoring scheme.

**Again note that MIPgen_scoring is only compatible with -score_method "svr"**

Run ./mipgen -doc to see the full documentation or read the text below

-----
MIPGEN PITFALLS
-----

-It is not recommended to include more than 8bp total of degenerate smMIP tags. Longer tags reduce specificity and complicate MIP rebalancing, leading to less complete coverage of targeted sites.

-Priority scores (used by default) modulate tiling behavior and density of tiling of low scoring regions. To reduce density, lower priority scores. In particular, for interrogation of very short loci (<10bp), priority scores should be set to 0 to prevent splitting the region into two MIP targets.

-Regions that do not map uniquely to the provided reference are not tiled by default. These and other bases that cannot be tiled will be reported in the <project_name>.coverage_failed.bed file.

-SVR scores below 1.4 and logistic scores below 0.7 rarely perform adequately. It is recommended to remove these probes from designed assays.

-The failure flags in the design files are a series of three bits. They represent, in order, non-uniqueness of the MIP target with respect to the provided genome sequence, inability to design SNP MIPs to any/all of the variants provided in the SNP file, and tandem repeats being detected in the MIP arms using Tandem Repeats Finder. SNP failures will be included in final probe sets. The others will not, assuming non-uniqueness is selected and a TRF executable is provided

-MIPgen does not remove candidate MIPs with overly high targeting arm copy numbers if there are no suitable alternatives to cover the desired regions. Be advised that copy numbers above 200 are likely to poisson your assay and should be removed unless you are specifically interrogating all such low complexity sites.

-----
OTHER NOTES
-----

-the subset of Boost required for compilation is included in the directory

-the makefile should compile the source without any external dependencies, but the executable that is generated makes system calls to SAMtools, BWA, Tabix, and, optionally, Tandem Repeats Finder

-the tools subdirectory has scripts for parsing fastqs prior to mapping and collapsing tag-defined read groups (TDRGs) and trimming targeting arms after mapping (Cython required)

-summary information is also printed as part of the collapsing process and requires NumPy and SciPy

-if there is at least 10bp of overlap between forward and reverse reads, it is recommended to merge read pairs into fr-reads using PEAR ( http://bioinformatics.oxfordjournals.org/content/early/2013/11/10/bioinformatics.btt593.long ) prior to processing and to utilize single end options

-see the HemoMIPs 2020 analysis pipeline for analyzing MIP sequencing data with html reports and full genotyping results: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007956
 
-----
GETTING STARTED
-----

Below is an example for getting MIPgen up and running and designing probes to EGFR, TERT and BRAF. 
The output generated from these commands is included in the repository (mipgen_example.tgz). 
A copy of the human reference genome (GRCh37, also known as hg19 in UCSC) is available at the 1000Genomes FTP site:

ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz

**Again, need to index the reference gene before running MIPgen or MIPgen_scoring:**

    bwa index human_g1k_v37.fasta

The file used in the example to pull out gene exons for the hg19 human reference (refGene.txt) can be downloaded here: 

http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/


First extract the tar file to make a MIPgen directory

    tar -xf mipgen.tgz

Enter the directory to compile the source code

    cd MIPGEN/

Use make to compile -- make sure you have a C++ compiler

    make

Make a directory for practicing designs

    mkdir ../mipgen_practice
    cd ../mipgen_practice

Create a new file and type a list of gene symbols you would like to tile
with MIPs

    vim practice_genes.txt

You can use one of the scripts in the tools subdirectory and a refGene text
file listing gene exons to pull out the genomic coordinates of your genes of
interest

    sh ../MIPGEN/tools/extract_coding_gene_exons.sh practice_genes.txt ../refGene.txt > practice_genes.bed

Check to make sure the coordinates are what you expect and that no errors
occurred

    less practice_genes.bed

Now you can perform designs with MIPgen!
Here is a design with very basic options.
Set the path to the genome reference and change the name of the FASTA file if yours is different.
Make sure you have the dependencies installed or accessible through a given
path (BWA, tabix, samtools, **for MIPgen_scoring: GSL**)!

**MIPgen_scoring: additional arguments are -initial_panel and -N_reads (described in more detail in the next section).**
**If no -initial_panel is supplied, MIPgen_scoring will just run MIPgen as created by Evan Boyle.**
**(dev goal: if no -initial_panel supplied, run MIPgen, then automatically use that as the initial panel for MIPgen_scoring)**

**To run MIPgen:**

    ../MIPGEN_scoring/mipgen -regions_to_scan practice_genes.bed -project_name practice_design -min_capture_size 162 -max_capture_size 162 -bwa_genome_index /<path to genome reference fasta file>/human_g1k_v37.fa

**To run MIPgen_scoring:**

    ../MIPGEN_scoring/mipgen -regions_to_scan practice_genes.bed -project_name practice_design -min_capture_size 162 -max_capture_size 162 -bwa_genome_index /<path to genome reference fasta file>/human_g1k_v37.fa -score_method "svr" -initial_panel practice_design_mipgen.csv -N_reads 50

**For MIPgen:** The final selection of MIPs is located in the picked MIPs file
review the scores to make sure the MIPs stand a good chance of success
(logistic scores below 0.6 are unlikely to provide usable data)

    less practice_design.picked_mips.txt

**Also for MIPgen:** By default all tested MIPs are output; this is a lot of output! (Turn it off with the silent_mode option)

    rm practice_design.all_mips.txt

**MIPgen_scoring:** **(in development; not in mipgen_example.tgz)** A summary of all tested panels are output in the all panels file
with the loss of uniform coverage and the random perturbation for each panel (the initial panel will have NA for random perturbation)

    less practice_design.all_panels.txt

**MIPgen_scoring:** **(in development; not in mipgen_example.tgz)** The final selection of MIPs is located in the final panel file

    less practice_design.final_panel.txt

sai, fq (and sam) files are not deleted automatically and can also take up
space

    rm -f *.sai *.fq

Generate a UCSC track with another tools script to visualize online

**(note: this is from MIPgen, but I had trouble getting this to work; still need to troubleshoot)**

    python ../MIPGEN/tools/generate_ucsc_track.py practice_design.picked_mips.txt practice_ucsc_track

Look at the other files in this directory to see designs for TERT, BRAF and
EGFR (we have not tested these designs experimentally so we cannot precisely
assess predicted performance for these probes)

-----
DESCRIPTION OF ALL OPTIONS
-----
**MIPgen_scoring additional options:**

    Required parameters to run MIPgen_scoring:

    -initial_panel                  path to csv file of an initial MIP panel for the iterative algorithm,
                                    described in more detail directly below;
								    if not provided, will default to running MIPgen as published by Boyle
    -score_method                   this is the same as MIPgen's, but MIPgen_scoring requires -score_method "svr"
                                    
    Other options:

    -N_reads                        desired number of reads for which we want to achieve uniform coverage
								    default is 100

The file for the initial panel should have one row per MIP with the following columns, in the specified order:
- **mip_start**: integer for the start position of the MIP scan sequence (between the extension and ligation arms) within the indexed genome.
- **mip_end**: integer for the end position of the MIP scan sequence (between the extension and ligation arms) within the indexed genome.
- **ext_len**: integer for the length of the extension arm of the MIP
- **lig_len**: integer for the length of the ligation arm of the MIP
- **probe_strand**: takes values + or -

See the mipgen_scoring_example for a sample .csv file that goes with Ellen's run of the MIPGEN example in their original documentation, as well as code to take an output .txt file of picked MIPs from MIPgen and convert to the required .csv format.

**MIPgen options (all of which are stil technically included in MIPgen_scoring):**

**Any required parameter in MIPgen is also a required parameter in MIPgen_scoring**

    usage: mipgen (<-parameter_name> <parameter_value>)*
    Created by Evan Boyle (boylee@uw.edu)
    Required parameters:

    -project_name                   prefix for output
                                    ex: /my/output/directory/batch_1
    -bwa_genome_index               samtools- and bwa-indexed genome reference file path
    -regions_to_scan                BED file of positions GUARANTEED to be captured
                                    regions will be merged as needed
    -min_capture_size               integer value for length of targeting arms plus insert region to be captured
                                    tested down to 120
    -max_capture_size               integer value for length of targeting arms plus insert region to be captured
                                    tested up to 250
    Oligo control:

    -arm_lengths                    manual setting of individual possible ligation and extension arm lengths
                                    format: <extension arm length 1>:<ligation arm length 1>,<extension arm length 2>:<ligation arm length 2>,...
                                    ex: 16:24,16:25,16:26
    -arm_length_sums                setting of arm lengths to include all pairs of arm lengths adding to a certain length
                                    default is 40,41,42,43,44,45
                                    format: <sum of extension and ligation arm lengths 1>,<sum of extension and ligation arm lengths 2>,...
                                    ex: 40,45
    -ext_min_length                 minimum length of extension arm
                                    default is 16
    -lig_min_length                 minimum length of ligation arm
                                    default is 18
    -tag_sizes                      specifies degenerate tag (Ns) to place in extension arm and ligation arm
                                    more than 8 bases severely reduces on-target rate of capture
                                    default is 5,0
                                    format: <extension arm tag size>,<ligation arm tag size>
                                    ex: 4,4
    -masked_arm_threshold           fraction of targeting arms having masked bases (floating point number from 0 to 1)
                                    to tolerate before optimizing with respect to masking (requires TRF executable)
                                    default is 0.5
    -target_arm_copy                threshold over which minimizing copy number to the reference takes priority
                                    default is 20
    -max_arm_copy_product           maximum permissible product of targeting arm copy number to the reference
                                    default is 75
    Tool dependencies:

    -tabix                          tabix build
                                    default is "tabix"
    -bwa                            bwa build
                                    default is "bwa"
    -trf                            tandem repeats finder executable
                                    default is "off"
                                    providing a TRF executable enables filtering of simple repeats in probe arms
    Input options:
    -genome_dir                     genome reference directory helpful for large jobs
                                    ex: /my/genome/split/by/chromosomes/
    -snp_file                       path to vcf file of snps to avoid during design
    -common_snps                    providing "off" will disable loading of common SNPs to avoid from NCBI using Tabix
                    DEPRECATED:ONLY SNP FILE OPTION NOW SUPPORTED
    -file_of_parameters             file containing any of the above parameters in the following format:
                                            -first_parameter_name first_parameter_value
                                            -second_parameter_name second_parameter_value
                                            <etc.>
                                            lines that do not start with '-' are ignored
    Tiling control:

    -feature_flank                  integer value for extending BED coordinates
                                    default is 0
    -capture_increment              integer step size for examining capture sizes
                                    default is 5
                                    starts at maximum and descends until minimum is reached or
                                    acceptable score found, whichever is first
    -logistic_heuristic             providing "off" will remove assumptions to reduce search space
    -max_mip_overlap                integer specifying the maximum number of nucleotides of overlap
                                    to test before selecting a MIP
                                    default is 30
    -starting_mip_overlap           integer specifying the starting number of nucleotides of overlap
                                    to use when linearly tiling with MIPs
                                    default is 0
    -check_copy_number              providing "off" enables targeting of non-unique sites in the reference
    -seal_both_strands              providing "on" disallows targeting arms occupying the same base
                                    positions on opposite strands
    -half_seal_both_strands         providing "on" disallows a second targeting arm from occupying more than half
                                    of the positions of a targeting arm placed on the other strand
    -double_tile_strand_unaware     providing "on" will ensure that all positions in targeted regions will be tiled twice
    -double_tile_strands_separately providing "on" will ensure that all positions in targeted regions will be tiled twice,
                                    at least once on each strand
    Scoring parameters:

    -score_method                   "logistic" will use the logistic model for scoring (default)
                                    "svr" will switch behavior to score with the svr model
                                    "mixed" will first filter with logistic scoring and finalize with the svr model
    -logistic_optimal_score         value of score metric at which to stop optimizing MIP score and
                                    accept current MIP, lower scores lead to fewer outliers
                                    default is 0.98
    -svr_optimal_score              analogous to logistic_optimal_score
                                    default is 2.2
    -logistic_priority_score        value of score metric below which MIPs are placed nonlinearly to
                                    optimize score and SNPs are tolerated
                                    default is 0.9,lower values enable more efficient linear tiling
    -svr_priotity_score             analogous to logistic_priority_score
                                    default is 1.5 for SVR
    Miscellaneous:

    -silent_mode                    providing "on" will reduce volume of text output
    -download_tabix_index           providing "on" will force redownload of common snp tbi file
                                        DEPRECATED
    -bwa_threads                    make use of BWA's multithreading option (-t _)
                                    default is 1

-----
TERMS OF USE AND LICENSING INFORMATION
-----

Copyright 2014 University of Washington.
 
Permission is hereby granted, to You to use, copy, modify, merge the MIPgen software and associated documentation files (the "Software") for the use in Your non-commercial academic and research activities. This includes electronic distribution for the use of client software, but only within Your organization. You may not otherwise distribute the Software, including any modified versions. You may not sublicense any rights under this license and You may not sell copies of the software.
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,  DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 
Contact license@uw.edu for any uses not covered under this license.

© University of Washington 2014
