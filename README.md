## Introduction

Leishmania major kDNA annotation is a python 3.7 package (modified from kDNA annotation) to annotate minicircles of mitochondrial DNA, also known as kinetoplast DNA (kDNA). It identifies guide RNAs (gRNAs) genes by alignment. Alignment of gRNAs with their cognate mRNAs are found. If gRNA transcriptomics data is available, the expression status of gRNAs are predicted (using kDNA annotation).

## Installation 

It is recommended to have anaconda3 installed. kDNA annotation uses Python3.7 packages pandas, numpy, biopython and yaml. kDNA annotation also requires the C libraries libyaml and the openMP library libgomp1. The package comes with binaries, but if you want to compile the C routines you will need gcc and libyaml-dev.

To clone the latest version of kDNA annotation use

    git clone https://github.com/Zedthedrifter/Leishmania_major_kDNA_annotation.git

Set up the PATH and PYTHONPATH environment variables in `.bash_profile`: 

    export PYTHONPATH=$HOME/Leishmania_major_kDNA_annotation:$PYTHONPATH
    export PATH=$HOME/Leishmania_major_kDNA_annotation/Leish_kDNA_annotation:$PATH


## Usage

The complexity and variablility of minicircles within and between species means that minicircle annotation of cassettes and gRNAs cannot be completely automated. User intervention is required throughout the annotation pipeline.

The pipeline commands are contained in the file `share/pipeline.py`. Alternatively, the pipeline can be done in the jupyter notebook `share/pipeline.ipynb`. This may be handier if the diagnostic output from the pipeline needs to be available. The pipeline is directed by the contents of the configuration file `config.yaml`, an example of which is in the `share/` directory. Some sequence data are also available in the '/share' directory.

kDNA annotation requires, at minimum, the following fasta files:
- minicircle sequences
- maxicircle sequence
- unedited mRNA sequences
- edited mRNA sequences

If transcriptomics of gRNAs are available, a gzipped SAM file of transcripts mapped to minicircles is needed.

#### Step 1

For each project set up the following directory structure:
```
    Project_name/
        pipeline.py
        config.yaml
        Work_files/
        In_files/
            minicircles file.fasta
            maxicircles file.fasta
            edited mRNA file.fasta
            unedited mRNA file.fasta
            transcripts file.sam.gz
```
The files `pipeline.py` and `config.yaml` should be copied from `kDNA_annotation/share/`. In the pipeline one additional directories will be created `Annotation/`.
to run pipeline.py, config file and work directory should be specified as inputs. This also allows you to analyse multiple strains within the same project, with individual config and work directory (to avoid overwritting)

#### Step 2

In `config.yaml` change
- the `project` parameter to the full path name of the project,
- the `minicircle fasta infile` parameter to the file name of the minicircle fasta file in the directory Infile,
- the `maxicircle fasta infile` parameter to the file name of the maxicircle fasta file in the directory Infile,
- the `unedited mRNA fasta infile` parameter to the file name of the unedited mRNA fasta file in the directory Infile,
- the `edited mRNA fasta infile` parameter to the file name of the edited mRNA fasta file in the directory Infile,
- if you have transcriptomics, the `transcriptomics infile` parameter to the file name of the transcriptomics sam file in the directory Infile,

#### Step 3

Run prepare(config_file) part of pipeline.py: clean_mini_and_maxicircles(),lm_mRNA_process(), align_maxi, align_mini

Clean the minicircle and maxicircle fasta files with `clean_mini_and_maxicircles()`.

Run `clean_mini_and_maxicircles()` in `pipeline.py` by un-commenting the call to `clean_mini_and_maxicircles()` and commenting out calls to all the other functions.

For kDNA annotation to work correctly all minicircles should be unique, aligned on conserved sequence block 1 (CSB-1) and have CSB-3 about 100nt downstream from CSB-1. Minicircles should also have the standard name "mO_XXX" where XXX runs from 001 upwards. 

The known CSB-1, 2 and 3 sequences as regular expressions are given in the config file. 

If there are no matches for CSB-1 in a minicircle's sequence this could because the minicircle has not been assembled completely, the regular expression is wrong or that a minor variant of the CSB-1 sequence exists. The program will output possible new variants that it finds at the start the minicircle sequences and then terminate. These variants can be added to the CSB-1 regular expression in the config file. On the other hand, if you want to automatically remove any minicircles without CSB-1 then set `remove minicircles without CSB1` in `config.yaml` to yes.

Minicircles containing CSB-1 but not aligned on it will be automatically re-aligned and a warning issued.

If CSB-3 does not exist in the correct position relative to CSB-1 a warning is given and the program terminates.

If any minicircle has an invalid name, all minicircles are re-named in the standard format mO_XXX (starting with 0). A log file is written out with the the old and new names of the minicircles. The name of the log file is given by `clean log text file` in the config file.

The maxicircle is also re-named `Maxicircle`.

The cleaned mini and maxicircle fasta files are re-written as `minicircle clean fasta file` and `maxicircle clean fasta file` respectively in the `working directory` as given in the config file.

#### Step 4

Process the edited and unedited mRNA sequences with `mRNA_process()`.

You may comment out `clean_mini_and_maxicircles()` and uncomment `mRNA_process()`.

The function `mRNA_process.py` will attempt to clean the inputted edited and unedited sequences given as fasta files in  `in directory` defined in the config file. It is important that the T-stripped edited and unedited sequences are identical otherwise it is usually impossible to reconstruct the insertions and deletions correctly. If this problem occurs, the offending sequences should be corrected. If the mis-match between the sequences occurs in a non-edited region this won't be a problem for gRNA identification. 

Several fasta files are outputted in `working directory`, including positions where deletions occur (`deletions mRNA text file`) and positions where insertions occur which are recorded as a lowercase "u" (`edited_mRNA_small_u.fasta`) or "t" (`edited_mRNA_small_t.fasta`).

#### Step 5

Align the maxicircle (`maxi_align`) and all the minicircles (`mini_align`) to all edited mRNA sequences. This code is written in parallelised C for speed. Both sense and anti-sense strands of the circles are aligned. 

The source code for these binaries are in the `src` directory of the package. If the source code is changed the binaries need re-compiling with the commands

    gcc -O3 -o ../kDNA_annotation/align_maxi align_maxi.c -fopenmp -lyaml
    gcc -O3 -o ../kDNA_annotation/align_mini align_mini.c -fopenmp -lyaml

These programs find all alignments at least as long as `minimum gRNA length` (defined in the config file). As well as finding all gRNAs, these programs also produce many false positives, duplicates and overlapping alignments. In addition, no checks are made for an anchor, nor whether an alignment contains insertions or deletions in the mRNA. These checks are all done at a later step. These programs output to two files `maxi_alignments.txt` and `mini_alignments.txt`. 

#### Step 6

Find high quality canonical gRNAs with `lm_hq_gRNAs(config_file)`.

As Leishmania does not have inverted repeats, we need to find the high quality canonical gRNAs; those that we can be quite sure are actual gRNAs. 

High quality gRNAs should be long (>= 35 nt) with a minimal number of mis-matches (<= 1), have a long anchor (>= 8 nt) of only Watson-Crick basepairs. These parameters are defined in the config file under `high quality gRNAs filter`.
The position of gRNA can be inferred and specified: 'start_position' 'end_position'. After specifying the expected region for gRNA, the function can be rerun with less stringent filter for high quality gRNA.

On running the function `hq_gRNAs()`, the alignments from Step 5 will be loaded and the high quality gRNAs extracted. The number of high quality gRNAs will be printed and the number of minicircles that contain these gRNAs. The gRNA lengths are plotted as a histogram. The histogram should be right-skewed with lengths up into the 50s. Also plotted as a histogram are the 5' positions of the high quality gRNAs on the minicircles. These should be clustered into one or more peaks which represent the cassettes. 

If no, or very few, high quality gRNAs are found then the filtering parameters in the config file should be weakened and  `hq_gRNAs()` re-run.

#### Step 7

Run 'lm_feature_id()' `annotate_minicircles()` to output a genbank file of all the minicircles and their annotated features in `genbank directory/genbank text file`and alignments of canonical gRNAs to edited mRNA in the `alignments directory`.  There is one file per mRNA.

#### Step 18

Analysis. TODO.

## Description of output files in Annotation directory

All positions of features on minicircles and mRNAs are 0-indexed in the following output files. However, in the genbank file, and in the names assigned to gRNAs,  positions are 1-indexed.

3' end positions of features are exclusive, i.e., the feature finishes at the position just before the recorded position.

In the following tables an asterisk denotes columns only occurring when transcriptomics data is available.

### cassettes_with_expression.txt and cassettes.txt

Contains information about each and every cassette identified and labelled on all minicircles.

column | description
:--- | :---
mO_name | Name of the minicircle the cassette is encoded on
forward_start | 5' position of forward inverted repeat on the minicircle
forward_end | 3' end of forward inverted repeat on the minicircle (exclusive)
forward_seq | The sequence of the forward repeat
reverse_start | 5' position of reverse inverted repeat on the minicircle
reverse_end | 3' end of reverse inverted repeat on the minicircle (exclusive)
reverse_seq | The sequence of the reverse repeat
cassette_label | The assigned cassette label (aka position)
expression* | **values = {expressed, non-expressed}**. The expression status of gRNAs in the cassette in available
type | **values = {canonical, non-canonical}**. Whether the cassette contains at least one canonical gRNA, or no canonical gRNAs

### gRNAs_with_expression.txt and gRNAs.txt

Contains information about each and every canonical gRNA found by alignment of minicircles and the maxicircle to edited mRNA. This includes all gRNAs in cassettes and orphan gRNAs outside of cassettes. A cassette may contain multiple gRNAs that edit different positions on mRNAs. Guide RNAs that edit the same region of different versions of the same mRNA gene are also included. 

column | description
:--- | :---
mO_name | Name of the minicircle (or maxicircle) the gRNA is encoded on
cassette_label | The cassette label or position the gRNA is encoded in. Equals "Maxi" if encoded on the maxicircle, or "Orphan" if not encoded in a cassette.
strand | **value = {template, coding}**. The strand the gRNA is encoded on.
length | The length of the gRNA from the start of the anchor to the end of the guiding region.
rel_start | For coding strand gRNAs: the distance from the 3' end of the forward repeat to the 5' end of the anchor. For template strand gRNAs: the distance from the 3' end of the anchor to the 5' end of the reverse repeat.
circle_start | 5' position of the gRNA on the minicircle
circle_end | 3' position of the gRNA on the minicircle (exclusive)
mRNA_name | The name of the mRNA the gRNA edits (includes mRNA version number if necessary)
product | The name of the mRNA the gRNA edits (does not include mRNA version number)
mRNA_start | 5' position on the edited mRNA the gRNA edits 
mRNA_end | 3' position on the edited mRNA the gRNA edits (exclusive)
mRNA_seq | The mRNA sequence edited by the gRNA (5' to 3')
gRNA_seq | The reversed gRNA sequence (i.e., 3' to 5') from 3' end of guiding region to 5' end of anchor
pairing | Base-pairing between the mRNA and gRNA ("\|" = Watson-Crick, ":" = GU wobble, "." = mismatch)
mismatches | Number of mismatches between the mRNA-gRNA duplex
name | The given name of the gRNA 
transcripts_total* | Total number of transcripts of the same strand as the gRNA that map to the gRNA on the minicircle. For cassette-encoded gRNAs this includes all transcripts whose 5' ends map within the cassette irrespective of whether they overlap with the gRNA.
transcripts_init_site* | Total number of transcripts whose 5' ends map within the initiation site (positions 30, 31 and 32 for *T. brucei*) of the gRNA's cassette. For orphan gRNAs this equals "transcripts_total"
p-value* | The probability that the number of transcripts mapping thr initiation site could have arisen by random chance under the null model that the gRNA is not expressed. Small p-values suggest rejection of the null model
expression* | **values = {expressed, non-expressed}**. Whether the gRNA is considered expressed or not depending on "p-value"
gene_rel_start* | For coding strand gRNAs in cassettes: the distance from the 3' end of the forward repeat to the 5' end of the initiation site. For template strand gRNAs in cassettes: the distance from the 3' end of the initiation site to the 5' end of the reverse repeat. Orphan gRNAs: distance to the 5' position on the minicircle with the most number of mapped cognate transcripts. Maxicircle gRNAs: set to "\<NA\>".
initiation_sequence* | 5 nucleotide initiation sequence of the gRNA gene. This starts at the 5' (3' for template strand gRNAs) position on the minicircle  with the most number of mapped cognate transcripts. Non-expressed gRNAs will have an initiation sequence if transcripts map to the initiation site, otherwise set to "NaN".
gene_rel_end* | For coding strand gRNAs in cassettes: the distance from the 3' end of the forward repeat to the 3' end of the gene. For template strand gRNAs in cassettes: the distance from the 5' end of the gRNA gene to the 5' end of the reverse repeat. For Orphan gRNAs this is relative to "gene_rel_start", i.e., the length of the gene. Maxicircle gRNAs: set to "\<NA\>".
rel_pos* | Distance from the 5' end of the initiation sequence to the 5' end of the anchor (3' for template gRNAs). Maxicircle gRNAs: set to zero.
anchor_type | **values = {unanchored, initiator, extenderA, extenderB}**. Unanchored means the 5' gRNA that creates the anchor for this gRNA is missing. Initiator means the anchor for this gRNA is fully within an unedited region. ExtenderA means that the 5' gRNA exists and creates the anchor region for this gRNA. ExtenderB means that the minimum anchor length (e.g. 6nt) is fully within an unedited region but that the full anchor includes insertions not covered by a previous gRNA. So gRNA maybe an initiator or an extender with missing 5' gRNA.
gene_mRNA_end* | The 3' position on the mRNA of the family this gRNA belongs to based on coincident mapping of 5' ends of genes
family_no | A number identifier for the family this gRNA belongs to. The number ascends for families running 5' to 3' along the mRNA (not used)
family_end | The 3' position on the mRNA of the family this gRNA belongs. Equals "gene_mRNA_end" is transcriptomics is available. otherwise equals the coincident ends of anchors.
family_id | The unique id of the family this gRNA belongs to



### genes_with_expression.txt

Only produced if transcriptomics is available. Gives information about each and every expressed gRNA gene (canonical and non-canonical, but not Maxicircle encoded gRNA genes) predicted by transcript mapping to minicircles. This includes all genes in cassettes and orphan genes outside of cassettes. A cassette may contain multiple genes that edit different positions on mRNAs. Genes of gRNAs that edit the same region of different versions of the same mRNA gene are also included. 


column | description
:--- | :---
mO_name | Name of the minicircle the gene is encoded on
cassette_label | The cassette label or position the gene is encoded in. Equals "Orphan" if not encoded in a cassette.
strand | **value = {template, coding}**. The strand the gRNA is encoded on.
transcripts_total | Total number of transcripts of the same strand as the gRNA that map to the gRNA on the minicircle. For cassette-encoded gRNAs this includes all transcripts whose 5' ends map within the cassette irrespective of whether they overlap with the gRNA.
transcripts_init_site | Total number of transcripts whose 5' ends map within the initiation site (positions 30, 31 and 32 for *T. brucei*) of the gRNA's cassette. For orphan gRNAs this equals "transcripts_total"
p-value | The probability that the number of transcripts mapping thr initiation site could have arisen by random chance under the null model that the gRNA is not expressed. Small p-values suggest rejection of the null model
rel_start | For coding strand gRNAs: the distance from the 3' end of the forward repeat to the 5' end of the initiation sequence. For template strand gRNAs: the distance from the 3' end of the initiation sequence to the 5' end of the reverse repeat.
initiation_sequence | 5 nucleotide initiation sequence of the gRNA gene. This starts at the 5' (3' for template strand gRNAs) position on the minicircle  with the most number of mapped cognate transcripts.
rel_end | For coding strand gRNAs in cassettes: the distance from the 3' end of the forward repeat to the 3' end of the gene. For template strand gRNAs in cassettes: the distance from the 5' end of the gRNA gene to the 5' end of the reverse repeat. For Orphan gRNAs this is relative to "rel_start", i.e., the length of the gene.
mRNA_name | The name of the mRNA the gRNA edits (includes mRNA version number if necessary). Equal to "NaN" for non-canonical gRNA genes.
circle_end | 3' position of the gene on the minicircle (exclusive)
circle_start | 5' position of the gene on the minicircle
gRNA_seq | The reversed gRNA gene sequence from 5' end of initiation sequence to 3' end of gene (3' to 5')
length | The length of the gene from the 5' end of the initiation sequence to the 3' end of the gene
mRNA_end | 3' position on the edited mRNA the gene edits (exclusive). Equal to "\<NA\>" for non-canonical gRNA genes
mRNA_seq | The mRNA sequence edited by the gRNA gene (5' to 3'). Equal to "NaN" for non-canonical gRNA genes.
mRNA_start | 5' position on the edited mRNA the gRNA edits. This includes the initiation sequence. Equal to "\<NA\>" for non-canonical gRNA genes
pairing | Base-pairing between the mRNA and gRNA ("\|" = Watson-Crick, ":" = GU wobble, "." = mismatch)
type | **values = {canonical, non-canonical}**. Whether the gene is is canonical or non-canonical

### Alignment files

Lines | Description
:--- | :---
1-4 | Position along mRNA (0-indexed)
5 | "M": Edited positions that are not covered by gRNAs (i.e., shows positions where gRNAs are missing)
6 | "@": Positions of the 5' ends of expressed canonical gRNA genes, i.e., from 5' end of initiation sequence (expressed genes are shown as lines of minus signs under each aligned gRNA).
7 | "X": Positions of the 5' ends of canonical gRNAs found by alignment to edited mRNA
8 | The number of deleted "U"s to the right of the base it is over
9 | Lowercase "u"s are insertions
10 | Protein sequence 

