# read2tree 

read2tree is a software tool that allows to obtain alignment matrices for tree inference. For this purpose it makes use of the OMA database and a set of reads. Its strength lies in the fact that it bipasses the several standard steps when obtaining such a matrix in regular analysis. These steps are read filtereing, assembly, gene prediction, gene annotation, all vs all comparison, orthology prediction, alignment and concatination. 




## Prerequisites

The following python packages are needed: [numpy](https://github.com/numpy/numpy), [scipy](https://github.com/scipy/scipy), [cython](https://github.com/cython/cython), [lxml](https://github.com/lxml/lxml), [tqdm](https://tqdm.github.io/docs/tqdm), [pysam](https://github.com/pysam-developers/pysam), [pyparsing](https://svn.code.sf.net/p/pyparsing/code/), [requests](http://python-requests.org), [filelock](https://github.com/benediktschmitt/py-filelock), [natsort](https://github.com/SethMMorton/natsort), [pyyaml](http://pyyaml.org/wiki/PyYAML), [biopython](https://github.com/biopython/biopython), [ete3](http://etetoolkit.org), [dendropy](http://packages.python.org/DendroPy/). You can install all of them using [conda](https://docs.conda.io/en/latest/miniconda.html).
```
conda install -c conda-forge biopython numpy Cython ete3 lxml tqdm scipy pyparsing requests natsort pyyaml
conda install -c bioconda dendropy 
```

Besides, you need softwares including [mafft](http://mafft.cbrc.jp/alignment/software/) (multiple sequence aligner), [iqtree](http://www.iqtree.org/) (phylogenomic inference), [ngmlr](https://github.com/philres/ngmlr), [ngm/nextgenmap](https://github.com/Cibiv/NextGenMap) (long and short read mappers), and [samtools](http://www.htslib.org/download/) which could be installed using conda.
```
conda install -c bioconda mafft  iqtree ngmlr nextgenmap  samtools
```

Finally, for installing [pyham](https://github.com/DessimozLab/pyham)(a library to work with HOGs), [pyoma](https://github.com/DessimozLab/pyoma)(library for retrieval of nucleotide sequences using OMA API) run `pip install pyham` `pip install pyoma` or alternatively:
```
git clone https://github.com/DessimozLab/pyham.git
python -m pip install -e ./pyham
git clone https://github.com/DessimozLab/pyoma.git
python -m pip install -e ./pyoma
```



## Installation

read2tree was built and tested with python 3.5.1. To set up read2tree on your local machine please follow the instructions below.

```
git clone https://github.com/DessimozLab/read2tree.git
cd read2tree
python setup.py install
```



## Run

To run read2tree two things are required as input:
1) The DNA sequencing reads as FASTQ file(s).
2) A set of reference orthologous groups, i.e. marker genes. 
This can be obtained from [OMA browser](https://omabrowser.org/oma/export_markers). 


```
read2tree --tree --standalone_path marker_genes/ --reads read_1.fq read_2.fq  --output_path output
```



## Test example


You can run the test example as follows:
```
cd tests
read2tree --standalone_path marker_genes/ --reads read_1.fq read_2.fq  --output_path output/
```
In the folder 'tests/output' you should be able to find the following folders:

| folder/file  | description           | 
| ------------- |-------------|
| 01_ref_ogs_aa | contains the selected OGs with amino acid data | 
| 01_ref_ogs_dna | contains the selected OGs with dna data |
| 02_ref_dna | contains the OGs reshuffeled by available species | 
| 03_align_aa | contains mafft alignment of aa data|
| 03_align_dna | contains codon replacement of aa alignments|
| 04_mapping_test_1b |contains the consensus sequences from the mapping|
| 05_ogs_map_test_1b_aa | contains the OGs with additional sequence test_1b|
| 05_ogs_map_test_1b_dna | contains the OGs with additional sequence test_1b|
| 06_align_test_1b_aa | contains the alignment with additional sequence test_1b|
| 06_align_test_1b_dna | contains the alignment with additional sequence test_1b|
| concat_test_1b_aa.phy | concatenated alignments from 06 amino acid folder|
| concat_test_1b_dna.phy| concatenated alignments from 06 dna folder|
| test_1b_all_cov.txt | summary of average numbers of reads used for selected sequences|
| test_1b_all_sc.txt | summary of average consensus length of reconstructed sequences|


Note that we consider species names as 5-letter codes e.g. AMPFI= Amphiura filiformis.

For running on clusters, you can run the first step of read2tree such that folders 01, 02 and 03 are computed (this allows for mapping). This can be done using the '--reference' option.  Since read2tree re-orders the OGs into the included species, it is possible to split the mapping step per species using multiple threads for the mapper. For this the '--single_mapping' option is available.


## Details of arguments

You can see the details of arguments of the read2tree package by running `read2tree -h`.
s
```
usage: read2tree [-h] [--version] [--threads THREADS] [--standalone_path STANDALONE_PATH] [--reads READS [READS ...]] [--read_type READ_TYPE]
                 [--split_reads] [--split_len SPLIT_LEN] [--split_overlap SPLIT_OVERLAP] [--split_min_read_len SPLIT_MIN_READ_LEN]
                 [--sample_reads] [--coverage COVERAGE] [--min_cons_coverage MIN_CONS_COVERAGE] [--genome_len GENOME_LEN]
                 [--output_path OUTPUT_PATH] [--dna_reference DNA_REFERENCE] [--ignore_species IGNORE_SPECIES] [--sc_threshold SC_THRESHOLD]
                 [--remove_species_mapping REMOVE_SPECIES_MAPPING] [--remove_species_ogs REMOVE_SPECIES_OGS]
                 [--ngmlr_parameters NGMLR_PARAMETERS] [--keep_all_ogs] [--check_mate_pairing] [--debug]
                 [--sequence_selection_mode SEQUENCE_SELECTION_MODE] [-s SPECIES_NAME] [--tree] [--merge_all_mappings] [-r]
                 [--min_species MIN_SPECIES] [--single_mapping SINGLE_MAPPING] [--ref_folder REF_FOLDER]

read2tree is a pipeline allowing to use read data in combination with an OMA standalone output run to produce high quality trees.

optional arguments:
  -h, --help            show this help message and exit
  --version             Show programme's version number and exit.
  --threads THREADS     [Default is 1] Number of threads for the mapping using ngm / ngmlr!
  --standalone_path STANDALONE_PATH
                        [Default is current directory] Path to oma standalone directory.
  --reads READS [READS ...]
                        [Default is none] Reads to be mapped to reference. If paired end add separated by space.
  --read_type READ_TYPE
                        [Default is short reads] Type of reads to use for mapping. Either ngm for short reads or ngmlr for long will be used.
  --split_reads         [Default is off] Splits reads as defined by split_len (200) and split_overlap (0) parameters.
  --split_len SPLIT_LEN
                        [Default is 200] Parameter for selection of read split length can only be used in combinationwith with long read option.
  --split_overlap SPLIT_OVERLAP
                        [Default is 0] Reads are split with an overlap defined by this argument.
  --split_min_read_len SPLIT_MIN_READ_LEN
                        [Default is 200] Reads longer than this value are cut into smaller values as defined by --split_len.
  --sample_reads        [Default is off] Splits reads as defined by split_len (200) and split_overlap (0) parameters.
  --coverage COVERAGE   [Default is 10] coverage in X.
  --min_cons_coverage MIN_CONS_COVERAGE
                        [Default is 1] Minimum number of nucleotides at column.
  --genome_len GENOME_LEN
                        [Default is 2000000] Genome size in bp.
  --output_path OUTPUT_PATH
                        [Default is current directory] Path to output directory.
  --dna_reference DNA_REFERENCE
                        [Default is None] Reference file that contains nucleotide sequences (fasta, hdf5). If not given it will usethe RESTapi
                        and retrieve sequences from http://omabrowser.org directly. NOTE: internet connection required!
  --ignore_species IGNORE_SPECIES
                        [Default is none] Ignores species part of the OMA standalone pipeline. Input is comma separated list without spaces, e.g.
                        XXX,YYY,AAA.
  --sc_threshold SC_THRESHOLD
                        [Default is 0.25; Range 0-1] Parameter for selection of sequences from mapping by completeness compared to its reference
                        sequence (number of ACGT basepairs vs length of sequence). By default, all sequences are selected.
  --remove_species_mapping REMOVE_SPECIES_MAPPING
                        [Default is none] Remove species present in data set after mapping step completed and only do analysis on subset. Input
                        is comma separated list without spaces, e.g. XXX,YYY,AAA.
  --remove_species_ogs REMOVE_SPECIES_OGS
                        [Default is none] Remove species present in data set after mapping step completed to build OGs. Input is comma separated
                        list without spaces, e.g. XXX,YYY,AAA.
  --ngmlr_parameters NGMLR_PARAMETERS
                        [Default is none] In case this parameters need to be changed all 3 values have to be changed [x,subread-length,R]. The
                        standard is: ont,256,0.25. Possibilities for these parameter can be found in the original documentation of ngmlr.
  --keep_all_ogs        [Default is on] Keep all orthologs after addition of mapped seq, which means also the OGs that have no mapped sequence.
                        Otherwise only OGs are used that have the mapped sequence for alignment and tree inference.
  --check_mate_pairing  Check whether in case of paired end reads we have consistent mate pairing. Setting this option will automatically select
                        the overlapping reads and do not consider single reads.
  --debug               [Default is off] Changes to debug mode: * bam files are saved!* reads are saved by mapping to OG
  --sequence_selection_mode SEQUENCE_SELECTION_MODE
                        [Default is sc] Possibilities are cov and cov_sc for mapped sequence.
  -s SPECIES_NAME, --species_name SPECIES_NAME
                        [Default is name of read 1st file] Name of species for mapped sequence.
  --tree                [Default is false] Compute tree, otherwise just output concatenated alignment!
  --merge_all_mappings  [Default is off] In case multiple species were mapped to the same reference this allows to merge this mappings and build
                        a tree with all included species!
  -r, --reference       [Default is off] Just generate the reference dataset for mapping.
  --min_species MIN_SPECIES
                        Min number of species in selected orthologous groups. If not selected it will be estimated such that around 1000 OGs are
                        available.
  --single_mapping SINGLE_MAPPING
                        [Default is none] Single species file allowing to map in a job array.
  --ref_folder REF_FOLDER
                        [Default is none] Folder containing reference files with sequences sorted by species.

read2tree (C) 2017-2022 David Dylus
```



## Change log


- version 0.2


- version 0.1

Adding covid analysis

- version 0.0

Initial work


## Authors

* [David Dylus](https://github.com/dvdylus), (main author)
* [Adrian Altenhoff](http://people.inf.ethz.ch/adriaal).
* [Sina Majidian](https://github.com/smajidian).

The authors would like to thank  Alex Warwick for help how to initiate such a package.

## License
This project is licensed under the MIT License.

