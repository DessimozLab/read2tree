# read2tree 

read2tree is a software tool that allows to obtain alignment matrices for tree inference. For this purpose it makes use of the OMA database and a set of reads. Its strength lies in the fact that it bipasses the several standard steps when obtaining such a matrix in regular analysis. These steps are read filtereing, assembly, gene prediction, gene annotation, all vs all comparison, orthology prediction, alignment and concatination. 

read2tree works in linux with  [![Python 3.10.8](https://img.shields.io/badge/python-3.10.8-blue.svg)](https://www.python.org/downloads/release/python-310/)


# Installation

There are three ways to install read2tree. You can choose either of them. 

### 1) Installation from source 

To set up read2tree on your local machine from source please follow the instructions below.


First, we need to create a fresh [conda](https://docs.conda.io/en/latest/miniconda.html) environment: 
```
conda create -n r2t python=3.10.8
```


#### Prerequisites

The following python packages are needed: [numpy](https://github.com/numpy/numpy), [scipy](https://github.com/scipy/scipy), [cython](https://github.com/cython/cython), [lxml](https://github.com/lxml/lxml), [tqdm](https://tqdm.github.io/docs/tqdm), [pysam](https://github.com/pysam-developers/pysam), [pyparsing](https://svn.code.sf.net/p/pyparsing/code/), [requests](http://python-requests.org), [filelock](https://github.com/benediktschmitt/py-filelock), [natsort](https://github.com/SethMMorton/natsort), [pyyaml](http://pyyaml.org/wiki/PyYAML), [biopython](https://github.com/biopython/biopython), [ete3](http://etetoolkit.org), [dendropy](http://packages.python.org/DendroPy/).  

You can install all of them using.
```
conda install -c conda-forge biopython numpy Cython ete3 lxml tqdm scipy pyparsing requests natsort pyyaml filelock
conda install -c bioconda dendropy pysam
```

Besides, you need softwares including [mafft](http://mafft.cbrc.jp/alignment/software/) (multiple sequence aligner), [iqtree](http://www.iqtree.org/) (phylogenomic inference), [ngmlr](https://github.com/philres/ngmlr), [ngm/nextgenmap](https://github.com/Cibiv/NextGenMap) (long and short read mappers), and [samtools](http://www.htslib.org/download/) which could be installed using conda.
```
conda install -c bioconda mafft iqtree ngmlr nextgenmap samtools
```

Then, you can install the read2tree package after downlaoding the package from this GitHub repo using

```
git clone https://github.com/DessimozLab/read2tree.git
cd read2tree
python setup.py install
```


### 2) Installation using Conda


```
conda create -n r2t python=3.10.8
conda install -c bioconda  read2tree 

```


### 3) Installation using Docker
The Dockerfile is also available in this repository. There is an example how to run in the [test example](#test-example) section.

A prebuild container can be loaded from dockerhub:
```
docker pull dessimozlab/read2tree:latest
```





## Run

To run read2tree two things are required as input:
1) The DNA sequencing reads as FASTQ file(s).
2) A set of reference orthologous groups, i.e. marker genes. 
In our wiki [page](https://github.com/DessimozLab/read2tree/wiki/obtaining-marker-genes), you may find information on how to obtain the marker genes using [OMA browser](https://omabrowser.org/oma/export_markers). 

Note that read2tree needs the Internet to download some additional files (cdna of OGs) from the oma database. 

### output 

The output of Read2Tree is the concatenated alignments as a fasta file where each record corresponds to one species. We also provide the option `--tree` for inferring the species tree using IQTREE as defualt.  


### Single species mode
```
read2tree --tree --standalone_path marker_genes/ --reads read_1.fastq read_2.fastq  --output_path output
```

### Multiple species mode
```
read2tree --standalone_path marker_genes/ --output_path output --reference  # this creates just the reference folder 01 - 03
read2tree --standalone_path marker_genes/ --output_path output --reads species1_R1.fastq species2_R2.fastq
read2tree --standalone_path marker_genes/ --output_path output --reads species2_R1.fastq species2_R2.fastq
read2tree --standalone_path marker_genes/ --output_path output --reads species3_R1.fastq species3_R2.fastq
read2tree --standalone_path marker_genes/ --output_path output --merge_all_mappings --tree
```


## Test example

The goal of this test example is to infer species tree for Mus musculus using its sequencing reads. You can download the full read data from from [SRR5171076](https://www.ncbi.nlm.nih.gov/sra/?term=SRR5171076) using [sra-tools](https://anaconda.org/bioconda/sra-tools). Alternatively, a small read dataset is provided in the `tests` folder. For this example, we consider five species including Mnemiopsis leidyi, Xenopus laevis, Homo sapiens, Gorilla gorilla, and Rattus norvegicus as the reference. Using [OMA browser](https://omabrowser.org/oma/export_markers), we downloaded 20 marker genes of these five species as the reference orthologous groups, located in the folder `tests/mareker_genes`. 

```
cd tests
read2tree --debug --tree --standalone_path marker_genes/ --reads sample_1.fastq sample_2.fastq --output_path output/
```


#### Run test example using docker

```
docker run --rm -i -v $PWD/tests:/input -v $PWD/tests/:/reads -v $PWD/output:/out -v $PWD/run:/run  dessimozlab/read2tree:latest  --tree --standalone_path /input/marker_genes --dna_reference /input/cds-marker_genes.fasta.gz --reads /reads/sample_1.fastq --output_path /out
```

### output files

You can check the inferred species tree for the sample and five reference species in Newick format:
```
$cat  output/tree_sample_1.nwk
(sample_1:0.0106979811,((HUMAN:0.0041202790,GORGO:0.0272785216):0.0433094119,(XENLA:0.1715052824,MNELE:0.9177670816):0.1141311779):0.0613339433,RATNO:0.0123413734);
```
For the full description of output files please check our wiki [page](https://github.com/DessimozLab/read2tree/wiki/output-files).  


Note that we consider species names as 5-letter codes e.g. XENLA = Xenopus laevis. If you want to rerun your analysis, make sure that you moved/deleted the files. Otherwise, read2tree continues the progress of previous analysis.  

For running on clusters, you can run the first step of read2tree such that folders 01, 02 and 03 are computed (this allows for mapping). This can be done using the '--reference' option.  Since read2tree re-orders the OGs into the included species, it is possible to split the mapping step per species using multiple threads for the mapper. For this the '--single_mapping' option is available.

Hint: As read2tree exploits the `progress` package, the user can benefit from continuing unfinished runs. However, if you want to conduct a new analysis with different inputs, you need to remove output of previous runs or change the `output_path`. 


### Details of arguments

To see the details of arguments, please take look at our wiki [page](https://github.com/DessimozLab/read2tree/wiki/Details-of-arguments)

## Possible issues

Installing on MAC sometimes drops this error:

```
raise ValueError, 'unknown locale: %s' % localename
ValueError: unknown locale: UTF-8
```

This can be mitigated using:

```
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
```

## Change log

- version 0.1.5:
  - fix issue with UnknownSeq being removed in Biopython>1.80
  - removing unused modeltester wrappers

- version 0.1.4:
   - allow reference folders not named marker_genes (#12)
   - update environment.yml file to contain all dependencies (#16)
   - documentation improvements
   - CI/CD pipeline

- version 0.1.3: 
   - improvements of documentation
   - adding support for docker
   - small bugfixes 

- version 0.1.2: packaging

- version 0.1.0: Adding covid analysis

- version 0.0: Initial work


## Authors

* [David Dylus](https://github.com/dvdylus), (main author)
* [Adrian Altenhoff](http://people.inf.ethz.ch/adriaal).


The authors would like to thank Alex Warwick for help how to initiate such a package.

## License
This project is licensed under the MIT License.

