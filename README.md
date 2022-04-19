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
python -m pip install -3 ./pyham
git clone https://github.com/DessimozLab/pyoma.git
python -m pip install -3 ./pyoma
```



## Installation

read2tree was build and tested with python 3.5.1. To set up read2tree on your local machine please follow the instructions below.

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
read2tree --standalone_path marker_genes/ --reads read_1.fq read_2.fq  --output_path output
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

