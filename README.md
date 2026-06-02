# read2tree 

read2tree is a software tool that allows obtaining alignment matrices for tree inference. For this purpose it makes use of the OMA database and a set of reads. Its strength lies in the fact that it bypasses the several standard steps when obtaining such a matrix in regular analysis. These steps are read filtering, assembly, gene prediction, gene annotation, all vs all comparison, orthology prediction, alignment and concatenation. 

read2tree works in linux with  [![Python 3.10.8](https://img.shields.io/badge/python-3.10.8-blue.svg)](https://www.python.org/downloads/release/python-310/)

## New release 
We are now releasing Read2Tree v2.0.0 with improved speed and logging. As the aligner, we are now using minimap2 [minimap2](https://github.com/lh3/minimap2). Also, MAFFT and IQtree are now using multiple threads. We would suggest running r2t with `--debug` which helps to debug later. Please note that arguments have slightly changed in this release(see below for details). 


## Read2Tree Talk:
You can watch David Dylus's presentation on Read2Tree as part of the SIB [in silico talks](https://www.sib.swiss/in-silico-talks/read2tree-inferring-phylogenetic-trees-from-raw-sequencing-data).

## Read2Tree publication
You can cite Read2Tree published in [Nature Biotechnology](https://doi.org/10.1038/s41587-023-01753-4): 
```
David Dylus, Adrian Altenhoff, Sina Majidian, Fritz J. Sedlazeck & Christophe Dessimoz.

Inference of phylogenetic trees directly from raw sequencing reads using Read2Tree.  Nat Biotechnol (2023). https://doi.org/10.1038/s41587-023-01753-4.
```


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
conda install -c conda-forge biopython numpy Cython ete3 lxml tqdm scipy pyparsing requests natsort pyyaml filelock libdeflate libcurl
conda install -c bioconda dendropy pysam
```

Besides, you need software packages including [mafft](http://mafft.cbrc.jp/alignment/software/) (multiple sequence aligner), [iqtree](http://www.iqtree.org/) (phylogenomic inference), [minimap2](https://github.com/lh3/minimap2) (long and short read mappers), and [samtools](http://www.htslib.org/download/) which could be installed using conda.
For this version, the `--read_type` argument accepts any minimap2 options string that defines how reads are aligned to the reference. For example, it could be `-ax sr`, `-ax map-hifi` or `-ax map-ont`. You can also pass `--threads 40` to be used with minimap2.
```
conda install -c bioconda mafft iqtree minimap2 samtools
```

For the coalescent species tree (step 4), [ASTER](https://github.com/chaoszhang/ASTER) (the C++ implementation of ASTRAL-III) is required. [ClipKIT](https://github.com/JLSteenwyk/ClipKIT) is optional and used only when `--trim` is passed.
```
conda install -c bioconda aster clipkit
```

Then, you can install the read2tree package after downlaoding the package from this GitHub repo using

```
git clone https://github.com/DessimozLab/read2tree.git
cd read2tree
python setup.py install
```




## Run

To run read2tree two things are required as input:
1) The DNA sequencing reads as FASTQ file(s).
2) A set of reference orthologous groups, i.e. marker genes. 
In our wiki [page](https://github.com/DessimozLab/read2tree/wiki/obtaining-marker-genes), you may find information on how to obtain the marker genes using [OMA browser](https://omabrowser.org/oma/export_markers). You can set the value of `Maximum nr of markers` as 200 or 400. Once you downloaded the tgz file, run this

```
tar xvzf  marker_genes_*.tgz 
ls marker_genes/*.fna | wc -l
cat marker_genes/*.fna > dna_ref.fa
``` 



### output 

The output of Read2Tree is the concatenated alignments as a fasta file where each record corresponds to one species. We also provide the option `--tree` for inferring the species tree using IQTREE as default (concatenation/supermatrix approach).

For a coalescent-based species tree that accounts for incomplete lineage sorting and differing gene tree histories, run the optional **step 4** after step 3 (see below).  


### Single species mode
```
read2tree --tree --standalone_path marker_genes/ --reads read_1.fastq read_2.fastq  --output_path output --dna_reference  dna_ref.fa 
```

### Multiple species mode

#### step1
```
read2tree  --step 1marker  --standalone_path marker_genes  --dna_reference dna_ref.fa --output_path output  --debug 
```

#### step2
The following could be run in parallel. 
```
read2tree --step 2map --standalone_path marker_genes  --dna_reference dna_ref.fa --reads species1_R1.fastq species2_R2.fastq  --output_path output --debug
read2tree --step 2map --standalone_path marker_genes  --dna_reference dna_ref.fa --reads species2_R1.fastq species2_R2.fastq  --output_path output --debug
read2tree --step 2map --standalone_path marker_genes  --dna_reference dna_ref.fa --reads species3_R1.fastq species3_R2.fastq  --output_path output  --debug
```

#### step3
```
read2tree  --step 3combine --standalone_path marker_genes  --dna_reference dna_ref.fa  --output_path output  --tree --debug
```

### Metagenomic mode

Pass `--meta` to enable metagenomic mode. In this mode, multiple consensus sequences per OG are kept (one per metagenomic species detected in the reads), instead of selecting a single best consensus.

```
read2tree --step 3combine --standalone_path marker_genes --dna_reference dna_ref.fa --output_path output --tree --meta --meta_min_markers 50 --meta_marker_fraction 0.5
```

Tunable filters (only active with `--meta`):
- `--meta_min_markers` (int, default `0`): minimum number of marker genes a metagenomic species must have to be retained.
- `--meta_marker_fraction` (float in 0-1, default `0.0`): minimum fraction of total marker genes a metagenomic species must have alignments in.

**Note on false positives:** metagenomic mode is permissive by design and may include species that share marker reads only by chance. We recommend tuning `--meta_min_markers` and `--meta_marker_fraction` to your dataset (e.g. `50` and `0.5` as a starting point for typical microbial communities) to reduce false positives.

#### step4 (optional: coalescent species tree)

Step 3 builds a supermatrix tree by concatenating all OG alignments. If you want a **coalescent-based species tree** instead, which better handles incomplete lineage sorting and the different evolutionary histories of individual genes, run step 4 after step 3:

```
read2tree --step 4astral --standalone_path marker_genes --dna_reference dna_ref.fa --output_path output --threads 24
```

Step 4 does the following automatically:
1. Filters the per-OG alignments from step 3 by taxon occupancy (`--min_samples`, default 10) and gap fraction (`--max_gap`, default 0.80).
2. Runs IQ-TREE on each passing alignment in parallel (`-m LG+F+G`, `-alrt 1000`, `--abayes`, `-fast`) to infer individual gene trees.
3. Collects all gene trees and passes them to [ASTER](https://github.com/chaoszhang/ASTER) to produce the final coalescent species tree.

**ASTER binary.** The ASTER suite provides several binaries, all installed by `conda install aster`. By default, step 4 auto-detects the first available in your PATH (`astral3`,`astral-pro3`,`astral-pro2`,`wastral`,`astral4`). Use `--astral_binary` to opt into a specific estimator:

| Binary | When to use |
|---|---|
| `astral3` | Default. Standard ASTRAL-III for single-copy orthologs. |
| `astral-pro3` | Gene trees with paralogs (multi-copy gene families). |
| `wastral` | Recommended for noisy gene trees. Weights each quartet by gene tree branch support, so poorly supported splits contribute less to the species tree. IQ-TREE's `--abayes` supports (already included in step 4) provide the weighting signal. |
| `astral4` | Large datasets with substantial missing taxa, or when substitution-rate branch lengths on internal nodes are needed for downstream rate analyses. |

```
# weighted ASTRAL - better accuracy when gene tree support is variable
read2tree --step 4astral --standalone_path marker_genes --dna_reference dna_ref.fa --output_path output --threads 24 --astral_binary wastral

# ASTRAL-IV - better robustness under missing data
read2tree --step 4astral --standalone_path marker_genes --dna_reference dna_ref.fa --output_path output --threads 24 --astral_binary astral4
```

**Choosing `--min_samples` for large datasets.** The default of 10 is intentionally permissive so the tool works out of the box for small test datasets. For studies with many samples, the occupancy threshold has a large effect on how many OGs survive filtering and on the quality of the resulting gene trees. A useful empirical guideline is to require at least 30–40% taxon occupancy. For example, `--min_samples 100` for a dataset of ~300 samples. In practice, applying a meaningful occupancy threshold together with `--max_gap 0.80` can reduce the number of OGs from tens of thousands to a few hundred; this is expected and desirable, as the surviving alignments are well-sampled across the tree and produce far more reliable gene trees for ASTRAL than a large set of sparse, gap-heavy alignments would. If the filtered set is very small (fewer than ~50 OGs), consider relaxing `--min_samples` slightly rather than `--max_gap`, since taxon occupancy drives gene tree resolution more than column-level gap content.

Optionally, pass `--trim` to run [ClipKIT](https://github.com/JLSteenwyk/ClipKIT) column-trimming on each alignment before gene tree inference:
```
read2tree --step 4astral --standalone_path marker_genes --dna_reference dna_ref.fa --output_path output --threads 24 --trim
```

**IQ-TREE and ASTER options.** Step 4 uses sensible defaults for per-gene tree inference (`-m LG+F+G -alrt 1000 --abayes -fast`). For cases where these need to be adjusted, three pass-through arguments are available:

| Argument | Default | Purpose |
|---|---|---|
| `--iqtree_model` | `LG+F+G` | Substitution model for per-gene IQ-TREE runs (e.g. `WAG+G`, `LG+G`, `TEST` for ModelFinder) |
| `--iqtree_args` | none | Extra flags appended verbatim to every per-gene IQ-TREE call (e.g. `"-B 1000"`) |
| `--astral_args` | none | Extra flags appended verbatim to the ASTER call; see the ASTER documentation for available options |
| `--no_fast` | off | Disable the `-fast` flag to run a full ML tree search per gene. Required when using bootstrap via `--iqtree_args` (e.g. `--iqtree_args "-B 1000"`), as `-fast` and bootstrap are incompatible in IQ-TREE |
| `--dna` | off | Use DNA alignments from `06_align_merge_dna` instead of amino acid alignments. Recommended for closely related species. Default model switches to `GTR+G` |

```
# Use WAG+G model instead of LG+F+G
read2tree --step 4astral --standalone_path marker_genes --dna_reference dna_ref.fa --output_path output --threads 24 --iqtree_model WAG+G

# Run ModelFinder per gene (much slower, but selects best-fit model for each OG)
read2tree --step 4astral --standalone_path marker_genes --dna_reference dna_ref.fa --output_path output --threads 24 --iqtree_model TEST

# Full ML search with ultrafast bootstrap (--no_fast required when using -B)
read2tree --step 4astral --standalone_path marker_genes --dna_reference dna_ref.fa --output_path output --threads 24 --no_fast --iqtree_args "-B 1000"

# DNA-based coalescent tree (recommended for closely related species)
read2tree --step 4astral --standalone_path marker_genes --dna_reference dna_ref.fa --output_path output --threads 24 --dna
```

Key output files written to `output/` (suffix is `aa` by default, `dna` when `--dna` is used):
- `07_astral_filtered_aa/` or `07_astral_filtered_dna/` - per-OG FASTA alignments that passed filtering
- `07_astral_trimmed_aa/` or `07_astral_trimmed_dna/` - ClipKIT-trimmed alignments (only when `--trim` is used)
- `08_gene_trees_aa/` or `08_gene_trees_dna/` - individual IQ-TREE gene tree files
- `gene_trees_merge_aa.nwk` or `gene_trees_merge_dna.nwk` - all gene trees concatenated into one file (input to ASTER)
- `astral_tree_merge_aa.nwk` or `astral_tree_merge_dna.nwk` - the final coalescent species tree in Newick format

### bootstraping

To have bootstrap values a metric for quality of internal nodes, you can run the following 
```
thread=20
iqtree -T ${thread} -s output/concat_*_aa.phy  -bb 1000 
```
The `.phy` file is either `concat_sample_aa.phy` or `concat_merge_aa.phy` corresponding to single- or multi-species mode. 

It is also possible to use [trimal](http://trimal.cgenomics.org/use_of_the_command_line_trimal_v1.2) for trimming msa `trimal -in <inputfile> -out <outputfile> -automated1` 

### DNA-mode tree inference
For closely related species, the user can infer tree using MSA of nucleotide sequences. 
```
thread=20
iqtree -T ${thread} -s output/concat_*_dna.phy 
```

## Test example

The goal of this test example is to infer species tree for Mus musculus using its sequencing reads. You can download the full read data from from [SRR5171076](https://www.ncbi.nlm.nih.gov/sra/?term=SRR5171076) using [sra-tools](https://anaconda.org/bioconda/sra-tools). Alternatively, a small read dataset is provided in the `tests` folder. For this example, we consider five species including Mnemiopsis leidyi, Xenopus laevis, Homo sapiens, Gorilla gorilla, and Rattus norvegicus as the reference. Using [OMA browser](https://omabrowser.org/oma/export_markers), we downloaded 20 marker genes of these five species as the reference orthologous groups, located in the folder `tests/mareker_genes`. 

```
cd tests
read2tree --tree --standalone_path marker_genes/ --reads sample_1.fastq sample_2.fastq  --output_path output --dna_reference  dna_ref.fa  
```


#### Run test example using docker
(to be updated  )
```
docker run --rm -i -v $PWD/tests:/input -v $PWD/tests/:/reads -v $PWD/outside_docker_out:/inside_docker_out -v $PWD/run:/run ${{ env.TEST_TAG }} --tree --standalone_path /input/marker_genes --dna_reference /input/dna_ref.fa --reads /reads/sample_1.fastq --output_path /inside_docker_out/output --debug --threads 1
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
It seems that minimap2 doesn't have the preset option  `-x sr` for short reads in older versions like 2.1. Then, you might get
```
Shell err: b"[E::main] unknown preset 'sr'\n"
```

Updating minimap2 to version 2.30 should fix the issue.


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

- version 2.0.1:
  - fixing a few bugs
- version 2.0.0:
  - pre-release: improve logging and make minimap2 as default aligner
- version 1.5:
  - using minimap2 as the read mapper in the minimap2 branch
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

* [David Dylus](https://github.com/dvdylus)
* [Adrian Altenhoff](http://people.inf.ethz.ch/adriaal).
* [Sina Majidian](https://sinamajidian.github.io/)


The authors would like to thank Alex Warwick for help how to initiate such a package.

## License
This project is licensed under the MIT License.

