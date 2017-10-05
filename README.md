# read2tree 

Given a set of reference OGs and sequencing reads this tool allows to automatically build phylogenetic trees of similar quality to a full blown genomics pipeline. 

## Getting Started

read2tree was build and tested with python 3.5.1. To set up read2tree on your local machine.
```
git clone https://github.com/dvdylus/read2tree.git
python setup.py install
```
### Prerequisites

read2tree integrates multiple software tools and allows to infer a phylogenetic tree skipping several steps of a usual pipeline.

* [mafft](http://mafft.cbrc.jp/alignment/software/) - Mutliple sequence alignment software
* [fasttree](http://www.microbesonline.org/fasttree/) - Dependency Management
* [ngmlr](https://github.com/philres/ngmlr) - Long read mapper for nanopore or PacBio read data
* [pyopa](...) - Implementation of Smith Waterman alignment algorithm in python
* [pyoma](...) - Library for retrieval of nucleotide sequences from oma run
* [pyham](...) - Library to work with HOGs
* [samtools](http://www.htslib.org/download/) - Set of programs to interact with high-throughput sequencing data

### Installing

For mafft, fasttree, ngmlr and samtools please follow the instructions provided by the individual packages. 

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Once successfully installed

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [pyCharm](https://www.jetbrains.com/pycharm) - Python IDE

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Alex Warwick for help how to initiate such a package.
__