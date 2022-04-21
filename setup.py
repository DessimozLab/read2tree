from setuptools import setup, find_packages

name = 'read2tree'

__version__ = None
with open('{:s}/__init__.py'.format(name), 'rt') as fp:
    for line in fp:
        if line.startswith('__version__'):
            exec(line.rstrip())

# conda install -c conda-forge biopython numpy Cython ete3 lxml tqdm scipy pyparsing requests natsort pyyaml
# conda install -c bioconda dendropy 

requirements = ['pysam', 'pyham', 'filelock']

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name=name,
    version=__version__,
    author='David Dylus and Fritz Sedlaczek',
    author_email='daviddylus@gmail.com',
    description='read2tree allows to build high quality phylogenetic trees '
                'using reads and a reference set of orthologous groups '
                '(DNA + Protein).',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dessimozlab/read2tree",
    packages=find_packages(".", exclude=["archive"]),
    include_package_data=True,
    package_data={
          'read2tree': ['logging/log.yaml']
      },
    install_requires=requirements,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Environment :: Console",
        "License :: OSI Approved :: MIT License",
    ],
    scripts=['bin/read2tree'],
    python_requires=">=3.5",
)
