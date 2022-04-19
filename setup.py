from setuptools import setup, find_packages

name = 'read2tree'
__version__ = None
with open('{:s}/__init__.py'.format(name), 'rt') as fp:
    for line in fp:
        if line.startswith('__version__'):
            exec(line.rstrip())

requirements = ['biopython', 'numpy', 'Cython', 'ete3', 'dendropy', 'lxml',
                'tqdm', 'scipy', 'pysam', 'pyham', 'pyparsing', 'requests',
                'filelock', 'natsort', 'pyyaml']



setup(
    name=name,
    version=__version__,
    author='David Dylus and Fritz Sedlaczek',
    author_email='daviddylus@gmail.com',
    description='read2tree allows to build high quality phylogenetic trees '
                'using reads and a reference set of orthologous groups '
                '(DNA + Protein).',
    packages=find_packages(".", exclude=["archive"]),
    include_package_data=True,
    package_data={
          'read2tree': ['logging/log.yaml']
      },
    install_requires=requirements,
    license='Proprietary',
    scripts=['bin/read2tree'])
