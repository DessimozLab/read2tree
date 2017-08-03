from setuptools import setup, find_packages

name = 'pore2tree'
__version__ = None
with open('{:s}/__init__.py'.format(name), 'rt') as fp:
    for line in fp:
        if line.startswith('__version__'):
            exec(line.rstrip())

requirements = ['biopython', 'ete3', 'dendropy', 'lxml', 'pandas', 'tqdm',
                'scipy', 'zoo']

setup(
    name=name,
    version=__version__,
    author='David Dylus and Fritz Sedlaczek',
    author_email='david.dylus@unil.ch',
    description='pore2tree allows to build high quality phylogenetic trees '
                'using reads and a reference set of orthologous groups (DNA + Protein).',
    packages=find_packages(),
    install_requires=requirements,
    license='Proprietary',
    scripts=['bin/pore2tree'])
