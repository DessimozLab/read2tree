import glob
import os
import re

from tqdm import tqdm
from Bio import SeqIO, Seq, SeqRecord
from pathlib import Path


OMA_STANDALONE_OUTPUT = 'Output'
OMA_MARKER_GENE_EXPORT = 'marker_genes'


class OMAOutputParser(object):

    def __init__(self, args):
        '''
        Initialise the OMAOutputParser
        '''
        self.args = args
        self.mode = ''
        self.oma_output_path = self._check_oma_output_path()
        self.num_selected_ogs = 0
        self.num_species = 0

        self.min_species = self._estimate_best_number_species()

        if self.args.ignore_species:
            self.ignore_species = self.args.ignore_species.split(",")
        else:
            self.ignore_species = []
        self.ogs = self._load_ogs_from_path()

    def _check_oma_output_path(self):
        """
        This method analysis how the standalone_path argument is used:
          1. OmaStandalone run, and standalone_path points to base folder
          2. OmaStandalone run, but standalone_path points to output folder
          3. Marker gene export, standalone_path points to the marker_genes folder
          4. Marker genes export, standalone_path points to a folder with marker genes, but named differently
             (maybe even different source)
          5. Marker genes export, standalone_path points to parent folder containing marker_genes subfolder

        All possibilities are analysed in this order and first one wins. The method sets then the `mode`,
        oma_output_path, og_fasta_path

        :return: None
        """
        standalone_path = Path(self.args.standalone_path)
        self.oma_output_path = None
        if (standalone_path / OMA_STANDALONE_OUTPUT / "OrthologousGroupsFasta").is_dir():
            self.mode = "standalone"
            self.oma_output_path = standalone_path / OMA_STANDALONE_OUTPUT
            self.og_fasta_path = standalone_path / "Output" / "OrthologousGroupsFasta"
        elif (standalone_path / "OrthologousGroupsFasta").is_dir():
            self.mode = "standalone"
            self.oma_output_path = standalone_path
            self.og_fasta_path = standalone_path / "OrthologousGroupsFasta"
        elif standalone_path.parts[-1] == OMA_MARKER_GENE_EXPORT:
            self.og_fasta_path = standalone_path
            self.mode = "marker_genes"
        else:
            self.mode = 'marker_genes'
            nr_fasta_files = len([name for name in os.listdir(standalone_path)
                                  if name.endswith(".fa") or name.endswith(".fasta")])
            if nr_fasta_files > 1:
                self.og_fasta_path = standalone_path
            elif (standalone_path / OMA_MARKER_GENE_EXPORT).is_dir():
                self.og_fasta_path = standalone_path / OMA_MARKER_GENE_EXPORT
            else:
                raise Exception("argument standalone_path seems not to point to a valid directory: {}"
                                .format(standalone_path))

    def _load_ogs_from_path(self):
        return self._filter_ogs_min_species()

    def _get_species_id(self, record):
        m = re.search(r"\[(?P<species>[^]]+)]", record.description)
        if m is not None:
            return m.group('species')
        else:
            return record.id[0:5]

    def _filter_ogs_min_species(self):
        """

        :return:
        """
        names_og = {}
        unique_species = set([])
        included_species = set([])
        print('--- Load OGs with min {} species from oma {} - mode = {} ---'
              .format(self.min_species, self.og_fasta_path, self.mode))
        files = (glob.glob(os.path.join(self.og_fasta_path, "*.fa")) or
                 glob.glob(os.path.join(self.og_fasta_path, "*.fasta")))
        for file in tqdm(files, desc='Loading files for pre-filter',
                         unit=' OGs'):
            name = os.path.splitext(os.path.basename(file))[0]
            name = name.replace("OMAGroup_", "OG")
            records = list(SeqIO.parse(file, 'fasta'))
            new_records = []
            seen_species = set([])
            for record in records:
                species = self._get_species_id(record)
                if species not in self.ignore_species:
                    if species in seen_species:
                        raise Exception("Invalid marker group: {} contains species '{}' more than once"
                                        .format(file, species))
                    seen_species.add(species)
                    new_id = record.id + "_" + name
                    record.id = new_id
                    new_records.append(record)

            unique_species.update(seen_species)
            if len(new_records) >= self.min_species:
                names_og[name] = new_records
                self.num_selected_ogs += 1
                included_species.update(seen_species)

        self.num_species = len(included_species)
        ignored = unique_species - included_species
        if len(ignored) > 0:
            print("{} species in marker genes were excluded as they never were present in sufficiently "
                  "complete marker genes".format(len(ignored)))
            for x in ignored:
                print(" - {}".format(x))
        return names_og

    def _estimate_best_number_species(self):
        """
        Estimate min number of species such that around 1000 OGs are selected
        :return:
        """
        min_species = 0
        if self.args.min_species != None:
            min_species = self.args.min_species

        return min_species
