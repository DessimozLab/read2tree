import glob
import os

from tqdm import tqdm
from Bio import SeqIO, Seq, SeqRecord
#from tables import *


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

        :return:
        """
        oma_output_path = ''
        if os.path.exists(os.path.join(self.args.standalone_path,
                                       "OrthologousGroupsFasta")):
            OMA_STANDALONE_OUTPUT = "."
        else:
            OMA_STANDALONE_OUTPUT = 'Output'

        if os.path.exists(os.path.join(self.args.standalone_path,
                                       OMA_STANDALONE_OUTPUT)):
            oma_output_path = os.path.join(self.args.standalone_path,
                                           OMA_STANDALONE_OUTPUT)
            if os.path.join(oma_output_path, "OrthologousGroups.orthoxml"):
                self.mode = 'standalone'
        # standalone path set to path up of marker_genes
        elif os.path.exists(os.path.join(self.args.standalone_path,
                                         OMA_MARKER_GENE_EXPORT)):
            oma_output_path = self.args.standalone_path
            self.mode = 'marker_genes'
        # standalone path set to marker_genes
        elif OMA_MARKER_GENE_EXPORT in self.args.standalone_path.split('/'):
            oma_output_path = os.path.join(self.args.standalone_path, "..")
            self.mode = 'marker_genes'
        return oma_output_path

    def _load_ogs_from_path(self):

        if self.mode == "marker_genes":
            ogs = self._filter_ogs_min_species_marker()
        else:
            ogs = self._filter_ogs_min_species()
        return ogs

    def _filter_ogs_min_species(self):
        """

        :return:
        """
        names_og = {}
        unique_species = []
        orthologous_groups_fasta = os.path.join(self.oma_output_path,
                                                "OrthologousGroupsFasta")
        print('--- Load OGs with min {} species from oma '
              'standalone! ---'.format(self.min_species))
        files = (glob.glob(os.path.join(orthologous_groups_fasta, "*.fa")) or
                 glob.glob(os.path.join(orthologous_groups_fasta, "*.fasta")))
        for file in tqdm(files, desc='Pre-filter files',
                         unit=' OGs'):
            name = file.split("/")[-1].split(".")[0]
            records = list(SeqIO.parse(file, 'fasta'))
            new_records = []
            for record in records:
                species = self._get_species_id(record)
                if species not in self.ignore_species:
                    new_records.append(record)
                    if species not in unique_species:
                        unique_species.append(species)
                        self.num_species += 1
                new_id = record.id + "_" + name
                record.id = new_id

            if len(new_records) >= self.min_species:
                names_og[name] = new_records
                self.num_selected_ogs += 1
        return names_og

    def _get_species_id(self, record):
        if '[' in record.description and ']' in record.description:
            return record.description[record.description.find(
                "[")+1:record.description.find("]")]
        else:
            return record.id[0:5]

    def _filter_ogs_min_species_marker(self):
        """

        :return:
        """
        names_og = {}
        unique_species = set([])
        included_species = set([])
        orthologous_groups_fasta = os.path.join(self.oma_output_path,
                                                "marker_genes")
        print('--- Load OGs with min {} species from oma '
              'marker gene export! ---'.format(self.min_species))
        files = (glob.glob(os.path.join(orthologous_groups_fasta, "*.fa")) or
                 glob.glob(os.path.join(orthologous_groups_fasta, "*.fasta")))
        for file in tqdm(files, desc='Loading files for pre-filter',
                         unit=' OGs'):
            name = file.split("/")[-1].split(".")[0].replace('OMAGroup_', 'OG')
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
