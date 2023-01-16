#!/usr/bin/env python
'''
    This class will check the current status of the pipeline computation and determine how to start once several steps have been finished and
    will save a lot of time when several steps are finished. Also it will allow to submit the mapping (part that takes longest time) as job array.

    -- David Dylus, July--XXX 2017
'''

import glob
import os
import logging


class Progress(object):

    def __init__(self, args, species_name=""):

        self.args = args

        self._reads = self.args.reads
        self._species_name = self.args.species_name

        self.logger = logging.getLogger(__name__)

        if species_name:
            self._species_name = species_name

        if self.args.remove_species_mapping:
            self.species_to_remove = self.args.remove_species_mapping.split(",")
        else:
            self.species_to_remove = []

        self._folder_ref_ogs_aa = os.path.join(self.args.output_path,
                                               '01_ref_ogs_aa')
        self._folder_ref_ogs_dna = os.path.join(self.args.output_path,
                                                '01_ref_ogs_dna')
        self._folder_ref_dna = os.path.join(self.args.output_path,
                                            '02_ref_dna')
        self._folder_mapping = os.path.join(self.args.output_path,
                                            '04_mapping_' + self._species_name)
        self._folder_align_aa = os.path.join(self.args.output_path,
                                            '03_align_aa')
        self._folder_align_dna = os.path.join(self.args.output_path,
                                            '03_align_dna')
        self._folder_append_og_aa = os.path.join(self.args.output_path,
                                             '05_ogs_map_' + self._species_name + '_aa')
        self._folder_append_og_dna = os.path.join(self.args.output_path,
                                              '05_ogs_map_' + self._species_name + '_dna')
        self._folder_align_append_aa = os.path.join(self.args.output_path,
                                             '06_align_' + self._species_name + '_aa')
        self._folder_align_append_dna = os.path.join(self.args.output_path,
                                              '06_align_' + self._species_name + '_dna')


        # holds the status of the computation
        self._num_species = self._get_number_of_references()
        self.ref_ogs_01 = self._get_og_set_status()
        self.ref_dna_02 = self._get_reference_status()
        self.ref_align_03 = self._get_alignment_status()
        self.mapping_04 = self._get_mapping_status()  # add here True for species removal test
        self.append_ogs_05 = self._get_append_og_set_status()
        self.align_06 = self._get_append_alignment_status()
        self.num_completed_mappings = len(self._get_finished_mapping_folders(self.args.output_path))
        self.tree = False

        # self.status_file = os.path.join(self.args.output_path, 'status.txt')

    def update_status(self):
        self._num_species = self._get_number_of_references()
        self.ref_ogs_01 = self._get_og_set_status()
        self.ref_dna_02 = self._get_reference_status()
        self.mapping_03 = self._get_mapping_status()
        self.append_ogs_04 = self._get_append_og_set_status()
        self.align_05 = self._get_alignment_status()
        self.tree = False

    def _extract_line_from_log(self, word, logfile):
        '''
        Extract relevant line from log
        Adapted from https://stackoverflow.com/questions/43177256/python-extract-single-line-from-file
        :param word: string to search
        :param logfile: log file (typically mplog.log)
        :return:
        '''

        try:
            with open(logfile, "r") as file:
                bestline = [line.split() for line in file if word in line]
                if bestline:
                    return bestline[-1]
                return None
        except FileNotFoundError:
            print('File {} not accessible'.format(logfile))
            return None

    def _get_number_of_OGs(self):
        '''
        Example log line:
        2018-11-23 12:13:53,688 - read2tree.OGSet - INFO - ass: Gathering of DNA seq for 5 OGs took 0.004261970520019531.
        :return: Number of OGs
        '''
        log_list = self._extract_line_from_log('Gathering', 'mplog.log')
        if log_list:
            return int(log_list[13])
        else:
            return 0

    def _get_number_of_appeneded_seq_to_OGs(self):
        '''
        Example log line:
        2018-11-23 12:13:53,688 - read2tree.OGSet - INFO - ass: Gathering of DNA seq for 5 OGs took 0.004261970520019531.
        :return: Number of OGs
        '''
        log_list = self._extract_line_from_log('Appending', 'mplog.log')
        if log_list:
            return int(log_list[9])
        else:
            return 0

    def _get_number_of_alignments(self):
        '''
        Example log line:
        2018-11-23 12:13:53,688 - read2tree.OGSet - INFO - ass: Gathering of DNA seq for 5 OGs took 0.004261970520019531.
        :return: Number of OGs
        '''
        log_list = self._extract_line_from_log('Alignment of', 'mplog.log')
        if log_list:
            return int(log_list[10])
        else:
            return 0

    def _get_number_of_references(self):
        '''
        Example log line:
        2018-11-23 12:13:53,691 - read2tree.ReferenceSet - INFO - ass: Extracted 6 reference species form 5 ogs took 0.0008709430694580078
        :return: Number of reference species
        '''
        log_list = self._extract_line_from_log('ReferenceSet', 'mplog.log')
        if log_list:
            return int(log_list[9])
        else:
            return 0

    def _count_files(self, path, ext):
        '''
        https://stackoverflow.com/questions/2632205/how-to-count-the-number-of-files-in-a-directory-using-python/16865840
        '''
        if len(os.listdir(path)) != 0:
            return len([f for f in glob.glob(os.path.join(path, ext)) if os.path.getsize(f) > 0])
        else:
            return 0

    def _get_og_set_status(self):
        '''
        Get OG status
        :return:
        '''
        num_ogs_expected = self._get_number_of_OGs()
        if os.path.exists(self._folder_ref_ogs_aa) and os.path.exists(self._folder_ref_ogs_dna):
            num_ogs_aa = self._count_files(self._folder_ref_ogs_aa, '*fa')
            num_ogs_dna = self._count_files(self._folder_ref_ogs_dna, '*fa')
            if (num_ogs_expected-num_ogs_aa) == 0 and (num_ogs_expected-num_ogs_dna) == 0:
                return True
            else:
                return False
        else:
            return False

    def _get_append_og_set_status(self):
        '''
        Get OG status
        :return:
        '''
        num_ogs_expected = self._get_number_of_appeneded_seq_to_OGs()
        if os.path.exists(self._folder_append_og_aa) and os.path.exists(self._folder_append_og_dna):
            num_ogs_aa = self._count_files(self._folder_append_og_aa, '*fa')
            num_ogs_dna = self._count_files(self._folder_append_og_dna, '*fa')
            if (num_ogs_expected-num_ogs_aa) <= 0 and (num_ogs_expected-num_ogs_dna) <= 0:
                return True
            else:
                return False
        else:
            return False

    def _get_reference_status(self):
        '''
        Get Reference status
        :return:
        '''
        num_ref_expected = self._get_number_of_references()
        if os.path.exists(self._folder_ref_dna):
            num_references = self._count_files(self._folder_ref_dna, '*fa')
            if (num_ref_expected-num_references) == 0:
                return True
            else:
                return False
        else:
            return False

    def _get_alignment_status(self):
        '''
        Get OG status
        :return:
        '''
        num_aligns_expected = self._get_number_of_alignments()
        if os.path.exists(self._folder_align_aa) and os.path.exists(self._folder_align_dna):
            num_align_aa = self._count_files(self._folder_align_aa, '*phy')
            num_align_dna = self._count_files(self._folder_align_dna, '*phy')
            if (num_aligns_expected-num_align_aa) == 0 and (num_aligns_expected-num_align_dna) == 0:
                return True
            else:
                return False
        else:
            return False

    def _get_append_alignment_status(self):
        '''
        Get OG status
        :return:
        '''
        num_aligns_expected = self._get_number_of_alignments()
        if os.path.exists(self._folder_align_append_aa) and os.path.exists(self._folder_align_append_dna):
            num_align_aa = self._count_files(self._folder_align_append_aa, '*phy')
            num_align_dna = self._count_files(self._folder_align_append_dna, '*phy')
            if (num_aligns_expected-num_align_aa) == 0 and (num_aligns_expected-num_align_dna) == 0:
                return True
            else:
                return False
        else:
            return False

    def _get_finished_mapping_folders(self, path):
        mapping_folders_finished = []
        num_expected_mappings = self._get_number_of_references()
        mapping_folders = [x for x in os.listdir(path) if '04' in x]
        for folder in mapping_folders:
            # NOTE: we are calculating the number of completed mappings as the number of existing cov files,
            # because these are written even if the mapping step did not find any reads to map to a particular reference
            computed_cov = [f for f in
                            glob.glob(os.path.join(self.args.output_path,
                                                   folder+'/*cov.txt'))]
            # it is finished if the number of generated coverage files is the same as the number of references
            if self.args.merge_all_mappings:
                if num_expected_mappings >= len(computed_cov) and len(computed_cov) >= 1:
                    mapping_folders_finished.append(folder)
            elif self._species_name in folder:
                if (num_expected_mappings - len(computed_cov)) == 0:
                    mapping_folders_finished.append(folder)
            # else:
        #         self.logger.debug(
        #             '{}: {} NOT completed'.format(self._species_name, path))
        # self.logger.info('{}: From {} mapping {} are completed'. format(self._species_name, len(mapping_folders_finished), len(mapping_folders)))
        return mapping_folders_finished

    def _get_mapping_status(self):
        mapping_folders = self._get_finished_mapping_folders(self.args.output_path)
        if mapping_folders:
            if len(mapping_folders) > 0:
                self.num_completed_mappings = len(mapping_folders)
                # self.logger.info('{}: Mapping completed!'.format(self._species_name))
                return True
            else:
                self.num_completed_mappings = 0
                # self.logger.info('{}: Mapping not completed!'.format(self._species_name))
                return False
        else:
            return False
