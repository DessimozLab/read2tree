#!/usr/bin/env python
'''
    This class will check the current status of the pipeline computation and determine how to start once several steps have been finished and
    will save a lot of time when several steps are finished. Also it will allow to submit the mapping (part that takes longest time) as job array.

    -- David Dylus, July--XXX 2017
'''

import mmap
import glob
import os

OMA_STANDALONE_OUTPUT = 'Output'
OMA_MARKER_GENE_EXPORT = 'marker_genes'

class Progress(object):

    def __init__(self, args):
        self.args = args

        if len(self.args.reads) == 2:
            self._reads = self.args.reads
            self._species_name = self._reads[0].split("/")[-1].split(".")[0]
        else:
            self._reads = self.args.reads[0]
            self._species_name = self._reads.split("/")[-1].split(".")[0]

        if self.args.remove_species:
            self.species_to_remove = self.args.remove_species.split(",")
        else:
            self.species_to_remove = []

        self._folder_ref_ogs_aa = os.path.join(self.args.output_path, "01_ref_ogs_aa")
        self._folder_ref_ogs_dna = os.path.join(self.args.output_path, "01_ref_ogs_dna")
        self._folder_ref_dna = os.path.join(self.args.output_path, '02_ref_dna')
        self._folder_mapping = os.path.join(self.args.output_path, "03_mapping_" + self._species_name)
        self._folder_ogs_map = os.path.join(self.args.output_path, "04_ogs_map" + self._species_name)

        self.status_file = os.path.join(self.args.output_path, 'status.txt')
        self.status = self._get_status()
        # self.oma_output_path = self.args.oma_output_path
        # self._num_ogs = oma_output.num_selected_ogs
        self._num_species = self._set_num_species()

        # self.status = self._determine_progress()

    def check_mapping(self):
        map_files = 0
        if os.path.exists(self._folder_mapping):
            for file in glob.glob(os.path.join(self._folder_mapping, "*cov.fa")):
                if os.path.getsize(file) > 0:
                    map_files += 1

        if map_files == self._num_species:
            return True
        else:
            return False

    def _set_num_species(self):
        """
        :return:
        """
        num_species = 0
        if self.status > 1:
            if self.args.remove_species:
                num_species = len(os.listdir(self._folder_ref_dna)) - len(self.args.remove_species)
            else:
                num_species = len(os.listdir(self._folder_ref_dna))
        else:
            num_species = 0
        return num_species

    def _get_status(self):
        status = 0
        if os.path.exists(self.status_file):
            self._find_last_completed_step()
            f = open(self.status_file, 'r', os.O_NONBLOCK)
            for line in f:
                # last_line = self._tail(self.status_file, 1)[-1].decode("utf-8")
                if '01_ref_ogs_aa: OK' in line:
                    status = 1
                elif '02_ref_dna: OK' in line:
                    status = 2
                elif '03_mapping_'+self._species_name+': OK' in line:
                    status = 3
                elif '04_ogs_map_'+self._species_name+': OK' in line:
                    status = 4
                elif '05_align_'+self._species_name+': OK' in line:
                    status = 5
            f.close()
        return status

    def set_status(self, status, ref=None):
        if not os.path.exists(self.status_file):
            to_append = self._write_header()
            with open(self.status_file, "w") as myfile:
                myfile.write(to_append)
        if status is 'ogs' and self.status < 1:
            status_text = '01_ref_ogs_dna: OK\n' \
                          '01_ref_ogs_aa: OK\n'
        elif status is 'ref' and self.status < 2:
            status_text = '02_ref_dna: OK\n'
        elif status is 'map' and self.status < 3:
            status_text = '03_mapping_'+self._species_name+': OK\n'
        # elif status is 'single_map' and ref is not None and self.status < 3:
        #     last_line = self._tail(self.status_file, 1)[-1].decode("utf-8")
        #     if 'OK' in last_line:
        #         status_text = '----- ' + self._species_name + ' -----\n'
        #         status_text += 'Mapping of ' + self._species_name + ' to ' + ref + '\n'
        #     else:
        #         status_text = 'Mapping of ' + self._species_name + ' to ' + ref + '\n'
        elif status is 're_ogs':
            status_text = '04_ogs_map_'+self._species_name+': OK\n'
        elif status is 'og_align':
            status_text = '05_align_'+self._species_name+': OK\n'
        self._append_status(status_text)

    def _tail(self, filename, n):
        """Returns last n lines from the filename.
        No exception handling
        https://stackoverflow.com/questions/136168/get-last-n-lines-of-a-file-with-python-similar-to-tail"""
        size = os.path.getsize(filename)
        with open(filename, "r") as f:
            # for Windows the mmap parameters are different
            fm = mmap.mmap(f.fileno(), 0, mmap.MAP_SHARED, mmap.PROT_READ)
            try:
                for i in range(size - 1, -1, -1):
                    if fm[i] == '\n':
                        n -= 1
                        if n == -1:
                            break
                return fm[i + 1 if i else 0:].splitlines()
            finally:
                fm.close()

    def _write_header(self):
        header = '--- Computation Status ---\n'
        return header

    def _append_status(self, status_text):
        to_append = status_text
        with open(self.status_file, "a") as myfile:
            myfile.write(to_append)

    def _find_last_completed_step(self):
        f = open(self.status_file, 'r')
        idx_last_completed_step = 0
        all_lines = []
        for i, line in enumerate(f):
            if 'OK' in line:
                idx_last_completed_step = i
            all_lines.append(line)
        open(self.status_file, 'w').writelines(all_lines[0:idx_last_completed_step+1])



