#!/usr/bin/env python
'''
    This class will check the current status of the pipeline computation and determine how to start once several steps have been finished and
    will save a lot of time when several steps are finished. Also it will allow to submit the mapping (part that takes longest time) as job array.

    -- David Dylus, July--XXX 2017
'''

import mmap
import glob
import os
import time
import subprocess

OMA_STANDALONE_OUTPUT = 'Output'
OMA_MARKER_GENE_EXPORT = 'marker_genes'

class Progress(object):

    def __init__(self, args, species_name=None):
        self.args = args

        # holds the status of the computation
        self.ref_ogs_01 = False
        self.ref_dna_02 = False
        self.mapping_03 = False
        self.append_ogs_04 = False
        self.align_05 = False
        self.tree = False

        if self.args.reads:
            if len(self.args.reads) == 2:
                self._reads = self.args.reads
                self._species_name = self._reads[0].split("/")[-1].split(".")[0]
            else:
                self._reads = self.args.reads[0]
                self._species_name = self._reads.split("/")[-1].split(".")[0]

        if self.args.species_name:
            self._species_name = self.args.species_name

        if not self.args.reads and not self.args.species_name:
            self._species_name = 'merged'

        if species_name:
            self._species_name = species_name

        if self.args.remove_species_mapping:
            self.species_to_remove = self.args.remove_species_mapping.split(",")
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
            map_files = len(glob.glob(os.path.join(self._folder_mapping, "*cov.txt")))

        if map_files is self._num_species:
            return True
        else:
            return False

    def _set_num_species(self):
        """
        :return:
        """
        if self.status > 1:
            if self.species_to_remove:
                num_species = len(os.listdir(self._folder_ref_dna)) - len(self.species_to_remove)
            else:
                num_species = len(os.listdir(self._folder_ref_dna))
        else:
            num_species = 0
        return num_species

    def _get_status(self, species_name=None):
        if not species_name:
            species_name = self._species_name
        status = 0
        if os.path.exists(self.status_file):
            self._wait_for_status_file()
            f = open(self.status_file, 'r')
            for line in f:
                if '01_ref_ogs_aa: OK' in line:
                    status = 1
                    self.ref_ogs_01 = True
                elif '02_ref_dna: OK' in line:
                    status = 2
                    self.ref_dna_02 = True
                elif '03_mapping_' + species_name + ': OK' in line:
                    status = 3
                    self.mapping_03 = True
                elif '04_ogs_map_' + species_name + ': OK' in line:
                    status = 4
                    self.append_ogs_04 = True
                elif '05_align_' + species_name + ': OK' in line:
                    status = 5
                    self.align_05 = True
            f.close()
        return status

    def set_status(self, status, ref=None):
        status_text = None
        if not os.path.exists(self.status_file):
            to_append = self._write_header()
            with open(self.status_file, "w") as myfile:
                myfile.write(to_append)
        if status is 'ogs' and self.ref_ogs_01 is False:
            status_text = '01_ref_ogs_dna: OK\n' \
                          '01_ref_ogs_aa: OK\n'
            self.ref_ogs_01 = True
        elif status is 'ref' and self.ref_dna_02 is False:
            status_text = '02_ref_dna: OK\n'
            self.ref_dna_02 = True
        elif status is 'map' and self.mapping_03 is False:
            status_text = '03_mapping_'+self._species_name+': OK\n'
            self.mapping_03 = True
        # elif status is 'single_map' and ref is not None and self.status < 3:
        #     last_line = self._tail(self.status_file, 1)[-1].decode("utf-8")
        #     if 'OK' in last_line:
        #         status_text = '----- ' + self._species_name + ' -----\n'
        #         status_text += 'Mapping of ' + self._species_name + ' to ' + ref + '\n'
        #     else:
        #         status_text = 'Mapping of ' + self._species_name + ' to ' + ref + '\n'
        elif status is 're_ogs' and self.append_ogs_04 is False:
            status_text = '04_ogs_map_'+self._species_name+': OK\n'
            self.append_ogs_04 = True
        elif status is 'og_align'and self.align_05 is False:
            status_text = '05_align_'+self._species_name+': OK\n'
            self.align_05 = True
        if status_text:
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
        with open(self.status_file) as f:
            contents = f.read()
        if status_text not in contents.split("\n"):
            with open(self.status_file, "a") as myfile:
                myfile.write(status_text)

    def _find_last_completed_step(self):
        return subprocess.check_output(['tail', '-1', self.status_file])

    def _is_locked(self):
        """
        Taken from: https://www.calazan.com/how-to-check-if-a-file-is-locked-in-python/
        Checks if a file is locked by opening it in append mode.
        If no exception thrown, then the file is not locked.
        """
        locked = None
        file_object = None
        if os.path.exists(self.status_file):
            try:
                buffer_size = 8
                # Opening file in append mode and read the first 8 characters.
                file_object = open(self.status_file, 'a', buffer_size)
                if file_object:
                    locked = False
            except IOError:
                locked = True
            finally:
                if file_object:
                    file_object.close()

        return locked

    def _wait_for_status_file(self):
        """
        Taken from: https://www.calazan.com/how-to-check-if-a-file-is-locked-in-python/
        Checks if the files are ready.
        For a file to be ready it must exist and can be opened in append
        mode.
        """
        wait_time = 5
        # If the file exists but locked, wait wait_time seconds and check
        # again until it's no longer locked by another process.
        while self._is_locked():
            print("WAITING FOR STATUS FILE!")
            time.sleep(wait_time)




