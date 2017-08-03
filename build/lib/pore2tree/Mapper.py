#!/usr/bin/env python
'''
    This file contains definitions of a class that allows to map reads to a reference file!
    Importantly the mapper function relies heavy on the following:
    ngm v0.5.4
    ngmlr v0.2.6
    samtools
    bcftools
    vcfutils.pl
    
    
    -- David Dylus, July--XXX 2017
'''

import pysam
import pyopa
import os
import tempfile
import subprocess
from tqdm import tqdm
from Bio import SeqIO, SeqRecord
from Bio.SeqIO.FastaIO import FastaWriter
from pore2tree.OGSet import OG
from pore2tree.ReferenceSet import Reference
from pore2tree.wrappers.read_mappers import NGM
from pore2tree.wrappers.read_mappers import NGMLR
from tables import *


class Mapper(object):
    """
    Structure for reference
    """

    def __init__(self, args, ref_set=None, og_set=None):
        print('--- Mapping of reads to reference sequences ---')
        self.args = args
        self._reads = args.reads

        # load pyopa related stuff
        self.defaults = pyopa.load_default_environments()
        self.envs = self.defaults['environments']
        self.env = self.envs[515]

        if ref_set is not None:
            self.mapped_records = self._map_reads_to_reference(ref_set)

        if self.mapped_records and og_set is not None:
            self.og_records = self._sort_by_og(og_set)

    def _map_reads_to_reference(self, reference):
        """
        Map all the reads to reference
        :param reference: 
        :return: dictionary with key og_name and value sequences mapped to each species
        """
        mapped_reads_species = {}
        output_folder = os.path.join(self.args.output_path, "03_mapping")
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        for species, value in tqdm(reference.items(), desc='Mapping reads to species', unit=' species'):
            # write reference into temporary file
            #ref_file_handle = tempfile.NamedTemporaryFile(mode='wt', delete=True)
            ref_file_handle = os.path.join(output_folder, species+'.fa')
            SeqIO.write(value.dna, ref_file_handle, 'fasta')
            # call the WRAPPER here
            if len(self._reads) == 2:
                ngm_wrapper = NGM(ref_file_handle, self._reads)
                ngm = ngm_wrapper()
            else:
                ngm_wrapper = NGMLR(ref_file_handle, self._reads)
                ngm = ngm_wrapper()
                sam_file = ngm['file']

            mapped_reads = list(SeqIO.parse(self._post_process_read_mapping(ref_file_handle, sam_file), 'fasta'))
            mapped_reads_species[species] = Reference()
            mapped_reads_species[species].dna = mapped_reads

        return mapped_reads_species

    def _post_process_read_mapping(self, ref_file, sam_file):
        """
        Function that will perform postprocessing of finished read mapping using the pysam functionality
        :param ngm: 
        :param reference_file_handle: 
        :return: 
        """
        outfile_name = ref_file.split(".")[0]+"_post"
        # with tempfile.NamedTemporaryFile(mode='wt') as file_handle:
        #     outfile_name = file_handle.name

        pysam.view("-bh", "-S", "-o", outfile_name + ".bam", sam_file, catch_stdout=False)
        pysam.sort("-o", outfile_name + "_sorted.bam", outfile_name + ".bam")
        pysam.index(outfile_name + "_sorted.bam")

        cmd = 'samtools mpileup -d 100000 -B -uf ' + ref_file + ' ' + outfile_name + '_sorted.bam | bcftools call -c | vcfutils.pl vcf2fq -d 1'

        with open(outfile_name + '_consensus_call.fq', "wb") as out:
            out.write(self._output_shell(cmd))

        fastq_records = list(SeqIO.parse(outfile_name + '_consensus_call.fq', 'fastq'))

        SeqIO.write(fastq_records, outfile_name + '_consensus_call.fa', 'fasta')


        return outfile_name+'_consensus_call.fa'

    def _output_shell(self, line):
        """
        Save output of shell line that has pipes
        taken from: https://stackoverflow.com/questions/7389662/link-several-popen-commands-with-pipes
        :param line: 
        :return: 
        """
        try:
            shell_command = subprocess.Popen(line, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        except OSError:
            return None
        except ValueError:
            return None

        (output, err) = shell_command.communicate()
        shell_command.wait()
        if shell_command.returncode != 0:
            print("Shell command failed to execute")
            return None

        return output

    def _sort_by_og(self, og_set):
        """
        Group the mapped sequences according to their OG name
        :return: dictionary with lists of records as values and og name as keys
        """
        og_records = {}
        for records in self.mapped_records.values():
            for record in records.dna:
                name = record.id.split("_")[-1]
                if name in og_records.keys():
                    og_records[name].dna.append(record)
                    aa = self._predict_best_protein_pyopa(record, og_set[name])
                    og_records[name].aa.append(aa)
                else:
                    og_records[name] = OG()
                    og_records[name].dna.append(record)
                    aa = self._predict_best_protein_pyopa(record, og_set[name])
                    og_records[name].aa.append(aa)

        return og_records

    def _predict_best_protein_position(self):
        raise NotImplementedError

    def _predict_best_protein_pyopa(self, record, og):
        """
        Given a list of sequences that are derived from mapped reads to multiple seq of a OG
        we find the best corresponding mapped seq by comparing it with a representative sequence of the original OG using
        pyopa local alignment and return the sequence with its highest score!
        :return: 
        """
        ref_og_seq = og.aa[0]
        s1 = pyopa.Sequence(str(ref_og_seq.seq))
        best_score = 0
        if 'X' not in record.seq and 'N' not in record.seq:
            frames = [record.seq[i:].translate(table='Standard', stop_symbol='X', to_stop=False, cds=False) for i
                      in range(3)]
            best_seq_idx = 0
            for i, seq in enumerate(frames):
                s2 = pyopa.Sequence(str(seq))
                # calculating local and global scores for the given sequences
                local_double = pyopa.align_double(s1, s2, self.env)
                # print('Local score: %f' % local_double[0])
                if local_double[0] > best_score:
                    best_score = local_double[0]
                    best_seq_idx = i
            best_translation = SeqRecord.SeqRecord(frames[best_seq_idx], id="UNKNOWN", description=record.description, name=record.name)
        return best_translation

    def write_by_og(self, output_folder):
        '''
        Write for each og all the mapped sequences into separate fasta files to a specified folder
        :param output_folder: folder where files should be stored
        '''
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        for key, value in tqdm(self.og_records.items(), desc="Writing DNA seq sorted by OG", unit=" OG"):
            handle = open(os.path.join(output_folder, 'mapped_'+key+'.fa'), "w")
            writer = FastaWriter(handle, wrap=None)
            writer.write_file(value)
            handle.close()
