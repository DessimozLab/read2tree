#!/usr/bin/env python
'''
    CoalescentInference: filter per-OG alignments, infer per-gene trees with IQ-TREE,
    and run ASTER to produce a coalescent species tree (step 4astral).

    -- Step 4 of the read2tree pipeline.
'''
import os
import glob
import time
import logging
from multiprocessing import Pool
from Bio import SeqIO, AlignIO

from read2tree.wrappers.treebuilders.iqtree import Iqtree, get_gene_tree_options, get_gene_tree_dna_options
from read2tree.wrappers.treebuilders.base_treebuilder import DataType
from read2tree.wrappers.treebuilders.aster import Aster
from read2tree.wrappers.options import StringOption
from read2tree.wrappers import WrapperError

logger = logging.getLogger(__name__)


def _run_gene_tree(task):
    """
    Module-level worker for multiprocessing pool.
    Runs IQ-TREE on a single alignment file and writes the treefile.
    Returns the Newick tree string, or None on failure.
    """
    alignment_file, gene_trees_folder, iqtree_model, iqtree_extra, no_fast, seq_type = task
    og_name = os.path.basename(alignment_file).rsplit('.', 1)[0]
    treefile = os.path.join(gene_trees_folder, og_name + '.treefile')
    if os.path.exists(treefile) and os.path.getsize(treefile) > 0:
        with open(treefile, 'r') as fh:
            return fh.read().strip()
    try:
        if seq_type == 'dna':
            iqtree_wrapper = Iqtree(alignment_file, datatype=DataType.DNA)
            iqtree_wrapper.options = get_gene_tree_dna_options()
        else:
            iqtree_wrapper = Iqtree(alignment_file, datatype=DataType.PROTEIN)
            iqtree_wrapper.options = get_gene_tree_options()
        if iqtree_model:
            iqtree_wrapper.options.options['-m'].set_value(iqtree_model)
        if no_fast:
            iqtree_wrapper.options.options['-fast'].active = False
        if iqtree_extra:
            iqtree_wrapper.options.options['_extra'] = StringOption('', iqtree_extra, active=True)
        tree = iqtree_wrapper()
        if tree:
            with open(treefile, 'w') as fh:
                fh.write(tree.strip() + '\n')
            return tree
    except Exception as e:
        logger.error('Gene tree failed for {}: {}'.format(og_name, e))
    return None


class CoalescentInference(object):
    """
    Orchestrates the coalescent species tree pipeline (step 4astral):

    1. Filter per-OG alignments from 06_align_merge_aa by gap fraction and
       taxon occupancy, writing clean FASTA to 07_astral_filtered_aa.
    2. Optionally trim filtered alignments with ClipKIT (--trim flag),
       writing results to 07_astral_trimmed_aa.
    3. Run IQ-TREE on each alignment in parallel using multiprocessing,
       writing individual gene treefiles to 08_gene_trees.
    4. Concatenate gene trees into a single treefile and run ASTER to
       produce the coalescent species tree.
    """

    def __init__(self, args):
        self.args = args
        self._species_name = 'merge'
        self._seq_type = 'dna' if getattr(args, 'dna', False) else 'aa'
        self._filtered_folder = self._make_output_path('07_astral_filtered_' + self._seq_type)
        self._trimmed_folder = self._make_output_path('07_astral_trimmed_' + self._seq_type) if args.trim else None
        self._gene_trees_folder = self._make_output_path('08_gene_trees_' + self._seq_type)
        self.elapsed_time = 0
        self.tree = None
        self._run()

    def _run(self):
        start = time.time()

        filtered_files = self._filter_alignments()
        if not filtered_files:
            logger.error('{}: No alignments passed filtering for step 4astral.'.format(self._species_name))
            return

        if self.args.trim:
            input_files = self._trim_alignments(filtered_files)
            if not input_files:
                logger.error('{}: No alignments remain after ClipKIT trimming.'.format(self._species_name))
                return
        else:
            input_files = filtered_files

        gene_tree_file = self._infer_gene_trees(input_files)
        self.tree = self._infer_species_tree(gene_tree_file)

        end = time.time()
        self.elapsed_time = end - start
        logger.info('{}: Step 4astral coalescent inference took {:.2f}s.'.format(
            self._species_name, self.elapsed_time))

    def _make_output_path(self, prefix):
        path = os.path.join(self.args.output_path, prefix)
        if not os.path.exists(path):
            os.makedirs(path)
        return path

    def _filter_alignments(self):
        """
        Filter per-OG alignments from 06_align_merge_aa by gap fraction and occupancy.

        The .fa files produced by step 3combine are phylip-relaxed format despite their
        extension; this method reads them correctly and writes passing alignments as
        standard FASTA to 07_astral_filtered_aa.

        :return: list of paths to filtered FASTA files
        """
        input_folder = os.path.join(self.args.output_path, '06_align_merge_' + self._seq_type)
        log_path = os.path.join(self._filtered_folder, 'filtering_summary.txt')
        filtered_files = []
        total = 0
        passed = 0
        dropped = 0

        with open(log_path, 'w') as log:
            log.write('OG_Name\tOriginal_Sequences\tSequences_Passing_Gap_Filter\tStatus\n')
            for filepath in sorted(glob.glob(os.path.join(input_folder, '*.fa'))):
                total += 1
                og_name = os.path.basename(filepath)
                orig_count = 0
                valid_records = []
                try:
                    alignment = AlignIO.read(filepath, 'phylip-relaxed')
                    for record in alignment:
                        orig_count += 1
                        seq_str = str(record.seq).upper()
                        gap_count = seq_str.count('-') + seq_str.count('X') + seq_str.count('N')
                        if len(seq_str) > 0 and (gap_count / len(seq_str)) <= self.args.max_gap:
                            valid_records.append(record)

                    if len(valid_records) >= self.args.min_samples:
                        out_path = os.path.join(self._filtered_folder,
                                                og_name.replace('.fa', '.fasta'))
                        SeqIO.write(valid_records, out_path, 'fasta')
                        filtered_files.append(out_path)
                        passed += 1
                        status = 'KEPT'
                    else:
                        dropped += 1
                        status = 'DROPPED (Low Occupancy)'
                except Exception as e:
                    dropped += 1
                    status = 'ERROR: {}'.format(e)

                log.write('{}\t{}\t{}\t{}\n'.format(og_name, orig_count, len(valid_records), status))

            log.write('\n=== TOTALS ===\n')
            log.write('Total evaluated: {}\nKept: {}\nDropped/failed: {}\n'.format(
                total, passed, dropped))
            if total > 0:
                log.write('Retention rate: {:.2f}%\n'.format((passed / total) * 100))

        logger.info('{}: Alignment filtering kept {} of {} OGs.'.format(
            self._species_name, passed, total))
        return filtered_files

    def _trim_alignments(self, filtered_files):
        """
        Run ClipKIT on each filtered alignment.

        :param filtered_files: list of FASTA alignment paths
        :return: list of trimmed FASTA paths that are non-empty after trimming
        """
        from read2tree.wrappers.aligners.clipkit import Clipkit
        trimmed_files = []
        for fasta_file in filtered_files:
            og_name = os.path.basename(fasta_file)
            trimmed_path = os.path.join(self._trimmed_folder, og_name)
            try:
                clipkit_wrapper = Clipkit(fasta_file, trimmed_path)
                result = clipkit_wrapper()
                if result:
                    trimmed_files.append(result)
                else:
                    logger.warning('{}: ClipKIT produced empty output for {}, skipping.'.format(
                        self._species_name, og_name))
            except WrapperError as e:
                logger.error('{}: ClipKIT failed for {}: {}'.format(self._species_name, og_name, e))
        logger.info('{}: ClipKIT trimming kept {} of {} alignments.'.format(
            self._species_name, len(trimmed_files), len(filtered_files)))
        return trimmed_files

    def _infer_gene_trees(self, alignment_files):
        """
        Run IQ-TREE on each alignment in parallel using multiprocessing.Pool.
        Collects all gene trees into a single Newick file for ASTER.

        :param alignment_files: list of FASTA alignment paths
        :return: path to the concatenated gene treefile
        """
        iqtree_model = getattr(self.args, 'iqtree_model', None)
        iqtree_extra = getattr(self.args, 'iqtree_args', None)
        no_fast = getattr(self.args, 'no_fast', False)
        tasks = [(f, self._gene_trees_folder, iqtree_model, iqtree_extra, no_fast, self._seq_type) for f in alignment_files]
        logger.info('{}: Running per-gene IQ-TREE on {} alignments with {} workers.'.format(
            self._species_name, len(tasks), self.args.threads))

        p = Pool(self.args.threads)
        results = p.map(_run_gene_tree, tasks)
        p.close()
        p.join()

        trees = [t for t in results if t is not None]
        gene_tree_file = os.path.join(self.args.output_path,
                                      'gene_trees_{}_{}.nwk'.format(self._species_name, self._seq_type))
        with open(gene_tree_file, 'w') as fh:
            for tree in trees:
                fh.write(tree.strip() + '\n')

        logger.info('{}: {} of {} gene trees successfully inferred.'.format(
            self._species_name, len(trees), len(tasks)))
        return gene_tree_file

    def _infer_species_tree(self, gene_tree_file):
        """
        Run ASTER on the collected gene trees to estimate a coalescent species tree.

        :param gene_tree_file: path to file containing one gene tree (Newick) per line
        :return: coalescent species tree in Newick format
        """
        species_tree_file = os.path.join(self.args.output_path,
                                         'astral_tree_{}_{}.nwk'.format(self._species_name, self._seq_type))
        aster_wrapper = Aster(gene_tree_file, species_tree_file,
                              binary=getattr(self.args, 'astral_binary', None))
        aster_wrapper.options.options['-t'].set_value(self.args.threads)
        astral_extra = getattr(self.args, 'astral_args', None)
        if astral_extra:
            aster_wrapper.options.options['_extra'] = StringOption('', astral_extra, active=True)
        tree = aster_wrapper()
        logger.info('{}: Coalescent species tree written to {}'.format(
            self._species_name, species_tree_file))
        return tree
