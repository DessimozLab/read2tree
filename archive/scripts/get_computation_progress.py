#!/usr/bin/env python

from __future__ import unicode_literals, division
import collections
import math
import os
import re
import multiprocessing
import logging
import itertools
logger = logging.getLogger(__name__)


def load_dataset(db_dir):
    genomes = [g[0:g.index('.')] for g in sorted(os.listdir(db_dir)) if g.endswith('.fa')]
    return genomes


def get_nr_protein(genome_db_path):
    with open(genome_db_path, 'rb') as fh:
        return sum(map(lambda line: line.startswith(b'<E>'), fh))


def get_batch_size(paramfile):
    pat = re.compile(b'AlignBatchSize\s*:=\s*(?P<size>[\d.e]+)[:;]')
    with open(paramfile, 'rb') as fh:
        for line in fh:
            m = pat.search(line)
            if m is not None:
                return int(float(m.group('size')))


def progress_of_pair(arg_tuple):
    g1, g2, allall, expected = arg_tuple
    compact = os.path.join(allall, g1, g2 + '.gz')
    if os.path.isfile(compact):
        is_exported = os.path.exists(os.path.join(allall, g1, g2 + '.sha2.gz'))
        return g1, g2, expected, 0, is_exported
    part_directory = os.path.join(allall, g1, g2)
    if not os.path.isdir(part_directory):
        return g1, g2, 0, 0, False
    files = os.listdir(part_directory)
    finished = sum(map(lambda x: x.endswith('.gz'), files))
    return g1, g2, finished, len(files) - finished, False


def build_genome_entrycount_dict(genomes, path):
    pool = multiprocessing.Pool()
    dbnames = [os.path.join(path, 'Cache', 'DB', g + '.db') for g in genomes]
    mapping = collections.OrderedDict(
        zip(genomes, pool.map(get_nr_protein, dbnames)))
    pool.close()
    pool.join()
    return mapping


class Dataset(object):
    def __init__(self, path, report_individual=False):
        self.path = path
        genomes = load_dataset(os.path.join(path, 'DB'))
        self.batch_size = get_batch_size(os.path.join(path, 'parameters.drw'))
        self.progress_results = []
        self.report_individual = report_individual
        self.genomes = build_genome_entrycount_dict(genomes, path)

    def expected_nr_chunks(self, g1, g2):
        if g1 == g2:
            nr_alignments = self.genomes[g1] * (self.genomes[g1] - 1) / 2
        else:
            nr_alignments = self.genomes[g1] * self.genomes[g2]
        return math.ceil(nr_alignments / self.batch_size)

    def genome_order(self, g1, g2):
        if (self.genomes[g1] < self.genomes[g2]) or (self.genomes[g1] == self.genomes[g2] and g1 < g2):
            return g1, g2
        else:
            return g2, g1

    def pairs_generator(self):
        allall_path = os.path.join(self.path, 'Cache', 'AllAll')
        for g1, g2 in itertools.combinations_with_replacement(self.genomes.keys(), 2):
            g1, g2 = self.genome_order(g1, g2)
            yield (g1, g2, allall_path, self.expected_nr_chunks(g1, g2))

    def compute_progress(self, nr_cpus=None):
        results = Reporter(self, self.report_individual)
        all_pairs = self.pairs_generator()
        pool = multiprocessing.Pool(nr_cpus)
        while True:
            pair_chunk = list(itertools.islice(all_pairs, 100))
            if len(pair_chunk) == 0:
                break
            stat = pool.map_async(progress_of_pair, pair_chunk, callback=results.handle_progress_result)
        pool.close()
        stat.wait()
        pool.join()
        results.print_report()
        return results


class Reporter(object):
    def __init__(self, dataset, report_individual_pairs):
        self.progress_results = []
        self.header_printed = False
        self.report_individual_pairs = report_individual_pairs
        self.dataset = dataset

    def handle_progress_result(self, result):
        logger.debug('handling result chunk of size {}'.format(len(result)))
        for pair in result:
            self.progress_results.append(pair)
            if self.report_individual_pairs:
                if not self.header_printed:
                    print(
                        "Genome1\tGenome2\tExported\tExpected chunks\tFinished chunks\tStarted chunks\tFraction of finished chunks")
                    self.header_printed = True
                expect = self.dataset.expected_nr_chunks(pair[0], pair[1])
                print("{}\t{}\t{}\t{}\t{}\t{}\t{:%}"
                      .format(pair[0], pair[1], pair[4], int(expect), int(pair[2]), int(pair[3]), pair[2] / expect))

    def print_report(self):
        tot = started = done = tot_no_export = done_no_export = 0
        for pair in self.progress_results:
            done += pair[2]
            started += pair[3]
            tot += self.dataset.expected_nr_chunks(pair[0], pair[1])
            if not pair[4]:
                tot_no_export += self.dataset.expected_nr_chunks(pair[0], pair[1])
                done_no_export += pair[2]
        print("\n\nSummary of OMA standalone All-vs-All computations:")
        print("--------------------------------------------------")
        print("Nr chunks started: {:d} ({:.2%})".format(int(started), started / tot))
        print("Nr chunks finished: {:d} ({:.2%})".format(int(done), done / tot))
        print("Nr chunks finished w/o exported genomes: {:d} ({:.2%})"
              .format(int(done_no_export), done_no_export / tot_no_export))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="""Compute overall status of OMA standalone AllAll computations. 
                    This script checks for all genome pairs to be computed the overall
                    progress. The reporting unit is the number of alignment chunks.""")
    parser.add_argument('root', nargs="?", default="./",
                        help="Path to the project root of the analysis, i.e. the directory"
                             "that contains the input genomes DB/ directory (default: %(default)s).")
    parser.add_argument('-n', '--nr-cpus', type=int,
                        help="use at most that many processes in parallel to calculate status. " +
                             "If not set, the number of available CPUs on the executing host will be used.")
    parser.add_argument('-i', '--individual', action="store_true", default=False,
                        help="Report a tsv-formatted line per genome pair with the current progress")
    parser.add_argument('-d', '--debug', action='store_true', help="increase output to debug level")
    conf = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG if conf.debug else logging.INFO,
                        format="%(asctime)-15s %(levelname)8s: %(message)s")
    d = Dataset(conf.root, report_individual=conf.individual)
    d.compute_progress(conf.nr_cpus)
