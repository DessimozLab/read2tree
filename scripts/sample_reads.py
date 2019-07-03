from tqdm import tqdm
import os, shutil, time, random, tempfile
import numpy as np
from read2tree.FastxReader import FastxReader
from math import ceil
import logging, yaml
from pkg_resources import resource_string
logging.getLogger(__name__).addHandler(logging.NullHandler())
conf = resource_string(__name__, '../read2tree/logging/log.yaml')

D = yaml.load(conf)
D.setdefault('version', 1)
logging.config.dictConfig(D)

logger = logging.getLogger(__name__)

def sample_from_reads(reads, genome_len, coverage):
    '''
    Main function taking in the reads of the object and processing it
    given the provided parameters
    :return: string that contains all the read sequences separated by '\n'
    '''
    sampled_reads = []
    start = time.time()
    if len(reads) == 1:
        idx_random = _get_vector_random_reads(reads, genome_len, coverage)
        if idx_random:
            sampled_reads = _sample_read_file(reads, idx_random)
        else:
            sampled_reads = reads
    elif len(reads) == 2:
        idx_random = _get_vector_random_reads(reads[0], genome_len, coverage)
        if idx_random:
            sampled_reads.append(_sample_read_file(reads[0], idx_random))
            sampled_reads.append(_sample_read_file(reads[1], idx_random))
            logger.info('Reads can be found at: {}'.format(sampled_reads))
            shutil.copy(sampled_reads[0], reads[0].split(".")[0]+".fa")
            shutil.copy(sampled_reads[1], reads[1].split(".")[0]+".fa")
        else:
            mate_pairs = _check_read_consistency(reads)
            if mate_pairs:
                sampled_reads = _select_mates_from_reads(reads, mate_pairs)
                shutil.copy(sampled_reads[0], reads[0].split(".")[0] + ".fa")
                shutil.copy(sampled_reads[1], reads[1].split(".")[0] + ".fa")
            else:
                sampled_reads = reads

    end = time.time()
    elapsed_time = end - start
    logger.info('Sampling of reads took {}.'.format(elapsed_time))

    return sampled_reads


def _check_read_consistency(reads):
    '''
    Function that checks whether all mate pairs are present and if not
    uses only reads for which mate pairs exist.
    '''
    print('--- Checking for consistent mate pairing ---')
    left_read = FastxReader(reads[0])
    right_read = FastxReader(reads[1])

    with left_read.open_fastx() as left_input:
        with right_read.open_fastx() as right_input:
            with_mate_pairs = set(left_read.readfq_id(left_input)) \
                .intersection(right_read.readfq_id(right_input))
            logger.info('Mate pairs have size: {}'.format(sys.getsizeof(with_mate_pairs)))

    with left_read.open_fastx() as left_input:
        len_left = len(set(left_read.readfq_id(left_input)))

    if len_left == len(with_mate_pairs):
        logger.info('Mate pairs are consistent.')
        return None
    else:
        logger.info('----> Mate pairing not consitent! ---')
        logger.info('Inconsistent number of mate pairs! '
                    'Will use only reads that have mate pair. '
                    'Consistent {} of {} total reads.'
                    .format(2*len(with_mate_pairs),
                            len_left+len_left))
        return with_mate_pairs

def _select_mates_from_reads(reads, mates):
    '''
    Main function taking in the reads of the object and processing it
    given the provided parameters
    :return: string that contains all the read sequences separated by '\n'
    '''
    sampled_reads = []
    start = time.time()
    sampled_reads.append(_sample_read_file(reads[0],
                                                mates))
    sampled_reads.append(_sample_read_file(reads[1],
                                                mates))

    end = time.time()
    elapsed_time = end - start
    logger.info('Selecting of reads with mates took {}.'
                .format(elapsed_time))
    return sampled_reads



def _get_num_reads_by_coverage(file, num_records, genome_len, coverage):
    read_len = _get_read_len(file, num_records)
    logger.info('Average read length estimated to {}.'.format(read_len))
    return int(ceil(genome_len * coverage / (2 * read_len)))


def _get_vector_random_reads(file, genome_len, coverage):
    total_records = _get_num_reads(file)
    num_reads_by_coverage = _get_num_reads_by_coverage(file, total_records, genome_len, coverage)
    logger.info('Sampling {} / {} reads for {}X coverage.'.format(num_reads_by_coverage, total_records, coverage))
    if num_reads_by_coverage > total_records:
        logger.info("Not enough reads available for sampling, using them all.")

        return None
    else:
        x = np.sort(np.random.choice((total_records + 1), num_reads_by_coverage, replace=False))
        return list(x)
        #old: return sorted(list(set(random.sample(range(total_records + 1), num_reads_by_coverage))))


def _get_num_reads(file):
    fastq_reader = FastxReader(file)
    num_lines = 0
    with fastq_reader.open_fastx() as f:
        for l in f: num_lines += 1
    return int(num_lines / 4)
        

def _get_read_len(file, num_records):
    fastq_reader = FastxReader(file)
    if num_records < 10000:
        samples = num_records
    else:
        samples = 10000
    lens = [0] * samples
    with fastq_reader.open_fastx() as f:
        for i, l in enumerate(fastq_reader.readfq(f)):
            if i < samples:
                lens[i] = len(l[1])
            else:
                break
        mean_len = np.mean(lens)
        median_len = np.median(lens)
    logger.info('The reads have a mean length of {} '
          'and a median length of {}.'
          .format(mean_len, median_len))
    return mean_len


def _get_2_line_fasta_string(read_id, seq, x=None):
    '''
    Transform 2 lines of read string to new read string providing the
    split information
    :param read_id: Read ID in the form of SRR00001
    :param read_num: Number of read usually after the read ID
    :param x: Numerical iterator
    :param seq: Sequence string
    :return: 2 lines that correspond to one read with adapted ID
    '''
    out = ''
    if x:
        new_name = ">" + read_id.replace('@','') + "_" + str(x) + ' length=' + \
            str(len(seq))
    else:
        new_name = ">" + read_id.replace('@','') + ' length=' + str(len(seq))
    out += new_name + "\n"
    out += seq + "\n"
    return out


def _sample_read_file(file, select_idx):
    initial_length = 0
    sampling_length = 0
    select_idx = list(select_idx)
    out_file = tempfile.NamedTemporaryFile(mode='at', suffix='.fa',
                                           delete=False)
    fastq_reader = FastxReader(file)
    with fastq_reader.open_fastx() as read_input:
        k = 0
        for i, f_input in enumerate(tqdm(fastq_reader.readfq(read_input),desc='Selecting reads from {}'.format(os.path.basename(file)),unit=' reads')):
            name = f_input[0]
            seq = f_input[1]
            initial_length = initial_length + len(seq)
            try:
                if i == select_idx[k]:
                    # logger.info(name)
                    seq_record_str = _get_2_line_fasta_string(name, seq)
                    out_file.write(seq_record_str)
                    sampling_length = sampling_length + len(seq)
                    k = k + 1
            except IndexError as e:
                logger.info('i={}, k={}, len={} and last select_idx={} has error:{}'.format(i,k,len(select_idx),select_idx[-1],e))
                pass
            if k == len(select_idx):
                break
    logger.info('Cummulative length of all reads {}bp. Cummulative length of sampled reads {}bp'.format( initial_length,sampling_length))
    out_file.close()
    return out_file.name

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="""Sample reads from files""")
    parser.add_argument('--oxml', default=None,
                        help='[Default is none] Remove species present '
                        'in data set after mapping step completed to '
                        'build OGs. Input is comma separated list '
                        'without spaces, e.g. XXX,YYY,AAA.')
    parser.add_argument('--coverage', type=float, default=10,
                            help='[Default is 10] coverage in X.')
    parser.add_argument('--genome_len', type=int, default=2000000,
                            help='[Default is 2000000] Genome size in bp.')
    parser.add_argument('--reads', nargs='+', default=None,
                            help='Reads to be mapped to reference. If paired '
                            'end add separated by space.')

    conf = parser.parse_args()
    logger.info('----- WE ARE SAMPLING THE READS -----')
    sample_from_reads(conf.reads, conf.genome_len, conf.coverage)
