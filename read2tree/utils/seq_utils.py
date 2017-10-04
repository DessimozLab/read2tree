#!/usr/bin/env python
'''
    Most of the functions here were taken from the zoo


    -- David Dylus, August--XXX 2017
'''

import os
import types
from collections import defaultdict

from enum import Enum
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq, UnknownSeq
from Bio.SeqRecord import SeqRecord

__all__ = ['is_dna', 'guess_datatype', 'identify_input']


def is_dna(record):
    """check whether a sequence is of type dna"""
    allchars = [char.upper() for char in str(record.seq) if char.upper() not in '-?X']
    return set(allchars).issubset(set(IUPAC.ambiguous_dna.letters))


def guess_datatype(alignment, from_filename=False):
    logger.warning("Guessing is not recommended - specify the sequence type with option datatype={DNA, PROTEIN}, be more confident")
    if from_filename:
        try:
            alignment = SeqIO.parse(alignment, 'fasta')
        except:
            alignment = SeqIO.parse(alignment, 'phylip-relaxed')
    return DataType.DNA if is_dna(alignment) else DataType.PROTEIN

AlignmentInput = Enum('AlignmentInput', 'OBJECT FILENAME')

def identify_input(alignment):
    """
    Work out if we're dealing with a Biopython object (return True), a file
    (return False), or invalid input (raise error)
    """
    try:
        if isinstance(alignment, (MultipleSeqAlignment, types.GeneratorType, list)):
            # `alignment` is a Biopython MultipleSequenceAlignment
            return AlignmentInput.OBJECT

        elif isinstance(alignment, str) and os.path.exists(alignment):
            # `alignment` is a filepath
            return AlignmentInput.FILENAME

    except:
        # `alignment` is some other thing we can't handle
        raise ValueError('{} is not an alignment object or a valid filename'.format(alignment))

def concatenate(alignments):
    """
    Concatenates a list of multiple sequence alignment objects.

    The alignments are concatenated based on their label, i.e. the
    sequences from the different alignments which have the same id/labels
    will become a single sequence. The order is preserved.

    If any sequences are missing in one or several alignments, these parts
    are padded with unknown data (:py:class:`Bio.Seq.UnknownSeq`).

    :param alignments: the list of alignments objects, i.e. list(:py:class:`Bio.Align.MultipleSeqAlignment`)
    :returns: a single :py:class:`Bio.Align.MultipleSeqAlignment`

    Example::

        >>> sequences = {'aln1': {'seq1': 'acgtca',
        ...                       'seq2': 'acgtt-',
        ...                       'seq3': 'ac-ta-'},
        ...              'aln2': {'seq2': 'ttg-cta',
        ...                       'seq3': 'tcgacta',
        ...                       'seq4': 'ttgacta'}}
        >>> alignments = [MultipleSeqAlignment([SeqRecord(Seq(sequence,
        ...                    alphabet=IUPAC.extended_dna), id=key)
        ...      for (key, sequence) in sequences[aln].items()])
        ...               for aln in ('aln1', 'aln2')]
        >>> con_alignment = concatenate(alignments)
        >>> con_alignment.sort()
        >>> print(con_alignment)
        ExtendedIUPACDNA() alignment with 4 rows and 13 columns
        acgtcaNNNNNNN seq1
        acgtt-ttg-cta seq2
        ac-ta-tcgacta seq3
        NNNNNNttgacta seq4

    :note:

       Limitations: any annotations in the sub-alignments are lost in
       the concatenated alignment.

    """

    # First check to see whether we're inputting filenames of alignments or the Biopython alignments
    # Assume that it's a biopython alignment if it's not a filename
    tmp_aligns = []
    for filename in alignments:
        if identify_input(filename).name == 'FILENAME':
            tmp_aligns.append(AlignIO.read(filename, "fasta"))
        else:
            tmp_aligns.append(filename)

    # Copy back to alignments
    alignments = tmp_aligns

    # Get the full set of labels (i.e. sequence ids) for all the alignments
    all_labels = set(seq.id for aln in alignments for seq in aln)

    # Make a dictionary to store info as we go along
    # (defaultdict is convenient -- asking for a missing key gives back an empty list)
    tmp = defaultdict(list)

    # Assume all alignments have same alphabet
    alphabet = alignments[0]._alphabet

    for aln in alignments:
        length = aln.get_alignment_length()

        # check if any labels are missing in the current alignment
        these_labels = set(rec.id for rec in aln)
        missing = all_labels - these_labels

        # if any are missing, create unknown data of the right length,
        # stuff the string representation into the tmp dict
        for label in missing:
            new_seq = UnknownSeq(length, alphabet=alphabet)
            tmp[label].append(str(new_seq))

        # else stuff the string representation into the tmp dict
        for rec in aln:
            tmp[rec.id].append(str(rec.seq))

    # Stitch all the substrings together using join (most efficient way),
    # and build the Biopython data structures Seq, SeqRecord and MultipleSeqAlignment
    return MultipleSeqAlignment(SeqRecord(Seq(''.join(v), alphabet=alphabet), id=k)
                                for (k, v) in tmp.items())