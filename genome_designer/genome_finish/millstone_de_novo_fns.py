import subprocess
import os

from django.conf import settings
import numpy as np
import pysam

from genome_finish import __path__ as gf_path_list
from genome_finish.insertion_placement_read_trkg import extract_left_and_right_clipped_read_dicts
from pipeline.read_alignment import _filter_out_interchromosome_reads
from pipeline.read_alignment_util import index_bam_file
from main.models import Dataset
from main.models import Variant
from main.models import VariantSet
from utils.bam_utils import clipping_stats
from variants.variant_sets import update_variant_in_set_memberships

GENOME_FINISH_PATH = gf_path_list[0]
VELVETH_BINARY = settings.TOOLS_DIR + '/velvet/velveth'
VELVETG_BINARY = settings.TOOLS_DIR + '/velvet/velvetg'


def get_altalign_reads(input_bam_path, output_bam_path, xs_threshold=None):
    input_af = pysam.AlignmentFile(input_bam_path, 'rb')
    output_af = pysam.AlignmentFile(output_bam_path, 'wb',
                template=input_af)

    for read in input_af:
        if read.has_tag('XS') and read.has_tag('AS'):
            if read.get_tag('AS') <= read.get_tag('XS'):
                output_af.write(read)

    output_af.close()
    input_af.close()



def get_piled_reads(input_bam_path, output_bam_path,
        clipping_threshold=None):
    """Creates bam of reads that have more than clipping_threshold bases
    of clipping and are stacked 3 standard deviations higher than the average
    pileup of clipped reads.  If no clipping_threshold specified, clipping
    stats for the alignment are calculated and the clipping_threshold is set
    to the mean + one stddev of the per read clipping of a sample of
    10000 reads.
    """
    if clipping_threshold is None:
        stats = clipping_stats(input_bam_path, sample_size=10000)
        clipping_threshold = int(stats['mean'] + stats['std'])

    input_af = pysam.AlignmentFile(input_bam_path, 'rb')
    output_af = pysam.AlignmentFile(output_bam_path, 'wb',
            template=input_af)
    lr_clipped = extract_left_and_right_clipped_read_dicts(
            input_af,
            clipping_threshold=clipping_threshold)
    input_af.close()

    for clipped_dict in [
            lr_clipped['left_clipped'], lr_clipped['right_clipped']]:
        stack_counts = map(len, clipped_dict.values())
        mean_stacking = np.mean(stack_counts)
        std_stacking = np.std(stack_counts)
        stacking_cutoff = mean_stacking + 3 * std_stacking
        for read_list in clipped_dict.values():
            if len(read_list) > stacking_cutoff:
                for read in read_list:
                    output_af.write(read)
    output_af.close()


def get_clipped_reads_smart(input_bam_path, output_bam_path,
        clipping_threshold=8, phred_encoding=None):

    """Gets reads not overlapping their adaptor with a terminal
    segment of clipping with average phred scores above the cutoff
    """

    phred_encoding_to_shift = {
        'Illumina 1.5': 31,
        'Sanger / Illumina 1.9': 0
    }

    CLIPPED_AVG_PHRED_CUTOFF = 20
    if (phred_encoding is not None and
        phred_encoding in phred_encoding_to_shift):

        CLIPPED_AVG_PHRED_CUTOFF += phred_encoding_to_shift[phred_encoding]

    SOFT_CLIP = 4
    HARD_CLIP = 5
    CLIP = [SOFT_CLIP, HARD_CLIP]

    input_af = pysam.AlignmentFile(input_bam_path, 'rb')
    output_af = pysam.AlignmentFile(output_bam_path, 'wb',
                template=input_af)

    for read in input_af:
        # If no cigartuples, i.e. unmapped, continue
        if read.cigartuples is None:
            continue

        if read.is_secondary or read.is_supplementary:
            continue

        # TODO: Account for template length
        # adapter_overlap = max(read.template_length - query_alignment_length, 0)

        # Determine left and right clipped counts
        left_clipping = (read.cigartuples[0][1]
                if read.cigartuples[0][0] in CLIP else 0)
        right_clipping = (read.cigartuples[-1][1]
                if read.cigartuples[-1][0] in CLIP else 0)

        # Write reads to file if clipped bases have average phred score
        # above cutoff
        if left_clipping > clipping_threshold:
            clipped_phred_scores = read.query_qualities[:left_clipping]
            if np.mean(clipped_phred_scores) > CLIPPED_AVG_PHRED_CUTOFF:
                output_af.write(read)
                continue
        if right_clipping > clipping_threshold:
            clipped_phred_scores = read.query_qualities[-right_clipping:]
            if np.mean(clipped_phred_scores) > CLIPPED_AVG_PHRED_CUTOFF:
                output_af.write(read)
                continue

    output_af.close()
    input_af.close()


def get_unmapped_reads(bam_filename, output_filename, avg_phred_cutoff=None):

    if avg_phred_cutoff is not None:
        intermediate_filename = '_unfiltered'.join(
                os.path.splitext(output_filename))
    else:
        intermediate_filename = output_filename

    cmd = '{samtools} view -h -b -f 0x4 {bam_filename}'.format(
            samtools=settings.SAMTOOLS_BINARY,
            bam_filename=bam_filename)
    with open(intermediate_filename, 'w') as output_fh:
        subprocess.call(
                cmd, stdout=output_fh, shell=True,
                executable=settings.BASH_PATH)

    if avg_phred_cutoff is not None:
        filter_low_qual_read_pairs(intermediate_filename, output_filename,
                avg_phred_cutoff)


<<<<<<< HEAD
def get_split_reads(bam_filename, output_filename):
    """Isolate split reads from a sample alignment.

    This uses a python script supplied with Lumppy, that is run as a
    separate process.

    NOTE THAT THIS SCRIPT ONLY WORKS WITH BWA MEM.
    """

    # Use lumpy bwa-mem split read script to pull out split reads.
    filter_split_reads = ' | '.join([
            # -F 0x100 option filters out secondary alignments
            '{samtools} view -h -F 0x100 {bam_filename}',
            'python {lumpy_bwa_mem_sr_script} -i stdin',
            '{samtools} view -Sb -']).format(
                    samtools=SAMTOOLS_BINARY,
                    bam_filename=bam_filename,
                    lumpy_bwa_mem_sr_script=
                            settings.LUMPY_EXTRACT_SPLIT_READS_BWA_MEM)

    try:
        with open(output_filename, 'w') as fh:
            subprocess.check_call(
                    filter_split_reads, stdout=fh, shell=True,
                    executable=BASH_PATH)

        # sort the split reads, overwrite the old file
        subprocess.check_call(
                [SAMTOOLS_BINARY, 'sort', output_filename,
                 os.path.splitext(output_filename)[0]])

    except subprocess.CalledProcessError:
        raise Exception('Exception caught in split reads generator, ' +
                        'perhaps due to no split reads')

def get_discordant_read_pairs(bam_filename, output_filename):
    """Isolate discordant pairs of reads from a sample alignment.
    """
    assert os.path.exists(bam_filename), "BAM file '%s' is missing." % (
            bam_filename)

    # NOTE: This assumes the index just adds at .bai, w/ same path otherwise
    # - will this always be true?
    if not os.path.exists(bam_filename+'.bai'):
        index_bam_file(bam_filename)

    # Use bam read alignment flags to pull out discordant pairs only
    filter_discordant = ' | '.join([
            '{samtools} view -u -F 0x0002 {bam_filename} ',
            '{samtools} view -u -F 0x0100 - ',
            '{samtools} view -u -F 0x0004 - ',
            '{samtools} view -u -F 0x0008 - ',
            '{samtools} view -b -F 0x0400 - ']).format(
                    samtools=settings.SAMTOOLS_BINARY,
                    bam_filename=bam_filename)

    try:
        with open(output_filename, 'w') as fh:
            subprocess.check_call(filter_discordant,
                    stdout=fh, shell=True, executable=settings.BASH_PATH)

        # sort the discordant reads, overwrite the old file
        subprocess.check_call([settings.SAMTOOLS_BINARY, 'sort', output_filename,
                os.path.splitext(output_filename)[0]])

        _filter_out_interchromosome_reads(output_filename)

    except subprocess.CalledProcessError:
        raise Exception('Exception caught in discordant reads generator, '+
                'perhaps due to no discordant reads')


def _parse_sam_line(line):
    parts = line.split()
    return {
        'read_id': parts[0],
        'flags': parts[1]
    }