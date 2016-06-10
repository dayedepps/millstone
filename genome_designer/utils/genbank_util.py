from Bio import SeqIO

def get_genbank_features_with_type(genbank_path, feature_type):

    chrom_intervals = {}
    with open(genbank_path, 'r') as fh:
        for seq_record in SeqIO.parse(fh, 'genbank'):
            interval_list = []
            for f in seq_record.features:
                if f.type == feature_type:
                    interval_list.append((f, f.extract(seq_record.seq)))

            chrom_intervals[seq_record.id] = interval_list

    return chrom_intervals


def generate_genbank_mobile_element_multifasta(genbank_path, output_fasta_path):
    """Extract mobile elements from a genbank-annotated genome into a multifasta
    record, for use with SV and mobile element-calling code.
    """
    me_features = get_genbank_features_with_type(
            genbank_path, 'mobile_element')
    me_sequences = {}
    for chrom, feature_seq_tuples in me_features.items():
        for feature, seq in feature_seq_tuples:
            if 'mobile_element_type' not in feature.qualifiers:
                continue
            me_type = feature.qualifiers['mobile_element_type'][0]
            me_type = me_type.replace(' ', '_')
            if me_type not in me_sequences:
                me_sequences[me_type] = seq

    with open(output_fasta_path, 'w') as fh:
        for me_type, seq in me_sequences.items():
            fh.write('>ME_%s\n' % me_type)
            fh.write(str(seq))
            fh.write('\n')


def get_overlapping_features(genbank_path, chromosome, interval):
    """
    TODO: This will be useful for annotating SVs and other quick
    gene-centric UI annotations.
    """
    raise NotImplementedError