import io
import os
import pandas as pd

def compare_seq(subject, query, type_blast='blastn'):
    tab_text = os.popen("%s -subject '%s' -query '%s' -outfmt 6" % (type_blast, subject, query)).read()
    data = pd.read_csv(io.StringIO(tab_text), sep='\t', index_col=0, header=None)
    header = ['qaccver',
              'saccver',
              'pident',
              'length',
              'mismatch',
              'gapopen',
              'qstart',
              'qend',
              'sstart',
              'send',
              'evalue',
              'bitscore']
    data.columns = header
    return data

