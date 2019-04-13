import os
import time

from Bio import Entrez

Entrez.email = 'l0404th@gmail.com'


def fetch(nid, start, end, output):
    record = ''
    while not record:
        try:
            record = Entrez.efetch(db='nuccore', id=nid, retmode='text', rettype='fasta', seq_start=start,
                                   seq_stop=end).read()
            record = record.split('\n')
            record = record[0] + '\n' + ''.join(record[1:])
        except:
            time.sleep(2)
    if os.path.isfile(output):
        if '(' in output:
            output_n = output.partition('(')[0]
            output_number = int(output.split('(')[-1].strip(')'))
            output = output_n + '(%s)' % str(output_number + 1)
        else:
            output = output + '(2)'
    with open(output, 'w') as f1:
        f1.write(record)


