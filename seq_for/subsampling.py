from __future__ import division
import random
import click
import gzip
from tqdm import tqdm


"""
# gz is ok
# fastq is ok

seqtk sample -s100 ~/data_bank/ZHJ_WES/raw_data/HUA_1.fastq.gz 10000 > ~/tools/Whole_pipelines/test_data/germline/s1_1.fq;
seqtk sample -s100 ~/data_bank/ZHJ_WES/raw_data/HUA_2.fastq.gz 10000 > ~/tools/Whole_pipelines/test_data/germline/s1_2.fq
"""

# todo: seqtk could do that....
def subsampling_fq(path):
    if path.endswith('.gz'):
        file_handle = gzip.open(path,mode='rt')
    else:
        file_handle = open(path,mode='r')


number_to_sample = 3000000
number_of_replicates = 10

with open("test.fastq") as input:
    num_lines = sum([1 for line in input])
total_records = int(num_lines / 4)
print("sampling " + str(number_to_sample) + " out of " + str(total_records) + " records")

output_files = []
output_sequence_sets = []
for i in range(number_of_replicates):
    output_files.append(open("sample.fastq." + str(i), "w"))
    output_sequence_sets.append(set(random.sample(xrange(total_records + 1), number_to_sample)))

record_number = 0
with open("test.fastq") as input:
        for line1 in input:
            line2 = input.next()
            line3 = input.next()
            line4 = input.next()
            for i, output in enumerate(output_files):
                if record_number in output_sequence_sets[i]:
                        output.write(line1)
                        output.write(line2)
                        output.write(line3)
                        output.write(line4)
                record_number += 1

for output in output_files:
    output.close()