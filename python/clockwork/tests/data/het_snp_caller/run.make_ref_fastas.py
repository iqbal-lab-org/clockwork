#!/usr/bin/env python3

import copy
import os
import random
import pyfastaq

random.seed(42)
ref_seq = [random.choice(['A', 'C', 'G', 'T']) for x in range(1000)]
ref_seq2 = [random.choice(['A', 'C', 'G', 'T']) for x in range(1000)]
mutated_seq = copy.copy(ref_seq)
mutated_seq[250] = 'G' if mutated_seq[250] == 'A' else 'A'

ref = pyfastaq.sequences.Fasta('ref1', ''.join(ref_seq))
mutated = pyfastaq.sequences.Fasta('mutated', ''.join(mutated_seq))


with open('run.ref_to_make_reads.fa', 'w') as f:
    print(ref, file=f)
    print(mutated, file=f)

with open('run.ref.fa', 'w') as f:
    ref2 = pyfastaq.sequences.Fasta('ref2', ''.join(ref_seq2))
    print(ref, file=f)
    print(ref2, file=f)
