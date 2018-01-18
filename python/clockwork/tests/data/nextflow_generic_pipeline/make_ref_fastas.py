#!/usr/bin/env python3

import copy
import os
import random
import pyfastaq

ref_seq = [random.choice(['A', 'C', 'G', 'T']) for x in range(1000)]
mutated_seq = copy.copy(ref_seq)
mutated_seq[250] = 'G' if mutated_seq[250] == 'A' else 'A'
mutated_seq2 = copy.copy(mutated_seq)
mutated_seq.insert(750, 'A')
mutated_seq.pop(500)

ref = pyfastaq.sequences.Fasta('ref', ''.join(ref_seq))
mutated = pyfastaq.sequences.Fasta('mutated', ''.join(mutated_seq))
mutated2 = pyfastaq.sequences.Fasta('mutated2', ''.join(mutated_seq2))

for seq in ref, mutated, mutated2:
    with open(os.path.join('Reference.tmp', seq.id + '.fa'), 'w') as f:
        print(seq, file=f)
