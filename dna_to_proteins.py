#!/usr/bin/env python
import sys

import Sequence
import Proteins

if __name__ == "__main__":
    seq = Sequence.Sequence(sys.argv[1])
    print(Proteins.get_protein(seq))

