#!/usr/bin/env python
import sys

import Sequence

if __name__ == "__main__":
    seq = Sequence.Sequence(sys.argv[1])
    print(seq.to_rna())
