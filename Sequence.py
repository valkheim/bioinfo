import os
import itertools

from dataclasses import dataclass
from collections import deque

import Base



class Sequence:
    def __init__(self, path):
        self.name = os.path.splitext(os.path.basename(path))[0]
        self.strand = ''
        self.length = 0
        self.bases = []
        self.distribution = {}
        self._parse_from_file(path)
        self.rna_strand = self.to_rna()

    def chunks(self, chunk_size=1):
        bases = iter(self.bases)
        while True:
            window = tuple(itertools.islice(bases, window_size))
            if not window:
                return

            yield window


    def walk_chunks(self, fptr, chunk_size=1):
        for chunk in self.chunks(chunk_size):
            fptr(chunk)


    def window(self, window_size=1):
        if window_size > self.length:
            window_size = self.length

        it = iter(self.bases)
        win = deque((next(it, None) for _ in range(window_size)), maxlen=window_size)
        yield win
        append = win.append
        for e in it:
            append(e)
            yield win


    def walk_window(self, fptr, window_size=1):
        for frame in self.window(window_size):
            fptr(frame)


    def to_rna(self):
        return self.strand.replace('T', 'U')




    def distance(self, sequence):
        """ Hamming distance """
        if len(self.strand) != len(sequence):
            raise ValueError("len mismatch")

        return sum(a != b for a, b in zip(self.strand, sequence))

    def _parse_from_file(self, path):
        idx = 0
        with open(path) as fh:
            while True:
                base = fh.read(1)
                if not base:
                    break

                base = base.upper()
                self.strand += base
                self.bases.append(Base.Base(base, base, idx))
                idx += 1
                try:
                    self.distribution[base]['total'] += 1
                except KeyError:
                    self.distribution[base] = {
                        'total': 1,
                        'amount': 0,
                        'frequency': 0
                    }

        self.length = idx
        for base in self.distribution.keys():
            self.distribution[base]['amount'] = self.distribution[base]['total'] / self.length
            self.distribution[base]['frequency'] = self.distribution[base]['amount'] / self.distribution[base]['total'] * 100

    def get_GC_AT_ratios(self):
        GC = self.distribution['G']['frequency'] + self.distribution['C']['frequency']
        AT = self.distribution['A']['frequency'] + self.distribution['T']['frequency']
        return (GC, AT)

    def __str__(self):
        dist = ' ; '.join([f"%{k} = {v['frequency']}" for k, v in self.distribution.items()])
        return  f"name:    {self.name}\n" \
                f"length:  {self.length}\n" \
                f"dist:    {dist}"
