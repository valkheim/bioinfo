import re
import os
import itertools

from typing import List
from dataclasses import dataclass
from collections import deque

from . import Base
from . import Data


class Sequence:
    def __init__(self, strand_str, name='sequence'):
        self.name = name
        self._init_strand(strand_str)
        self.rna_strand = self.to_rna()


    def append(self, bases: str):
        for base in bases:
            self.strand += base

        self._init_strand(self.strand)


    def chunks(self, chunk_size=1):
        bases = iter(self.bases)
        while True:
            window = tuple(itertools.islice(bases, chunk_size))
            if not window:
                return

            yield window


    def walk_chunks(self, fptr, chunk_size=1):
        for chunk in self.chunks(chunk_size):
            fptr(chunk)


    def walk_bases(self, fptr):
        for base in self.bases:
            fptr(base)


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
            yield fptr(frame)

    def get_complement_base(self, base):
        # Adenine - Thymine
        # Cytosine - Guanin
        if base.letter == 'A': return 'T'
        if base.letter == 'T': return 'A'
        if base.letter == 'C': return 'G'
        if base.letter == 'G': return 'C'


    def complement(self):
        complement_str = ''
        for base in self.bases:
            complement_str += self.get_complement_base(base)
        return Sequence(complement_str)

    def invert(self):
        return Sequence(self.strand[::-1])


    def coding_regions(self):  # aka coding sequence (CDS)
        coding_regions = []
        coding_region = ''
        coding = False
        for bases in self.chunks(Data.CODON_LENGTH):
            codon = ''.join([ base.letter for base in bases ])
            if codon in Data.STOP:
                if coding_region:
                    coding_region += codon
                    coding_region_length = len(coding_region)
                    if coding_region_length > 6:  # start + stop
                        coding_regions.append(coding_region)

                coding_region = ''
                coding = False

            if codon == Data.START[0]:
                coding_region = Data.START[0]
                coding = True

            if not coding:
                continue

            if codon == Data.START[0] or codon in Data.STOP:
                continue

            coding_region += codon


        return coding_regions


    def to_rna(self):
        return self.strand.replace('T', 'U')


    def distance(self, sequence):
        """ Hamming distance """
        if len(self.strand) != len(sequence):
            raise ValueError("len mismatch")

        return sum(a != b for a, b in zip(self.strand, sequence))


    def walk_coding_regions(self):
        bits = re.split('|'.join(Data.STOP), self.strand)
        bits = []
        for bit in re.split('|'.join(Data.STOP), self.strand):
            candidate = re.search(f'^{Data.START}', bit)
            if not candidate:
                continue

            print(candidate.string)
            length = len(candidate.string)
            if length <= 3 or length % 3 != 0:
                continue

            bits.append(bit[3:])

        return bits


    def _init_strand(self, strand_str):
        idx = 0
        self.strand = ''
        self.strand = ''
        self.length = 0
        self.bases = []
        self.distribution = {}
        for base in strand_str:
            base = base.upper()
            self.strand += base
            self.bases.append(Data.Base(base, base, idx))
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
