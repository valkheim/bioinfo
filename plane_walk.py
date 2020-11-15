#!/usr/bin/env python

import copy
import sys
import Sequence

import matplotlib.pyplot as plt

def plane_walk(seq, compression=1):
    """
    Creates a 2D representation of a sequence using the following coordinates system

          A

          |
     G ---+--- C
          |

          T
    """
    point = [0, 0]
    points = [ point ]
    fig = plt.figure()
    def get_next_point(bases):
        next_base = bases[-1].name
        if next_base == 'C':
            point[0] += 1
        elif next_base == 'G':
            point[0] -= 1

        if next_base == 'A':
            point[1] += 1
        elif next_base == 'T':
            point[1] -= 1

        points.append(copy.deepcopy(point))


    seq.walk_window(get_next_point, compression)
    for i in range(1, len(points)):
        a, b = points[i-1], points[i]
        plt.plot([a[0], b[0]], [a[1], b[1]], color='blue')

    figname = f'{seq.name}-{compression}.png'
    fig.savefig(figname)
    print(f'[+] plot written to {figname}')


if __name__ == "__main__":
    seq = Sequence.Sequence(sys.argv[1])
    plane_walk(seq, int(seq.length * .01) or 1)
