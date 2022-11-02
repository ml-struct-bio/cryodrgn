"""Plot the learning curve"""

import argparse

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from cryodrgn import analysis


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", nargs="+", help="Input run.log file(s)")
    parser.add_argument("-o", help="Output PNG")
    return parser.parse_args()


def main(args):
    cmap = matplotlib.cm.get_cmap("jet")
    i = 0
    cs = np.arange(len(args.input)) / len(args.input)
    for f in args.input:
        loss = analysis.parse_loss(f)
        c = cmap(cs[i])
        plt.plot(loss, label=f, c=c)
        print(f)
        print(loss)
        i += 1
    plt.xlabel("epoch")
    plt.ylabel("loss")
    plt.legend(loc="best")
    if args.o:
        plt.savefig(args.o)
    else:
        plt.show()


if __name__ == "__main__":
    main(parse_args())
