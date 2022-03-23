import sys
import math


def get_n0(n: int):
    n0 = int(math.log2(n) * 100)
    return n0


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: %s <N>' % sys.argv[0])
    try:
        n = int(sys.argv[-1])
    except:
        print('A positive integer should be provided.', file=sys.stderr)
        exit(-1)
    n0 = get_n0(n)
    print('%d' % (n0))
