import sys
import math


def get_n0n1(n: int, policy: str):
    if policy == 'policy1':
        n0 = math.ceil(math.sqrt(n))
        n1 = math.ceil(math.log2(n))
    elif policy == 'policy2':
        n0 = min(n, max(2 * math.sqrt(n), 1024))
        n1 = min(max(1, math.ceil(n / n0)), 5 + math.ceil(math.log2(n)))
    else:
        print('Invalid argument policy: %s', policy)
        exit(-1)
    return n0, n1


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: %s <N>' % sys.argv[0])
    try:
        n = int(sys.argv[-1])
    except:
        print('A positive integer should be provided.', file=sys.stderr)
        exit(-1)
    n0, n1 = get_n0n1(n, 'policy2')
    print('%d %d' % (n0, n1))
