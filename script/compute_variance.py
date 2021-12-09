import sys
import math
import argparse


def read_record(input: str, input_format: str, n: int):
    data = []
    with open(input) as f:
        lines = f.readlines()
    if n and len(lines) <= n:
        raise ValueError('the length of the file (%d) is less than n (%d)' % (len(lines), n))
    for i in range(n):
        if input_format == 'simple':
            p = float(lines[i])
        else:
            p = float(lines[i].split()[1])
    return data


def compute_var(data):
    mean = sum(data) / len(data)
    var = sum(map(lambda x: x ** 2, map(lambda x: x - mean, data))) / (len(data) - 1)
    return var


def main():
    parser = argparse.ArgumentParser(description='Computing variance of a series.')
    parser.add_argument('--input-format', type=str, choices=['simple', 'multirecord'], default='multirecord')
    parser.add_argument('--input', type=str, required=True)
    parser.add_argument('--result-type', type=str, choices=['var', 'std'], default='std')
    parser.add_argument('--n', type=int)
    args = parser.parse_args(sys.argv[1:])

    data = read_record(args.input, args.input_format, args.n)
    var = compute_var(data)
    if args.result_type == 'std':
        print(math.sqrt(var))
    else:
        print(var)


if __name__ == '__main__':
    main()
