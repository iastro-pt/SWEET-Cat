#!/usr/bin/python

import argparse

def _parse():
    '''Calculate the surface gravity from M and R'''
    p = argparse.ArgumentParser(description='Remove new line from last line of file.')
    p.add_argument('filein', help='file to edit last line', type=str)
    return p.parse_args()


def remove_newline_last_line(filein):
    """Remove new line from last line of file."""
    f = open(filein, 'r')
    d = f.read()
    f.close()
    d = d[:-1]
    open(filein, 'w').write(d)
    return

def main():
    args = _parse()
    remove_newline_last_line(args.filein)

if __name__ == '__main__':
    main()
