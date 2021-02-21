import sys
import numpy as np
import matplotlib.pyplot as plt
import argparse
import re

def parse_header(header_line):
    p = re.compile('\S+ \[\S*\]')
    header = []
    while (p.search(header_line)):
        begin, end = p.search(header_line).span()
        header.append(p.search(header_line).group())
        header_line = header_line[end:]
    return header

def read_data(filename, head=False):
    """ Reads data formated in one line header then data.
    The header gives the variables and units in : 
    variable [units] format """
    with open(filename, 'r') as fp:
        header = parse_header(fp.readline().strip())
        data = np.loadtxt(fp)
        if len(data[0, :]) != len(header):
            raise ValueError("Different values of header length and data length")
        indices_sort = np.argsort(data[:, 0])
        for j in range(len(header)):
            data[:, j] = np.array([data[i, j] for i in indices_sort])

    if head:
        return header, data
    else:
        return data

def plot_data(data, figname):
    fig, ax = plt.subplots()
    ax.plot(data[:, 0], data[:, 1], label='$y_1$')
    ax.plot(data[:, 0], data[:, 2], label='$y_2$')
    ax.plot(data[:, 0], data[:, 3], label='$y_3$')
    ax.legend()
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim([1e-10, 2])
    ax.grid(True)
    fig.savefig(figname, bbox_inches='tight')

if __name__ == '__main__':
    fig_dir = 'figures/'
    fn = 'cvRobert_dns.dat'
    data = read_data(fn)
    plot_data(data, f'{fig_dir}cvRobert_dns')