#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Sat Mar  2 12:13:59 2019

@author: fanlizhou

Plot selected gene expression data against total expression data

Usage: data_plot.py [-h] [--name NAME] [--label LABEL] data_files [data_files ...] 

Options:
--name           SP(Default) or LP
--label          Define the label of out-put files. Default="top"
data_files       Path to the input data files

"""

import io, os, argparse
import matplotlib.pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser(description=
            'Plot histograms of gene expression data\n')
    parser.add_argument('data_files', nargs = '+',
                        help = 'log_mean.txt gene.txt')
    parser.add_argument('--name', '-n', required = False, default = 'SP',
                        help = 'SP(Default) or LP')
    parser.add_argument('--label', '-l', 
                        type = str, required = False, default = 'top', 
                        help = 'Define the label of out-put files. Default="top"\n')
    
    args = parser.parse_args()
    
    paths = list(args.data_files)
    
    for path in paths:
        if not os.path.isfile(path):
            parser.error('File "%s" cannot be found.' % (path))
    return args


def read_data(filename):
    
    data = []
    file = io.open(filename)
    file.readline()
    file.readline()
    
    for line in file:
        
        split_line = line.split('\t')
        data.append([split_line[-1]])
        for i in range(0, len(split_line)-1):
            data[-1].append(float(split_line[i]))
            
    file.close()
    return data

def data_plot(data_all, data_sel, name, label):
    
    fig, _ = plt.subplots(nrows = 2, ncols = 2, figsize = (10, 10))
    sample = ['wtC','wtH','gcn2C','gcn2H']
    
    for i in range(1, 5):
        # plot complete expression data
        fig.axes[i-1].hist([gene[i] for gene in data_all], bins = 'auto')
        # plot selected gene expression data
        fig.axes[i-1].hist([gene[i] for gene in data_sel], 
                            bins = 'auto', color = 'red')
        
        title = 'Histogram: ' + name + ' ' + label + ' ' + sample[i-1]
        fig.axes[i-1].set_title(title)
        fig.axes[i-1].set_xlabel('log expression')
        fig.axes[i-1].set_ylabel('# of genes')

    plt.savefig(f'../results/{name}_{label}.png')           

# main flow
args = parse_args()
data_all = read_data(args.data_files[0])
data_sel = read_data(args.data_files[1])
data_plot(data_all, data_sel, args.name, args.label)
        