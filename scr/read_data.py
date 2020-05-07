#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Feb 28 13:19:27 2019

@author: fanlizhou

Extract SP and LP data from 'TL_expr_gcrma_filtered_20181128.csv'.
Apply log transformation and get mean values for 'wtC','wtH','gcn2C', 
and 'gcn2H' R1_wtC1 and R3_wtH were not included
Write data into 'SP_log_mean.txt' and 'LP_log_mean.txt'

Usage: src/read_data.py data_file 

Options:

data_file  Path to the input data file

"""
import csv, io, os, math, argparse
import matplotlib.pyplot as plt


def parse_args():    
    parser = argparse.ArgumentParser(description = 
            'Read and write data\n')
    parser.add_argument('data_file', help = 'one input data file\n')

    args = parser.parse_args()

    if not os.path.isfile(args.data_file):
        parser.error('File "%s" cannot be found.\n' % (args.data_file))
    
    return args


def read_data(file_path):  
    # expression data
    # index 0 is gene name, expression data starts at index 1
    sp_data = []
    lp_data = []
    # column names    
    sample = []     
    
    with open(file_path, 'r') as csvfile:
        file_reader = csv.reader(csvfile, delimiter= ',', quotechar = '"' )
        first_line = file_reader.__next__()
        
        # save column names
        for n in range(0, len(first_line)):
            if((n % 9 in [4,5,6] and (n not in [4, 15]))):
                sample.append(first_line[n].split('_')[1])
        
        # save data, each line is data from one gene
        for line in file_reader:
            # save gene name
            sp_data.append([line[-1]])
            lp_data.append([line[-1]])
            for n in range(1, len(line)):
                
                # save SP data for each gene
                if((n % 9 in [4,5,6]) and (n not in [4, 15])):
                    sp_data[-1].append(float(line[n]))
                
                # save LP data for each gene
                elif(n % 9 in [7,8,0] and (n not in [7, 18])):
                    lp_data[-1].append(float(line[n]))

    return (sample, sp_data, lp_data)


def data_log_mean(data, sample):   
    log_mean = []
    
    for gene in data:
        # log of expression data
        exp = math.log10(gene[1])
        count = 1
        # gene name
        log_mean.append([gene[0]])
        
        for i in range(2, len(gene)):
            # if current column has the same genotype as previous column
            # sample index starts from 0, 
            # gene index starts from 1 for data, index 0 is gene name
            if(sample[i-2] == sample[i-1]):
                exp += math.log10(gene[i])
                count += 1
           
            # else, calculate the log mean value of previous genotype
            else:
                log_mean[-1].append(exp / count)
                count = 1
                exp = math.log10(gene[i])
                
        # calculate the log mean value of the last genotype
        log_mean[-1].append(exp / count)
        
    return log_mean


def data_plot(data, name):   
    fig, _ = plt.subplots(nrows=2, ncols=2, figsize=(10,10))
    sample = ['wtC', 'wtH', 'gcn2C', 'gcn2H']
    
    # histograms
    for i in range(1, len(data[0])):
        fig.axes[i-1].hist([gene[i] for gene in data], bins='auto')
        title = 'Histogram: ' + name + ' '+sample[i-1]
        fig.axes[i-1].set_title(title)
        fig.axes[i-1].set_xlabel('log expression')
        fig.axes[i-1].set_ylabel('# of genes')


def write_data(data, name):    
    file = io.open('../data/clean_data' + name + '_log_mean.txt', 'w')
    
    file.write('{0}{1}\n'.format(name,'_log_mean'))
    
    file.write('%-20s\t%-20s\t%-20s\t%-20s\t%s\n' %
               ('wtC', 'wtH', 'gcn2C', 'gcn2H', 'name'))
    for gene in data:
        file.write('%-20s\t%-20s\t%-20s\t%-20s\t%s\n' %
               (gene[1], gene[2], gene[3], gene[4], gene[0]))
        
    file.close()
    
    
# main flow
args = parse_args()
(sample, sp_data, lp_data) = read_data(args.data_file)

sp_log_mean = data_log_mean(sp_data, sample)
lp_log_mean = data_log_mean(lp_data, sample)

data_plot(sp_log_mean, 'SP')
data_plot(lp_log_mean, 'LP')

write_data(sp_log_mean, 'SP')
write_data(lp_log_mean, 'LP')

