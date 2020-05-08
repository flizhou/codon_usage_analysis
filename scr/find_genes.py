#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Feb 28 17:04:55 2019

@author: fanlizhou

Read data from 'SP_log_mean.txt' or 'LP_log_mean.txt'
Use quick partition to find top/bottom (cut_off) genes in
'wtC','wtH','gcn2C', or 'gcn2H'
Find genes intersections in 'wtC','wtH','gcn2C', and 'gcn2H'
Write those gene expression data into 'SP_gene.txt' or 'LP_gene.txt'

Usage: find_genes.py [-h] [--top TOP] [--cut_off CUT_OFF] data_file 

Options:
--top            Select top genes (1) or bottom genes (0). Default=1
--cut_off        Define the cut_off to select genes. Default=0.1
data_file        Path to the input data file

"""

import io, os, argparse


def parse_args():
    parser = argparse.ArgumentParser(description = 
            'Find top or bottom genes\n')
    parser.add_argument('data_file', help = 'one input data file\n')
    parser.add_argument('--top', '-t', 
                        type = int, required = False, default = 1, 
                        help = 'Select top genes (1) or bottom genes (0). Default=1\n')
    parser.add_argument('--cut_off', '-c', 
                        type = float, required = False, default = 0.1, 
                        help = 'Define the cut_off to select genes. Default=0.1\n')

    args = parser.parse_args()

    if not os.path.isfile(args.data_file):
        parser.error('File "%s" cannot be found.\n' % (args.data_file))
    
    return args


def read_data(filename):   
    # list of expression data lists for each genotype
    data = [[], [], [], []]
    # gene names
    genes = []
    # index for each gene
    index = []
    count = 0
    
    file = io.open(filename)
    file.readline()
    file.readline()
    
    for line in file:
        split_line = line.split('\t')
        
        for i in range(0, 4):
            data[i-1].append(float(split_line[i]))
        
        genes.append(split_line[-1])
        index.append(count)
        count += 1
    
    return (data, genes, index)


# helper function for quick_partition
# take a float list and use quick partition to find (gene_num) genes
# with highest/lowest expression
# data is not changed, only index is modified
def partition(data, index, start, size, gene_num, top):    
    left = start + 1
    right = start+size - 1
    pivot = data[index[start]]
    
    # if top is True, move numbers (>pivot) to the right, 
    # and move numbers (<pivot) to the left
    # if top is False, move numbers the other way around
    while(left < right):

        while(left != start + size and 
              ((top and data[index[left]] <= pivot) or 
               (not top and data[index[left]] >= pivot))):
            left += 1
            
        while(right != start and 
              ((top and data[index[right]] >= pivot) or
               (not top and data[index[right]] <= pivot))):
            right -= 1
         
        if(left < right and 
           ((top and data[index[left]] > data[index[right]]) or
            (not top and data[index[left]] < data[index[right]]))):
            index[right], index[left] = index[left], index[right]
            left += 1
            right -= 1
            
    # if top is True, get genes with highest expression levels       
    if(top and data[index[right]] < pivot):    
        index[start], index[right] = index[right], index[start]

    # if top is False, get genes with lowest expression levels       
    if (not top and data[index[right]] > pivot):    
        index[start], index[right] = index[right], index[start]
        
    # continue sorting left of the pivot if less than (gene_num) genes 
    # are in the right part of the pivot and the current left part is not empty
    if(len(data) - left < gene_num and right > start + 1):
        partition(data, index, start, right - start, gene_num, top)
    
    # continue sorting right of the pivot if more than (gene_num) genes 
    # are in the right part of the pivot and the current right part is not empty
    elif(len(data) - left > gene_num and left < start + size - 1):
        partition(data, index, left, start + size - left, gene_num, top)
    

def quick_partition(data, index, cut_off, top):
    # total number of genes to save
    gene_num = cut_off * len(data)
    
    partition(data, index, 0, len(data), gene_num, top)
    
    # target gene indice are save at the right end of index list
    return set(index[len(data) - int(gene_num) : len(data)]) 
         

def find_intersections(gene_sets):
    return gene_sets[0] & gene_sets[1]\
            & gene_sets[2] & gene_sets[3]


def write_data(data, genes, gene_index, name, top):
    top_or_bottom = '_top'
    if not top:
        top_or_bottom = '_bot'
        
    filename = name[:2] + top_or_bottom + '_gene.txt'
    file = io.open(filename, 'w')
    
    file.write('%s\n' % (filename))
    file.write('%-20s\t%-20s\t%-20s\t%-20s\t%s\n' %
               ('wtC','wtH','gcn2C','gcn2H','name'))
    
    for index in gene_index:
        file.write('%-20s\t%-20s\t%-20s\t%-20s\t%s' %
                   (data[0][index], data[1][index], 
                    data[2][index], data[3][index], genes[index]))
    
    file.close()

  
# main flow
args = parse_args()
(data, genes, index) = read_data(args.data_file)
gene_sets = []
for i in range(4):
    gene_sets.append(quick_partition(data[i], index.copy(), args.cut_off, args.top))

gene_index = find_intersections(gene_sets)    

write_data(data, genes, gene_index, args.data_file, args.top)
    
    
            
        
        