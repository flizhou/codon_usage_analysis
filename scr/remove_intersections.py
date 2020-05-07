#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wed Mar 13 15:31:40 2019

@author: fanlizhou

Remove gene intersections from 'SP_gene.txt' and 'LP_gene.txt'
Sort genes based on gene expression in descending order
Write gene names into 'SP_(label)_genes_no_overlaps.txt' 
and 'LP_(label)_genes_no_overlaps.txt'

Usage: remove_intersections.py [-h] [--label LABEL] sp_file lp_file 

Options:
--label          Define the label of out-put files. Default="top"
sp_file          Path to the SP data files
lp_file          Path to the LP data files
"""

import io, os, argparse

def parse_args():    
    parser = argparse.ArgumentParser(description = 
            'Remove intersections and sort genes\n')
    parser.add_argument('sp_file', help = 'one input SP data file\n')
    parser.add_argument('lp_file', help = 'one input LP data file\n')
    parser.add_argument('--label', '-l', 
                        type = str, required = False, default = 'top', 
                        help = 'Define the label of out-put files. Default="top"\n')

    args = parser.parse_args()
       
    for path in [args.sp_file, args.lp_file]:
        if not os.path.isfile(path):
            parser.error('File "{0}" cannot be found.' % (path))
    
    return args


def read_data(filename):    
    # Set of gene names
    gene_set = set()
    # Dictict key: gene name -> mean expression
    gene_dict = {}
    
    file = io.open(filename)
    file.readline()
    file.readline()
    
    for line in file:
        split_line = line.split()
        mean = 0.0
        for i in range(0, 4):
            mean += float(split_line[0]) / 4
        
        for i in range(4, len(split_line)):
            gene_set.add(split_line[i])
            gene_dict[split_line[i]]= mean
            
    file.close()
     
    return gene_set, gene_dict


def remove_intersections(sp_set, lp_set):    
    sp_len = len(sp_set)
    lp_len = len(lp_set)
    intersection = len(sp_set & lp_set)
    
    # Remove intersections from sp genes
    sp_genes = sp_set - lp_set
    # Remove intersections from lp genes
    lp_genes = lp_set - sp_set
    
    print('Groups\t\tNumber of genes\nSP\t\t%d\nLP\t\t%d\nIntersection\t%d\n'
          % (sp_len, lp_len, intersection))
    
    return (sp_genes, lp_genes)


# Helper function for sort_data
# Output is the index to insert target
def binary_insert(values, target):    
    left = 0
    right = len(values)
    
    while left < right:
        mid = left + (right - left) // 2
        if target > values[mid]:
            left = mid + 1
        else:
            right = mid
            
    return right


# Sort data based on mean expression in descending order
def sort_data(gene_set, gene_dict):    
    # List of gene names sorted by corresponding gene expression
    genes = []
    # List of gene expression sorted by gene expression
    exp_val = []
    
    for gene in gene_set:
        val = gene_dict[gene]
        # Find the index to insert based on gene expression
        index = binary_insert(exp_val, val)
        exp_val.insert(index, val)
        genes.insert(index, gene)
    
    return genes


def write_data(gene_set, gene_dict, name, label):         

    file = io.open(f'../results/{name}_{label}_gene_dist.txt', 'w')

    genes = sort_data(gene_set, gene_dict)

    # Write sorted gene names to file 
    for gene in genes:
        if(gene != 'NA'):            
            file.write(gene)            
            file.write('\n')
        
    file.close()
    
    
# Main Flow
args = parse_args()
sp_set, sp_dict = read_data(args.sp_file)
lp_set, lp_dict = read_data(args.lp_file)

sp_genes, lp_genes = remove_intersections(sp_set, lp_set)  
  
write_data(sp_genes, sp_dict, 'SP', args.label)  
write_data(lp_genes, lp_dict, 'LP', args.label) 
 