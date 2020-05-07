#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wed Mar 13 16:08:13 2019

@author: fanlizhou

Read gene coding sequence from 'Araport11_genes.201606.txt'
Find gene sequence for all genes in 'SP_gene_dist.txt' and 'LP_gene_dist.txt'
Write sequence data into 'SP_seq.txt' and 'LP_seq.txt'

Usage: get_gene_seq.py [-h] [--label LABEL] sp_file lp_file 

Options:
--label          Define the label of out-put files. Default="top"
sp_file          Path to the SP data files
lp_file          Path to the LP data files

"""

import io, os, argparse, collections

def parse_args():    
    parser = argparse.ArgumentParser(description = 
            'Find gene sequences for SP and LP\n')
    parser.add_argument('sp_file', help = 'one input SP data file\n')
    parser.add_argument('lp_file', help = 'one input LP data file\n')
    parser.add_argument('--label', '-l', 
                        type = str, required = False, default = 'top', 
                        help = 'Define the label of out-put files. Default="top"\n')

    args = parser.parse_args()
    
    for path in [args.sp_file, args.lp_file]:
        if not os.path.isfile(path):
            parser.error('File "%s" cannot be found.' % (path))
 
    return args


# The Gene_Bank class stores all gene sequences from the input genome
class Gene_Bank:
    
    def __init__(self, filename):
        # gene_bank structure:
        # 1st dict key: chrom -> 2nd dict key: gene location//1000 
        #                        -> 3rd dict key: gene_name
        #                           -> [(gene info, gene sequence)] 
        self.gene_bank = collections.defaultdict(
                                lambda:collections.defaultdict(
                                      lambda:collections.defaultdict(
                                             list)))
        
        self.build_gene_bank(filename)
        
        
    def build_gene_bank(self, filename):
        file = io.open(filename)
        seq = ''
        name = ''
        gene_info = ''
        count_all = 0
        count_non_triple = 0
        
        for line in file:
            # Read a gene inofrmation line
            if line[0] == '>':
                # If has read a gene already, save the gene sequence
                if name != '':
                    # Save chromosome name as 1st bucket key
                    chrom = name[1 : 5]
                    # Save gene location // 1000 as 2nd bucket key
                    loc = int(name[5 : ]) // 1000
                    # Save gene sequence
                    self.gene_bank[chrom][loc][name[1 : ]].append((gene_info, seq))
                    
                    count_all += 1
                    if len(seq) % 3:
                        count_non_triple += 1
                
                # Read the gene information
                name = line.split()[0].split('.')[0]
                gene_info = line.strip()
                seq = ''
            
            # Read a sequence line
            else:
                seq += line.strip()
        
        print('%d genes added\n%d are non-triple\n' %
                          (count_all, count_non_triple))
        
    # Search a gene by chromsome and location
    # Return gene sequence     
    def search_genes(self, name):
        chrom = name[:4]
        loc = int(name[4:]) // 1000
        return self.gene_bank[chrom][loc][name]
                          

def read_data(filename):
    file = io.open(filename)
    # Target gene names
    genes = []
    
    for line in file:
        genes.append(line.strip())
    
    file.close()
    return genes
    

def get_sequence(gene_bank, genes):
    # List of gene sequences
    gene_seq = []
    
    for gene in genes:
        gene_seq.extend(gene_bank.search_genes(gene))
    
    return gene_seq


def write_data(gene_seq, name, label):

    file = io.open(f'../results/{name}_{label}_seq.txt', 'w')
    
    for gene in gene_seq:
        # Write gene name to file
        file.write('%s\n' % (gene[0]))
        
        # Write gene sequence to file, 70 AAs per line
        n = 70
        while n < len(gene[1]):
            file.write('%s\n' % (gene[1][n-70 : n]))
            n += 70
        file.write('%s\n' % (gene[1][n-70 : ]))
    
    file.close()
  

# Main Flow
args = parse_args()       
 
gene_bank = Gene_Bank('Araport11_genes.201606.txt')
sp_genes = read_data(args.sp_file)          
lp_genes = read_data(args.lp_file)          

sp_gene_seq = get_sequence(gene_bank, sp_genes)             
write_data(sp_gene_seq, 'SP', args.label)                

lp_gene_seq = get_sequence(gene_bank, lp_genes)             
write_data(lp_gene_seq, 'LP', args.label)   

               