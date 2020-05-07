#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wed Mar 13 17:34:32 2019

@author: fanlizhou

Analyze codon usage of sequence from 'SP_gene_seq.txt' and 'LP_gene_seq.txt'
Plot heatmap of amino acid usage and codon usage
Plot codon usage in each gene for each amino acid. Genes were arranged so that
the gene expression of SP decrease from 0 to 50 (x-axis) and the gene expression
of LP increase from 51 to 100 (x-axis)

Usage: codon_usage.py [-h] [--label LABEL] sp_file lp_file 

Options:
--label          Define the label of out-put files. Default="top"
sp_file          Path to the SP data files
lp_file          Path to the LP data files

"""

import io, os, argparse, collections
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description=
            'Analyze codon usage of SP and LP\n')
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


# A Codon_Usage class to store codon usage information for each genotype
class Codon_Usage:
    
    def __init__(self, filename):   
        self.seq, self.gene_num = self.get_seq(filename)
        
        
    def get_seq(self, filename):      
        file = io.open(filename)
        # List of selected gene sequences, excluded genes that are non-triple
        all_seq = []
        gene_seq = ''
        count_all = 0
        count_non_triple = 0
    
        for line in file:
            # Read a gene information line
            if line[0]=='>':
                count_all += 1
                
                # If a gene has been read, then append it to all_seq if the
                # sequence is triple
                if gene_seq!='':                    
                    if len(gene_seq)%3:
                        count_non_triple += 1
                    else:
                        all_seq.append(gene_seq)
                        
                gene_seq = ''
                
            # Read a gene sequence line    
            else:
                gene_seq += line.strip()
                
    
        file.close()   
        print('%s:\n%d genes added\n%d are non-triple\n'%
                          (filename[:2],count_all, count_non_triple))
        
        return (all_seq, count_all - count_non_triple)
        

    def get_AA(self, codon):
        # Dict key: codon -> AA
        codon_map = {
        'TTT':'Phe', 'TTC':'Phe', 'TTA':'Leu',  'TTG':'Leu',
        'TCT':'Ser', 'TCC':'Ser', 'TCA':'Ser',  'TCG':'Ser',
        'TAT':'Tyr', 'TAC':'Tyr', 'TAA':'STOP', 'TAG':'STOP',
        'TGT':'Cys', 'TGC':'Cys', 'TGA':'STOP', 'TGG':'Trp',
        'CTT':'Leu', 'CTC':'Leu', 'CTA':'Leu',  'CTG':'Leu',
        'CCT':'Pro', 'CCC':'Pro', 'CCA':'Pro',  'CCG':'Pro',
        'CAT':'His', 'CAC':'His', 'CAA':'Gln',  'CAG':'Gln',
        'CGT':'Arg', 'CGC':'Arg', 'CGA':'Arg',  'CGG':'Arg',
        'ATT':'Ile', 'ATC':'Ile', 'ATA':'Ile',  'ATG':'Met',
        'ACT':'Thr', 'ACC':'Thr', 'ACA':'Thr',  'ACG':'Thr',
        'AAT':'Asn', 'AAC':'Asn', 'AAA':'Lys',  'AAG':'Lys',
        'AGT':'Ser', 'AGC':'Ser', 'AGA':'Arg',  'AGG':'Arg',
        'GTT':'Val', 'GTC':'Val', 'GTA':'Val',  'GTG':'Val',
        'GCT':'Ala', 'GCC':'Ala', 'GCA':'Ala',  'GCG':'Ala',
        'GAT':'Asp', 'GAC':'Asp', 'GAA':'Glu',  'GAG':'Glu',
        'GGT':'Gly', 'GGC':'Gly', 'GGA':'Gly',  'GGG':'Gly'}

        if codon in codon_map:
            return codon_map[codon] 
        else:
            return ''
    
        
    def get_usage_dict(self, seq):
        # usage_dict structure:
        # dict key: AA -> [
        #                        dict key: codon -> 
        #                                     [codon_count,
        #                                      codon_count/AA_count]
        #                        AA_count
        #                 ]  
        usage_dict = \
            collections.defaultdict(lambda: 
                                    [
                                     collections.defaultdict(
                                             lambda: [0, 0]), 
                                     0
                                    ])
        # Save AAs usage information
        for index in range(0, len(seq), 3):
            codon = seq[index:index+3]
            AA = self.get_AA(codon)
            if AA:
                # Count how many times the AA appears
                usage_dict[AA][1] += 1
                # Count how many times the codon is used
                usage_dict[AA][0][codon][0] += 1
        
        # Calculate the codon usage percentage for an AA
        for AA in usage_dict:
            for codon in usage_dict[AA][0]:
                usage_dict[AA][0][codon][1] = \
                    usage_dict[AA][0][codon][0]/usage_dict[AA][1]

        return usage_dict


    def get_AA_dict(self):      
        # AA_dict structure:
        # 1st dict key: AA -> 2nd dict key: codon -> a list of codon usage 
        #                                           percentage of each gene                                 
        AA_dict = \
            collections.defaultdict(
                    lambda:collections.defaultdict(list))
        
        # Dict key: AA -> codon list
        AA_map = {
                    'Phe':['TTT', 'TTC'],
                    'Leu':['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
                    'Ser':['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 
                    'Tyr':['TAT', 'TAC'], 
                    'STOP':['TAA', 'TAG', 'TGA'],
                    'Cys':['TGT', 'TGC'],  
                    'Trp':['TGG'],
                    'Pro':['CCT', 'CCC', 'CCA', 'CCG'],
                    'His':['CAT', 'CAC'], 
                    'Gln':['CAA', 'CAG'],
                    'Arg':['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
                    'Ile':['ATT', 'ATC', 'ATA'], 
                    'Met':['ATG'],
                    'Thr':['ACT', 'ACC', 'ACA', 'ACG'],
                    'Asn':['AAT', 'AAC'], 
                    'Lys':['AAA', 'AAG'],
                    'Val':['GTT', 'GTC', 'GTA', 'GTG'],
                    'Ala':['GCT', 'GCC', 'GCA', 'GCG'],
                    'Asp':['GAT', 'GAC'], 
                    'Glu':['GAA', 'GAG'],
                    'Gly':['GGT', 'GGC', 'GGA', 'GGG']
                  }
        
        # list of codon usage for each gene
        usage_dict_list = []
        
        # Get codon usage information for each gene
        for seq in self.seq:
            usage_dict_list.append(self.get_usage_dict(seq))
        
        # Get the list of codon usage percentage from each gene 
        for AA in list(AA_map.keys()):
            for codon in AA_map[AA]:
                # Get codon usage information from each gene
                for usage_dict in usage_dict_list:
                    # Append codon usage percentage in the gene
                    AA_dict[AA][codon].append(
                            usage_dict[AA][0][codon][1])
            
        return AA_dict  
    

def heatmap_SP_LP(sp_AA_dict, lp_AA_dict, label):    
    # List of Chi-Square test results
    AA_chisquare = []
    # AA plotting annotation information
    AA_text = []
    
    # List of student's t-test results
    codon_ttest = []
    # Codon plotting annotaion information
    codon_text = []
    
    i = 0
    j = 0
    # Number of genes analyzed
    count_all = 0
    # Number of genes that show significant results
    count_sig = 0
    
    for AA in list(sp_AA_dict.keys()):        
        # Mean values of codon usage for each AA
        sp_codon_mean = []
        lp_codon_mean = []     
        
        for codon in sp_AA_dict[AA]:
            # Calculate ttest results             
            p_val = stats.ttest_ind(sp_AA_dict[AA][codon],
                                    lp_AA_dict[AA][codon],
                                    equal_var = False)[1]
            
            # Display eight codons in a row
            if not i % 8:
                codon_ttest.append([])
                codon_text.append([])
            i += 1
            
            # Handle NULL values
            if np.isnan(p_val):
                codon_ttest[-1].append(0)
                codon_text[-1].append(codon + '\n NA')
            # Save ttest p-values and annotation information    
            else:                
                codon_ttest[-1].append(p_val)
                codon_text[-1].append(codon + '\n' + str(round(p_val, 2)))
                count_all += 1
                if p_val < 0.5:
                    count_sig += 1
            
            sp_codon_mean.append(np.mean(sp_AA_dict[AA][codon]))
            lp_codon_mean.append(np.mean(lp_AA_dict[AA][codon]))       
        
        # Get Chi-Square test results of each AA
        p_val = stats.chisquare(np.array([sp_codon_mean, lp_codon_mean]), 
                                axis = None)[1]
        
        # Display six AA in a row
        if not j % 6:
                AA_chisquare.append([])
                AA_text.append([])
        j += 1
        
        # Handle Null values
        if np.isnan(p_val):            
            AA_chisquare[-1].append(0)
            AA_text[-1].append(AA + '\n NA')
        # Save Chi-Square test p-values and annotation information
        else:                                
            AA_chisquare[-1].append(p_val)
            AA_text[-1].append(AA + '\n' + str(round(p_val, 2)))
    
    # Handle empty cells
    for n in range(j % 6, 6):
        AA_chisquare[-1].append(0)
        AA_text[-1].append('')
    
    # Get list of AAs that show significant difference between SP and LP groups
    AAs = choose_codons(codon_ttest, codon_text)    

    AA_chisquare = np.array(AA_chisquare)
    codon_ttest = np.array(codon_ttest)
    
    AA_text = np.array(AA_text)
    codon_text = np.array(codon_text)

    print('%d out of %d codon show significant usage difference \
          between SP and LP genes (p_value < 0.5)\n' % 
          (count_sig, count_all))
    plot_heatmap(AA_chisquare, AA_text, 'AAs_ChiSquare', label)
    plot_heatmap(codon_ttest, codon_text, 'Codons_ttest', label)
    
    return AAs


def plot_heatmap(data, text, cbarlabel, label):
    
    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (10, 5))

    im, cbar = heatmap(data, ax, 'YlGn', cbarlabel)
    
    annotate_heatmap(im, text)

    fig.tight_layout()
    plt.show
    plt.savefig(f'../results/{cbarlabel}_{label}.png')     
    
def heatmap(data, ax, cmap, cbarlabel):
    
    if not ax:
        ax = plt.gca()
        
    im = ax.imshow(data, cmap)
    
    cbar = ax.figure.colorbar(im, ax=ax)

    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    ax.set_xticklabels(range(data.shape[1]))
    ax.set_yticklabels(range(data.shape[0]))

    ax.tick_params(top=False, bottom=True,
                   labeltop=False, labelbottom=True)

    # Draw white space between squares
    for edge, spine in ax.spines.items():
        spine.set_visible(False)
    
    ax.set_xticks(np.arange(data.shape[1] + 1) - 0.5, minor = True)
    ax.set_yticks(np.arange(data.shape[0] + 1) - 0.5, minor = True)
    ax.grid(which = 'minor', color = 'w', linestyle = '-', linewidth = 3)
    ax.tick_params(which = 'minor', bottom = False, left = False)    
    cbar.ax.set_ylabel(cbarlabel, va = 'top')

    return im, cbar


def annotate_heatmap(im, text_label):
    textcolors = ['black','white']

    data = im.get_array()
    # Set threshold to decide color
    threshold = im.norm(data.max()) / 2
        
    kw = dict(horizontalalignment = 'center',
              verticalalignment = 'center')
    
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color = textcolors[im.norm(data[i,j]) > threshold])
            im.axes.text(j, i, text_label[i,j], **kw)


def choose_codons(ttest, text):           
    # Dict key: AA -> codon
    # Only contains AAs with only two codon choices    
    codon_map = {
        'TTT':'Phe', 'TTC':'Phe', 'TAT':'Tyr', 'TAC':'Tyr',
        'TGT':'Cys', 'TGC':'Cys', 'CAT':'His', 'CAC':'His', 
        'CAA':'Gln', 'CAG':'Gln', 'AAT':'Asn', 'AAC':'Asn', 
        'AAA':'Lys', 'AAG':'Lys', 'GAT':'Asp', 'GAC':'Asp', 
        'GAA':'Glu', 'GAG':'Glu'}    
    
    codon_dict = collections.defaultdict(list)
    for i in range(len(ttest)):
        for j in range(len(ttest[i])):
            if ttest[i][j] < 0.01:
                codon = text[i][j][:3]
                if codon in codon_map:
                    codon_dict[codon_map[codon]].append(codon)
        
    file = io.open('AAs_to_compare.txt', 'w')    
    file.write('Compare following AAs\n')
    # AAs that have only two codon choices and show significant 
    # codon usage difference between SP and LP groups
    AAs = []
    
    for AA in codon_dict.keys():
        AAs.append(AA)       
        if len(codon_dict[AA]) == 2:
            file.write('%s: %s, %s\n' % 
                       (AA, codon_dict[AA][0], codon_dict[AA][1]))
        else:
            file.write('%s: %s\n' % (AA, codon_dict[AA][0]))
            
    file.close()
    
    return AAs
        

def plot_SP_LP(sp_AA_dict, lp_AA_dict):
    # Plot each AA
    for AA in list(sp_AA_dict.keys()):        
        # List of codon usage information
        codon_data = []
        # List of codon names
        codons = []
        
        for codon in sp_AA_dict[AA]: 
            # LP group data is displayed from lowest expressed genes 
            # to highest expressed genes
            lp_AA_dict[AA][codon].reverse()
            
            codons.append(codon)            
            codon_data.append([])
            # Display SP group data first and then LP group data
            codon_data[-1].append(sp_AA_dict[AA][codon])            
            codon_data[-1].append(lp_AA_dict[AA][codon])
        
        # Plot usage curves    
        codon_usage_plot(codon_data, AA, codons)

    
def codon_usage_plot(data, AA, codons):
    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (15,5))
    
    for i in range(len(data)):
        # 0-50 shows SP group data
        x_sp = np.linspace(0, 50, len(data[i][0]))
        # 50-100 shows LP group data
        x_lp = np.linspace(50, 100, len(data[i][1]))
        
        ax.plot(x_sp, data[i][0], label = 'sp_' + codons[i])
        ax.plot(x_lp, data[i][1], label = 'lp_' + codons[i])
        ax.legend(loc = 1)
        ax.set_title(AA)

 
def plot_distribution(sp_dict, lp_dict, AA):
    fig, axes = plt.subplots(nrows = 2, ncols =1, figsize = (40, 20))

    for codon in sp_dict[AA]:
        x = np.arange(len(sp_dict[AA][codon]))
        sp_y = np.array(sp_dict[AA][codon])
        lp_y = np.array(lp_dict[AA][codon])
        
        axes[0].plot(x, sp_y)
        axes[1].plot(x, lp_y)
        
    plt.show


def get_skellam_distribution(sp_dict, lp_dict, AA):    
    sp_mu = {}
    lp_mu = {}
    codons = []
    
    # Get mean values
    for codon in sp_dict[AA]:
        codons.append(codon)
        sp_mu[codon] = np.mean(sp_dict[AA][codon])
        lp_mu[codon] = np.mean(lp_dict[AA][codon])
    
    skellam_plot(sp_mu[codons[0]], sp_mu[codons[1]], 'SP-' + AA)
    skellam_plot(lp_mu[codons[0]], lp_mu[codons[1]], 'LP-' + AA)


def skellam_plot(mu1, mu2, name):    
    print(mu1,' ', mu2, ' ', mu1-mu2, ' ', name)

    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (5, 5))        
    x = np.arange(stats.skellam.ppf(0.01, mu1, mu2), 
                  stats.skellam.ppf(0.99, mu1, mu2))
    ax.plot(x, stats.skellam.pmf(x, mu1, mu2), marker = 'o', label = name)
    ax.legend(loc = 1)
    
    plt.show
      
    
# Main Flow
args = parse_args()
sp_codon_usage = Codon_Usage(args.sp_file)
lp_codon_usage = Codon_Usage(args.lp_file)

sp_AA_dict = sp_codon_usage.get_AA_dict()    
lp_AA_dict = lp_codon_usage.get_AA_dict()

print("Analyzing SP and LP %s group data\n" % (args.label))
    
AAs = heatmap_SP_LP(sp_AA_dict, lp_AA_dict, args.label)
plot_SP_LP(sp_AA_dict, lp_AA_dict)

# Optional
# Get Skellam distributions of AAs that have only two codon choices 
# and show distictive usage between SP and LP
'''
sp_all_codon_usage = Codon_Usage('SP_all_gene_seq.txt')
lp_all_codon_usage = Codon_Usage('LP_all_gene_seq.txt')

sp_all_AA_dict = sp_all_codon_usage.get_AA_dict()    
lp_all_AA_dict = lp_all_codon_usage.get_AA_dict()

for AA in AAs:
    plot_distribution(sp_all_AA_dict, lp_all_AA_dict, AA)
    get_skellam_distribution(sp_all_AA_dict, lp_all_AA_dict, AA)
'''