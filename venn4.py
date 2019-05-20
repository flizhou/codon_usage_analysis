#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  1 10:33:35 2019

@author: fanlizhou

4-way Venn diagram
Modified from ksahlin/pyinfor

"""

import io, os, argparse
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

def parse_args():    
    parser = argparse.ArgumentParser(description=
            'python3 venn4.py SP_top.txt SP_bot.txt LP_top.txt LP_bot.txt\n')
    parser.add_argument('sp_top', help = 'one input SP_top data file\n')
    parser.add_argument('sp_bot', help = 'one input SP_bot data file\n')
    parser.add_argument('lp_top', help = 'one input LP_top data file\n')
    parser.add_argument('lp_bot', help = 'one input LP_bot data file\n')
   
    args = parser.parse_args()
        
    for path in [args.sp_top, args.sp_bot, args.lp_top, args.lp_bot]:
        if not os.path.isfile(path):
            parser.error('File "%s" cannot be found.' % (path))
    
    return args


def read_data(filename):    
    file = io.open(filename)
    
    gene_set = set()
    
    for line in file:
        gene_set.add(line.strip())
    
    file.close()
    
    return gene_set


def venn4(data, groups):    
    fig, ax = plt.subplots(1, 1, figsize = (10, 10))
    colors = ['r', 'g', 'b', 'c']
    
    patches = []
    width, height = 170, 110
    
    # Make ellipse for each group
    patches.append(Ellipse((170, 170), width, height, - 45, 
                   color = colors[0], alpha = 0.5))
    patches.append(Ellipse((200, 200), width, height, - 45, 
                   color = colors[1], alpha = 0.5))
    patches.append(Ellipse((200, 200), width, height, - 135, 
                   color = colors[2], alpha = 0.5))
    patches.append(Ellipse((230, 170), width, height, - 135, 
                   color = colors[3], alpha = 0.5))
    
    for ellipse in patches:
        ax.add_patch(ellipse)
    
    ax.set_xlim(80, 320)
    ax.set_ylim(80, 320)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal")
    
    labels = set_labels(data)
    
    alignment = {'horizontalalignment' : 'center', 
                 'verticalalignment' : 'baseline'}
    
    # LP-top
    plt.text(120, 200, labels['1000'], **alignment)
    plt.text(280, 200, labels['0100'], **alignment)
    plt.text(155, 250, labels['0010'], **alignment)
    plt.text(245, 250, labels['0001'], **alignment)

    # LP-bot
    plt.text(200, 115, labels['1100'], **alignment)
    plt.text(140, 225, labels['1010'], **alignment)
    plt.text(145, 155, labels['1001'], **alignment)
    plt.text(255, 155, labels['0110'], **alignment)
    plt.text(260, 225, labels['0101'], **alignment)
    plt.text(200, 240, labels['0011'], **alignment)

    # SP-top
    plt.text(235, 205, labels['0111'], **alignment)
    plt.text(165, 205, labels['1011'], **alignment)
    plt.text(225, 135, labels['1101'], **alignment)
    plt.text(175, 135, labels['1110'], **alignment)

    # SP-bot
    plt.text(200, 175, labels['1111'], **alignment)

    # names of different groups
    plt.text(110, 110, groups[0], fontsize = 16, **alignment)
    plt.text(290, 110, groups[1], fontsize = 16, **alignment)
    plt.text(130, 275, groups[2], fontsize = 16, **alignment)
    plt.text(270, 275, groups[3], fontsize = 16, **alignment)

    leg = ax.legend(groups, loc = 'best', fancybox = True)
    leg.get_frame().set_alpha(0.5)
    
    plt.show()

  
def set_labels(data):    
    all_data = set()
    for gene_set in data:
        all_data |= gene_set
        
    results = {}
    
    # Calculate gene numbers in each area
    for i in range(1, 2 ** 4):
        key = bin(i).split('0b')[-1].zfill(4)
        value = all_data.copy()
        
        intersection = [data[i] for i in range(4) if key[i] == '1']
        difference = [data[i] for i in range(4) if key[i] == '0']
        
        for gene_set in intersection:
            value &= gene_set
        for gene_set in difference:
            value -= gene_set
            
        results[key] = value
            
    labels = {k : ("%d"  % (len(results[k]))) for k in results}
    
    return labels


# Main flow
args = parse_args()
sp_top = read_data(args.sp_top)
sp_bot = read_data(args.sp_bot)
lp_top = read_data(args.lp_top)
lp_bot = read_data(args.lp_bot)

venn4([sp_top, sp_bot, lp_top, lp_bot], 
      ["sp_top", "sp_bot", "lp_top", "lp_bot"])