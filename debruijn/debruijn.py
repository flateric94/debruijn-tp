#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import textwrap
import matplotlib.pyplot as plt
matplotlib.use("Agg")
from collections import OrderedDict

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file (default contigs.fasta)")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as an image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    with open (fastq_file, "r") as f:
        for i in f:
            yield next(f).strip("\n")
            next(f)
            next(f)
    pass


def cut_kmer(read, kmer_size):
    for i,_ in enumerate(read[:len(read) - kmer_size+1]):
        yield read[i:i+kmer_size]
    pass


def build_kmer_dict(fastq_file, kmer_size):
    dict = {}
    for read in read_fastq(fastq_file):
        for k_mer in cut_kmer(read, kmer_size):
            if k_mer in dict:
                dict[k_mer] += 1
            else:
                dict[k_mer] = 1
    return dict
    pass


def build_graph(kmer_dict):
    digraph = nx.DiGraph()
    for k_mer, occurence in kmer_dict.items():
        digraph.add_edge(k_mer[:-1], k_mer[1:], weight = occurence)
    return digraph
    pass


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    for path in path_list:
        if delete_entry_node and delete_sink_node:
            # supprime tous les noeuds du chemin
            graph.remove_nodes_from(path)
        if delete_entry_node and not delete_sink_node:
            # supprime tous les noeuds du chemin, sauf le dernier
            graph.remove_nodes_from(path[:len(path)-1])
        if not delete_entry_node and delete_sink_node:
            # supprime tout les noeuds du chemin, sauf le premier
            graph.remove_nodes_from(path[1:len(path)])
        if not delete_entry_node and not delete_sink_node:
            # supprime tout les noeuds du chemin, sauf le premier et le dernier
            graph.remove_nodes_from(path[1:len(path)-1])
    return graph
    pass


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    ecart_type_weight = statistics.stdev(weight_avg_list)
    ecart_type_lenght = statistics.stdev(path_length)
    if (ecart_type_weight > 0):
        best_path = weight_avg_list.index(max(weight_avg_list))
        path_list.pop(best_path)
        remove_paths(graph, path_list, delete_entry_node, delete_sink_node)
    if (ecart_type_weight == 0):
        if (ecart_type_lenght > 0):
            best_path = path_length.index(max(path_length))
            path_list.pop(best_path)
            remove_paths(graph, path_list, delete_entry_node, delete_sink_node)
        if (ecart_type_lenght == 0):
            best_path = randint(0,len(path_length))
            path_list.pop(best_path)
            remove_paths(graph, path_list, delete_entry_node, delete_sink_node)
    return graph
    pass

def path_average_weight(graph, path):
    """Compute the weight of a path"""
    return statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])

def solve_bubble(graph, ancestor_node, descendant_node):
    path_list = []
    path_length = []
    weight_avg_list = [] 
    for path in nx.all_simple_paths(graph, ancestor_node, descendant_node):
        path_list.append(path)
        path_length.append(len(path))
        weight_avg_list.append(path_average_weight(graph,path))
    graph = select_best_path(graph, path_list, path_length, weight_avg_list)
    return graph
    pass

def simplify_bubbles(graph):
    pass

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    roots = [v for v, d in graph.in_degree() if d == 0]
    return roots
    pass

def get_sink_nodes(graph):
    leaves = [v for v, d in graph.out_degree() if d == 0]
    return leaves
    pass

def get_contigs(graph, starting_nodes, ending_nodes):
    result = []
    for i in starting_nodes:
        for e in ending_nodes:
            for paths in nx.all_simple_paths(graph, i, e):
                contig=""
                contig += paths[0][0]
                if nx.has_path(graph,i,e):
                    for node in paths:
                        contig += node[-1]     
                    result.append([contig, len(contig)])
    return result
    pass

def save_contigs(contigs_list, output_file):
    with open(output_file, "w") as f:
        j=0
        for i in contigs_list:
            f.write(f">contig_{j} len={i[1]}\n")
            f.write(f"{textwrap.fill(i[0],width=80)}\n")
            j+=1
    pass


def draw_graph(graph, graphimg_file):
    """Draw the graph
    """                                    
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # graph = nx.DiGraph()
    # graph.add_edges_from([(("AG", "TC"), ("CA", "GT")), (("AC", "TG"), ("CA", "GT")), (("CA", "GT"), ("AG", "TC")),
    #     (("AG", "TC"), ("CG", "GC")), (("CG", "GC"), ("CG", "GC")), (("CG", "GC"), ("CT", "GA")), (("CT", "GA"), ("AT", "TC")),
    #     (("CT", "GA"), ("AA", "TT"))])
    # contig_list = get_contigs(graph, ["TC", "AC"], ["AT" , "AA"])
    # print(contig_list)
    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)


if __name__ == '__main__':
    main()
