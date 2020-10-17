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

__author__ = "Opale Rambaud"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Opale Rambaud"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Opale Rambaud"
__email__ = "opale.rambaud@gmail.com"
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
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()

def read_fastq(fastq_file):

    """Read a fastq file and return by yield sequences
    Parameters : fastq_file : a fastq_file given in argument
    Returns : sequences (yield)  
    """

    with open(fastq_file, "r") as filin:
        for line in filin:
            yield next(filin).strip()
            next(filin)
            next(filin)
            
def cut_kmer(seq, kmer_size):
    """
    Read all sequences and return kmers
    Parameters : seq : a sequence to cut
    		kmer_size : a given size of kmers in argument
    Returns : kmers
    """
    for i in range(len(seq)-kmer_size+1):
        yield seq[i:i+kmer_size]            

    
def build_kmer_dict(fastq_file, kmer_size):
    """
    Build a dictionnary of kmers for all sequences of the fastq file and the occurence of each kmers
    Parameters : fastq_file : a fastq_file given in argument
    		kmer_size: a given size of kmers in argument
    Returns: dictionnary of kmers with their occurences 
    """
    l_seq = []
    dic_kmer = {}
    for i in read_fastq(fastq_file):
        l_seq.append(i)
    for j in range(len(l_seq)):
        l_kmer = []
        for i in cut_kmer(l_seq[j], kmer_size):
            l_kmer.append(i)
        for i in range(len(l_kmer)):
            if l_kmer[i] not in dic_kmer:
                dic_kmer[l_kmer[i]] = l_kmer.count(l_kmer[i])
                
    return dic_kmer        
    
def build_graph(dic_kmer):
    """
    Build a tree of prefixes and suffixes kmers
    Parameters : dic_kmer : dictionnary of kmers with their occurences
    Returns : tree of prefixes and suffixes kmers
    """
    tree_graph = nx.DiGraph()
    for key, val in dic_kmer.items():
        tree_graph.add_edge(key[:-1], key[1:], weight=val)
        
    return tree_graph
    
        
def get_starting_nodes(tree_graph):
    """
    Construct a list of starting nodes from the tree_graph
    Parameters : tree_graph : tree of prefixes and suffixes kmers
    Returns : list of starting nodes 
    """
    start = []
    for node in tree_graph.nodes:
        if len(list(tree_graph.predecessors(node))) == 0:
            start.append(node)

    return start


def get_sink_nodes(tree_graph):
    """
    Construct a list of sink nodes from the tree_graph
    Parameters : tree_graph : tree of prefixes and suffixes kmers
    Returns : list of sink nodes 
    """
    sink = []
    for node in tree_graph.nodes:
        if len(list(tree_graph.successors(node))) == 0:
            sink.append(node)

    return sink
    
    
def get_contigs(tree_graph, start, sink):
    """
    Get the different possible contigs from list of starting nodes, list of sink nodes and tree_graph
    Paramters : tree_graph :tree of prefixes and suffixes kmers
    		start: list of starting nodes
    		sink: list of sink nodes 
    Returns : A list of tuples with contig and length of contig
    """  
    contigs = []
    for start_node in start:
        for sink_node in sink:
            for path in nx.all_simple_paths(tree_graph, start_node, sink_node):
                contig = path[0]
                for i in range(1,len(path)):
                    contig += path[i][-1]
                contigs.append((contig, len(contig)))
                
    return contigs    
def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def save_contigs(contig_tuple, output_file):
    """
    Save the contigs in a text file with fasta form
    Parameters : contig_tuple :tuples with contig and length of contig
    		output_file : the name of the output file given in arguments
    Returns : a text file with fasta form with contigs inside
    """
    with open(output_file, "w+") as filout: 
        for i, contig in enumerate(contig_tuple):
            filout.write(">contig_"+str(i)+" len="+str(contig[1])+"\n")
            filout.write(fill(contig[0])+"\n")
        filout.close()
    

def std(val):
    """
    Calcul std
    Parameters : val: liste of integer
    Returns: std of the list of integer
    """
    return statistics.stdev(val)

def path_average_weight(tree_graphe, path):
    """
    Caclul the average weights for all paths 
    Parameters : tree_graph :tree of prefixes and suffixes kmers
    		paths : list of paths
    Returns : Mean of the weight for each paths
    """
    avg_weights = []
    for edge in tree_graphe.subgraph(path).edges(data=True):
        avg_weights.append(edge['weight'])
    return statistics.mean(avg_weights)

def remove_paths(tree_graphe, paths, delete_entry_node=False, delete_sink_node=False):
    """
    Remove a path containing starting or ending nodes
    Parameters : tree_graph :tree of prefixes and suffixes kmers
    		paths : list of paths
    Returns : tree_graphe without paths containings starting or ending nodes
    """
    for path in paths:
        if delete_entry_node and delete_sink_node is not True:
            tree_graphe.remove_nodes_from(path[:-1])
        elif delete_entry_node and delete_sink_node:
            tree_graphe.remove_nodes_from(path)
        elif delete_entry_node is not True and delete_sink_node:
            tree_graphe.remove_nodes_from(path[1:])            
        elif delete_entry_node is not True and delete_sink_node is not True:
            tree_graphe.remove_nodes_from(path[1:-1])
          
    return tree_graphe

def select_best_path(tree_graph, paths, path_len_m, path_weight_m, delete_entry_node=False, delete_sink_node=False):
    """
    Select the best path between all paths
    Parameters : tree_graph :tree of prefixes and suffixes kmers
    		paths : list of paths to remove
    		path_len_m: list of paths lenght
    		path_weight_m: list of paths weights
    Returns : graph with only best path
    """
    random.seed(9001) 
      
    #weight :
   
    max_weight = max(path_weight_m)   
    weight_ind = [i for i, j in enumerate(path_weight_m) if j == max_weight]
    best_paths = [paths[i] for i in weight_ind]
    
    #lenght :
    
    path_len_m = [path_len_m[i] for i in idx]
    max_len = max(path_len_m)
    len_ind = [i for i, j in enumerate(path_len_m) if j == max_len]
    best_paths = [best_paths[i] for i in len_ind]
    
    
    if len(best_paths) > 1:
        best_paths = best_paths[random.randint(0, len(best_paths))]
    paths.remove(best_paths[0])
    tree_graph = remove_paths(tree_graph, paths, delete_entry_node, delete_sink_node)
    
    return tree_graph
    
def solve_bubble(tree_graph, ancestor_node, new_node):
    """
    Remove one bubble in tree_graph
    Parameters :tree_graph :tree of prefixes and suffixes kmers
    		ancestor_node 
    		new_node
    """
    paths_l = list(nx.all_simple_paths(tree_graph, ancestor_node, new_node))
    path_len = []
    path_weight = []
    for path in paths_l:
        path_len.append(len(path))
        path_weight.append(path_average_weight(tree_graph, path))
    graph = select_best_path(tree_graph, paths, path_len, path_weight)
    
    return tree_graph
    
    
def simplify_bubbles(tree_graphe):
    """
    Remove bubbles in tree_graph
    Parameters : tree_graphe: tree_graph :tree of prefixes and suffixes kmers
    Returns: tree_graph without bubbles
    """  
    bubbles = []
    for node in tree_graphe.nodes():
        predecessors_list = list(tree_graphe.predecessors(node))
        if len(predecessors_list) > 1:
            ancestor = nx.lowest_common_ancestor(tree_graphe, predecessors_list[0], predecessors_list[1])
            bubbles.append([ancestor, node])
    for i in range(len(bubbles)):
        tree_graphe = solve_bubble(tree_graphe, bubbles[i][0], bubbles[i][1])
        
    return tree_graphe      

#==============================================================
# Main program
#==============================================================
def main():

    """Main function"""

    #Lecture du fichier et construction du graphe :
    args = get_arguments()
    dic_kmer = build_kmer_dict(args.fastq_file, args.kmer_size)
    tree_graph = build_graph(dic_kmer)

    #Résolution des bulles : 
    tree_graph = simplify_bubbles(tree_graph)    

    #Ecriture des contigs : 
    strt_nodes = get_starting_nodes(tree_graph)
    s_nodes = get_sink_nodes(tree_graph)
    contigs = get_contigs(tree_graph, strt_nodes, s_nodes)
    save_contigs(contigs, args.output_file)
    

if __name__ == '__main__':
    main()
