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

    """Read all sequences and return kmers
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
    for i in read_fastq(fastq_file):
        l_seq.append(i)

    dic_kmer = {}
    for j in range(len(l_seq)):
        l_kmer = []
        for i in cut_kmer(l_seq[j], kmer_size):
            l_kmer.append(i)

        for i in range(len(l_kmer)):
            if l_kmer[i] not in dic_kmer:
                dic_kmer[l_kmer[i]] = l_kmer.count(l_kmer[i])
    return dic_kmer        
    

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

if __name__ == '__main__':
    main()
