#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2021-04-16 18:16

import argparse


def assembly_to_groups(assembly, len_cutoff):
    ctg_dict = dict()
    cluster_list = list()
    small_frag = set()
    with open(assembly) as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.split()
            if line.startswith('>'):
                ctg_dict[cols[1]] = cols[0][1:]
                if int(cols[2]) < len_cutoff:
                    small_frag.add(cols[1])
            else:
                cluster_list.append([num.strip('-') for num in cols if num.strip('-') not in small_frag])
    return ctg_dict, cluster_list


def output_clusters(ctg_dict, cluster_list):
    with open('prunning.clusters.txt', 'w') as f:
        f.write('#Group\tnContigs\tContigs\n')
        ngroup = len(cluster_list)
        for n, nums in enumerate(cluster_list, 1):
            f.write('{0}g{1}\t{2}\t{3}\n'.format(ngroup, n, len(nums), ' '.join([ctg_dict[num] for num in nums])))


def output_counts(ctg_dict, counts):
    with open(counts) as fin, open('sub.'+counts, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                fout.write(line)
            else:
                cols = line.split()
                if cols[0] in ctg_dict.values():
                    fout.write(line)


def output_fasta(ctg_dict, fasta):
    output = False
    with open(fasta) as fin, open('sub.'+fasta, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                if line.split()[0][1:] in ctg_dict.values():
                    output = True
                    fout.write(line)
                else:
                    output = False
            elif output:
                fout.write(line)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('assembly', help='*.review.assembly (output file of juicebox manual grouping), used to generate new prunning.clusters.txt')
    parser.add_argument('--fasta', default=None, help='input fasta file of contigs, this parameter will remove contigs not in .review.assembly, optional')
    parser.add_argument('--counts', default=None, help='input prunning.counts_RE.txt, this parameter will remove contigs not in .review.assembly, optional')
    parser.add_argument('--len_cutoff', default=100, type=float, help='length cutoff, default: %(default)s Kbp')
    args = parser.parse_args()

    ctg_dict, cluster_list = assembly_to_groups(args.assembly, args.len_cutoff*1000)
    output_clusters(ctg_dict, cluster_list)
    if args.fasta:
        output_fasta(ctg_dict, args.fasta)
    if args.counts:
        output_counts(ctg_dict, args.counts)


if __name__ == '__main__':
    main()

