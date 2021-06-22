#!/usr/bin/env python
import sys
import os
import argparse
import multiprocessing
import pysam


def get_opt():
	group = argparse.ArgumentParser()
	group.add_argument('-r', '--ref', help='reference contig level assembly', required=True)
	group.add_argument('-g', '--alleletable', help='Allele.ctg.table', required=True)
	group.add_argument('-b', '--bam', help='bam file, default: prunning.bam', default='prunning.bam')
	group.add_argument('-d', '--workdir', help='work directory, default: wrk_dir', default='wrk_dir')
	group.add_argument('-t', '--thread', help='threads, default: 10', type=int, default=10)
	return group.parse_args()


def read_fasta(in_fa):
	fa_db = {}
	with open(in_fa, 'r') as fin:
		for line in fin:
			if line[0] == '>':
				id = line.strip().split()[0][1:]
				fa_db[id] = []
			else:
				fa_db[id].append(line.strip())
	for id in fa_db:
		fa_db[id] = ''.join(fa_db[id])
	
	return fa_db


def load_allele(allele_table):
	ctg_on_chr = {}
	chr_contain_ctg = {}
	with open(allele_table, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			chrn = data[0]
			if chrn.startswith('tig') or chrn.startswith('scaffold') or chrn.startswith('utg') or chrn.startswith('ctg'):
				continue
			for ctg in data[2:]:
				if ctg not in ctg_on_chr:
					ctg_on_chr[ctg] = {}
				if chrn not in ctg_on_chr[ctg]:
					ctg_on_chr[ctg][chrn] = 0
				ctg_on_chr[ctg][chrn] += 1
	for ctg in ctg_on_chr:
		max_chr = ""
		max_cnt = 0
		for chrn in ctg_on_chr[ctg]:
			if ctg_on_chr[ctg][chrn] > max_cnt:
				max_cnt = ctg_on_chr[ctg][chrn]
				max_chr = chrn
		ctg_on_chr[ctg] = max_chr
		if max_chr not in chr_contain_ctg:
			chr_contain_ctg[max_chr] = {}
		chr_contain_ctg[max_chr][ctg] = 1
	return ctg_on_chr, chr_contain_ctg


def split_files(chrn, allele_table, ref, bam_file, wrkdir):
	wrk_dir = os.path.join(wrkdir, chrn)
	if not os.path.exists(wrk_dir):
		os.mkdir(wrk_dir)
	
	print("\tDealing %s"%chrn)
	ctg_on_chr, chr_contain_ctg = load_allele(allele_table)
	fa_db = read_fasta(ref)

	sub_bam = os.path.join(wrk_dir, chrn+'.bam')
	sub_fa = os.path.join(wrk_dir, chrn+'.fa')
	with open(sub_fa, 'w') as fout:
		for ctg in chr_contain_ctg[chrn]:
			fout.write(">%s\n%s\n"%(ctg, fa_db[ctg]))

	with pysam.AlignmentFile(bam_file, 'rb') as fin:
		with pysam.AlignmentFile(sub_bam, 'wb', template=fin) as fout:
			for ctg in chr_contain_ctg[chrn]:
				for line in fin.fetch(contig=ctg):
					if line.next_reference_name and line.next_reference_name in ctg_on_chr and ctg_on_chr[line.next_reference_name]==chrn:
						fout.write(line)
	

def partition_gmap(ref, allele_table, bam, wrkdir, threads):
	if not os.path.exists(wrkdir):
		os.mkdir(wrkdir)
	
	print("Getting groups")
	chrn_db = {}
	with open(allele_table, 'r') as fin:
		for line in fin:
			chrn_db[line.strip().split()[0]] = 1

	bai = bam+'.bai'
	if not os.path.exists(bai):
		print("BAI file not found, starting index...")
		ret = os.system('samtools index %s'%bam)
		if ret==0:
			print("Index success")
		else:
			print("Fatal: bam file must be sorted")
			sys.exit(-1)

	print("Splitting files")
	if len(chrn_db) < threads:
		threads = len(chrn_db)
	pool = multiprocessing.Pool(processes=threads)
	for chrn in chrn_db:
		pool.apply_async(split_files, (chrn, allele_table, ref, bam, wrkdir,))
	pool.close()
	pool.join()
	print("Notice: If you got errors of \"Length mismatch\" during allhic extract, it is normal because we split bam with the same header, it will not effect the result")
	print("Finished")


if __name__ == '__main__':
	opts = get_opt()
	ref = opts.ref
	allele_table = opts.alleletable
	bam = opts.bam
	wrkdir = opts.workdir
	threads = opts.thread
	partition_gmap(ref, allele_table, bam, wrkdir, threads)

