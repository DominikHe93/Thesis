#!/usr/bin/env python
import argparse
import sys
import os
import gzip
from Bio import SeqIO

version = "0.2.0 (12.06.17)"
name = "run_minimus2.py"
parser=argparse.ArgumentParser(description = "Will extract the 16S genes from genbank files and write them to a multi fasta file. version = " + version)
parser.add_argument('-1', '--in1', action = "store", dest = "input1", required = True, help = "reference_fasta file")
parser.add_argument('-2', '--in2', action = "store", dest = "input2", required = True, help = "subject_fasta file")
parser.add_argument('-o', '--out', action = "store", dest = "output", default = "minimus2_output", help = "Output-basename")
parser.add_argument('--minoverlap', action = "store", type = int, dest = "minoverlap", default = 500, help = "Minimum overlap for minimus2 (default=500)")
parser.add_argument('--minident', action = "store", dest = "minident", default = "94", help = "minimum identity for overlap")
parser.add_argument('--maxconserr', action = "store", dest = "maxconserr", default = "0.06", help = "Maximum consensus error (0..1) (Default 0.06)")
parser.add_argument('--maxtrim', action = "store", dest = "maxtrim", default = "20", help = "Maximum sequence trimming length (Default 20bp)")
parser.add_argument('--minpath', action = "store", dest = "minimus2_path", default = "", help = "Path to minimus2-executables") #ignoring this for now, just assuming minimus2 in PATH
parser.add_argument('-V', '--version', action="version", version=name + " version " + version) 
args = parser.parse_args()

def DELME_read_fasta(input_fasta):
	if input_fasta.endswith(".gz"):
		try:
			readfile = gzip.open(input_fasta, 'r')
		except Exception, ex:
			sys.stderr.write(ex.__class__.__name__ + " : " + str(ex) + "\n")
			sys.stderr.write("\nInput file is expected to be gzip-compressed based on file extension. But it isn't!\n")
			raise IOError
	else:
		readfile = open(input_fasta, "r")
	
	fasta_list = list(SeqIO.parse(readfile, "fasta"))
	readfile.close()
	
	read_dict = { fasta_list[x].id : x + 1 for x in range(0, len(fasta_list)) }
	
	return read_dict

def read_fasta(input_fasta):
	if input_fasta.endswith(".gz"):
		try:
			readfile = gzip.open(input_fasta, 'r')
		except Exception, ex:
			sys.stderr.write(ex.__class__.__name__ + " : " + str(ex) + "\n")
			sys.stderr.write("\nInput file is expected to be gzip-compressed based on file extension. But it isn't!\n")
			raise IOError
	else:
		readfile = open(input_fasta, "r")
	
	fasta_list = list(SeqIO.parse(readfile, "fasta"))
	readfile.close()
	
	return fasta_list

def combine_input(input1, input2, output_basename):
	input1_count = len(input1)
	input2_count = len(input2)
	merged_input = input1 + input2
	combined_inputname = output_basename + ".inputseq"
	sys.stderr.write("\nmerging input fastas to {}\n".format(combined_inputname))
	merged_file = SeqIO.write(merged_input, combined_inputname, "fasta")
	
	sys.stderr.write("\ncreating contig-name hashtable\n")
	scg_dict = { input1[x].id : x + 1 for x in range(0, input1_count) }
	mg_dict = { input2[x].id : x + 1 + input1_count for x in range(0, input2_count) }
	print "scg_dict = {} entries".format(len(scg_dict))
	print "mg_dict = {} entries".format(len(mg_dict))
	return scg_dict, mg_dict

def call_minimus2(refcount, output_basename):
	from subprocess import call
	toamos_cmd = "{cmd} -s {inname} -o {outname}".format(cmd = os.path.join(args.minimus2_path, "toAmos"), inname = output_basename + ".inputseq", outname = output_basename + ".afg")
	minimus2_cmd = "{cmd} {outname} -D REFCOUNT={refcount} -D OVERLAP={overlap} -D CONSERR{conserr} -D MINID={minid} -D MAXTRIM={maxtrim}".format(cmd = os.path.join(args.minimus2_path,"minimus2"), outname = output_basename, refcount = refcount, overlap = args.minoverlap, conserr = args.maxconserr, minid = args.minident, maxtrim = args.maxtrim)
	
	sys.stderr.write("\nrunning: {}\n".format(toamos_cmd))
	call(toamos_cmd.split())
	sys.stderr.write("\nrunning: {}\n".format(minimus2_cmd))
	call(minimus2_cmd.split())
	return output_basename + ".contig"

class merged_SCG(object):
	def __init__(self, line_tokens):
		self.index = line_tokens[0].lstrip('#')
		self.length = line_tokens[2]
		self.merged_SCG_names = []
		self.merged_MG_names = []
		self.merged_SCG_lengths = []
		self.merged_MG_lengths = []
		
	
	def merge(self, linetokens, scg_dict, mg_dict):
		import re
		ctoken = linetokens[0]
		#print "------"
		#print ctoken
		#print re.search(r"\(\d+\) ", ctoken).start()
		#print "-----"
		contigname = ctoken[:re.search(r"\(\d+\)$", ctoken).start()].lstrip('#')
		length = linetokens[2]
		if contigname in scg_dict:
			#print "-->{} in scg dict".format(contigname)
			self.merged_SCG_names.append(contigname)
			self.merged_SCG_lengths.append(length)
		else:
			#print "###-->{} in mg dict".format(contigname)
			self.merged_MG_names.append(contigname)
			self.merged_MG_lengths.append(length)
	
	def outnamestring(self):
		return "{name}\t{length}\t{scglengths}\t{mglengths}\t{scgnames}\t{mgnames}\n".format(name = self.index, length = self.length, scglengths = ";".join(self.merged_SCG_lengths), mglengths = ";".join(self.merged_MG_lengths), scgnames = ";".join(self.merged_SCG_names), mgnames = ";".join(self.merged_MG_names))
	
	def outindex_string(self, scg_dict, mg_dict):
		return "{name}\t{length}\t{scglengths}\t{mglengths}\t{scgnames}\t{mgnames}\n".format(name = self.index, length = self.length, scglengths = ";".join(self.merged_SCG_lengths), mglengths = ";".join(self.merged_MG_lengths), scgnames = ";".join([str(scg_dict[i]) for i in self.merged_SCG_names]), mgnames = ";".join([str(mg_dict[i]) for i in self.merged_MG_names]))

def read_minimus_file(filename, scg_dict, mg_dict):
	merged_list = []
	infile = open(filename, "r")
	for line in infile:
		if line.startswith("#"):
			if line.startswith("##"):
				merged_list.append(merged_SCG(line.split()))
				continue
			merged_list[-1].merge(line.split(), scg_dict, mg_dict)
	infile.close()
	return merged_list

def write_output(merged_list, ref_dict, subj_dict, output_basename):
	header = "merged_contig\tlength\tref-contig_lengths\tsubj-contig_lengths\tref-contigs\tsubjcontigs\n"
	sys.stderr.write("\nwriting contig names table\n")
	outnames = open(output_basename + "_merged_contignames.tab", "w")
	outnames.write(header)
	for m in merged_list:
		outnames.write(m.outnamestring())
	outnames.close()
	
	sys.stderr.write("\nwriting contig indices table\n")
	outindices = open(output_basename + "_merged_contigindices.tab", "w")
	outindices.write(header)
	for m in merged_list:
		outindices.write(m.outindex_string(ref_dict, subj_dict))
	outindices.close()

def main():
	reference = read_fasta(args.input1)
	subject = read_fasta(args.input2)
	ref_dict, subj_dict = combine_input(reference, subject, args.output)
	minimus2_result = call_minimus2(len(reference), args.output)
	merged_list = read_minimus_file(minimus2_result, ref_dict, subj_dict)
	write_output(merged_list, ref_dict, subj_dict, args.output)

main()

def test():
	minimus2_result = "premerge01_marine03.contig"
	print "reading reference"
	ref_dict = DELME_read_fasta(args.input1)
	print "reading subject"
	subj_dict = DELME_read_fasta(args.input2)
	print "merging list"
	merged_list = read_minimus_file(minimus2_result, ref_dict, subj_dict)
	write_output(merged_list, ref_dict, subj_dict, args.output)
#test()
