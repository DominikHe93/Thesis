#!/usr/bin/env python
import argparse
import os.path
import sys
import re
import gzip
import bz2

version = "0.1.2 (22.08.2018)"
name = os.path.basename(sys.argv[0]) #get scriptname from actual script filename that was called
parser=argparse.ArgumentParser(description = "Prepare BLAST table for contig-classification using KRONA-tools. \
For each protein of each contig, the Blast hits are filtered to remove all hits with a score lower than X% of the best score for that protein.\
Then protein names are truncated to/replaced with the corresponding contig names and all hits are sorted based on contig name and (descending) score.\
The Y% best scoring hits for each contig wil be kept and used for subsequent LCA-classification with KRONA-tools")
parser.add_argument('-i','--input', action = "store", dest = "input_table", required = True, help = "Input file must be tabular blast output ('outfmt 6'). Query proteins must be based on PRODIGAL predictions")
parser.add_argument('-o','--output', action = "store", dest = "output_table", required = True, help = "output file (for use as input in KRONA-TOOLS)")
parser.add_argument("-sc", "--score_cutoff", action = "store", dest = "score_cutoff", type = float, default = 0.5, help = "[FLOAT] fraction of the best observed bit-score of each contig, below which new hits will not be recorded. (set between 0.0-1.0, default = 0.5) set to 0 in order to consider ALL hits")
parser.add_argument("--maxhits", action = "store", dest = "maxhits", type = float, default = 0.1, help = "maximum number of hits to consider per contig. If set between 0-1: fraction of the total hits per contig; if set larger than 1: absolute count; default = 0.1 (10%% best hits). Set to 1 in order to consider ALL hits")
parser.add_argument("--lookup_file", action = "store", dest = "lookup_file", default = None, help = "If proteins are NOT named according to prodigal standards (<Contig-name>_<ORF-Nr>), provide either a .gff, .gbk or a .tab (tab seperated lookup_table) file linking each protein to a contig. If a .gbk or .gff is supplied, protein_names will be assumed to be locus_tags") #TODO make option for lookuptable
parser.add_argument('-V', '--version', action = "version", version = name + " version " + version)
args = parser.parse_args()

def openfile(infile):
	assert os.path.exists(infile), "Could not find file %s" %infile
	try:
		if infile.endswith(".gz") or infile.endswith(".gzip"):
			print "COMPRESSED"
			readfile = gzip.open(infile, 'r')
		elif infile.endswith(".bz2") or infile.endswith(".bzip2"):
			readfile = bz2.BZ2File(infile, 'r')
		else:
			readfile = open(infile, 'r')
	except Exception, ex:
		print ex.__class__.__name__ + " : " + str(ex)
		return None
	else:
		return readfile

def get_lookup_dict(lookup_file):
	sys.stderr.write("\nopening supplied lookup_file\n")
	for i in [".gb.gz", ".gb", ".gbk.gz", ".gbk", ".gbff.gz", ".gbff"]:
		if lookup_file.endswith(i):
			return lookup_from_gb(lookup_file)
	for i in [".gff", ".gff.gz"]:
		if lookup_file.endswith(i):
			return lookup_from_gff(lookup_file)
	for i in [".tab", ".tab.gz", ".tsv", ".tsv.gz", ".csv", ".csv.gz"]:
		if lookup_file.endswith(i):
			return lookup_from_tab(lookup_file)
	raise IOError("\nERROR: do not recognize file type from suffix of \"{}\"\n".format(lookup_file))

def lookup_from_gb(lookup_file):
	sys.stderr.write("\n\t--> determined as \"genbank\" based on suffix\n")
	from Bio import SeqIO
	lookup_dict = {}
	infile = openfile(lookup_file)
	for record in SeqIO.parse(infile, "genbank"):
		contig = record.id
		for feature in record.features:
			if feature.type == "CDS":
				identifier = feature.qualifiers["locus_tag"][0]
				lookup_dict[identifier] = contig
	infile.close()
	return lookup_dict

def lookup_from_gff(lookup_file):
	sys.stderr.write("\n\t--> determined as \"gff\" based on suffix\n")
	lookup_dict = {}
	infile = openfile(lookup_file)
	for line in infile:
		if line.startswith("#"):
			if line.startswith("##FASTA"):
				break #stop reading file if sequence is reached
			continue #ignore comment and header lines
		linetokens = line.split()
		#import pdb; pdb.set_trace()
		ncbicompliant_contigname = linetokens[0]
		contigname = ncbicompliant_contigname.split("|")[-1] #prokka may create ncbi-compliant output that has this scheme "<seqmeth>|<center>|contigname". This should take ony the actual contigname
		featuretype = linetokens[2]
		comment = linetokens[8]
		if featuretype == "CDS": #ignore other feature types (such as "gene")
			#WARNING: better make this regex-based in the future! This is just a quick-and-dirty solution
			comment_tokens = comment.split(";")
			for ct in comment_tokens:
				if ct.startswith("ID="): #double checking, if ID is listed in nondefault position in comments
					locus_tag = ct[3:].strip() #remove prefixing "ID="-substring and remove line ending (in case id is at end of comment-line)
					break
				raise RuntimeError("Could not find ID for CDS in comment of this gff line:\n{}\n".format(line))
			lookup_dict[locus_tag] = contigname
	infile.close()
	return lookup_dict

def lookup_from_tab(lookup_file):
	sys.stderr.write("\n\t--> determined as tab-seperated table based on suffix\n")
	lookup_dict = {}
	infile = openfile(lookup_file)
	for line in infile:
		if line.startswith("#") or len(line.split("\t")) < 2:
			continue
		tokens = line.strip().split("\t")
		lookup_dict[tokens[0]] = tokens[1]
	infile.close()
	return lookup_dict

def add_to_dict_or_not(contig_dict, contigname, columns): #TODO: this can be sped up by setting a switch that automatically discards all further hits od the same protein, if the previous was already below score-cutoff
	#print columns
	bitscore = columns[11]
	protein = columns[0]
	#print "{}\n{}\n{}\n".format(bitscore, protein, contigname)
	addkept, addprot = 1, 1
	if contigname in contig_dict:
		if protein in contig_dict[contigname]:
			addprot = 0
			#print "{} >= {}??".format(bitscore, contig_dict[contigname][protein][0][10] * args.score_cutoff)
			if bitscore >= contig_dict[contigname][protein][0][10] * args.score_cutoff: # "contig_dict[contigname][protein][0][11]"corresponds to the bit-score of the first (best) hit of that protein
				#print "---> YAY!"
				contig_dict[contigname][protein].append(columns[1:])
			else:
				#print "---> NOOOO!"
				addkept = 0
		else:
			contig_dict[contigname][protein] = [columns[1:]]
	else:
		contig_dict[contigname] = {protein : [columns[1:]]}
	return contig_dict, addkept, addprot

def readfile(input_table, lookup_dict = None):
	sys.stderr.write("\nReading input data...\n")
	if not (os.path.exists(args.input_table) and os.path.isfile(input_table)): #todo: add support for compressed input-tables
		raise IOError("\n{} does not exist or is not a file\n".format(input_table))
	intable = []
	infile = openfile(input_table)
	pattern = re.compile("_\d+$")
	contig_dict = {}
	contig_count, protein_count, totalhitcount, kepthitcount = 0, 0, 0, 0
	for line in infile:
		if line.startswith("#"):#headers MAY be added as comments
			continue
		totalhitcount += 1
		columns = line.rstrip().split("\t")
		columns[11] = float(columns[11])# convert score value from string to float for correct sorting later on
		intable.append(columns[0:12]) #use only first 12 columns and discard any additional custom fields (Most Downstream tools gets confused if it finds custom fields)
		if not lookup_dict:
			assert re.search(pattern, columns[0]), "\nERROR: protein {} does not seem to be named according to prodigal standards. You should supply a lookup-table\n".format(columns[0])
			contigname = re.sub(pattern, "", columns[0])
		else:
			contigname = lookup_dict[columns[0]]
		contig_dict, addkept, addprot = add_to_dict_or_not(contig_dict, contigname, columns)
		kepthitcount += addkept
		protein_count += addprot
	sys.stderr.write("\nkept {} of {} hits to {} proteins on {} contigs\n".format(kepthitcount, totalhitcount, protein_count, len(contig_dict)))
	return contig_dict

def get_maxhits(contig, contig_counts):
	if args.maxhits > 1:
		return args.maxhits
	else:
		return contig_counts[contig] * args.maxhits

def process_contig_dict(contig_dict):
	sys.stderr.write("\nprocessing data\n")
	from operator import itemgetter
	final_list = []
	for contig in contig_dict:
		contig_list = []
		for prot in contig_dict[contig]:
			contig_list.extend([ [contig] + x for x in contig_dict[contig][prot] ])
		totalcontighits = len(contig_list)
		if args.maxhits > 0 and args.maxhits <= 1:
			hitcount_cutoff = int(round(totalcontighits * args.maxhits))
			if hitcount_cutoff < 3: #always use at least the three best hits per contig (if available)
				hitcount_cutoff = 3
		else:
			hitcount_cutoff = int(args.maxhits)
		#print contig_list
		final_list.extend(sorted(contig_list, key=itemgetter(11), reverse=True)[:hitcount_cutoff])
	return final_list

def writefile(final_list, outfile_name):
	sys.stderr.write("\nwriting to output_file \"{}\"\n".format(outfile_name))
	outfile=open(outfile_name, "w")
	for column_list in final_list:
		outfile.write("{}\n".format("\t".join([str(x) for x in column_list])))

def main():
	assert args.score_cutoff > 0 and args.score_cutoff <= 1, "\nERROR: score_cutoff must be set between 0 - 1\n"
	if args.lookup_file != None:
		lookup_dict = get_lookup_dict(args.lookup_file)
	else:
		lookup_dict = None
	contig_dict = readfile(args.input_table, lookup_dict)
	final_list = process_contig_dict(contig_dict)
	writefile(final_list, args.output_table)

main()
