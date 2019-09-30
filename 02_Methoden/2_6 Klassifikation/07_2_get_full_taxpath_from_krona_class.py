#!/usr/bin/env python
import argparse
import sys
import os
import re
from Bio import SeqIO
from operator import itemgetter
# TODO: adjust cutoffs for minimum hits to accept "most-abundant" annotation for each taxon level
# TODO: differentiate between "most likely taxon" (supported by at least 10% of the contigs/assembly) and "most closely related taxon" supported by less

version = "0.2b (20.04.19)"
name = os.path.basename(sys.argv[0]) #get scriptname from actual script filename that was called
parser=argparse.ArgumentParser(description="get full taxonomic paths from KRONA LCA classifications")
#input_group = parser.add_mutually_exclusive_group(required=True)
parser.add_argument('-i','--input_classif', action = "store", dest = "input_tables", nargs = "+", required = True, help = "Input Classification(s) (Output_tables from \"ktClassifyBLAST\" OR rnammer/SINA). If multiple files are given, they will be considered as hierarchically ranked from highest to lowest significance (e.g. 16S -> 23S -> marker prots -> total prots). Each contig wil be annoated according to the highes ranking classification (referring to the next level if classification = \"unknown\")")
#input_group.add_argument('-hi','--multiple_input', action = "store", nargs = "+", dest = "hierarchical_inputs", help = "multipe classification files from SINA AND ktClassifyBLAST. MUST BE ORDERED FROM HIGHEST TO LOWEST SIGNIFICANCE (e.g. 16SrRNA -> 23S rRNA -> marker-prots -> total prots)! Will classify each contig according to the highest ranking classification (or the next highest if the highest turns out to be \"unknown\")")
parser.add_argument('-f', '--fasta', action = "store", dest = "fasta", default = None, help= "fatafile (or at least file with fasta_headers) belonging to targen genome. Required for excluding annotations from other genomes (e.g. combined SINA) and/or for filtering out contaminations based on taxonomy)")
parser.add_argument('--taxfile', action = "store", dest = "taxfile", default = None, help = "Full path to the KRONA taxonomy-file (including the actual filename). Script will attempt to find it itself if not specified here)")
parser.add_argument('--filter', action = "store", dest = "taxfilter", choices = ["p", "phylum", "c", "class", "o", "order", "f", "family", "g", "genus", "s", "species", "None"], default = "none", help = "Filter out all contigs with differing annotations (but ignoring \"unknown\") at the specified level. Default : None (= do not filter)")
parser.add_argument("--filter_cutoff", action = "store", dest = "filter_cutoff", type = float, default = 40.0, help = "Minimum average BLAST identity for considering contig-classification for filtering. default = \"40.0\". classifications with confidence below this value will not be filtered. (Set to \"0.0\" if ALL classifications schould be considered)")
parser.add_argument("--unofficials", action = "store_true", dest = "include_unofficials", default = False, help = "Also include \"unofficial\" Candidate phyla in Taxon counts (default False)")
parser.add_argument('-o','--output_basename', action = "store", dest = "output_prefix", default = None, help = "Prefix for outputfiles (default : \"<First_Input_filename>\")")
parser.add_argument('-V', '--version', action = "version", version = name + " version " + version)
args=parser.parse_args()

#TODO: add option to use read-coverage as additional weight

#Notes:
#Tax_table-format:
# <Tax-id> <depth> <parent> <rank> <name>
#some things to know:
# candidate phyla also have a rank that is called "phylum"
#    but they have an indermediate rank (at depth 4) called "Bacteria candidate phyla" (--> can be used to automatically recognize them if ever needed)
official_ranks = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
krona_table = "pathtotable.tab"

def openfile(infile):
	assert os.path.exists(infile), "Could not find file %s" %infile
	try:
		if infile.endswith(".gz") or infile.endswith(".gzip"):
			import gzip
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

def read_infasta(infile):
	readfile = openfile(infile)
	return { record.id.split("|")[-1]: record for record in SeqIO.parse(readfile, "fasta")} #remove ncbi stype stuff from prokka contignames

def get_taxonomyfile_path():#TODO: currently this may only work for conda-installations of KRONA. Find out where other installation forms store their taxonomy files
	import time
	for path in os.environ["PATH"].split(os.pathsep):
		searchterm=os.path.join(path, "ktImportTaxonomy")
		if os.path.exists(searchterm):
			if os.path.islink(searchterm):
				scriptdir = os.path.dirname(os.path.realpath(os.path.join(path, os.readlink(searchterm))))
				taxonomydir = os.path.join(scriptdir, "../taxonomy")
				if os.path.exists(taxonomydir):
					sys.stderr.write("\n\t... located most probable taxonomy path at \"{}\"\n".format(taxonomydir))
					taxonomyfile = os.path.join(taxonomydir, "taxonomy.tab")
					if os.path.exists(taxonomyfile) and os.path.isfile(taxonomyfile):
						sys.stderr.write("\t... located taxonomy-file at \"{}\"\n".format(taxonomyfile))
						file_age_in_days = (time.time() - os.path.getmtime(taxonomyfile)) / (60*60*24)
						sys.stderr.write("\t  --> taxonomy file is {:.1f} days old".format(file_age_in_days))
						if file_age_in_days >= 180:
							sys.stderr.write("\n\n\t=== WARNING: taxonomy file is OLDER THAN 6 MONTHS! ===\n\t=== You should consider updating it! \n\n")
						else:
							sys.stderr.write("  --> OK!\n")
						return taxonomyfile
	raise IOException("\nERROR: could not locate the KRONA \"taxonomy.tab\". please update/create it or provide a precise path to it\n")

def open_taxtable(infilename):
	sys.stderr.write("\nLoading Taxonomy...\n")
	#official_ranks = ["species", "genus", "family", "order", "class", "phylum"] #wanted to record only official ranks at first, but ditch that. record ALL
	taxIDdict = {}
	taxNAMEdict = {}
	infile = openfile(infilename)
	for line in infile:
		tokens = line.strip().split("\t")
		taxid = int(tokens[0]) #integers use less memory than strings
		depth = int(tokens[1])
		parent = int(tokens[2])
		rank = tokens[3]
		name = tokens[4]
		taxIDdict[taxid] = {"depth" : depth, "parent" : parent, "rank" : rank, "name" : name}
		taxNAMEdict[name] = {"taxid" : taxid, "rank" : rank, "depth" : depth} #ignoring unfortunate double entry of Actinobacteria (class AND phylum) here
	infile.close()
	#print taxdict[1224]
	return taxIDdict, taxNAMEdict

def get_full_taxpath(taxid, taxdict):
	#print "testing"
	#print taxdict[1224]
	prefixdict = {"species": "s", "genus": "g", "family" : "f", "order" : "o", "class" : "c", "phylum" : "p", "superkingdom" : "d"}
	official_rank_dict = { rank : None for rank in official_ranks }
	taxpath = [ "" ] * (taxdict[taxid]["depth"])
	backup_phylum = None
	is_candidate_phylum = False
	itercount = 0
	currtax = taxid
	safelock_maxiter = 100 #safeguard to make sure there are no "dead-end" taxons that do not lead to root# spoiler: there seem to be none
	while taxdict[currtax]["depth"] != 0:
		itercount += 1
		if itercount > safelock_maxiter:
			raise RunTimeError("\nERROR: to many iterations for taxon \"{}\"\n".format(taxid))
		currdepth = taxdict[currtax]["depth"]
		currparent = taxdict[currtax]["parent"]
		currrank = taxdict[currtax]["rank"]
		if currrank in prefixdict:
			prefix = prefixdict[currrank]
			official_rank_dict[currrank] = taxdict[currtax]["name"]
		else:
			prefix = "nr"
			if taxdict[currtax]["name"].startswith("Candidatus"):
				backup_phylum = taxdict[currtax]["name"]
			elif taxdict[currtax]["name"] == "Bacteria candidate phyla":
				is_candidate_phylum = True
		currname = taxdict[currtax]["name"]
		taxpath[taxdict[currtax]["depth"]-1] = "{}_{}".format(prefix, currname.replace(" ", "_"))
		currtax = currparent
	if args.include_unofficials and official_rank_dict["phylum"] == None and is_candidate_phylum: #annotate also lesser known candidate phyla
		official_rank_dict["phylum"] = "Unofficial_{}".format(backup_phylum)
	return {"taxpath" : taxpath, "ranks" : official_rank_dict}

def parsetable(intablename, taxIDdict, taxNAMEdict, fastadict=None):
	sys.stderr.write("\ntrying to recognize input type (SINA or KronaClass) by headerline of input file (will fail if header was removed or modified)\n")
	infile = openfile(intablename)
	headerline = infile.readline()
	if headerline.startswith("#queryID\ttaxID\tAvg. % identity"):
		return parse_kronatable(intablename, taxIDdict, fastadict)
	elif headerline.startswith("\"job_id\";\"sequence_number\";\"sequence_identifier\""):
		return parse_sinatable_csv(intablename, taxIDdict, taxNAMEdict, fastadict) # TODO: add additional check for "tax_slv" to make sure file is actually SINA classified
	elif headerline.startswith("job_id\tsequence_number\tsequence_identifier\t"):
		return parse_sinatable_tsv(intablename, taxIDdict, taxNAMEdict, fastadict)
	else:
		raise IOError("\nERROR: cannot determine classification type of \"{}\"\n".format(intablename))

def parse_sinatable_csv(intablename, taxIDdict, taxNAMEdict, fastadict):
	sys.stderr.write("\n  parsing {} as SINA_csv\n".format(intablename))
	column_pattern = re.compile("\"[^\"]*\"")
	infile = openfile(intablename)
	outdict = {}
	for line in infile:
		tokens = [ x.strip("\"'") for x in re.findall(column_pattern, line) ]
		if tokens[0] == "job_id":
			continue #skip header
		rnammerID = tokens[2]
		originalID = get_original_seqID(rnammerID)
		if fastadict and (originalID not in fastadict):
			continue
		ident = float(tokens[8])
		sinatax = tokens[-1]
		taxid = sina2taxid(sinatax, taxIDdict, taxNAMEdict)
		if taxid == None:
			continue
		taxpath = get_full_taxpath(taxid, taxIDdict)
		outdict[originalID] = {"taxID":taxid, "taxpath": taxpath, "confidence": ident}
	#import pdb; pdb.set_trace()
	return outdict

def parse_sinatable_tsv(intablename, taxIDdict, taxNAMEdict, fastadict):
	sys.stderr.write("\n  parsing {} as SINA_tsv\n".format(intablename))
	column_pattern = re.compile("\"")
	infile = openfile(intablename)
	outdict = {}
	for line in infile:
		tokens = line.strip().split("\t")
		if tokens[0] == "job_id":
			continue #skip header
		rnammerID = tokens[2]
		originalID = get_original_seqID(rnammerID)
		if fastadict and (originalID not in fastadict):
			continue
		ident = float(tokens[8])
		sinatax = tokens[-1]
		taxid = sina2taxid(sinatax, taxIDdict, taxNAMEdict)
		if taxid == None:
			continue
		taxpath = get_full_taxpath(taxid, taxIDdict)
		outdict[originalID] = {"taxID":taxid, "taxpath": taxpath, "confidence": ident}
	#import pdb; pdb.set_trace()
	return outdict

def get_original_seqID(rnammerID): #CAUTION: original Contignames/SeqIDs must not contain "_rRNA_"! #TODO: maybe find a safer solution?
	prefix = re.compile("^.+_rRNA_")
	suffix = re.compile("_\d+-\d+_DIR[+-]$")
	return re.sub(suffix, "", re.sub(prefix, "", rnammerID))

def sina2taxid(sinatax, taxIDdict, taxNAMEdict):
	if sinatax.startswith("Unclassified") or sinatax.startswith("-"):
		return None
	taxonomy = sinatax.strip(";").split(";")
	for index in reversed(range(len(taxonomy))):
		name = taxonomy[index]
		#import pdb; pdb.set_trace()
		if name == "Actinobacteria": #resolve the unfortunate double-level "Actinobacteria" problem
			#import pdb; pdb.set_trace()
			if taxonomy[index -1 ] == "Actinobacteria":
				return 1760 #taxid of the CLASS "Actinobacteria"
			else:
				return 201174 #taxid of the PHYLUM "Actinobacteria"
		if name in taxNAMEdict:
			return taxNAMEdict[name]["taxid"]
		if "Candidatus {}".format(name) in taxNAMEdict:
			return taxNAMEdict["Candidatus {}".format(name)]["taxid"]
		if "{} group".format(name) in taxNAMEdict:
			return taxNAMEdict["{} group".format(name)]["taxid"]
	return None #Return None if classification not fount in KONA/NCBI-lookuptable on any level

def parse_kronatable(intablename, taxdict, fastadict):
	sys.stderr.write("\n  parsing {} as KRONA table\n".format(intablename))
	infile = openfile(intablename)
	#outfile = open(outtablename, "w")
	#header = "identifier\ttax_id\ttax_path\tconfidence\n"
	#outfile.write(header)
	outdict = {}
	for line in infile:
		if line.startswith("#"):
			continue
		tokens = line.strip().split("\t")
		contig = tokens[0]
		if fastadict and (contig not in fastadict):
			continue
		taxid = int(tokens[1])
		confidence = tokens[2]
		taxpath = get_full_taxpath(taxid, taxdict)
		outdict[contig] = {"taxID":taxid, "taxpath": taxpath, "confidence": confidence}
		#outline = "{c}\t{ti}\t{op}\t{conf}\n".format(c = contig, ti = taxid, op = ";".join(taxpath), conf = confidence)
		#outfile.write(outline)
	infile.close()
	#outfile.close()
	#sys.stderr.write("\n--> finished! Wrote results to {}\n".format(outtablename))
	return outdict

def get_taxon_count(taxon_counts, taxinfo, weight): #TODO add actual count for cases where contiglengths are used as wheigts
	if weight == None:
		weight = 1
	tsuperkingdom=taxinfo["taxpath"]["ranks"]["superkingdom"]
	tphylum = taxinfo["taxpath"]["ranks"]["phylum"]
	tclass = taxinfo["taxpath"]["ranks"]["class"]
	torder = taxinfo["taxpath"]["ranks"]["order"]
	tfamily = taxinfo["taxpath"]["ranks"]["family"]
	tgenus = taxinfo["taxpath"]["ranks"]["genus"]
	tspecies = taxinfo["taxpath"]["ranks"]["species"]
	if tsuperkingdom in taxon_counts:
		taxon_counts[tsuperkingdom]["totalweight"] += weight
		taxon_counts[tsuperkingdom]["totalcount"] += 1
	else:
		taxon_counts[tsuperkingdom] = {"totalweight" : weight, "totalcount" : 1}
	if tphylum in taxon_counts[tsuperkingdom]:
		taxon_counts[tsuperkingdom][tphylum]["totalweight"] += weight
		taxon_counts[tsuperkingdom][tphylum]["totalcount"] += 1
	else:
		taxon_counts[tsuperkingdom][tphylum] = {"totalweight" : weight, "totalcount" : 1}
	if tclass == None:
		return taxon_counts
	if tclass in taxon_counts[tsuperkingdom][tphylum]:
		taxon_counts[tsuperkingdom][tphylum][tclass]["totalweight"] += weight
		taxon_counts[tsuperkingdom][tphylum][tclass]["totalcount"] += 1
	else:
		taxon_counts[tsuperkingdom][tphylum][tclass] = {"totalweight" : weight, "totalcount" : 1}
	if torder == None:
		return taxon_counts
	if torder in taxon_counts[tsuperkingdom][tphylum][tclass]:
		taxon_counts[tsuperkingdom][tphylum][tclass][torder]["totalweight"] += weight
		taxon_counts[tsuperkingdom][tphylum][tclass][torder]["totalcount"] += 1
	else:
		taxon_counts[tsuperkingdom][tphylum][tclass][torder] = {"totalweight" : weight, "totalcount" : 1}
	if tfamily == None:
		return taxon_counts
	if tfamily in taxon_counts[tsuperkingdom][tphylum][tclass][torder]:
		taxon_counts[tsuperkingdom][tphylum][tclass][torder][tfamily]["totalweight"] += weight
		taxon_counts[tsuperkingdom][tphylum][tclass][torder][tfamily]["totalcount"] += 1
	else:
		taxon_counts[tsuperkingdom][tphylum][tclass][torder][tfamily] = {"totalweight" : weight, "totalcount" : 1}
	if tgenus in taxon_counts[tsuperkingdom][tphylum][tclass][torder][tfamily]:
		taxon_counts[tsuperkingdom][tphylum][tclass][torder][tfamily][tgenus]["totalweight"] += weight
		taxon_counts[tsuperkingdom][tphylum][tclass][torder][tfamily][tgenus]["totalcount"] += 1
	else:
		taxon_counts[tsuperkingdom][tphylum][tclass][torder][tfamily][tgenus] = {"totalweight" : weight, "totalcount" : 1}
	if tspecies == None:
		return taxon_counts
	if tspecies in taxon_counts[tsuperkingdom][tphylum][tclass][torder][tfamily][tgenus]:
		taxon_counts[tsuperkingdom][tphylum][tclass][torder][tfamily][tgenus][tspecies]["totalweight"] += weight
		taxon_counts[tsuperkingdom][tphylum][tclass][torder][tfamily][tgenus][tspecies]["totalcount"] += 1
	else:
		taxon_counts[tsuperkingdom][tphylum][tclass][torder][tfamily][tgenus][tspecies] = {"totalweight" : weight, "totalcount" : 1}
	return taxon_counts

def combine_taxclass(outdicts, fastadict):#TODO: add option to use read-coverage as additional weight
	done_ids = set()
	taxon_counts = {} #should be nested {domain : {"totalweight"=nd, "phyla" = {phylum1: {"totalweight":np}, "classes": {class1: etc...}}} for domain in anotation }
	combined_dict = {}
	if fastadict:
		contigs = fastadict.keys()
	else:
		contigs = set()
		for od in outdicts:
			for contig in od.keys():
				contigs.add(contig)
		contigs = list(contigs)
	#isroot = True
	#import pdb; pdb.set_trace()
	for contig in contigs:
		#import pdb; pdb.set_trace()
		#print "\n--------\n{}".format(contig)
		taxinfo = None
		if fastadict:
			length = len(fastadict[contig])
		else:
			length = None
		for od in range(len(outdicts)):
			#print "-->{}".format(od)
			if contig in outdicts[od]:
				#print "IT EXISTS"
				taxinfo = outdicts[od][contig]
				#print taxinfo
				#print contig
				if taxinfo["taxpath"]["ranks"]["phylum"] != None: #Currently: accepting highest level annotation IF at least annotated to phylum. otherwise looks in next lower level
					taxon_counts = get_taxon_count(taxon_counts, taxinfo, length)
					#print "annotated at least to phylum"
					break
		if taxinfo:
			outinfo = {"taxID" : taxinfo["taxID"], "ranks" : taxinfo["taxpath"]["ranks"], "taxpath" : taxinfo["taxpath"]["taxpath"], "confidence" : taxinfo["confidence"], "inference" : "{}/{}".format(od + 1, len(outdicts))}
		else:
			#print "Never heard of it"
			outinfo = {"taxID" : -1, "ranks" : { rank : None for rank in official_ranks}, "taxpath" : ["unclassified"], "confidence" : 0, "inference" : "-"}				
		if fastadict:
			outinfo["length"] = length
		combined_dict[contig] = outinfo
	return combined_dict, taxon_counts

def write_combined_taxclass(combined_dict, outfilename):
	outfile = open(outfilename, "w")
	header = "contig\ttaxID\ttaxpath\tav._ident\tinference\n"
	if "length" in combined_dict[combined_dict.keys()[0]]:
		header = "contig\ttaxID\ttaxpath\tav._ident\tinference\tlength\n"
	outfile.write(header)
	for contig in combined_dict:
		taxinfo = combined_dict[contig]
		if "length" in taxinfo:
			length = "\t{}".format(taxinfo["length"])
		else:
			length = ""
		#print "--------------- {} -----------".format(contig)
		#print taxinfo
		outline = "{c}\t{ti}\t{op}\t{conf}\t{inf}{l}\n".format(c = contig, ti = taxinfo["taxID"], op = ";".join(taxinfo["taxpath"]), conf = taxinfo["confidence"], inf = taxinfo["inference"], l = length)
		#print outline
		outfile.write(outline)
	outfile.close()

def get_predominant_taxa(taxon_counts, total_weights, total_count):
	testdict = taxon_counts #todo: add checks at every level that there IS a next level
	
	#max_sk, max_p, max_c, max_o, max_f, max_g, max_s = None, None, None, None, None, None, None
	predom_dict = {taxon : None for taxon in official_ranks}
	predomabund_dict = {taxon : 0 for taxon in official_ranks}
	predomcount_dict = {taxon : 0 for taxon in official_ranks}
	if args.fasta:
		countcutoff = 2000 #combined length must be more than 2 kbp (in order to exclude filtering based on spurious annotation on just one or to 500bp contigs)
		#maybe change this to make up at least 10% of total size?
	else:
		countcutoff = 2 # if no length-weights: at least two contigs must be annotated as this soecific taxon in order to filter for it
		#maybe change this to require at least 10% of total contigs?
	for r in official_ranks:
		#print "-----------------------------"
		#print r
		#print testdict
		sortedbydom = sorted([(x, testdict[x].get("totalweight"), testdict[x].get("totalcount")) for x in testdict if (x != "totalweight" and x != "totalcount") ], key=itemgetter(1), reverse = True)
		#sortedbydom = sorted(testdict.iteritems(), key=lambda x:x[1].get("totalweight") if type(x[1]) == dict else None, reverse = True) # <-- this one was troublesome (I'm no good with lambda functions)
		#print sortedbydom
		#~ print "********"
		#~ print r
		#~ print sortedbydom
		#~ print "********"
		#~ print "-----"
		#~ print sortedbydom[0]
		#~ print "------"
		if len(sortedbydom) > 0: #apparently there can be empty lists in some cases --> toDO: FIND OUT WHY!
			mostdom = sortedbydom[0][0]
			mostdom_abund = sortedbydom[0][1]
			mostdom_count = sortedbydom[0][2]
			if len(sortedbydom) > 1:
				secondmostdom = sortedbydom[1][0]
				secondmostdom_abund = sortedbydom[1][1]
				secondmostdom_count = sortedbydom[1][2]
				if mostdom != None and secondmostdom != None:#there is only really a conflict if none of the two equally high ranking taxons are "None"
					if mostdom_abund == secondmostdom_abund:
						sys.stderr.write("\nWARNING: two taxa with equal abundance at taxlevel \"{}\":\n  \"{}\"  and  \"{}\"\n  --> not able to determine most abundand taxon for this level and above\n!".format(r, mostdom, secondmostdom))
						break
				elif mostdom == None and secondmostdom != None: #prefer non-None annotations IF they fulfill minimum criteria
					if secondmostdom_abund >= countcutoff:# TODO: Add third level --> what is if most abundant is "unclassified" and second&third are equally abundant?
						mostdom = secondmostdom
						mostdom_abund = secondmostdom_abund
						mostdom_count = secondmostdom_count
					else:
						break
			elif mostdom == None or mostdom_abund < countcutoff: #also do not accept mostdom if it does not fulfil minimum criteria
				break
		else:
			mostdom = None
			mostdom_abund = 0
			mostdom_count = 0
		predom_dict[r] = mostdom
		predomabund_dict[r] = mostdom_abund
		predomcount_dict[r] = mostdom_count
		if mostdom not in testdict or not testdict[mostdom] or len(testdict[mostdom])<=2: #must have at least "totalcounts" key AND a "totalweights" key AND a taxon key
			break
		debugdict = testdict.copy()
		testdict = testdict[mostdom]
	#~ max_sk = max(taxoncounts.iterkeys(), key=(lamdba key: taxoncounts[key]["totalweight"]))
	#~ max_p = max(taxoncounts[max_sk].iterkeys(), key=(lambda key: taxoncounts[max_sk][key]["totalweight"]))
	#~ max_c = max(taxoncounts[max_sk][max_p].iterkeys(), key=(lambda key: taxoncounts[max_sk][max_p][key]["totalweight"]))
	#~ max_o = max(taxoncounts[max_sk][max_p][max_c].iterkeys(), key=(lambda key: taxoncounts[max_sk][max_p][max_c][key]["totalweight"]))
	#~ max_f = max(taxoncounts[max_sk][max_p][max_c][max_o].iterkeys(), key=(lambda key: taxoncounts[max_sk][max_p][max_c][max_o][key]["totalweight"]))
	#~ max_g = max(taxoncounts[max_sk][max_p][max_c][max_o][max_f].iterkeys(), key=(lambda key: taxoncounts[max_sk][max_p][max_c][max_o][max_f][key]["totalweight"]))
	#~ max_s = max(taxoncounts[max_sk][max_p][max_c][max_o][max_f].iterkeys(), key=(lambda key: taxoncounts[max_sk][max_p][max_c][max_o][max_f][key]["totalweight"]))
	currentorg = args.output_prefix
	#~ print predom_dict
	#~ print "------------------"
	#~ print predomabund_dict
	#~ print "------------------"
	#~ print predomcount_dict
	#~ print "-------------------"
	sys.stdout.write("\npredominant taxa\n{}\n{}\n".format("\t".join(["bin/genome"] + official_ranks), "\t".join([currentorg] + ["{}[:w{:.2f}%;c{:.2f}%]".format(str(predom_dict.get(r)), (float(predomabund_dict.get(r))/total_weights) * 100,(float(predomcount_dict.get(r))/total_count) * 100)  for r in official_ranks] )))
	sys.stdout.flush()
	#print predom_dict
	return predom_dict

def filter_fasta(combined_dict, predomdict, fastadict, outfastaname, outtabname, max_rank):
	sys.stderr.write("\nfiltering fasta\n\n")
	outfile = open(outfastaname, "w")
	poscount, negcount = 0, 0
	for contig in fastadict:
		will_write = True
		for r in official_ranks:
			if predomdict[r] == None:
				break
			if combined_dict[contig]["ranks"][r] != None and combined_dict[contig]["ranks"][r] != predomdict[r]:
				#print "{} : {} != {}".format(contig, combined_dict[contig]["ranks"][r], predomdict[r])
				#print "\t-->confidence {} >= {}?".format(combined_dict[contig]["confidence"], args.filter_cutoff)
				if float(combined_dict[contig]["confidence"]) >= args.filter_cutoff:
					#print "\t\tYES --> DELETING!"
					#print type(combined_dict[contig]["confidence"])
					#print type(args.filter_cutoff)
					negcount += 1
					will_write = False
				#else:
				#	print "--> WOULD delete {} \n     because \"{}\" != \"{}\". BUT confidence is too low ({})".format(contig, combined_dict[contig]["ranks"][r], predomdict[r], combined_dict[contig]["confidence"])
				break
			if r == max_rank:
				break
		if will_write:
			poscount += 1
			SeqIO.write([fastadict[contig]], outfile, "fasta")
			if sum([poscount, negcount]) % 100 == 0:
				sys.stderr.write("\rread {} contigs, kept {}, removed {}".format(sum([poscount, negcount]), poscount, negcount))
		else:
			combined_dict.pop(contig)
	sys.stderr.write("\rread {} contigs, kept {}, removed {}\n".format(sum([poscount, negcount]), poscount, negcount))
	sys.stderr.write("now writing filtered table\n")
	write_combined_taxclass(combined_dict, outtabname)
	outfile.close()



def main():
	for intable in args.input_tables:
		assert os.path.exists(intable), "\nERROR: could not find \"{}\"".format(intable)
	if args.fasta:
		fastadict = read_infasta(args.fasta)
	else:
		fastadict = None
	if args.output_prefix == None:
		args.output_prefix = args.input_tables[0]
	if args.taxfile:
		if not os.path.exists(args.taxfile):
			raise IOError("\nERROR: could not find KRONA taxonomy file \"{}\"\n".format(args.taxfile))
		else:
			taxfile = args.taxfile
	else:
		taxfile = get_taxonomyfile_path()
	taxIDdict, taxNAMEdict = open_taxtable(taxfile)
	outdicts = [ parsetable(it, taxIDdict, taxNAMEdict, fastadict) for it in args.input_tables ]
	#import pdb; pdb.set_trace()
	combined_dict, taxon_counts = combine_taxclass(outdicts, fastadict)
	write_combined_taxclass(combined_dict, "{}_FULLTAXPATH.tab".format(args.output_prefix))
	if args.fasta:
		total_weights = sum([ len(x) for x in fastadict.values()])
	else:
		total_weights = len(combined_dict)
	predom_dict = get_predominant_taxa(taxon_counts, total_weights, len(combined_dict))
	if args.taxfilter != "none":
		filter_fasta(combined_dict, predom_dict, fastadict, "{}_FILTERED_CONTIGS.fasta".format(args.output_prefix), "{}_FULLTAXPATH_FILTERED.tab".format(args.output_prefix),args.taxfilter)

main()


