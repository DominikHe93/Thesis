#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""
Spyder Editor

"""
import json
import urllib2
import pickle
import traceback
import sys

def load_pickle(pickelfilename):
	print "loading picklefiles"
	prokaryotes_dict = {}
	picklefile = open(pickelfilename, "rb")
	load_dict = pickle.load(picklefile)
	picklefile.close()
	for i in load_dict:
		prokaryotes_dict[load_dict[i]["name"]] = None
	print "done loading picke files"
	return prokaryotes_dict
##############################################
	

def taxoncounts(taxonlist, prokaryotes_dict, target, strangetaxadict, rowdata):  #strangetaxadict listet Eukaryoten und nicht in der NCBI Liste vorkommenden Prokaryoten auf
	rowdata["total_count"] = "0" 
	rowdata["target_count"] = "0"
	rowdata["prokaryote_count"] = "0"
	debuglistofnonprokaryotes = []
	
	for te in taxonlist:
		rowdata["total_count"] = str(int(rowdata["total_count"]) + te[1])
		
		if te[0].startswith("unclassified"):
			if te[0] == "unclassified (derived from Bacteria)":
				te[0] = "unclassified Bacteria"
				rowdata["prokaryote_count"] =str(int(rowdata["prokaryote_count"]) + te[1])
			else:
				continue
		elif te[0] in prokaryotes_dict:
			rowdata["prokaryote_count"] = str(int(rowdata["prokaryote_count"]) + te[1])		
			if te[0] == target:
				rowdata["target_count"] = str(int(rowdata["target_count"]) + te[1])
		elif te[0].startswith("Candidatus "):
			if te[0][11:] in prokaryotes_dict:
				rowdata["prokaryote_count"] = str(int(rowdata["prokaryote_count"]) + te[1])		
			if te[0][11:] == target:
				rowdata["target_count"] = str(int(rowdata["target_count"]) + te[1])		
		else:
			strangetaxadict[te[0]] = None
			
		if int(rowdata["total_count"]) == 0:
			rowdata["perc_target_total"] = None
		else: 
			rowdata["perc_target_total"]= str((float(rowdata["target_count"]))/int(rowdata["total_count"]) * 100)		
		
		if int(rowdata["prokaryote_count"]) == 0:
			rowdata["perc_target_prokaryote"] = None
		else: 
			rowdata ["perc_target_prokaryote"] = str((float(rowdata["target_count"]))/int(rowdata["prokaryote_count"]) * 100)
	
	return rowdata, strangetaxadict
	
##############################################	

#Zugriff auf MG-Rast Daten

#Pfad: data, instance, mixs
def get_rowdata(i, tempid):
	
	#Reihenfolg der Daten
	
	#for i in tdatadict["data"]:
		
	#Eintraege auf "-" setzen
	rowdata={}
	for k in order:
		rowdata[k]= "-"
	rowdata["id"] = tempid	
	
	if "name" in i:
		rowdata["name"] = i["name"]
		
	if "project_metagenomes" in i:
		rowdata["project_metagenomes"] = str(len(i["project_metagenomes"]))
		
		
	#Pfad data, instance, mixs
	if "mixs" in i:
		if "country" in i["mixs"]:
			rowdata["country"] = i["mixs"]["country"]
		if "feature" in i["mixs"]:
			rowdata["feature"] = i["mixs"]["feature"]
		if "project_name" in i["mixs"]:
			rowdata["project_name"] = i["mixs"]["project_name"]
		if "seq_method" in i["mixs"]:
			rowdata["seq_method"] = i["mixs"]["seq_method"]
		if "latitude" in i["mixs"]:
			rowdata["latitude"] = str(i["mixs"]["latitude"])
		if "longitude" in i["mixs"]:
			rowdata["longitude"] = str(i["mixs"]["longitude"])
		if "collection_date" in i["mixs"]:
			rowdata["collection_date"] = i["mixs"]["collection_date"]
		if "material" in i["mixs"]:
			rowdata["material"] = i["mixs"]["material"] 
		if "sequence_type" in i["mixs"]:
			rowdata["sequence_type"] = i["mixs"]["sequence_type"]
		if "env_package_type" in i["mixs"]:
			rowdata["env_package_type_mixs"] = i["mixs"]["env_package_type"]
		if "biome" in i["mixs"]:
			rowdata["biome"] = i["mixs"]["biome"]
		if "project_id" in i["mixs"]:
			rowdata["project_id"] = i["mixs"]["project_id"]
		if "location" in i["mixs"]:
			rowdata["location"] = i["mixs"]["location"]
		
	#Pfad data, instance, metadata, library
	if "metadata" in i:
		if "library" in i["metadata"] and "data" in i["metadata"]["library"] and "seq_meth" in i["metadata"]["library"]["data"]:
			rowdata["seq_meth"] = i["metadata"]["library"]["data"]["seq_meth"]
		
	
	#Pfad data, instance, metadata, project
		if "project" in i["metadata"]:
			if "name" in i["metadata"]["project"]:
				rowdata["name"] = i["metadata"]["project"]["name"].replace("\n" , ";")
			if "id" in i["metadata"]["project"]:
				rowdata["project_id"] = i["metadata"]["project"]["id"].replace("\n" , ";")
			if "public" in i["metadata"]["project"]:
				rowdata["public"] = i["metadata"]["project"]["public"].replace("\n" , ";")
			if "project_description" in i["metadata"]["project"]["data"]:
				rowdata["project_description"] = i["metadata"]["project"]["data"]["project_description"].replace("\n" , ";").replace("\r" , "")
			if  "study_abstract" in i["metadata"]["project"]["data"]:
				rowdata["study_abstract"] = i["metadata"]["project"]["data"]["study_abstract"].replace("\n" , ";").replace("\r" , "")
			if  "study_description" in i["metadata"]["project"]["data"]:
				rowdata["study_description"] = i["metadata"]["project"]["data"]["study_description"].replace("\n" , ";").replace("\r" , "")
			if  "study_title" in i["metadata"]["project"]["data"]:
				rowdata["study_title"] = i["metadata"]["project"]["data"]["study_title"].replace("\n" , ";").replace("\r" , "")
		if "env_package" in i["metadata"] and "type" in i["metadata"]["env_package"]:
			rowdata["type"] = i["metadata"]["env_package"]["type"].replace("\n" , ";").replace("\r" , "")
		
	#Pfad data, instance, statistics, sequence_stats
	if "statistics" in i:
		if "sequence_count_raw" in  i["statistics"]["sequence_stats"]:
			rowdata["sequence_count_raw"] = str(i["statistics"]["sequence_stats"]["sequence_count_raw"])
		if "bp_count_raw" in  i["statistics"]["sequence_stats"]:
			rowdata["bp_count_raw"] = str(i["statistics"]["sequence_stats"]["bp_count_raw"])
	return rowdata
		
##############################################
	
#Veraenderbare Parameter
taxlevel = "phylum"
picklefile = "/data/ibg5_dominik/db/ncbi_taxonomy/{taxlevel}.pickle".format(taxlevel=taxlevel)

order=["id", "name", "biome","collection_date","env_package_type_mixs","feature","material",\
	"country", "location", "latitude","longitude","prokaryote_count", "perc_target_prokaryote","perc_target_total", "bp_count_raw", \
	"project_metagenomes", "sequence_type","seq_meth","seq_method","sequence_count_raw",\
	"target_count","total_count","type","project_name","project_id","public","project_description",\
	"study_abstract", "study_description", "study_title"]	

mainoutfilename = "mg_rast_data.tab"
strangefilename = "strangetaxa.list"
##############################################

def main():
	print "now in main()"
	strangetaxadict = {}
	offset = 0
	endpoint = 28000
	wantedtypes = ["wgs", "WGS", "shotgun metagenome"]
	portionsize = 1000
	strangedict = {}
	strangeseqtypedict = {}
	prokaryotes_dict = load_pickle(picklefile)
	#print prokaryotes_dict
	#Headerline der Tabelle
	print "will now try to get data..."
	headerline="\t".join(order)
	mainoutfile = open(mainoutfilename, "w")
	mainoutfile.write(headerline + "\n")
	while offset <= endpoint: 
		
		tdatap=urllib2.urlopen("http://api.mg-rast.org/metagenome?limit={ps}&order=name&offset={offs}&verbosity=full".format(ps=portionsize, offs=offset)).read()
	
		offset=offset + portionsize
		tdatadict=json.loads(tdatap)
		#print "mgid\ttotal_count\tprokaryote_count\tperc_target_total"
		for i in tdatadict["data"]:
			tempid = i["id"]
			print tempid
			try:
				tempseqtype = i["sequence_type"]
				if not tempseqtype in wantedtypes:
					strangeseqtypedict[tempseqtype]=None
					#print "{} is of unwanted seqtype '{}'".format(tempid, tempseqtype)
					continue
				if taxlevel in i["statistics"]["taxonomy"]:
					if len(i["statistics"]["taxonomy"][taxlevel]) != 0:
						taxonlist = i["statistics"]["taxonomy"][taxlevel]
				
				rowdata = get_rowdata(i, tempid)
				rowdata, strangetaxtadict = taxoncounts(taxonlist, prokaryotes_dict, "Chloroflexi", strangetaxadict, rowdata)
				outline_list=[]
				for k in order:
					outline_list.append(rowdata[k])
					outline="\t".join(outline_list) + "\n"
				#print outline
				mainoutfile.write(outline.encode("utf8"))
			except Exception as e:
				print "oops there was an error with \"{}\"".format(tempid)
				print e
				print traceback.format_exc()
				sys.stdout.flush()
				continue
			#print strangedict
			#outfile=open("mg_id.tab" , "a")
			#outfile.write(i["id"] + "\n")
	mainoutfile.close()
	strangefile = open(strangefilename, "w")
	strangefile.write("\n".join(strangedict.keys()))
	strangefile.close()
	strangeseqtypefile = open("strangeseqtypes.list", "w")
	strangeseqtypefile.write("\n".join(strangeseqtypedict.keys()))
	strangeseqtypefile.close()

print "hello"
main()
