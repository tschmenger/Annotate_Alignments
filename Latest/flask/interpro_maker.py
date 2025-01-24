import sys
import urllib.request
import re
nums = re.compile(r"[+-]?\d+(?:\.\d+)?")
whitespace_killer=re.compile(r"\s+")
import time
import ast
import logging
import os
import gzip
###############################################################################
def interprodownloader(identif):
    BASE_URL 	= "https://www.ebi.ac.uk/interpro/api/protein/UniProt/"+identif+"/?residues&page_size=200"
    req 		= urllib.request.urlopen(BASE_URL)
    the_page 	= req.read().decode('utf-8')
    interpro 	= ast.literal_eval(the_page)
    interpro_processed = {}
    for k in interpro:
        for loca in interpro[k]["locations"]:
            descr = ""
            for b in loca["description"]:	###              'locations': [{'description': 'GEF (guanine nucleotide exchange factor) interaction site',
                descr = descr+b
            if "(" in descr:
                front_descr = descr.split("(")[0]
                back_descr = descr.split(")")[1]
                descr = front_descr+back_descr
            for categ in loca["fragments"]: ### 		  categ are dictionaries, again, because the nesting never ends here
                for element in categ:
                    starting = categ["start"]
                    ending = categ["end"]
                    if int(starting) == int(ending):
                        residue = int(starting)
                    if identif not in interpro_processed:
                        interpro_processed[identif]={}
                        interpro_processed[identif][descr]=[residue]
                    elif descr not in interpro_processed[identif]:
                        interpro_processed[identif][descr]=[residue]
                    else:
                        interpro_processed[identif][descr].append(residue)
    return	interpro_processed
###############################################################################
try:
    os.system("rm donnus.txt")
except:
    pass
os.system("cat interpro_data_V0_20250116.txt |cut -f1|uniq >> donnus.txt")
protaccessions = []
progress = []
with open("donnus.txt") as progressfile:
    for line in progressfile:
        unip = line.split("\t")[0].replace(" ","")
        if unip not in progress:
            progress.append(unip)
del progress[-1] #when I interrupt the previous run I cannot be sure that the last uniprot entry has completed so I need to start the new run with this last item
humancontainer = []
with gzip.open("Uniprot_Acc_20250104_V2.txt.gz") as content:
    for oldline in content:
        line = oldline.decode('utf-8')
        if line.strip() and "None" not in line:
            primary_acc = line.split("\t")[0]
            if primary_acc not in protaccessions and primary_acc not in progress:
                protaccessions.append(primary_acc.replace("\n","").replace(" ",""))
                if "HUMAN" in line:
                    humancontainer.append(primary_acc.replace("\n","").replace(" ",""))

try:
    for entry in humancontainer:
        feature_dict = {}
        try:
            feature_dict = interprodownloader(entry)
            for k in feature_dict:
                for v in feature_dict[k]:
                    print(k,"\t",v,"\t",feature_dict[k][v])
        except:
            #logging.exception("message")
            pass
except:
    #logging.exception("message")
    pass

for entry in protaccessions:
    feature_dict = {}
    try:
        feature_dict = interprodownloader(entry)
        for k in feature_dict:
            for v in feature_dict[k]:
                print(k,"\t",v,"\t",feature_dict[k][v])
    except:
        #logging.exception("message")
        pass
