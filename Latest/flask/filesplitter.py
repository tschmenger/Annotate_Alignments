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
donnus = []
dirs_done = []
# with gzip.open("Data_V0_20250104.txt.gz") as content:
#     for oldline in content:
#         line = oldline.decode('utf-8')
#         if line.strip() and "None" not in line:
#             primary_acc = line.split("\t")[0].replace(" ","")
#             subdirname  = primary_acc[0:4]
#             if subdirname not in dirs_done:
#                 os.system("mkdir static/uniprotdata/"+str(subdirname))
#                 filenamus = primary_acc+".txt"
#                 with open("static/uniprotdata/"+str(subdirname)+"/"+filenamus, "a") as outfl:
#                     outfl.write(line)
#                 outfl.close()
#                 dirs_done.append(subdirname)
#             else:
#                 filenamus = primary_acc+".txt"
#                 with open("static/uniprotdata/"+str(subdirname)+"/"+filenamus, "a") as outfl:
#                     outfl.write(line)
#                 outfl.close()

with open("interpro_data_V0_20250116.txt") as content:
    for line in content:
        if line.strip() and "None" not in line:
            primary_acc = line.split("\t")[0].replace(" ","")
            subdirname  = primary_acc[0:4]
            if subdirname not in dirs_done:
                os.system("mkdir interprodata/"+str(subdirname))
                filenamus = primary_acc+"_interpro.txt"
                with open("interprodata/"+str(subdirname)+"/"+filenamus, "a") as outfl:
                    outfl.write(line)
                outfl.close()
                dirs_done.append(subdirname)
            else:
                filenamus = primary_acc+"_interpro.txt"
                with open("interprodata/"+str(subdirname)+"/"+filenamus, "a") as outfl:
                    outfl.write(line)
                outfl.close()
