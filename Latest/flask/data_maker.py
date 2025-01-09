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

inputdictionary_two = {}
translator = {}
with gzip.open("uniprot_sprot.dat.gz","rb") as content:
#with open("RHOA_test.txt","r") as content:
    for oldline in content:
        #print(oldline)
        line = oldline.decode('utf-8')
        newline = whitespace_killer.sub(" ",line).replace("\n","")
        try:
            if "ID" in newline and "Reviewed;" in newline and "AA." in newline:	### ID   RHOA_HUMAN              Reviewed;         193 AA.
                checktype = ""
                checker = "FALSE"
                position = ""
                comment = ""
                idcollector = []
                langname = newline.split(" ")[1].replace(" ","").replace("\n","")
            if "AC" in newline.split(" ")[0]:
                done = "false"
                for k in newline.split(" ")[1:]:
                    a = k.replace(";","").replace(" ","")
                    if done == "false":
                        if "-" in a and len(a)>1:
                            idcollector.append(a.split("-")[0])
                            if a.split("-")[0] not in translator: ### k is uniprot ID
                                translator[a.split("-")[0]]=[langname,[]]
                                done = str(a)
                        else:
                            if len(a)>1:
                                idcollector.append(a)
                                if a not in translator: ### k is uniprot ID
                                    translator[a]=[langname,[]]
                                    done = str(a)
                    else:
                        if "-" in a and len(a)>1:
                            translator[done][1].append(a.split("-")[0])
                        else:
                            if len(a)>1:
                                translator[done][1].append(a)
            #print(translator)
            #print(idcollector)
            if "FT" in newline.split(" ")[0] and "VARIANT" in newline:
                if checker == "FALSE":	### this will happen when I get to such a line:			FT   VARIANT         47
                    comment = ""
                    position = newline.split(" ")[2]
                    checker = "TRUE"
                    checktype = "VARIANT"
            elif "FT" in newline.split(" ")[0] and "MUTAGEN" in newline:
                if checker == "FALSE":	### this will happen when I get to such a line:			FT   MUTAGEN         226
                    comment = ""
                    position = newline.split(" ")[2]
                    checker = "TRUE"
                    checktype = "MUTAGEN"
            elif "BINDING" in newline:	### this will happen when I get to such a line:			FT   BINDING         10..17
                if "FT" in newline.split(" ")[0]:
                    if checker == "FALSE":	### this will happen when I get to such a line:			FT   MUTAGEN         226
                        comment = ""
                        position = newline.split(" ")[2]
                        checker = "TRUE"
                        checktype = "BINDING"
            elif "FT" in newline.split(" ")[0] and "MOD_RES" in newline:	### this will happen when I get to such a line:			FT   BINDING         38
                if checker == "FALSE":	### this will happen when I get to such a line:			FT   MUTAGEN         226
                    comment = ""
                    position = newline.split(" ")[2]
                    checker = "TRUE"
                    checktype = "MOD_RES"
            elif "FT" in newline.split(" ")[0] and "ACT_SITE" in newline:
                if checker == "FALSE":	### this will happen when I get to such a line:			FT   ACT_SITE        108
                    comment = ""
                    position = newline.split(" ")[2]
                    checker = "TRUE"
                    checktype = "ACT_SITE"
            else:
                if checker == "TRUE":	### this should happen right on the next line after setting checker to TRUE
                    moddedline = whitespace_killer.sub(" ",newline).replace("\n","").replace(" ","_")
                #	FT_/note="G->V:_Increased_Rho_protein_signal_transduction.
                #	FT_Constitutively_ active."
                #	FT_/evidence="ECO:0000269|PubMed:19948726,

                    if "FT_/note=" in moddedline:		### additional failsafe
                        if checktype == "MUTAGEN":
                            #FT   MUTAGEN         34
                            #FT                   /note="Y->A: Abolishes interaction with DGKQ."
                            #FT                   /evidence="ECO:0000269|PubMed:10066731,
                            if "Missing" in moddedline:
                                pass
                            else:
                                res_one = moddedline.split(":")[0].split("=\"")[1].split("->")[0]
                                res_two = moddedline.split(":")[0].split("=\"")[1].split("->")[1]
                                comment = moddedline.split(":")[1]
                        elif 	checktype == "VARIANT":
                            #FT   VARIANT         212
                            #FT                   /note="D -> V (in dbSNP:rs56143363)"
                            #FT                   /evidence="ECO:0000269|PubMed:17344846"
                            #FT                   /id="VAR_040389"
                            if "Missing" in newline:
                                pass
                            else:
                                res_one = moddedline.split("_")[1].split("=\"")[1]
                                res_two = moddedline.split("_")[3]
                                comment_raw = moddedline.split("_")[4:]
                                for item in comment_raw:
                                    comment = comment + item + "_"
                                comment = comment + "/"
                        elif checktype == "MOD_RES":
                            #FT   MOD_RES         689
                            #FT                   /note="Phosphoserine"
                            #FT                   /evidence="ECO:0007744|PubMed:23186163"
                            if "Missing" in moddedline:
                                pass
                            else:
                                binder = moddedline.split("=\"")[1]
                                comment = binder[:-1].replace("\"","")

                                checker = "FALSE"
                                homol_mutat = position
                                kommentar = checktype+"/"+comment
                                for k in idcollector:
                                    if k not in inputdictionary_two:
                                        inputdictionary_two[k]={}
                                        inputdictionary_two[k]["MOD_RES"]=[str(homol_mutat)]
                                    elif "MOD_RES" not in inputdictionary_two[k]:
                                        inputdictionary_two[k]["MOD_RES"]=[str(homol_mutat)]
                                    else:
                                        inputdictionary_two[k]["MOD_RES"].append(str(homol_mutat))
                                checker = "FALSE"
                                checktype = ""
                        elif checktype == "ACT_SITE":
                            #FT   ACT_SITE        108
                            #FT                   /note="Proton acceptor"
                            #FT                   /evidence="ECO:0000250"
                            if "Missing" in moddedline:
                                pass
                            else:
                                binder = moddedline.split("=\"")[1]
                                comment = binder[:-1].replace("\"","")
                                checker = "FALSE"
                                homol_mutat = position
                                if ".." in homol_mutat:	#FT   BINDING         10..17
                                    bindingposone = homol_mutat.split("..")[0]
                                    bindingpostwo = homol_mutat.split("..")[1]
                                    for bindpos in range(int(bindingposone),int(bindingpostwo)+1):
                                        for k in idcollector:
                                            if k not in inputdictionary_two:
                                                inputdictionary_two[k]={}
                                                inputdictionary_two[k]["ACT_SITE"]=[str(bindpos)]
                                            elif "ACT_SITE" not in inputdictionary_two[k]:
                                                inputdictionary_two[k]["ACT_SITE"]=[str(bindpos)]
                                            else:
                                                inputdictionary_two[k]["ACT_SITE"].append(str(bindpos))
                                else:
                                    for k in idcollector:
                                        if k not in inputdictionary_two:
                                            inputdictionary_two[k]={}
                                            inputdictionary_two[k]["ACT_SITE"]=[str(homol_mutat)]
                                        elif "ACT_SITE" not in inputdictionary_two[k]:
                                            inputdictionary_two[k]["ACT_SITE"]=[str(homol_mutat)]
                                        else:
                                            inputdictionary_two[k]["ACT_SITE"].append(str(homol_mutat))
                                checker = "FALSE"
                                checktype = ""
                        else:
                            pass
                    elif "/ligand=" in newline:
                        if checktype == "BINDING":
                        #FT   BINDING         79
                        #FT   BINDING         10..17
                        #FT                   /note="ATP"
                        #or	FT                   /note="Fatty acid"
                        #FT                   /evidence="ECO:0000255|PROSITE-ProRule:PRU00159"
                            if "Missing" in newline:
                                pass
                            else:
                                binder = newline.split("=\"")[1]
                                comment = binder[:-1].replace("\"","")
                                checker = "FALSE"
                                kommentar = checktype+"/"+comment
                                if ".." in position:	#FT   BINDING         10..17
                                    bindingposone = position.split("..")[0]
                                    bindingpostwo = position.split("..")[1]
                                    for bindpos in range(int(bindingposone),int(bindingpostwo)+1):
                                        for k in idcollector:
                                            if k not in inputdictionary_two:
                                                inputdictionary_two[k]={}
                                                inputdictionary_two[k]["BINDING"]=[str(bindpos)]
                                            elif "BINDING" not in inputdictionary_two[k]:
                                                inputdictionary_two[k]["BINDING"]=[str(bindpos)]
                                            else:
                                                inputdictionary_two[k]["BINDING"].append(str(bindpos))
                                else:
                                    homol_mutat = position
                                    for k in idcollector:
                                        if k not in inputdictionary_two:
                                            inputdictionary_two[k]={}
                                            inputdictionary_two[k]["BINDING"]=[str(homol_mutat)]
                                        elif "BINDING" not in inputdictionary_two[k]:
                                            inputdictionary_two[k]["BINDING"]=[str(homol_mutat)]
                                        else:
                                            inputdictionary_two[k]["BINDING"].append(str(homol_mutat))
                                checker = "FALSE"
                                checktype = ""
                    elif "/evidence" in newline:
                        # this typically marks the end of one entry
                        ### now I need to build the info I wanna actually store
                        #FT   MOD_RES         180
                        #FT                   /note="Cysteine methyl ester"
                        #FT                   /evidence="ECO:0000305|PubMed:8424780"
                        if checktype == "MOD_RES":
                            pass
                        elif checktype == "BINDING":
                            pass
                        elif checktype == "ACT_SITE":
                            pass
                        elif checktype == "VARIANT":
                            homol_mutat = position
                            if ".." in homol_mutat:	#FT   BINDING         10..17
                                    bindingposone = homol_mutat.split("..")[0]
                                    bindingpostwo = homol_mutat.split("..")[1]
                                    for bindpos in range(int(bindingposone),int(bindingpostwo)+1):
                                        for k in idcollector:
                                            if k not in inputdictionary_two:
                                                inputdictionary_two[k]={}
                                                inputdictionary_two[k]["VARIANT"]=[str(bindpos)]
                                            elif "VARIANT" not in inputdictionary_two[k]:
                                                inputdictionary_two[k]["VARIANT"]=[str(bindpos)]
                                            else:
                                                inputdictionary_two[k]["VARIANT"].append(str(bindpos))
                            else:
                                for k in idcollector:
                                    if k not in inputdictionary_two:
                                        inputdictionary_two[k]={}
                                        inputdictionary_two[k]["VARIANT"]=[str(homol_mutat)]
                                    elif "VARIANT" not in inputdictionary_two[k]:
                                        inputdictionary_two[k]["VARIANT"]=[str(homol_mutat)]
                                    else:
                                        inputdictionary_two[k]["VARIANT"].append(str(homol_mutat))
                            checker = "FALSE"	### set checker to False, so the loop can repeat for another instance
                            checktype = ""
                        elif checktype == "MUTAGEN":
                            homol_mutat = position
                            if ".." in homol_mutat:	#FT   BINDING         10..17
                                    bindingposone = homol_mutat.split("..")[0]
                                    bindingpostwo = homol_mutat.split("..")[1]
                                    for bindpos in range(int(bindingposone),int(bindingpostwo)+1):
                                        for k in idcollector:
                                            if k not in inputdictionary_two:
                                                inputdictionary_two[k]={}
                                                inputdictionary_two[k]["MUTAGEN"]=[str(bindpos)]
                                            elif "MUTAGEN" not in inputdictionary_two[k]:
                                                inputdictionary_two[k]["MUTAGEN"]=[str(bindpos)]
                                            else:
                                                inputdictionary_two[k]["MUTAGEN"].append(str(bindpos))
                            else:
                                for k in idcollector:
                                    if k not in inputdictionary_two:
                                        inputdictionary_two[k]={}
                                        inputdictionary_two[k]["MUTAGEN"]=[str(homol_mutat)]
                                    elif "MUTAGEN" not in inputdictionary_two[k]:
                                        inputdictionary_two[k]["MUTAGEN"]=[str(homol_mutat)]
                                    else:
                                        inputdictionary_two[k]["MUTAGEN"].append(str(homol_mutat))
                            checker = "FALSE"	### set checker to False, so the loop can repeat for another instance
                            checktype = ""
                        else:
                            pass
                    else:
                        #FT   MUTAGEN         14
                        #FT                   /note="G->V: Increased Rho protein signal transduction.
                        #FT                   Constitutively active."
                        #FT                   /evidence="ECO:0000269|PubMed:19948726,
                        #FT                   ECO:0000269|PubMed:31570889"
                        comment_add = newline.split("_")[1:]
                        for item in comment_add:
                            comment = comment + "/" + item
        except:
            print(logging.exception("message"))
            #print(newline)
            #print(moddedline)
            pass



for k in inputdictionary_two:
    for v in inputdictionary_two[k]:
        for i in inputdictionary_two[k][v]:
            print(k, "\t", translator[k][0],"\t",str(translator[k][1]).replace("[","").replace("]","").replace("'",""),"\t", v,"\t",i)


### 17511 human proteins, 292k proteins in total (all reviewed) as of 04. Jan 2025
