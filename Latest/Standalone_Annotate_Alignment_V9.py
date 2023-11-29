#	# ###### ###### #    # ###### ######
#	# #	 #      ##   # #      #
#	# #      ###### # #  # #####  ######
#	# #      #      #  # #      # #  
#	# #      #      #   ##      # #
######	# ###### ###### #    # ###### ######
# 	This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, 
#	either version 3 of the License, or (at your option) any later version.
# 	This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#	See the GNU General Public License for more details.
# 	You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 
#########################################################################################################################################################################################################################
#!/usr/bin/env/Python 3.6.8
import sys
import urllib.request
import re
nums = re.compile(r"[+-]?\d+(?:\.\d+)?")
whitespace_killer=re.compile(r"\s+")
import time
import ast
import logging
import os
#########################################################################################################################################################################################################################
positions = {}
translator = {}
idcontainer = []
# https://www.jalview.org/help/html/colourSchemes/clustal.html
Clustalcolors = {"A":"hydrophobic",
		"I":"hydrophobic",
		"L":"hydrophobic",
		"M":"hydrophobic",
		"F":"hydrophobic",
		"W":"hydrophobic",
		"V":"hydrophobic",
		"C":"hydrophobic",
		"K":"positive",
		"R":"positive",
		"E":"negative",
		"D":"negative",
		"N":"polar",
		"Q":"polar",
		"S":"polar",
		"T":"polar",
		"G":"glycine",
		"P":"proline",
		"H":"aromatic",
		"Y":"aromatic"}
clustaltypes = {"hydrophobic":"blue",			
		"positive":"red",
		"negative":"magenta",
		"polar":"green",
		"glycine":"black",
		"proline":"orange",
		"aromatic":"cyan"}
#########################################################################################################################################################################################################################
def idgetter(alignmentfile):
	with open(alignmentfile,"r") as idsource:
		for line in idsource:	
### 	CLUSTAL O(1.2.4) multiple sequence alignment
###
###
###	Q9BYZ6      ------------------------------------------------------------	0
			if "CLUSTAL" not in line:
				newline = whitespace_killer.sub(" ",line).replace("\n","")
				try:
					identifier = line.split(" ")[0]
					if "." in identifier:
						newidentifier = identifier.split(".")[0]
						if newidentifier not in idcontainer:
							idcontainer.append(newidentifier)
					elif "|" in identifier:
						newidentifier = identifier.split("|")[1]
						if newidentifier not in idcontainer:
							idcontainer.append(newidentifier)
					else:
						if identifier not in idcontainer:
							idcontainer.append(identifier)
				except:	
					pass
	idcontainer.remove("\n")
	idcontainer.remove("")
	return idcontainer
# ---------------------------------------------------------------------------------------------------------------------------------------------
def HOMOL_UNIP_GETTER(idc, inputdictionary_two):
	for k in idc:
			content=urllib.request.urlopen("https://rest.uniprot.org/uniprotkb/" + k + ".txt")
			checktype = ""
			checker = "FALSE"
			position = ""
			comment = ""
			for oldline in content:
				line = oldline.decode('utf-8')
				try:
					if "ID" in line and "Reviewed;" in line:	### ID   RHOA_HUMAN              Reviewed;         193 AA.
						newline = whitespace_killer.sub(" ",line).replace("\n","")
						langname = newline.split(" ")[1].replace(" ","").replace("\n","")
						if k not in translator:
							translator[k]=langname
					if "FT" and "VARIANT" in line:
						if checker == "FALSE":	### this will happen when I get to such a line:			FT   VARIANT         47
							comment = ""				
							newline = whitespace_killer.sub(" ",line).replace("\n","")
							position = newline.split(" ")[2]
							checker = "TRUE"
							checktype = "VARIANT"
					elif "FT" and "MUTAGEN" in line:
						if checker == "FALSE":	### this will happen when I get to such a line:			FT   MUTAGEN         226
							comment = ""	
							newline = whitespace_killer.sub(" ",line).replace("\n","")
							position = newline.split(" ")[2]
							checker = "TRUE"
							checktype = "MUTAGEN"
					elif "BINDING" in line:	### this will happen when I get to such a line:			FT   BINDING         10..17
						if "FT" in line:
							if checker == "FALSE":	### this will happen when I get to such a line:			FT   MUTAGEN         226
								comment = ""	
								newline = whitespace_killer.sub(" ",line).replace("\n","")
								position = newline.split(" ")[2]
								checker = "TRUE"
								checktype = "BINDING"
					elif "FT" and "MOD_RES" in line:	### this will happen when I get to such a line:			FT   BINDING         38
						if checker == "FALSE":	### this will happen when I get to such a line:			FT   MUTAGEN         226
							comment = ""	
							newline = whitespace_killer.sub(" ",line).replace("\n","")
							position = newline.split(" ")[2]
							checker = "TRUE"
							checktype = "MOD_RES"
					elif "FT" and "ACT_SITE" in line:
						if checker == "FALSE":	### this will happen when I get to such a line:			FT   ACT_SITE        108
							comment = ""	
							newline = whitespace_killer.sub(" ",line).replace("\n","")
							position = newline.split(" ")[2]
							checker = "TRUE"
							checktype = "ACT_SITE"
					else:
						if checker == "TRUE":	### this should happen right on the next line after setting checker to TRUE
							newline = whitespace_killer.sub(" ",line).replace("\n","").replace(" ","_")
						#	FT_/note="G->V:_Increased_Rho_protein_signal_transduction.
						#	FT_Constitutively_ active."
						#	FT_/evidence="ECO:0000269|PubMed:19948726,
	
							if "FT_/note=" in newline:		### additional failsafe
								if checktype == "MUTAGEN":
									#FT   MUTAGEN         34
									#FT                   /note="Y->A: Abolishes interaction with DGKQ."
									#FT                   /evidence="ECO:0000269|PubMed:10066731,
									if "Missing" in newline:
										pass
									else:
										res_one = newline.split(":")[0].split("=\"")[1].split("->")[0]
										res_two = newline.split(":")[0].split("=\"")[1].split("->")[1]
										comment = newline.split(":")[1]
								elif 	checktype == "VARIANT":
									#FT   VARIANT         212
									#FT                   /note="D -> V (in dbSNP:rs56143363)"
									#FT                   /evidence="ECO:0000269|PubMed:17344846"
									#FT                   /id="VAR_040389"
									if "Missing" in newline:
										pass
									else:
										res_one = newline.split("_")[1].split("=\"")[1]
										res_two = newline.split("_")[3]
										comment_raw = newline.split("_")[4:]
										for item in comment_raw:
											comment = comment + item + "_"
										comment = comment + "/"	
								elif checktype == "MOD_RES":
									#FT   MOD_RES         689
									#FT                   /note="Phosphoserine"
									#FT                   /evidence="ECO:0007744|PubMed:23186163"
									if "Missing" in newline:
										pass
									else:
										binder = newline.split("=\"")[1]						
										comment = binder[:-1].replace("\"","")

										checker = "FALSE"
										homol_mutat = position
										kommentar = checktype+"/"+comment
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
									if "Missing" in newline:
										pass
									else:
										binder = newline.split("=\"")[1]						
										comment = binder[:-1].replace("\"","")
										checker = "FALSE"
										homol_mutat = position
										if ".." in homol_mutat:	#FT   BINDING         10..17
											bindingposone = homol_mutat.split("..")[0]
											bindingpostwo = homol_mutat.split("..")[1]
											for bindpos in range(int(bindingposone),int(bindingpostwo)+1):
												if k not in inputdictionary_two:
													inputdictionary_two[k]={}
													inputdictionary_two[k]["ACT_SITE"]=[str(bindpos)]
												elif "ACT_SITE" not in inputdictionary_two[k]:
													inputdictionary_two[k]["ACT_SITE"]=[str(bindpos)]
												else:
													inputdictionary_two[k]["ACT_SITE"].append(str(bindpos))
										else:
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
												if k not in inputdictionary_two:
													inputdictionary_two[k]={}
													inputdictionary_two[k]["BINDING"]=[str(bindpos)]
												elif "BINDING" not in inputdictionary_two[k]:
													inputdictionary_two[k]["BINDING"]=[str(bindpos)]
												else:
													inputdictionary_two[k]["BINDING"].append(str(bindpos))
										else:
											homol_mutat = position
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
												if k not in inputdictionary_two:
													inputdictionary_two[k]={}
													inputdictionary_two[k]["VARIANT"]=[str(bindpos)]
												elif "VARIANT" not in inputdictionary_two[k]:
													inputdictionary_two[k]["VARIANT"]=[str(bindpos)]
												else:
													inputdictionary_two[k]["VARIANT"].append(str(bindpos))
									else:
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
												if k not in inputdictionary_two:
													inputdictionary_two[k]={}
													inputdictionary_two[k]["MUTAGEN"]=[str(bindpos)]
												elif "MUTAGEN" not in inputdictionary_two[k]:
													inputdictionary_two[k]["MUTAGEN"]=[str(bindpos)]
												else:
													inputdictionary_two[k]["MUTAGEN"].append(str(bindpos))
									else:
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
					pass
	return inputdictionary_two, translator
#########################################################################################################################################################################################################################
feature_dict = {}
gapletters = [".","-"]
translator = {}
beginnerdict = {}
def CUSTOM_ALIGN(targetfile):
	alignments_to_keep = {}
	with open(targetfile,"r") as alignfile:
            for line in alignfile:
                if "CLUSTAL" not in line:
                    newline = whitespace_killer.sub(" ",line).replace("\n","")
                    found_seqname = "no"
                    if len(newline.strip())!= 0:
                        try:
                            seq_name = newline.split(" ")[0].split("_")[0]
                            if "." in seq_name:
                                seq_name = newline.split(" ")[0].split("_")[0].split(".")[0]
                                found_seqname = "yes"
                            elif "|" in seq_name:
                                seq_name = newline.split(" ")[0].split("_")[0].split("|")[1]
                                found_seqname = "yes"
                            else:
                                found_seqname = "yes"   
                        except:
                            seq_name = newline.split(" ")[0]

                        if found_seqname == "no":
                            try:
                                seq_name = newline.split(" ")[0].split(".")[0]
                                found_seqname = "yes"
                            except:
                                seq_name = newline.split(" ")[0]

                        if found_seqname == "no":
                            try:
                                seq_name = newline.split(" ")[0].split("|")[1]
                                found_seqname = "yes"
                            except:
                                seq_name = newline.split(" ")[0]

                        seq = newline.split(" ")[1]
                        if seq_name != "":
                            if seq_name not in alignments_to_keep:
                                alignments_to_keep[seq_name]=str(seq)
                            else:
                                alignments_to_keep[seq_name]=alignments_to_keep[seq_name]+str(seq)
#	print(alignments_to_keep)
	return alignments_to_keep
# ---------------------------------------------------------------------------------------------------------------------------------------------
def conservation_checker(identifier, seqdict, relevantpositions):
	for protein in seqdict:
		sequence = seqdict[protein]
		if identifier in protein:
			truesequenzler = sequence
	sequencelength = len(truesequenzler) ### including gaps, meaning this is the alignment length
	positionalcounter = 1
	conservational_dictionary = {}
	theforbiddenalignpos = []
	for i in range(0,sequencelength):
		keepITin = "false"
		### i corresponds to the alignment position!
		identitycontainer = []
		for ident in seqdict:
			if identifier in ident:
				orires = seqdict[ident][i]
				identitycontainer.append(seqdict[ident][i].upper())
			else:
				identitycontainer.append(seqdict[ident][i].upper())
		identitypercentage = float(identitycontainer.count(orires.upper()))/float(len(identitycontainer))	### so far this also includes "-" as the original truesequence residue, be cautious
		oritype = "none"

		if orires not in gapletters:
			#print positionalcounter,"\t", identitycontainer, "\t", orires,"\t",identitypercentage,"\t",
			if int(positionalcounter) not in conservational_dictionary:
				conservational_dictionary[int(positionalcounter)] = [float(identitypercentage), orires]
			positionalcounter+=1
		elif orires in gapletters:
			if identitypercentage >= 0.90:
				theforbiddenalignpos.append(i+1)		
		else:
			pass
	return conservational_dictionary, theforbiddenalignpos	
# ---------------------------------------------------------------------------------------------------------------------------------------------
def SHOWORDER(seqs, doi, starti, endi, goi):
	# dictionary of interest, residue of interest, windowsize, gene of interest
	showtime = {}
	for k in doi:	### uniprot ID = k
		featurecount = []
		sequenzler = seqs[k]
		residue = 0
		for i, letter in enumerate(sequenzler,start = 1):
			if letter not in gapletters:
				residue += 1
				if i >= starti:
					if i <= endi:
						for v in doi[k]: ### categories, i.e. VARIANT = v
							for vv in doi[k][v]:	### residue number = vv
								try:
									if int(vv) == residue:
										if int(vv) not in featurecount:
											featurecount.append(vv)
								except:
									pass
		if k != goi:							
			if k not in showtime:
				showtime[k]=len(featurecount)
	#print showtime
	ranking = sorted(showtime, key=lambda x: (-showtime[x], x))
	ranking.insert(0,goi)
	#print ranking
	#for ranked in ranking:
	#	print ranked, "\t", showtime[ranked]
	return ranking
# ---------------------------------------------------------------------------------------------------------------------------------------------
def interprodownloader(identif):
	import ast
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
					if descr not in interpro_processed:
						interpro_processed[descr]=[residue]
					else:
						if residue not in interpro_processed[descr]:
							interpro_processed[descr].append(residue)
	return	interpro_processed
# ---------------------------------------------------------------------------------------------------------------------------------------------
def create_svg(sequences_dict, positions, colordict, startposition, windowsize, poi, forbidden, proteinfeatures):
    heatmapper = {}
    startposition_checker = startposition
    lengeforplotting = len(coloringcategories)
    if lengeforplotting < 4:
        lengeforplotting = 4
    Heatmapstart = 60-(len(coloringcategories)+1)*10
    Konservierungsypsilon = Heatmapstart - 20
    Categoryypsilon = Heatmapstart - 130	
    #### do this when havng constructed the dictionary with interesting positions
    #### here it is supplied as is, but needs to be further modified
    for item in positions:
        for categ in colordict:
            if categ not in positions[item]:
                positions[item][categ]=[]
    if startposition == "none":
        startposition = 1
    try:
        filename = translator[poi]+"_Position"+str(startposition)+"_Windowsize"+str(windowsize)+".svg"
    except:
        filename = poi+"_Position"+str(startposition)+"_Windowsize"+str(windowsize)+".svg"
    with open(filename,"w+") as dwg:
    #dwg = open(filename,"a")
        x = 50
        y = 80
        sequence_of_interest = sequences_dict[poi]
        non_minus_count = 0
        distance_end = len(sequence_of_interest)+100	### to make sure it gets weeded out below, if none of the if statements directly below trigger
        distance_start = 0					### to make sure it gets weeded out below, if none of the if statements directly below trigger
        for i, letter in enumerate(sequence_of_interest,start = 1):
            if letter not in gapletters:
                non_minus_count += 1
                if non_minus_count == startposition:
                    startpos = i	### this is the alignment position that corresponds to the residue of interest. alignment position includes "-"
                if non_minus_count == startposition+windowsize:
                    distance_end = i
                if non_minus_count == startposition-windowsize:
                    distance_start = i
        maxcharactercnt = non_minus_count		### should capture the true length of the sequence of interest
        ### make sure the windowsize does not conflict with positions close to the start or end of the sequence
        if distance_start <= 0:
            distance_start = 1
        if distance_end > len(sequence_of_interest):
            distance_end = len(sequence_of_interest)
        roworder = SHOWORDER(sequences_dict, positions, distance_start, distance_end, poi)
    #    for uniprot in sequences_dict:
    #	if uniprot not in roworder:
    #		roworder.append(uniprot)

        maximumdistance = distance_end - distance_start
        viewboxcounter = 1
        all_x_vals = []
        highlightingID = 0
        highlightsaver = {}
        print("Drawing alignment...")
        #print(roworder)
        for uniprot in roworder:
            seq 	= sequences_dict[uniprot]
            namus 	= uniprot
            ###print(namus,"\t",seq)
            startingpoint = startposition - windowsize	### this is required for the correct labeling according to the sequence of interest
            try:
                drawname = translator[namus]
            except:
                drawname = namus


            if poi.replace(" ","") in namus.replace(" ",""):
                #print(poi,"\t",drawname,"true")
                old_x = x
                old_y = y
                x = 50
                y = 60
                dwg.write("<rect fill='yellow' height='14' width='90' x='"+str(x-100)+"' y='"+str(y)+"' />")
                #dwg.add(dwg.rect((x-100, y), (90, 14), fill="yellow"))
                if len(drawname) < 8:
                    dwg.write("<text dominant-baseline='central' fill='black' font-family='Arial' font-size='10px' font-weight='bold' text-anchor='end' x='"+str(x-45)+"' y='"+str(y+7)+"'>"+str(drawname)+"</text>")
                    #dwg.add(dwg.text(drawname, insert = (x-45,y+7), text_anchor='end', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill='black'))
                    dwg.write("""<text dominant-baseline="central" fill="black" font-family="Arial" font-size="10px" font-weight="bold" text-anchor="end" x='"""+str(x-45)+"""' y='"""+str(Konservierungsypsilon+5)+"""'>"""+drawname+"""</text>""")
                    #dwg.add(dwg.text(drawname, insert = (x-45,Konservierungsypsilon+5), text_anchor='end', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill='black'))
                else:
                    dwg.write("""<text dominant-baseline="central" fill="black" font-family="Arial" font-size="7px" font-weight="bold" text-anchor="end" x='"""+str(x-45)+"""' y='"""+str(y+7)+"""'>"""+drawname+"""</text>""")
                    #dwg.add(dwg.text(drawname, insert = (x-45,y+7), text_anchor='end', dominant_baseline='central', font_size='7px', font_family='Arial', font_weight='bold', fill='black'))
                    dwg.write("""<text dominant-baseline="central" fill="black" font-family="Arial" font-size="7px" font-weight="bold" text-anchor="end" x='"""+str(x-45)+"""' y='"""+str(Konservierungsypsilon+5)+"""'>"""+drawname+"""</text>""")
                    #dwg.add(dwg.text(drawname, insert = (x-45,Konservierungsypsilon+5), text_anchor='end', dominant_baseline='central', font_size='7px', font_family='Arial', font_weight='bold', fill='black'))
            else:
                #print(poi,"\t",drawname,"false")
                if len(drawname) < 8:
                    dwg.write("""<text dominant-baseline="central" fill="black" font-family="Arial" font-size="10px" font-weight="bold" text-anchor="end" x='"""+str(x-45)+"""' y='"""+str(y+7)+"""'>"""+drawname+"""</text>""")
                    #dwg.add(dwg.text(drawname, insert = (x-45,y+7), text_anchor='end', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill='black'))
                else:
                    dwg.write("""<text dominant-baseline="central" fill="black" font-family="Arial" font-size="7px" font-weight="bold" text-anchor="end" x='"""+str(x-45)+"""' y='"""+str(y+7)+"""'>"""+drawname+"""</text>""")
                    #dwg.add(dwg.text(drawname, insert = (x-45,y+7), text_anchor='end', dominant_baseline='central', font_size='7px', font_family='Arial', font_weight='bold', fill='black'))
            #charactercount = 0
            totalcount = 0
            if startingpoint <= 0:
                startnumberlabel = 1
            elif startingpoint >= maxcharactercnt:
                startnumberlabel = maxcharactercnt
            else:
                startnumberlabel = startingpoint 
            charactercount = 0
            tempfeat = {}
            featcount = 0
            firstdone = "false"
            lastdone = "false"
            forbidden_start = "false"
            forbidden_end = "false"
            gapcounter = 0
            for i, letter in enumerate(seq, start=1):
                totalcount += 1		#### gives the alignment position, including gaps
                letter = seq[i-1]
                if x not in all_x_vals:
                    all_x_vals.append(x)
                if letter not in gapletters:
                    charactercount += 1
                    if totalcount <= distance_end:	### distance_end refers to the last alignment position that will be considered, which is +windowsize non-gap residues from the input position
                        endcounter = charactercount
                        testlenge = int(distance_end)-int(totalcount)
                        if testlenge <= maximumdistance:	### checks that we still operate around the position of interest +/- residues only
                            if totalcount >= distance_start:
                                if firstdone == "false":
                                    forbidden_start = "true"
                                    startcounter = charactercount
                                    dwg.write("""<text dominant-baseline="central" fill="black" font-family="Arial" font-size="10px" font-weight="bold" text-anchor="end" x="35" y='"""+str(y+8)+"""'>"""+str(startcounter)+"""</text>""")
                                    #dwg.add(dwg.text(startcounter, insert=(35, y+8), text_anchor='end', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill="black"))					
                                    firstdone = "true"
                                if totalcount not in forbidden: ### totalcount is int and forbidden is a list of ints
                                    if poi in namus:
                                        konserv_val = Konserve[startnumberlabel][0]
                                        if float(konserv_val)>= 0.7:
                                            dwg.write("""<rect fill='"""+(clustaltypes[Clustalcolors[letter.upper()]])+"""' height='"""+str(len(roworder)*20)+"""' opacity="0.2" width="10" x='"""+str(x)+"""' y='"""+str(y)+"""' />""")
                                            #dwg.add(dwg.rect((x,y),(10,len(roworder)*20), fill= clustaltypes[Clustalcolors[letter.upper()]], opacity=0.2))
                                        viewboxcounter += 1

                                        if int(startnumberlabel) == int(startposition):
                                                position_interest_x = x
                                                position_interest_y = y
                                                dwg.write("""<rect fill="black" height="14" width="10" x='"""+str(x)+"""' y='"""+str(Konservierungsypsilon)+"""' />""")
                                                #dwg.add(dwg.rect((x, Konservierungsypsilon), (10, 14), fill="black"))
                                                toproof = 1-float(konserv_val)
                                                dwg.write("""<rect fill="white" height='"""+str(14*toproof)+"""' width="10" x='"""+str(x)+"""' y='"""+str(Konservierungsypsilon)+"""' />""")
                                                #dwg.add(dwg.rect((x, Konservierungsypsilon), (10, 14*toproof), fill="white"))
                                                if startposition_checker != "none":
                                                    dwg.write("""<text dominant-baseline="central" fill="red" font-family="Arial" font-size="8px" font-weight="bold" text-anchor="middle" x='"""+str(x+5)+"""' y='"""+str(Konservierungsypsilon-3)+"""'>"""+str(startnumberlabel)+"""</text>""")
                                                    #dwg.add(dwg.text(str(startnumberlabel), insert=(x+5, Konservierungsypsilon-3), text_anchor='middle', dominant_baseline='central', font_size='8px', font_family='Arial', font_weight='bold', fill='red'))			
                                        else:
                                                if int(startnumberlabel)>= int(startingpoint):
                                                    dwg.write("""<rect fill="black" height="14" width="10" x='"""+str(x)+"""' y='"""+str(Konservierungsypsilon)+"""' />""")
                                                    #dwg.add(dwg.rect((x, Konservierungsypsilon), (10, 14), fill="black"))
                                                    toproof = 1-float(konserv_val)
                                                    dwg.write("""<rect fill="white" height='"""+str(14*toproof)+"""' width="10" x='"""+str(x)+"""' y='"""+str(Konservierungsypsilon)+"""' />""")
                                                    #dwg.add(dwg.rect((x, Konservierungsypsilon), (10, 14*toproof), fill="white"))
                                                    if int(startnumberlabel)%10 == False:
                                                        dwg.write("""<text dominant-baseline="central" fill="black" font-family="Arial" font-size="8px" font-weight="bold" text-anchor="middle" x='"""+str(x+5)+"""' y='"""+str(Konservierungsypsilon-3)+"""'>"""+str(startnumberlabel)+"""</text>""")
                                                        #dwg.add(dwg.text(str(startnumberlabel), insert=(x+5, Konservierungsypsilon-3), text_anchor='middle', dominant_baseline='central', font_size='8px', font_family='Arial', font_weight='bold', fill='black'))
                                        for feat in proteinfeatures:
                                                if startnumberlabel in proteinfeatures[feat]:
                                                    if feat not in tempfeat:
                                                        tempfeat[feat]=[featurecolors[featcount],featcount]
                                                        featcount+=1
                                                    elevator = tempfeat[feat][1]
                                                    elevator_floor = 0
                                                    if elevator >= 8:
                                                        if elevator_floor <= 10:
                                                            elevator = elevator_floor
                                                            elevator_floor += 1
                                                        else:
                                                            elevator_floor = 0
                                                            elevator = elevator_floor
                                                    y_level = -45 + (elevator*3)
                                                    y_level_text = -95 + (elevator*5)
                                                    dwg.write("""<rect fill='"""+tempfeat[feat][0]+"""' height="2" width="10" x='"""+str(x)+"""' y='"""+str(y_level)+"""' />""")
                                                    #dwg.add(dwg.rect((x, y_level), (10, 2), fill=tempfeat[feat][0]))
                                                    if "done" not in tempfeat[feat]:
                                                        dwg.write("""<text dominant-baseline="central" fill='"""+tempfeat[feat][0]+"""' font-family="Arial" font-size="6px" font-weight="bold" text-anchor="start" x='"""+str(x+15)+"""' y='"""+str(y_level_text)+"""'>"""+feat+"""</text>""")
                                                        #dwg.add(dwg.text(str(feat), insert=(x+15, y_level_text), text_anchor='start', dominant_baseline='central', font_size='6px', font_family='Arial', font_weight='bold', fill=tempfeat[feat][0]))
                                                        tempfeat[feat].append("done")			
                                                            


                                        startnumberlabel+=1
                                    try:
                                        drawn = 0
                                        radius = 7
                                        for colorcateg in coloringcategories:
                                            if str(charactercount) in positions[namus][colorcateg]:
                                                if drawn != 1:
                                                    highlightingID += 1
                                                    hightlightstring = namus+"/"+letter+str(charactercount) + "|" + colorcateg
                                                else:
                                                    hightlightstring = hightlightstring + "}" + colorcateg
                                                dwg.write("""<circle cx='"""+str(x+5)+"""' cy='"""+str(y+7.5)+"""' fill='"""+colordict[colorcateg]+"""' r='"""+str(radius)+"""' />""")
                                                #dwg.add(dwg.circle((x+5, y+7.5), (radius), fill=colordict[colorcateg]))
                                                dwg.write("""<text dominant-baseline="central" fill="black" font-family="Arial" font-size="10px" font-weight="bold" text-anchor="middle" x='"""+str(x+5)+"""' y='"""+str(y+8)+"""'>"""+letter+"""</text>""")
                                                #dwg.add(dwg.text(letter, insert=(x+5, y+8), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill="black"))
                                                drawn = 1
                                                if x not in heatmapper:
                                                    heatmapper[x]={}
                                                    heatmapper[x][colorcateg]=1
                                                elif colorcateg not in heatmapper[x]:
                                                    heatmapper[x][colorcateg]=1
                                                else:
                                                    heatmapper[x][colorcateg]+=1
                                            radius -= 1
                                        if drawn == 1:
                                            if str(highlightingID) not in highlightsaver:						
                                                highlightsaver[str(highlightingID)]=[x+5,y+7.5,hightlightstring]
                                        if drawn == 0:
                                            dwg.write("""<text dominant-baseline="central" fill="black" font-family="Arial" font-size="10px" font-weight="bold" text-anchor="middle" x='"""+str(x+5)+"""' y='"""+str(y+8)+"""'>"""+letter+"""</text>""")
                                            #dwg.add(dwg.text(letter, insert=(x+5, y+8), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill="black"))
                                    except:
                                        dwg.write("""<text dominant-baseline="central" fill="black" font-family="Arial" font-size="10px" font-weight="bold" text-anchor="middle" x='"""+str(x+5)+"""' y='"""+str(y+8)+"""'>"""+letter+"""</text>""")
                                        #dwg.add(dwg.text(letter, insert=(x+5, y+8), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill="black"))	
                                    x += 10
                            else:
                                gapcounter += 1
                else:	### will draw just a "-" for a gap in the alignment
                    if totalcount >= distance_start:
                        if totalcount <= distance_end:
                            if totalcount not in forbidden:
                                x += 10
            viewboxcounter = x
            lastx = x
            lasty = y
            finalresidue = startcounter+gapcounter+(2*windowsize)
            dwg.write("""<text dominant-baseline="central" fill="black" font-family="Arial" font-size="10px" font-weight="bold" text-anchor="middle" x='"""+str(lastx+20)+"""' y='"""+str(y+8)+"""'>"""+str(endcounter)+"""</text>""")
            #dwg.add(dwg.text(endcounter, insert=(lastx+20, y+8), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill="black"))	
                            
            if poi in namus:
                dwg.write("""<rect fill="none" stroke="black" stroke_width="1" height="14" width='"""+str(x+60)+"""' x="-50" y='"""+str(y)+"""' />""")
                #dwg.add(dwg.rect((-50,y),(x+60,14), fill="none",stroke="black",stroke_width=1))	
                x = 50
                y = old_y
            else:
                x = 50
                y += 20
        dwg.flush()
        dwg.close()
    with open(filename,"r") as dwg_file:
        existing_content = dwg_file.read()
    dwg = open(filename,"w+")
    viewboxwidth = (viewboxcounter+140)
    viewboxheight = len(roworder)*20+100+(200-Categoryypsilon)
    viewbox_dimensions = str(str(-50)+','+str(Categoryypsilon-80)+','+str(viewboxwidth)+','+str(viewboxheight))
    combined_content = """<?xml version="1.0" encoding="utf-8" ?><svg baseProfile="full" height="100%" version="1.1" viewBox='"""+viewbox_dimensions+"""' width="100%" xmlns="http://www.w3.org/2000/svg" xmlns:ev="http://www.w3.org/2001/xml-events" xmlns:xlink="http://www.w3.org/1999/xlink"><defs />""" + existing_content
    dwg.write(combined_content)
    #dwg.viewbox(-50, Categoryypsilon-80,viewboxwidth,viewboxheight)

    if startposition_checker != "none":
        dwg.write("""<rect fill="none" stroke="black" stroke_width="1" height='"""+str(len(roworder)*20)+"""' width="10" x='"""+str(position_interest_x)+"""' y='"""+str(position_interest_y)+"""' />""")
    	#dwg.add(dwg.rect((position_interest_x, position_interest_y), (10, len(roworder)*20),fill="none",stroke="black",stroke_width=1))

    x = 50
    y = 0

    maxfinder = {}

    for xval in heatmapper:
        for category in colors:
             if category not in heatmapper[xval]:
                 heatmapper[xval][category]=0
        for categ in heatmapper[xval]:
            if categ not in maxfinder:
                maxfinder[categ]=[int(heatmapper[xval][categ])]
            else:
                maxfinder[categ].append(int(heatmapper[xval][categ]))
    for allxval in all_x_vals:
        if allxval not in heatmapper:
            heatmapper[allxval]={}
            for category in colors:
                if category not in heatmapper[allxval]:
                    heatmapper[allxval][category]=0.0
    mapx = 40
    mapy = Heatmapstart
    #print heatmapper
    catcounter = 0
    for category in colors:
        catcounter += 1
        try:
            heatmap_maximum = max(maxfinder[category])
        except:
            heatmap_maximum = 1
        dwg.write("""<text dominant-baseline="central" fill="black" font-family="Arial" font-size="10px" font-weight="bold" text-anchor="end" x='20' y='"""+str(mapy+5)+"""'>"""+str(category)+"""</text>""")
        #dwg.add(dwg.text(category, insert=(20, mapy+5), text_anchor='end', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill='black'))
        for xval in heatmapper:
		#print xval, "\t", mapy
            try:
                opac = float(heatmapper[xval][category])/float(heatmap_maximum)
            except:
                opac = 0.0
            if float(opac) == 0.0:
                dwg.write("""<rect fill="lightblue" opacity='"""+str(opac)+"""' height="10" width="10" x='"""+str(xval)+"""' y='"""+str(mapy)+"""' />""")
                #dwg.add(dwg.rect((xval, mapy), (10, 10), fill="lightblue", opacity = 0.15 ))
            else:
                dwg.write("""<rect fill='"""+colors[category]+"""' opacity='"""+str(opac)+"""' height="10" width="10" x='"""+str(xval)+"""' y='"""+str(mapy)+"""' />""")
                #dwg.add(dwg.rect((xval, mapy), (10, 10), fill=colors[category], opacity = opac ))
            if mapy == 20:
                pass
        dwg.write("""<rect fill="none" stroke="black" stroke_width="0.5" height="10" width='"""+str(lastx-mapx-10)+"""' x='"""+str(50)+"""' y='"""+str(mapy)+"""' />""")
        #dwg.add(dwg.rect((50, mapy), (lastx-mapx-10, 10),fill="none",stroke="black",stroke_width=0.5))	### <<<<
        mapy += 10

    for i in range(50,lastx-10,10):
        correct_height_to_draw = catcounter*10
        dwg.write("""<rect fill="none" stroke="black" stroke_width="0.5" height='"""+str(correct_height_to_draw)+"""' width="10" x='"""+str(i)+"""' y='"""+str(Heatmapstart)+"""' />""")
        #dwg.add(dwg.rect((i, Heatmapstart), (10, 40),fill="none",stroke="black",stroke_width=0.5))

    x = 50
    y = 0
    for category in colors:
        dwg.write("""<rect fill='"""+colors[category]+"""' height="10" width="60" x='"""+str(x-30)+"""' y='"""+str(Categoryypsilon)+"""' />""")
        #dwg.add(dwg.rect((x-30, Categoryypsilon), (60, 10), fill=colors[category]))
        dwg.write("""<text dominant-baseline="central" fill="black" font-family="Arial" font-size="10px" font-weight="bold" text-anchor="middle" x='"""+str(x)+"""' y='"""+str(Categoryypsilon+5)+"""'>"""+category+"""</text>""")
        #dwg.add(dwg.text(category, insert=(x, Categoryypsilon+5), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill='black'))
        x += 60

    dwg.write("</svg>")
    dwg.close()
    
    styletext = """<style>
   		<![CDATA[
    		text.moo {
         		font-family: "arial";
         		fill: black;
         		font-size: 100%;
    			}
    		rect.hiss {
         		fill:white;
    			}
   			]]>
   		svg text.moo {display: none;}
   		svg rect.hiss {display: none;}
   		svg g:hover text {display: block;}
   		svg g:hover rect {display: block;}
 		</style>"""

    imagefile = open(filename,"r")
    data= imagefile.read()
    data = data.replace("</svg>", styletext+"</svg>")
    imagefile.close()
    writeFile = open(filename, "w")
    writeFile.write(data)
    writeFile.close()

    circletext = ""
    for hlid in highlightsaver:
        cx = highlightsaver[hlid][0]
        cy = highlightsaver[hlid][1]
        txt = highlightsaver[hlid][2]

        uppertext = txt.split("|")[0]
        lowertext = txt.split("|")[1]
        delty = len(lowertext.split("}"))*10
        tspany = cy+15
        whiteboxheight = len(lowertext.split("}"))*20+30
        tspanner = ""
        for showfeature in lowertext.split("}"):
            tspanner = tspanner + """<text class="moo" x='"""+str(cx)+"""' y='"""+str(tspany-28-delty)+"""'><tspan class="text">"""+str(showfeature)+"""</tspan></text>"""
            tspany += 15

        circletext = circletext+"""<g xmlns="http://www.w3.org/2000/svg">
          <circle xmlns="http://www.w3.org/2000/svg" cx='"""+str(cx)+"""' cy='"""+str(cy)+"""' r="7" style="fill:transparent;stroke:transparent;stroke-width:0.5;fill-opacity:0.25;stroke-opacity:0.25"/>      
          <rect class="hiss" x='"""+str(cx-5)+"""' y='"""+str(cy-40-delty)+"""' height='"""+str(whiteboxheight)+"""' width='"""+str(len(uppertext)+90)+"""'></rect>
          <text class="moo" x='"""+str(cx)+"""' y='"""+str(cy-28-delty)+"""'><tspan class="text">"""+uppertext+"""</tspan></text>"""+tspanner+"""</g>"""

    imagefile = open(filename,"r")
    imagefile.seek(0)

    data = imagefile.read()
    imagefile.close()
    data = data.replace("</svg>", circletext+"</svg>")

    writeFile = open(filename, "w")
    writeFile.write(data)
    writeFile.close()
    
    return filename
#########################################################################################################################################################################################################################
###	python3 Annotate_Alignment_V6.py P61586 34 30 RHOA_BlastpExample_ClustalMSA.clustal RHOA_Blastp_info.txt Features_RHOA.txt
if len(sys.argv) != 7:
	print("""Please provide the correct input parameters.\nThe command should look like this:\npython3 Standalone_Annotate_Alignment_V8.py Uniprot_Accession Position Viewing_Range Alignmentfile.clustal Annotations[default=none] Structural_Features[default=none]""")	
	exit()

alignmentfile 		= sys.argv[4]
ids 			= idgetter(alignmentfile)
sequences 		= CUSTOM_ALIGN(alignmentfile)
protein_of_interest 	= sys.argv[1]
try:
    position_of_interest = int(sys.argv[2])
except:
    position_of_interest = str(sys.argv[2])
window = int(sys.argv[3])
if position_of_interest == "none":
    window = 30000

print("Your process is running. This may take a few moments...")
if str(sys.argv[5]) == "none":
	try:
                print("Fetching information from Uniprot...")
                positions, translator 	= HOMOL_UNIP_GETTER(ids,positions)
	except TimeoutError:
		print("Uniprot timeout. Will try again in 10 secs.")
		time.sleep(10)
		print("Attempting to fetch information from Uniprot...")
		positions, translator 	= HOMOL_UNIP_GETTER(ids,positions)
	except:
		#logging.exception("message")
		pass
else:
	with open(sys.argv[5]) as f:
    		data_align = f.read()
	positions = ast.literal_eval(data_align)
############
Konserve, TheForbiddenPositions 	= conservation_checker(protein_of_interest,sequences, positions)
############
try:
    if str(sys.argv[6]) == "none":
        print("Fetching structural information from Interpro...")
        feature_dict = interprodownloader(protein_of_interest)
    else:
        with open(sys.argv[6]) as ff:
            data_feat = ff.read()
        feature_dict = ast.literal_eval(data_feat)
except:
    #logging.exception("message")
    pass
############
positioncolors = ["lightgreen","salmon","yellow","orchid","lightblue"]
colors = {}
coloringcategories = []
counter = 0
############
for k in positions:
    for v in positions[k]:
        if v not in colors:
            if v not in coloringcategories:
                coloringcategories.append(v)
                colors[v]=positioncolors[counter]
                counter+=1

for seqident in sequences:
    if seqident not in positions:
        positions[seqident]={}
        for colcateg in colors:
            positions[seqident][colcateg]=[]
featurecolors = ["firebrick","tomato","orange","olive","palegreen","teal","dodgerblue","blueviolet","deeppink",
		"firebrick","tomato","orange","olive","palegreen","teal","dodgerblue","blueviolet","deeppink",
		"firebrick","tomato","orange","olive","palegreen","teal","dodgerblue","blueviolet","deeppink",
		"firebrick","tomato","orange","olive","palegreen","teal","dodgerblue","blueviolet","deeppink",
		"firebrick","tomato","orange","olive","palegreen","teal","dodgerblue","blueviolet","deeppink",
		"firebrick","tomato","orange","olive","palegreen","teal","dodgerblue","blueviolet","deeppink",
		"firebrick","tomato","orange","olive","palegreen","teal","dodgerblue","blueviolet","deeppink",
		"firebrick","tomato","orange","olive","palegreen","teal","dodgerblue","blueviolet","deeppink",
		"firebrick","tomato","orange","olive","palegreen","teal","dodgerblue","blueviolet","deeppink"]
############
resultfile = create_svg(sequences, positions, colors, position_of_interest, window, protein_of_interest, TheForbiddenPositions, feature_dict)

if str(sys.argv[5]) == "none":
    position_dict_file = resultfile.replace("Position"+str(position_of_interest)+"_Windowsize"+str(window)+".svg","positional_info.txt")
    with open(position_dict_file,'w') as data: 
        data.write(str(positions))
    print("You will find the collected positional information here:",position_dict_file)
if str(sys.argv[6]) == "none":
    structure_dict_file = resultfile.replace("Position"+str(position_of_interest)+"_Windowsize"+str(window)+".svg","structural_info.txt")
    with open(structure_dict_file,'w') as data: 
        data.write(str(feature_dict))
    print("You will find the collected structural information here:",structure_dict_file)
print("Complete! --> ",resultfile)
### python3 Standalone_Annotate_Alignment_V7.py P61586 34 30 ../RHOA_BlastpExample_ClustalMSA.clustal none none





















