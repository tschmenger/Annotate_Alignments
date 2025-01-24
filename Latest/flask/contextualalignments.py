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
import gzip
from flask import Flask, request, jsonify, render_template, send_from_directory
import random
import shutil
app = Flask(__name__)
#########################################################################################################################################################################################################################

idcontainer = []
positioncolors = ["lightgreen","salmon","yellow","orchid","lightblue"]
colors = {}
coloringcategories = []
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
    idsource = alignmentfile.splitlines()
    for line in idsource:
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
    try:
        idcontainer.remove("\n")
    except:
        pass
    try:
        idcontainer.remove("")
    except:
        pass
    return idcontainer
# ---------------------------------------------------------------------------------------------------------------------------------------------
def HOMOL_UNIP_GETTER(idc, inputdictionary_two):
	for k in idc:
		subdirname = k[0:4]
		try:
			filepathname = "uniprotdata/"+subdirname+"/"+str(k)+".txt.gz"
			with gzip.open(filepathname,"r") as infl:
				for oldline in infl:
					line = oldline.decode('utf-8')
					primary_acc = line.split("\t")[0]
					secondary_acc = line.split("\t")[2]
					genus   = line.split("\t")[1].split("_")[0]
					translator[k]=genus.replace(" ","")
					mod_type = line.split("\t")[3].replace(" ","")
					affectedposition = line.split("\t")[4].replace(" ","").replace("\n","")
					if k not in inputdictionary_two:
						inputdictionary_two[k]={}
						inputdictionary_two[k][mod_type]=[affectedposition]
					elif mod_type not in inputdictionary_two[k]:
						inputdictionary_two[k][mod_type]=[affectedposition]
					else:
						inputdictionary_two[k][mod_type].append(affectedposition)
		except:
			#logging.exception("message")
			pass

	return inputdictionary_two, translator
#########################################################################################################################################################################################################################
feature_dict = {}
gapletters = [".","-"]
translator = {}
beginnerdict = {}
def CUSTOM_ALIGN(targetfile):
    alignments_to_keep = {}
    alignfile = targetfile.splitlines()
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
def SHOWORDER(seqs, doi, starti, endi, goi, maxnumber):
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
	raw_ranking = sorted(showtime, key=lambda x: (-showtime[x], x))
	raw_ranking.insert(0,goi)
	try:
		ranking = raw_ranking[0:int(maxnumber)]
	except:
		ranking = raw_ranking
	return ranking
# ---------------------------------------------------------------------------------------------------------------------------------------------
def interprodownloader(identif):
	interpro_processed = {}
	try:
		subdirname = identif[0:4]
		filepathname = "interprodata/"+subdirname+"/"+str(identif)+"_interpro.txt.gz"
		with gzip.open(filepathname,"r") as infl:
			for oldline in infl:
				line = oldline.decode('utf-8')	### Q197B6   BINDING: ATP.   [62, 62, 62, 213, 213, 213]
				unipid = line.split("\t")[0].replace(" ","")
				if unipid == identif:
					categ = line.split("\t")[1].replace(" ","")
					rawpossinterest = line.split("\t")[2].replace(" ","")
					possinterest = str(rawpossinterest).replace("[","").replace("]","").replace(" ","").replace("\n","").split(",")
					if categ not in interpro_processed:
						interpro_processed[categ]=[]
						for stelle in possinterest:
							if int(stelle) not in interpro_processed[categ]:
								interpro_processed[categ].append(int(stelle))
	except:
		#logging.exception("message")
		pass
	return	interpro_processed
# ---------------------------------------------------------------------------------------------------------------------------------------------
def create_svg(sequences_dict, positions, colordict, startposition, windowsize, poi, forbidden, proteinfeatures, Konserve, featurecolors, translator, topgun):
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

    filename 		= poi+"_Position"+str(startposition)+"_Windowsize"+str(windowsize)+"_Sequences"+str(topgun)+".svg"
    filename_print 	= poi+"_Position"+str(startposition)+"_Windowsize"+str(windowsize)+"_Sequences"+str(topgun)+"_print.svg"
    with open(filename_print,"w+") as dwg:
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
        roworder = SHOWORDER(sequences_dict, positions, distance_start, distance_end, poi, topgun)
    #    for uniprot in sequences_dict:
    #	if uniprot not in roworder:
    #		roworder.append(uniprot)

        maximumdistance = distance_end - distance_start
        viewboxcounter = 1
        all_x_vals = []
        highlightingID = 0
        highlightsaver = {}
        konservation_text = []
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
                                                savelist = [x, Konservierungsypsilon, konserv_val]
                                                konservation_text.append(savelist)
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
                                                    savelist = [x, Konservierungsypsilon, konserv_val]
                                                    konservation_text.append(savelist)
                                                    #dwg.add(dwg.rect((x, Konservierungsypsilon), (10, 14), fill="black"))
                                                    toproof = 1-float(konserv_val)
                                                    dwg.write("""<rect fill="white" height='"""+str(14*toproof)+"""' width="10" x='"""+str(x)+"""' y='"""+str(Konservierungsypsilon)+"""' />""")
                                                    #dwg.add(dwg.rect((x, Konservierungsypsilon), (10, 14*toproof), fill="white"))
                                                    if int(startnumberlabel)%10 == False:
                                                        dwg.write("""<text dominant-baseline="central" fill="black" font-family="Arial" font-size="8px" font-weight="bold" text-anchor="middle" x='"""+str(x+5)+"""' y='"""+str(Konservierungsypsilon-3)+"""'>"""+str(startnumberlabel)+"""</text>""")
                                                        #dwg.add(dwg.text(str(startnumberlabel), insert=(x+5, Konservierungsypsilon-3), text_anchor='middle', dominant_baseline='central', font_size='8px', font_family='Arial', font_weight='bold', fill='black'))
                                        elevator_testval = "no"
                                        failcounter = 0
                                        elevator_floor = 0
                                        for feat in proteinfeatures:
                                                if startnumberlabel in proteinfeatures[feat]:
                                                    if feat not in tempfeat:
                                                        tempfeat[feat]=[featurecolors[featcount],featcount]
                                                        featcount+=1
                                                    elevator = tempfeat[feat][1] ## this will be the color
                                                    #elevator_floor = 0
                                                    if elevator >= 8: ## I have a range of 9 colors here
                                                        if elevator_floor <= 10:
                                                            elevator = elevator_floor
                                                            elevator_floor += 1
                                                        else:
                                                            if elevator_testval  == "no":
                                                                elevator_floor = 0
                                                                elevator = elevator_floor
                                                                elevator_testval="yes"

                                                            else:
                                                                elevator = elevator_floor
                                                                elevator_floor += 1
                                                                failcounter += 1
                                                                if failcounter >= 8:
                                                                    elevator_testval = "no"
                                                    #print(elevator, elevator_floor, tempfeat[feat][0], feat, elevator_testval)
                                                    # 2 0 orange GTP/Mg2+ binding site
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
    with open(filename_print,"r") as dwg_file:
        existing_content = dwg_file.read()
    dwg = open(filename_print,"w+")
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

    imagefile = open(filename_print,"r")
    data= imagefile.read()
    data = data.replace("</svg>", styletext+"</svg>")
    imagefile.close()
    writeFile = open(filename, "w")
    writeFile.write(data)
    writeFile.close()

    circletext = ""
    circletext_two = ""
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
    for hltxt in konservation_text:
        cx = hltxt[0]
        cy = hltxt[1]
        txt = str(int(hltxt[2]*100))+"%"
        whiteboxheight = 20
        circletext_two = circletext_two+"""<g xmlns="http://www.w3.org/2000/svg">
          <circle xmlns="http://www.w3.org/2000/svg" cx='"""+str(cx)+"""' cy='"""+str(cy)+"""' r="7" style="fill:transparent;stroke:transparent;stroke-width:0.5;fill-opacity:0.25;stroke-opacity:0.25"/>
          <rect class="hiss" x='"""+str(cx)+"""' y='"""+str(cy-20)+"""' height='"""+str(whiteboxheight)+"""' width='"""+str(30)+"""'></rect>
          <text class="moo" x='"""+str(cx)+"""' y='"""+str(cy)+"""'><tspan class="text">"""+txt+"""</tspan></text></g>"""
    #print(konservation_text)
    imagefile = open(filename,"r")
    imagefile.seek(0)

    data = imagefile.read()
    imagefile.close()
    data = data.replace("</svg>", circletext+"</svg>")
    data = data.replace("</svg>", circletext_two+"</svg>")

    writeFile = open(filename, "w")
    writeFile.write(data)
    writeFile.close()

    return filename
#########################################################################################################################################################################################################################
########################################################################################################################
@app.route('/')
def index():
    return render_template('index.html')  # Stellt die index.html bereit
RESULT_DIR = os.getcwd()
@app.route('/generate', methods=['POST'])
def generate_result():
    # Benutzereingaben aus der Anfrage abrufen
    positions = {}
    translator = {}
    data = request.json
    input1 = data.get("input1")
    input2 = data.get("input2")
    slider = data.get("slider")
    topguns = data.get("topguns")
    alignmentfile = data.get("clustalInput")

    if not input1 or not input2 or not alignmentfile or not topguns:
        return jsonify({"error": "All fields are required!"}), 400


    # Dateinamen erstellen
    filename = f"{input1}_Position{input2}_Windowsize{slider}_Sequences{topguns}.svg"
    filename_print = f"{input1}_Position{input2}_Windowsize{slider}_Sequences{topguns}_print.svg"
    ids 			         = idgetter(alignmentfile)
    sequences 		        = CUSTOM_ALIGN(alignmentfile)
    protein_of_interest 	= input1
    try:
        position_of_interest    = int(input2)
    except Exception as e:
        return jsonify({'Please check your input.': str(e)}), 500
    try:
        window                  = int(slider)
    except Exception as e:
        return jsonify({'Please check your input.': str(e)}), 500
    try:
        topgunner 				= int(topguns)
    except Exception as e:
        return jsonify({'Please check your input.': str(e)}), 500
    positions, translator 	= HOMOL_UNIP_GETTER(ids,positions)
    ################################################################################################################################################
    Konserve, TheForbiddenPositions 	= conservation_checker(protein_of_interest,sequences, positions)
    ############
    try:
        feature_dict = interprodownloader(protein_of_interest)
    except:
        feature_dict = {}
    #print(feature_dict)
    ############

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
    resultfile = create_svg(sequences, positions, colors, position_of_interest, window, protein_of_interest, TheForbiddenPositions, feature_dict, Konserve, featurecolors, translator, topguns)

    RESULT_DIR_two = os.path.join(os.getcwd(), 'static', 'results')

    # move the normal file
    temp_file_path = os.path.join(os.getcwd(), filename)  # Assume initially created in working directory
    resultfile_path = os.path.join(RESULT_DIR_two, filename)
    shutil.move(temp_file_path, resultfile_path)

    # Same process for the print version
    temp_file_path_print = os.path.join(os.getcwd(), filename_print)
    resultfile_print_path = os.path.join(RESULT_DIR_two, filename_print)
    shutil.move(temp_file_path_print, resultfile_print_path)

    # Dateinamen an das Frontend zurückgeben
    return jsonify({"filename": filename, "filename_print": filename_print})


@app.route('/download/<filename>')
def download_file(filename):
    return send_from_directory(RESULT_DIR, filename, as_attachment=True)
@app.route('/example_data', methods=['GET'])
def get_example_data():
    # Define the file paths for the example data
    example_clustal_path = './static/example/RHOA_BlastpExample_ClustalMSA.clustal'

    try:
        # Read the Clustal alignment from the file
        with open(example_clustal_path, 'r') as file:
            clustal_data = file.read()

        # Return example data as JSON
        return jsonify({
            'input1': 'P61586',  # Example Uniprot accession
            'input2': '13',  # Example position
            'slider': '10',  # Example slider value
            'clustalInput': clustal_data,  # Example Clustal alignment
            'svgPath': '/static/example/P61586_Position13_Windowsize10.svg',
            'svgPrintPath': '/static/example/P61586_Position13_Windowsize10_print.svg'
        })

    except Exception as e:
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    app.run(debug=True, use_reloader=False)
