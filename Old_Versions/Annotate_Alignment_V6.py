#!/usr/bin/env/Python 3.6.8
import sys
from Bio import AlignIO
sys.path.append('/home/bq_tschmenger/svgwrite/')	### https://pypi.org/project/svgwrite/
import svgwrite	### preferably install it using pip
#from Bio import AlignIO
import ast
import logging
import re
nums = re.compile(r"[+-]?\d+(?:\.\d+)?")
whitespace_killer=re.compile(r"\s+")
import math
import heapq
### THIS VERSION WORKS FOR PYTHON3

feature_dict = {}
gapletters = [".","-"]
translator = {}
beginnerdict = {}
def CUSTOM_ALIGN(targetfile):
	alignments_to_keep = {}
	with open(targetfile,"r") as alignfile:
		for record in AlignIO.read(alignfile,"clustal"):
			try:
				seq_name = record.id.split("_")[1]
			except:
				seq_name = record.id
			try:
				seq_name = record.id.split(".")[0]
			except:
				seq_name = record.id

			seq = record.seq
			if seq_name not in alignments_to_keep:
				alignments_to_keep[seq_name]=str(seq)
	
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
								if int(vv) == residue:
									if int(vv) not in featurecount:
										featurecount.append(vv)
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
# ---------------------------------------------------------------------------------------------------------------------------------------------
def create_svg(sequences_dict, positions, colordict, startposition, windowsize, poi, forbidden, proteinfeatures):
    heatmapper = {}
    startposition_checker = startposition
    Heatmapstart = 60-(len(coloringcategories)+1)*10
    Konservierungsypsilon = Heatmapstart - 20
    Categoryypsilon = Heatmapstart - 130	
    #### do this when havng constructed the dictionary with interesting positions
    #### here it is supplied as is, but needs to be further modified
    for item in positions:
        for categ in colordict:
            if categ not in positions[item]:
                positions[item][categ]=[]
    ### adjust the proteinfeatures dictionary
    for feat in proteinfeatures:
        oldstuff = proteinfeatures[feat]
        maxi = max(proteinfeatures[feat])
        mini = min(proteinfeatures[feat])
        for i in range(mini,maxi+1):
            proteinfeatures[feat].append(i)
    if startposition == "none":
        startposition = 1
    try:
        filename = translator[poi]+"_Position"+str(startposition)+"_Windowsize"+str(windowsize)+".svg"
    except:
        filename = poi+"_Position"+str(startposition)+"_Windowsize"+str(windowsize)+".svg"
    dwg = svgwrite.Drawing(filename, profile='full')
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
    for uniprot in roworder:
        seq 	= sequences_dict[uniprot]
        namus 	= uniprot
        startingpoint = startposition - windowsize	### this is required for the correct labeling according to the sequence of interest
        try:
            drawname = translator[namus]
        except:
            drawname = namus

	#print drawname
        if poi in namus:	##### this if/else conditional can probably be put in yet another function to reduce the amount of code being used here
            old_x = x
            old_y = y
            x = 50
            y = 60
            dwg.add(dwg.rect((x-90, y), (80, 14), fill="yellow"))
            if len(drawname) < 8:
                dwg.add(dwg.text(drawname, insert = (x-60,y+7), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill='black'))
                dwg.add(dwg.text(drawname, insert = (x-60,Konservierungsypsilon+5), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill='black'))
            else:
                dwg.add(dwg.text(drawname, insert = (x-60,y+7), text_anchor='middle', dominant_baseline='central', font_size='7px', font_family='Arial', font_weight='bold', fill='black'))
                dwg.add(dwg.text(drawname, insert = (x-60,Konservierungsypsilon+5), text_anchor='middle', dominant_baseline='central', font_size='7px', font_family='Arial', font_weight='bold', fill='black'))
        else:
            if len(drawname) < 8:
                dwg.add(dwg.text(drawname, insert = (x-60,y+7), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill='black'))
            else:
                dwg.add(dwg.text(drawname, insert = (x-60,y+7), text_anchor='middle', dominant_baseline='central', font_size='7px', font_family='Arial', font_weight='bold', fill='black'))
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
                                dwg.add(dwg.text(startcounter, insert=(35, y+8), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill="black"))					
                                firstdone = "true"
                            if totalcount not in forbidden: ### totalcount is int and forbidden is a list of ints
                                if poi in namus:
                                    konserv_val = Konserve[startnumberlabel][0]
                                    if float(konserv_val)>= 0.7:
                                        dwg.add(dwg.rect((x,y),(10,len(roworder)*20), fill= clustaltypes[Clustalcolors[letter.upper()]], opacity=0.2))
                                    viewboxcounter += 1

                                    if int(startnumberlabel) == int(startposition):
                                            position_interest_x = x
                                            position_interest_y = y
                                            dwg.add(dwg.rect((x, Konservierungsypsilon), (10, 14), fill="black"))
                                            toproof = 1-float(konserv_val)
                                            dwg.add(dwg.rect((x, Konservierungsypsilon), (10, 14*toproof), fill="white"))
                                            if startposition_checker != "none":
                                                dwg.add(dwg.text(str(startnumberlabel), insert=(x+5, Konservierungsypsilon-3), text_anchor='middle', dominant_baseline='central', font_size='8px', font_family='Arial', font_weight='bold', fill='red'))			
                                    else:
                                            if int(startnumberlabel)>= int(startingpoint):
                                                dwg.add(dwg.rect((x, Konservierungsypsilon), (10, 14), fill="black"))
                                                toproof = 1-float(konserv_val)
                                                dwg.add(dwg.rect((x, Konservierungsypsilon), (10, 14*toproof), fill="white"))
                                                if int(startnumberlabel)%10 == False:
                                                    dwg.add(dwg.text(str(startnumberlabel), insert=(x+5, Konservierungsypsilon-3), text_anchor='middle', dominant_baseline='central', font_size='8px', font_family='Arial', font_weight='bold', fill='black'))
                                    for feat in proteinfeatures:
                                            if startnumberlabel in proteinfeatures[feat]:
                                                if feat not in tempfeat:
                                                    tempfeat[feat]=[featurecolors[featcount],featcount]
                                                    featcount+=1
                                                elevator = tempfeat[feat][1]
                                                elevator_floor = 0
                                                if elevator >= 10:
                                                    if elevator_floor <= 10:
                                                        elevator = elevator_floor
                                                        elevator_floor += 1
                                                    else:
                                                        elevator_floor = 0
                                                        elevator = elevator_floor
                                                y_level = -35 + (elevator*3)
                                                y_level_text = -95 + (elevator*4.5)
                                                dwg.add(dwg.rect((x, y_level), (10, 2), fill=tempfeat[feat][0]))
                                                if "done" not in tempfeat[feat]:
                                                    dwg.add(dwg.text(str(feat), insert=(x+15, y_level_text), text_anchor='middle', dominant_baseline='central', font_size='8px', font_family='Arial', font_weight='bold', fill=tempfeat[feat][0]))
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
                                            dwg.add(dwg.circle((x+5, y+7.5), (radius), fill=colordict[colorcateg]))
                                            dwg.add(dwg.text(letter, insert=(x+5, y+8), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill="black"))
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
                                        dwg.add(dwg.text(letter, insert=(x+5, y+8), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill="black"))
                                except:
                                    dwg.add(dwg.text(letter, insert=(x+5, y+8), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill="black"))	
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
        dwg.add(dwg.text(endcounter, insert=(lastx+20, y+8), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill="black"))	
			
        if poi in namus:
            dwg.add(dwg.rect((-40,y),(x+50,14), fill="none",stroke="black",stroke_width=1))	
            x = 50
            y = old_y
        else:
            x = 50
            y += 20
    viewboxwidth = (viewboxcounter+140)
    viewboxheight = len(roworder)*20+100+(200-Categoryypsilon)
    dwg.viewbox(-50, Categoryypsilon-80,viewboxwidth,viewboxheight)

    if startposition_checker != "none":
    	dwg.add(dwg.rect((position_interest_x, position_interest_y), (10, len(roworder)*20),fill="none",stroke="black",stroke_width=1))

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
    for category in colors:
        try:
            heatmap_maximum = max(maxfinder[category])
        except:
            heatmap_maximum = 1
        dwg.add(dwg.text(category, insert=(20, mapy+5), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill='black'))
        for xval in heatmapper:
		#print xval, "\t", mapy
            try:
                opac = float(heatmapper[xval][category])/float(heatmap_maximum)
            except:
                opac = 0.0
            if float(opac) == 0.0:
                dwg.add(dwg.rect((xval, mapy), (10, 10), fill="lightblue", opacity = 0.15 ))
            else:
                dwg.add(dwg.rect((xval, mapy), (10, 10), fill=colors[category], opacity = opac ))
            if mapy == 20:
                pass		
        dwg.add(dwg.rect((50, mapy), (lastx-mapx-10, 10),fill="none",stroke="black",stroke_width=0.5))	### <<<<
        mapy += 10

    for i in range(50,lastx-10,10):
        dwg.add(dwg.rect((i, Heatmapstart), (10, 40),fill="none",stroke="black",stroke_width=0.5))

    x = 50
    y = 0
    for category in colors:
        dwg.add(dwg.rect((x-30, Categoryypsilon), (60, 10), fill=colors[category]))
        dwg.add(dwg.text(category, insert=(x, Categoryypsilon+5), text_anchor='middle', dominant_baseline='central', font_size='10px', font_family='Arial', font_weight='bold', fill='black'))
        x += 60
    dwg.save()
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
# ---------------------------------------------------------------------------------------------------------------------------------------------
alignmentfile = sys.argv[4]	#### change this to the location of the alignmentfile.
###
sequences 	= CUSTOM_ALIGN(alignmentfile)
protein_of_interest = sys.argv[1]
try:
    position_of_interest = int(sys.argv[2])
except:
    position_of_interest = str(sys.argv[2])
window = int(sys.argv[3])

if position_of_interest == "none":
    window = 30000

with open(sys.argv[5]) as f:
    data_align = f.read()
positions = ast.literal_eval(data_align)

Konserve, TheForbiddenPositions 	= conservation_checker(protein_of_interest,sequences, positions)

### note: This is the real CONNECTORalpha dictionary I could fetch from uniprot on 31.01.2023


try:
    with open(sys.argv[6]) as ff:
        data_feat = ff.read()
    feature_dict = ast.literal_eval(data_feat)
    for k in feature_dict:
        start 	= feature_dict[k][0]
        end	= feature_dict[k][1]
        for i in range(start, end+1):
            if i not in feature_dict[k]:
                feature_dict[k].append(i)
except:
    logging.exception("message")
    pass
### to define the annotation colors we want to use

positioncolors = ["lightgreen","salmon","yellow","orchid","lightblue"]
colors = {}
coloringcategories = []
counter = 0
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
featurecolors = ["firebrick","tomato","orange","olive","palegreen","teal","dodgerblue","blueviolet","deeppink"]
###
create_svg(sequences, positions, colors, position_of_interest, window, protein_of_interest, TheForbiddenPositions, feature_dict)	### def create_svg(sequences_dict, positions, colors, startposition, windowsize, poi):


