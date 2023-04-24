import sys
import urllib
import re
nums = re.compile(r"[+-]?\d+(?:\.\d+)?")
whitespace_killer=re.compile(r"\s+")
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
idcontainer = []
def idgetter(alignmentfile):
	with open(alignmentfile,"r") as idsource:
		for line in idsource:	### >P61585.1:1-193 RecName: Full=Transforming protein RhoA; AltName: Full=Gb; AltName: Full=p21; Flags: Precursor [Bos taurus]
			if ">" in line:
				try:
					identifier = line.split(">")[1].split(".")[0].replace(" ","")
				except:	
					identifier = line.split(">")[1].replace(" ","").replace("\n","")
				if identifier not in idcontainer:
					idcontainer.append(identifier)
	return idcontainer

infodictionary = {}
translatedictionary = {}
def HOMOL_UNIP_GETTER(idc, inputdictionary_two):
	for k in idc:
			content=urllib.urlopen("https://www.uniprot.org/uniprot/" + k + ".txt")
			checktype = ""
			checker = "FALSE"
			position = ""
			comment = ""
			for line in content:
				try:
					if "ID" and "Reviewed;" in line:	### ID   RHOA_HUMAN              Reviewed;         193 AA.
						newline = whitespace_killer.sub(" ",line).replace("\n","")
						langname = newline.split(" ")[1].replace(" ","").replace("\n","")
						if translatedictionary.has_key(k)==False:
							translatedictionary[k]=langname
					if "FT   VARIANT" in line:
						if checker == "FALSE":	### this will happen when I get to such a line:			FT   VARIANT         47
							comment = ""				
							newline = whitespace_killer.sub(" ",line).replace("\n","")
							position = newline.split(" ")[2]
							checker = "TRUE"
							checktype = "VARIANT"
					elif "FT   MUTAGEN" in line:
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
					elif "FT   MOD_RES" in line:	### this will happen when I get to such a line:			FT   BINDING         38
						if checker == "FALSE":	### this will happen when I get to such a line:			FT   MUTAGEN         226
							comment = ""	
							newline = whitespace_killer.sub(" ",line).replace("\n","")
							position = newline.split(" ")[2]
							checker = "TRUE"
							checktype = "MOD_RES"
					elif "FT   ACT_SITE" in line:
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
										if inputdictionary_two.has_key(k)==False:
											inputdictionary_two[k]={}
											inputdictionary_two[k]["MOD_RES"]=[str(homol_mutat)]
										elif inputdictionary_two[k].has_key("MOD_RES")==False:
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
												if inputdictionary_two.has_key(k)==False:
													inputdictionary_two[k]={}
													inputdictionary_two[k]["ACT_SITE"]=[str(bindpos)]
												elif inputdictionary_two[k].has_key("ACT_SITE")==False:
													inputdictionary_two[k]["ACT_SITE"]=[str(bindpos)]
												else:
													inputdictionary_two[k]["ACT_SITE"].append(str(bindpos))
										else:
											if inputdictionary_two.has_key(k)==False:
												inputdictionary_two[k]={}
												inputdictionary_two[k]["ACT_SITE"]=[str(homol_mutat)]
											elif inputdictionary_two[k].has_key("ACT_SITE")==False:
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
												if inputdictionary_two.has_key(k)==False:
													inputdictionary_two[k]={}
													inputdictionary_two[k]["BINDING"]=[str(bindpos)]
												elif inputdictionary_two[k].has_key("BINDING")==False:
													inputdictionary_two[k]["BINDING"]=[str(bindpos)]
												else:
													inputdictionary_two[k]["BINDING"].append(str(bindpos))
										else:
											homol_mutat = position
											if inputdictionary_two.has_key(k)==False:
												inputdictionary_two[k]={}
												inputdictionary_two[k]["BINDING"]=[str(homol_mutat)]
											elif inputdictionary_two[k].has_key("BINDING")==False:
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
												if inputdictionary_two.has_key(k)==False:
													inputdictionary_two[k]={}
													inputdictionary_two[k]["VARIANT"]=[str(bindpos)]
												elif inputdictionary_two[k].has_key("VARIANT")==False:
													inputdictionary_two[k]["VARIANT"]=[str(bindpos)]
												else:
													inputdictionary_two[k]["VARIANT"].append(str(bindpos))
									else:
											if inputdictionary_two.has_key(k)==False:
												inputdictionary_two[k]={}
												inputdictionary_two[k]["VARIANT"]=[str(homol_mutat)]
											elif inputdictionary_two[k].has_key("VARIANT")==False:
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
												if inputdictionary_two.has_key(k)==False:
													inputdictionary_two[k]={}
													inputdictionary_two[k]["MUTAGEN"]=[str(bindpos)]
												elif inputdictionary_two[k].has_key("MUTAGEN")==False:
													inputdictionary_two[k]["MUTAGEN"]=[str(bindpos)]
												else:
													inputdictionary_two[k]["MUTAGEN"].append(str(bindpos))
									else:
											if inputdictionary_two.has_key(k)==False:
												inputdictionary_two[k]={}
												inputdictionary_two[k]["MUTAGEN"]=[str(homol_mutat)]
											elif inputdictionary_two[k].has_key("MUTAGEN")==False:
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
	return inputdictionary_two, translatedictionary
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
alignmentfile 			= sys.argv[1]
ids 				= idgetter(alignmentfile)
infodictionary, translated 	= HOMOL_UNIP_GETTER(ids,infodictionary)

print infodictionary



