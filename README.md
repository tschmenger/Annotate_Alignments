___
# Introduction
## This script
The file **Latest/create_svg_20230327.py** comes with several python2 functions that can annotate or highlight positions in a Clustal-formatted alignment.

___
## General use case
In general, this script needs two things: 
* an alignment in clustal format 
* a python dictionary formatted as {Protein:{Feature:\[Residues]}}, 
* an **optional** protein feature dictionary, formatted as (and based on protein of interest) {Feature:[Startposition, Endposition]}.

to work. The dictionary can be created elsewhere and could contain different features than the one I included here, so it is **versatile**.

**For example** using the files in /Latest:

`python create_svg_20230327.py P61586 34 30 RHOA_BlastpExample_ClustalMSA.clustal RHOA_Blastp_info.txt`

___
## Preparing alignment
This applies to users who do not already have an alignment. The following steps can used to create an alignment using blastp and clustal omega.

- Step 1: Download the fasta sequence of your protein of interest. For RHOA you could do this via Uniprot, like [this](https://rest.uniprot.org/uniprotkb/P61586.fasta). Note: You can easily build this url using **https://rest.uniprot.org/uniprotkb/** + uniprotID +**.fasta**


- Step 2: Use the downloaded fasta sequence to perform a blast search for similar sequences [blastp](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins). For more information on how to use BLAST please see [Blast Help](https://blast.ncbi.nlm.nih.gov/doc/blast-help/). **Make sure to select "Swissprot" as a database, or otherwise make sure that the retained accessions will be UniprotIDs.**

- Step 3: Select the sequences you prefer and download them.

![BlastP Example](https://github.com/tschmenger/Annotate_Alignments/blob/cd9645b4be6740b617e480c256a3d558ac4fb6a1/manual_blastp.png?sanitize=true)|
|:--:| 
| *Make sure to download the complete sequences.* |

- Step 4: Add the protein of interests fasta manually to the top of the just downloaded file. Change the formatting to roughly mimick the formatting of the remaining entries.

![fasta example](https://github.com/tschmenger/Annotate_Alignments/blob/f462d36028dde94573a3918502b88d342cef4f3e/manual_FastaAdded.png?sanitize=true)|
|:--:| 
| *This is how your prepared fasta sequences should look like.* |

- Step 5: Perform a multiple sequence alignment using [Clustal Omega](https://www.ebi.ac.uk/Tools/msa/clustalo/). Make sure to select **Protein**. Download & save the alignment file for usage with this script. Input the sequences via copy & paste or upload a file.

![Clustal Example](https://github.com/tschmenger/Annotate_Alignments/blob/202367f70da00a851817482dfaf57b9c5146c3c7/manual_clustal.png?sanitize=true)

Make sure you download the complete MSA (including the clustal version, followed by 2 empty lines, followed by the MSA).

- Step 6: Annotate the alignment manually or programmatically with whatever information you want, following the aforementioned format. **Recommended** to use the **Create_Information.py** script (see next section).

___
## Required libraries/software (side script|Create_Information.py)

Python 2.7.

import sys

import urllib

import re

___
## Prepare information
Simply run for example **general** `python Create_Information.py Sequences.txt >> Info_Output.txt`

or using **example** `python Create_Information.py RHOA_BlastP_sequences.txt >> RHOA_Blastp_info.txt`.

This will produce a python dictionary file that is suitable for usage with the main script, collecting UniProt information of the **BINDING**, **ACT_SITE**, **MUTAGEN** or **VARIANT** categories.

**You should now have both the alignment file and an annotation file and can proceed to use the main script below.**

___
## Optional: Retrieving features
We can also retrieve structural features (domains, motifs) from [Interpro](https://www.ebi.ac.uk/interpro/) following this link format:

`https://www.ebi.ac.uk/interpro/api/protein/UniProt/UniprotID/?residues&page_size=200`

that will be, in our example [https://www.ebi.ac.uk/interpro/api/protein/UniProt/P61586/?residues&page_size=200](https://www.ebi.ac.uk/interpro/api/protein/UniProt/P61586/?residues&page_size=200).

There, simply copy & paste the python dictionary into a new file that you can use together with the main script (see below for how to do that).

![Interpro Dictionary Start](https://github.com/tschmenger/Annotate_Alignments/blob/c286bb1273af2be4edef650cb0db541de08d4866/manual_interpro.png)|
|:--:| 
| *The beginning of the InterPro python dictionary.* |

___
## Required libraries/software (main script|Create_svg.py)

Python 2.7.

import svgwrite

from Bio import AlignIO

import ast

import sys
___
## Features
- **New** Added transparent rectangles to highlight a sequence conservation (= identity) over >= 70 %, based on the sequence of interest. The colors for this are taken from CLUSTAL/Jalview.
- **New** Change highlighting to circles. Circle radius can later be adjusted based on evidence.
- **New** Added basic heatmapping above the alignment, showing how many highlights per position & per category we have.
- **New** Added start and end positions for each displayed sequence.
- Command line functionality. 

To use the script we can now execute the following command:
`python create_svg_20230327.py P61586 34 30 RHOA_BlastpExample_ClustalMSA.clustal RHOA_Blastp_info.txt Features_RHOA.txt` 

This command has several fields after calling the script:

| Field        | Example           | Description  |
| ------------- |:-------------:| -----:|
| 0     | P61586 | The uniprot ID of the protein we are interested in |
| 1     | 34 | The position to be highlighted |
| 2     | 30 | The Windowsize, we show +/- the windowsize of residues around the highlighted position|
| 3     | RHOA_BlastpExample_ClustalMSA | The alignment file |
| 4     | RHOA_Blastp_info.txt | The file containing positional information |
| 5     | Features_RHOA.txt | A file containing structural/domain features, numbering based on **protein of interest** This is **optional** and can easily be retrieved from InterPro. |

**Note**: The script allows for a little hack here. If you want a (large) .svg containing the whole alignment just give a big number in field 2, for example 20000. The script will then produce a complete alignment view. **New** Giving "none" instead of a position to be highlighted (field 1) works the same + it removed the position specific rectangle.

- Named Output files. The resultfile will already be named depending on the input settings, so one can easily try different settings. The name follows this format: 
`poi+"_Position"+str(startposition)+"_Windowsize"+str(windowsize)+".svg"`

- Conservation: Gives a black rectangle as an indicator of sequence identity (top) for the POI residue at that position.

- Residue Numbering: Gives the residue number (every 10 residues), based on sequence of interest, and highlights the input residue in red.

- Feature Annotation: Gives a colored background based on type of annotation (taken from uniprot) to the respective residue.

- Sorted Sequences: Sequences with fewer uniprot annotations are sorted to the bottom of the alignment.

- GAPs removed: Gaps are printed with white color (i.e. invisible on a white background). Additionally, columns with more than 90 % GAPs are removed from the alignment. Sequences affected by this (i.e. the up to 10 % of sequences that did not have a gap at that position) **are kept and not removed**. 

- Highlighting protein features, *here* for example p-loop, Switch I and the Effector region of RHOA. We currently support the displaying of up to 9 features (dependent on the given colors in *featurecolors* on line 2518 of this example script).

___
## The most recent type of results
The result of executing `python create_svg_20230314_kinasesV2.py P61586 34 30 RHOA_Alignment.clustal_num AlignInfos_RHOA.txt Features_RHOA.txt`.

<img src="https://github.com/tschmenger/Annotate_Alignments/blob/41455e4ea5adc7696164e70ba5a98bec6bb21215/Version4/P61586_Position34_Windowsize30.svg?sanitize=true">

The result of executing `python create_svg_20230327.py P61586 34 30 RHOA_BlastpExample_ClustalMSA.clustal RHOA_Blastp_info.txt`.

<img src="https://github.com/tschmenger/Annotate_Alignments/blob/cd9645b4be6740b617e480c256a3d558ac4fb6a1/Latest/P61586_Position34_Windowsize30.svg?sanitize=true">
