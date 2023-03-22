___
# Introduction
## This script
The file **create_svg_20230314_kinases_V2.py** comes with several python2 functions that can annotate or highlight positions in a Clustal-formatted alignment.

___
## General use case
In general, this script needs two things: 
* an alignment in clustal format 
* a python dictionary formatted as {Protein:{Feature:\[Residues]}}, **also see Version3/AlignmentAnnotationdDictionary_Example.txt**
* an **optional** protein feature dictionary, formatted as (and based on protein of interest) {Feature:[Startposition, Endposition]}, **also see Version3/RegionAnnotationsDictionary_Example.txt**

to work. The dictionary can be created elsewhere and could contain different features than the one I included here, so it is **versatile**.

___
## Required libraries/software

Python 2.7. (sorry folks)

import svgwrite

from Bio import AlignIO

import ast

imposrt sys
___
## Features
- **New** Added transparent rectangles to highlight a sequence conservation (= identity) over >= 70 %, based on the sequence of interest. The colors for this are taken from CLUSTAL/Jalview.
- **New** Change highlighting to circles. Circle radius can later be adjusted based on evidence.
- **New** Added basic heatmapping above the alignment, showing how many highlights per position & per category we have.
- **New** Added start and end positions for each displayed sequence.
- Command line functionality. 
To use the script we can now execute the following command:
`python create_svg_20230314_kinasesV2.py P61586 34 30 RHOA_Alignment.clustal_num AlignInfos_RHOA.txt Features_RHOA.txt` 

This command has several fields after calling the script:

| Field        | Example           | Description  |
| ------------- |:-------------:| -----:|
| 0     | P61586 | The uniprot ID of the protein we are interested in |
| 1     | 34 | The position to be highlighted |
| 2     | 30 | The Windowsize, we show +/- the windowsize of residues around the highlighted position|
| 3     | RHOA_Alignment.clustal_num | The alignment file |
| 4     | AlignInfos_RHOA.txt | The file containing positional information |
| 5     | Features_RHOA.txt | A file containing structural/domain features, numbering based on **protein of interest** |

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
## What is next?
- Make circle size adjustable by evidence.

___
## The most recent type of results
The result of executing `python create_svg_20230314_kinasesV2.py P61586 34 30 RHOA_Alignment.clustal_num AlignInfos_RHOA.txt Features_RHOA.txt`.
<img src="https://github.com/russelllab/kinaseResistance/blob/ac8fb82c5fbf26a116d23f3b84c61e7c543108b2/Create_SVG/Version_K(inases)/MAP2K3_Position84_Windowsize30.svg?sanitize=true">
