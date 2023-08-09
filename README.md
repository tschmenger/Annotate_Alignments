___
# Introduction
The file *Latest/Standalone_Annotate_Alignment_V8.py* comes with several python3.6.8 functions that can annotate and highlight positions in a Clustal-formatted alignment. The resultfile is a .svg that is fully compatible for further editing in Inkscape/Illustrator or for inspection via your browser.

___
## Citation
Please cite 

`Gurdeep Singh (1)*+, Torsten Schmenger (1)*+, Juan-Carlos Gonzalez-Sanchez (1), Anastasiia Kutkina (1), Nina Bremec (1), Gaurav Diwan (1), Cristina Lopez (2), Rocio Sotillo (3), Robert B Russell (1)+; "Identify activating, deactivating and resistance variants in protein kinases"; 2023` 

when using this work.
___
## General use case
In general, this script needs only one thing to work: 
* an alignment in clustal format **mandatory**
* **optional** a python dictionary formatted as {Protein:{Feature:\[Residues]}}, 
* **optional** a protein feature dictionary, formatted as (and based on protein of interest) {Feature:[Startposition, Endposition]}.

**Note**: Both dictionaries can either be provided by the user or will be automatically generated by the script. In that case, the script returns both dictionaries for future use.


**For example** using the files in /Latest:

**User prodived dictionaries**

`python3 Standalone_Annotate_Alignment_V8.py P61586 34 30 RHOA_BlastpExample_ClustalMSA.clustal RHOA_Blastp_info.txt Features_RHOA.txt`

**Minimal example**

`python3 Standalone_Annotate_Alignment_V8.py P61586 34 30 RHOA_BlastpExample_ClustalMSA.clustal none none`

___
## Preparing alignment
This applies to users who do not already have an alignment. The following steps can used to create an alignment using blastp and clustal omega.

- Step 1: Download the fasta sequence of your protein of interest. For RHOA you could do this via Uniprot, like [this](https://rest.uniprot.org/uniprotkb/P61586.fasta). Note: You can easily build this url using **https://rest.uniprot.org/uniprotkb/** + uniprotID +**.fasta**


- Step 2: Use the downloaded fasta sequence to perform a blast search for similar sequences [blastp](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins). For more information on how to use BLAST please see [Blast Help](https://blast.ncbi.nlm.nih.gov/doc/blast-help/). **Make sure to select "Swissprot" as a database, or otherwise make sure that the retained accessions will be UniprotIDs.**

- Step 3: Select the sequences you prefer and download them.

![BlastP Example](https://github.com/tschmenger/Annotate_Alignments/blob/cd9645b4be6740b617e480c256a3d558ac4fb6a1/manual_blastp.png?sanitize=true)|
|:--:| 
| *Make sure to download the complete sequences.* |

- Step 4: Add the protein of interests fasta manually to the top of the just downloaded file. Change the formatting to roughly mimic the formatting of the remaining entries.

![fasta example](https://github.com/tschmenger/Annotate_Alignments/blob/f462d36028dde94573a3918502b88d342cef4f3e/manual_FastaAdded.png?sanitize=true)|
|:--:| 
| *This is how your prepared fasta sequences should look like.* |

- Step 5: Perform a multiple sequence alignment using [Clustal Omega](https://www.ebi.ac.uk/Tools/msa/clustalo/). Make sure to select **Protein**. Download & save the alignment file for usage with this script. Input the sequences via copy & paste or upload a file.

![Clustal Example](https://github.com/tschmenger/Annotate_Alignments/blob/202367f70da00a851817482dfaf57b9c5146c3c7/manual_clustal.png?sanitize=true)

Make sure you download the complete MSA (including the clustal version, followed by 2 empty lines, followed by the MSA).

- Step 6: Annotate the alignment manually or programmatically with whatever information you want, following the aforementioned format. **Recommended** to use the **Create_Information.py** script (see next section).


___
## Required libraries/software (main script)

Python 3.6.8+
___
## Features
- **NEW** Returns dictionaries of collected information for re-usage by the user in future runs.
- **NEW** Fused preparation and main script into a single script.
- **NEW** Upgraded script to python 3 plus added some cosmetic changes.
- **New** Added tooltips that show up on mouseover events. Works on residues with functional information.
- **New** Added transparent rectangles to highlight a sequence conservation (= identity) over >= 70 %, based on the sequence of interest. The colors for this are taken from CLUSTAL/Jalview.
- **New** Change highlighting to circles. Circle radius can later be adjusted based on evidence.
- **New** Added basic heatmapping above the alignment, showing how many highlights per position & per category we have.
- **New** Added start and end positions for each displayed sequence.
- Command line functionality. 

To use the script we can now execute the following command:
`python Standalone_Annotate_Alignment_V8.py P61586 34 30 RHOA_BlastpExample_ClustalMSA.clustal RHOA_Blastp_info.txt Features_RHOA.txt` 

This command has several fields after calling the script:

| Field        | Example           | Description  |
| ------------- |:-------------:| -----:|
| 0     | P61586 | The uniprot ID of the protein we are interested in |
| 1     | 34 | The position to be highlighted |
| 2     | 30 | The Windowsize, we show +/- the windowsize of residues around the highlighted position|
| 3     | RHOA_BlastpExample_ClustalMSA | The alignment file.  |
| 4     | RHOA_Blastp_info.txt | The file containing positional information. IDs must match to those present in the alignment file. **Can be left empty by writing "none". The script will then automatically populate the dictionary via Uniprot**|
| 5     | Features_RHOA.txt | A file containing structural/domain features, numbering based on **protein of interest.** IDs must match those present in the alignment file. **Can be left empty by writing "none". The script will then automatically populate the dictionary via Interpro** |

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
The result of executing `python3 Standalone_Annotate_Alignment_V7.py P61586 34 30 RHOA_BlastpExample_ClustalMSA.clustal RHOA_Blastp_info.txt Features_RHOA.txt`.

<img src="https://github.com/tschmenger/Annotate_Alignments/blob/e9a1434d2f515a950cbb3734d9430d7020339324/Latest/RHOA_HUMAN_Position34_Windowsize30.svg?sanitize=true">

___
# License
tschmenger/Annotate_Alignments is licensed under the **GNU General Public License v3.0**

# Acknowledgements
Thanks to @gurdeep330 for providing feedback and suggestions.
