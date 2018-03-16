Mod-seeker
==========

Mod-seeker is an analysis pipeline for Mod-seq data.  Mod-seeker requires fastq sequence files, a target sequence built as a bowtie index, and gene annotations as input files. Mod-seeker returns sites of significant modification on each annotated gene as output, along with odds-ratios and results of significance testing at each site.

Mod-seeker is implemented in Python. It is comprised of two programs. Mod-seeker-map.py is a wrapper program that calls on other open source software to analyze sequence data throught the following steps:

1) Trim 5’ and 3’ adapter sequences (using cutadapt).
2) Align modification-stop sequences (no 5’ adapter) to a target sequence (using bowtie)
3) Remove untemplated nucleotides added by reverse transcriptase from resulting alignments
4) Map aligned sequence reads to locations on annotated genes (using samtools and bedtools).
5) Count the number of modification stop reads at each position in each gene from each sample.

Mod-seeker-stats.py can be run on the output files from Mod-seeker-map.py (CountMod*tab) and returns sites of significant modification in the treated sample, as compared to the control. If no replicates are present, the program usese the Chi-square test. Any number of replicates can be added, with analysis utilizing the Cochran-Mantel-Haenzsel test.

Dependencies
============

Mod-seeker requires Python 2.6 or later. In addition, you must have the following installed and in your PATH:
1) cutadapt (https://code.google.com/p/cutadapt/).
2) bowtie (http://bowtie-bio.sourceforge.net/index.shtml).
3) samtools (http://samtools.sourceforge.net).
4) bedtools (https://github.com/arq5x/bedtools2).

Usage
=====

1) Place your fastq files in the “InputFiles” directory.
2) Check the settings in map_settings.txt and then type “python Mod-seeker-map.py” in the command line.
3) Copy or move the CountMod*.tab files into the Mod-seeker directory.
4) Check the settings in stats_settings.txt and then type “python Mod-seeker-stats.py” in the command line.

Settings files are used to set the parameters for Mod-seeker-map.py (map_settings.txt) and Mod-seeker-stats.py (stats_settings.txt). Examples are given below.

map_settings.txt:

ref_path = /Applications/bowtie-1.0.0/indexes/  # Set the path to the bowtie target genome (from bowtie-build)
ref_genome = Scer # The name of the bowtie target genome (from bowtie-build)
GeneAnnotationFile = ~/Desktop/Mod-seq/SGD_SacCer3_filtered.gff # Set the path to the annotation file (in gff or bed format only).
adapter5 = ^ATCGTAGGCACCTGAAA # The 5’ adapter sequence used to make Mod-seq libraries
adapter3 = CTGTAGGCACCATCAAT # The 3’ adapter sequence used to make Mod-seq libraries.
output_path = default # Mod-seeker-map.py will create a “Results” directory
bowtie_alignment_threads = 2 # number of threads (-p #) for bowtie alignment
max_python_threads = 4 # number of threads for Mod-seeker-map.py.

stats_settings.txt:

FDR = 0.05 # The FDR correction level (here 5%)
Odds_Ratio_Threshold = 1.5 # The minimal fold-enrichment in treatment vs. control comparisons
Number_Of_Replicates = 2

File_List:
Treated_1 = CountMod_Scer_Scer_WT_Rep1_100mM_small.tab
Control_1 = CountMod_Scer_Scer_WT_Rep1_0mM_small.tab
Treated_2 = CountMod_Scer_Scer_WT_Rep2_100mM_small.tab
Control_2 = CountMod_Scer_Scer_WT_Rep2_0mM_small.tab

Test Data
=========

The pipeline is provided with short test data files in the InputFiles directory.  These files contain 10,000 sequence reads from each of two treatment / control replicate samples from S. cerevisiae.  To use these files, first ensure you have all necessary software installed (see dependencies), have built a bowtie target index for S. cerevisiae, and have a gene annotation file specified in map-settings.txt (in bed or gff format).

Reference
=========
Talkish J, May GE, Lin Y, Woolford JL. Jr, McManus CJ. Mod-seq: High-throughput sequencing for chemical probing of RNA structure. (accepted to RNA, 2014).

License
=======
Mod-seeker is provided under the Gnu General Public License (v3).

Copyright (C) 2014  Yizhu Lin
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Changes
=======

v0.1
———-
Initial release
