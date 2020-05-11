# Nampum_Rice_Analysis
Pipelines to analyze rearrangement and indel mutations around CRISPR cut sites

code runs with python 3.7.7

**important package to be installed**   
[biopython 1.73](https://biopython.org/wiki/Download), [bwa](https://github.com/lh3/bwa), and their dependencies

**command**  
To run rearrangement analysis
```
python master.py files.txt
```
The script returns an excel (or csv) with columns of plant name and rows of different mutations detected

To run indel analysis
```
python autoscan.py
```
and then
```
python analysisindel.py
```
**Settings:**  
Default average length per seq is 150bp, default length setting for fragment is 25bp, and default data folder is 'raw_data'.  
A demo for analyzing indels for one file is available in ```handle_sam.ipynb```
