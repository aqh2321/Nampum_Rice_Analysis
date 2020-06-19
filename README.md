# Nampum_Rice_Analysis
Pipelines to analyze rearrangement and indel mutations around CRISPR cut sites

code runs with python 3.7.7

**important package to be installed**   
[biopython 1.73](https://biopython.org/wiki/Download) its dependencies
This project also used [bwa](https://github.com/lh3/bwa). This tool is already set up in the module.

**command**  
To open the help menu:
```
python master.py -h
```
To run the entire project:
```
python master.py raw_data
```
To run rearrangement analysis alone:
```
python master.py raw_data --m=1
```
To run indel analysis alone:
The script returns 2 csv files with the same counts for indel mutations but transposed column and row names.
```
python master.py raw_data --m=2
```
To run indel analysis with high coverage over sequences:
```
python master.py raw_data -c --m=2
```
To run diversity analysis, enter ```y``` to continue, ```n``` to exit.

**Settings:**  
Default average length per seq is 150bp, default length setting for fragment is 25bp, and default data folder is 'raw_data'.  
A demo for analyzing indels for one file is available in ```handle_sam.ipynb```
