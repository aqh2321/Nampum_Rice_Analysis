# Nampum_Rice_Analysis
Pipelines to analyze rearrangement and indel mutations around CRISPR cut sites

code runs with python 3.7.7

**important package to be installed**   
[biopython 1.73](https://biopython.org/wiki/Download), [bwa](https://github.com/lh3/bwa) and their dependencies

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
To add labels to rearrangement and indel data, please enter ```y``` to continue, or ```n``` to exit.
A csv file containing different general types of mutations on different protospacers will also be generated if proceed to ```y```.

**Settings:**  
Default average length per seq is 150bp, default length setting for fragment is 25bp, and default data folder is 'raw_data'.  
A demo for analyzing indels for one file is available in ```handle_sam.ipynb```
A demo for combining rearrangement and indel data together, add labels, and calculate the frequencies is avai in ```add_del_inv_label.ipynb```

**Additional Information**
Some scripts in folder analyze2 allows manual adjustments for more functions. Information are noted at the beginning of each script.
For example, in autoscan.py and antoscan2.py, there're instructions to generate tdf files from sam files if needed.
Also, in analysisindel.py and analysisindel2.py, there're options to turn on a threshold filter to filter out low mutation reads.
