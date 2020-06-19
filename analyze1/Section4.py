'''Section4 writing output'''

import xlsxwriter
from analyze1.Supporter import get_pair, listsort
from collections import defaultdict
import os, sys

sortprep = []
types = []
finaldictlist = []
name = []
files = open(sys.argv[1])#filenames

for file in files:
    with open('%s_output.txt'%file[:-1]) as counting:
        name.append(file[:-1])
        analysis = {}
        analysis = defaultdict(int)
        for line in counting:
            k, v = get_pair(line)
            if int(k[0]) not in sortprep:
                sortprep.append(int(k[0]))
            if k not in types:
                types.append(k)
            analysis[k] += 1
        sortprep.sort()
        finaldictlist.append(analysis)

files.close()

finallist = listsort(types, sortprep)

workbook = xlsxwriter.Workbook('report.xlsx')
worksheet = workbook.add_worksheet()

col = 0
row = 0

worksheet.write(0,0, 'Type')
for item in finallist:
    row += 1
    worksheet.write(row, col, item)

col = 0
row = 0

for i in range(len(name)):
    worksheet.write(0,i+1, name[i])
    dictlist = finaldictlist[i]
    for item in finallist:
        row += 1
        for k,v in dictlist.items():
            if k == item:
                worksheet.write(row, col + 1, v)
        #if not find match, create empty cells, add 0 in excel afterwards
    row = 0
    col += 1


workbook.close()

print("======finish writing all the documents======")
