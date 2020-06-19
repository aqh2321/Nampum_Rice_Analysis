'''
SECTION_1
building library class
can do individual object check
'''
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

class library:
    def __init__(self,record):
        self.record = record
        self.origin = record[0].seq


    def LF(self):
        lf = {}
        for i in range(1,8):
                spacer = self.record[i].seq
                spacerrev = spacer.reverse_complement()
                lenth = len(spacer)
                for m in range(0, len(self.origin)-lenth+1):
                    comparison = self.origin[m:lenth]
                    if spacer == comparison:
                        lf[self.record[i].id, "LF"] = self.origin[(m-30):(m-20)]
                    elif spacerrev == comparison:
                        lf[self.record[i].id, "LF"] = self.origin[(m-30):(m-20)]
                    else:
                        m += 1
                        lenth += 1
        del lf['Protospacer_6', 'LF']
        lf['Protospacer_5', 'MF'] = 'CTCGCTAACC'
        return lf

    def RF(self):
        rf = {}
        for i in range(1,8):
                spacer = self.record[i].seq
                spacerrev = spacer.reverse_complement()
                lenth = len(spacer)
                for m in range(0, len(self.origin)-lenth+1):
                    comparison = self.origin[m:lenth]
                    if spacer == comparison:
                        rf[self.record[i].id, "RF"] = self.origin[(m-30):(m-20)]
                    elif spacerrev == comparison:
                        rf[self.record[i].id, "RF"] = self.origin[(m-30):(m-20)]
                    else:
                        m += 1
                        lenth += 1
        del rf['Protospacer_5', 'RF']
        return rf

    def LR(self):
        lr = {}
        for i in range(1,8):
                spacer = self.record[i].seq
                spacerrev = spacer.reverse_complement()
                lenth = len(spacer)
                for m in range(0, len(self.origin)-lenth+1):
                    comparison = self.origin[m:lenth]
                    if spacer == comparison:
                        lr[self.record[i].id, "LR"] = self.origin[(m-30):(m-20)].reverse_complement()
                    elif spacerrev == comparison:
                        lr[self.record[i].id, "LR"] = self.origin[(m-30):(m-20)].reverse_complement()
                    else:
                        m += 1
                        lenth += 1
        del lr['Protospacer_6', 'LR']
        lr['Protospacer_5', 'MR'] = 'GGTTAGCGAG'
        return lr

    def RR(self):
        rr = {}
        for i in range(1,8):
                spacer = self.record[i].seq
                spacerrev = spacer.reverse_complement()
                lenth = len(spacer)
                for m in range(0, len(self.origin)-lenth+1):
                    comparison = self.origin[m:lenth]
                    if spacer == comparison:
                        rr[self.record[i].id, "RR"] = self.origin[(m-30):(m-20)].reverse_complement()
                    elif spacerrev == comparison:
                        rr[self.record[i].id, "RR"] = self.origin[(m-30):(m-20)].reverse_complement()
                    else:
                        m += 1
                        lenth += 1
        del rr['Protospacer_5', 'RR']
        return rr

def get_pair(line):
    k, sep, v = line.partition("!")
    return (k, v)

def sumlist(target):
    if len(target) < 2:
        return target[0]
    else:
        return (target[0]+sumlist(target[1:]))

def flipFR(sequence):
    x = list(sumlist(sequence))
    if x[2] == 'F' and x[5] == 'R':
        x[2] = 'R'
        x[5] = 'F'
        newsequence = sumlist(x)
        return newsequence
    elif x[2] == 'R' and x[5] == 'F':
        x[2] = 'F'
        x[5] = 'R'
        newsequence = sumlist(x)
        return newsequence

def simplify(sequence):
    x = list(sumlist(sequence))
    if 'F' not in x:
        lenth = len(x)
        for m in range(lenth//3):
            del x[(lenth//3 - m)*3 - 1]
            newsequence = sumlist(x)
        return newsequence
    while 'F' in x:
        x.remove('F')
        newsequence = sumlist(x)
    return newsequence

def flipnum(sequence):
    sequence.reverse()
    return sequence

def flip(test):
    if test[0][2] == 'F' and test[1][2] == 'F':
        return simplify(test) #default
    elif test[0][2] == 'R' and test[1][2] == 'R':
        return simplify(flipnum(test))
    else:
        if (int(test[0][0]) - int(test[1][0])) < 0:
            return sumlist(test) #default, small number in the front
        elif (int(test[0][0]) - int(test[1][0])) == 0:
            if test[0][1] != test[1][1]:
                if test[0][2] == 'F': #default specific around same protospacer, first position ends in F
                    return sumlist(test)
                else:
                    return flipFR(flipnum(test))
            else:
                return sumlist(test) #head to head or tail to tail duplication
        else:
            return flipFR(flipnum(test))

def flipmultiple(test):
    letter = []
    x = list(sumlist(test))
    lenth = len(x)
    for m in range(lenth//3):
        letter.append(x[(m+1)*3-1])
    check = 'F'
    checkrev = 'R'
    bEqual = True
    brevEqual = True
    for item in letter:
        if item != check:
            bEqual = False
        elif item != checkrev:
            brevEqual = False

    if bEqual:
        return simplify(test)
    elif brevEqual:
        return simplify(flipnum(test))
    else:
        return sumlist(test)

def listsort(list, sort):
    newlist = []
    for i in sort:
        for item in list:
            if int(item[0]) == i:
                newlist.append(item)
    return newlist
