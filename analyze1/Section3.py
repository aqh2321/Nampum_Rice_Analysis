'''Section4'''

# In[ ]:
import os,sys
from analyze1.Supporter import get_pair,sumlist,flipFR,simplify,flipnum,flip,flipmultiple
file = sys.argv[1]

#SECTION_4
#a) wt
with open('%s.txt'%file) as f:
    with open('%s_output.txt'%file, "w") as output:
        names = []
        ids = []
        position = []
        location = []

        for line in f:
            k, v = get_pair(line)
            ids.append(k)
            position.append(v[1:4])
            location.append(v[4:])

        for m in ids:
            if m not in names:
                names.append(m)

        for name in names:
            actual =[]
            #distance = []
            for n in range(len(ids)):
                if ids[n] == name:
                    actual.append(position[n])
                    #distance.append(location[n])
            if len(actual) == 2:
                writing = flip(actual)
                #writing_distance = abs(int(distance[1])-int(distance[0]))
                output.write("%s ! %s\n" % (str(writing),name))
            elif len(actual) > 2:
                writing = flipmultiple(actual)
                output.write("%s ! %s\n" % (str(writing), name))
        print('\n$_$ finish section3: rearrangement searching for %s' %file)
