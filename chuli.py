import os
import sys
filename=sys.argv[1]
outfile=open(sys.argv[2],'w')
for line in open(filename):
    newline=line.rstrip().split('\t')
    if len(newline)==8:
        #newline.append(' ')
        #outfile.write('\t'.join(newline)+'\n')
        continue
    else:
        outfile.write('\t'.join(newline)+'\n')
outfile.close()
