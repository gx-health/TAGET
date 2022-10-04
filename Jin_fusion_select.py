import sys

statfile = open(sys.argv[1])
longfile = open(sys.argv[2])
shortfile = open(sys.argv[3])
fafile = open(sys.argv[4])

absdiff = 50
ratiodiff = 0.01

lengthdict = {}
for line in longfile.readlines():
    info = line.rstrip().split('\t')
    if info[0] == '0':
        continue
    thisname = info[3]
    thislength = int(info[7]) - int(info[6])
    if thisname not in lengthdict:
        lengthdict[thisname] = [0, 0, 0]
    lengthdict[thisname][0] += thislength

for line in shortfile.readlines():
    info = line.rstrip().split('\t')
    if info[0] == '0':
        continue
    thisname = info[3]
    thislength = int(info[7]) - int(info[6])
    if thisname not in lengthdict:
        lengthdict[thisname] = [0, 0, 0]
    lengthdict[thisname][1] += thislength

for line in fafile.readlines():
    if line.startswith('>'):
        aa = line.split(' ')
        thisname = aa[0][1:]
        info = aa[1].split(';')
        thislength = int(info[1].split('=')[-1])
        if thisname not in lengthdict:
            lengthdict[thisname] = [0, 0, 0]
        lengthdict[thisname][2] = thislength


for line in statfile.readlines():
    if line.startswith('ID'):
        continue
    info = line.rstrip().split('\t')
    if line.find('FUSION') == -1:
        continue
    longratio = lengthdict[info[0]][0] / lengthdict[info[0]][2]
    shortratio = lengthdict[info[0]][1] / lengthdict[info[0]][2]
    if lengthdict[info[0]][1] - lengthdict[info[0]][0] >= absdiff or shortratio - longratio >= ratiodiff:
        print(info[0] + '\tShort')
    else:
        print(info[0] + '\tLong')
