import sys

inputfile = open(sys.argv[1])

for line in inputfile.readlines():
 if line.startswith('ID'):
  continue
 if line.find('FUSION') == -1:
  continue
 info = line.split('\t')
 if len(info) > 19:
  print(line, end='')