"""
Author: Zhang Chengsheng, @2019.03.21
Release: 2019.05.27
"""
__VERSION__ = '0.1'


import sys


def overlapSplit(seq, readName, overlap, aveLength, MIN_READ_LENGTH):
    res = ''
    for idx, i in enumerate(range(0, len(seq), aveLength - overlap)):
        left = i
        right = left + aveLength
        left = 0 if left < 0 else left
        seqReal = seq[left:right]
        if len(seqReal) < MIN_READ_LENGTH:
            res = res.rstrip('\n') + seqReal + '\n'
        else:
            res += '{}_{}\n{}\n'.format(readName,idx,seqReal)
        if right > len(seq):
            break
    return res


def main(fileIn,fileOut,READ_LENGTH,OVERLAP_LENGTH,MIN_READ_LENGTH):
    with open(fileIn,'r') as f, open(fileOut,'w') as o:
        initial = 1
        for i in f.readlines():
            if initial:
                if i.startswith('>'):
                    initial = 0
                    readHeader = i.strip('\n')
                    readName = readHeader.split(' ')[0]
                    seq = ''
                continue
            if i.startswith('>'):
                res = overlapSplit(seq,readName,OVERLAP_LENGTH,READ_LENGTH,MIN_READ_LENGTH)
                o.write(res)
                readHeader = i.strip('\n')
                readName = readHeader.split(' ')[0]
                seq = ''
            else:
                seq += i.strip()
        else:
            res = overlapSplit(seq, readName, OVERLAP_LENGTH, READ_LENGTH,MIN_READ_LENGTH)
            o.write(res)
            pass


def options(argv):
    from argparse import ArgumentParser as AP
    usages = "python3 {} -i file_in -o file_out".format(argv[0])
    p = AP(usage=usages)
    p.add_argument("-i", dest="file_in", metavar="file_in", help="fasta file.")
    p.add_argument("-o", dest="file_out", metavar="file_out", help="split fasta file.")
    p.add_argument("-l", dest="length", metavar="[int]", help="splitted read length, default: 100", type=int, default=100)
    p.add_argument("-v", dest="overlap", metavar="[int]", help="splitted read overlap length, default: 80", type=int, default=80)
    p.add_argument("-m", dest="min", metavar="[int]",help="min read length, default: 30", default=30)
    if len(argv) == 1:
        p.print_help()
        exit(1)
    return p.parse_args(argv[1:])


if __name__ == '__main__':
    args = options(sys.argv)
    main(args.file_in,args.file_out,args.length,args.overlap,args.min)