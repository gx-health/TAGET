"""
Author: Zhang Chengsheng, @2020.06.01
输入包含FL count的fa、软件的结果前缀，输出不同水平的三代测序TPM列表
"""

import os,sys


def exp(fa,out):
    """提取fa文件中的count, 格式：>transcript/0 full_length_coverage=4;length=9621;num_subreads=22"""
    with open(fa,'r') as f, open(out,'w') as o:
        o.write('#ID\tFLC\tLEN\tSRC\n')
        for i in f.readlines():
            if not i.strip():
                continue
            if i.startswith('>'):
                line = i.strip().split(' ')
                ID = line[0].lstrip('>')
                FLC = line[1].split(';')[0].split('=')[1]
                LEN = line[1].split(';')[1].split('=')[1]
                #SRC = line[1].split(';')[2].split('=')[1]
                #o.write('{}\t{}\t{}\t{}\n'.format(ID,FLC,LEN,SRC))
                o.write('{}\t{}\t{}\n'.format(ID,FLC,LEN))
    return out


def expParse(expFile,FL=1):
    """FL指提取的列"""
    expD = {}
    with open(expFile,'r') as f:
        for i in f.readlines():
            if not i.strip() or i.startswith('#'):
                continue
            line = i.strip().split('\t')
            if line[0] not in expD:
                expD[line[0]] = int(line[FL])
            else:
                expD[line[0]] += int(line[FL])
    return expD


def isoformExpCalc(stat,exp,out):
    TotalExp = sum([exp[i] for i in exp])

    with open(stat,'r',encoding='utf8') as f, open(out,'w') as o:
        for i in f.readlines():
            if not i.strip():
                continue
            if i.startswith('ID'):
                o.write('\t'.join(i.strip().split('\t')[:9])+'\tFLC\tTPM\n')
                continue
            line = i.strip().split('\t')
            if line[0] in exp:
                FLC = exp[line[0]]
                TPM = float(FLC)/TotalExp*1000000  #TPM FUNCTION
                o.write('{}\t{}\t{}\n'.format('\t'.join(line[:9]),FLC,TPM))


def transcriptExpCalc(cluster,exp, out):
    TotalExp = sum([exp[i] for i in exp])

    with open(cluster,'r') as f, open(out,'w') as o:
        for i in f.readlines():
            if not i.strip():
                continue
            line = i.strip().split('\t')
            FLC = 0
            for id in line[-1].split(','):
                if id in exp:
                    FLC += exp[id]
            TPM = float(FLC) / TotalExp * 1000000  # TPM FUNCTION
            o.write('{}\t{}\t{}\n'.format(i.strip(),FLC,TPM))


def geneExpCalc(clusters,exp,out):
    TotalExp = sum([exp[i] for i in exp])

    geneDict = {}
    for file in clusters:
        with open(file,'r') as f:
            for i in f.readlines():
                if not i.strip():
                    continue
                line = i.strip().split('\t')
                gene = line[0]
                ids = line[-1].split(',')
                FLC = 0
                for id in ids:
                    if id in exp:
                        FLC += exp[id]
                if gene not in geneDict:
                    geneDict[gene] = FLC
                else:
                    geneDict[gene] += FLC
    with open(out,'w') as o:
        o.write('Genes\tFLC\tTPM\n')
        for id in geneDict:
            TPM = float(geneDict[id]) / TotalExp * 1000000  # TPM FUNCTION
            o.write('{}\t{}\t{}\n'.format(id,geneDict[id],TPM))


def main(fa,prefix,outputdir):
    stat = prefix + '.stat'
    cluster = [prefix+'.SM.cluster',prefix+'.NC.cluster',prefix+'.FUSION.cluster',prefix+'.GENIC.cluster']
    if not os.path.exists(stat):
        print('{} not found !'.format(stat))
        exit(11)
    if not os.path.exists(fa):
        print('{} not found !'.format(stat))
        exit(11)
    for i in cluster:
        if not os.path.exists(i):
            print('{} not found !'.format(i))
            exit(11)
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    expFile = exp(fa,os.path.join(outputdir,os.path.basename(prefix+'.readsExp.tsv')))
    expD = expParse(expFile,FL=1)
    isoformExpCalc(stat,expD,os.path.join(outputdir,os.path.basename(prefix+'.IsoformExp.tsv')))
    for i in cluster:
        transcriptExpCalc(i,expD,os.path.join(outputdir,os.path.basename(i.replace('.cluster','')+'.transcriptExp.tsv')))
    geneExpCalc(cluster,expD,os.path.join(outputdir,os.path.basename(prefix+'.geneExp.tsv')))


def option(argv):
    from argparse import ArgumentParser as AP
    usages = "python3 {} -f fa -i prefix -o outputdir".format(argv[0])
    p = AP(usage=usages)
    p.add_argument("-f", dest="fa", metavar="[fa]", help="fasta file.",required=True)
    p.add_argument("-i", dest="prefix", metavar="[prefix]", help="prefix of *.NC.cluster",required=True)
    p.add_argument("-o", dest="out", metavar="[outputdir]", help="outputdir path", required=True)

    if len(argv) == 1:
        p.print_help()
        exit(1)
    return p.parse_args(argv[1:])


if __name__ == '__main__':
    args = option(sys.argv)
    main(args.fa,args.prefix,args.out)
