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


def isoformExpCalc(stat,exp,out,out_taget):
    TotalExp = sum([exp[i] for i in exp])

    with open(os.path.join(out_taget,stat),'r',encoding='utf8') as f, open(out,'w') as o:
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


def transcriptExpCalc(cluster,exp, out,out_taget,stat,outputdir):
    TotalExp = sum([exp[i] for i in exp])
    dic={}
    with open(os.path.join(out_taget,cluster),'r') as f, open(out,'w') as o:
        for i in f.readlines():
            if not i.strip():
                continue
            if i.startswith('Chrom'):
                o.write(i.rstrip()+'\tTPM\tFLC\n')
            else:
                line = i.strip().split('\t')
                FLC = 0
                for id in line[-1].split(','):
                    if id in exp:
                        FLC += exp[id]
                TPM = float(FLC) / TotalExp * 1000000  # TPM FUNCTION
                o.write('{}\t{}\t{}\n'.format(i.strip(),TPM,FLC))
                trans_name=line[4].split(',')
                for tmp in trans_name:
                    if tmp not in dic:
                        dic[tmp]=[str(FLC),str(TPM)]
                    else:
                        print('error')
                        continue
    outfile=open(os.path.join(outputdir,stat+'1'),'w')     
    for line in open(os.path.join(out_taget,stat)):
        if line.startswith('ID'):
            outfile.write(line.rstrip()+'\tFLC\tTPM\n')
        else:
            newline=line.rstrip().split('\t')
            if newline[0] in dic:
                outfile.write(line.rstrip()+'\t'+'\t'.join(dic[newline[0]])+'\n')
            else:
                outfile.write(line.rstrip()+'\tNA\tNA\n')
                #print(newline)
                #print('key error')

def geneExpCalc(clusters,exp,out,out_taget):
    TotalExp = sum([exp[i] for i in exp])

    geneDict = {}
    for file in clusters:
        with open(os.path.join(out_taget,file),'r') as f:
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


def main(fa,prefix,outputdir,out_taget):
    stat = prefix + '.stat'
    #cluster = [prefix+'.SM.cluster',prefix+'.NC.cluster',prefix+'.FUSION.cluster',prefix+'.GENIC.cluster']
    cluster=[prefix+'.cluster.transcript']
    if not os.path.exists(os.path.join(out_taget,stat)):
        print('{} not found !'.format(stat))
        exit(11)
    if not os.path.exists(fa):
        print('{} not found !'.format(stat))
        exit(11)
    for i in cluster:
        if not os.path.exists(os.path.join(out_taget,i)):
            print('{} not found !'.format(i))
            exit(11)
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    expFile = exp(fa,os.path.join(outputdir,os.path.basename(prefix+'.readsExp.tsv')))
    expD = expParse(expFile,FL=1)
    isoformExpCalc(stat,expD,os.path.join(outputdir,os.path.basename(prefix+'.IsoformExp.tsv')),out_taget)
    for i in cluster:
        transcriptExpCalc(i,expD,os.path.join(outputdir,os.path.basename(i+'.exp')),out_taget,stat,outputdir)
    geneExpCalc(cluster,expD,os.path.join(outputdir,os.path.basename(prefix+'.geneExp.tsv')),out_taget)
    #merge(prefix, outputdir)
def merge(prefix,outputdir):
    dic={}
    for line in open(os.path.join(outputdir,prefix+'.SM.transcriptExp.tsv')):
        newline=line.rstrip().split('\t')
        t = newline[-1]
        newline[-1] = newline[-2]
        newline[-2] = t
        if newline[0] not in dic:
            dic[newline[0]]=[newline]
        else:
            dic[newline[0]].append(newline)
    for line in open(os.path.join(outputdir,prefix+'.NC.transcriptExp.tsv')):
        newline=line.rstrip().split('\t')
        t = newline[-1]
        newline[-1] = newline[-2]
        newline[-2] = t
        if newline[0] not in dic:
            dic[newline[0]]=[newline]
        else:
            dic[newline[0]].append(newline)
    for line in open(os.path.join(outputdir,prefix+'.GENIC.transcriptExp.tsv')):
        newline=line.rstrip().split('\t')
        t=newline[-1]
        newline[-1]=newline[-2]
        newline[-2]=t
        if newline[0] not in dic:
            dic[newline[0]]=[newline]
        else:
            dic[newline[0]].append(newline)
    dic1={}
    for line in open(os.path.join(outputdir,prefix+'.IsoformExp.tsv')):
        newline=line.rstrip().split('\t')
        if newline[3] not in dic1:
            dic1[newline[3]]=newline[5]
        else:
            continue
    outfile=open(os.path.join(outputdir,prefix+'.transcript.exp'),'w')
    outfile.write('Chrom\tGene\tTranscript\tReadsCount\tReadsID\tTPM\tFLC\n')
    for line in dic:
        for tmp in dic[line]:
            if line in dic1:
                outfile.write(dic1[line]+'\t'+'\t'.join(tmp)+'\n')
            else:
                chr1=line.split(':')
                outfile.write(chr1[0] + '\t' + '\t'.join(tmp) + '\n')
    outfile.close()


def option(argv):
    from argparse import ArgumentParser as AP
    usages = "python3 {} -f fa -i prefix -o outputdir".format(argv[0])
    p = AP(usage=usages)
    p.add_argument("-f", dest="fa", metavar="[fa]", help="fasta file.",required=True)
    p.add_argument("-i", dest="prefix", metavar="[prefix]", help="prefix of .annot.cluster.transcript",required=True)
    p.add_argument("-o", dest="out", metavar="[outputdir]", help="outputdir path", required=True)
    p.add_argument("-p", dest="out_taget", metavar="[outputdir_taget]", help="outputdir of taget path", required=True)
    if len(argv) == 1:
        p.print_help()
        exit(1)
    return p.parse_args(argv[1:])


if __name__ == '__main__':
    args = option(sys.argv)
    main(args.fa,args.prefix,args.out,args.out_taget)
