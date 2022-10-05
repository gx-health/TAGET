"""
Author: Zhang Chengsheng, @2020.03.06
建库，老脚本没了，被迫改一个更老的
"""

import os,sys
import re
import pickle


def dict_make(gtf,gtf_out=0):
    global db_dict
    db_dict = {}
    with open(gtf,'r') as f:
        for line in f.readlines():
            if line.startswith('#'):
                continue
            context_parse(line)
    db_patch()
    if not gtf_out:
        gtf_out = gtf + '.pickle'
    pickle_make(db_dict,gtf_out)


def pickle_make(dict_in,file_out):
    with open(file_out, 'wb') as o:
        pickle.dump(dict_in, o)


def load_pickle(file_in):
    return pickle.load(open(file_in, 'rb'))


def context_parse(line):
    global db_dict
    header = ['chr', 'db', 'record', 'start', 'end', 'tmp', 'strand', 'tmp', 'info']
    keywords = ['gene_id','transcript_id','gene_name','transcript_name','exon_number','gene_biotype']
    context = line.strip('\n').split('\t')
    record = context[2]
    gene_biotype = 'unknown'
    chr = context[0]
    start, end = int(context[3]), int(context[4])
    start, end = min(start, end), max(start, end)
    strand = context[6]  # +/-
    #strings = context[8]
    if chr not in db_dict:
        db_dict[chr] = {}
    if record == 'gene':
        gene_id = re_find_keyword(line, 'gene_name')  # gene_id  # for ENSEMBL & NCBI & GENCODE
        #gene_biotype = re_find_keyword(line, 'gene_biotype')  # for ENSEMBL
        if gene_id not in db_dict[chr]:
            db_dict[chr][gene_id] = [gene_biotype,strand,[start,end],{}]
        else:
            if not db_dict[chr][gene_id][2]:
                db_dict[chr][gene_id][2] = [start,end]
            if not db_dict[chr][gene_id][1]:
                db_dict[chr][gene_id][1] = strand
    if record == 'transcript':
        gene_id = re_find_keyword(line, 'gene_name')  # gene_id  # for ENSEMBL & NCBI & GENCODE
        #gene_biotype = re_find_keyword(line, 'gene_biotype')  # for ENSEMBL
        transcript_id = re_find_keyword(line, 'transcript_name')  # transcript_id transcript_name  for ENSEMBL & GENCODE
        #transcript_id = re_find_keyword(line, 'transcript_id')  # transcript_id transcript_name  for NCBI
        if gene_id not in db_dict[chr]:
            db_dict[chr][gene_id] = [gene_biotype,strand,[],{}]
        if transcript_id not in db_dict[chr][gene_id][3]:
            db_dict[chr][gene_id][3][transcript_id] = [strand,[start,end],{}]
        else:
            if not db_dict[chr][gene_id][3][transcript_id][1]:
                db_dict[chr][gene_id][3][transcript_id][1] = [start,end]
    if record == 'exon':
        gene_id = re_find_keyword(line, 'gene_name')  # gene_id  # for ENSEMBL & NCBI & GENCODE
        #gene_biotype = re_find_keyword(line, 'gene_biotype')  # for ENSEMBL
        transcript_id = re_find_keyword(line, 'transcript_name')  # transcript_id transcript_name  for ENSEMBL & GENCODE
        #transcript_id = re_find_keyword(line, 'transcript_id')  # transcript_id transcript_name   for NCBI
        exon_number = int(re_find_keyword(line, 'exon_number'))  # for ENSEMBL
        #exon_number = int(re_find_keyword_wtf_gencode(line))  # for GENCODE
        if gene_id not in db_dict[chr]:
            db_dict[chr][gene_id] = [gene_biotype,strand,[],{}]
        if transcript_id not in db_dict[chr][gene_id][3]:
            db_dict[chr][gene_id][3][transcript_id] = [strand,[],{}]
        if not list(db_dict[chr][gene_id][3][transcript_id][2]):
            exon_number = 1
        else:
            exon_number = max(list(db_dict[chr][gene_id][3][transcript_id][2]))+1
        if exon_number not in db_dict[chr][gene_id][3][transcript_id][2]:
            db_dict[chr][gene_id][3][transcript_id][2][exon_number] = [start,end]


def re_find_keyword_wtf_gencode(strings,keyword='exon_number'):
    return re.findall(keyword + ' (.*?);', strings)[0]


def re_find_keyword(strings,keyword='gene_name'):
    return re.findall(keyword + ' "(.*?)"', strings)[0]


def db_patch():
    def _get_transcript_region(dict):
        p = []
        for exon in dict:
            p.append(dict[exon][0])
            p.append(dict[exon][1])
        return [min(p),max(p)]

    def _get_gene_region(dict):
        p = []
        for transcript in dict:
            p.append(dict[transcript][1][0])
            p.append(dict[transcript][1][1])
        return [min(p),max(p)]

    global db_dict
    for chr in db_dict:
        for gene in db_dict[chr]:
            gene_exon_idx = {}  # 4
            gene_exon = []
            gene_ss = []  # 5
            gene_exon_num = []  # 6
            gene_transcipt_exon_idx = []  # 7
            for transcript in db_dict[chr][gene][3]:
                if not db_dict[chr][gene][3][transcript][1]:
                    db_dict[chr][gene][3][transcript][1] = _get_transcript_region(db_dict[chr][gene][3][transcript][2])
                for exon in db_dict[chr][gene][3][transcript][2]:
                    e1,e2 = db_dict[chr][gene][3][transcript][2][exon]
                    if db_dict[chr][gene][3][transcript][2][exon] not in gene_exon:
                        gene_exon.append(db_dict[chr][gene][3][transcript][2][exon])
                    if e1 not in gene_ss:
                        gene_ss.append(e1)
                    if e2 not in gene_ss:
                        gene_ss.append(e2)
                exon_num = len(db_dict[chr][gene][3][transcript][2])
                gene_exon_num.append(exon_num)
            gene_exon = sorted(gene_exon,key=lambda x:(x[0],x[1]))
            for _idx,i in enumerate(gene_exon):
                gene_exon_idx[_idx] = i
                for transcript in db_dict[chr][gene][3]:
                    transcript_exon_idx = []
                    for exon in db_dict[chr][gene][3][transcript][2]:
                        transcript_exon_idx.append(gene_exon.index(db_dict[chr][gene][3][transcript][2][exon]))
                    db_dict[chr][gene][3][transcript].append(transcript_exon_idx)
                    gene_transcipt_exon_idx.append(transcript_exon_idx)
            if not db_dict[chr][gene][2]:
                db_dict[chr][gene][2] = _get_gene_region(db_dict[chr][gene][3])
            db_dict[chr][gene].append(gene_exon_idx)
            db_dict[chr][gene].append(gene_ss)
            db_dict[chr][gene].append(gene_exon_num)
            db_dict[chr][gene].append(gene_transcipt_exon_idx)


def main():
    gtf = r'/public/source/share/zcs/data/GTF/hg38.ensembl.gtf'
    gtf_db = r'/public/source/share/zcs/data/GTF/hg38.ensembl.v20200306.1.pickle'
    dict_make(gtf, gtf_out=gtf_db)
    db_patch()


if __name__ == '__main__':
    main()