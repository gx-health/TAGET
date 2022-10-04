"""

"""

import os,sys,pickle


def load_pickle(file_in):
    return pickle.load(open(file_in, 'rb'))


def pickle_make(dict_in,file_out):
    with open(file_out, 'wb') as o:
        pickle.dump(dict_in, o)


def bed_read(bed):
    isoforms_bed = {}
    with open(bed,'r') as f:
        for idx,line in enumerate(f.readlines()):
            if not line.strip():
                continue
            id = line.split('\t')[3]
            if not idx:
                id_temp = id
                bed_string = line
                continue
            if id == id_temp:
                bed_string += line
            else:
                if id_temp in isoforms_bed:
                    print('error: id repeat in bed files: {}'.format(id_temp))
                isoforms_bed[id_temp] = bed_string
                id_temp = id
                bed_string = line
        else:
            isoforms_bed[id_temp] = bed_string
    return isoforms_bed


def check_overlap(ref, gene, length=0):
    flag = 0
    overlap_length = 0
    if ref[0] <= gene[0]:
        if ref[1] <= gene[0]:
            flag = 0  # 不挨着ref在左
        elif ref[1] >= gene[1]:
            flag = 2  # ref 包含gene
            overlap_length = abs(gene[1] - gene[0] + 1)
        else:
            flag = 1  # overlap
            overlap_length = abs(ref[1] - gene[0] + 1)
    else:
        if ref[0] >= gene[1]:
            flag = 0  # 不挨着ref在右
        elif ref[1] <= gene[1]:
            flag = 2  # gene包含ref
            overlap_length = abs(ref[1] - ref[0] + 1)
        else:
            flag = 1  # overlap
            overlap_length = abs(gene[1] - ref[0] + 1)
    if length:
        return overlap_length
    return flag


def trans(seq, transform=True,reverse=False):
    dict1 = {"A":"T",
             "T":"A",
             "G":"C",
             "C":"G",
             "a":"t",
             "t":"a",
             "c":"g",
             "g":"c",
             "N":"N",
             "n":"n"}
    new_seq = ""
    if transform:
        for i in seq:
            if i not in dict1:
                new_seq += i
            else:
                new_seq += dict1[i]
    else:
        new_seq = seq
    if reverse:
        new_seq = new_seq[::-1]
    return new_seq


def refseq_extract(site_chr,site_site,strand,fbuffer,faidxDict):
    faidx_dict = faidxDict
    start,end = int(site_site[0]),int(site_site[1])
    if site_chr in faidx_dict:
        chrom = site_chr
    elif site_chr.lstrip('chr') in faidx_dict:
        chrom = site_chr.lstrip('chr')
    elif 'chr'+site_chr in faidx_dict:
        chrom = 'chr'+site_chr
    else:
        return 0
    offset = int(faidx_dict[chrom][1])
    line = int(faidx_dict[chrom][2])
    size = int(faidx_dict[chrom][3])
    location = offset + int(start / line) + start - 1
    length = (int(end / line) - int(start / line)) * (size - line) + end - start + 1
    fbuffer.seek(location, 0)
    sequence = fbuffer.read(length)
    if sequence.startswith('\n'):
        fbuffer.seek(location-1,0)
        s1 = fbuffer.read(1)
        sequence = s1+sequence
    sequence = sequence.upper().replace('\n','')
    if strand not in ['1',1,'+']:
        sequence = trans(sequence,transform=True,reverse=True)
    return sequence


def refseqIdx(faidx):
    faidx_dict = {}
    with open(faidx, 'r') as o:
        for i in o.readlines():
            faidx_dict[i.strip().split('\t')[0]] = i.strip().split('\t')[1:]
    return faidx_dict


def ss_compare(ref,target):
    MAX_DRAFT = 10
    target_start, target_end = [i[0] for i in target], [i[1] for i in target]
    target_range = [min(target_start), max(target_end)]
    if not ref:
        return 0,target,target_range
    ref_start,ref_end = [i[0] for i in ref],[i[1] for i in ref]
    ref_range = [min(ref_start),max(ref_end)]
    tk,rk=0,0
    match = []
    for t_idx in range(len(target_start)):
        start_match,end_match = 0,0
        if not check_overlap(ref_range,target[t_idx]):
            continue
        for j in range(len(ref_start[rk:])):
            r_idx = rk+j
            if abs(target_start[t_idx]-ref_start[r_idx]) < MAX_DRAFT:
                start_match = 1
            if abs(target_end[t_idx]-ref_end[r_idx]) < MAX_DRAFT:
                end_match = 1
            if start_match or end_match:
                rk += j
                break
        if not start_match+end_match:
            match = []
            break
        match.append(start_match+end_match)
    same = 0
    if len(match) > 1:
        if 1 not in match[1:-1]:
            same = 1
            if abs(ref_range[-1]-ref_range[0]) > abs(target_range[-1]-target_range[0]):
                target = ref
                target_range = ref_range
    return same,target,target_range


def runSplit(c):
    import Long2Short
    faSplit = os.path.join(c['OUTPUT_DIR'],c['SAMPLE_UNIQUE_NAME']+c['FA_SPLIT_SUFFIX'])
    if os.path.exists(faSplit):
        print('find fasplit file, use it!')
        return
    Long2Short.main(c['FASTA'],faSplit,int(c['READ_LENGTH']),int(c['READ_OVERLAP']),int(c['MIN_READ_LENGTH']))


def runSam2Bed(sam,bed,Type,process,readLength,readOverlap):
    import Sam2Bed
    Sam2Bed.main(sam,bed,type=Type,process=int(process),READ_LENGTH=int(readLength),OVERLAP_LENGTH=int(readOverlap))


def runMinimap2(c):
    print('run minimap2')
    sam = os.path.join(c['OUTPUT_DIR'],c['SAMPLE_UNIQUE_NAME']+c['SAM_SUFFIX_MINIMAP2'])
    samSort = os.path.join(c['OUTPUT_DIR'],c['SAMPLE_UNIQUE_NAME']+c['SAM_SORT_SUFFIX_MINIMAP2'])
    bed = os.path.join(c['OUTPUT_DIR'], c['SAMPLE_UNIQUE_NAME'] + c['BED_SUFFIX_MINIMAP2'])
    if os.path.exists(bed):
        print('find minimap2 result, use it!')
        return bed
    if os.path.exists(samSort):
        print('find minimap2 sorted sam file, use it!')
        runSam2Bed(samSort, bed, 'long', c['PROCESS'], c['READ_LENGTH'], c['READ_OVERLAP'])
        return bed
    cmd1 = '{} -ax splice -t {} -uf --secondary=no -C5 {} {} > {}'.format(c['MINIMAP2'],c['PROCESS'],c['GENOME_FA'],c['FASTA'],sam)
    cmd2 = '{} sort -n -o {} {}'.format(c['SAMTOOLS'],samSort,sam)
    os.system(cmd1)
    os.system(cmd2)
    os.remove(sam)
    runSam2Bed(samSort,bed,'long',c['PROCESS'],c['READ_LENGTH'],c['READ_OVERLAP'])
    return bed


def runHisat2(c):
    print('run hisat2')
    faSplit = os.path.join(c['OUTPUT_DIR'],c['SAMPLE_UNIQUE_NAME']+c['FA_SPLIT_SUFFIX'])
    sam = os.path.join(c['OUTPUT_DIR'],c['SAMPLE_UNIQUE_NAME']+c['SAM_SUFFIX_HISAT2'])
    samSort = os.path.join(c['OUTPUT_DIR'],c['SAMPLE_UNIQUE_NAME']+c['SAM_SORT_SUFFIX_HISAT2'])
    bed = os.path.join(c['OUTPUT_DIR'], c['SAMPLE_UNIQUE_NAME'] + c['BED_SUFFIX_HISAT2'])
    if os.path.exists(bed):
        print('find hisat2 result, use it!')
        return bed
    if os.path.exists(samSort):
        print('find hisat2 sorted sam file, use it!')
        runSam2Bed(samSort, bed, 'long', c['PROCESS'], c['READ_LENGTH'], c['READ_OVERLAP'])
        return bed
    cmd1 = '{} -p {} -f --score-min L,0,-0.8 -x {} -U {} -S {}'.format(c['HISAT2'],c['PROCESS'],c['HISAT2_INDEX'],faSplit,sam)
    cmd2 = '{} sort -n -o {} {}'.format(c['SAMTOOLS'],samSort,sam)
    os.system(cmd1)
    os.system(cmd2)
    os.remove(sam)
    runSam2Bed(samSort, bed, 'short', c['PROCESS'], c['READ_LENGTH'], c['READ_OVERLAP'])
    return bed


def runGmap(c):
    print('run gmap')
    sam = os.path.join(c['OUTPUT_DIR'], c['SAMPLE_UNIQUE_NAME'] + c['SAM_SUFFIX_GMAP'])
    samSort = os.path.join(c['OUTPUT_DIR'], c['SAMPLE_UNIQUE_NAME'] + c['SAM_SORT_SUFFIX_GMAP'])
    bed = os.path.join(c['OUTPUT_DIR'], c['SAMPLE_UNIQUE_NAME'] + c['BED_SUFFIX_GMAP'])
    if os.path.exists(bed):
        print('find gmap result, use it!')
        return bed
    if os.path.exists(samSort):
        print('find gmap sorted sam file, use it!')
        runSam2Bed(samSort, bed, 'long', c['PROCESS'], c['READ_LENGTH'], c['READ_OVERLAP'])
        return bed
    cmd1 = ''
    cmd2 = '{} sort -n -o {} {}'.format(c['SAMTOOLS'], samSort, sam)
    #os.system(cmd1)
    #os.system(cmd2)
    #os.remove(sam)
    #runSam2Bed(samSort, bed, 'long', c['PROCESS'], c['READ_LENGTH'], c['READ_OVERLAP'])
    return bed

