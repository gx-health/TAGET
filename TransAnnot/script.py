# -*- coding: utf-8 -*-
import os
import sys

def get_file(inputdir1,filename1):
    file1_stat = os.path.join(inputdir1,str(filename1 + '.annot.stat'))
    file1_junc = os.path.join(inputdir1,str(filename1 + '.annot.junction'))
    file1_anno = os.path.join(inputdir1,str(filename1 + '.annot.multiAnno'))
    file1_anno_file = os.path.join(inputdir1,str(filename1 + '.annot.multiAnno.filter'))
    trans_out_exp1 = os.path.join(inputdir1,str('merge\\' + filename1 + '.transcript.exp'))
    mkdir(os.path.join(inputdir1,str('merge')))
    
    return file1_stat,file1_junc,file1_anno,file1_anno_file,trans_out_exp1


def mkdir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)

def get_stat_info(file1_stat):
    dic_exp_T2 = {}
    dic_tpm_T2 = {}
    for row in open(file1_stat,encoding = 'utf-8'):
        line = row[:-1].split('\t')
        dic_exp_T2[line[0]] = line[-2]
        dic_tpm_T2[line[0]] = line[-1]
    
    return dic_exp_T2,dic_tpm_T2

def get_junc_info(file1_junc):
    dic_T2_con = {}
    for row in open(file1_junc):
        line = row.split('\t')
        if line[0] in dic_T2_con:
            dic_T2_con[line[0]].append(line[10])
        else:
            dic_T2_con[line[0]] = [line[10]]
    
    return dic_T2_con

def get_anno_info(file1_anno,file1_anno_file,dic_T2_con):
    list_T2_trans = []
    fo = open(file1_anno_file,'w')
    for row in open(file1_anno):
        line = row[:-1].split('\t')
        if '|' in line[3]:
            continue
        
        if ':' in line[1]:
            continue
        
        if line[2] == 'NA':
            if line[0] in dic_T2_con:
                if set(dic_T2_con[line[0]]) == {'canonical'}:
                    fo.write(row)
                    list_T2_trans.append(line[0])
                    continue
                else:
                    continue
            else:
                continue
        
        fo.write(row)
        list_T2_trans.append(line[0])
    
    fo.close()
    return list_T2_trans

def write_new_exp(file2_exp_file,out_exp1,list_T2_trans,dic_exp_T2,dic_tpm_T2):
    fo = open(out_exp1,'w')
    for row in open(file2_exp_file):
        line = row[:-1].split('\t')
        if line[0] == 'Chrom':
            fo.write(row)
            continue
        
        if line[4] == '':
            fo.write(row)
            continue
        
        list_new = []
        exp = 0
        tpm = 0
        for i in line[4].split(','):
            if i in list_T2_trans:
                exp = exp + int(dic_exp_T2[i])
                tpm = tpm + float(dic_tpm_T2[i])
                list_new.append(i)
        
        fo.write('\t'.join(line[0:3]) + '\t' + str(len(list_new)) + '\t' + ','.join(list_new) + '\t' + str(tpm) + '\t' + str(exp) + '\n')
    
    fo.close()

def runtime_file(inputdir2,filename2,file2_exp_file):
    file2_stat,file2_junc,file2_anno,file2_anno_file,trans_out_exp2 = get_file(inputdir2,filename2)
    dic_exp_T2,dic_tpm_T2 = get_stat_info(file2_stat)
    dic_T2_con = get_junc_info(file2_junc)
    list_T2_trans = get_anno_info(file2_anno,file2_anno_file,dic_T2_con)
    write_new_exp(file2_exp_file,trans_out_exp2,list_T2_trans,dic_exp_T2,dic_tpm_T2)
    return trans_out_exp2

def get_info(trans_out_exp2,dic_gene,list_trans):
    dic_1_t = {}
    for row in open(trans_out_exp2):
        if row.startswith('Chrom'):
            continue
        line = row[:-1].split('\t')
        if line[0] not in list_trans:
            list_trans.append(line[2])
        
        dic_1_t[line[2]] = int(line[6])
        dic_gene[line[2]] = line[0] + '\t' + line[1]
    
    return dic_gene,list_trans,dic_1_t

def write_merge(trans_out_exp1,trans_out_exp2,trans_outfile,filename1,filename2):
    dic_gene = {}
    list_trans = []
    dic_gene,list_trans,dic_1_t = get_info(trans_out_exp1,dic_gene,list_trans)
    dic_gene,list_trans,dic_2_t = get_info(trans_out_exp2,dic_gene,list_trans)
    
    fo = open(trans_outfile,'w')
    fo.write('Chrom\tGene\tTranscript\t' + filename1 + '\t' + filename2 + '\n')
    for i in sorted(list(set(list_trans))):
        fo.write(dic_gene[i] + '\t' + i + '\t')
        t1 = 0
        t2 = 0
        if i in dic_1_t:
            t1 = dic_1_t[i]
        
        if i in dic_2_t:
            t2 = dic_2_t[i]
        
        fo.write(str(t1) + '\t' + str(t2) + '\n')
    
    fo.close()


def write_merge_known(trans_outfile):
    outfile_known = trans_outfile + '.knwon'
    fo = open(outfile_known,'w')
    for row in open(trans_outfile):
        if row.startswith('Chrom'):
            fo.write(row)
            continue
        
        if ':' in row:
            continue
        
        if 'Genic' in row:
            continue
        
        if '\t0\t0' in row:
            continue
        
        if 'NNC' in row:
            continue
        
        if 'NIC' in row:
            continue
        
        fo.write(row)
    
    fo.close()
    return outfile_known


def runtime_gene(gtf_file,trans_outfile_known,filename1,filename2,gene_outfile):
    dic_gtf_site,gtf_site,dic_gene_trans = get_gtf_info(gtf_file)
    list_site = get_trans_merge_info(trans_outfile_known,dic_gtf_site)
    try:
        write_gene_merge(gene_outfile,trans_outfile_known,filename1,filename2)
        print('Gene合并文件成功生成。')
    except:
        print('Gene合并文件未生成。')

def get_gtf_info(gtf_file):
    dic_gtf_site = {}
    gtf_site = []
    dic_gene_trans = {}
    for row in open(gtf_file):
        line = row[:-1].split('\t')
        if line[0] == 'Chr':
            continue
        if line[4] not in dic_gtf_site:
            dic_gtf_site[line[4]] = [int(line[1])]
            dic_gtf_site[line[4]].append(int(line[2]))
        else:
            dic_gtf_site[line[4]].append(int(line[1]))
            dic_gtf_site[line[4]].append(int(line[2]))
        
        gene = '-'.join(line[4].split('-')[:-1])
        if gene not in dic_gene_trans:
            dic_gene_trans[gene] = [line[4]]
        else:
            dic_gene_trans[gene].append(line[4])
    
    for key in dic_gtf_site:
        dic_gtf_site[key] = sorted(dic_gtf_site[key])[1:-1]
        gtf_site.append(dic_gtf_site[key])
    
    return dic_gtf_site,gtf_site,dic_gene_trans


def get_trans_merge_info(trans_outfile_known,dic_gtf_site):
    list_site = {}
    for row in open(trans_outfile_known):
        line = row[:-1].split('\t')
        if 'Chrom' in row:
            continue
        
        if line[1] not in list_site:
            list_site[line[1]] = [dic_gtf_site[line[2]][1:-1]]
        else:
            if dic_gtf_site[line[2]][1:-1] in list_site[line[1]]:
                list_site[line[1]].append(dic_gtf_site[line[2]][1:-1])
            else:
                list_site[line[1]].append(dic_gtf_site[line[2]][1:-1])
    
    return list_site


def write_gene_other(gene_outfile):
    file1 = gene_outfile + '.nozero'
    file2 = gene_outfile + '.count2'
    fo = open(file1,'w')
    for row in open(gene_outfile):
        if '\t0' in row:
            continue
        
        fo.write(row)
    
    fo.close()
    
    fo = open(file2,'w')
    for row in open(gene_outfile):
        line = row[:-1].split('\t')
        if line[0] == 'Gene':
            fo.write(row)
            continue
        
        if ((int(line[1]) >= 2)&(int(line[2]) >= 2)):
            fo.write(row)
    
    fo.close()
    return file1,file2


def write_gene_merge(gene_outfile,trans_outfile_known,filename1,filename2):
    list_gene = []
    dic_gene_1 = {}
    dic_gene_2 = {}
    fo = open(gene_outfile,'w')
    for row in open(trans_outfile_known):
        line = row[:-1].split('\t')
        if line[0] == 'Chrom':
            continue
        
        if line[1] not in list_gene:
            list_gene.append(line[1])
            dic_gene_1[line[1]] = int(line[3])
            dic_gene_2[line[1]] = int(line[4])
        else:
            dic_gene_1[line[1]] = dic_gene_1[line[1]] + int(line[3])
            dic_gene_2[line[1]] = dic_gene_2[line[1]] + int(line[4])
        
    fo.write('Gene\t' + filename1 + '\t' + filename2 + '\n')
    for i in list_gene:
        fo.write(i + '\t' + str(dic_gene_1[i]) + '\t' + str(dic_gene_2[i]) + '\n')
    
    fo.close()

def get_info_from_config(config_file):
    for row in open(config_file,encoding = 'utf-8'):
        if row.startswith('#'):
            continue
        
        line = row.strip().split('\t')
        if line[0] == 'inputdir1':
            inputdir1 = line[1]
        elif line[0] == 'filename1':
            filename1 = line[1]
        elif line[0] == 'inputexp1':
            file1_exp_file = line[1]
        elif line[0] == 'inputdir2':
            inputdir2 = line[1]
        elif line[0] == 'filename2':
            filename2 = line[1]
        elif line[0] == 'inputexp2':
            file2_exp_file = line[1]
        elif line[0] == 'exon_gtf':
            gtf_file = line[1]
        elif line[0] == 'out_dir':
            trans_out_dir = line[1]
        elif line[0] == 'out_sample':
            out_sample = line[1]
        elif line[0] == 'diff_script_file':
            diff_script_file = line[1]
        elif line[0] == 'diff':
            DEG_runtime = line[1]
    
    diff_dir = os.path.join(trans_out_dir,'diff')
    trans_outfile_pre = os.path.join(trans_out_dir,out_sample)
    trans_outfile = trans_outfile_pre + '.trans.exp'
    gene_outfile = trans_outfile_pre + '.gene.exp.known'
    return inputdir1,filename1,file1_exp_file,inputdir2,filename2,file2_exp_file,gtf_file,out_sample,diff_script_file,diff_dir,trans_outfile,gene_outfile,DEG_runtime



def main():
    config_file = sys.argv[1]
    inputdir1,filename1,file1_exp_file,inputdir2,filename2,file2_exp_file,gtf_file,out_sample,diff_script_file,diff_dir,trans_outfile,gene_outfile,DEG_runtime = get_info_from_config(config_file)
    trans_out_exp2 = runtime_file(inputdir2,filename2,file2_exp_file)
    trans_out_exp1 = runtime_file(inputdir1,filename1,file1_exp_file)
    try:
        write_merge(trans_out_exp1,trans_out_exp2,trans_outfile,filename1,filename2)
        print('transcript合并文件成功生成。')
    except:
        print('transcript合并文件未生成。')
    try:
        trans_outfile_known = write_merge_known(trans_outfile)
        print('transcript合并过滤文件成功生成。')
        if DEG_runtime:
            trans_outfile1_diff = os.path.join(diff_dir, str(out_sample + '.trans.exp.knwon'))
            mkdir(diff_dir)
            mkdir(trans_outfile1_diff)
            try:
                print('transcript合并过滤文件差异分析开始运行.')
                diff(trans_outfile_known,trans_outfile1_diff,diff_script_file)
            except:
                print('transcript合并过滤文件差异分析未运行.')
            
    except:
        print('transcript合并过滤文件未生成。')
    runtime_gene(gtf_file,trans_outfile_known,filename1,filename2,gene_outfile)
    try :
        gene_out_nozero,gene_out_count2 =  write_gene_other(gene_outfile)
        print('Gene合并过滤文件（去0、count>=2）成功生成。')
        if DEG_runtime:
            gene_outfile1_diff = os.path.join(diff_dir, str(out_sample + '.gene.exp.known.nozero'))
            gene_outfile2_diff = os.path.join(diff_dir, str(out_sample + '.gene.exp.known.count2'))
            mkdir(diff_dir)
            mkdir(gene_outfile1_diff)
            mkdir(gene_outfile2_diff)
            try:
                print('Gene合并过滤文件（去0）差异分析开始运行.')
                diff(gene_out_nozero,gene_outfile1_diff,diff_script_file)
            except:
                print('Gene合并过滤文件（去0）差异分析未运行.')
            try:
                print('Gene合并过滤文件（count>=2）差异分析开始运行.')
                diff(gene_out_count2,gene_outfile2_diff,diff_script_file)
            except:
                print('Gene合并过滤文件（count>=2）差异分析未运行.')
    except:
        print('Gene合并过滤文件（去0、count>=2）未生成。')

def diff(inputfile,outfile,diff_script_file):
    if 'trans' in inputfile:
        diff1 = diff_script_file.replace('.r','_trans.r')
        cmd = 'Rscript "' + diff1 + '" "' + inputfile + '" "' + outfile + '"'
    else:
        cmd = 'Rscript "' + diff_script_file + '" "' + inputfile + '" "' + outfile + '"'
    os.system(cmd)

if __name__ == '__main__' :
    main()
    print ("Finish!")


