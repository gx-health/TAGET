# -*- coding: utf-8 -*-
"""
Spyder Editor
https://blog.csdn.net/zhu_si_tao/article/details/71079842?utm_medium=distribute.pc_relevant.none-task-blog-baidujs-2
This is a temporary script file.
"""

import pandas as pd
import os
import sys
import numpy as np
from scipy.stats import chi2_contingency,fisher_exact
#from statsmodels.sandbox.stats.multicomp import multipletests
#from FisherExact import fisher_exact
#import FisherExact
import argparse
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
from rpy2.robjects.vectors import FloatVector
import time
stats=importr('stats')
stats_mod=importr('statmod')
def read_expression(filename,filter_list={},filter_expression=10):
    dic={}
    i=0
    j=0
    print('filter expression:',filter_expression)
    if filter_list!={}:
        for line in open(filename):
            if not line.startswith('Gene') and not line.startswith('Chrom'):
                newline=line.rstrip().split('\t')
                #print(newline)
                i=i+1
                #print(i)
                if len(newline)>4:
                    if newline[1] in filter_list:
                        tmp=[float(t) for t in newline[3:]]
                    else:
                        tmp=[0]
                else:
                    if newline[0] in filter_list:
                        tmp=[float(t) for t in newline[2:]]
                    else:
                        tmp=[0]
                #print(tmp)
                if max(tmp)<filter_expression:
                    #i=i+1
                    #print(i,newline)
                    continue
                else:
                    j+=1
                   
                    if len(newline)>4:
                        if newline[1] not in dic:

                            dic[newline[1]]=[[newline[2]]+[float(t) for t in newline[3:]]]
                        else:
                            dic[newline[1]].append([newline[2]]+[float(t) for t in newline[3:]])
                    else:
                        if newline[0] not in dic:
                            dic[newline[0]]=[[newline[1]]+[float(t) for t in newline[2:]]]
                        else:
                            dic[newline[0]].append([newline[1]]+[float(t) for t in newline[2:]])
            else:
                continue
    else:
        for line in open(filename):
            if not line.startswith('Gene') and not line.startswith('Chrom'):
                newline=line.rstrip().split('\t')
                if len(newline)>4:
                    tmp=[float(t) for t in newline[3:]]
                else:
                    tmp=[float(t) for t in newline[2:]]
                if max(tmp)<filter_expression:
                    continue
                else:
                
                    if newline[0] not in dic:
                        dic[newline[0]]=[[newline[1]]+[float(t) for t in newline[2:]]]
                    else:
                        dic[newline[0]].append([newline[1]]+[float(t) for t in newline[2:]])
            else:
                continue
    print(len(dic.keys()),' numbers of gene read')
    outfile=open('test.txt','w')
    i=0
    for key in dic:
        for line in dic[key]:
            #i=i+1
            #print(i)
            outfile.write(key+'\t'+line[0]+'\t'+str(line[1])+'\t'+str(line[2])+'\n') 
    outfile.close()   
    return dic
def read_gene_expression(filename):
    dic={}

    for line in open(filename):
        if not line.startswith('Gene'):
            newline=line.rstrip().split('\t')
            if newline[0] not in dic:
                dic[newline[0]]=[float(t) for t in newline[1:]]
            else:
                new_list=[i for i in range(len(newline[1:]))]
                
                #print('repeat name:' ,newline[0])
                for i in range(len(newline[1:])):
                    tmp=newline[1:]
                    #print(tmp,i,tmp[i])
                    new_list[i]=dic[newline[0]][i]+float(tmp[i])
                dic[newline[0]]=new_list
    count=0
    for key in dic:
        count+=1
    print('total gene number:',count)
    return dic

def filter_low_express_gene(dic,low_express=1,mode=1):
    dic_new={}
    print(low_express)
    if mode==1:
        tmp_key_list=[]
        for key in dic:
            if max(dic[key])<=low_express:
                tmp_key_list.append(key)
            else:
                continue
        for key in dic:
            if key not in tmp_key_list:
                dic_new[key]=dic[key]
            
    elif mode==2:
        tmp_key_list=[]
        for key in tmp_key_list:
            if max(dic[key])<=low_express and min(dic[key])>0:
                tmp_key_list.append(key)
            else:
                continue
        for key in dic:
            if key not in tmp_key_list:
                dic_new[key]=dic[key]
    else:
        pass
    count=0
    for key in dic_new:
        count+=1
    print('total gene number after filtering:',count,len(tmp_key_list))
    return dic_new
        
        
            

def gene_dataframe(dic,low_express_transcript_filter=10,filter_genic='TRUE'):
    gene_dic={}
    print(len(dic.keys()),' number of genes will be process')
    for key in dic:
        column=['tumor','normal']
        index_name=[]
        data_list=[]
        if filter_genic=='FALSE':
            for tmp in dic[key]:
                if max(int(tmp[1]),int(tmp2))<=low_express_transcript_filter:
                    continue
                else:
                    index_name.append(tmp[0])
                    data_list.append(tmp[1:])
                    #print(index_name)
                    #print(data_list)
        else:
            for tmp in dic[key]:
                if 'Genic' in tmp[0]:
                    #print('done')
                    continue
                else:
                   
                    if max(int(tmp[1]),int(tmp[2]))<=low_express_transcript_filter:
                        continue
                    else:
                        index_name.append(tmp[0])
                        data_list.append(tmp[1:])
        num_df=pd.DataFrame(data=data_list,index=index_name,columns=column)
        #print(num_df)
        if len(index_name)>1:
            gene_dic[key]=num_df
            num_df.to_csv('pandas.txt', header=None,sep='\t', mode='a')
        #else:
        #    print(num_df)
   
        
            
    print(len(gene_dic.keys()),' numbers of genes were generated data.frame ')
    return gene_dic
def chi2(gene_dic,eps=0.00001,method='fisher'):
    dic={}
    flag=0
    wrong=0
    for key in gene_dic:
        #print(gene_dic[key])
        a=gene_dic[key].replace(0,eps)
        #print(gene_dic[key])
        if method=='fisher':
            #odds_ratio,pvalue=fisher_exact(a)
            a=gene_dic[key]
            try:
                pvalue=fisher_exact_test(a)
                #pvalue1=FisherExact.fisher_exact(a)
                dic[key]=pvalue[0]
                if pvalue[1]==1:
                    flag+=1
                #print(pvalue,pvalue1)
            except:
                wrong+=1
                continue
                #print('error')
                #print(key,a)
        else:
            try:
                p=chi2_contingency(a)
                dic[key]=p[1]
            except:
                print(key,a)
            
    print('fisher:',flag)
    gene_list=[]
    p_value=[]
    for key in dic:
        gene_list.append(key)
        p_value.append(dic[key])
    #print(p_value)
    FDR=p_adjust_bh(p_value)
    FDR1 = stats.p_adjust(FloatVector(p_value), method = 'bonferroni')
    FDR2= stats.p_adjust(FloatVector(p_value), method = 'BY')
    #p_adjusted = multipletests(p_value, method='bonferroni')
    #print(FDR)
    dic=dict(zip(gene_list,FDR))
    dic_result={}
    outfile=open('adjust_p_value.txt','w')
    #print(len(gene_list),len(p_value),len(FDR),len(p_adjusted[1]))
    for i in range(len(gene_list)):
        outfile.write(gene_list[i]+'\t'+str(p_value[i])+'\t'+str(FDR[i])+'\t'+str(FDR1[i])+'\t'+str(FDR2[i])+'\n')
        dic_result[gene_list[i]]=[p_value[i],FDR[i]]
    dic_tmp={}
    #for key in dic:
    #    dic_tmp[key]=[diic[key]]+dic_result[key]
    diff=0
    same=0
    for key in dic:
        if dic[key]<=0.05:
            diff+=1
        else:
            same+=1
    print(same,diff,diff/(same+diff))
    switch={}
    not_switch={}
    switch_diff=0
    switch_df={}
    for key in dic:
        df=gene_dic[key]
        trans1=df['tumor'].idxmax()
        trans2=df['normal'].idxmax()
        if trans1==trans2:
            switch[key]=[dic[key],trans1,trans2]
        else:
            not_switch[key]=dic[key]
            if float(dic[key])<=0.05:
                if 'Genic' not in trans1 or 'Genic' not in trans2: 
                #if ('NIC' in trans1 or 'NNC' in trans1) and ('NIC' not in trans1 and 'NNC' not in trans1):
                    switch_diff+=1
                    switch_df[key]=[df,trans1,trans2]
                    #switch[key]=[dic[key],trans1,trans2]
                else:
                    switch[key]=[dic[key],trans1,trans2]
                #print(df)
            else:
                switch[key]=[dic[key],trans1,trans2]
                
    print(switch_diff,switch_diff/(same+diff))
    print(len(switch.keys()),' number of gene had been analysis')
    return [dic[key],switch,not_switch,switch_df,dic_result]
def total_usage(dic):
    usage={}
    for key in dic:
        df=dic[key]
        b=np.sum(df,axis=0)
        b=b.replace(0,0.0001)
        ratio=df/b
        total=np.sum(abs(ratio['tumor']-ratio['normal']))*100*0.5
        usage[key]=total
    return usage

def gene_DE(filename):
    dic={}
    for line in open(filename):
        if not line.startswith('"Gene'):
            newline=line.rstrip().split('\t')
            if newline[0] not in dic:
                #print(newline)
                dic[newline[0]]=[newline[4],newline[6]]
    return dic
def plot_vacinno(DE_dic,total_usage_dic,change_dic):
    result={}
    for key in total_usage_dic:
        if key in DE_dic:
            if key in change_dic:
                result[key]=[str(total_usage_dic[key]),str(DE_dic[key][0]),'1']
            else:
                result[key]=[str(total_usage_dic[key]),str(DE_dic[key][0]),'0']
    outfile=open('test.txt','w')
    for key in result:
        outfile.write(key+'\t'+'\t'.join(result[key])+'\n')
    outfile.close()
def p_adjust_bh(p):
    p=np.asfarray(p)
    by_descend=p.argsort()[::-1]
    by_orig=by_descend.argsort()
    steps=float(len(p))/np.arange(len(p),0,-1)
    q=np.minimum(1,np.minimum.accumulate(steps*p[by_descend]))
    return q[by_orig]

def outfile(p,usage_dic,DE_dic,output):
    change_dic=p[-2]
    switch=p[1]
    outfile=open(output+'_DIU_result.txt','w')
    outfile.write('Gene\tp_value\tp_adjust\ttotal_usage\tchange major isoform\ttumor_major\tnormal_major\tlogFC\tDE_pvalue\n')
    for key in p[-1]:
        if key in DE_dic:
            outfile.write(key+'\t')
            for line in p[-1][key]:
                outfile.write(str(line)+'\t')
            outfile.write(str(usage_dic[key])+'\t')
            if key in change_dic:
                outfile.write('1'+'\t'+change_dic[key][-2]+'\t'+change_dic[key][-1]+'\t')
            else:
                outfile.write('0'+'\t'+switch[key][-2]+'\t'+switch[key][1]+'\t')
            outfile.write(str(DE_dic[key][0])+'\t'+str(DE_dic[key][1])+'\n')
        else:
            outfile.write(key+'\t')
            for line in p[-1][key]:
                outfile.write(str(line)+'\t')
            outfile.write(str(usage_dic[key])+'\t')
            if key in change_dic:
                outfile.write('1'+'\t'+change_dic[key][-2]+'\t'+change_dic[key][-1]+'\t')
            else:
                outfile.write('0'+'\t'+switch[key][-2]+'\t'+switch[key][1]+'\t')
            outfile.write(str(0)+'\t'+str(1)+'\n')
    outfile.close()

def filter_low_transcript(gene_dic,filter_tpm_ratio=0.1,low_express=5):
    gene_dic_filter={}
    print(len(gene_dic.keys()),' numbers of gene will be fiter',low_express)
    for key in gene_dic:
        df=gene_dic[key]
        #if max(np.sum(df,axis=0))<=low_express:
        #    continue
        #else:
        ratio=df/np.sum(df,axis=0).replace(0,0.0001)
        df2=ratio.copy()
        df2.drop(df2[(df2.tumor < filter_tpm_ratio) & (df2.normal<filter_tpm_ratio)].index,inplace=True)
        row_index=[]
        for index, row in df2.iterrows():
            row_index.append(index)
        df_new=df.loc[row_index]
        if np.shape(df_new)[0]>1:
            gene_dic_filter[key]=df_new
    print(len(gene_dic_filter.keys()),' numbers of gene after fitering',low_express)
    return gene_dic_filter            
def fisher_exact_test(df):
    a=FloatVector(df['tumor'])
    b=FloatVector(df['normal'])
    c=robjects.r['cbind'](a,b)
    tag=0
    #print(df)
    try:
        d=robjects.r['fisher.test'](c)
        pvalue=d[0][0]
        #fisher+=1
        tag=1
        #print('fisher:',fisher)
    except:
        #print('Fisher exact fail,exe chiq2')
        a=df.replace(0,0.00001)
        p=chi2_contingency(a)
        pvalue=p[1]
        pvalue=1
        #print('chisq:done')
        
    return [pvalue,tag]        
def compute_gene_express(filename,output):
    dic={}
    for line in open(filename):
        if not line.startswith('Gene') and not line.startswith('Chrom'):
            newline=line.rstrip().split('\t')
            if len(newline)>4:
                if newline[1] not in dic:
                    dic[newline[1]]=[float(newline[3]),float(newline[4])]
                else:
                    dic[newline[1]][0]=dic[newline[1]][0]+float(newline[3])
                    dic[newline[1]][1]=dic[newline[1]][1]+float(newline[4])
            else:
                if newline[0] not in dic:
                    dic[newline[0]]=[float(newline[2]),float(newline[3])]
                else:
                    dic[newline[0]][0]=dic[newline[0]][0]+float(newline[2])
                    dic[newline[0]][1]=dic[newline[0]][1]+float(newline[3])
    #outfile=open('test1.txt','w')
    #for key in dic:
    #    outfile.write(key+'\t'+str(dic[key][0])+'\t'+str(dic[key][1])+'\n')
    #outfile.close()
    column=['tumor','normal']
    index_name=[]
    data_list=[]
    for key in dic:
        index_name.append(key)
        data_list.append(dic[key])
    num_df=pd.DataFrame(data=data_list,index=index_name,columns=column)
    a=FloatVector(num_df['tumor'])
    b=FloatVector(num_df['normal'])
    c=robjects.r['sage.test'](a,b,sum(num_df['tumor']),sum(num_df['normal']))
    FDR1 = stats.p_adjust(FloatVector(c), method = 'fdr')
    num_df['p-value']=c
    num_df['FDR']=FDR1
    num_df.to_csv(output+'_DGE_result.txt', sep='\t')
        
    return dic
        
def compute_DTE(filename,output,filter_expression=1):
    dic={}
    for line in open(filename):
        if not line.startswith('Gene') and not line.startswith('Chrom'):
            newline=line.rstrip().split('\t')
            if len(newline)>4:
                tmp=[float(t) for t in newline[3:]]
            else:
                tmp=[float(t) for t in newline[2:]]
            if max(tmp)<filter_expression:
                    continue
            else:

                if newline[0] not in dic:
                    dic[newline[0]]=[[newline[1]]+[float(t) for t in newline[2:]]]
                else:
                    dic[newline[0]].append([newline[1]]+[float(t) for t in newline[2:]])
        else:
            continue

    column=['transcript','tumor','normal']
    index_name=[]
    data_list=[]
    for key in dic:
        column=['transcript','tumor','normal']
        #index_name=[]
        #data_list=[]
        for line in dic[key]:

            if 'NIC' in line[0] or 'NNC' in line[0] or 'Genic' in line[0]:
                pass
            else:
                index_name.append(key)
                data_list.append(line)
    num_df=pd.DataFrame(data=data_list,index=index_name,columns=column)
    a=FloatVector(num_df['tumor'])
    b=FloatVector(num_df['normal'])
    c=robjects.r['sage.test'](a,b,sum(num_df['tumor']),sum(num_df['normal']))
    FDR1 = stats.p_adjust(FloatVector(c), method = 'fdr')
    num_df['p-value']=c
    num_df['FDR']=FDR1
    num_df.to_csv(output+'_DTE_result.txt', sep='\t')
    
if __name__ == '__main__':
    start = time.time()
    parser = argparse.ArgumentParser()
    parser.description='please enter three files transcirpt.exp, gene.exp and DE express '
    parser.add_argument("-t", "--transcript", help="this is parameter transcript.exp")
    parser.add_argument("-g", "--gene", help="this is parameter gene.exp",default='')
    parser.add_argument("-d", "--DE", help="this is parameter diffierential expression",default='')
    parser.add_argument("-o", "--output_priex", help="this is parameter output_preix",default='./test')
    parser.add_argument("-r", "--ratio", help="this is paramete of filter low express transcirpt",type=float, default="0.05")
    parser.add_argument("-p", "--gene_ratio", help="this is paramete of filter low express gene",type=float, default="50")
    parser.add_argument("-i", "--filter_transcirpt_expression", help="this is paramete of filter low express transcript",type=float, default="10")
    args = parser.parse_args()
    
    #filename=r'G:\express\new_expression\LL_T_transcript.exp'
    #gene_filename=r'G:\express\new_expression\LL_T_gene.exp'
    #DE_filename=r'G:\express\DEGs\LL_1\LL_1_DEG_50.txt'
    filename=args.transcript
    gene_filename=args.gene
    DE_filename=args.DE
    outfile_prex=args.output_priex
    low_express_transcript_ratio=float(args.ratio)
    low_express_gene_filter=float(args.gene_ratio)
    low_express_transcript_filter=float(args.filter_transcirpt_expression)
    print('parameter: ',low_express_transcript_ratio,low_express_gene_filter)
    print('gene_express is doing')
    if gene_filename!='':
        gene_express=read_gene_expression(gene_filename,outfile_prex)
    else:
        gene_express=compute_gene_express(filename,outfile_prex)
    print('gene_express is done')
    print('filter_low_express_gene is doing')
    tmp_dic=filter_low_express_gene(gene_express,low_express=low_express_gene_filter,mode=1)
    print('filter_low_express_gene is done')
    print('transcript_express is doing')
    dic= read_expression(filename,tmp_dic,filter_expression=low_express_transcript_filter)
    print('transcript_express is done')
    print('generate dataframe is doing')
    gene_dic=gene_dataframe(dic,low_express_transcript_filter,'TRUE')
    print('generate dataframe is done')
    if low_express_transcript_ratio>0:
        print('filter low transcirpt is doing')
        gene_dic_filter=filter_low_transcript(gene_dic,filter_tpm_ratio=low_express_transcript_ratio,low_express=low_express_transcript_filter) 
        print('filter low transcirpt is done')
    else:
        gene_dic_filter=gene_dic.copy()
    print('fisher exact is doing')
    p=chi2(gene_dic_filter,method='fisher')
    print('fisher exact is done')
    usage_dic=total_usage(gene_dic_filter)
    if DE_filename=='':
        DE_dic={}
    else:
        DE_dic=gene_DE(DE_filename)
    #plot_vacinno(DE_dic,usage_dic,p[-2])
    print('DTE is doing')
    compute_DTE(filename,outfile_prex,filter_expression=low_express_transcript_filter)
    outfile(p,usage_dic,DE_dic,outfile_prex)
    elapsed = (time.time() - start)
    print("Time used:",elapsed)
