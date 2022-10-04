"""
Author: Zhang Chengsheng, @2020.08.04
不用pandas numpy还是不舒服，所以这个就用了
"""

import os,sys
import pandas as pd
import numpy as np
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from scipy import interpolate

def tableRead(file):
    return pd.read_table(file)


def statMerge(dfStat,dfRead):
    dfStatTrim = dfStat[['ID','Classification','Subtype','Seq_length','Seq_exon_num','exon_miss_to_transcript_start','exon_miss_to_transcript_end']]
    dfnew = pd.merge(left=dfStatTrim,right=dfRead,left_on='ID',right_on='ReadID',how='outer')
    dfnew.drop('ReadID',axis=1,inplace=True)
    return dfnew
    #dfnew.to_csv(r'D:\zcs-genex\SCRIPTS\ISOzcs\TransAnnot\report\t.txt',sep='\t',index=False)


def dfLineCount(df,key,target):
    return df.groupby(key,as_index=False)[target].count()


def dfLineSum(df,key,target):
    return df.groupby(key, as_index=False)[target].sum()


class DF2filter:
    def __init__(self,df):
        self.df = df
        self.types = ['FSM', 'ISM', 'NIC', 'NNC', 'Genic', 'Intergenic', 'FUSION', 'UNKNOWN']

    def typeCount(self,value1):
        key1 = 'Classification'
        df = self.df[self.df[key1]==value1]
        reads = len(df)
        ts = len(dfLineCount(df, 'Transcript', key1))
        gene = len(dfLineCount(df, 'Gene', key1))
        return reads,ts,gene

    def allTypeCount(self):
        reads = len(self.df)
        ts = len(dfLineCount(self.df,'Transcript','Classification'))
        gene = len(dfLineCount(self.df, 'Gene', 'Classification'))
        return reads,ts,gene

    def sheetClassification(self):
        cla = {}
        for i in self.types:
            cla[i] = list(self.typeCount(i))
        #cla['Total'] = self.allTypeCount()
        df = pd.DataFrame(cla,index=['Isoform','Transcript','Gene'],columns=self.types)
        return df.T

    def lengthClassification(self):
        res = {}
        for i in self.types:
            df = self.df[self.df['Classification'] == i]
            res[i] = df['Seq_length'].values
        return res


class ReportPdf:
    def __init__(self):
        self.colors = {'FSM':'mediumblue',
                       'ISM':'lightskyblue',
                       'NIC':'green',
                       'NNC':'red',
                       'Genic':'darkorange',
                       'Intergenic':'tan',
                       'FUSION':'purple',
                       'UNKNOWN':'grey'}

    def page0(self,sampleName):
        fig = plt.figure(figsize=(8,8))
        fig.text(0.5, 0.6,'TransAnnot Report',fontsize=40,ha='center',va='bottom')
        fig.text(0.5, 0.3, 'Sample: {}'.format(sampleName), fontsize=24, ha='center', va='bottom')
        fig.text(0.5,0.1,'V0.0.2 test',fontsize=24,ha='center',va='bottom')
        #fig.show()
        return fig

    def page1(self,df,df2):
        """
        Total count & Structural Classification
        参数1: A.sheetClassification()
        参数2：A.allTypeCount()
        """
        fig = plt.figure(figsize=(8,8))
        #ax1
        ax = fig.add_axes([0.1,0.05,1,0.4])
        cols = ['lightblue']
        table = ax.table(cellText=df.values,colWidths=[0.17]*3,rowLabels=df.index.tolist(),rowColours=cols*8,colLabels=df.columns.tolist(),colColours=cols*3,loc='bottom',bbox=[0.2,0.2,0.55,0.7])
        table.auto_set_font_size(False)
        table.set_fontsize(14)
        table.scale(1,1)
        ax.axis('off')
        #ax2
        ax2 = fig.add_axes([0.1,0.35,1,0.3])
        ax2.text(0.4,0.4,'Structural Classification',fontsize=20,verticalalignment='center',horizontalalignment='center')
        ax2.axis('off')
        #ax3
        ax3 = fig.add_axes([0.1,0.5,1,0.3])
        table2 = ax3.table(cellText=[[df2[0]],[df2[1]],[df2[2]]],colWidths=[1],rowLabels=['Unique Isoforms:','Unique Transcripts:','Unique Genes:'],loc='bottom',bbox=[0.5,0.6,0.5,0.5],edges='open')
        table2.auto_set_font_size(False)
        table2.set_fontsize(18)
        table2.scale(5, 5)
        ax3.axis('off')
        #fig.show()
        return fig

    def page2(self,df):
        """
        Isoform length
        参数：statMerge(dfStat,dfRead)['Seq_length']
        """
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_axes([0.12, 0.1, 0.8, 0.8])
        num = min(len(df),80)
        c,bins,d = ax.hist(df.values,bins=num,alpha=0.7)
        bins = (bins[1:]+bins[:-1])/2
        xnew = np.linspace(start=bins[0],stop=bins[-1],num=min(50,10*num))
        func = interpolate.interp1d(bins,c,kind='cubic')
        ynew = func(xnew)
        ax.plot(xnew, ynew)
        ax.set_xlabel('Isoform Length',fontsize=14)
        ax.set_ylabel('Isoform Count',fontsize=15)
        ax.xaxis.set_tick_params(labelsize=13)
        ax.yaxis.set_tick_params(labelsize=12)
        ax.set_title('Distribution of Isoform Lengths',fontsize=20)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        #fig.show()
        return fig

    def page3(self,df):
        """
        Isofrom length in classification
        参数：A.lengthClassification()
        """
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        for i in df:
            num = min(len(df[i]), 80)
            c, bins, d = ax.hist(df[i], bins=num, alpha=0)
            bins = (bins[1:] + bins[:-1]) / 2
            xnew = np.linspace(start=bins[0], stop=bins[-1], num=min(300, 10 * num))
            func = interpolate.interp1d(bins, c, kind=3)
            ynew = func(xnew)
            ax.plot(xnew, ynew,color=self.colors[i])
        legend = ax.legend(list(df),fontsize=12,markerscale=10,loc=7)
        for i in legend.get_lines():
            i.set_linewidth(3)
        ax.set_xlabel('Isoform Length', fontsize=14)
        ax.set_ylabel('Isoform Count', fontsize=15)
        ax.xaxis.set_tick_params(labelsize=13)
        ax.yaxis.set_tick_params(labelsize=12)
        ax.set_title('Distribution of Isoform Lengths', fontsize=20)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        #fig.show()
        return fig

    def page4(self,df):
        """
        Isoform Classification barplot
        参数: A.sheetClassification()
        """
        text = ['{:.2f}%'.format(i) for i in (df['Isoform']/df['Isoform'].sum()*100).values]
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_axes([0.12, 0.12, 0.8, 0.8])
        ax.bar(x=list(self.colors),height=df['Isoform'],color=self.colors.values(),edgecolor='black')
        for x,y,s in zip(list(self.colors),df['Isoform'],text):
            ax.text(x,y,s,fontsize=13,ha='center',va='bottom')

        ax.set_ylabel('Isoforms Count', fontsize=15)
        ax.xaxis.set_tick_params(labelsize=13,rotation=45)
        ax.yaxis.set_tick_params(labelsize=12)
        ax.set_title('Classification Distribution of Isoforms', fontsize=20)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        #fig.show()
        return fig

    def page5(self,df):
        """

        """
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_axes([0.12, 0.12, 0.8, 0.8])
        #ax.boxplot()


def test():
    sampleName = '759133C'
    stat = r'D:\zcs-genex\SCRIPTS\ISOzcs\TransAnnot\report\759133C.annot.stat'
    junctionStat = r'D:\zcs-genex\SCRIPTS\ISOzcs\TransAnnot\report\759133C.annot.junction'
    tsClu = r'D:\zcs-genex\SCRIPTS\ISOzcs\TransAnnot\report\759133C.annot.cluster.transcript'
    geneClu = r'D:\zcs-genex\SCRIPTS\ISOzcs\TransAnnot\report\759133C.annot.cluster.gene'
    readClu = r'D:\zcs-genex\SCRIPTS\ISOzcs\TransAnnot\report\759133C.annot.cluster.reads'
    outputPDF = r'D:\zcs-genex\SCRIPTS\ISOzcs\TransAnnot\report\t1.pdf'
    dfStat = tableRead(stat)
    dfRead = tableRead(readClu)
    df1 = statMerge(dfStat,dfRead)
    A = DF2filter(df1)
    df2 = A.sheetClassification()  # classification count
    df3 = A.allTypeCount()  # total count
    df4 = A.lengthClassification()  # 各分类的reads长度列表

    B = ReportPdf()
    f0 = B.page0(sampleName)
    f1 = B.page1(df2,df3)
    #B.page2(df1['Seq_exon_num'])
    f2 = B.page2(df1['Seq_length'])
    f3 = B.page3(df4)
    f4 = B.page4(df2)

    with PdfPages(outputPDF) as pdf:
        pdf.savefig(f0)
        pdf.savefig(f1)
        pdf.savefig(f2)
        pdf.savefig(f3)
        pdf.savefig(f4)




if __name__ == '__main__':
    test()