"""
Author: Zhang Chengsheng, @2020.06.15
"""

__version__ = '0.0.1'
from tkinter import *
import pickle
from matplotlib import pyplot as plt
import matplotlib.lines as lines
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


class TAGUI:
    def __init__(self,db):
        self.db = db['gene2transcript']
        self.fusion = db['fusion']
        self.info = db['reads']
        self.boxplot = boxplot2(db)
        self.UI = Tk()
        self._initWindow()

    def _initBox(self,row,col,height=15,width=15,rowspan=1,columnspan=1,func=Listbox,mode=BROWSE):
        Frame = LabelFrame(self.UI, text='')
        Frame.grid(row=row, column=col, rowspan=rowspan,columnspan=columnspan,padx=(20, 0), pady=(0, 10), sticky=W + N + E + S)
        ScrollX = Scrollbar(Frame, width=15, orient=HORIZONTAL)
        ScrollX.pack(side=BOTTOM, fill=X)
        ScrollY = Scrollbar(Frame, width=15, orient=VERTICAL)
        ScrollY.pack(side=RIGHT, fill=Y)
        if func == Listbox:
            listbox = func(Frame, width=width, height=height, selectmode=mode,exportselection=False,yscrollcommand=ScrollY.set, xscrollcommand=ScrollX.set)
        else:
            listbox = func(Frame,width=width, height=height,yscrollcommand=ScrollY.set, xscrollcommand=ScrollX.set)
        ScrollX.config(command=listbox.xview)
        ScrollY.config(command=listbox.yview)
        listbox.pack(fill=BOTH)
        return Frame,listbox

    def _initWindow(self):
        rowSearchLable,rowSearch,rowLab1,rowBox1,rowLab2,rowBox2,rowShowFigure,rowInfoLable,rowInfo=0,1,2,3,4,5,6,7,8
        col1,col2,col3,col4,col5=[0,2,4,6,8]
        self.UI.title("TransAnnot Viewer V0.0.1")
        self.UI.geometry('858x681+200+200')

        self.sampleLable = Label(self.UI, text="Sample")
        self.sampleLable.grid(row=rowLab1, column=col1)
        self.chromLable = Label(self.UI, text="Chromosome")
        self.chromLable.grid(row=rowLab1, column=col2)
        self.geneLable = Label(self.UI,text="Gene Level")
        self.geneLable.grid(row=rowLab1,column=col3)
        self.transcriptLable = Label(self.UI,text="Transcript Level")
        self.transcriptLable.grid(row=rowLab1,column=col4)
        self.fusionTranscriptLable = Label(self.UI, text="Fusion Transcript")
        self.fusionTranscriptLable.grid(row=rowLab2, column=col4)
        self.readLable = Label(self.UI,text="Read Level")
        self.readLable.grid(row=rowLab1,column=col5)
        self.fusionreadLable = Label(self.UI, text="Fusion Read")
        self.fusionreadLable.grid(row=rowLab2, column=col5)
        self.infoReadLable = Label(self.UI, text="Read infomation")
        self.infoReadLable.grid(row=rowInfoLable,column=col1)
        self.infoRefLable = Label(self.UI, text="Ref information")
        self.infoRefLable.grid(row=rowInfoLable,column=col4)

        # box
        checkInput = self.UI.register(self._inputCheck)
        self.searchLable = Label(self.UI, text="Gene Keyword Search")
        self.searchLable.grid(row=rowSearchLable, column=col1,padx=(20, 0))
        self.searchBox = Entry(self.UI, width=65, textvariable=StringVar(), validate='key',validatecommand=(checkInput,'%P','%v','%W'))
        self.searchBox.grid(row=rowSearch, column=col1, columnspan=15, padx=(20, 0), pady=(2, 2), sticky=W + N + S + E)
        self.genes = {}
        self.samples = {}
        self.f2g = {}
        self.g2f = {}
        self.sampleFrame,self.sampleList = self._initBox(rowBox1,col1,height=18,width=15,rowspan=3,mode=EXTENDED)
        self.chromFrame, self.chromList = self._initBox(rowBox1, col2, height=18, width=10, rowspan=3,mode=EXTENDED)
        self.geneFrame, self.geneList = self._initBox(rowBox1, col3, height=18, width=15, rowspan=3)
        self.transcriptFrame, self.transcriptList = self._initBox(rowBox1, col4, height=10, width=18,mode=EXTENDED)
        self.fTranscriptFrame, self.fTranscriptList = self._initBox(rowBox2, col4, height=5, width=18)
        self.readFrame, self.readList = self._initBox(rowBox1, col5, height=10, width=30,mode=EXTENDED)
        self.fReadFrame, self.fReadList = self._initBox(rowBox2, col5, height=5, width=30,mode=EXTENDED)
        self.readInfoFrame, self.readInofList = self._initBox(rowInfo, col1, height=12, width=50, columnspan=6,func=Text)
        self.refInfoFrame, self.refInofList = self._initBox(rowInfo, col4, height=12, width=50, columnspan=6, func=Text)
        self.readInofList.config(stat=DISABLED)  # NORMAL
        self.refInofList.config(stat=DISABLED)  # NORMAL
        for c in sorted(self.db):
            self.chromList.insert(END,c)
            for g in self.db[c]:
                self.geneList.insert(END,g)
                self.genes[g] = [c,[]]
                for t in self.db[c][g]:
                    for s in self.db[c][g][t]:
                        if s not in self.genes[g][1]:
                            self.genes[g][1].append(s)
                        if s not in self.samples:
                            self.samples[s] = {"gene":{g:0},"chrom":{c:0}}
                            self.sampleList.insert(END,s)
                        else:
                            self.samples[s]["gene"][g] = 0
                            self.samples[s]["chrom"][c] = 0
        self.chromList.insert(END, 'FUSION')
        for g in self.fusion:
            gs = g.split('|')
            for t in self.fusion[g]:
                for i in gs:
                    if i not in self.g2f:
                        self.g2f[i] = {}
                    self.g2f[i][t] = self.fusion[g][t]
                if g not in self.f2g:
                    self.f2g[t] = g

        ## 列表绑定函数
        self.geneList.bind('<<ListboxSelect>>',self._clickGeneList)
        self.transcriptList.bind('<<ListboxSelect>>', self._clickTranscriptList)
        self.chromList.bind('<<ListboxSelect>>', self._clickChromList)
        self.sampleList.bind('<<ListboxSelect>>',self._clickSample)
        self.fTranscriptList.bind('<<ListboxSelect>>',self._clickfTranscriptList)
        self.fReadList.bind('<<ListboxSelect>>',self._clickfReadList)
        self.readList.bind('<<ListboxSelect>>',self._clickReadList)
        ## 按钮
        self.drawButton = Button(self.UI,text='Show Figure',command=self._draw)
        self.drawButton.grid(row=rowShowFigure,column=col1,columnspan=15,padx=(20,0),sticky=EW)

        #self.searchBox.bind('<key-input>',print(self.KEYWORD))

    def _inputCheck(self,content,reason,name):
        self._refreshGenebox(content)
        #print(content,reason,name)
        return True

    def _refreshGenebox(self,keyword):
        self.geneList.delete(0,END)
        #self.chromList.select_clear(0, END)
        #self.sampleList.select_clear(0, END)
        samples = [self.sampleList.get(i) for i in self.sampleList.curselection()]
        samples = samples if samples else self.sampleList.get(0,END)
        chroms = [self.chromList.get(i) for i in self.chromList.curselection()]
        chroms = chroms if chroms else self.chromList.get(0,END)
        genes = []

        if keyword == '':
            for sample in samples:
                for g in self.samples[sample]["gene"]:
                    if self.genes[g][0] in chroms:
                        if g not in genes:
                            genes.append(g)
                            self.geneList.insert(END, g)
        else:
            for sample in samples:
                for g in self.samples[sample]["gene"]:
                    if keyword.upper() in g.upper() and self.genes[g][0] in chroms:
                        if g not in genes:
                            genes.append(g)
                            self.geneList.insert(END,g)

    def _clickGeneList(self,event):
        self.transcriptList.delete(0,END)
        self.readList.delete(0, END)
        self.fReadList.delete(0,END)
        self.fTranscriptList.delete(0,END)
        g = self.geneList.get(self.geneList.curselection())
        samples = [self.sampleList.get(i) for i in self.sampleList.curselection()]
        samples = samples if samples else self.samples
        if g in self.genes:
            c,ss = self.genes[g]
            for t in self.db[c][g]:
                for s in samples:
                    if s in self.db[c][g][t]:
                        self.transcriptList.insert(END, t)
                        break
                for s in self.db[c][g][t]:
                    if s in samples:
                        for r in self.db[c][g][t][s]:
                            self.readList.insert(END, '{}    {}'.format(r, s))
        else:
            for t in self.fusion[g]:
                for s in samples:
                    if s in self.fusion[g][t]:
                        self.fTranscriptList.insert(END,t)
                        break
                for s in self.fusion[g][t]:
                    if s in samples:
                        for r in self.fusion[g][t][s]:
                            self.fReadList.insert(END,'{}    {}'.format(r,s))
        if g in self.g2f:
            for t in self.g2f[g]:
                for s in samples:
                    if s in self.g2f[g][t]:
                        self.fTranscriptList.insert(END, t)
                        break
                for s in self.g2f[g][t]:
                    if s in samples:
                        for r in self.g2f[g][t][s]:
                            self.fReadList.insert(END,'{}    {}'.format(r,s))

    def _clickReadList(self,event):
        self.fTranscriptList.select_clear(0,END)
        self.fReadList.select_clear(0,END)

        reads = [self.readList.get(i) for i in self.readList.curselection()]
        if len(reads) == 1:
            self.readInofList.config(stat=NORMAL)
            self.refInofList.config(stat=NORMAL)
            self.readInofList.delete(1.0, END)
            self.refInofList.delete(1.0, END)
            r,s = reads[0].split('    ')
            info = self.info[s][r]
            txt = ''
            for idx in range(len(info['chrom'])):
                txt += '{}\t{}\t{}\t{}\t{}\n'.format(info['chrom'][idx],info['strand'][idx],info['start'][idx],info['end'][idx],info['anno'][idx])
            self.readInofList.insert(END,txt)
            self.refInofList.config(stat=DISABLED)
            self.readInofList.config(stat=DISABLED)

    def _clickfReadList(self,event):
        self.transcriptList.select_clear(0,END)
        self.readList.select_clear(0,END)

    def _clickfTranscriptList(self,event):
        self.fReadList.delete(0,END)
        self.transcriptList.select_clear(0,END)
        self.readList.select_clear(0,END)
        samples = [self.sampleList.get(i) for i in self.sampleList.curselection()]
        samples = samples if samples else self.samples
        t = self.fTranscriptList.get(self.fTranscriptList.curselection())
        g = self.f2g[t]
        for s in self.fusion[g][t]:
            if s in samples:
                for r in self.fusion[g][t][s]:
                    self.fReadList.insert(END,'{}    {}'.format(r,s))

    def _clickTranscriptList(self,event):
        self.fTranscriptList.select_clear(0, END)
        self.fReadList.select_clear(0, END)
        self.readList.delete(0,END)
        g = self.geneList.get(self.geneList.curselection())
        c, ss = self.genes[g]
        samples = [self.sampleList.get(i) for i in self.sampleList.curselection()]
        samples = samples if samples else self.samples
        ts = [self.transcriptList.get(i) for i in self.transcriptList.curselection()]
        ts = ts if ts else self.transcriptList.get(0,END)
        for t in ts:
            for s in self.db[c][g][t]:
                if s in samples:
                    for r in self.db[c][g][t][s]:
                        self.readList.insert(END,'{}    {}'.format(r,s))

    def _clickChromList(self,event):
        self.geneList.delete(0,END)
        self.transcriptList.delete(0,END)
        self.readList.delete(0,END)
        samples = [self.sampleList.get(i) for i in self.sampleList.curselection()]
        genes = []
        if self.chromList.curselection():
            for i in self.chromList.curselection():
                c = self.chromList.get(i)
                if c == 'FUSION':
                    for fu in self.fusion:
                        self.geneList.insert(END,fu)
                    continue
                samples = samples if samples else self.samples
                for s in samples:
                    for g in self.samples[s]['gene']:
                        if g not in genes and self.genes[g][0] == c:
                            genes.append(g)
                            self.geneList.insert(END, g)
        else:
            for g in self.genes:
                if g not in genes:
                    genes.append(g)
                    self.geneList.insert(END, g)

    def _clickSample(self,event):
        self.transcriptList.delete(0, END)
        self.readList.delete(0, END)
        self.chromList.delete(0, END)
        self.geneList.delete(0, END)
        if self.sampleList.curselection():
            for i in self.sampleList.curselection():
                s = self.sampleList.get(i)
                for c in self.samples[s]['chrom']:
                    self.chromList.insert(END, c)
                for g in self.samples[s]['gene']:
                    self.geneList.insert(END, g)
        else:
            cs,gs = [],[]
            for s in self.samples:
                for c in self.samples[s]['chrom']:
                    if c not in cs:
                        cs.append(c)
                        self.chromList.insert(END, c)
                for g in self.samples[s]['gene']:
                    if g not in gs:
                        gs.append(g)
                        self.geneList.insert(END, g)
        if self.fusion:
            self.chromList.insert(END, 'FUSION')

    def _draw(self):
        g = self.geneList.get(self.geneList.curselection())
        c = self.genes[g][0]
        #t = [self.transcriptList.get(i) for i in self.transcriptList.curselection()]
        s = [self.sampleList.get(i) for i in self.sampleList.curselection()]
        r = [self.readList.get(i) for i in self.readList.curselection()]
        r = r if r else self.readList.get(0,END)
        plt = self.boxplot.draw_gene(s,g,r)
        x,y = plt.get_size_inches()
        draw = Toplevel()
        draw.geometry('1000x500+10+10')
        #draw.geometry('{}x{}+10+10'.format(int(x*100),int(y*100)))
        title = g
        draw.title(title)
        p = FigureCanvasTkAgg(plt,draw)
        p.draw()
        p.get_tk_widget().pack(side=TOP,fill=BOTH,expand=1)

    def show(self):
        self.UI.mainloop()


class boxplot2:
    def __init__(self,db):
        self.db = db

    def draw_gene(self,sample,gene,read):
        reads, anno, g2t = self.db['reads'], self.db['annotation'], self.db['gene2transcript']
        rs = [i.split('    ') for i in read]
        test_dict = {}
        for r,s in rs:
            test_dict[r] = reads[s][r]
        if not test_dict:
            print("?????????????????????")
        anno_db = anno[gene] if gene in anno else {}
        fig = self.normal_boxplot(test_dict, anno_db=anno_db, save=0, TAUI=1)
        return fig

    def normal_boxplot(self,draw_db, anno_db, save=0, fusion=0, TAUI=0):
        y_num = len(draw_db) + len(anno_db)
        # path = outputdir
        width = 0.5
        colTable = {0: 'gray', 10: '#3B3D42', 1: '#1975C4', 2: '#E68E29', 3: '#199F15', 4: '#9F3241', 5: '#18D6CB',
                    6: '#D670D2', 7: '#8BB6D6', 8: '#C6D611', 9: '#9425D6'}
        col_id_list, col_anno_list = [], []
        col_idx = 0
        colDict = {}
        all_lr = []
        for id in draw_db:
            col_id_list.append(draw_db[id]['info'][1])
            all_lr += draw_db[id]['start']
            all_lr += draw_db[id]['end']
        left, right = min(all_lr), max(all_lr)
        for transcript in anno_db:
            col_anno_list.append(transcript)
            all_lr += anno_db[transcript]['start']
            all_lr += anno_db[transcript]['end']
        for i in col_anno_list:
            if i not in colDict:
                colDict[i] = 10
        for i in col_id_list:
            if i in colDict:
                if colDict[i] == 10:
                    col_idx += 1
                    colDict[i] = min(col_idx, 10)
            else:
                colDict[i] = 0
        r = max(all_lr) - min(all_lr)
        multi = min(int(r / 150000) + 1, 6)
        fig = plt.figure(figsize=(multi * 10, (1 + (multi - 1) * 0.5) * 5 * (int(y_num / 20) + 1)))
        ax = fig.add_axes([0.2, 0.1, 0.75, 0.85])
        y_ticks = [idx for idx, i in enumerate(draw_db)]
        y_tickslabels = [i.split('_')[-1] for i in draw_db]
        for idx1, id in enumerate(draw_db):
            this_lr = []
            this_lr += draw_db[id]['start']
            this_lr += draw_db[id]['end']
            colBase = colTable[colDict[draw_db[id]['info'][1]]]
            for idx2, i in enumerate(draw_db[id]['start']):
                col = colBase
                if not idx2 or idx2 == len(draw_db[id]['anno']) - 1:
                    if draw_db[id]['anno'][idx2] not in ['左边缘外显子', '右边缘外显子']:
                        col = 'red'
                else:  # idx2 and idx2 != len(draw_db[id]['anno']) -1:
                    if fusion:
                        if draw_db[id]['anno'][idx2] not in ['', '左边缘外显子', '右边缘外显子']:
                            col = 'red'
                    elif draw_db[id]['anno'][idx2]:
                        col = 'red'
                length = abs(int(draw_db[id]['start'][idx2]) - int(draw_db[id]['end'][idx2]))
                rect = plt.Rectangle((int(i), idx1 - width / 2), length, width, color=col, joinstyle='miter',
                                     capstyle='butt', alpha=1)
                ax.add_patch(rect)
            ax.add_line(lines.Line2D((min(this_lr), max(this_lr)), (idx1, idx1), linewidth=0.5, solid_capstyle='butt',
                                     solid_joinstyle='miter', color='black', alpha=0.5))
        for idx1, transcript in enumerate(sorted(anno_db)):
            start = anno_db[transcript]['start']
            end = anno_db[transcript]['end']
            colBase = 'black' if fusion else colTable[colDict[transcript]]
            r = start + end
            s = idx1 + len(draw_db) + 1
            y_ticks.append(s)
            y_tickslabels.append(transcript)
            ax.add_line(
                lines.Line2D((min(r), max(r)), (s, s), linewidth=0.5, solid_capstyle='butt', solid_joinstyle='miter',
                             color='black', antialiased=False, alpha=0.5))
            for idx2, k in enumerate(start):
                length = abs(start[idx2] - end[idx2])
                rect = plt.Rectangle((min(start[idx2], end[idx2]), s - width / 2), length, width, color=colBase,
                                     joinstyle='miter', capstyle='butt', alpha=1)
                ax.add_patch(rect)

        ax.set_xticks(
            [min(all_lr), int(0.75 * min(all_lr) + 0.25 * max(all_lr)), int(0.5 * (min(all_lr) + max(all_lr))),
             int(0.25 * min(all_lr) + 0.75 * max(all_lr)), max(all_lr)])
        ax.set_xticklabels(
            [min(all_lr), int(0.75 * min(all_lr) + 0.25 * max(all_lr)), int(0.5 * (min(all_lr) + max(all_lr))),
             int(0.25 * min(all_lr) + 0.75 * max(all_lr)), max(all_lr)])
        ax.set_xlim([min(all_lr), max(all_lr)])
        ax.set_yticks(y_ticks)
        ax.set_yticklabels(y_tickslabels)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(axis=u'y', which=u'both', length=0)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize((1 + (multi - 1) * 0.5) * 10)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize((1 + (multi - 1) * 0.5) * 10)
        ax.plot()
        # fig.show()
        if save:
            try:
                fig.savefig(save, bbox_inches='tight')
            except:
                pass
        else:
            if TAUI:
                return fig
            else:
                fig.show()


def load_pickle(file_in):
    return pickle.load(open(file_in, 'rb'))


def showGUI():
    #dbfile = r'D:\zcs-genex\SCRIPTS\ISOzcs\TransAnnot\TAUI\1\776415C.annot.db.pickle'
    #dbfile = r'D:\zcs-genex\SCRIPTS\ISOzcs\TransAnnot\TAUI\1\merge.db.pickle'
    dbfile = r'D:\zcs-genex\SCRIPTS\ISOzcs\workflow\20200730\merge.db.pickle'
    #dbfile = r'D:\zcs-genex\SCRIPTS\ISOzcs\TransAnnot\report\759133C.annot.db.pickle'
    db = load_pickle(dbfile)
    A = TAGUI(db)
    A.show()


if __name__ == '__main__':
    showGUI()
