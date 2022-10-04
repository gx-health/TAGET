"""

"""


import os,sys
from matplotlib import pyplot as plt
import matplotlib.lines as lines
import pickle


def load_pickle(file_in):
    return pickle.load(open(file_in, 'rb'))


def check_overlap(ref, gene):
    flag = 0
    if ref[0] <= gene[0]:
        if ref[1] <= gene[0]:
            flag = 0  # 不挨着ref在左
        elif ref[1] >= gene[1]:
            flag = 2  # ref 包含gene
        else:
            flag = 1  # overlap
    else:
        if ref[0] >= gene[1]:
            flag = 0  # 不挨着ref在右
        elif ref[1] <= gene[1]:
            flag = 2  # gene包含ref
        else:
            flag = 1  # overlap
    return flag


def get_gene_structure(chr,region):
    genes = []
    chr = chr.lstrip('chr')
    if chr == 'MT':
        return genes
    for d in [d1,d2,d3]:
        for gene in d[chr]:
            gene_start, gene_end = d[chr][gene][2]
            overlap = check_overlap([int(region[0]), int(region[1])], [gene_start, gene_end])
            if overlap:
                genes.append(d[chr][gene])
    return genes


def check_gene_split(list1):
    def get_split(type=1):
        if type == 1:
            res_all = []
            lock = 0
            for idx,[ref1,ref2,chr1,strand,anno] in enumerate(list1):
                if not idx:
                    chr_temp = chr1
                    res = [[ref1,ref2,chr1,strand,anno]]
                    continue
                if chr1 != chr_temp:
                    if lock:
                        return 0
                    res_all.append(res)
                    res = [[ref1,ref2,chr1,strand,anno]]
                    lock = 1
                    chr_temp = chr1
                else:
                    res.append([ref1,ref2,chr1,strand,anno])
            else:
                res_all.append(res)

        elif type == 2:
            res_all = []
            lock = 0
            for idx, [ref1, ref2, chr1, strand, anno] in enumerate(list1):
                if not idx:
                    site_temp = int(ref1),int(ref2)
                    res = [[ref1, ref2, chr1, strand, anno]]
                    continue
                if abs(int(ref1) - site_temp[0]) > 500000:
                    if lock:
                        return 0
                    res_all.append(res)
                    res = [[ref1, ref2, chr1, strand, anno]]
                    lock = 1
                    site_temp = int(ref1),int(ref2)
                else:
                    res.append([ref1, ref2, chr1, strand, anno])
            else:
                res_all.append(res)

        else:
            res_all = 0
        return res_all

    chr_list,sites = [],[]
    for ref1,ref2,chr1,strand,anno in list1:
        chr_list.append(str(chr1))
        sites.append(int(ref1))
        sites.append(int(ref2))
    if len(set(chr_list)) > 1:
        res = get_split(type=1)
        if not res:
            return -1
        return res
    else:
        sites = sorted(sites)
        for idx,i in enumerate(sites):
            if not idx:
                continue
            if i - sites[idx-1] > 500000:
                res = get_split(type=2)
                if not res:
                    return -1
                return res
    return 0


def bed_read(file_in,id_specific=0):
    transcripts = {}
    genes = {}
    with open(file_in,'r') as f:
        for i in f.readlines():
            if not i.strip():
                continue
            chrom,start,end,transcript,p,strand,seq_start,seq_end = i.strip().split('\t')[0:8]
            try:
                anno = i.strip().split('\t')[8]
            except:
                anno = ''
            if chrom == 0 or chrom == '0':
                continue
            if id_specific and transcript != id_specific:
                continue
            if transcript not in transcripts:
                transcripts[transcript] = [[start,end,chrom,strand,anno]]
            else:
                transcripts[transcript].append([start,end,chrom,strand,anno])
            if ':' in anno:
                gene = anno.split(':')[0]
                if gene in genes:
                    if transcript not in genes[gene]:
                        genes[gene].append(transcript)
                else:
                    genes[gene] = [transcript]
    return transcripts,genes


def fusion_inter_box(data,t_id,outputdir):
    def get_range(list1):
        r = []
        for ref1,ref2,chr1,strand,anno in list1:
            r.append(int(ref1))
            r.append(int(ref2))
        return chr1,[min(r),max(r)]

    def loc_trans(n):
        return n-min(all_lr2)+max(all_lr1)+space

    #data[0],data[1] = data[1],data[0]
    t_id = t_id.split('_')[-1]
    fig_name = t_id.replace('/','_') + '.pdf'
    fig_name = os.path.join(outputdir,fig_name)
    r1 = get_range(data[0])
    r2 = get_range(data[1])
    genes1 = get_gene_structure(*r1)
    if not genes1:
        return 0
    genes2 = get_gene_structure(*r2)

    y_num_1 = sum([len(i[3]) for i in genes1])
    y_num_2 = sum([len(i[3]) for i in genes2])
    y_num = max(y_num_1,y_num_2) + 1
    left1,right1 = r1[1]
    left2,right2 = r2[1]
    all_lr1 = r1[1]
    all_lr2 = r2[1]
    for idx1, i in enumerate(genes1):
        gene_type, strand, lr, transcripts = i
        all_lr1.append(lr[0])
        all_lr1.append(lr[1])
    rr1 = max(all_lr1)-min(all_lr1)
    for idx1, i in enumerate(genes2):
        gene_type, strand, lr, transcripts = i
        all_lr2.append(lr[0])
        all_lr2.append(lr[1])
    rr2 = max(all_lr2)-min(all_lr2)
    space = int((rr1+rr2)/10)
    space_s = (rr1+0.5*space)/(rr1+rr2+space)

    multi = min(int((rr1+rr2+space) / 150000) + 1, 6)
    fig = plt.figure(figsize=(multi * 20, (1 + (multi - 1) * 0.5) * 5 * (int(y_num / 20) + 1)))
    ax = fig.add_axes([0.25, 0.1, 0.65, 0.75])
    y_ticks1 = [0]
    y_tickslabels1 = [t_id]
    y_ticks2 = [-1,0]
    y_tickslabels2 = ['',t_id]

    width = 0.5

    ## fig1
    ax.add_line(lines.Line2D((left1,right1), (0, 0), linewidth=0.5, solid_capstyle='butt', solid_joinstyle='miter', color='black',antialiased=False, alpha=0.5))
    for idx,i in enumerate(data[0]):
        ref1, ref2, chr1, strand, anno = i
        length = int(ref2)-int(ref1)
        rect = plt.Rectangle((int(ref1), 0-width/2), length, width, color='r',joinstyle='miter',capstyle='butt', alpha=1)
        ax.add_patch(rect)
        ## add junction line
        if not idx:
            pass
        else:
            pass
    for idx1,i in enumerate(genes1):
        gene_type, strand, lr, transcripts = i
        for idx2,j in enumerate(transcripts):
            #print(j, transcripts[j])
            exons = transcripts[j][2]
            r = transcripts[j][1]
            s = sum([len(genes1[i][3]) for i in range(max(0,idx1))]) + idx2 + 1
            y_ticks1.append(s)
            y_tickslabels1.append(j)
            ax.add_line(lines.Line2D((r[0], r[1]), (s, s), linewidth=0.5, solid_capstyle='butt', solid_joinstyle='miter',color='black', antialiased=False,alpha=0.5))
            for k in exons:
                region = transcripts[j][2][k]
                length = abs(region[1]-region[0])
                rect = plt.Rectangle((min(region), s - width / 2), length, width, color='black',joinstyle='miter',capstyle='butt', alpha=1)
                ax.add_patch(rect)
    else:
        y_ticks1.append(max(y_ticks1)+1)
        y_tickslabels1.append('')

    # fig2
    ax.add_line(lines.Line2D((loc_trans(left2), loc_trans(right2)), (0, 0), linewidth=0.5, solid_capstyle='butt', solid_joinstyle='miter',
                             color='black', antialiased=False, alpha=0.5))
    for idx, i in enumerate(data[1]):
        ref1, ref2, chr1, strand, anno = i
        length = int(ref2) - int(ref1)
        rect = plt.Rectangle((int(loc_trans(int(ref1))), 0 - width / 2), length, width, color='b', joinstyle='miter', capstyle='butt',alpha=1)
        ax.add_patch(rect)
        ## add junction line
        if not idx:
            pass
        else:
            pass
    for idx1, i in enumerate(genes2):
        gene_type, strand, lr, transcripts = i
        for idx2, j in enumerate(transcripts):
            # print(j, transcripts[j])
            exons = transcripts[j][2]
            r = transcripts[j][1]
            s = sum([len(genes2[i][3]) for i in range(max(0, idx1))]) + idx2 + 1
            y_ticks2.append(s)
            y_tickslabels2.append(j)
            ax.add_line(lines.Line2D((loc_trans(r[0]), loc_trans(r[1])), (s, s), linewidth=0.5, solid_capstyle='butt', solid_joinstyle='miter',color='black', antialiased=False, alpha=0.5))
            for k in exons:
                region = transcripts[j][2][k]
                length = abs(region[1] - region[0])
                rect = plt.Rectangle((loc_trans(min(region)), s - width / 2), length, width, color='black', joinstyle='miter',capstyle='butt', alpha=1)
                ax.add_patch(rect)
    else:
        y_ticks2.append(max(y_ticks2)+1)
        y_tickslabels2.append('')

    xtickslabel = []
    xticks = []
    xticks1,xticks2 = [],[]
    xtickslabel1,xtickslabel2 = [],[]


    ticks_num1 = 2 if rr1*5 < rr2 else 5
    for i in range(ticks_num1):
        ratio = i*1/(ticks_num1-1)
        tick = min(all_lr1)*(1-ratio)+max(all_lr1)*ratio
        xticks.append(int(tick))
        xtickslabel.append(int(tick))

        xticks1.append(int(tick))
        xtickslabel1.append(int(tick))
    ticks_num2 = 2 if rr2 * 5 < rr1 else 5
    for i in range(ticks_num2):
        ratio = i * 1 / (ticks_num2 - 1)
        tick = min(all_lr2) * (1 - ratio) + max(all_lr2) * ratio
        xticks.append(loc_trans(int(tick)))
        xtickslabel.append(int(tick))

        xticks2.append(loc_trans(int(tick)))
        xtickslabel2.append(int(tick))

    ax.set_xticks(xticks)
    ax.set_xticklabels(xtickslabel,size=(1 + (multi - 1) * 0.5) * 10)
    ax.set_xlim(min(xticks),max(xticks))

    a = 0.02
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((space_s, space_s), (-a, a), **kwargs)

    ax.set_yticks(y_ticks1)
    ax.set_yticklabels(y_tickslabels1,size=(1 + (multi - 1) * 0.5) * 10)
    ax.set_ylim(-1,max(max(y_ticks1),max(y_ticks2)))
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis=u'y', which=u'both', length=0)

    ax2 = ax.twinx()
    ax2.set_yticks(y_ticks2)
    ax2.set_yticklabels(y_tickslabels2,size=(1 + (multi - 1) * 0.5) * 10)
    ax2.set_ylim(-1, max(max(y_ticks1), max(y_ticks2)))
    ax2.spines['top'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.tick_params(axis=u'y', which=u'both', length=0)
    #for tick in ax.xaxis.get_major_ticks():
    #    tick.label.set_fontsize((1 + (multi - 1) * 0.5) * 10)
    #for tick in ax.yaxis.get_major_ticks():
    #    tick.label.set_fontsize((1 + (multi - 1) * 0.5) * 10)
    #for tick in ax2.xaxis.get_major_ticks():
    #    tick.label.set_fontsize((1 + (multi - 1) * 0.5) * 10)
    #for tick in ax2.yaxis.get_major_ticks():
    #    tick.label.set_fontsize((1 + (multi - 1) * 0.5) * 10)
    ax.plot()
    fig.show()
    #fig.savefig(fig_name, bbox_inches='tight')


def normal_boxplot(genes,gene,tss,outputdir,save=1):
    y_num_1 = sum([len(i[3]) for i in genes])
    y_num_2 = len(tss)
    y_num = y_num_1 + y_num_2
    path = outputdir
    width = 0.5
    fig_name = os.path.join(path,gene.split('/')[-1]+'.pdf')
    all_lr = []
    for t in tss:
        for i in tss[t]:
            all_lr.append(int(i[0]))
            all_lr.append(int(i[1]))
    left,right = min(all_lr),max(all_lr)
    for idx1, i in enumerate(genes):
        gene_type, strand, lr, transcripts = i[0:4]
        all_lr.append(lr[0])
        all_lr.append(lr[1])
    r = max(all_lr)-min(all_lr)
    multi = min(int(r/150000)+1,6)
    fig = plt.figure(figsize=(multi*20,(1+(multi-1)*0.5)*5*(int(y_num/20)+1)))
    ax = fig.add_axes([0.3,0.1,0.65,0.8])
    y_ticks = [idx for idx,i in enumerate(tss)]
    y_tickslabels = [i.split('_')[-1] for i in tss]
    cols = ['red','blue']
    col_idx = 0
    col_tab = {}
    for idx1,t in enumerate(tss):
        this_lr = []
        for i in tss[t]:
            anno = i[4]
            this_lr.append(int(i[0]))
            this_lr.append(int(i[1]))
            if ':' in anno:
                gene_name = anno.split(':')[0]
                if gene_name not in col_tab:
                    col_tab[gene_name] = cols[col_idx]
                    col_idx = 1
                col = col_tab[gene_name]
            else:
                col = 'grey'
            length = abs(int(i[1])-int(i[0]))
            rect = plt.Rectangle((int(i[0]), idx1 - width / 2), length, width, color=col, joinstyle='miter',capstyle='butt', alpha=1)
            ax.add_patch(rect)
        ax.add_line(lines.Line2D((min(this_lr),max(this_lr)), (idx1, idx1), linewidth=0.5, solid_capstyle='butt', solid_joinstyle='miter', color='black',alpha=0.5))
    for idx1,i in enumerate(genes):
        gene_type, strand, lr, transcripts = i[0:4]
        for idx2,j in enumerate(sorted(transcripts,key=lambda x:x)):
            if idx2 > 20:
                continue
            exons = transcripts[j][2]
            r = transcripts[j][1]
            s = sum([len(genes[i][3]) for i in range(max(0,idx1))]) + idx2 + y_num_2
            y_ticks.append(s)
            y_tickslabels.append(j)
            ax.add_line(lines.Line2D((r[0], r[1]), (s, s), linewidth=0.5, solid_capstyle='butt', solid_joinstyle='miter',color='black', antialiased=False,alpha=0.5))
            for k in exons:
                region = transcripts[j][2][k]
                length = abs(region[1] - region[0])
                rect = plt.Rectangle((min(region), s - width / 2), length, width, color='black', joinstyle='miter',capstyle='butt', alpha=1)
                ax.add_patch(rect)

    ax.set_xticks([min(all_lr),int(0.75*min(all_lr)+0.25*max(all_lr)),int(0.5*(min(all_lr)+max(all_lr))),int(0.25*min(all_lr)+0.75*max(all_lr)),max(all_lr)])
    ax.set_xticklabels([min(all_lr),int(0.75*min(all_lr)+0.25*max(all_lr)),int(0.5*(min(all_lr)+max(all_lr))),int(0.25*min(all_lr)+0.75*max(all_lr)),max(all_lr)])
    ax.set_xlim([min(all_lr),max(all_lr)])
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_tickslabels)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis=u'y', which=u'both', length=0)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize((1+(multi-1)*0.5)*10)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize((1+(multi-1)*0.5)*10)
    ax.plot()
    #fig.show()
    if save:
        fig.savefig(fig_name,bbox_inches='tight')
    else:
        fig.show()


def main():
    file_in = r'D:\zcs-genex\SCRIPTS\ISOzcs\script_v191115\anno\draw_20200117\759133C.drawfusion.ensembl.bed'
    outputdir = r'D:\zcs-genex\SCRIPTS\ISOzcs\script_v191115\anno\draw_20200117\759133C_all'
    transcripts,genes = bed_read(file_in)
    t_drop = []
    #genes_i_need = ['TYMP','ODF3B']
    #genes_i_need = ['ASB16', 'RPS6KL1','BCL2L12']
    #genes_i_need = ['MAG']
    for i in transcripts:
        a = check_gene_split(transcripts[i])
        if a:
            t_drop.append(i)
            print('fusion:',i)
            if a != -1:
                fusion_inter_box(a,i,outputdir)
    for gene in genes:
        #if gene not in genes_i_need:
            #continue
        print(gene)
        j = []
        ts = genes[gene]
        tss = {}
        for t in ts:
            if t in t_drop:
                continue
            tss[t]= transcripts[t]
            for j1, j2, chr1, strand, anno in transcripts[t]:
                j.append(j1)
                j.append(j2)
        else:
            if not tss:
                continue
        region = [min(j), max(j)]
        if int(region[1]) - int(region[0]) > 1000000:
            continue
        genes_anno = get_gene_structure(chr1, region)
        normal_boxplot(genes_anno, gene, tss, outputdir)


def iso_anno(file_in,db):
    def get_gene_structure_iso(chr, region,db):
        genes = []
        chr = chr.lstrip('chr')
        if chr == 'MT':
            return genes
        for d in [db]:
            for gene in d[chr]:
                gene_start, gene_end = d[chr][gene][2]
                overlap = check_overlap([int(region[0]), int(region[1])], [gene_start, gene_end])
                if overlap:
                    genes.append(d[chr][gene])
        return genes
    outputdir = ''
    transcripts, genes = bed_read(file_in)
    for i in transcripts:
        for j in transcripts[i]:
            print(j)
        break
    for i in transcripts:
        a = check_gene_split(transcripts[i])
        if a:
            #t_drop.append(i)
            print('fusion:', i)
            if a != -1:
                fusion_inter_box(a, i, outputdir)

    for gene in genes:
        j = []
        ts = genes[gene]
        tss = {}
        for t in ts:
            tss[t] = transcripts[t]
            for j1, j2, chr1, strand, anno in transcripts[t]:
                j.append(j1)
                j.append(j2)
        else:
            if not tss:
                continue
        region = [min(j), max(j)]
        if int(region[1]) - int(region[0]) > 1000000:
            continue
        genes_anno = get_gene_structure_iso(chr1, region,db)
        normal_boxplot(genes_anno, gene, tss, outputdir,save=0)


if __name__ == '__main__':
    db_ensembl = r'D:\zcs-genex\SCRIPTS\SVzcs\split_mapping\dev0725\over_1_0829\test191224\alldata\hg38.ensembl.gtf.pickle'
    db_ncbi = r'D:\zcs-genex\SCRIPTS\SVzcs\split_mapping\dev0725\over_1_0829\test191224\alldata\hg38.ncbi.gtf.pickle'
    db_gencode = r'D:\zcs-genex\SCRIPTS\SVzcs\split_mapping\dev0725\over_1_0829\test191224\alldata\hg38.gencode.gtf.pickle'
    d1 = load_pickle(db_ensembl)
    d2 = load_pickle(db_ncbi)
    d3 = load_pickle(db_gencode)
    main()
