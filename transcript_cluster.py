"""
Author: Zhang Chengsheng, @2020.03.20
"""

import os,sys


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


def ss_compare(ref,target):
    #ref = [[10000,10050],[10100,10200],[10500,10800],[10900,12000],[12500,15000],[15100,15800],[16000,16800],[19850,20000],[200100,25000]]
    #target = [[10600,10800],[10900,12000],[12500,15000],[15100,15800],[16000,16800],[19850,20000]]
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


class structure_save:
    def __init__(self):
        self.d = {}

    def add_id(self,id,structure):
        if id not in self.d:
            self.d[id] = structure

    def saveReadsinfo(self,A):
        id = A.id
        info = A.reads_info()
        self.d[id] = info

    def saveGene2transcript(self,A):
        gene = A.best_gene
        transcript = A.best_transcript
        classification = A.classification
        id = A.id
        if gene in self.d:
            if transcript in self.d[gene]:
                self.d[gene][transcript][id] = classification
            else:
                self.d[gene][transcript] = {id:classification}
        else:
            self.d[gene] = {transcript:{id:classification}}

    def addGene2transcript(self,db):
        for gene in db:
            if gene not in self.d:
                self.d[gene] = db[gene]
            else:
                for transcript in db[gene]:
                    if transcript in self.d[gene]:
                        for read in db[gene][transcript]:
                            if read not in self.d[gene][transcript]:
                                self.d[gene][transcript][read] = db[gene][transcript][read]
                    else:
                        self.d[gene][transcript] = db[gene][transcript]

    def saveAnnoDB(self,A):
        id = A.id
        anno = A.db_used
        for gene in anno:
            if gene not in self.d:
                self.d[gene] = anno[gene]

    def saveClassification(self,A):
        id = A.id
        classification = A.classification
        subtype = A.subtype
        if classification in self.d:
            self.d[classification][id] = subtype
        else:
            self.d[classification] = {id:subtype}


class tcluster:
    def __init__(self):
        #self.type = classification
        self.tree = {}
        self.NSM = {}
        self.FUSION = {}

    def init_NSM_node(self,region,ref,t_id,gene,classificaion):
        return {'node':[region,[t_id],{0:[region,ref,[t_id],[gene],[classificaion]]},[gene]],'up':0,'down':0}

    def add_NSM_node(self,chrom,region,exons,t_id,gene,classificaion,node=0):
        if not node:  # 初次检索
            if chrom in self.tree:
                if check_overlap(self.tree[chrom]['node'][0],region):
                    self.tree[chrom]['node'][1].append(t_id)
                    if gene not in self.tree[chrom]['node'][3]:
                        self.tree[chrom]['node'][3].append(gene)
                    for i in self.tree[chrom]['node'][2]:
                        if check_overlap(self.tree[chrom]['node'][2][i][0],region):
                            same, new_ref, new_range = ss_compare(self.tree[chrom]['node'][2][i][1],exons)
                            if same:
                                self.tree[chrom]['node'][2][i][2].append(t_id)
                                self.tree[chrom]['node'][2][i][1] = new_ref
                                self.tree[chrom]['node'][2][i][0] = [min(self.tree[chrom]['node'][2][i][0]+region),max(self.tree[chrom]['node'][2][i][0]+region)]
                                self.tree[chrom]['node'][2][i][4].append(classificaion)
                                if gene not in self.tree[chrom]['node'][2][i][3]:
                                    self.tree[chrom]['node'][2][i][3].append(gene)
                                break
                    else:
                        self.tree[chrom]['node'][2][len(self.tree[chrom]['node'][2])] = [region,exons,[t_id],[gene],[classificaion]]
                        self.tree[chrom]['node'][0] = [min(self.tree[chrom]['node'][0]+region), max(self.tree[chrom]['node'][0]+region)]
                elif region[0] < self.tree[chrom]['node'][0][0]:  # left
                    if self.tree[chrom]['up']:
                        self.add_NSM_node(chrom,region,exons,t_id,gene,classificaion,node=self.tree[chrom]['up'])
                    else:
                        self.tree[chrom]['up'] = self.init_NSM_node(region,exons,t_id,gene,classificaion)
                else:  # right
                    if self.tree[chrom]['down']:
                        self.add_NSM_node(chrom,region,exons,t_id,gene,classificaion,node=self.tree[chrom]['down'])
                    else:
                        self.tree[chrom]['down'] = self.init_NSM_node(region,exons,t_id,gene,classificaion)
            else:  # 新建chrom
                self.tree[chrom] = self.init_NSM_node(region,exons,t_id,gene,classificaion)
        else:  # 递归添加
            if check_overlap(node['node'][0],region):
                node['node'][1].append(t_id)
                if gene not in node['node'][3]:
                    node['node'][3].append(gene)
                for i in node['node'][2]:
                    if check_overlap(node['node'][2][i][0], region):
                        same, new_ref, new_range = ss_compare(node['node'][2][i][1], exons)
                        if same:
                            node['node'][2][i][2].append(t_id)
                            node['node'][2][i][1] = new_ref
                            node['node'][2][i][0] = [min(new_range+node['node'][2][i][0]),max(new_range+node['node'][2][i][0])]
                            node['node'][2][i][4].append(classificaion)
                            if gene not in node['node'][2][i][3]:
                                node['node'][2][i][3].append(gene)
                            break
                else:
                    node['node'][2][len(self.tree[chrom]['node'][2])] = [region, exons, [t_id],[gene],[classificaion]]
                node['node'][0] = [min(region+node['node'][0]),max(node['node'][0])]
            elif region[0] < node['node'][0][0]:  # left
                if node['up']:
                    self.add_NSM_node(chrom,region,exons,t_id,gene,classificaion,node=node['up'])
                else:
                    node['up'] = self.init_NSM_node(region,exons,t_id,gene,classificaion)
            else:  # right
                if node['down']:
                    self.add_NSM_node(chrom,region,exons,t_id,gene,classificaion,node=node['down'])
                else:
                    node['down'] = self.init_NSM_node(region,exons,t_id,gene,classificaion)

    def add_SM_node(self,chrom,gene,transcript,t_id):
        if chrom not in self.tree:
            self.tree[chrom] = {}
        if gene not in self.tree[chrom]:
            self.tree[chrom][gene] = [[t_id],{}]
        else:
            self.tree[chrom][gene][0].append(t_id)
        if transcript not in self.tree[chrom][gene][1]:
            self.tree[chrom][gene][1][transcript] = [t_id]
        else:
            self.tree[chrom][gene][1][transcript].append(t_id)

    def init_FUSION_node(self,gene,t_id,structure,part_idx):
        genes = sorted(gene.split('|'))
        return [genes,[t_id],{1:[genes,[t_id],structure,[part_idx]]}]

    def add_FUSION_node(self,A):
        refStart = A.REF_START
        refEnd = A.REF_END
        part_idx = A.parts_idx
        gene = A.best_gene
        t_id = A.id
        structure = [[refStart[i],refEnd[i]] for i in range(len(refStart))]
        if not self.tree:
            self.tree[0] = self.init_FUSION_node(gene,t_id,structure,part_idx)
        else:
            for node in self.tree:
                genes_anno = self.tree[node][0]
                genes = sorted(gene.split('|'))
                gene_not_in = []
                in_flag = 0
                for i in genes:
                    if i in genes_anno:
                        in_flag = 1
                    else:
                        gene_not_in.append(i)
                if in_flag:
                    in_flag = 0
                    self.tree[node][0] = self.tree[node][0] + gene_not_in
                    self.tree[node][1].append(t_id)
                    for i in self.tree[node][2]:
                        genes_anno_1 = self.tree[node][2][i][0]
                        structure_anno = self.tree[node][2][i][2]
                        gene_not_in1_1 = [i for i in genes if i not in genes_anno_1]
                        if len(genes_anno_1) == len(genes) and not gene_not_in1_1:
                            same,structure_anno1,structure_range = ss_compare(structure_anno,structure)
                            if same:
                                self.tree[node][2][i][2] = structure_anno1
                                self.tree[node][2][i][1].append(t_id)
                                self.tree[node][2][i][3].append(part_idx)
                                break
                    else:
                        self.tree[node][2][max(list(self.tree[node][2]))+1] = [genes,[t_id],structure,[part_idx]]
                    break
            else:
                self.tree[max(list(self.tree))+1] = self.init_FUSION_node(gene,t_id,structure,part_idx)

    def tree_FUSION_view(self):
        txt = ''
        suffix = '_fusion_'
        for node in self.tree:
            name1 = self.tree[node][0]
            id1 = '|'.join(name1)
            suffix_idx = 0
            for node2 in self.tree[node][2]:
                name2 = self.tree[node][2][node2][0]
                id2 = '|'.join(name2)
                if id2 not in self.FUSION:
                    self.FUSION[id2] = {}
                while '{}{}{}'.format(id2,suffix,suffix_idx) in self.FUSION[id2]:
                    suffix_idx += 1
                class_t = '{}{}{}'.format(id2,suffix,suffix_idx)
                suffix_idx = 0
                self.FUSION[id2][class_t] = {}
                id3 = ','.join(self.tree[node][2][node2][1])
                txt += '{}\t{}\t{}\n'.format(id1,id2,id3)
                for i in range(len(self.tree[node][2][node2][1])):
                    self.FUSION[id2][class_t][self.tree[node][2][node2][1][i]] = self.tree[node][2][node2][3][i]
        return txt

    def tree_NSM_view(self,tree):
        def nodes_Traversal(nodes):
            if nodes['up']:
                nodes_Traversal(nodes['up'])
            res.append(nodes['node'])
            if nodes['down']:
                nodes_Traversal(nodes['down'])
        res = []
        for chrom in tree:
            nodes_Traversal(tree[chrom])
        txt = ''
        for i in res:
            txt += self.stat_NSM_format(i)
        return txt

    def stat_NSM_format(self,info):
        tab1 = '|'.join(info[3])
        txt = ''
        for _idx,i in enumerate(info[2]):
            ids = ','.join(info[2][i][2])
            tab2 = '|'.join(info[2][i][3])
            tab_novel = "{}_novel_{}".format(tab2,_idx)
            txt += '{}\t{}\t{}\t{}\t{}\n'.format(tab1,tab2,tab_novel,len(info[2][i][2]),ids)

            if tab2 in self.NSM:
                self.NSM[tab2][tab_novel] = {j:info[2][i][4][_idx] for _idx,j in enumerate(info[2][i][2])}
            else:
                if tab1 in self.NSM:
                    self.NSM[tab1][tab_novel] = {j:info[2][i][4][_idx] for _idx,j in enumerate(info[2][i][2])}
                else:
                    self.NSM[tab1] = {tab_novel:{j:info[2][i][4][_idx] for _idx,j in enumerate(info[2][i][2])}}
        return txt

    def tree_SM_view(self):
        txt = ''
        for chrom in self.tree:
            for gene in self.tree[chrom]:
                for transcript in self.tree[chrom][gene][1]:
                    ids = ','.join(self.tree[chrom][gene][1][transcript])
                    txt += '{}\t{}\t{}\t{}\n'.format(gene,transcript,len(self.tree[chrom][gene][1][transcript]),ids)
        return txt

    def txt_write(self,txt,out):
        with open(out,'w') as o:
            o.write(txt)


if __name__ == '__main__':
    ss_compare(0,0)