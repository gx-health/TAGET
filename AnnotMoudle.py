"""
Author: Zhang Chengsheng, @2019.12.19
release: 2020.03.25
"""
import utils
import math


def ss_db_patch(exon_db):
    ss_start_db = []
    ss_end_db = []
    for idx in exon_db:
        i,j = exon_db[idx]
        if i not in ss_start_db:
            ss_start_db.append(i)
        if j not in ss_end_db:
            ss_end_db.append(j)
    return sorted(ss_start_db),sorted(ss_end_db)


class Isoform_anno:
    def __init__(self,id,string,fabuffer,faidxDict):
        self.id = id
        self.string = string
        self.fabuffer = fabuffer
        self.faidxDict = faidxDict
        self.SJtxt = ''
        self.chrom_mask = []
        self._origin_string_parse()
        self.isoform_length = sum(self.exons_length)
        self._bias_define()
        self._bias_calc()
        self.merge = 0
        self._pure_info()
        self.anno = {}
        self.db_used = {}
        self.parts_idx = {}
        self.subtypes = {1: '',
                        11: '右边缘外显子',
                        12: '左边缘外显子',
                        18: '左边缘外显子融合',
                        19: '右边缘外显子融合',
                        20: '两端剪切点都已知的新外显子',
                        21: '内含子保留',
                        31: '左端新剪切点',
                        32: '右端新剪切点',
                        35: '左端不对齐延伸外显子',
                        36: '右端不对齐延伸外显子',
                        40: '双端新剪切点',
                        41: '双端新剪切点基因内无注释的新外显子',
                        42: '左边缘外显子新剪切点',
                        43: '右边缘外显子新剪切点',
                        44: '单外显子转录本双端不对齐',
                        50: '边缘双端不对齐', # 已经没有这个了
                        60: '双端新剪切点基因外无注释的新外显子',
                        100: 'BUG_100',
                        101: 'BUG_101',
                        200: '无注释',
                        400: 'BUG_NA'
                        }

    def _origin_string_parse(self):
        """解析并返回原始真实bed分段"""
        THRESHOLD_EXONS_BIAS = 3
        lines = self.string.rstrip().split('\n')
        self.chroms = []
        self.ref_start = []
        self.ref_end = []
        self.seq_start = []
        self.seq_end = []
        self.strands = []
        self.unmapped_length = 0
        self.mask_length = 0
        self.exons_length = []
        self.bias_start = []
        self.bias_end = []
        for line in lines:
            mask = 1
            if mask:
                line = self._0_mask(line)
            s = line.rstrip().split('\t')
            if s[0] == '0':  #0:unmapped,-1:masked
                self.unmapped_length += abs(int(s[7])-int(s[6]))
            elif s[0] == '-1':
                self.mask_length += abs(int(s[7])-int(s[6]))
            self.chroms.append(s[0])
            self.ref_start.append(int(s[1]))
            self.ref_end.append(int(s[2]))
            self.strands.append(s[5])
            self.seq_start.append(int(s[6]))
            self.seq_end.append(int(s[7]))
            self.exons_length.append(int(s[7])-int(s[6])+1)
            self.bias_start.append(0)
            self.bias_end.append(0)

    def _0_mask(self,line):
        """对单外显子短于阈值的做mask处理"""
        MIN_EXON_LENGHT = 10
        s = line.rstrip().split('\t')
        if abs(int(s[7])-int(s[6])) <= MIN_EXON_LENGHT:
            return '-1\t0\t0\t{}\t*\t0\t{}\t{}'.format(s[3],s[6],s[7])
        else:
            return line

    def _bias_define(self):
        """根据exons mask情况确定Splice Junction(SJ)的可变范围"""
        for idx,i in enumerate(self.exons_length):
            if self.chroms[idx] == '-1':
                if idx-1>=0:
                    if self.chroms[idx-1] not in ['-1','0']:
                        if self.strands[idx-1] == '1':
                            self.bias_end[idx-1] = i
                        else:
                            self.bias_start[idx-1] = i
                if idx+1<len(self.exons_length):
                    if self.chroms[idx+1] not in ['-1','0']:
                        if self.strands[idx+1] == '1':
                            self.bias_start[idx+1] = i
                        else:
                            self.bias_end[idx+1] = i

    def _bias_calc(self):
        """计算返回外显子SJ的可信区间"""
        MAX_NATURAL_EXONS_BIAS_THRESHOLD = 0  # 剪切位点误差范围
        self.bias_start_region = []
        self.bias_end_region = []
        for idx, i in enumerate(self.bias_start):
            self.bias_start_region.append([self.ref_start[idx]-i-MAX_NATURAL_EXONS_BIAS_THRESHOLD,self.ref_start[idx]+MAX_NATURAL_EXONS_BIAS_THRESHOLD])
        for idx, i in enumerate(self.bias_end):
            self.bias_end_region.append([self.ref_end[idx]-MAX_NATURAL_EXONS_BIAS_THRESHOLD,self.ref_end[idx]+i+MAX_NATURAL_EXONS_BIAS_THRESHOLD])

    def _pure_info(self,chrom_mask=[]):
        MIN_INTRON = 50
        self.CHROMS = []
        self.REF_START = []
        self.REF_END = []
        self.SEQ_START = []
        self.SEQ_END = []
        self.STRANDS = []
        self.EXONS_LENGTH = []
        self.BIAS_START = []
        self.BIAS_END = []
        self.BIAS_START_REGION = []
        self.BIAS_END_REGION = []
        for idx,i in enumerate(self.chroms):
            if i not in ['0','-1'] and i not in chrom_mask:
                self.CHROMS.append(self.chroms[idx])
                self.REF_START.append(self.ref_start[idx])
                self.REF_END.append(self.ref_end[idx])
                self.SEQ_START.append(self.seq_start[idx])
                self.SEQ_END.append(self.seq_end[idx])
                self.STRANDS.append(self.strands[idx])
                self.EXONS_LENGTH.append(self.exons_length[idx])
                self.BIAS_START.append(self.bias_start[idx])
                self.BIAS_END.append(self.bias_end[idx])
                self.BIAS_START_REGION.append(self.bias_start_region[idx])
                self.BIAS_END_REGION.append(self.bias_end_region[idx])
        delete_idx = []
        for idx,i in enumerate(self.CHROMS):
            if not idx:
                continue
            if self.CHROMS[idx] == self.CHROMS[idx-1] and self.STRANDS[idx] == self.STRANDS[idx-1]:
                strand = self.STRANDS[idx]
                if strand in [-1,'-1']:
                    gap1 = abs(self.REF_START[idx-1] - self.REF_END[idx])
                    gap2 = abs(gap1 - abs(self.SEQ_START[idx] - self.SEQ_END[idx-1]))
                else:
                    gap1 = abs(self.REF_START[idx] - self.REF_END[idx-1])
                    gap2 = abs(gap1 - abs(self.SEQ_START[idx] - self.SEQ_END[idx-1]))
                if gap1 <= MIN_INTRON or gap2 <= MIN_INTRON:
                    new_ref_start = min(self.REF_START[idx],self.REF_START[idx-1],self.REF_END[idx],self.REF_END[idx-1])
                    new_ref_end = max(self.REF_START[idx],self.REF_START[idx-1],self.REF_END[idx],self.REF_END[idx-1])
                    new_seq_start = min(self.SEQ_START[idx],self.SEQ_START[idx-1],self.SEQ_END[idx],self.SEQ_END[idx-1])
                    new_seq_end = max(self.SEQ_START[idx],self.SEQ_START[idx-1],self.SEQ_END[idx],self.SEQ_END[idx-1])
                    self.REF_START[idx-1] = new_ref_start
                    self.REF_END[idx-1] = new_ref_end
                    self.SEQ_START[idx-1] = new_seq_start
                    self.SEQ_END[idx-1] = new_seq_end
                    if new_ref_start == self.REF_START[idx]:
                        self.BIAS_START[idx-1] = self.BIAS_START[idx]
                        self.BIAS_START_REGION[idx-1] = self.BIAS_START_REGION[idx]
                    else:
                        self.BIAS_END[idx-1] = self.BIAS_END[idx]
                        self.BIAS_END_REGION[idx-1] = self.BIAS_END_REGION[idx]
                    delete_idx.append(idx)
        for idx in sorted(delete_idx,reverse=1):
            self._info_refresh(idx)
        if delete_idx: self.merge = 1
        #####
        delete_idx = []
        for idx, i in enumerate(self.CHROMS):
            if not idx or idx+1 == len(self.CHROMS):
                continue
            idx1,idx2 = idx-1,idx+1
            if self.CHROMS[idx1] == self.CHROMS[idx2] and self.STRANDS[idx1] == self.STRANDS[idx2]:
                strand = self.STRANDS[idx1]
                if strand in [1,'1']:
                    gap1 = abs(self.REF_START[idx2] - self.REF_END[idx1])
                    gap2 = abs(gap1 - abs(self.SEQ_START[idx2] - self.SEQ_END[idx1]))
                else:
                    gap1 = abs(self.REF_START[idx1] - self.REF_END[idx2])
                    gap2 = abs(gap1 - abs(self.SEQ_START[idx2] - self.SEQ_END[idx1]))
                if gap1 <= MIN_INTRON or gap2 <= MIN_INTRON:
                    new_ref_start = min(self.REF_START[idx1], self.REF_START[idx2], self.REF_END[idx1],self.REF_END[idx2])
                    new_ref_end = max(self.REF_START[idx1], self.REF_START[idx2], self.REF_END[idx1],self.REF_END[idx2])
                    new_seq_start = min(self.SEQ_START[idx1], self.SEQ_START[idx2], self.SEQ_END[idx1],self.SEQ_END[idx2])
                    new_seq_end = max(self.SEQ_START[idx1], self.SEQ_START[idx2], self.SEQ_END[idx1],self.SEQ_END[idx2])
                    self.REF_START[idx2] = new_ref_start
                    self.REF_END[idx2] = new_ref_end
                    self.SEQ_START[idx2] = new_seq_start
                    self.SEQ_END[idx2] = new_seq_end
                    if new_ref_start == self.REF_START[idx1]:
                        self.BIAS_START[idx2] = self.BIAS_START[idx1]
                        self.BIAS_START_REGION[idx2] = self.BIAS_START_REGION[idx1]
                    else:
                        self.BIAS_END[idx2] = self.BIAS_END[idx1]
                        self.BIAS_END_REGION[idx2] = self.BIAS_END_REGION[idx1]
                    delete_idx.append(idx)
                    delete_idx.append(idx2)
        for idx in sorted(list(set(delete_idx)), reverse=1):
            self._info_refresh(idx)
        if delete_idx: self.merge = 1

    def _info_refresh(self,idx):
        """一次一个删除更新列表，不可逆"""
        self.mask_length += self.seq_end[idx] - self.seq_start[idx] + 1
        del self.CHROMS[idx]
        del self.REF_START[idx]
        del self.REF_END[idx]
        del self.SEQ_START[idx]
        del self.SEQ_END[idx]
        del self.STRANDS[idx]
        del self.EXONS_LENGTH[idx]
        del self.BIAS_START[idx]
        del self.BIAS_END[idx]
        del self.BIAS_START_REGION[idx]
        del self.BIAS_END_REGION[idx]

    def _break_point_cluster(self):
        """isoform根据断点(异染色体、大段gap、不同链性、染色体重叠，等)分类"""
        classification = 1
        if not self.CHROMS:
            return classification
        if len(self.CHROMS) > 1:  # mulit-exon
            chrom_num = self._multi_valid_chroms_check()
            gap_num, self.parts_idx = self._gap_check()
            if not self.CHROMS:  # 排除一种分成几段但是哪段都不行的情况，p.s.但是把单染色体映射短于200的也过滤掉了
                return 1
            if chrom_num > 1:  # 多染色体映射
                if gap_num > 4:
                    return -5  # 复杂的融合
                return -2  # 跨染色体融合
            else:  # 单染色体映射
                if (self.mask_length + self.unmapped_length) / self.isoform_length > 0.5:
                    return 2  # 无效区域过长
                if gap_num > 4:
                    return 3  # 复杂的新转录本
                elif gap_num:  # 分成两段及以上
                    #  疑似同染色体远距离融合||新转录本
                    return -1
                elif len(self.CHROMS) == 1:  # 单外显子
                    return 10
                else:  # 仅一段
                    # 正常转录本
                    return 0

        ##################### 0 ######################
        else:  # mono-exon
            chrom_num = self._multi_valid_chroms_check()
            gap_num, self.parts_idx = self._gap_check()
            return 10

    def _debug_run(func):
        def wrap(self,*args):
            try:
                res = func(self,*args)
                return res
            except:
                print(self.id)
                return [0,0,0,0,0,0,0]
        return wrap

    def _multi_valid_chroms_check(self):
        """检查是否存在多染色体映射，同时检查其长度 将低于阈值的分段移除"""
        if len(set(self.CHROMS)) == 1:
            return 1
        FUSION_MIN_LENGTH_THRESHOLD = 200
        self.CHROMS_CLUSTER = {}
        for _idx,chrom in enumerate(self.CHROMS):
            if chrom not in self.CHROMS_CLUSTER:
                self.CHROMS_CLUSTER[chrom] = self.EXONS_LENGTH[_idx]
            else:
                self.CHROMS_CLUSTER[chrom] += self.EXONS_LENGTH[_idx]
        chrom_num = 0
        for chrom in self.CHROMS_CLUSTER:
            if self.CHROMS_CLUSTER[chrom] > FUSION_MIN_LENGTH_THRESHOLD:
                chrom_num += 1
            else:
                self.chrom_mask.append(chrom)
        self._pure_info(chrom_mask=self.chrom_mask)
        return chrom_num

    def _calc_length(self,dict_in):
        length = 0
        for i in dict_in:
            length += abs(dict_in[i][-1] - dict_in[i][0])+1
        return length

    def _gap_check(self):
        start = self.REF_START
        end = self.REF_END
        strand = self.STRANDS
        chrom = self.CHROMS
        MAX_GAP_LENGTH = 1000000
        MIN_ISOFORM_LENGTH = 200
        res_idx = {0:[]}
        indicator = 0
        for _idx,i in enumerate(start):
            if not i:
                pass
            if not _idx:
                res_idx[indicator].append(_idx)
                continue
            gap1 = abs(start[_idx] - start[_idx-1])
            gap2 = abs(end[_idx] - end[_idx - 1])
            if gap1 > MAX_GAP_LENGTH or gap2 > MAX_GAP_LENGTH or strand[_idx] != strand[_idx-1] or chrom[_idx] != chrom[_idx-1]:
                indicator += 1
                res_idx[indicator] = [_idx]
            else:
                res_idx[indicator].append(_idx)
        ###### 段内长度检查 ######
        remove_list = []
        idx_list = []
        for i in res_idx:
            length = 0
            for idx in res_idx[i]:
                length += abs(end[idx] - start[idx])
            if length < MIN_ISOFORM_LENGTH:
                for idx in res_idx[i]:
                    idx_list.append(idx)
                remove_list.append(i)
        for i in remove_list:
            res_idx.pop(i)
        for i in sorted(idx_list,reverse=1):
            self._info_refresh(i)
        if idx_list:
            gap_num, res_idx = self._gap_check()
        else:
            gap_num = len(res_idx)-1
        return gap_num,res_idx

    def anno_used_add(self,db,chrom):
        """准备后续画图用的注释库"""
        for gene in db:
            if gene not in self.db_used:
                gene_type = db[gene][0]
                strand = db[gene][1]
                ts = db[gene][3]
                self.db_used[gene] = {}
                for t in ts:
                    if t not in self.db_used[gene]:
                        self.db_used[gene][t] = {}
                        exons = ts[t][2]
                        chroms = [chrom] * len(exons)
                        strands = [strand] * len(exons)
                        left = [exons[i][0] for i in exons]
                        right = [exons[i][1] for i in exons]
                        anno = [''] * len(exons)
                        self.db_used[gene][t]['info'] = {}
                        self.db_used[gene][t]['info'] = [gene, t, gene_type, '']
                        self.db_used[gene][t]['chrom'] = chroms
                        self.db_used[gene][t]['strand'] = strands
                        self.db_used[gene][t]['start'] = left
                        self.db_used[gene][t]['end'] = right
                        self.db_used[gene][t]['anno'] = anno

    def reads_info(self):
        info = {}
        info['info'] = [self.best_gene,self.best_transcript,self.classification,self.subtype]
        info['chrom'] = self.CHROMS
        info['strand'] = self.STRANDS
        info['start'] = self.REF_START
        info['end'] = self.REF_END
        info['anno'] = self.exon_anno
        info['parts_idx'] = self.parts_idx
        return info

    def init_parts_anno(self,flag=0):
        self.EXONS_TYPE = [200 for i in self.REF_START]
        self.correct_REF_START = [-1 for i in self.REF_START]
        self.correct_REF_END = [-1 for i in self.REF_START]
        self.GTAG_START = [-1 for i in self.REF_START]
        self.GTAG_END = [-1 for i in self.REF_START]
        self.SJknown = [-1 for i in range(len(self.REF_START)-1)]
        self.SJcanonical = [-1 for i in range(len(self.REF_START)-1)]
        self.gene_left_right_edge_flag = [0 for i in self.REF_START]
        self.best_anno_info = ['NA','NA','NA','NA','NA','NA','NA','NA']
        self.classification = 'NA'
        self.subtype = 'NA'
        self.best_gene = 'NA'
        self.best_transcript = 'NA'
        self.multiAnno = {}
        self.multiAnnoTXT = ''
        self.anno_chrom = 'NA'
        self.anno_strand = 'NA'
        self.diff_to_gene_start = 'NA'
        self.diff_to_gene_end = 'NA'
        self._5_ref_draft = 'NA'
        self._3_ref_draft = 'NA'
        self._5_seq_draft = 'NA'
        self._3_seq_draft = 'NA'
        self.exon_anno = ['' for i in self.REF_START]
        self.parts_anno = {}
        self.part_best = {}
        if flag <= 0 or flag == 10:
            for part in self.parts_idx:
                self.parts_anno[part] = {}
                self.part_best[part] = []
                for idx in self.parts_idx[part]:
                    self.parts_anno[part][idx] = {'exon_type': 200,
                                                'correct_REF_START': -1,
                                                'correct_REF_END': -1,
                                                'ss_start': 0,
                                                'ss_end': 0,
                                                'exon_miss_5': 'NA',
                                                'exon_miss_3': 'NA',
                                                'GTAG_START': 'NA',
                                                'GTAG_END': 'NA',
                                                }

    def anno_func_v1(self,REF_START, REF_END, CHROMS,STRANDS, BIAS_START_REGION, BIAS_END_REGION, anno_db,part,idxs,mono=0):
        best_gene = 'NA'
        best_transcript = 'NA'
        if not anno_db:
            best_gene = '{}:{}-{}'.format(CHROMS[0], min(REF_START), max(REF_END))
            for IDX in idxs:
                self.parts_anno[part][IDX]['exon_type'] = 60
            return best_gene, [best_transcript,0,['?'],[1],['NA','NA','NA','NA','NA','NA','NA','NA']]
        elif len(anno_db) > 5:
            ## 基因注释过多，可能会卡，则在此设限
            pass

        def sj_db_mk(transcripts):
            sj_db = []
            for transcript in transcripts:
                exon_idx = transcripts[transcript][3]
                for idx in range(len(exon_idx)):
                    if not idx:
                        continue
                    sj = '{}-{}'.format(exon_idx[idx-1],exon_idx[idx])
                    if sj not in sj_db:
                        sj_db.append(sj)
            return sj_db

        def sj_check(sj_db,sjs):
            sj_flag = []
            wtf = []
            flag = 1
            for idx in range(len(sjs)):
                if not isinstance(sjs[list(sjs)[idx]], int):
                    if 'U' in sjs[list(sjs)[idx]] or 'NA' in sjs[list(sjs)[idx]]:
                        sj_flag.append(0)
                        if idx:
                            if '-' in str(sjs[list(sjs)[idx]]):
                                sj = '{}-{}'.format(sjs[list(sjs)[idx-1]],sjs[list(sjs)[idx]].split('-')[0])
                                if sj not in sj_db:
                                    wtf.append(0)
                                else:
                                    wtf.append(1)
                if not idx:
                    continue

                st = sjs[list(sjs)[idx - 1]].split('-')[-1] if '-' in str(sjs[list(sjs)[idx-1]]) else sjs[list(sjs)[idx]]
                sj = '{}-{}'.format(sjs[list(sjs)[idx-1]],st)
                if sj not in sj_db:
                    wtf.append(0)
                    sj_flag.append(0)
                else:
                    wtf.append(1)
                    sj_flag.append(1)
            if sj_flag and 0 in sj_flag:
                flag = 0
            return flag,wtf

        def ss_in_or_not(ss_db_in,point_start,point_end):
            """检查断点是否在目标基因区段内"""
            a,b = min(ss_db_in),max(ss_db_in)
            flag_start = 1 if a < point_start[0] < b or a < point_start[-1] < b else 0
            flag_end = 1 if a < point_end[0] < b or a < point_end[-1] < b else 0
            flag_end = 1 if min(point_start) < a < max(point_end) or min(point_start) < b < max(point_end) else flag_end
            return flag_start,flag_end

        def gene_left_right_edge_flag_make(gene_left_right_edge_flag,gene_exon_in):
            left_idx,right_idx = 0,0
            for _idx,IDX in enumerate(gene_left_right_edge_flag):
                if gene_exon_in[_idx] != 0:
                    left_idx = IDX
                    break
            for _idx, IDX in enumerate(list(gene_left_right_edge_flag)[::-1]):
                if gene_exon_in[::-1][_idx] != 0:
                    right_idx = IDX
                    break
            if STRANDS[0] in [1,'1']:
                gene_left_right_edge_flag[left_idx] = -1 if mono != 10 else 0
                gene_left_right_edge_flag[right_idx] = 1 if mono != 10 else 0
            else:
                gene_left_right_edge_flag[left_idx] = 1 if mono != 10 else 0
                gene_left_right_edge_flag[right_idx] = -1 if mono != 10 else 0
            return gene_left_right_edge_flag

        def exon_check(exon_db,start,end,start_flag=0):
            """用于检查每个外显子与注释的比较情况"""
            EXON_DRAFT_1 = 10
            idx_start_list = {}
            idx_end_list = {}
            ss_start,ss_end = 0,0  # 剪切点是否已知,0:novel 1:known
            for _x,idx in enumerate(exon_db):
                ref_start, ref_end = exon_db[idx]
                if start[0] <= ref_start <= start[-1]:
                    start_draft = 0
                else:
                    start_draft = min(start[0] - ref_start, start[-1] - ref_start) if start[0] - ref_start > 0 else max(start[0] - ref_start, start[-1] - ref_start)
                if end[0] <= ref_end <= end[-1]:
                    end_draft = 0
                else:
                    end_draft = min(end[0] - ref_end, end[-1] - ref_end) if end[0] - ref_end > 0 else max(end[0] - ref_end, end[-1] - ref_end)
                flag_start_1 = 1 if abs(start_draft) <= EXON_DRAFT_1 else 0
                flag_end_1 = 1 if abs(end_draft) <= EXON_DRAFT_1 else 0
                if not _x:
                    start_draft_min = start_draft
                    start_point_most = exon_db[idx][0]
                    end_draft_min = end_draft
                    end_point_most = exon_db[idx][-1]
                else:
                    start_draft_min, start_point_most = [start_draft, exon_db[idx][0]] if abs(start_draft) <= abs(start_draft_min) else [start_draft_min, start_point_most]
                    end_draft_min, end_point_most = [end_draft, exon_db[idx][-1]] if abs(end_draft) <= abs(end_draft_min) else [end_draft_min, end_point_most]
                    if end_point_most < start_point_most:
                        end_draft_min, end_point_most = [end_draft, exon_db[idx][-1]]
                if start_flag == -1:
                    if flag_end_1:
                        idx_end_list[idx] = end_draft
                        idx_start_list[idx] = start_draft
                    continue
                elif start_flag == 1:
                    if flag_start_1:
                        idx_start_list[idx] = start_draft
                        idx_end_list[idx] = end_draft
                    continue
                if flag_start_1:
                    idx_start_list[idx] = start_draft
                if flag_end_1:
                    idx_end_list[idx] = end_draft
            exon_id_anno = {}  # {exon_idx:[left_draft,rigth_draft,subtype]}
            ## subtype: 0:正常，1:双端对齐，2:单端对齐,3:双端不对齐
            default_score = 1000
            if idx_start_list and idx_end_list:
                ss_start,ss_end = 1,1
                for m in idx_start_list:
                    for n in idx_end_list:
                        id,subtype = [m,0] if m == n else ['{}-{}'.format(m,n),1]
                        exon_id_anno[id] = [idx_start_list[m],idx_end_list[n],subtype]
            elif idx_start_list:
                ss_start,ss_end = 1,0
                for m in idx_start_list:
                    id = '{}-SU'.format(m)
                    exon_id_anno[id] = [idx_start_list[m],default_score,2]
            elif idx_end_list:
                ss_start,ss_end = 0,1
                for m in idx_end_list:
                    id = 'SU-{}'.format(m)
                    exon_id_anno[id] = [default_score,idx_end_list[m],2]
            else:
                ss_start,ss_end = 0,0
                id = 'DU'
                exon_id_anno[id] = [default_score,default_score,3]
            return exon_id_anno,ss_start,ss_end

        def transcript_check(transcripts,gene_exon_idx_anno,sj_db):
            best_res = []
            for transcript in transcripts:
                exon_idxs = transcripts[transcript][3]
                t_start,t_end = transcripts[transcript][1]
                t_map_list = []  # 外显子与注释的匹配情况，-2超范围，正常id：外显子编号，-1没匹配上
                for _idx,IDX in enumerate(gene_exon_idx_anno):
                    ov = self.check_overlap([t_start,t_end],[self.REF_START[IDX],self.REF_END[IDX]])
                    if not ov:
                        t_map_list.append(-2)
                        continue
                    else:
                        for id in gene_exon_idx_anno[IDX]:
                            if id in exon_idxs:
                                t_map_list.append(id)
                                break
                        else:
                            t_map_list.append(-1)
                if self.list_intersection(t_map_list, exon_idxs):  # 该转录本可作为候选最佳注释
                    ref_exon_num = len(exon_idxs)
                    e_idxs = exon_idxs.index(t_map_list[0]),exon_idxs.index(t_map_list[-1])
                    exon_miss_5,exon_miss_3 = min(e_idxs),len(exon_idxs)-max(e_idxs)-1
                    score = sum([(abs(gene_exon_idx_anno[IDX][t_map_list[idx]][0]) + abs(gene_exon_idx_anno[IDX][t_map_list[idx]][1])) for idx,IDX in enumerate(list(gene_exon_idx_anno))])
                    miss_exon_num = ref_exon_num - len(REF_START)
                    reflength = self._calc_length(transcripts[transcript][2])
                    diff_to_transcript_end = max(REF_END) - t_end
                    diff_to_transcript_start = t_start - min(REF_START)
                    strand = STRANDS[0]
                    diff_to_transcript_start, diff_to_transcript_end = [diff_to_transcript_start,diff_to_transcript_end] if strand == 1 or strand == '1' else [diff_to_transcript_end, diff_to_transcript_start]
                    best_res.append([transcript,miss_exon_num,score,t_map_list,ref_exon_num,reflength,diff_to_transcript_start, diff_to_transcript_end,exon_miss_5,exon_miss_3])
            best_transcript = 'NA'
            exon_type_dict = {}
            best_gene_exon_idx_anno = {}
            if best_res:
                best_res = sorted(best_res,key=lambda x:(x[1],x[2]))
                self.multiAnnot(gene=gene,res=best_res,mode=1)  # add in V0.0.2
                best_transcript_anno = best_res[0]
                transcript, miss_exon_num, score, t_map_list, ref_exon_num, reflength, diff_to_transcript_start, diff_to_transcript_end,exon_miss_5,exon_miss_3 = best_transcript_anno
                anno_info = [ref_exon_num, reflength, diff_to_transcript_start, diff_to_transcript_end,exon_miss_5,exon_miss_3]
                for _idx,IDX in enumerate(gene_exon_idx_anno):
                    best_gene_exon_idx_anno[IDX] = [t_map_list[_idx],gene_exon_idx_anno[IDX][t_map_list[_idx]]]
                best_transcript = best_transcript_anno[0]
                best_idx = best_transcript_anno[3]
                for IDX in gene_exon_idx_anno:
                    if gene_left_right_edge_flag[IDX]:
                        edge_add = 0 if gene_left_right_edge_flag[IDX] == 1 else 1
                        overlap_exons = []
                        for _idx1 in transcripts[best_transcript][2]:
                            if self.check_overlap(transcripts[best_transcript][2][_idx1], [self.REF_START[IDX],self.REF_END[IDX]]):
                                overlap_exons.append(_idx1)
                        if len(overlap_exons) > 1:
                            exon_type_dict[IDX] = 18+edge_add  # 边缘外显子融合  18左端，19右端
                        else:
                            exon_type_dict[IDX] = 11+edge_add  # 正常边缘外显子  11左端，12右端
                    else:
                        exon_type_dict[IDX] = 1
            else:
                miss_exon_num = 10000
                best_idx = dict([(IDX,sorted(gene_exon_idx_anno[IDX], key=lambda x: (gene_exon_idx_anno[IDX][x][2],abs(gene_exon_idx_anno[IDX][x][1]) + abs(gene_exon_idx_anno[IDX][x][0])))[0]) for IDX in gene_exon_idx_anno])
                score = sum([abs(gene_exon_idx_anno[IDX][best_idx[IDX]][0])+abs(gene_exon_idx_anno[IDX][best_idx[IDX]][1]) for IDX in gene_exon_idx_anno])
                anno_info = ['NA','NA','NA','NA','NA','NA']
                t_map_list = [best_idx[i] for i in best_idx]
                for _idx,IDX in enumerate(gene_exon_idx_anno):
                    best_gene_exon_idx_anno[IDX] = [t_map_list[_idx],gene_exon_idx_anno[IDX][t_map_list[_idx]]]
                for IDX in best_idx:
                    if isinstance(best_idx[IDX],int):  # 正常外显子
                        exon_type = 1
                        if gene_left_right_edge_flag[IDX]:
                            edge_add = 0 if gene_left_right_edge_flag[IDX] == 1 else 1
                            exon_type = 11 + edge_add  # 正常边缘外显子  11左端，12右端
                    elif 'U' in best_idx[IDX]:  # 双端不对齐+右端不对齐+左端不对齐
                        overlap_exons = []
                        for _idx1 in exon_db:
                            if self.check_overlap(exon_db[_idx1], [self.REF_START[IDX], self.REF_END[IDX]]):
                                overlap_exons.append(_idx1)
                        if 'DU' in best_idx[IDX]:  # 双端不对齐
                            exon_type = 40
                            if gene_left_right_edge_flag[IDX]:
                                edge_add = 0 if gene_left_right_edge_flag[IDX] == 1 else 1
                                exon_type = 42 + edge_add
                            elif not overlap_exons:
                                exon_type = 41  # 平白无故多出来的外显子
                            if mono == 10:
                                exon_type = 44  # 单外显子转录本双端不对齐
                        elif 'SU' in best_idx[IDX]:  # 单端不对齐
                            bdq_add = 0 if 'SU-' in best_idx[IDX] else 1  # 左端不对齐+0，右端不对齐+1
                            if not overlap_exons:
                                exon_type = 41  # 平白无故多出来的外显子，当然这种情况是不会发生的
                            else:
                                exon_max_length = 0
                                for _idx2 in overlap_exons:
                                    exon_max_length = abs(exon_db[_idx2][-1] - exon_db[_idx2][0]) if abs(exon_db[_idx2][-1] - exon_db[_idx2][0]) > exon_max_length else exon_max_length
                                if abs(self.REF_END[IDX] - self.REF_START[IDX]) - exon_max_length > 50:
                                    exon_type = 35+bdq_add  ## exon_extension
                                else:
                                    exon_type = 31+bdq_add  ## 正常单端对齐外显子，短一截，后续需要讨论exon二合一问题
                        else:  # 不存在这种情况
                            exon_type = 101
                    elif '-' in best_idx[IDX]:  # 两段端都对齐但分属不同外显子
                        exon_type = 20  # 两端都对齐但分属不同外显子，又不是外显子保留
                        se_idx = [int(i) for i in best_idx[IDX].split('-')]
                        for transcript in transcripts:
                            tsc_idx = transcripts[transcript][3]
                            if se_idx[0] in tsc_idx and se_idx[1] in tsc_idx:
                                exon_type = 21  # 内含子保留
                                break
                    else:  # 不知道是啥
                        exon_type = 100
                    exon_type_dict[IDX] = exon_type
            sj_flag,sjs = 0,[]
            if [i for i in exon_type_dict if exon_type_dict[i] > 19]:
                sj_flag,sjs = sj_check(sj_db, best_idx)
            return best_transcript,score,sj_flag,sjs,exon_type_dict,best_gene_exon_idx_anno,anno_info,miss_exon_num

        anno_score = {}
        for gene in anno_db:
            anno_score[gene] = {}
            gene_region = anno_db[gene][2]
            transcripts = anno_db[gene][3]
            ss_db = anno_db[gene][5]
            exon_db = anno_db[gene][4]
            sj_db = sj_db_mk(transcripts)
            gene_exon_start_in,gene_exon_end_in = {},{}  # 外显子是否在参考注释区段内{idx:flag} 1:in,0:out
            gene_exon_id_anno = {}  # 外显子注释参考{idx:{exon_idx:[left_draft,right_draft,subtype]}} draft:seq-ref，更多注释在exon_check()
            gene_ss_start,gene_ss_end = {},{}  # 外显子剪切点是否已知，用于判定NNC。{idx:flag} 1:known,0:novel
            sj_flag = 0  # 是否存在新的SJ，1:存在，0:不存在
            gene_left_right_edge_flag = {}  # 是否为边缘外显子，-1:左端，1:右端,0：中间
            for _idx,IDX in enumerate(idxs):
                flag_start, flag_end = ss_in_or_not(ss_db, BIAS_START_REGION[_idx], BIAS_END_REGION[_idx])
                gene_exon_start_in[IDX] = flag_start
                gene_exon_end_in[IDX] = flag_end
                gene_left_right_edge_flag[IDX] = 0
            gene_exon_in = [gene_exon_start_in[i] + gene_exon_end_in[i] for i in gene_exon_start_in]
            gene_left_right_edge_flag = gene_left_right_edge_flag_make(gene_left_right_edge_flag, gene_exon_in)
            for _idx,IDX in enumerate(idxs):
                exon_id_anno,ss_start,ss_end = exon_check(exon_db, BIAS_START_REGION[_idx], BIAS_END_REGION[_idx], start_flag=gene_left_right_edge_flag[IDX])
                gene_exon_id_anno[IDX] = exon_id_anno
                gene_ss_start[IDX] = ss_start
                gene_ss_end[IDX] = ss_end
            best_transcript, score, sj_flag,sjs, exon_type_dict, best_gene_exon_idx_anno, anno_info,miss_exon_num = transcript_check(transcripts, gene_exon_id_anno,sj_db)
            diff_to_gene_start = max(REF_END) - max(gene_region)
            diff_to_gene_end = min(gene_region) - min(REF_START)
            diff_to_gene_start, diff_to_gene_end = [diff_to_gene_start, diff_to_gene_end] if STRANDS[0] == 1 else [diff_to_gene_end, diff_to_gene_start]
            anno_info = anno_info + [diff_to_gene_start, diff_to_gene_end]
            anno_score[gene] = [best_transcript,score,best_gene_exon_idx_anno,sj_flag,sjs,gene_ss_start,gene_ss_end,exon_type_dict,gene_exon_start_in,gene_exon_end_in,anno_info,gene_left_right_edge_flag,miss_exon_num]

        for idx,gene in enumerate(sorted(anno_score,key=lambda x:(anno_score[x][-1],anno_score[x][1]))):
            if anno_score[gene][0] != 'NA':
                best_gene = gene
                break
            else:
                if not idx:
                    best_gene = gene
                else:
                    if anno_score[gene][1] == 'NA':
                        continue
                    elif anno_score[best_gene][1] == 'NA':
                        best_gene = gene
                    elif anno_score[gene][1] < anno_score[best_gene][1]:
                        best_gene = gene
        best_transcript, score, best_gene_exon_idx_anno,sj_flag,sjs, gene_ss_start, gene_ss_end, exon_type_dict, gene_exon_start_in, gene_exon_end_in, anno_info,gene_left_right_edge_flag,miss_exon_num = anno_score[best_gene]
        gene_exon_in = [gene_exon_start_in[i] + gene_exon_end_in[i] for i in gene_exon_start_in]
        if not sum(gene_exon_in):
            best_gene = '{}:{}-{}'.format(CHROMS[0], min(REF_START), max(REF_END))
            return best_gene, ['NA',0,['?'],[1],['NA','NA','NA','NA','NA','NA','NA','NA']]
        for IDX in best_gene_exon_idx_anno:
            self.EXONS_TYPE[IDX] = exon_type_dict[IDX]
            correct_REF_START,correct_REF_END = [anno_db[best_gene][4][best_gene_exon_idx_anno[IDX][0]][0],anno_db[best_gene][4][best_gene_exon_idx_anno[IDX][0]][1]] if isinstance(best_gene_exon_idx_anno[IDX][0],int) else [-1,-1]
            self.correct_REF_START[IDX],self.correct_REF_END[IDX] = correct_REF_START,correct_REF_END
            self.GTAG_START = -1  # TODO
            self.GTAG_END = -1  # TODO
            self.parts_anno[part][IDX] = {'exon_type': exon_type_dict[IDX],
                                          'correct_REF_START': correct_REF_START,
                                          'correct_REF_END': correct_REF_END,
                                          'ss_start': gene_ss_start[IDX],
                                          'ss_end': gene_ss_end[IDX],
                                          'exon_miss_5': 'NA',
                                          'exon_miss_3': 'NA',
                                          'GTAG_START': 'NA',
                                          'GTAG_END': 'NA',
                                          }
        return best_gene,[best_transcript,sj_flag,sjs,gene_exon_in,anno_info]

    #@_debug_run
    def annotation(self,db,flag,debug=0):
        if debug:
            #transcript_boxplot.iso_anno(bed_test,db)
            pass
        self.init_parts_anno(flag)
        MAX_RANGE = 1000000  # 单部分最大跨度
        _5_flank, _3_flank = 'NA', 'NA'
        best_gene, best_transcript, reflength, refexon_num = 'NA','NA','NA','NA'
        best_anno_info = ['NA','NA','NA','NA','NA','NA','NA','NA']
        self.correct_REF_START = self.REF_START.copy()
        self.correct_REF_END = self.REF_END.copy()

        def _0_in_ss_flag_in(ss_flag_in,idxs):
            res = []
            idx = []
            for _idx, i in enumerate(ss_flag_in):
                if i == 0:
                    res.append(_idx)
                    idx.append(idxs[_idx])
            return res,idx

        if flag <= 0 or flag == 10:  # 正常转录本
            anno_res = []
            if (flag >= 0 and len(self.parts_idx) > 1) or not self.parts_idx:
                self.classification = 'UNKNOWN'
                self.subtype = ' UNKNOWN'
                return 0
            for part in self.parts_idx:
                idxs = self.parts_idx[part]
                REF_START = self._get_parts_by_idx(self.REF_START, idxs)
                REF_END = self._get_parts_by_idx(self.REF_END, idxs)
                CHROMS = self._get_parts_by_idx(self.CHROMS, idxs)
                STRANDS = self._get_parts_by_idx(self.STRANDS, idxs)
                BIAS_START_REGION = self._get_parts_by_idx(self.BIAS_START_REGION, idxs)
                BIAS_END_REGION = self._get_parts_by_idx(self.BIAS_END_REGION, idxs)
                if max(REF_END) - min(REF_START) > MAX_RANGE:
                    print(self.id, '注释区域长度超标！')
                    self.classification = 'UNKNOWN'
                    self.subtype = ' UNKNOWN'
                    return 0
                anno_db = self.db_search_by_region(db, CHROMS[0], min(REF_START), max(REF_END))
                self.anno_used_add(anno_db, CHROMS[0])
                ## 可在此加上判断注释基因数量模块
                ## 对注释中的每个基因进行判定
                best_gene, best_anno = self.anno_func_v1(REF_START, REF_END, CHROMS, STRANDS, BIAS_START_REGION,BIAS_END_REGION, anno_db, part, idxs, mono=flag)
                best_transcript, sj_flag,sjs, gene_exon_in, anno_info = best_anno
                anno_res.append([best_gene, best_transcript, sj_flag,sjs, gene_exon_in, anno_info,part,CHROMS[0],STRANDS[0]])
                if gene_exon_in == [1]:
                    juncIDX = idxs
                else:
                    juncIDX = [idxs[i] for i in range(len(idxs)) if gene_exon_in[i]]
                for _i in range(len(juncIDX)):
                    if not _i:
                        continue
                    if juncIDX[_i] - juncIDX[_i - 1] == 1:
                        if not sjs:
                            self.SJknown[juncIDX[_i - 1]] = 1
                        elif '?' in sjs:
                            self.SJknown[juncIDX[_i - 1]] = 0
                        else:
                            self.SJknown[juncIDX[_i - 1]] = sjs[_i - 1]

                flag_in_idx,flag_idx = _0_in_ss_flag_in(gene_exon_in,idxs)  # 数位表示外显子，1为在范围内，0为不在
                cc = 0
                pop_list = [best_gene]
                if not flag_in_idx:  # 没啥融合的转录本
                    pass
                while flag_in_idx:
                    REF_START = self._get_parts_by_idx(REF_START, flag_in_idx)
                    REF_END = self._get_parts_by_idx(REF_END, flag_in_idx)
                    CHROMS = self._get_parts_by_idx(CHROMS, flag_in_idx)
                    STRANDS = self._get_parts_by_idx(STRANDS, flag_in_idx)
                    BIAS_START_REGION = self._get_parts_by_idx(BIAS_START_REGION, flag_in_idx)
                    BIAS_END_REGION = self._get_parts_by_idx(BIAS_END_REGION, flag_in_idx)
                    anno_db_new = self.db_search_by_region(db, CHROMS[0], min(REF_START), max(REF_END))
                    for i in pop_list:
                        if i in anno_db_new:
                            anno_db_new.pop(i)
                    self.anno_used_add(anno_db_new, CHROMS[0])
                    best_gene, best_anno = self.anno_func_v1(REF_START, REF_END, CHROMS, STRANDS, BIAS_START_REGION,BIAS_END_REGION, anno_db_new, part, flag_idx, mono=flag)
                    pop_list.append(best_gene)
                    best_transcript, sj_flag,sjs, gene_exon_in, anno_info = best_anno
                    anno_res.append([best_gene, best_transcript, sj_flag,sjs, gene_exon_in, anno_info,part,CHROMS[0],STRANDS[0]])
                    if len(flag_in_idx) != len(gene_exon_in):
                        # 这段代码是因为novel intergenic transcript返回的gene_exon_in 为[1]与flag_in_idx不等长，为了省事这么写了
                        juncIDX = flag_in_idx
                    else:
                        juncIDX = [flag_in_idx[i] for i in range(len(flag_in_idx)) if gene_exon_in[i]]
                    for _i in range(len(juncIDX)):
                        if not _i:
                            continue
                        if juncIDX[_i] - juncIDX[_i - 1] == 1:
                            if not sjs:
                                self.SJknown[juncIDX[_i - 1]] = 1
                            elif '?' in sjs:
                                self.SJknown[juncIDX[_i - 1]] = 0
                            else:
                                self.SJknown[juncIDX[_i - 1]] = sjs[_i - 1]
                    flag_in_idx,flag_idx = _0_in_ss_flag_in(gene_exon_in,idxs)
                    cc += 1
                    if cc > 5:
                        print(self.id, '循环超限！')
                        break
                for IDX in self.parts_anno[part]:
                    self.exon_anno[IDX] = self.subtypes[self.parts_anno[part][IDX]['exon_type']]
            SM_NC_EXON_NUM = 1
            SM_NC_MINUES_LENGTH = -300
            SM_NC_ADD_LENGTH = 100
            chroms_all = [str(i[7]) for i in anno_res]
            strand_all = [str(i[8]) for i in anno_res]
            best_gene_all = [str(i[0]) for i in anno_res]
            best_chrom = '|'.join(chroms_all)
            best_strand = '|'.join(strand_all)
            best_gene = '|'.join(best_gene_all)
            self.anno_chrom = best_chrom
            self.anno_strand = best_strand
            if flag < 0:
                classification = 'FUSION'
                subtype = 'INNER_CHROM_FUSION' if len(set(chroms_all)) == 1 else 'INTER_CHROM_FUSION'
            elif flag == 10:
                subtype = 'mono_exon'
                classification = 'Intergenic'
                for i in anno_res:
                    if ':' not in i[0]:
                        classification = 'Genic'
                if len(anno_res) == 1 and best_transcript not in ['NA','',0]:
                    ref_exon_num, reflength, diff_to_transcript_start, diff_to_transcript_end, exon_miss_5, exon_miss_3, diff_to_gene_start, diff_to_gene_end = anno_res[0][5]
                    best_anno_info = [ref_exon_num, reflength, diff_to_transcript_start, diff_to_transcript_end,exon_miss_5, exon_miss_3, diff_to_gene_start, diff_to_gene_end]
                    ISM_5_EXON = 1 if exon_miss_5 != 'NA' and exon_miss_5 >= SM_NC_EXON_NUM else 0
                    ISM_3_EXON = 1 if exon_miss_3 != 'NA' and exon_miss_3 >= SM_NC_EXON_NUM else 0
                    classification = 'ISM'
                    if ISM_5_EXON and ISM_3_EXON:
                        subtype = 'INTERNAL_ISM_EXON'
                    elif ISM_3_EXON:
                        subtype = '3_ISM_EXON'
                    elif ISM_5_EXON:
                        subtype = '5_ISM_EXON'
                    elif diff_to_transcript_start > SM_NC_ADD_LENGTH and diff_to_transcript_end > SM_NC_ADD_LENGTH:
                        subtype = '53_FLANK'
                    elif diff_to_transcript_start > SM_NC_ADD_LENGTH:
                        subtype = '5_FLANK'
                    elif diff_to_transcript_end > SM_NC_ADD_LENGTH:
                        subtype = '3_FLANK'
                    elif diff_to_transcript_start <= SM_NC_MINUES_LENGTH and diff_to_transcript_end <= SM_NC_MINUES_LENGTH:
                        subtype = 'INTERNAL_ISM'
                    elif diff_to_transcript_start <= SM_NC_MINUES_LENGTH:
                        subtype = '5_ISM'
                    elif diff_to_transcript_end <= SM_NC_MINUES_LENGTH:
                        subtype = '3_ISM'
                    else:
                        classification = 'FSM'
                        subtype = 'FSM'
            elif len(anno_res) > 1:
                ## 注释结果数大于1时应判断其余是否为intergenic
                for i in anno_res[1:]:
                    #print(i)
                    if ':' not in i[0]:
                        classification = 'FUSION'
                        subtype = 'NEIGHBOUR_FUSION'
                        best_gene = '|'.join([i[0] for i in anno_res])
                        break
                else:
                    best_gene = anno_res[0][0]
                    classification = 'NNC'
                    subtype = 'EXTEND_EXON_REGION_NNC'
            else:
                best_res = anno_res[0]
                ref_exon_num, reflength, diff_to_transcript_start, diff_to_transcript_end, exon_miss_5, exon_miss_3, diff_to_gene_start, diff_to_gene_end = best_res[5]
                best_anno_info[6],best_anno_info[7] = diff_to_gene_start, diff_to_gene_end
                flags = list(set([self.parts_anno[0][i]['exon_type'] for i in self.parts_anno[0] if self.parts_anno[0][i]['exon_type'] > 13]))
                if flags:
                    classification = 'NNC' if max(flags) > 29 else 'NIC'
                    subtype = '|'.join([self.subtypes[i] for i in flags])
                    if classification == 'NIC':
                        if anno_res[0][2]:
                            subtype += '|combination_of_known_junction'
                        else:
                            subtype += '|combination_of_known_splicesite'
                elif isinstance(diff_to_transcript_start,int):
                    best_anno_info = [ref_exon_num, reflength, diff_to_transcript_start, diff_to_transcript_end,exon_miss_5, exon_miss_3, diff_to_gene_start, diff_to_gene_end]
                    ISM_5_EXON = 1 if exon_miss_5 != 'NA' and exon_miss_5 >= SM_NC_EXON_NUM else 0
                    ISM_3_EXON = 1 if exon_miss_3 != 'NA' and exon_miss_3 >= SM_NC_EXON_NUM else 0
                    classification = 'ISM'
                    if ISM_5_EXON and ISM_3_EXON:
                        subtype = 'INTERNAL_ISM_EXON'
                    elif ISM_3_EXON:
                        subtype = '3_ISM_EXON'
                    elif ISM_5_EXON:
                        subtype = '5_ISM_EXON'
                    elif diff_to_transcript_start > SM_NC_ADD_LENGTH and diff_to_transcript_end > SM_NC_ADD_LENGTH:
                        subtype = '53_FLANK'
                    elif diff_to_transcript_start > SM_NC_ADD_LENGTH:
                        subtype = '5_FLANK'
                    elif diff_to_transcript_end > SM_NC_ADD_LENGTH:
                        subtype = '3_FLANK'
                    elif diff_to_transcript_start <= SM_NC_MINUES_LENGTH and diff_to_transcript_end <= SM_NC_MINUES_LENGTH:
                        subtype = 'INTERNAL_ISM'
                    elif diff_to_transcript_start <= SM_NC_MINUES_LENGTH:
                        subtype = '5_ISM'
                    elif diff_to_transcript_end <= SM_NC_MINUES_LENGTH:
                        subtype = '3_ISM'
                    else:
                        classification = 'FSM'
                        subtype = 'FSM'
                elif not flag:
                    classification = 'NIC'
                    if anno_res[0][2]:
                        subtype = 'combination_of_known_junction'
                    else:
                        subtype = 'combination_of_known_splicesite'
                else:
                    classification = 'BUG'
                    subtype = 'BUG'
        else:  ## flag > 0  不正常转录本
            self.classification = 'UNKNOWN'
            self.subtype = 'UNKNOWN'
            return 0
        if debug:
            print(self.id)
            print(classification,subtype)
            print(best_gene,best_transcript,self.CHROMS[0],self.STRANDS[0])
            print(self.EXONS_TYPE)
            print(self.exon_anno)
        self.classification = classification
        self.subtype = subtype
        self.best_gene = best_gene
        self.best_transcript = best_transcript
        self.best_anno_info = best_anno_info
        txt = self.sjCanonicalCalc()
        self.SJtxt += txt
        self.multiAnnot(mode=2)
        return 1

    def db_search_by_region(self,db,chrom,start,end):
        res = {}
        if chrom not in db:
            return res
        for gene in db[chrom]:
            ref_start,ref_end = db[chrom][gene][2]
            if self.check_overlap([start,end],[min(ref_start,ref_end),max(ref_start,ref_end)]):
                res[gene] = db[chrom][gene]
        return res

    def multiAnnot(self,gene='',res=[],mode=1):
        """
        TransAnnot V0.0.2新增外挂函数，用于处理FSM/ISM注释时可能出现的多注释情况。
        为了避免频繁内嵌莫名其妙代码的情况（虽然已经嵌入很多了），特此将实现本功能的函数全部集中于此，通过mode=123隔开。
        mode1: 通过一个类内的全局变量self.multiAnno标识多注释（可能有BUG）
        mode2: 通过最终分类识别FSM/ISM来唤醒多注释，并生成txt
        """
        if mode == 1:
            for b in res:
                self.multiAnno[b[0]] = [gene, b[0],b[1],b[2]]
        if mode == 2:
            score = self.Zsoftmax()
            mm = ''
            if self.classification in ['FSM','ISM'] and self.multiAnno:
                m = []
                #for t in sorted(self.multiAnno,key=lambda x:(self.multiAnno[x][2],self.multiAnno[x][3])):
                #    m.append('{},{},{},{}'.format(self.multiAnno[t][0],self.multiAnno[t][1],self.multiAnno[t][2],self.multiAnno[t][3]))
                for t in score:
                    m.append('{}:{:.2%}'.format(t,score[t]))
                mm = '|'.join(m)
            self.multiAnnoTXT += '{}\t{}\t{}\t{}\n'.format(self.id,self.best_gene,self.best_transcript,mm)

    def Zsoftmax(self):
        keys = []
        exons = []
        score = []
        if not self.multiAnno:
            return {}
        for i in sorted(self.multiAnno,key=lambda x:(self.multiAnno[x][2],self.multiAnno[x][3])):
            keys.append(i)
            exons.append(self.multiAnno[i][2]+1)
            score.append(self.multiAnno[i][3]+1)
        if len(keys) == 1:
            return {keys[0]:1}
        new_score = []
        for idx in range(len(exons)):
            if not idx:
                new_score.append(score[idx])
                sTemp1 = score[idx]
                sTemp2 = 0
                eTemp1 = exons[idx]
                continue
            if exons[idx] != eTemp1:
                sTemp2 += sTemp1
                sTemp1 = score[idx]
                eTemp1 = exons[idx]
                s = score[idx] + sTemp2
                new_score.append(s)
            else:
                new_score.append(score[idx] + sTemp2)
                sTemp1 = score[idx]
        ave = sum(new_score) / len(new_score)
        a = [2 * ave - i for i in new_score]
        d = [i / sum(a) for i in a]
        res = {keys[idx]:d[idx] for idx in range(len(keys))}
        return res

    def sjCanonicalCalc(self):
        res = ''
        knwonDict = {1:'known',0:'novel',-1:'NA'}
        canonicalDict = {1:'canonical',0:'non-canonical',-1:'NA'}
        #start,end = self.REF_START,self.REF_END
        start,end = self.correct_REF_START,self.correct_REF_END
        if not self.CHROMS or len(self.CHROMS) == 1: return ''
        leftMagicNum = [31,35,40,41,43,60,200]
        rightMagicNum = [32,36,40,41,42,60,200]
        for i in range(len(self.CHROMS)):
            if not i:
                continue
            startKnown,endKnown = 0,0
            chr1,chr2 = self.CHROMS[i-1],self.CHROMS[i]
            strand1,strand2 = self.STRANDS[i-1],self.STRANDS[i]
            type1,type2 = self.EXONS_TYPE[i-1],self.EXONS_TYPE[i]
            SJknown = self.SJknown[i-1]
            if chr1 != chr2 or strand1 != strand2:
                return ''
            if strand1 in ['1',1,'+']:
                site = [end[i-1] if end[i-1] not in [-1] else self.REF_END[i-1],start[i] if start[i] not in [-1] else self.REF_START[i]]
                s,e = site[0],site[1]
                if SJknown == 1:
                    startKnown,endKnown = 1,1
                elif SJknown == 0:
                    if type1 not in rightMagicNum:
                        startKnown = 1
                    if type2 not in leftMagicNum:
                        endKnown = 1
                else:
                    startKnown,endKnown = -1,-1
            else:
                site = [end[i] if end[i] not in [-1] else self.REF_END[i],start[i-1] if start[i-1] not in [-1] else self.REF_START[i-1]]
                s,e = site[1],site[0]
                if SJknown == 1:
                    startKnown, endKnown = 1, 1
                elif SJknown == 0:
                    if type1 not in leftMagicNum:
                        startKnown = 1
                    if type2 not in rightMagicNum:
                        endKnown = 1
                else:
                    startKnown, endKnown = -1, -1
            if SJknown == -1:
                startCanonical, endCanonical, SJcanonical, JunctionSeq = -1,-1,-1,'NA'
            else:
                startCanonical, endCanonical, SJcanonical, JunctionSeq = self.canonical(chr1,site,strand1)
            txt = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(self.id,self.best_gene,self.best_transcript,chr1,strand1,s,e,JunctionSeq,canonicalDict[startCanonical],canonicalDict[endCanonical],canonicalDict[SJcanonical],knwonDict[startKnown],knwonDict[endKnown],knwonDict[SJknown])
            res += txt
        return res

    def canonical(self,chrom,site,strand):
        canonicalList = ['GTAG','GCAG','ATAC']
        donorList = ['GT','GC','AT']
        acceptorList = ['AG','AC']
        left,right = site
        left2 = [left+1,left+2] if strand in ['1',1,'+'] else [left+1,left+2]
        right2 = [right-2,right-1] if strand in ['1',1,'+'] else [right-2,right-1]
        leftSeq = utils.refseq_extract(chrom, left2, strand, self.fabuffer, self.faidxDict)
        rightSeq = utils.refseq_extract(chrom, right2, strand, self.fabuffer, self.faidxDict)
        JunctionSeq = leftSeq + rightSeq if strand in ['1',1,'+'] else rightSeq+leftSeq
        if strand in ['1', 1, '+']:
            leftFlag = 1 if leftSeq in donorList else 0
            rightFlag = 1 if rightSeq in acceptorList else 0
            junctionFlag = 1 if JunctionSeq in canonicalList else 0
            startFlag, endFlag = leftFlag, rightFlag
        else:
            leftFlag = 1 if leftSeq in acceptorList else 0
            rightFlag = 1 if rightSeq in donorList else 0
            junctionFlag = 1 if JunctionSeq in canonicalList else 0
            startFlag, endFlag = rightFlag, leftFlag
        return startFlag,endFlag,junctionFlag,JunctionSeq

    def check_overlap(self,ref,gene,length=0):
        flag = 0
        overlap_length = 0
        if ref[0] <= gene[0]:
            if ref[1] <= gene[0]:
                flag = 0  # 不挨着ref在左
            elif ref[1] >= gene[1]:
                flag = 2  # ref 包含gene
                overlap_length = abs(gene[1]-gene[0]+1)
            else:
                flag = 1  # overlap
                overlap_length = abs(ref[1] - gene[0] + 1)
        else:
            if ref[0] >= gene[1]:
                flag = 0  # 不挨着ref在右
            elif ref[1] <= gene[1]:
                flag = 2  # gene包含ref
                overlap_length = abs(ref[1]-ref[0]+1)
            else:
                flag = 1  # overlap
                overlap_length = abs(gene[1]-ref[0]+1)
        if length:
            return overlap_length
        return flag

    def _get_parts_by_idx(self,parts,idx):
        return [parts[i] for i in idx]

    def list_intersection(self,list_test,list_base):
        """检查目标与注释外显子列表的连续一致性"""
        a = list(set(list_base).intersection(set(list_test)))
        try:
            if sorted(a) == sorted(list_test):
                if abs(list_base.index(list_test[-1]) - list_base.index(list_test[0]))+1 == len(a):
                    return 1
        except:
            pass
        return 0

    def loc_format(self):
        res = []
        for i in range(len(self.CHROMS)):
            res.append([self.REF_START[i],self.REF_END[i]])
        return res

    def bed_format(self,t=0):
        txt = ''
        for x,i in enumerate(self.CHROMS):
            if t == 1:
                if self.EXONS_TYPE[x] not in [1]:
                    s1,s2 = 0,0
                else:
                    s1 = (self.correct_REF_START[x] - self.REF_START[x]) if self.correct_REF_START[x] > 0 else 0
                    s2 = (self.correct_REF_END[x] - self.REF_END[x]) if self.correct_REF_END[x] > 0 else 0
                txt += '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(self.CHROMS[x], self.REF_START[x], self.REF_END[x],self.id, '*', self.STRANDS[x], self.SEQ_START[x],self.SEQ_END[x], s1,s2,self.exon_anno[x])
            elif t == 10:  # 矫正断点
                if not x or x+1 == len(self.CHROMS):
                    s1,s2 = self.REF_START[x],self.REF_END[x]
                else:
                    s1 = self.correct_REF_START[x] if self.correct_REF_START[x] > 0 else self.REF_START[x]
                    s2 = self.correct_REF_END[x] if self.correct_REF_END[x] > 0 else self.REF_END[x]
                txt += '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(self.CHROMS[x], s1, s2,self.id, '*', self.STRANDS[x], self.SEQ_START[x],self.SEQ_END[x], self.exon_anno[x])
            else:
                txt += '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(self.CHROMS[x],self.REF_START[x],self.REF_END[x],self.id,'*',self.STRANDS[x],self.SEQ_START[x],self.SEQ_END[x],self.exon_anno[x])
        return txt

    def stat_format(self,exp=0):
        ref_exon_num, reflength, diff_to_transcript_start, diff_to_transcript_end, exon_miss_5, exon_miss_3, diff_to_gene_start, diff_to_gene_end = self.best_anno_info
        if not exp:
            txt = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(self.id,self.classification,self.subtype,self.best_gene,self.best_transcript,self.anno_chrom,self.anno_strand,self.isoform_length,len(self.SEQ_START),reflength,ref_exon_num,diff_to_gene_start,diff_to_gene_end,diff_to_transcript_start,diff_to_transcript_end,exon_miss_5,exon_miss_3)
        else:
            txt = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(self.id,self.classification,self.subtype,self.best_gene,self.best_transcript,self.anno_chrom,self.anno_strand,self.isoform_length,len(self.SEQ_START),reflength,ref_exon_num,diff_to_gene_start,diff_to_gene_end,diff_to_transcript_start,diff_to_transcript_end,exon_miss_5,exon_miss_3,exp[0],exp[1])
        return txt


if __name__ == '__main__':
    import utils
    bed = r'D:\zcs-genex\SCRIPTS\ISOzcs\workflow\20200723\759133C.minimap2.bed'
    #bed = r'D:\zcs-genex\SCRIPTS\ISOzcs\workflow\20200723\test1.txt'
    dbFile = r'D:\zcs-genex\SCRIPTS\ISOzcs\script_v191115\dev_test_0325\hg38.ensembl.v20200306.1.pickle'
    fa = r'D:\zcs-genex\SCRIPTS\Scripts_for_work\lijuan_refseq\GRCh38_latest_genomic.fna.filter.fa'
    fai = r'D:\zcs-genex\SCRIPTS\Scripts_for_work\lijuan_refseq\GRCh38_latest_genomic.fna.filter.fa.fai'
    fabuffer = open(fa,'r')
    faidx_dict = utils.refseqIdx(fai)
    db = utils.load_pickle(dbFile)
    d = utils.bed_read(bed)
    for id in d:
        Cell = Isoform_anno(id, d[id],fabuffer,faidx_dict)
        flag = Cell._break_point_cluster()
        Cell.annotation(db, flag,debug=0)
        #print(Cell.SJtxt)
        print(Cell.multiAnnoTXT)
    fabuffer.close()