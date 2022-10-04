"""
Author: Zhang Chengsheng
"""

import utils
import numpy


class Isoform_anno:
    def __init__(self, id, string, fabuffer, faidxDict):
        self.id = id
        self.string = string
        self.fabuffer = fabuffer
        self.faidxDict = faidxDict
        self.chrom_mask = []  # 多染色体映射时因长度过短被屏蔽掉的染色体存放于此列表中
        self._origin_string_parse()
        self._bias_calc()
        self._pure_info()
        self.merge = 0
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
                         50: '边缘双端不对齐',  # 已经没有这个了
                         60: '双端新剪切点基因外无注释的新外显子',
                         100: 'BUG_100',
                         101: 'BUG_101',
                         200: '无注释',
                         400: 'BUG_NA'
                         }

    def _origin_string_parse(self):
        """解析原始BED文件，返回原始分段信息"""
        mask = True  # 是否做0mask处理的开关
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
            if mask:
                line = self._0_mask(line)
            s = line.rstrip().split('\t')
            if s[0] == '0':  # 0:unmapped,-1:masked
                self.unmapped_length += abs(int(s[7]) - int(s[6]))
            elif s[0] == '-1':
                self.mask_length += abs(int(s[7]) - int(s[6]))
            self.chroms.append(s[0])
            self.ref_start.append(int(s[1]))
            self.ref_end.append(int(s[2]))
            self.strands.append(s[5])
            self.seq_start.append(int(s[6]))
            self.seq_end.append(int(s[7]))
            self.exons_length.append(int(s[7]) - int(s[6]) + 1)
            self.bias_start.append(0)
            self.bias_end.append(0)
        self.isoformLength = sum(self.exons_length)

    def _0_mask(self,line):
        """对单外显子短于阈值的分段做0mask处理"""
        MIN_EXON_LENGHT = 10
        s = line.rstrip().split('\t')
        if abs(int(s[7]) - int(s[6])) <= MIN_EXON_LENGHT:
            return '-1\t0\t0\t{}\t*\t0\t{}\t{}'.format(s[3], s[6], s[7])
        else:
            return line

    def _bias_define(self):
        """根据exons mask情况确定SJ bias的范围"""
        for idx, i in enumerate(self.exons_length):
            if self.chroms[idx] == '-1':
                if idx - 1 >= 0:
                    if self.chroms[idx - 1] not in ['-1', '0']:
                        if self.strands[idx - 1] in [1,'1','+']:
                            self.bias_end[idx - 1] = i
                        else:
                            self.bias_start[idx - 1] = i
                if idx + 1 < len(self.exons_length):
                    if self.chroms[idx + 1] not in ['-1', '0']:
                        if self.strands[idx + 1] in [1,'1','+']:
                            self.bias_start[idx + 1] = i
                        else:
                            self.bias_end[idx + 1] = i

    def _bias_calc(self):
        """初始化外显子SJ的可信区间，留坑而已，目前没啥用"""
        MAX_NATURAL_EXONS_BIAS_THRESHOLD = 0  # 剪切位点误差范围
        self.bias_start_region = []
        self.bias_end_region = []
        for idx, i in enumerate(self.bias_start):
            self.bias_start_region.append([self.ref_start[idx]-i-MAX_NATURAL_EXONS_BIAS_THRESHOLD,self.ref_start[idx]+MAX_NATURAL_EXONS_BIAS_THRESHOLD])
        for idx, i in enumerate(self.bias_end):
            self.bias_end_region.append([self.ref_end[idx]-MAX_NATURAL_EXONS_BIAS_THRESHOLD,self.ref_end[idx]+i+MAX_NATURAL_EXONS_BIAS_THRESHOLD])

    def _pure_info(self, chrom_mask=[]):
        MIN_INTRON = 50  #最小intron长度，低于此长度的intron将导致其前后exon被此段代码合并。
        self.CHROMS = [self.chroms[i] for i in range(len(self.chroms)) if (self.chroms[i] not in ['0', '-1'] and self.chroms[i] not in chrom_mask)]
        self.REF_START = [self.ref_start[i] for i in range(len(self.chroms)) if (self.chroms[i] not in ['0', '-1'] and self.chroms[i] not in chrom_mask)]
        self.REF_END = [self.ref_end[i] for i in range(len(self.chroms)) if (self.chroms[i] not in ['0', '-1'] and self.chroms[i] not in chrom_mask)]
        self.SEQ_START = [self.seq_start[i] for i in range(len(self.chroms)) if (self.chroms[i] not in ['0', '-1'] and self.chroms[i] not in chrom_mask)]
        self.SEQ_END = [self.seq_end[i] for i in range(len(self.chroms)) if (self.chroms[i] not in ['0', '-1'] and self.chroms[i] not in chrom_mask)]
        self.STRANDS = [self.strands[i] for i in range(len(self.chroms)) if (self.chroms[i] not in ['0', '-1'] and self.chroms[i] not in chrom_mask)]
        self.EXONS_LENGTH = [self.exons_length[i] for i in range(len(self.chroms)) if (self.chroms[i] not in ['0', '-1'] and self.chroms[i] not in chrom_mask)]
        self.BIAS_START = [self.bias_start[i] for i in range(len(self.chroms)) if (self.chroms[i] not in ['0', '-1'] and self.chroms[i] not in chrom_mask)]
        self.BIAS_END = [self.bias_end[i] for i in range(len(self.chroms)) if (self.chroms[i] not in ['0', '-1'] and self.chroms[i] not in chrom_mask)]
        self.BIAS_START_REGION = [self.bias_start_region[i] for i in range(len(self.chroms)) if (self.chroms[i] not in ['0', '-1'] and self.chroms[i] not in chrom_mask)]
        self.BIAS_END_REGION = [self.bias_end_region[i] for i in range(len(self.chroms)) if (self.chroms[i] not in ['0', '-1'] and self.chroms[i] not in chrom_mask)]
        delete_idx = []
        for idx, i in enumerate(self.CHROMS):
            if not idx: continue
            if self.CHROMS[idx] == self.CHROMS[idx - 1] and self.STRANDS[idx] == self.STRANDS[idx - 1]:
                strand = self.STRANDS[idx]
                if strand in [-1, '-1','-']:
                    gap1 = abs(self.REF_START[idx - 1] - self.REF_END[idx])
                    gap2 = abs(gap1 - abs(self.SEQ_START[idx] - self.SEQ_END[idx - 1]))
                else:
                    gap1 = abs(self.REF_START[idx] - self.REF_END[idx - 1])
                    gap2 = abs(gap1 - abs(self.SEQ_START[idx] - self.SEQ_END[idx - 1]))
                if gap1 <= MIN_INTRON or gap2 <= MIN_INTRON:
                    new_ref_start = min(self.REF_START[idx], self.REF_START[idx - 1], self.REF_END[idx],self.REF_END[idx - 1])
                    new_ref_end = max(self.REF_START[idx], self.REF_START[idx - 1], self.REF_END[idx],self.REF_END[idx - 1])
                    new_seq_start = min(self.SEQ_START[idx], self.SEQ_START[idx - 1], self.SEQ_END[idx], self.SEQ_END[idx - 1])
                    new_seq_end = max(self.SEQ_START[idx], self.SEQ_START[idx - 1], self.SEQ_END[idx],self.SEQ_END[idx - 1])
                    self.REF_START[idx - 1] = new_ref_start
                    self.REF_END[idx - 1] = new_ref_end
                    self.SEQ_START[idx - 1] = new_seq_start
                    self.SEQ_END[idx - 1] = new_seq_end
                    if new_ref_start == self.REF_START[idx]:
                        self.BIAS_START[idx - 1] = self.BIAS_START[idx]
                        self.BIAS_START_REGION[idx - 1] = self.BIAS_START_REGION[idx]
                    else:
                        self.BIAS_END[idx - 1] = self.BIAS_END[idx]
                        self.BIAS_END_REGION[idx - 1] = self.BIAS_END_REGION[idx]
                    delete_idx.append(idx)
        for idx in sorted(delete_idx, reverse=True): self._info_refresh(idx)
        if delete_idx: self.merge = 1
        #####
        delete_idx = []
        for idx, i in enumerate(self.CHROMS):
            if not idx or idx + 1 == len(self.CHROMS):
                continue
            idx1, idx2 = idx - 1, idx + 1
            if self.CHROMS[idx1] == self.CHROMS[idx2] and self.STRANDS[idx1] == self.STRANDS[idx2]:
                strand = self.STRANDS[idx1]
                if strand in [1, '1']:
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
        for idx in sorted(list(set(delete_idx)), reverse=True):
            self._info_refresh(idx)
        if delete_idx: self.merge = 1

    def _info_refresh(self, idx):
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

    def _multi_valid_chroms_check(self):
        """检查是否存在多染色体映射，同时检查其长度 将低于阈值的分段移除"""
        if len(set(self.CHROMS)) == 1:
            return 1
        FUSION_MIN_LENGTH_THRESHOLD = 200
        self.CHROMS_CLUSTER = {}
        for _idx, chrom in enumerate(self.CHROMS):
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

    def _calc_length(self, dict_in):
        length = 0
        for i in dict_in:
            length += abs(dict_in[i][-1] - dict_in[i][0]) + 1
        return length

    def _gap_check(self):
        start = self.REF_START
        end = self.REF_END
        strand = self.STRANDS
        chrom = self.CHROMS
        MAX_GAP_LENGTH = 1000000
        MIN_ISOFORM_LENGTH = 200
        res_idx = {0: []}
        indicator = 0
        for _idx, i in enumerate(start):
            if not i:
                pass
            if not _idx:
                res_idx[indicator].append(_idx)
                continue
            gap1 = abs(start[_idx] - start[_idx - 1])
            gap2 = abs(end[_idx] - end[_idx - 1])
            if gap1 > MAX_GAP_LENGTH or gap2 > MAX_GAP_LENGTH or strand[_idx] != strand[_idx - 1] or chrom[_idx] != chrom[_idx - 1]:
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
        for i in sorted(idx_list, reverse=1):
            self._info_refresh(i)
        if idx_list:
            gap_num, res_idx = self._gap_check()
        else:
            gap_num = len(res_idx) - 1
        return gap_num, res_idx


def debug():
    bed = r'D:\zcs-genex\SCRIPTS\ISOzcs\workflow\20200723\759133C.minimap2.bed'
    d = utils.bed_read(bed)
    fabuffer = ''
    faidx_dict = ''
    for id in d:
        Cell = Isoform_anno(id, d[id], fabuffer, faidx_dict)
        exit(111)


if __name__ == '__main__':
    debug()
