"""
python3 script
Author: Zhang Chengsheng, @2019.07.01
dev 200310
"""

import sys
import re
import multiprocessing
import utils


def sam_read(sam_in,type='long'):
    samDict = {}
    with open(sam_in,'r') as f:
        for i in f.readlines():
            if i.startswith('@'):
                continue
            line = i.strip().split('\t')
            if type == 'long':
                id = line[0]
            else:
                id = '_'.join(line[0].split('_')[:-1])
            if id not in samDict:
                samDict[id] = [i.strip('\n')]
            else:
                samDict[id].append(i.strip('\n'))
    return samDict


def long_sam_parse(cluster,c=0,report=0):
    if c % 100 == 1 and report:
        print('\r{} isoforms finished !'.format(c-1),end='',flush=True)
    res = {}
    for i in cluster:
        line = i.strip().split('\t')
        read = line
        loc, strand1, chrom1 = read_location(read,0)
        for j in loc:
            key = '|'.join([str(i) for i in loc[j]])
            if j not in res:
                if loc[j] == [0,0,0]:
                    res[j] = {}
                    continue
                else:
                    res[j] = {key:1}
            else:
                if loc[j] == [0,0,0]:
                    continue
                if key not in res[j]:
                    res[j][key] = 1
                else:
                    res[j][key] += 1
    return res


def short_sam_parse(cluster,c=0,OVERLAP_LENGTH=80,READ_LENGTH=100,report=0):
    if c % 100 == 1 and report:
        print('\r{} read finished...'.format(c),end='',flush=True)
    dict_cluster = short_sam_cluster_trans(cluster)
    res = {}
    for i in dict_cluster:
        for ii in dict_cluster[i]:
            loc, strand1, chrom1 = read_location(ii, i,OVERLAP_LENGTH,READ_LENGTH)
            for j in loc:
                key = '|'.join([str(i) for i in loc[j]])
                if j not in res:
                    if loc[j] == [0, 0, 0]:
                        res[j] = {}
                        continue
                    else:
                        res[j] = {key: 1}
                else:
                    if loc[j] == [0, 0, 0]:
                        continue
                    if key not in res[j]:
                        res[j][key] = 1
                    else:
                        res[j][key] += 1
    return res


def short_sam_cluster_trans(cluster):
    cluster = sorted(cluster, key=lambda x: int(x.split('\t')[0].split('_')[-1]))
    initial = True
    dictCluster = {}
    for i in cluster:
        line = i.strip().split('\t')
        tid = int(line[0].split('_')[-1])
        if initial:
            initial = False
            ttid = tid
            reads = [line]
            continue
        if tid == ttid:
            ttid = tid
            reads.append(line)
        else:
            dictCluster[ttid] = reads
            ttid = tid
            reads = [line]
    else:
        dictCluster[ttid] = reads
    return dictCluster


def read_location(line,tidx,OVERLAP_LENGTH = 80,READ_LENGTH = 100):
    CIGAR, start, strandString,chrom = line[5],int(line[3]),line[1],line[2]
    start = int(start)
    strand = (-1 if bin(int(strandString))[-5] == '1' else 1) if len(bin(int(strandString))) > 6 else 1
    miss_match = (0 if bin(int(strandString))[-3] == '1' else 1) if len(bin(int(strandString))) > 4 else 1
    res = {}
    c = re.findall(re.compile('\D'), CIGAR)
    n = re.findall(re.compile('\d+'), CIGAR)
    cn = [[c[i], n[i]] for i in range(len(c))] if CIGAR not in ['*',''] else []

    tLen = 0
    if not cn:
        tLen = len(line[9])
    else:
        for i, j in cn:
            if i not in ['D', 'N']:
                tLen += int(j)

    site = [(READ_LENGTH-OVERLAP_LENGTH) * tidx + i + 1 for i in range(tLen)] if strand == 1 else [(READ_LENGTH-OVERLAP_LENGTH) * tidx + tLen - i for i in range(tLen)]
    seq_idx, ref_idx = 0, 0
    if not miss_match:
        for i in range(tLen):
            res[site[i]] = [0,0,0]
        return res,0,'*'

    for i,j in cn:
        if i in ['M','=']:
            for k in range(int(j)):
                res[site[seq_idx]] = [start + ref_idx,strand,chrom]
                seq_idx += 1
                ref_idx += 1
        if i in ['D','N']:
            for k in range(int(j)):
                ref_idx += 1
        if i == 'I':
            for k in range(int(j)):
                res[site[seq_idx]] = [-1,0,0]  # start + ref_idx + 0.5 insertion
                seq_idx += 1
        if i in ['S','H']:
            for k in range(int(j)):
                res[site[seq_idx]] = [0,0,0]  # deletion
                seq_idx += 1
    return res,strand,chrom


def loc_parse(loc,oped=False):
    loc_idx = sorted(loc)
    start, end = loc[loc_idx[0]][0], loc[loc_idx[-1]][0]
    if oped:
        return start,end
    MAX_GAP = 10
    loc_seq,loc_ref = [],[]
    loc_strand,loc_chrom = [],[]
    for idx,i in enumerate(loc_idx):
        if not idx:
            temp = i
            miss_temp = 0
        else:
            if loc[loc_idx[idx]][0] == -1:  # or loc[loc_idx[idx-1]][0] == -1:
                miss_temp += 1
                continue
            elif abs(loc[loc_idx[idx]][0] - loc[loc_idx[idx-1-miss_temp]][0]) > MAX_GAP or loc[loc_idx[idx]][1] != loc[loc_idx[idx-1-miss_temp]][1] or loc[loc_idx[idx]][2] != loc[loc_idx[idx-1-miss_temp]][2]:
                loc_seq.append([temp,i-1-miss_temp])
                loc_ref.append([loc[temp][0],loc[loc_idx[idx-1-miss_temp]][0]])
                loc_strand.append(loc[temp][1])
                loc_chrom.append(loc[temp][2])
                if miss_temp:
                    loc_seq.append([i-miss_temp,i-1])
                    loc_ref.append([0,0])
                    loc_strand.append(0)
                    loc_chrom.append(0)
                temp = i
                miss_temp = 0
            else:
                miss_temp = 0
    else:
        loc_seq.append([temp, i])
        loc_ref.append([loc[temp][0], loc[i][0]])
        loc_strand.append(loc[temp][1])
        loc_chrom.append(loc[temp][2])
    return loc_seq,loc_ref,loc_strand,loc_chrom


def res_combine(short,long,debug=0):
    weight_long = 5
    if not long:
        return short
    else:
        for idx,i in enumerate(sorted(short)):
            if i in long and long[i]:
                if debug:
                    print(i,long[i],'-----||||-----',short[i])
                for j in long[i]:
                    if j in short[i]:
                        short[i][j] += weight_long*long[i][j]
                    else:
                        short[i][j] = weight_long*long[i][j]
    return short


def loc_vote(res):
    def get_dict_max(idx):
        dict = res[idx]
        base_weight = 10  # 基础权重
        continus_weight = 1  # 连续性权重
        max_point,max_score = [],0
        for i in dict:
            loc,strand,chr = i.split('|')
            loc,strand,chr = int(loc),int(strand),str(chr)
            weight = base_weight * dict[i]
            if idx-1 and loc:
                if abs(loc - res[idx-1][0]) <= 2 and strand == int(res[idx-1][1]):
                    weight += 5
            if len(res)-idx-1 > 0 and loc:
                for j in res[idx+1]:
                    if abs(loc-int(j.split('|')[0])) <= 2 and strand == int(j.split('|')[1]):
                        weight += continus_weight*res[idx+1][j]
                        break
            if weight > max_score:
                max_score = weight
                max_point = [loc,strand,chr]
        else:
            pass
        return max_point

    for idx in sorted(res):
        if not res[idx]:
            res[idx] = [0,0,0]
        elif len(res[idx]) == 1:
            data = sorted(res[idx])[0].split('|')
            data[0],data[1] = int(data[0]),int(data[1])
            res[idx] = data
        else:
            data = get_dict_max(idx)  # 投票表决，考虑上下游情况，保持连续性
            res[idx] = data
    return res


def _wrap(func):
    def wrap(self,*args, **kewargs):
        try:
            res = func(self)
            return res
        except Exception as e:
            print('ERROR: ' + self.t_id)
            print(e)
            return 404
    return wrap


class vote:
    def __init__(self,res,t_id):
        self.res = res
        self.t_id = t_id

    @_wrap
    def travers(self):
        res_length = max(self.res)
        region = self._find_multi_region()
        for idx,i in enumerate(region[2]):
            res = self.get_res_region(i)
            res = self.two_region_parse(res)  # 做区域分析返回解集

            if len(res) > 10:
                return 0
            left,right = max(0,i[0]-1),min(res_length,i[1]+1)
            upstream = list(self.get_res_region([left,left])[left])[0] if left else 0
            downstream = list(self.get_res_region([right,right])[right])[0] if res_length-right else 0
            res1 = self.parts_score(res, i,upstream,downstream,debug=idx)  # 解集寻路
            if not res1:
                return 0
            # 区域最优解替换原始解
            for i in res1:
                self.res[i] = res1[i]
        return self.res

    def get_res_region(self,region):
        res = {}
        for i in range(region[0],region[1]+1,1):
            res[i] = self.res[i]
        return res

    def _find_multi_region(self):
        region = {0: [], 1: [], 2: []}
        flag = 0
        temp = [0, 0]
        MIN_ONE_LENGTH_THRESHOLD = 50
        for i in sorted(self.res):
            if len(self.res[i]) == 1:
                if flag == 1:
                    pass
                else:
                    temp[-1] = int(i-1)
                    region[flag].append(temp)
                    temp = [i,0]
                    flag = 1
            else:  # len(self.res[i]) > 1:
                if flag == 2:
                    pass
                else:
                    temp[-1] = int(i - 1)
                    region[flag].append(temp)
                    temp = [i, 0]
                    flag = 2
        else:
            temp[-1] = int(i)
            region[flag].append(temp)
        #### 合并 1 中连续长度低于MIN_ONE_LENGTH_THRESHOLD 的分段
        mini_combo = {}
        mini_combo_idx = {}
        for idx,i in enumerate(region[1]):
            if abs(i[-1]-i[0]) < MIN_ONE_LENGTH_THRESHOLD:
                mini_combo[idx] = [i[0]-1,i[-1]+1]
        for i in mini_combo:
            k1,k2 = 0,0
            idx1,idx2 = 0,0
            for idx,j in enumerate(region[2]):
                if mini_combo[i][0] in j:
                    k1 = 1
                    idx1 = idx
                if mini_combo[i][1] in j:
                    k2 = 1
                    idx2 = idx
            else:
                if k1 and k2:
                    mini_combo_idx[i] = [idx1,idx2]
        if mini_combo_idx:
            for i in sorted(mini_combo_idx,key=lambda x:max(mini_combo_idx[x]),reverse=1):
                region[2][min(mini_combo_idx[i])] = [region[2][min(mini_combo_idx[i])][0],region[2][max(mini_combo_idx[i])][1]]
                region[2].pop(max(mini_combo_idx[i]))
                region[1].pop(i)

        return region

    def two_region_parse(self,two_region):
        MAX_GAP_THRESHOLD = 10
        res = []
        names = locals()
        roads = []
        empty = {}
        for loc in two_region:
            road_num = len(two_region[loc])
            if not road_num:
                empty[loc] = {}
            if len(roads) < road_num:  # 初始化road
                for i in range(len(roads),road_num,1):
                    names['line_%s' % str(i)] = {}
                    roads.append(names['line_%s' % str(i)])

            for idx,road in enumerate(roads):  # 更新road
                if not road:
                    continue
                r_key = list(road[loc-1])[0]
                r_ref, r_strand, r_chr = r_key.split('|')
                for key in list(two_region[loc]):
                    ref,strand,chr = key.split('|')
                    if chr != r_chr or strand != r_strand:  # road上传
                        continue
                    elif abs(int(ref)-int(r_ref)) > MAX_GAP_THRESHOLD:  # road上传
                        continue
                    else:  # road更新
                        value = two_region[loc].pop(key)
                        roads[idx][loc] = {key:value}
                        break
                else:
                    res.append(road)
                    roads[idx] = {}

            else:  # 新建road
                for key in list(two_region[loc]):
                    for idx,road in enumerate(roads):
                        if not road:
                            value = two_region[loc].pop(key)
                            roads[idx][loc] = {key:value}
                            break
        else:
            for idx,road in enumerate(roads):
                if road:
                    res.append(road)
        if empty:
            res.append(empty)
        return res

    def parts_score(self,parts,region,up,down,debug=0):
        MAX_GAP = 10

        def _get_range_and_score(idx,part):
            """
            input: parts_dict:  [loc:{ref:score} ...]
            idx,score,strand,chr,[loc_start,loc_end],[ref_start,ref_end]
            """
            loc = [min(list(part)),max(list(part))]
            score = 0
            site_min,site_max = 0,0
            for i in part:
                if not part[i]:
                    return idx,loc[-1]-loc[0]+1,0,0,loc,[0,0]
                key = list(part[i])[0]
                value = part[i][key]
                site,strand,chr = key.split('|')
                score += value
                site = int(site)
                strand = int(strand)
                site_min = site if not site_min else (site if site and site < site_min else site_min)
                site_max = site if site > site_max else site_max
            return idx,score,strand,chr,loc,[site_min,site_max]

        def find_the_path(parts_info,region):
            """
            input: sorted parts_info, region
            output: res
            """
            OVERLAP_THRESHOLD = 50
            parts_num = len(parts_info)
            parts_idx = range(parts_num)
            endless = region[-1]
            pathway = []
            res = []

            def iter_circle(end,p_idx,start=0):
                ss = 0 if start else p_idx[-1]+1
                temp_p_idx = p_idx.copy()
                if 0 < endless - end < OVERLAP_THRESHOLD:
                    if temp_p_idx not in res:
                        res.append(temp_p_idx)
                for idx in parts_idx[ss:]:
                    start1,end1 = parts_info[idx][4]
                    num = parts_info[idx][0]
                    if start:
                        if start1 == end:
                            if end1 == endless:
                                temp_p_idx.append(idx)
                                res.append(temp_p_idx)
                            else:
                                temp_p_idx.append(idx)
                                iter_circle(end1,temp_p_idx)
                    elif -1 <= int(end - start1) <= OVERLAP_THRESHOLD:
                        if end1 == endless:
                            temp_p_idx.append(idx)
                            res.append(temp_p_idx)
                        else:
                            temp_p_idx.append(idx)
                            iter_circle(end1,temp_p_idx)
                    else:
                        pass
                    temp_p_idx = p_idx.copy()

            iter_circle(region[0],pathway,start=1)
            return res

        def _get_path_score(path):
            """
            input: path
            output: score
            """

            def _overlap_jump_score(info,before_info):
                b_i, b_score, b_strand, b_chr, b_seq, b_ref = before_info
                a_i, a_score, a_strand, a_chr, a_seq, a_ref = info
                a_len = a_seq[-1] - a_seq[0] + 1
                b_len = b_seq[-1] - b_seq[0] + 1
                # 判断seq overlap， 打分， 选出最优
                bias = 0
                if a_seq[0] - b_seq[1] == 1:
                    pass
                elif a_seq[0] - b_seq[1] < 1:
                    overlap_len = abs(a_seq[0] - b_seq[1])+1
                    bias += min(overlap_len*0.01,0.1)
                else:
                    pass

                if a_chr != b_chr and a_chr and b_chr:
                   bias += 0.4

                a_end = a_ref[0] if a_strand else a_ref[1]
                b_end = b_ref[0] if b_strand else b_ref[1]

                if abs(a_end-b_end) > 100000 and a_end and b_end:
                    bias += 0.1
                elif a_strand != b_strand:
                    if a_chr and b_chr:
                        bias += 0.1
                elif (a_end - b_end) * a_strand < 0:
                    bias += min(abs((a_end - b_end))*0.01,0.1)
                bias = min(0.4,bias)
                return bias

            def _start_end_score(chr,strand,ref_site,type='start'):
                MAX_GAP_THRESHOLD = 100000
                multi = 0
                _site, _strand, _chr = (up.split('|') if up else [0,0,0]) if type == 'start' else (down.split('|') if down else [0,0,0])
                if chr != _chr:
                    multi += 0.2
                if abs(int(_site) - int(ref_site)) > MAX_GAP_THRESHOLD:
                    multi += 0.2
                else:
                    multi += 0.2*max(0,abs(int(_site) - int(ref_site))/MAX_GAP_THRESHOLD*0.05)
                if strand != _strand:
                    multi += 0.05
                multi = min(0.2,multi)
                return multi

            total_score = 0
            num = len(path)
            base_score = sum([i[1] for i in parts_info])
            base_length = region[-1]-region[0]+1
            base_average = int(base_score/base_length)

            part_self_score_minus = 0
            parts_link_score_minus = 0
            parts_link_score_minus_list = [0]
            parts_up_down_minus = 0
            part_number_minus = 0
            for idx,i in enumerate(path):
                info = parts_info[i]  # (7, 174, 1, '7', [2249, 2300], [75398916, 75398967])
                id = info[0]
                strand = int(info[2])
                chr = info[3]
                seq_start,seq_end = info[4]
                seq_length = seq_end - seq_start + 1
                if seq_length < 10:
                    part_self_score_minus += 0.05
                ref_start,ref_end = info[5][int((strand-1)/2)],info[5][int((-strand-1)/2)]
                part_score = info[1]  # SCORE_1
                total_score += part_score
                dict_info = parts[id]
                if not idx:
                    multi_start = _start_end_score(chr,strand,ref_start,type='start')  # 计算和up的差异评分
                if num-idx == 1:
                    multi_end = _start_end_score(chr,strand,ref_end,type='end')  # 计算和down的差异评分
                if not idx:
                    before_info = info
                    continue
                else:
                    parts_link_score_minus_list.append(_overlap_jump_score(info,before_info))
                    before_info = info
            else:

                parts_up_down_minus = (multi_start+multi_end)*base_score
                part_self_score_minus = max(parts_link_score_minus_list)*base_score
                total_score -= parts_up_down_minus
                total_score -= parts_link_score_minus
                total_score -= part_self_score_minus

            return total_score

        up_site, up_strand, up_chr = up.split('|') if up else [0, 0, 0]
        down_site, down_strand, down_chr = down.split('|') if down else [0, 0, 0]
        if up_site and down_site and up_strand == down_strand and up_chr == down_chr:
            if abs(abs(int(up_site) - int(down_site)) - abs(region[-1] - region[0])) < MAX_GAP:
                # 触发自动补全,未判断长度
                res = {}
                direction = 1 if int(down_site) > int(up_site) else -1
                for _idx, i in enumerate(range(region[0], region[1] + 1, 1)):
                    site = int(up_site) + (_idx + 1) * direction if _idx + 1 < abs(int(down_site) - int(up_site)) else int(down_site)-direction
                    res[i] = {'{}|{}|{}'.format(site, up_strand, up_chr): 1} if site else {}
                return res
        parts_info = []
        for idx,i in enumerate(parts):
            info = _get_range_and_score(idx,i)
            parts_info.append(info)
        parts_info = sorted(parts_info, key=lambda x: (x[4][0], x[4][1]))
        if len(parts_info) > 50:
            return 0
        res = find_the_path(parts_info,region)
        if not res or len(res) > 1000:
            return 0
        res_score_list = []
        for i in res:
            score = _get_path_score(i)
            res_score_list.append([score,i])
        res_score_list = sorted(res_score_list,key=lambda x:x[0],reverse=True)
        best_path = res_score_list[0][1]
        best_res = {}
        for idx in best_path:
            dict1 = parts[parts_info[idx][0]]
            for idx1 in dict1:
                best_res[idx1] = dict1[idx1]
        return best_res


def bed_write(chrom,name,seq,res,strand):
    context = ''
    if len(seq) == len(res) == len(strand):
        for i in range(len(seq)):
            start,end = min(res[i][0],res[i][1]),max(res[i][0],res[i][1])
            context += "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom[i],start,end,name,'*',strand[i],seq[i][0],seq[i][1])
    return context


def worker(cluster, idx, t_id, type, OVERLAP_LENGTH=80, READ_LENGTH=100):
    if type == 'long':
        res = long_sam_parse(cluster, idx, report=0)
    elif type == 'short':
        res = short_sam_parse(cluster, idx,OVERLAP_LENGTH,READ_LENGTH)
    else:
        print('!!!!!!!!!')
        exit(1)
    A = vote(res, t_id)
    res = A.travers()
    if res in [0,404]:
        return res
    res = loc_vote(res)  # 投票出唯一位点
    loc_seq, loc_ref, loc_strand, loc_chrom = loc_parse(res)
    context = bed_write(loc_chrom, t_id, loc_seq, loc_ref, loc_strand)
    return context


def main(sam,output,process=1,type='long',OVERLAP_LENGTH=80,READ_LENGTH=100):
    long = sam_read(sam, type=type)
    p = multiprocessing.Pool(processes=process)
    multiprocess_result = []
    o = open(output, 'w')
    #o1 = open(output + '.undo', 'w')
    done,undo = 0,0
    for idx,t_id in enumerate(long):
        cluster = long[t_id]
        if process > 1:
            re1 = p.apply_async(worker,args=(cluster,idx,t_id,type,OVERLAP_LENGTH,READ_LENGTH,))
            multiprocess_result.append([t_id,re1])
        else:
            context = worker(cluster,idx,t_id,type,OVERLAP_LENGTH,READ_LENGTH)
            if context in [0,404]:
                undo += 1
                #o1.write('{}\t{}\n'.format(t_id,context))
            else:
                done += 1
                o.write(context)
    if process > 1:
        p.close()
        p.join()
        for t_id,re1 in multiprocess_result:
            if re1.get() in [0,404]:
                undo += 1
                #o1.write('{}\t{}\n'.format(t_id,re1.get()))
            else:
                done += 1
                o.write(re1.get())
    #print('Sam2Bed: {} reads done, {} reads undo'.format(done,undo))
    #o1.close()
    o.close()


def option(argv):
    from argparse import ArgumentParser as AP
    usages = "python3 {} -i sam -o bed_out -t [type]".format(argv[0])
    p = AP(usage=usages)
    p.add_argument("-i", dest="sam", metavar="sam", help="nsorted sam file.")
    p.add_argument("-o", dest="bed_out", metavar="bed_out", help="bed file out.")
    p.add_argument("-t", dest="type", metavar="[long/short]", help="type of read, long or short. default: short", default='short')
    p.add_argument("-l", dest="read_length", metavar="[int]", help="short read length, default: 100", type=int, default=100)
    p.add_argument("-v", dest="overlap_length", metavar="[int]", help="overlab length, default: 80", type=int, default=80)
    p.add_argument("-p", dest="process", metavar="[int]", help="process, default: 1", type=int,default=1)

    if len(argv) == 1:
        p.print_help()
        exit(1)
    return p.parse_args(argv[1:])


if __name__ == '__main__':
    args = option(sys.argv)
    main(args.sam,args.bed_out,type=args.type,process=args.process,READ_LENGTH=args.read_length,OVERLAP_LENGTH=args.overlap_length)


