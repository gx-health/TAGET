# -*- coding: utf-8 -*-
"""
Author: Zhang Chengsheng, @2020.03.16
临时性脚本，用于从二代三代独立注释结果中选出最优解。是二三代混装注释最终版的临时替代品。
不排除后期跳票扶正的可能性
"""

import os,sys,pickle
import ISO_anno_V02
import transcript_cluster


sys.setrecursionlimit(100000000)


def pickle_make(dict_in,file_out):
    with open(file_out, 'wb') as o:
        pickle.dump(dict_in, o)


def main(bed_long,bed_short,output):
    long = ISO_anno_V02.bed_read(bed_long) if bed_long else {}
    short = ISO_anno_V02.bed_read(bed_short) if bed_short else {}
    bed_out = output + '.bed'
    stat_out = output + '.stat'
    classification_q = {
                        'FSM':0,
                        'ISM':1,
                        'NIC':2,
                        'NNC':3,
                        'FUSION': 4,
                        'Genic':5,
                        'Intergenic':6,
                        'UNKNOWN':10,
                        }
    StructureDB = transcript_cluster.structure_save()
    Anno_DB = transcript_cluster.structure_save()
    Gene2transciptDB = transcript_cluster.structure_save()
    ClassificationDB = transcript_cluster.structure_save()
    ClusterSM = transcript_cluster.tcluster()
    ClusterNC = transcript_cluster.tcluster()
    ClusterGENIC = transcript_cluster.tcluster()
    ClusterFUSION = transcript_cluster.tcluster()

    def add_CLUSTER(A):
        StructureDB.saveReadsinfo(A)
        ClassificationDB.saveClassification(A)
        Anno_DB.saveAnnoDB(A)
        if A.classification in ['FSM', 'ISM']:
            ClusterSM.add_SM_node(A.CHROMS[0], A.best_gene, A.best_transcript, id)
            Gene2transciptDB.saveGene2transcript(A)
        elif A.classification in ['NIC', 'NNC']:
            ClusterNC.add_NSM_node(A.CHROMS[0], [min(A.REF_START), max(A.REF_END)], A.loc_format(), id, A.best_gene,A.classification,node=0)
        elif A.classification in ['Genic', 'Intergenic']:
            ClusterGENIC.add_NSM_node(A.CHROMS[0], [min(A.REF_START), max(A.REF_END)], A.loc_format(), id, A.best_gene,A.classification,node=0)
        elif A.classification in ['FUSION']:
            ClusterFUSION.add_FUSION_node(A)

    with open(bed_out,'w') as bed, open(stat_out,'w') as stat:
        stat.write('ID\tType\tSubtype\tGene\tTranscript\tChrom\tStrand\tSeq_length\tSeq_exon_num\tRef_length\tRef_exon_num\tdiff_to_gene_start\tdiff_to_gene_end\tdiff_to_transcript_start\tdiff_to_transcript_end\n')
        count = 0
        t = 10
        for id in sorted(list(set(list(long) + list(short)))):
            count += 1
            print('\r{}\t{}'.format(count,id),end='',flush=1),
            if id in long:
                A = ISO_anno_V02.Isoform_bed_read(id, long[id])
                A1 = A._break_point_cluster()
                A.annotation(db, A1)
                classificationA = A.classification
                _23A = 1 if [i for i in A.EXONS_TYPE if i < 40] else 0
            else:
                classificationA = 0
                _23A = 0
            if id in short:
                B = ISO_anno_V02.Isoform_bed_read(id, short[id])
                B1 = B._break_point_cluster()
                B.annotation(db, B1)
                classificationB = B.classification
                _23B = 1 if [i for i in B.EXONS_TYPE if i < 40] else 0
            else:
                classificationB = 0
                _23B = 0
            if classificationA and classificationB:
                fA = 8 if classification_q[classificationA] == 3 and not _23A else classification_q[classificationA]
                fB = 8 if classification_q[classificationB] == 3 and not _23B else classification_q[classificationB]
                C,C1,_23 = [A,B,2] if fA < fB else [B,A,3]
                add_CLUSTER(C)
                txt1 = C.stat_format()
                txt2 = C.bed_format(t=t)
                bed.write(txt2)
                if classification_q[classificationA] != classification_q[classificationB]:
                    stat.write(txt1.strip('\n')+'\t{}\t{}\t{}\n'.format(_23,C1.classification,C1.subtype))
                else:
                    stat.write(txt1)
            elif classificationA:
                add_CLUSTER(A)
                txt1 = A.stat_format()
                txt2 = A.bed_format(t=t)
                bed.write(txt2)
                stat.write(txt1)
            elif classificationB:
                add_CLUSTER(B)
                txt1 = B.stat_format()
                txt2 = B.bed_format(t=t)
                bed.write(txt2)
                stat.write(txt1)
            else:
                #info = [0,0,0,0,0,0,0,0,0,0,0,0,0]
                txt1 = B.stat_format()
                txt2 = B.bed_format(t=t)
                bed.write(txt2)
                stat.write(txt1)
    txt3 = ClusterFUSION.tree_FUSION_view()
    ClusterFUSION.txt_write(txt3,output+'.FUSION.cluster')
    txt0 = ClusterSM.tree_SM_view()
    ClusterSM.txt_write(txt0,output+'.SM.cluster')
    txt1 = ClusterNC.tree_NSM_view(ClusterNC.tree)
    ClusterNC.txt_write(txt1,output+'.NC.cluster')
    txt2 = ClusterGENIC.tree_NSM_view(ClusterGENIC.tree)
    ClusterGENIC.txt_write(txt2, output + '.GENIC.cluster')
    Gene2transciptDB.addGene2transcript(ClusterNC.NSM)
    Gene2transciptDB.addGene2transcript(ClusterGENIC.NSM)
    PaintDB = {'reads':StructureDB.d,'gene2transcript':Gene2transciptDB.d,'annotation':Anno_DB.d,'classification':ClassificationDB.d,'fusion':ClusterFUSION.FUSION}
    pickle_make(PaintDB,output+'.db.pickle')


def option(argv):
    from argparse import ArgumentParser as AP
    usages = "python3 {} -d db_in -2 short_reads_bed -3 long_reads_bed -o output".format(argv[0])
    p = AP(usage=usages)
    p.add_argument("-d", dest="db_in", metavar="sam_dict", help="gtf_dict file.")
    p.add_argument("-2", dest="short", metavar="short_reads_bed", help="short reads bed format file",default=0)
    p.add_argument("-3", dest="long", metavar="long_reads_bed", help="long reads bed format file",default=0)
    p.add_argument("-o", dest="out_preffix", metavar="output_file", help="bed&stat file out.")
    p.add_argument("-p", dest="process", metavar="[int]", help="process, default: 1", type=int,default=1)

    if len(argv) == 1:
        p.print_help()
        exit(1)
    return p.parse_args(argv[1:])


if __name__ == '__main__':
    #gtf = r'/public/source/share/zcs/data/GTF/hg38.ensembl.v20200306.1.pickle'
    args = option(sys.argv)
    gtf = args.db_in
    db = ISO_anno_V02.load_pickle(gtf)
    bed_long = args.long
    bed_short = args.short
    output = args.out_preffix
    main(bed_long, bed_short, output)
