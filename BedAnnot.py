"""

"""

import os,sys
import AnnotMoudle
import utils
from TransAnnotMerge import TreeCluster
import gtf2db


def main(hisat2Bed,minimap2Bed,gmapBed,baseOutDir,sample,suffix1,db,fa,cnn=0,TPM=0,report=0):
    if not os.path.exists(baseOutDir):
        os.makedirs(baseOutDir)
    fabuffer = open(fa, 'r')
    faidx_dict = utils.refseqIdx(fa+'.fai')
    hisat2 = utils.bed_read(hisat2Bed) if hisat2Bed else {}
    minimap2 = utils.bed_read(minimap2Bed) if minimap2Bed else {}
    gmap = utils.bed_read(gmapBed) if gmapBed else {}
    bedOut = os.path.join(baseOutDir,sample+suffix1+'.bed')
    statOut = os.path.join(baseOutDir,sample+suffix1+'.stat')
    SJout = os.path.join(baseOutDir,sample+suffix1+'.junction')
    DIUout = os.path.join(baseOutDir,sample+suffix1+'.multiAnno')
    tpmDict = tpmRead(TPM) if TPM else {}
    classification_q = {
        'FSM': 0,
        'ISM': 1,
        'NIC': 2,
        'NNC': 3,
        'FUSION': 4,
        'Genic': 5,
        'Intergenic': 6,
        'UNKNOWN': 10,
    }
    ReadsDB = StructureDB()
    AnnotationDB = StructureDB()
    ClassificationDB = StructureDB()
    BigTree = TreeCluster(TPM=TPM)
    cnn = 1 if cnn in [1,'1',True] else 0

    def annotBlue(d):
        Cell = AnnotMoudle.Isoform_anno(id, d[id],fabuffer,faidx_dict)
        flag = Cell._break_point_cluster()
        Cell.annotation(db, flag)
        classification = Cell.classification
        _23rd = 1 if [i for i in Cell.EXONS_TYPE if i < 40] else 0
        f = 8 if classification_q[classification] == 3 and not _23rd else classification_q[classification]
        return Cell,classification,_23rd,f

    def dbRed(Annot):
        ReadsDB.saveReadsinfo(Annot,sample)
        AnnotationDB.saveAnnoDB(Annot)
        ClassificationDB.saveClassification(Annot,sample)
        exons = [[Annot.REF_START[i],Annot.REF_END[i]] for i in range(len(Annot.REF_START))]
        tpm = float(tpmDict[Annot.id][1]) if TPM and Annot.id in tpmDict else 0
        flc = int(tpmDict[Annot.id][0]) if TPM and Annot.id in tpmDict else 0
        try:
            chrom = Annot.CHROMS[0]
            BigTree.add(sample, Annot.id, chrom, Annot.best_gene, Annot.best_transcript,[min(Annot.REF_START), max(Annot.REF_END)], exons, Annot.classification, part_idx=Annot.parts_idx, tpm=tpm, flc=flc)
        except Exception as e:
            return

    with open(bedOut,'w') as bed, open(statOut,'w') as stat,open(SJout,'w') as sj,open(DIUout,'w') as DIU:
        if TPM:
            stat.write('ID\tClassification\tSubtype\tGene\tTranscript\tChrom\tStrand\tSeq_length\tSeq_exon_num\tRef_length\tRef_exon_num\tdiff_to_gene_start\tdiff_to_gene_end\tdiff_to_transcript_start\tdiff_to_transcript_end\texon_miss_to_transcript_start\texon_miss_to_transcript_end\tFLC\tTPM\n')
        else:
            stat.write('ID\tClassification\tSubtype\tGene\tTranscript\tChrom\tStrand\tSeq_length\tSeq_exon_num\tRef_length\tRef_exon_num\tdiff_to_gene_start\tdiff_to_gene_end\tdiff_to_transcript_start\tdiff_to_transcript_end\texon_miss_to_transcript_start\texon_miss_to_transcript_end\n')
        sj.write('ID\tGene\tTranscript\tChrom\tStrand\tStart\tEnd\tSequence\tStartCanonical\tEndCanonical\tJunctionCanonical\tStartKnown\tEndKnown\tJunctionKnown\n')
        DIU.write('ID\tGene\tTranscript\tAnnotation\n')
        count = 0
        for id in sorted(list(set(list(hisat2)+list(minimap2)+list(gmap)))):
            count += 1
            if report: print('\r{}\t{}  '.format(count, id), end='', flush=True),
            H, Hc, H23,Hf = annotBlue(hisat2) if id in hisat2 else [0,'UNKNOWN',0,10]
            M, Mc, M23,Mf = annotBlue(minimap2) if id in minimap2 else [0,'UNKNOWN',0,10]
            G, Gc, G23,Gf = annotBlue(gmap) if id in gmap else [0,'UNKNOWN',0,10]
            txt1,txt2,txt3,txt4 = '','','',''
            if min(Hf,Mf,Gf) == Mf:
                if not M: continue
                dbRed(M)
                txt1 = M.stat_format(exp=tpmDict[id] if id in tpmDict else 0)
                txt2 = M.bed_format(t=cnn)
                txt3 = M.SJtxt
                txt4 = M.multiAnnoTXT
            elif min(Hf,Mf,Gf) == Gf:
                if not G: continue
                dbRed(G)
                txt1 = G.stat_format(exp=tpmDict[id] if id in tpmDict else 0)
                txt2 = G.bed_format(t=cnn)
                txt3 = G.SJtxt
                txt4 = G.multiAnnoTXT
            else:
                if not H: continue
                dbRed(H)
                txt1 = H.stat_format(exp=tpmDict[id] if id in tpmDict else 0)
                txt2 = H.bed_format(t=cnn)
                txt3 = H.SJtxt
                txt4 = H.multiAnnoTXT
            bed.write(txt2)
            stat.write(txt1)
            sj.write(txt3)
            DIU.write(txt4)
        BigTree.forest()
        g1, g2, g3 = BigTree.DBView(sample,showNone=True)
        f1, f2, f3, FusionDB = BigTree.DBviewFusion(sample,showNone=True)
        out1 = os.path.join(baseOutDir,sample+suffix1+'.cluster.gene')
        out2 = os.path.join(baseOutDir,sample+suffix1+'.cluster.transcript')
        out3 = os.path.join(baseOutDir,sample+suffix1+'.cluster.reads')
        with open(out1, 'w') as o1, open(out2, 'w') as o2, open(out3, 'w') as o3:
            o1.write(g1)
            o2.write(g2)
            o3.write(g3)
            o1.write(f1)
            o2.write(f2)
            o3.write(f3)
    PaintDB = {'reads': ReadsDB.d, 'gene2transcript': BigTree.db, 'annotation': AnnotationDB.d,'classification': ClassificationDB.d, 'fusion': FusionDB}
    utils.pickle_make(PaintDB, os.path.join(baseOutDir,sample+suffix1+'.db.pickle'))


class StructureDB:
    def __init__(self):
        self.d = {}

    def add_id(self,id,structure):
        if id not in self.d:
            self.d[id] = structure

    def saveReadsinfo(self,A,sample):
        id = A.id
        info = A.reads_info()
        if sample not in self.d:
            self.d[sample] = {}
        self.d[sample][id] = info

    def saveAnnoDB(self,A):
        id = A.id
        anno = A.db_used
        for gene in anno:
            if gene not in self.d:
                self.d[gene] = anno[gene]

    def saveClassification(self,A,sample):
        id = A.id
        classification = A.classification
        subtype = A.subtype
        if classification in self.d:
            if sample in self.d[classification]:
                self.d[classification][sample][id] = subtype
            else:
                self.d[classification][sample] = {id:subtype}
        else:
            self.d[classification] = {sample:{id:subtype}}


def tpmRead(fileIn):
    res = {}
    try:
        with open(fileIn,'r') as f:
            for i in f.readlines():
                if not i.strip() or i.startswith('#'):
                    continue
                line = i.strip().split('\t')
                id = line[0]
                tpm = line[2]
                flc = line[1]
                res[id] = [flc,tpm]
    except Exception as e:
        print('ERROR: TPM file is invalid !',file=sys.stderr)
        return 0
    return res


def option(argv):
    from argparse import ArgumentParser as AP
    usages = "python3 {} -f genome -g gtf -2 short_reads_bed -3 long_reads_bed -o output".format(argv[0])
    p = AP(usage=usages)
    p.add_argument("-g", dest="gtf", metavar="[gtf]", help="gtf annotation file",required=True)
    p.add_argument("-f", dest="fa", metavar="[fa]", help="reference genome fasta file", required=True)
    p.add_argument("--hisat2", dest="hisat2", metavar="[hisat2 bed file]", help="hisat2 bed file")
    p.add_argument("--minimap2", dest="minimap2", metavar="[minimap2 bed file]", help="minimap2 bed file")
    p.add_argument("--gmap", dest="gmap", metavar="[gmap bed file]", help="gmap bed file")
    p.add_argument("--tpm", dest="tpm", metavar="[TPM exp file]", help="TPM exp file")
    p.add_argument("--cnn", dest="cnn", action='store_true',help="Use machine learning model to help correct the mapping reslut")
    p.add_argument("--report", dest="report", action="store_true", help="report when running")
    p.add_argument("-o", dest="output", metavar="[output_dir]", help="output directory of bed&stat file out.",required=True)
    p.add_argument("-n", dest="sample", metavar="[sample name]", help="sample unique name",required=True)
    p.add_argument("-s", dest="suffix", metavar="[sample suffix name]", help="sample name suffix [.annot]", default='.annot')
    p.add_argument("-p", dest="process", metavar="[int]", help="process [1]", type=int,default=1)

    if len(argv) == 1:
        p.print_help()
        exit(1)
    return p.parse_args(argv[1:])


if __name__ == '__main__':
    args = option(sys.argv)
    main(args.hisat2, args.minimap2, args.gmap, args.output,args.sample, args.suffix, gtf2db.dict_make(args.gtf),args.fa, cnn=args.cnn, TPM=args.tpm, report=args.report)