"""
Author: Zhang Chengsheng, @2020.06.02
"""

import utils
import os,sys


def isoformExpParse(stat):
    res = {}
    if not stat:
        return res
    with open(stat,'r',encoding='utf8') as f:
        for i in f.readlines():
            if not i.strip() or i.startswith('#') or i.startswith('ID'):
                continue
            line = i.strip().split('\t')
            t_id,classification,gene,transcript,chrom,FLC,TPM = line[0],line[1],line[3],line[4],line[5],line[17],line[18]
            if '|' in chrom:
                chrom = chrom.split('|')[0]
            if t_id not in res:
                res[t_id] = [t_id,classification,gene,transcript,chrom,FLC,TPM]
    return res


def configRead(config):
    """sampleID stat bed db"""
    c = {}
    with open(config,'r') as f:
        for i in f.readlines():
            if not i.strip()or i.startswith('#'):
                continue
            line = i.strip().split('\t')
            sample,stat,bed,db = line[0],line[1],line[2],line[3]
            c[sample] = [stat,bed,db]
    return c


def bedParse(bed):
    res = {}
    with open(bed,'r',encoding='utf8') as f:
        for i in f.readlines():
            if not i.strip() or i.startswith('#'):
                continue
            line = i.strip().split('\t')
            chrom,start,end,t_id = line[0],line[1],line[2],line[3]
            if chrom in [0,'0'] or start in ['0'] or end in ['0']:
                continue
            if t_id not in res:
                res[t_id] = [[int(start),int(end)]]
            else:
                res[t_id].append([int(start),int(end)])
    return res


def calcRegion(exons):
    return [min([i[0] for i in exons]),max([i[1] for i in exons])]


def easyMerge(config,outputDir,mode=''):
    samples = configRead(config)
    geneExpOut = os.path.join(outputDir,'gene.exp')
    transcriptExpOut = os.path.join(outputDir,'transcript.exp')
    db_out = os.path.join(outputDir,'merge.db.pickle')

    bigTree = TreeCluster(TPM=mode)
    dbr = DBmerge()
    dba = DBmerge()
    dbc = DBmerge()
    for i in samples:
        exp,bed,db = samples[i]
        expRes,bedRes,db = isoformExpParse(exp),bedParse(bed),utils.load_pickle(db)
        dbr.readsMerge(db)
        dba.annotMerge(db)
        dbc.classMerge(db)
        if expRes:
            for j in expRes:
                t_id, classification, gene, transcript, chrom, FLC, TPM = expRes[j]
                if t_id not in bedRes:
                    continue
                exons = bedRes[t_id]
                region = calcRegion(exons)
                bigTree.add(i,t_id,chrom,gene,transcript,region,exons,classification,part_idx=db['reads'][i][t_id]['parts_idx'],tpm=TPM,flc=FLC)
        else:
            for t_id in bedRes:
                exons = bedRes[t_id]
                region = calcRegion(exons)
                chrom = db['reads'][i][t_id]['chrom'][0]
                gene, transcript,classification = db['reads'][i][t_id]['info'][:3]
                bigTree.add(i, t_id, chrom, gene, transcript, region, exons, classification,part_idx=db['reads'][i][t_id]['parts_idx'], tpm=0, flc=0)
    bigTree.forest()
    bigTree.fusionTreeToDB()
    t1,t2 = bigTree.expMatrix(list(samples),mode=mode)  ## mode = TPM/FLC
    with open(geneExpOut,'w') as o1, open(transcriptExpOut,'w') as o2:
        o1.write(t1)
        o2.write(t2)
    for i in samples:
        g1,g2,g3=bigTree.DBView(i)
        out1 = os.path.join(outputDir,'{}.gene.exp'.format(i))
        out2 = os.path.join(outputDir,'{}.transcript.exp'.format(i))
        out3 = os.path.join(outputDir,'{}.reads.exp'.format(i))
        with open(out1,'w') as o1,open(out2,'w') as o2, open(out3,'w') as o3:
            o1.write(g1)
            o2.write(g2)
            o3.write(g3)
    PaintDB = {'reads': dbr.db, 'gene2transcript': bigTree.db, 'annotation': dba.db,'classification': dbc.db, 'fusion': bigTree.fusionDB}
    utils.pickle_make(PaintDB, db_out)


class TreeCluster:
    def __init__(self,TPM=''):
        self.TPM= 1 if TPM in ['FLC','TPM','tpm','flc'] else 0
        self.treeNIC = {}
        self.treeNNC = {}
        self.treeGenic = {}
        self.treeIntergenic = {}
        self.treeFusion = {}
        self.treeGenic = {}
        self.db = {}
        self.fusionDB = {}

    def addSMnode(self,sample,t_id,chrom,gene,transcript,classification,tpm=0,flc=0):
        if chrom not in self.db:
            self.db[chrom] = {}
        if gene not in self.db[chrom]:
            self.db[chrom][gene] = {}
        if transcript not in self.db[chrom][gene]:
            self.db[chrom][gene][transcript] = {}
        if sample not in self.db[chrom][gene][transcript]:
            self.db[chrom][gene][transcript][sample] = {}
        if t_id not in self.db[chrom][gene][transcript][sample]:
            self.db[chrom][gene][transcript][sample][t_id] = {'TPM':tpm,'FLC':flc,'classification':classification}

    def _initNCnode(self,sample,t_id,region,exons,classification,tpm=0,flc=0):
        return {'zcs':{'idx':region,'cluster':{1:{'idx':region,'ref':exons,'reads':{sample:{t_id:{'TPM':tpm,'FLC':flc,'classification':classification}}}}}},'up':0,'down':0}

    def addNode(self,tree,sample,t_id,chrom,gene,region,exons,classification,tpm=0,flc=0,node=0):
        if ':' not in gene:
            if chrom not in tree:
                tree[chrom] = {}
            if gene not in tree[chrom]:
                tree[chrom][gene] = {}
                novelName = '{}_{}_{}'.format(gene,classification,1)
                tree[chrom][gene][novelName] = {'idx':region,'name':novelName,'ref':exons,'reads':{sample:{t_id:{'TPM':tpm,'FLC':flc,'classification':classification}}}}
            else:
                for cluster in tree[chrom][gene]:
                    if utils.check_overlap(tree[chrom][gene][cluster]['idx'], region):
                        same, new_ref, new_range = utils.ss_compare(tree[chrom][gene][cluster]['ref'], exons)
                        if same or classification in ['Genic','Intergenic']:
                            tree[chrom][gene][cluster]['idx'] = [min(tree[chrom][gene][cluster]['idx']+region),max(tree[chrom][gene][cluster]['idx']+region)]
                            tree[chrom][gene][cluster]['ref'] = new_ref
                            if sample in tree[chrom][gene][cluster]['reads']:
                                tree[chrom][gene][cluster]['reads'][sample][t_id] = {'TPM': tpm, 'FLC': flc,'classification': classification}
                            else:
                                tree[chrom][gene][cluster]['reads'][sample] = {t_id: {'TPM': tpm, 'FLC': flc, 'classification': classification}}
                            break
                else:
                    novelName = '{}_{}_{}'.format(gene, classification, len(tree[chrom][gene]) + 1)
                    tree[chrom][gene][novelName] = {'idx': region, 'name': novelName, 'ref': exons, 'reads': {sample: {t_id: {'TPM': tpm, 'FLC': flc, 'classification': classification}}}}
        else:
            if chrom not in tree:
                tree[chrom] = {}
            if 'unknown' not in tree[chrom]:
                tree[chrom]['unknown'] = self._initNCnode(sample,t_id,region,exons,classification,tpm=tpm,flc=flc)
            else:
                if not node: node = tree[chrom]['unknown']
                if utils.check_overlap(node['zcs']['idx'],region):  # zcs
                    node['zcs']['idx'] = [min(node['zcs']['idx']+region),max(node['zcs']['idx']+region)]
                    for nidx in node['zcs']['cluster']:
                        same, new_ref, new_range = utils.ss_compare(node['zcs']['cluster'][nidx]['ref'], exons)
                        if same or classification in ['Genic','Intergenic']:
                            node['zcs']['cluster'][nidx]['ref'] = new_ref
                            node['zcs']['cluster'][nidx]['idx'] = [min(node['zcs']['cluster'][nidx]['idx']+region),max(node['zcs']['cluster'][nidx]['idx']+region)]
                            if sample in node['zcs']['cluster'][nidx]['reads']:
                                node['zcs']['cluster'][nidx]['reads'][sample][t_id] = {'TPM':tpm,'FLC':flc,'classification':classification}
                            else:
                                node['zcs']['cluster'][nidx]['reads'][sample] = {t_id:{'TPM':tpm,'FLC':flc,'classification':classification}}
                            break
                    else:  # add new ns
                        nidx = len(node['zcs']['cluster']) + 1
                        node['zcs']['cluster'][nidx] = {'idx':region,'ref':exons,'reads':{sample:{t_id:{'TPM':tpm,'FLC':flc,'classification':classification}}}}
                elif region[0] < node['zcs']['idx'][0]:  # left
                    if node['up']:
                        self.addNode(tree,sample, t_id, chrom, gene, region, exons, classification, tpm=tpm, flc=flc,node=node['up'])
                    else:
                        node['up'] = self._initNCnode(sample,t_id,region,exons,classification,tpm=tpm,flc=flc)
                else:  # right
                    if node['down']:
                        self.addNode(tree,sample, t_id, chrom, gene, region, exons, classification, tpm=tpm, flc=flc,node=node['down'])
                    else:
                        node['down'] = self._initNCnode(sample,t_id,region,exons,classification,tpm=tpm,flc=flc)

    def _initFUSIONnode(self,sample,gene,t_id,classification,exons,part_idx,tpm=0,flc=0):
        genes = sorted(gene.split('|'))
        return {'genes':genes,'cluster':{1:{'genes':genes,'ref':exons,'reads':{sample:{t_id:{'TPM':tpm,'FLC':flc,'classification':classification,'pidx':part_idx}}}}}}

    def addFUSIONnode(self,sample,t_id,classification,gene,exons,part_idx,tpm=0,flc=0):
        if not self.treeFusion:
            self.treeFusion[0] = self._initFUSIONnode(sample,gene,t_id,classification,exons,part_idx,tpm=tpm,flc=flc)
        else:
            for node in self.treeFusion:
                genesAnnot = self.treeFusion[node]["genes"]
                genes = sorted(gene.split('|'))
                geneNotIn = []
                in_flag = 0
                for i in genes:
                    if i in genesAnnot:
                        in_flag = 1
                    else:
                        geneNotIn.append(i)
                if in_flag:
                    in_flag = 0
                    self.treeFusion[node]['genes'] = self.treeFusion[node]['genes'] + geneNotIn
                    for i in self.treeFusion[node]['cluster']:
                        genesAnnotSub = self.treeFusion[node]['cluster'][i]['genes']
                        exonsRef = self.treeFusion[node]['cluster'][i]['ref']
                        geneNotInSub = [i for i in genes if i not in genesAnnotSub]
                        if len(genesAnnotSub) == len(genes) and not geneNotInSub:
                            same,exonRefNew,exonRange = utils.ss_compare(exonsRef,exons)
                            if same:
                                self.treeFusion[node]['cluster'][i]['ref'] = exonRefNew
                                if sample in self.treeFusion[node]['cluster'][i]['reads']:
                                    self.treeFusion[node]['cluster'][i]['reads'][sample][t_id] = {'TPM':tpm,'FLC':flc,'classification':classification,'pidx':part_idx}
                                else:
                                    self.treeFusion[node]['cluster'][i]['reads'][sample] = {t_id:{'TPM':tpm,'FLC':flc,'classification':classification,'pidx':part_idx}}
                                break
                    else:
                        self.treeFusion[node]['cluster'][max(list(self.treeFusion[node]['cluster']))+1] = {'genes':genes,'ref':exons,'reads':{sample:{t_id:{'TPM':tpm,'FLC':flc,'classification':classification,'pidx':part_idx}}}}
                    break
            else:
                self.treeFusion[max(list(self.treeFusion))+1] = self._initFUSIONnode(sample,gene,t_id,classification,exons,part_idx,tpm=tpm,flc=flc)

    def add(self, sample, t_id, chrom, gene, transcript, region, exons, classification, part_idx=0, tpm=0, flc=0, node=0):
        if classification == 'NNC':
            self.addNode(self.treeNNC, sample, t_id, chrom, gene, region, exons, classification, tpm=tpm, flc=flc, node=node)
        elif classification == 'NIC':
            self.addNode(self.treeNIC, sample, t_id, chrom, gene, region, exons, classification, tpm=tpm, flc=flc, node=node)
        elif classification == 'Genic':
            self.addNode(self.treeGenic, sample, t_id, chrom, gene, region, exons, classification, tpm=tpm, flc=flc, node=node)
        elif classification == 'Intergenic':
            self.addNode(self.treeIntergenic, sample, t_id, chrom, gene, region, exons, classification, tpm=tpm, flc=flc, node=node)
        elif classification in ['FSM','ISM']:
            self.addSMnode(sample,t_id,chrom,gene,transcript,classification,tpm=tpm,flc=flc)
        elif classification in ['FUSION']:
            self.addFUSIONnode(sample, t_id, classification, gene, exons, part_idx=part_idx, tpm=tpm, flc=flc)

    def DBviewFusion(self,sample,showNone=True):
        res = {}
        isoformTXT = ''
        isoformTXTblank = '-\t{}\t{}\t{}\t{}\t{}\n' if self.TPM else '-\t{}\t{}\t{}\n'
        transcriptTXT = ''
        transcriptTXTblank = '-\t{}\t{}\t{}\t{}\t{}\t{}\n' if self.TPM else '-\t{}\t{}\t{}\t{}\n'
        geneTXT = ''
        geneTXTblank = '-\t{}\t{}\t{}\t{}\t{}\n' if self.TPM else '-\t{}\t{}\t{}\n'
        for node in self.treeFusion:
            gene = '|'.join(self.treeFusion[node]["genes"])
            if gene not in res: res[gene] = {}
            tsCount, gTPM, gFLC = 0, 0, 0
            NoneFlag = 1  # 1:empty, 0:at least 1
            ts_ids = []
            for cluster in self.treeFusion[node]["cluster"]:
                transcript = '|'.join(self.treeFusion[node]["cluster"][cluster]["genes"])
                if transcript not in res[gene]: res[gene][transcript] = {}
                if sample in self.treeFusion[node]["cluster"][cluster]['reads']:
                    if sample not in res[gene][transcript]: res[gene][transcript][sample] = {}
                    NoneFlag = 0
                    isoCount = len(self.treeFusion[node]["cluster"][cluster]['reads'][sample])
                    tsCount += isoCount
                    tTPM, tFLC = 0, 0
                    t_ids = []
                    for t_id in self.treeFusion[node]["cluster"][cluster]['reads'][sample]:
                        if t_id not in res[gene][transcript][sample]: res[gene][transcript][sample][t_id] = self.treeFusion[node]["cluster"][cluster]['reads'][sample][t_id]['pidx']
                        t_ids.append(t_id)
                        iTPM = float(self.treeFusion[node]["cluster"][cluster]['reads'][sample][t_id]['TPM'])
                        tTPM += iTPM
                        iFLC = int(self.treeFusion[node]["cluster"][cluster]['reads'][sample][t_id]['FLC'])
                        tFLC += iFLC
                        isoformTXT += isoformTXTblank.format(gene, transcript, t_id, iTPM, iFLC)
                    transcriptTXT += transcriptTXTblank.format(gene, transcript, isoCount, ','.join(t_ids), tTPM,tFLC)
                    ts_ids += t_ids
                    if self.TPM:
                        gTPM += tTPM
                        gFLC += tFLC
                elif showNone:
                    transcriptTXT += transcriptTXTblank.format(gene, transcript, 0, '', 0, 0)
            if not NoneFlag:
                geneTXT += geneTXTblank.format(gene, tsCount, ','.join(ts_ids), gTPM, gFLC)
            elif showNone:
                geneTXT += geneTXTblank.format(gene, 0, '', 0, 0)
        return geneTXT, transcriptTXT, isoformTXT, res

    def DBView(self,sample,showNone=True):
        isoformTXT = 'Chrom\tGene\tTranscript\tReadID\tTPM\tFLC\n' if self.TPM else 'Chrom\tGene\tTranscript\tReadID\n'
        isoformTXTblank = '{}\t{}\t{}\t{}\t{}\t{}\n' if self.TPM else '{}\t{}\t{}\t{}\n'
        transcriptTXT = 'Chrom\tGene\tTranscript\tReadsCount\tReadsID\tTPM\tFLC\n' if self.TPM else 'Chrom\tGene\tTranscript\tReadsCount\tReadsID\n'
        transcriptTXTblank = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n' if self.TPM else '{}\t{}\t{}\t{}\t{}\n'
        geneTXT = 'Chrom\tGene\tReadsCount\tReadsID\tTPM\tFLC\n' if self.TPM else 'Chrom\tGene\tReadsCount\tReadsID\n'
        geneTXTblank = '{}\t{}\t{}\t{}\t{}\t{}\n' if self.TPM else '{}\t{}\t{}\t{}\n'
        for chrom in self.db:
            for gene in self.db[chrom]:
                tsCount, gTPM, gFLC = 0, 0, 0
                NoneFlag = 1  # 1:empty, 0:at least 1
                ts_ids = []
                for transcript in self.db[chrom][gene]:
                    if sample in self.db[chrom][gene][transcript]:
                        NoneFlag = 0
                        isoCount = len(self.db[chrom][gene][transcript][sample])
                        tsCount += isoCount
                        tTPM,tFLC = 0,0
                        t_ids = []
                        for t_id in self.db[chrom][gene][transcript][sample]:
                            t_ids.append(t_id)
                            iTPM = float(self.db[chrom][gene][transcript][sample][t_id]['TPM'])
                            tTPM += iTPM
                            iFLC = int(self.db[chrom][gene][transcript][sample][t_id]['FLC'])
                            tFLC += iFLC
                            isoformTXT += isoformTXTblank.format(chrom, gene, transcript, t_id, iTPM, iFLC)
                        transcriptTXT += transcriptTXTblank.format(chrom, gene, transcript, isoCount,','.join(t_ids), tTPM,tFLC)
                        ts_ids += t_ids
                        if self.TPM:
                            gTPM += tTPM
                            gFLC += tFLC
                    elif showNone:
                        transcriptTXT += transcriptTXTblank.format(chrom, gene, transcript, 0,'', 0, 0)
                if not NoneFlag:
                    geneTXT += geneTXTblank.format(chrom,gene,tsCount,','.join(ts_ids),gTPM,gFLC)
                elif showNone:
                    geneTXT += geneTXTblank.format(chrom, gene, 0, '', 0, 0)
        return geneTXT, transcriptTXT, isoformTXT

    def fusionTreeToDB(self):
        for node in self.treeFusion:
            gene = '|'.join(self.treeFusion[node]["genes"])
            if gene not in self.fusionDB: self.fusionDB[gene] = {}
            for cluster in self.treeFusion[node]["cluster"]:
                transcript = '|'.join(self.treeFusion[node]["cluster"][cluster]["genes"])
                if transcript not in self.fusionDB[gene]: self.fusionDB[gene][transcript] = {}
                for sample in self.treeFusion[node]["cluster"][cluster]["reads"]:
                    if sample not in self.fusionDB[gene][transcript]: self.fusionDB[gene][transcript][sample] = {}
                    for t_id in self.treeFusion[node]["cluster"][cluster]["reads"][sample]:
                        if t_id not in self.fusionDB[gene][transcript][sample]:
                            self.fusionDB[gene][transcript][sample][t_id] = self.treeFusion[node]["cluster"][cluster]["reads"][sample][t_id]

    def treeToDB(self,tree):
        def nodes_Traversal(nodes):
            if nodes['up']:
                nodes_Traversal(nodes['up'])
            geneName = '{}:{}-{}'.format(chrom,min(nodes['zcs']['idx']),max(nodes['zcs']['idx']))
            if geneName not in self.db[chrom]:
                self.db[chrom][geneName] = {}
            for cluster in nodes['zcs']['cluster']:
                transcriptName = '{}:{}-{}'.format(chrom,min(nodes['zcs']['cluster'][cluster]['idx']),max(nodes['zcs']['cluster'][cluster]['idx']))
                self.db[chrom][geneName][transcriptName] = nodes['zcs']['cluster'][cluster]['reads']
            if nodes['down']:
                nodes_Traversal(nodes['down'])

        for chrom in tree:
            if chrom not in self.db: self.db[chrom] = {}
            for gene in tree[chrom]:
                if gene not in ['unknown']:
                    if gene not in self.db[chrom]: self.db[chrom][gene] = {}
                    for cluster in tree[chrom][gene]:
                        transcript = tree[chrom][gene][cluster]['name']
                        self.db[chrom][gene][transcript] = tree[chrom][gene][cluster]['reads']
                else:
                    nodes_Traversal(tree[chrom]['unknown'])

    def forest(self):
        self.treeToDB(self.treeNIC)
        self.treeToDB(self.treeNNC)
        self.treeToDB(self.treeGenic)
        self.treeToDB(self.treeIntergenic)

    def expMatrix(self,samples,mode='TPM'):
        geneTXT= 'Gene\t{}\n'.format('\t'.join(samples))
        transcriptTXT = 'Gene\tTranscript\t{}\n'.format('\t'.join(samples))
        for chrom in self.db:
            for gene in self.db[chrom]:
                tpmGene = []
                for transcript in self.db[chrom][gene]:
                    tpmTranscript = []
                    for sample in samples:
                        tpmSample = 0
                        if sample in self.db[chrom][gene][transcript]:
                            for t_id in self.db[chrom][gene][transcript][sample]:
                                tpmSample += float(self.db[chrom][gene][transcript][sample][t_id][mode])
                        tpmTranscript.append(tpmSample)
                    transcriptTXT += '{}\t{}\t{}\n'.format(gene,transcript,'\t'.join([str(i) for i in tpmTranscript]))
                    tpmGene = [tpmGene[i]+tpmTranscript[i] for i in range(len(tpmGene))] if tpmGene else tpmTranscript
                geneTXT += '{}\t{}\n'.format(gene,'\t'.join([str(i) for i in tpmGene]))
        return geneTXT,transcriptTXT


class DBmerge:
    def __init__(self):
        self.db = {}

    def readsMerge(self,db):
        for sample in db['reads']:
            if sample not in self.db:
                self.db[sample] = db['reads'][sample]

    def annotMerge(self,db):
        for gene in db['annotation']:
            if gene not in self.db:
                self.db[gene] = db['annotation'][gene]

    def classMerge(self,db):
        for c in db['classification']:
            if c not in self.db:
                self.db[c] = db['classification'][c]
            for sample in db['classification'][c]:
                if sample not in self.db[c]:
                    self.db[c][sample] = db['classification'][c][sample]


def option(argv):
    from argparse import ArgumentParser as AP
    usages = "python3 {} -c config -o outputdir".format(argv[0])
    p = AP(usage=usages)
    p.add_argument("-c", dest="config", metavar="[config]", help="config file", required=True)
    p.add_argument("-o", dest="output", metavar="[output_dir]", help="output directory",required=True)
    p.add_argument("-m", dest="mode", metavar="[TPM/FLC]", help="chose TPM or FLC or None, default: None", required=False,default='')

    if len(argv) == 1:
        p.print_help()
        exit(1)
    return p.parse_args(argv[1:])


if __name__ == '__main__':
    args = option(sys.argv)
    easyMerge(args.config,args.output,args.mode)