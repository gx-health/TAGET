"""
Author: Zhang Chengsheng
python3 script
"""

__version__ = '0.0.1'
__email__ = 'zhangchengsheng1@qq.com'
__codex__ = 'https://github.com/captorr'


import os,sys

if sys.version_info[0] < 3:
    sys.stderr.write("ERROR: python version > 3.0.0 required ! Now is {}.{}.{}\n".format(sys.version_info[0], sys.version_info[1],sys.version_info[2]))
    exit()

import InitConfig
import utils
import BedAnnot
import gtf2db


baseConfig = os.path.join(os.path.dirname(sys.argv[0]),'TransAnnot.Config')


def optinalArgs():
    from argparse import ArgumentParser as AP
    usages = "python3 {} [optinal]".format(sys.argv[0])
    p = AP(usage=usages)
    p.add_argument("-f", dest="fa", metavar="[fasta]", help="input file in fasta format")
    p.add_argument("-g", dest="genome", metavar="[genome]", help="genome file in fasta format")
    p.add_argument("-o", dest="output", metavar="[output]", help="output directory path")
    p.add_argument("-n", dest="name", metavar="[sample name]", help="sample unique name")
    p.add_argument("-a", dest="gtf", metavar="[gtf]", help="reference annotation file in gtf format")
    p.add_argument("-p", dest="process", metavar="[int]", help="Number of process used. [1]",type=int)
    p.add_argument("-c", dest="config", metavar="[config]", help="set the other args by DIY config file")

    p.add_argument("--use_minimap2", dest="minimap2",metavar='[1/0]', help="Use minimap2 as aligner [0]")
    p.add_argument("--use_hisat2", dest="hisat2", metavar="[hisat2 index/1/0]", help="Use hisat2 as aligner and supply the prefix of index created by hisat2-build [0]")
    p.add_argument("--use_gmap", dest="gmap", metavar="[gmap index/1/0]", help="Use gmap as aligner and supply the prefix of index created by gmap_build [0]")
    p.add_argument("--tpm", dest="tpm", metavar="[tpm file]",help="Supply a tpm/flc file.")
    p.add_argument("--length", dest="length", metavar="[int]", help="when hisat2 was used, set the length of short reads [100]", type=int)
    p.add_argument("--overlap", dest="overlap", metavar="[int]", help="when hisat2 was used, set the overlap length of short reads [80]", type=int)
    p.add_argument("--minLength", dest="minLength", metavar="[int]",help="when hisat2 was used, set the min length of short reads [30]", type=int)
    p.add_argument("--cnn", dest="cnn",action='store_true',help="Use machine learning model to help correct the mapping reslut")
    p.add_argument("--version", dest="version",action='store_true',help="show version")
    p.add_argument("--report", dest="report", metavar="[0/1]", help="Print information when running[0]")
    p.add_argument("-d", dest="help",action='store_true', help="show more detail")

    if len(sys.argv) == 1:
        p.print_help()
        exit(1)
    return p.parse_args(sys.argv[1:])


def main():
    args = optinalArgs()
    try:
        config = InitConfig.configMake(args)
    except Exception as e:
        print(e)
        InitConfig.warning()

    if not os.path.exists(config['OUTPUT_DIR']):
        os.makedirs(config['OUTPUT_DIR'])
    if config['USE_HISAT2'] == '1':
        utils.runSplit(config)

    minimap2Bed = utils.runMinimap2(config) if config['USE_MINIMAP2'] == '1' else None
    gmapBed = utils.runGmap(config) if config['USE_GMAP'] == '1' else None
    hisat2Bed = utils.runHisat2(config) if config['USE_HISAT2'] == '1' else None
    print('load gtf annotation')
    DB = gtf2db.dict_make(config['GTF_ANNOTATION'])
    print('Annotate transcript')
    BedAnnot.main(hisat2Bed, minimap2Bed, gmapBed, config['OUTPUT_DIR'],config['SAMPLE_UNIQUE_NAME'], config['ANNO_SUFFIX_1'], DB, config['GENOME_FA'],TPM=config['TPM_LIST'], report=config['REPORT_IN_RUNNING'])
    cmd=config['PYTHON']+' '+config['TAGET_DIR']+ ' /ISO_anno_combine.py -2 '+os.path.join(config['OUTPUT_DIR'],config['SAMPLE_UNIQUE_NAME'])+'.hisat2.bed -3 '+os.path.join(config['OUTPUT_DIR'],config['SAMPLE_UNIQUE_NAME'])+'.minimap2.bed -d '+ config['GTF_DB']+ ' -o '+os.path.join(config['OUTPUT_DIR'],config['SAMPLE_UNIQUE_NAME'])+'.anno.tmp'
    #print(cmd)
    os.system(cmd)
    print('\nFinished !\n')


if __name__ == '__main__':
    main()
