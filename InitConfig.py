"""

"""

__version__ = '0.0.1'
__BASECONFIG__ = 'https://github.com/captorr'
__BASELENGTH__ = 31

import os,sys
import shutil
import TransAnnot


def configRead(options_dict,config,basecheck=0):
    with open(config,'r') as f:
        for i in f.readlines():
            if i.startswith('#') or not i.strip():
                continue
            line = i.strip('\n').split('=')
            if basecheck: options_dict[line[0].strip()] = line[1].strip()
            elif line[1].strip(): options_dict[line[0].strip()] = line[1].strip()
    #if basecheck and len(options_dict) != __BASELENGTH__:
    #    warning()
    return options_dict


def configMake(iArgs):
    if iArgs.help:
        help()
        exit(1)
    if iArgs.version:
        print("Version: {}".format(TransAnnot.__version__))
        exit(1)
    baseConfig = TransAnnot.baseConfig
    if not os.path.exists(baseConfig):
        warning()
    config = configRead({},baseConfig,basecheck=1)
    if iArgs.config: config = configRead(config,iArgs.config,basecheck=0)
    if iArgs.minimap2: config['USE_MINIMAP2'] = iArgs.minimap2
    if iArgs.hisat2:
        if iArgs.hisat2 in ['0','1']:
            config['USE_HISAT2'] = iArgs.hisat2
        else:
            config['HISAT2_INDEX'] = iArgs.hisat2
    if iArgs.gmap:
        if iArgs.gmap in ['0','1']:
            config['USE_GMAP'] = iArgs.gmap
        else:
            config['GMAP_INDEX'] = iArgs.gmap
    if iArgs.fa: config['FASTA'] = iArgs.fa
    if iArgs.output: config['OUTPUT_DIR'] = iArgs.output
    if iArgs.genome: config['GENOME_FA'] = iArgs.genome
    if iArgs.gtf: config['GTF_ANNOTATION'] = iArgs.gtf
    if iArgs.process: config['PROCESS'] = iArgs.process
    if iArgs.length: config['READ_LENGTH'] = iArgs.length
    if iArgs.overlap: config['READ_OVERLAP'] = iArgs.overlap
    if iArgs.minLength: config['MIN_READ_LENGTH'] = iArgs.minLength
    if iArgs.tpm: config['TPM_LIST'] = iArgs.tpm
    if iArgs.report: config['REPORT_IN_RUNNING'] = iArgs.report
    if iArgs.name: config['SAMPLE_UNIQUE_NAME'] = iArgs.name
    argsCheck(config)
    if not config['SAMPLE_UNIQUE_NAME']: config['SAMPLE_UNIQUE_NAME'] = os.path.basename(config['FASTA']).split('.')[0]
    if config['USE_MINIMAP2'] == '0' and config['USE_MINIMAP2'] == '0' and config['USE_MINIMAP2'] == '0':
        print('ERROR: Choose at least one aligner!',file=sys.stderr)
        exit(1)
    return config


def argsCheck(c):
    fileCheck('fasta input',c['FASTA'])
    fileCheck('genome fasta',c['GENOME_FA'])
    fileCheck('gtf annotation',c['GTF_ANNOTATION'])
    samtoolsCheck(c)
    if c['TPM_LIST']:
        fileCheck('TPM file',c['TPM_LIST'])
    if c['USE_HISAT2'] == '1':
        softwareCheck('hisat2',c['HISAT2'])
        fileCheck('hisat2 index',c['HISAT2_INDEX'],basedir=1)
    elif c['USE_HISAT2'] == '0':
        pass
    else:
        print('ERROR: invalid flag in "USE_HISAT2" ! (need 0 or 1)',file=sys.stderr)
        exit()
    if c['USE_MINIMAP2'] == '1':
        softwareCheck('minimap2',c['MINIMAP2'])
    elif c['USE_MINIMAP2'] == '0':
        pass
    else:
        print('ERROR: invalid flag in [USE_MINIMAP2] ! (need 0 or 1)',file=sys.stderr)
        exit()
    if c['USE_GMAP'] == '1':
        softwareCheck('gmap',c['GMAP'])
        fileCheck('gmap index',c['GMAP_INDEX'],basedir=1)
    elif c['USE_GMAP'] == '0':
        pass
    else:
        print('ERROR: invalid flag in "USE_GMAP" ! (need 0 or 1)',file=sys.stderr)
        exit()
    try:
        rl = int(c['READ_LENGTH'])
        ro = int(c['READ_OVERLAP'])
        mrl = int(c['MIN_READ_LENGTH'])
        if ro >= rl or mrl >= rl:
            print("ERROR: [READ_OVERLAP]\[MIN_READ_LENGTH] shuld be less than ['READ_LENGTH'] !",file=sys.stderr)
            exit()
    except:
        print(r'ERROR: [READ_LENGTH]\[READ_OVERLAP]\[MIN_READ_LENGTH] should be int',file=sys.stderr)
        exit()


def softwareCheck(name,cmd):
    if not shutil.which(cmd):
        print("ERROR: You choose {} as aligner, but {} is not an environment variable!".format(name,cmd),file=sys.stderr)
        exit()


def samtoolsCheck(c):
    if not shutil.which(c['SAMTOOLS']):
        print("ERROR: samtools is not an environment variable!")
        exit()
    if not os.path.exists(c["GENOME_FA"]+'.fai'):
        print("WARNING: can't find genome fasta index file, try to make it.")
        cmd = '{} faidx {}'.format(c['SAMTOOLS'],c["GENOME_FA"])
        os.system(cmd)
        print('genome fa index finished!\n')


def fileCheck(name,path,basedir=0):
    flag = os.path.exists(path) if not basedir else os.path.exists(os.path.dirname(path))
    if not flag:
        print("ERROR: {} doesn't exist!".format(name),file=sys.stderr)
        exit()


def warning():
    print("ERROR: Invalid base config file !",file=sys.stderr)
    print("================================")
    print("Please make sure the base config file is the latest complete version and fit to your software !")
    print("You can download latest base config file by command:\n\twget {}".format(__BASECONFIG__))
    print("Or download latest software by command:\n\twget {}".format(TransAnnot.__codex__))
    print("================================")
    exit()


def help():
    print('todo')
    exit(1)


def chelp():
    print("Thank you for using this script, you can download the latest base config file at:\n{}".format(__BASECONFIG__))
    exit(1)


if __name__ == '__main__':
    chelp()