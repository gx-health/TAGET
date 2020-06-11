# TransAnnot 使用手册

TransAnnot 是一款用于转录本结构注释的软件。

## 软件安装

	git clone https://github....

## 环境依赖

 * HISAT2/MINIMAP2/GMAP中至少一款软件
 * samtools
 * python3 （仅使用基础模块）

## FAST RUN

	python main.py -f [fasta] -g [genome fasta] -o [output directory] -a[annot gtf] -p [process] --use_minimap2 [1] --use_hisat2 [hisat2 index]

## 软件运行

1. 软件所需参数可由以下3种方式传入，序号高的参数会覆盖低序参数：
	1. 软件目录下的基础配置文件TransAnnot.Config
	2. 参数 `-c [config]` 导入的临时配置文件
	3. 直接导入参数

2. 首次运行建议先调整基础配置文件[TransAnnot.Config](), 主要设置以后运行时不需要频繁改动的参数：
	* 相关软件的环境变量
	* 默认使用的映射软件、映射软件的INDEX文件路径
	* 参考基因组(FASTA)、参考基因组注释文件（GTF）、默认PROCESS数

3. 如有设置好的基础配置文件，运行某一样本时仅需提供**输入fasta文件**和**输出文件夹路径**即可。可由`-c config`设置这两行，也可由`-f [fatsa] -o [output]`设置

4. 软件支持reads\transcript\gene水平的表达量统计，需要通过`--tpm`参数或`-c config`输入基于**reads ID**的表达量列表

5. 软件支持多样本的联合分析\表达量统计，可通过单独运行[TransAnnotMerge.py]()实现

## 运行结果

如未改动基础配置文件[TransAnnot.Config]()中的后缀参数，输出结果将有：

* [{sample_id}.annot.bed]() 注释过reads的bed格式结果

* [{sample_id}.annot.stat]()  所有reads的注释统计结果

* [{sample_id}.annot.db.pickle]()  分析数据，用于后续可视化分析

* [{sample_id}.annot.cluster.gene]()  注释过reads的gene水平聚类结果

* [{sample_id}.annot.cluster.transcript]()  注释过reads的transcript水平聚类结果

* [{sample_id}.annot.cluster.reads]()  注释过reads的信息

### {sample_id}.annot.stat 各列信息：

* `ID`： reads ID
* `Classification`： classificatioin of reads
* `Subtype`: subtype of reads
* `Gene`: gene annotation or region in genome[chr1:100000-100500]
* `Transcript`: transcript annotation
* `Chrom`: chromosome
* `Strand`: strand
* `Seq_length`: reads length
* `Seq_exon`: exon number of reads
* `Ref_length`: length of transcript annotation
* `Ref_exon_num`: exon number of transcript annotatioin
* `diff_to_gene_start`: 5` site difference of reads and annotation gene in reference genome
* `diff_to_gene_end`: 3` site difference of reads and annotation gene in reference genome
* `diff_to_transcript_start`: 5` site difference of reads and annotation transcript in reference genome
* `diff_to_transcript_end`: 3` site difference of reads and annotation transcript in reference genome
* `exon_miss_to_transcript_start`: number of exon missed in 5` site between reads and transcript annotation
* `exon_miss_to_transcript_end`: number of exon missed in 3` site between reads and transcript annotation


## TransAnnot的转录本分类方式

基于SQANTI的分类做了一些改动，目前有以下几类：

### 转录本水平分类

* `FSM`: full splice site match
* `ISM`: incomplete splice site match
* `NIC`: novel in catalog
* `NNC`: novel not in catalog
* `GENIC`: genic
* `INTERGENIC`: intergenic
* `FUSION`: fusion
* `UNKNOWN`： unknown

### 外显子水平分类

* `KE`: known exon
* `LEKE`: left end known exon
* `REKE`: right end known exon
* `NEKSLE`: novel exon with known splice site in left end exon and has the unique region overlap with at least two known exons
* `NEKSRE`: novel exon known splice site in right end exon and has the unique region overlap with at least two known exons
* `IE`: intron retention: two known splice sites from the same transcript's sequential exon
* `NEDT`: novel exon with two known splice sites from different transcript
* `NELS`: novel exon with novel left splice site
* `NERS`: novel exon with novel right splice site
* `LEE`: left exon_extension： the novel splice site in the left end of the exon which is longer than any exons overlap with it
* `REE`: right exon_extension: the novel splice site in the right end of the exon which is longer than any exons overlap with it
* `NEDS`: novel exon:double novel splice sites overlap with at least one known exon
* `NEIG`: novel exon inner-gene：novel exon inside the gene and without any overlap with known exon
* `NEOG`: novel exon inter-gene：novel exon outside the gene
* `NELE`: novel exon with novel splice site in the far left exon
* `NERE`: novel exon with novel splice site in the far right exon
* `MDNS`: monoexon with double novel splice sites
