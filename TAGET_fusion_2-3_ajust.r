#! /usr/bin/Rscript
#====================================================================
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stringr"))
#====================================================================
option_list <- list(
  make_option(c("-l", "--long"), type="character", default= FALSE, action = "store",
              help="bed files from Minimap2"),
  make_option(c("-s", "--short"), type="character", default= FALSE, action = "store",
              help="bed files from Hisat2"),
  make_option(c("-a", "--stat"), type="character", default= FALSE, action = "store",
              help="STAT files from TAGET"),
  make_option(c("-f", "--fa"), type="character", default= FALSE, action = "store",
              help="fa files from TAGET"),
  make_option(c("-j", "--py"), type="character", default= FALSE, action = "store",
              help="python script for selecting"),
  make_option(c("-e", "--select"), type="character", default= FALSE, action = "store",
              help="python script for selecting stat"),
  make_option(c("-t", "--gtf"), type="character", default= FALSE, action = "store",
              help="gtf files of genome"),
  make_option(c("-n", "--name"), type="character", default= FALSE, action = "store",
              help="name of output files"),
  make_option(c("-o", "--output"), type="character", default= "Fusion_annotation_TAGET", 
              help="To specify output path"),
  make_option(c("-c", "--chu"), type="character", default= FALSE, action = "store",
              help="python script for chuli"))

parser <- OptionParser(usage = "%prog [options]", option_list=option_list,
                       description = "\nThe script is to annotate breakpoint of fusions from TAGET bed results (Hisat2 or Minimap2).",
                       epilogue="Yuhao Dong (yhdong18@fudan.edu.cn),  0104, 2022.\n"
)

arguments <- parse_args(parser)

longfile=arguments$long
shortfile=arguments$short
statfile=arguments$stat
fafile=arguments$fa

pyfile=arguments$py
pystatfile=arguments$select
pychuli=arguments$chu
gtffile=arguments$gtf
path=arguments$output
name=arguments$name

if(!dir.exists(path)){
system(paste0("mkdir ", path))
}
if(!dir.exists(paste0(path,"/data"))){
system(paste0("mkdir ",path, "/data"))
}


#
options(stringsAsFactors = FALSE)

#gene count function
gene_count <- function(bed){
    bed$gene_count <- 1
    for(v in 1:(nrow(bed)-1)){
        if((bed[v, 'Transcript'] == bed[v+1, 'Transcript'])&(bed[v, 'id'] != bed[v+1, 'id'])){
            bed[v+1, 'gene_count'] <- bed[v, 'gene_count'] + 1}
        
        if((bed[v, 'Transcript'] == bed[v+1, 'Transcript'])&(bed[v, 'id'] == bed[v+1, 'id'])){
            bed[v+1, 'gene_count'] <- bed[v, 'gene_count']}
    }
    return(bed)
}

#unlist function
agg_unlist <- function(y){
  y1 <- gsub("[c(]","",gsub("[\"]","", y))
  y2 <- gsub(" ","",gsub("[)]","", y1))
  y3 <- gsub(",","|",gsub("[\n]","", y2))
  return(y3)
}


#gtf file to bed file
if(!file.exists(paste0(path, "/genome_gene.bed"))){
gtf <- read.table(gtffile, sep="\t", header=F)
gene <- subset(gtf, V3 == "gene")
gene$name <- str_split_fixed(gene$V9, "[;]", n=5)[,3]
gene$gene_name <-gsub("gene_name ", "", gene$name)
write.table(gene[c(1,4,5,11)], paste0(path,"/genome_gene.bed"), sep="\t", row.names=F, col.names=F, quote=F)
}

if(!file.exists(paste0(path, "/genome_exon.bed"))){
gtf <- read.table(gtffile, sep="\t", header=F)
exon <- subset(gtf, V3 == "exon")
exon$exon_number <- gsub(" ","",str_split_fixed(exon$V9, "[;]", n=17)[,5])
exon$tran <- str_split_fixed(exon$V9, "[;]", n=17)[,9]
exon$name <- paste(exon$tran, exon$exon, sep="_")
exon$name <- gsub("transcript_name ","",exon$name)
exon$name <- gsub("_number","",exon$name)
write.table(exon[c(1,4,5,12)], paste0(path,"/genome_exon.bed"), sep="\t", row.names=F, col.names=F, quote=F)
}

#fusion bed file selecting

system(paste0("python ", pyfile," ",statfile," ",longfile," ",shortfile," ",fafile," > ", path,"/data/",name,".selecting.bed"))

system(paste0("python ", pystatfile, " ",statfile, " > ", path,"/data/",name,".select.stat"))

print("selecting complete")

#bed merge
s <- read.csv(paste0(path,"/data/",name,".selecting.bed"), sep="\t",header=F)

long <- read.csv(longfile, sep="\t",header=F)
short <- read.csv(shortfile, sep="\t",header=F)

sl <- subset(s, V2 == "Long")
ss <- subset(s, V2 == "Short")

b <- rbind(long[long$V4 %in% sl$V1,], short[short$V4 %in% ss$V1,])

#bedtools annotation
b <- subset(b, V1 > 0)

write.table(b, paste0(path,"/data/",name,".data.bed"), sep="\t",row.names=F, col.names=F,quote=F)
system(paste0("bedtools intersect -a ", path,"/data/",name,".data.bed", " -b ", path, "/genome_gene.bed -loj | bedtools groupby -i - -g 1-8 -c 12 -o collapse > ",path,"/data/",name,".annotation.bed"))
#print(paste0("python ",pychuli," ",path,"/data/",name,".annotation.bed",path,"/data/",name,".annotation1.bed"))
system(paste0("python ",pychuli," ",path,"/data/",name,".annotation.bed"," ",path,"/data/",name,".annotation1.bed"))
system(paste0("bedtools intersect -a ", path,"/data/",name,".annotation1.bed", " -b ", path, "/genome_exon.bed -loj | bedtools groupby -i - -g 1-9 -c 13 -o collapse > ",path,"/data/",name,".annotation_exon.bed"))

print("annotation complete")

#read annotation bed files
dat <- read.csv(paste0(path,"/data/",name,".annotation_exon.bed"), sep="\t", header=F)[-c(5)]
colnames(dat) = c('Chr','Start','End','Transcript','Strand','len_start','len_end','all_gene','all_exon')
if(nrow(dat[dat$all_exon == ".",]) > 0){
    dat[dat$all_exon == ".",]$all_exon <- "non-exon"
}

dat$gene <- str_split_fixed(dat$all_gene, "[,]", n=2)[,1]
dat$length <- dat$End-dat$Start

###cut length
dat <- subset(dat, length >= 20)

dat$coor <- paste(dat$Chr, dat$Start, dat$End,sep="_")
dat[dat$gene == ".",]$gene <- dat[dat$gene == ".",]$coor
dat$id <- paste(dat$Transcript, dat$gene, sep="_")


#filter fusion transcripts
uni <- dat[!duplicated(dat$id),]
tran <- as.data.frame(table(uni$Transcript))

##gene
#ge <- subset(tran, Freq == 1)
#datt_gene <- dat[dat$Transcript %in% ge$Var1,]
#datt_u_gene <- subset(datt_gene[!duplicated(datt_gene$id),], gene != ".")
#write.table(datt_u_gene, paste0(path,"/data/",name,".gene_annotation.txt"), sep="\t", row.names = F, col.names = T, quote = F)

##fusion
fu <- subset(tran, Freq > 1)
datt <- dat[dat$Transcript %in% fu$Var1,]

#filter overlap genes
t <- datt[grep("[,]",datt$all_gene),]
dt <- as.data.frame(table(t$Transcript))
dt <- subset(dt, Freq > 0)
ov <- as.data.frame(0)
for(i in 1:nrow(dt)){
    v <- dt[i,'Var1']
    v <- as.character(v)
    vt <- subset(datt, Transcript == v)
    x <- str_split(vt$all_gene, ",")

    ge <- as.data.frame(0)
    colnames(ge) <- "Var1"
    for(j in 1:length(x)){
        t <- as.data.frame(table(x[[j]]))
        ge <- rbind(ge, t[1])
    }

    if(max(table(ge[-1,])) == nrow(vt)){
        ov <- rbind(ov, v)
    }
}
ov <- ov[-1,]
datt <- datt[!datt$Transcript %in% ov,]

write.table(b[b$V4 %in% datt$Transcript,], paste0(path,"/data/",name,".fusion_filter.bed"), sep="\t",row.names=F, col.names=F,quote=F)
print("filter complete")

#length merge
gene_length <- aggregate(datt[,c('length')], list(datt[,'id']), sum)
colnames(gene_length) <- c("id", "length_merge")



#filter middle dupgenes
datt$number_gene <- 1

for(i in row.names(datt[grep(",",datt$all_gene),])){
    st <- as.data.frame(table(str_split(datt[i,'all_gene'],',')))
    datt[i,'number_gene'] <- nrow(st)
}

datt$id <- paste(datt$Transcript, datt$gene, sep="_")

datt$id_all <- paste(datt$Transcript, datt$all_gene, sep="_")
datt$idd <- paste(datt$Start, datt$id, sep="_")
datt$rank <- row.names(datt)

#remove middle
first <- row.names(unique(datt['id_all'], fromLast =T))
last <- row.names(unique(datt['id_all'], fromLast =F))

datt <- datt[c(first,last),]

data <- datt[!duplicated(datt$idd),]
data$rank <- as.numeric(data$rank)
data <- data[order(data$rank, decreasing = F),]

#combine-find min genes
row.names(data) <- 1:nrow(data)
ta <- data
for(i in 1:(nrow(data))){
    id <- data[i, 'idd']
    row.names(ta) <- 1:nrow(ta)
    rank <- as.numeric(row.names(ta[ta$idd %in% id,]))

    if(data[i,'number_gene']==1){
        ta[rank,'all_gene'] <- ta[rank,'gene']
        data[i,'all_gene'] <- data[i,'gene']
        next
        }

    tg <- str_split(data[i,'all_gene'], ',')

    same_before <- if(rank == 1){same_before = F}else{(ta[rank-1, 'Strand'] == data[i, 'Strand']) & (ta[rank-1, 'Transcript'] == data[i, 'Transcript'])}
    same_after <- if(rank == nrow(ta)){same_after = F}else{(ta[rank+1, 'Strand'] == data[i, 'Strand']) & (ta[rank+1, 'Transcript'] == data[i, 'Transcript'])}

    belong_before <- if(rank == 1){belong_before = F}else{all(str_split(ta[rank-1, 'all_gene'], ',')[[1]] %in% tg[[1]])}
    belong_after <- if(rank == nrow(ta)){belong_after = F}else{all(str_split(ta[rank+1, 'all_gene'], ',')[[1]] %in% tg[[1]])}

    
    if(belong_before & same_before){
        ta[rank-1, 'Start'] <- min(ta[rank-1, 'Start'], ta[rank, 'Start'])
        ta[rank-1, 'End'] <- max(ta[rank-1, 'End'], ta[rank, 'End'])
        ta <- ta[-rank,]

    }else if(belong_after & same_after){
        ta[rank+1, 'Start'] <- min(ta[rank+1, 'Start'], ta[rank, 'Start'])
        ta[rank+1, 'End'] <- max(ta[rank+1, 'End'], ta[rank, 'End'])
        ta <- ta[-rank,]

    }


}
ta <- na.omit(ta)


#remove middle exons
first <- row.names(unique(ta['id'], fromLast =T))
last <- row.names(unique(ta['id'], fromLast =F))

ta <- ta[c(first,last),]

data <- ta[!duplicated(ta$idd),]
data$rank <- as.numeric(data$rank)
data <- data[order(data$rank, decreasing = F),]

#deduplicate gene merge
dup <- as.data.frame(0)
for(i in 1:(nrow(data)-1)){
    if((data[i, 'id'] == data[i+1, 'id']) & (data[i, 'Strand'] == data[i+1, 'Strand'])){
        data[i+1, 'Start'] <- min(data[i+1, 'Start'], data[i, 'Start'])
        data[i+1, 'End'] <- max(data[i+1, 'End'], data[i, 'End'])

        e1 <- unlist(str_split(data[i, 'all_exon'], ","))
        e2 <- unlist(str_split(data[i+1, 'all_exon'], ","))

        e <- sort(unique(c(e1, e2)))
        ex <- 0
        for(x in 1:length(e)){ex <- paste(ex, e[x],sep=",")}

        data[i+1, 'all_exon'] <- gsub("0,", "", ex)
        dup <- rbind(dup, i)
    }
}
if(nrow(dup) == 1){
    bed <- data
}else{bed <- as.data.frame(data[-dup[-1,],])}
row.names(bed) <- 1:nrow(bed)

print("deduplicate complete")

#break
bed <- gene_count(bed)#gene_count function
bed$gene_count <- as.numeric(bed$gene_count)


bed$breakpoint <- 0
bed$breakpoint_length <- 0
bed$y <- 0

for(y in 1:nrow(bed)){
    ybed <- bed[bed$Transcript %in% bed[y, 'Transcript'],]
    ymax <- max(ybed$gene_count)
    if(bed[y, 'gene_count'] == 1){bed[y, 'y'] <- "first"}
    if(bed[y, 'gene_count'] < ymax & bed[y, 'gene_count'] > 1){bed[y, 'y'] <- "mid"}
    if(bed[y, 'gene_count'] == ymax){bed[y, 'y'] <- "last"}
    bed[y, 'gene_count'] <- ymax
}

bed$b <- paste(bed$Strand, bed$y, sep="_")

bed$start <- paste(bed$Chr, bed$Start, sep=":")
bed$end <- paste(bed$Chr, bed$End, sep=":")


for(x in 1:nrow(bed)){
    bed[x, 'breakpoint'] <- switch(bed[x, 'b'], 
        '-1_first' = bed[x, 'start'], 
        '-1_last' = bed[x, 'end'], 
        '1_first' = bed[x, 'end'], 
        '1_last' = bed[x, 'start'],
        '1_mid' = paste(bed[x, 'start'], bed[x, 'end'], sep="|"),
        '-1_mid' = paste(bed[x, 'start'], bed[x, 'end'], sep="|")
        )
    
    bed[x, 'breakpoint_length'] <- bed[x,'len_end']
}


#merge length
length <- merge(bed, gene_length, by="id",all=F)
bed <- length


#transcript vesion speculate
bed$Trans_version <- 0
for(i in 1:nrow(bed)){
    t <- str_split_fixed(unlist(str_split(bed[i, 'all_exon'], ",")), "[_]",n=2)[,1]
    tt <- as.data.frame((table(t)))
    tran <- tt[tt$Freq == max(tt$Freq),][1]
    tx <- 0
    for(x in 1:nrow(tran)){tx <- paste(tx, tran[x,],sep=";")}
    bed[i, 'Trans_version'] <- gsub("0;", "", tx)
}

xx <- bed[grep("[,]",bed$all_gene),]
x <- xx[grep("[.]",xx$all_gene),]

for(j in 1:nrow(x)){
    c <- str_split_fixed(x[j,]$all_gene, "[,]", n=4)
    x[j,]$gene <- c[-grep("[.]", c)][1]
}

bed[bed$idd %in% x$idd,]$gene <- x$gene

write.table(bed, paste0(path,"/data/",name,".fusion_annotation_raw.txt"), sep="\t", row.names = F, col.names = T, quote = F)
print("breakpoint complete")


#final combine
#t <- subset(raw, Subtype == "NEIGHBOUR_FUSION" & length_merge < 100)

test <- aggregate(bed[,c('all_exon')], list(bed[,'Transcript']), paste)
for(i in 1:nrow(test)){
    if(all(test$x[[i]] %in% "non-exon")){
        tr <- test[i,'Group.1']
        bed <- bed[!bed$Transcript %in% tr,]
    }

}

fusion <- aggregate(bed[,c('gene', 'Strand','breakpoint','breakpoint_length','length_merge','Trans_version')], list(bed[,'Transcript']), paste)
colnames(fusion)[1:3] <- c('Transcript', "FUSION","Strand")

row <- row.names(unique(bed['Transcript'], fromLast =T))
trans <- bed[row, c('Transcript', 'gene_count')]

m <- merge(trans, fusion, by="Transcript")

pa <- function(F){
d <- as.data.frame(0)
for(i in 1:nrow(F)){
    d <- rbind(d, paste(F[i,1], F[i,2], sep="|"))
}
return(d[-1,])
}

if(all(m$gene_count == 2)){
    F <- as.data.frame(m$FUSION)
    S <- as.data.frame(m$Strand)
    B <- as.data.frame(m$breakpoint)
    BL <- as.data.frame(m$breakpoint_length)
    L <- as.data.frame(m$length_merge)
    T <- as.data.frame(m$Trans_version)
    m$FUSION <- pa(F)
    m$Strand <- pa(S)
    m$breakpoint <- pa(B)
    m$breakpoint_length <- pa(BL)
    m$length_merge <- pa(L)
    m$Trans_version <- pa(T)
}else{
    m$FUSION <- agg_unlist(m$FUSION)
    m$Strand <- agg_unlist(m$Strand)
    m$breakpoint <- agg_unlist(m$breakpoint)
    m$breakpoint_length <- agg_unlist(m$breakpoint_length)
    m$length_merge <- agg_unlist(m$length_merge)
    m$Trans_version <- agg_unlist(m$Trans_version)
}


final <- m
final$FUSION <- gsub(" ", "", final$FUSION)
final$Trans_version <- gsub(" ", "", final$Trans_version)

write.table(final, paste0(path, "/",name,".fusion_annotation_merge.txt"), sep="\t", row.names=FALSE, quote=FALSE)
#write.csv(final, paste0(path, "/",name,".fusion_annotation_merge.txt"), sep="\t", row.names=F)

print(paste0(name, " all complete, ", nrow(final), " fusions reported"))








#new old fusion merge
#merge stat
options(stringsAsFactors = FALSE)

  c <- final
  colnames(c) <- paste0("new_", colnames(c))
  colnames(c)[1] <- "V1"
  
  s_stat <- read.table(paste0(path,"/data/",name,".select.stat"),sep="\t", fill=TRUE, header=FALSE,fileEncoding="utf8")[c(1,18:20)]
  #s_stat <- s_stat[grep("transcript", s_stat$V1),]
  
  stat <- read.table(statfile, sep="\t",fill=TRUE, header=FALSE, fileEncoding="utf8")[-1,1:4]
  #stat <- stat[grep("transcript", stat$V1),]
  
  s <- merge(stat, s_stat, by="V1",all=TRUE)
  
  data <- merge(c, s, by="V1", all=TRUE)
  data <- subset(data, V2 == "FUSION"| V19 == "FUSION"| new_gene_count > 0)
  data[is.na(data)] <- 0
  
  colnames(data)[1] <- 'Transcript'
  colnames(data)[8:13] <- c('old_Type','old_SubType','old_FUSION','alternative_FUSION','alternative_Type','alternative_SubType')
  
  #write.csv(data, paste0("new_fusions/merge/",i,"_clean.merge.csv"),row.names = F)
  
  write.table(data, paste0(path,"/",name,".new_old_FUSION_merge.txt"),quote=FALSE,sep="\t",row.names=FALSE)
  
  
#merge
plot <- merge(dat, final, by="Transcript", all=FALSE)
write.table(plot, paste0(path,"/",name,".fusion_exons_merge.txt"),quote=FALSE,sep="\t",row.names=FALSE)

