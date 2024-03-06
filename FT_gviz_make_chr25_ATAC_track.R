#aims
#itrack: from UCSC directly
#gtrack: depend on others
#atrack: cpgIslands, UCSC cpg.bed to cpg.gr (chr25 only)
#grtrack: genemodel.bed from customized GTF chr25 only

library(Gviz)

###make genomic backbone###
#load GTF or bed and transfer to gr 
setwd("/home/miller/Desktop/FT_for_Gviz")
getwd()

chr25_gtf <- "galGal6_with_newKRT_RefSeq_20210129_onlyRegChr25.sorted.gtf"
chr25_gr <- rtracklayer::import(chr25_gtf)
head(chr25_gr)
unique(chr25_gr@seqnames)

chr25_CpG <- "galGal6_UCSC_chr25_CpG.bed"
CpG_gr <- rtracklayer::import(chr25_CpG)
seqlevelsStyle(CpG_gr) <- "Ensembl"
unique(CpG_gr@seqnames)

rm(chr25_gtf,chr25_CpG)


#Annotation track, title ="CpG"
atrack <- AnnotationTrack(CpG_gr, name = "CpG",collapse = TRUE)
plotTracks(atrack)

## genomic coordinates
gtrack <- GenomeAxisTrack()
plotTracks(list(gtrack, atrack))

#Chromosme name : "25"
chr <- as.character(unique(seqnames(CpG_gr)))
chr
#Ideogram track (fetches data from UCSC)
itrack <- IdeogramTrack(genome = "galGal6", chromosome = '25', name='Chr25')
class(itrack)
#remove UCSC "chr"
levels(itrack@bandTable$chrom) <- sub("^chr", "", levels(itrack@bandTable$chrom), ignore.case=T)
View(itrack)

plotTracks(list(itrack, gtrack, atrack),
           from = 100000, to = 110000)


#Add gene model
#Load gene model data (GG6a_ens104_chr25_wt_customizedKRT_exon.bed)
geneModels <- read.table(file.choose(),sep="\t",header=T) 
head(geneModels)

grtrack <- GeneRegionTrack(geneModels, genome = "galGal6",
                           chromosome = chr, name = "Genes",
                           showId = TRUE)


grtrack <- GeneRegionTrack(geneModels, genome = gen,
                           chromosome = chr, name = "Gene Model", 
                           transcriptAnnotation = "symbol",
                           background.panel = "#FFFEDB",
                           background.title = "darkblue")

plotTracks(list(itrack, gtrack, atrack, grtrack),from = 2140000
, to = 2200000)


###Add data track###
E14F <- 'AT_E14F_stepB_dedup_opt.bigwig'
wingE <- 'AT_wingE_stepB_dedup_opt.bigwig'



E14F_track <- DataTrack(E14F,
                        genome = "galGal6",
                        chromosome = "chr25",
                        window = -1,
                        col.mountain='blue',	
                        fill.mountain=c("black", "blue"),
                        col.axis="black",
                        name = "E14F")

wingE_track <- DataTrack(wingE,
                        genome = "galGal6",
                        chromosome = "chr25",
                        window = -1,
                        col.mountain='blue',	
                        fill.mountain=c("black", "blue"),
                        col.axis="black",
                        name = "wingE")


#plot all tracks
plotTracks(list(itrack, E14F_track, wingE_track,atrack,gtrack, grtrack), 
           sizes=c(1,2,2,0.75,1,12), #track hight
           type="polygon", from = 2035000, to =2450000, 
           ylim = c(0,0.4), #ATAC track scale
           #lwd.border.title=0.2,
           cex = 4,cex.axis=2.5,cex.title = NULL, title.width = 0.35,fontcolor= "black",fontsize=2)




