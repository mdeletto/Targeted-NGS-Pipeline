#!/usr/bin/Rscript
library(dplyr)
library("ggplot2")
library("stringr")
library("DNAcopy")
library(rtracklayer)
library(GenomicRanges)
library(Homo.sapiens)
library(biovizBase)
library(biomaRt)
library(jsonlite)
library(RMySQL)
library(org.Hs.eg.db)
#library(ggrepel)
#library(BSgenome.Hsapiens.UCSC.hg19)
data("ideoCyto")

base_output <- "MP18-142"
diffCoverage <- "/home/michael/Downloads/diffCoverage.seg"
cnvseg <- "/home/michael/Downloads/test.cnv.seg.detailed.tsv"

targets <- "/home/michael/YNHH/Reference_Files/OCP/AmpliSeq_OCP/OCP.20150630.designed.bed"
targets <- "/home/michael/YNHH/Reference_Files/OCP/v3/OCAv3.20170110.designed.bed"
#targets <- "/home/michael/YNHH/Reference_Files/CHPv2/CHP2.Somatic.20131001.designed.bed"
#targets <- "/home/michael/YNHH/Reference_Files/CCP/CCP.20131001.designed.bed"
targets <- "/home/michael/YNHH/Reference_Files/AmpliSeqWhEx/AmpliSeqExome.20141113.designed.bed"


args = commandArgs(trailingOnly=TRUE)
base_output <- args[1]
diffCoverage <- args[2]
cnvseg <- args[3]
targets <- args[4]
  
  
  
  
########################################
### CONNECT TO BIOMART FOR GENE OVERLAP COMPARISONS
### NOTE: gene_symbols frochrNames <- as.vector(unique(CNA.segmented.df$chrName))m biovizBase is also used (legacy)
########################################

grch37 = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)

########################################
### READ GENE_IDs FROM AMPLISEQ BED FILE
########################################
target_regions <- read.table(targets, skip=1)
colnames(target_regions) <- c("chr","begin","finish")
target_regions.genomicrange <- transformDfToGr(target_regions, seqnames = "chr", start = "begin", end="finish")
chrNames <- as.vector(unique(target_regions$chr))
target_regions.genomicrange <- keepSeqlevels(target_regions.genomicrange, chrNames)
seqlengths(target_regions.genomicrange) <- as.numeric(seqlengths(ideoCyto$hg19)[chrNames])

########################################
### PULL IN DATABASE INFO FOR GENE ALIASES
########################################

# use sql to get alias table and gene_info table (contains the symbols)
# first open the database connection
dbCon <- org.Hs.eg_dbconn()
# write your SQL query
sqlQuery <- 'SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'
# execute the query on the database
aliasSymbol <- dbGetQuery(dbCon, sqlQuery)

checkGeneAlias <- function(gene){
  if (!(gene %in% aliasSymbol$symbol)){
    index <- which(aliasSymbol[,2] == gene)
    symbol <- aliasSymbol[index,5]
    return(symbol)
  }
}

########################################
### PATTERN MATCHING FOR GENE_ID CHANGES BASED ON TARGET BED
########################################

if (grepl("CCPHSMV2", targets)){
  index <- 5
  target_gene_list <- unique(target_regions[,index])
  target_gene_list.cnv_only <- target_gene_list
  confidence_cutoff <- 20
} else if (grepl("CHP2", targets)){
  index <- which(apply(target_regions, 2, function(x) any(grepl("GENE_ID", x))))
  target_gene_list <- unique(str_match(target_regions[,index], "GENE_ID=([:graph:]+)")[,2])
  target_gene_list.cnv_only <- target_gene_list
  confidence_cutoff <- 20
} else if (grepl("OCP", targets)){
  index <- which(apply(target_regions, 2, function(x) any(grepl("GENE_ID", x))))
  target_regions.cnv_only <- subset(target_regions, grepl("CNV_HS=1", target_regions[,index]))
  target_gene_list.cnv_only <- unique(str_match(target_regions.cnv_only[,index], "GENE_ID=([:graph:]+?)(\n|;)")[,2])
  target_gene_list <- unique(str_match(target_regions[,index], "GENE_ID=([:graph:]+?)(\n|;)")[,2])
  #target_gene_list <- c(target_gene_list, "MYCL1")
  confidence_cutoff <- 20
} else {
  index <- which(apply(target_regions, 2, function(x) any(grepl("GENE_ID", x))))
  target_gene_list <- unique(str_match(target_regions[,index], "GENE_ID=([:graph:]+?)(\n|;)")[,2])
  target_gene_list.cnv_only <- target_gene_list
  confidence_cutoff <- 20
  #target_gene_list <- c(target_gene_list, "MYCL1")
}

# Subset those target genes that require replacement
target_gene_list.aliases <- subset(target_gene_list, !(target_gene_list %in% aliasSymbol$symbol))
target_gene_list <- unlist(replace(target_gene_list, 
                                    which(target_gene_list %in% target_gene_list.aliases), 
                                    values=lapply(target_gene_list.aliases, FUN=checkGeneAlias)))


target_gene_list <- unlist(target_gene_list)

########################################
### LOAD DATA AND RESTRUCTURE
########################################

# Load data
data(ideoCyto, package = "biovizBase")
data(genesymbol, package = "biovizBase")

# Read segment files
diffCoverage <- read.table(file = diffCoverage, sep = '\t', header = TRUE)
cnvseg <- read.table(file = cnvseg, sep = '\t', header=TRUE)
colnames(cnvseg) <- c("chrName","begin","end","len","ploidy","numtiles","confidence","precision")

# Rename 'end' columns to 'finish' to not cause problems with GenomicRanges objects
names(diffCoverage)[names(diffCoverage)=="end"] <- "finish"
names(cnvseg)[names(cnvseg)=="end"] <- "finish"
# Reorder files with specific chromosomal order
diffCoverage$chrName <- as.character(diffCoverage$chrName)
cnvseg$chrName <- as.character(cnvseg$chrName)
#Then turn it back into an ordered factor

chrNames <- c(paste("chr",1:22,sep=""),"chrX","chrY","chrM")

diffCoverage$chrName <- factor(diffCoverage$chrName, levels=chrNames
)

cnvseg$chrName <- factor(cnvseg$chrName, levels=chrNames
)

diffCoverage <- diffCoverage[order(diffCoverage$chrName),]
cnvseg <- cnvseg[order(cnvseg$chrName),]


diffCoverage$ploidycutoff <- cut(diffCoverage$difference,
                                 breaks = c(-Inf, 4, Inf),
                                 labels = c("<=4",">4")
)


################################################
### APPLY NEW CBS TO BINS CALLED BY IONREPORTER
################################################

# Create CNA object
CNA.object <- CNA(genomdat=log2(diffCoverage$difference),
                  chrom=diffCoverage$chrName,
                  maploc=diffCoverage$begin,
                  data.type='logratio')

# Smooth bins
CNA.smoothed <- smooth.CNA(CNA.object)

# Perform CBS
CNA.segmented <- segment(CNA.smoothed, verbose=0, min.width=2)

# Add p-value and confidence limits to segments
CNA.segmented.df.pvalue <- segments.p(CNA.segmented, ngrid=100, tol=1e-6, alpha=0.05, search.range=100, nperm=1000)

# Convert seg.mean log2ratio back to ploidy
CNA.segmented.df.pvalue$ploidy <- 2 ^ CNA.segmented.df.pvalue$seg.mean
CNA.segmented.df.pvalue$ID <- base_output

CNA.segmented.df.pvalue.nochrom <- CNA.segmented.df.pvalue
CNA.segmented.df.pvalue.nochrom$chrom <- gsub(pattern = "chr", replacement = "", CNA.segmented.df.pvalue$chrom)

# Write CBS segments in .tsv format (.seg)

write.table(CNA.segmented.df.pvalue.nochrom[c("ID","chrom","loc.start","loc.end","num.mark","seg.mean")],
            file = paste0(base_output, ".cbs.seg", collapse = ""),
            sep = "\t",
            quote = F,
            row.names = F)

####################################################
### FIND OVERLAPPING GENES FOR SEGMENTS FROM ENSEMBL
####################################################

# change "chrom" column to "chrName"
# We'll need this for plotting later
colnames(CNA.segmented.df.pvalue)[2] <- "chrName"
CNA.segmented.df.pvalue$chrName <- as.character(CNA.segmented.df.pvalue$chrName)

# Transform CNA dataframe into GenomicRanges object
CNA.segmented.df.pvalue.genomicrange <- transformDfToGr(CNA.segmented.df.pvalue, seqnames = "chrName", start = "loc.start", end="loc.end")
chrNames <- as.vector(unique(CNA.segmented.df.pvalue$chrName))
CNA.segmented.df.pvalue.genomicrange <- keepSeqlevels(CNA.segmented.df.pvalue.genomicrange, chrNames)
seqlengths(CNA.segmented.df.pvalue.genomicrange) <- as.numeric(seqlengths(ideoCyto$hg19)[chrNames])


getGeneBoundaries <- function(gene){
  query <- getBM(attributes = c('chromosome_name', 'start_position', 'end_position'),
                 filters = c('hgnc_symbol'), 
                 values = list(gene),
                 mart = grch37)
  
  return(query)
  
}

getGenes <- function(chr, start, end){
  query <- getBM(filters = c('chromosome_name', 'start', 'end'),
                 attributes = c('hgnc_symbol'), 
                 values = list(chr, start, end),
                 mart = grch37)
  
  return(query)
  
}

apply_BM_overlap.all <- function(chr, start, end){
  query <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'hgnc_symbol','band'),
                 filters = c('chromosome_name', 'start', 'end'), 
                 values = c(as.list(sub("^chr", "", chr)), start, end),
                 mart = grch37)
  # Select unique genes
  gene_symbols.all <- query
  gene_symbols.all.unique <- unique(gene_symbols.all$hgnc_symbol)
  # Remove blank symbols
  remove_symbols <- c("")
  gene_symbols.all.unique <- gene_symbols.all.unique[! gene_symbols.all.unique %in% remove_symbols]
  # Subset genes only included in targeted assay
  gene_symbols.targeted <- subset(gene_symbols.all.unique, gene_symbols.all.unique %in% target_gene_list)
  gene_symbols.targeted.unique <- unique(gene_symbols.targeted) 
  # Get the cytoband information
  cytoband.min <- head(query$band, n=1)
  cytoband.max <- tail(query$band, n=1)
  
  return_values <- list("chrName" = chr,
                        "loc.start" = start,
                        "loc.end" = end,
                        "all.genes" = paste(gene_symbols.all.unique, collapse="|"),
                        "targeted.genes" = paste(gene_symbols.targeted.unique, collapse="|"),
                        "cytoband.start" = cytoband.min,
                        "cytoband.end" = cytoband.max
  )
  
  return(return_values)
}

# Apply apply_BM_overlap

gene_symbol_and_cytobands <- mapply(apply_BM_overlap.all, CNA.segmented.df.pvalue$chrName, CNA.segmented.df.pvalue$loc.start, CNA.segmented.df.pvalue$loc.end)
gene_symbol_and_cytobands.transposed <- data.frame(t(gene_symbol_and_cytobands), stringsAsFactors=FALSE)

# Merge gene_symbols_and_cytobands information with segments
CNA.segmented.df.pvalue.annotated <- merge(CNA.segmented.df.pvalue, gene_symbol_and_cytobands.transposed, by=c("chrName", "loc.start", "loc.end"), all=TRUE)

####################################################
### PREPARE DATA AND WRITE TO FILE (JSON)
####################################################

# Change "ID" column to be the base_output input
CNA.segmented.df.pvalue.annotated$ID <- base_output
CNA.segmented.df.pvalue.annotated$ploidy.rounded <- round(CNA.segmented.df.pvalue.annotated$ploidy, digits = 0)
CNA.segmented.df.pvalue.annotated$length <- as.integer(CNA.segmented.df.pvalue.annotated$loc.end) - as.integer(CNA.segmented.df.pvalue.annotated$loc.start)

# Write ISCN nomenclature
CNA.segmented.df.pvalue.annotated$iscn <- paste(sub("^chr", "", CNA.segmented.df.pvalue.annotated$chrName),
                                                CNA.segmented.df.pvalue.annotated$cytoband.start,
                                                CNA.segmented.df.pvalue.annotated$cytoband.end,
                                                "(",
                                                CNA.segmented.df.pvalue.annotated$loc.start,
                                                "-",
                                                CNA.segmented.df.pvalue.annotated$loc.end,
                                                ")x",
                                                CNA.segmented.df.pvalue.annotated$ploidy.rounded,
                                                sep = "")

# covert type to character so it prints nicely in the JSON
CNA.segmented.df.pvalue.annotated %>% mutate_if(is.list, as.character) -> CNA.segmented.df.pvalue.annotated

# Clean up colnames to not include "." character.  This is often used for nesting in JSON
colClean <- function(x){ colnames(x) <- gsub("\\.", "_", colnames(x)); x } 
CNA.segmented.df.pvalue.annotated <- colClean(CNA.segmented.df.pvalue.annotated)

# Write segments out to file
write(toJSON(CNA.segmented.df.pvalue.annotated), file = paste0(base_output, ".cbs.seg.json", collapse = ""))

####################################################
### PLOT CBS
####################################################

# Annotate CNA object with genesymbol overlap, then convert back to dataframe for plotting
wh <- sort(subsetByOverlaps(genesymbol, CNA.segmented.df.pvalue.genomicrange, ignore.strand=TRUE), ignore.strand=TRUE)
if (length(target_gene_list) > 0 & !is.na(target_gene_list)){
  wh <- subset(wh, wh$symbol %in% target_gene_list)
}
# Collapse overlapping transcripts for each gene
wh <- unlist(range(split(wh, ~symbol)))
wh <- data.frame(chrName=seqnames(wh),
                 start=start(wh),
                 ends=end(wh),
                 symbol=wh@ranges@NAMES)
# Create gene midpoint
wh$mid <- as.integer(round((wh$start + wh$end)/2))

# Create CBS plot 
plot_cbs <- function(CNA.segmented.df.pvalue.annotated, diffCoverage, wh){
  
  cbs_plot <- ggplot() +
    geom_segment(aes(x=as.integer(loc.start),
                     xend=as.integer(loc.end),
                     y=as.integer(round(ploidy, digits=0)),
                     yend=as.integer(round(ploidy, digits=0)), 
                     color="blue"), CNA.segmented.df.pvalue.annotated) +
    geom_point(aes(as.integer(begin),
                   as.numeric(difference),
                   color=ploidycutoff), diffCoverage) + 
    facet_wrap(~ chrName, 
               nrow = 4, 
               scales = "free_x") + 
    ylim(-3, max(diffCoverage$difference)+1) +
    scale_color_manual(values = c("black", "blue","red")) +
    ylab("ploidy") +
    theme(legend.position="none",
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank())
  
  # if number of genes > 500, then plot might look messy, so don't plot gene names
  if (nrow(wh) < 500){
    cbs_plot <- cbs_plot + geom_text(data = wh, 
                                     aes(x = mid,
                                         y = -2,
                                         label = symbol, 
                                         angle = 90),
                                     size=2, 
                                     position=position_dodge(width=0.5) )
  }
  
  return(cbs_plot)
}

CNA.segmented.df.pvalue$chrName <- factor(CNA.segmented.df.pvalue$chrName, levels=chrNames)

cbs_plot <- plot_cbs(CNA.segmented.df.pvalue, diffCoverage, wh)

ggsave(filename= paste0(base_output, ".grandlinear.cbs.segments.png", collapse = ""))

###########################################################
### PLOT GRANDLINEAR (only covered regions of chromosomes)
###########################################################


plot_hmm <- function(cnvseg, diffCoverage, wh){
  
  hmm_plot <- ggplot() +
    geom_segment(aes(x=as.integer(begin),
                     xend=as.integer(finish),
                     y=as.integer(round(ploidy, digits=0)),
                     yend=as.integer(round(ploidy, digits=0)), 
                     color="blue"), cnvseg) +
    geom_point(aes(as.integer(begin),
                   as.numeric(difference),
                   color=ploidycutoff), diffCoverage) + 
    facet_wrap(~ chrName, 
               nrow = 4, 
               scales = "free_x") + 
    ylim(-3, max(diffCoverage$difference)+1) +
    scale_color_manual(values = c("black", "blue","red")) +
    ylab("ploidy") +
    theme(legend.position="none",
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank())
  
  # if number of genes > 500, then plot might look messy, so don't plot gene names
  if (nrow(wh) < 500){
    hmm_plot <- hmm_plot + geom_text(data = subset(wh, wh$symbol %in% target_gene_list.cnv_only), 
                                     aes(x = mid,
                                         y = -2,
                                         label = symbol, 
                                         angle = 90),
                                     size=2, 
                                     position=position_dodge(width=0.5) )
  }
  
  return(hmm_plot)
}

hmm_plot <- plot_hmm(cnvseg, diffCoverage, wh)

ggsave(filename= paste0(base_output, ".grandlinear.hmm.segments.png", collapse = ""))



#################################################
### CONVERT IONREPORTER SEGMENT DATA TO GENOMICRANGES OBJECT
#################################################

# Load additional libraries
# ggbio may override ggplot2 
#detach("package:biomaRt")
library(ggbio)


# Transform dataframe to GenomicRange object
diffCoverage.genomicrange <- transformDfToGr(diffCoverage, seqnames = "chrName", start = "begin", end="finish")
cnvseg.genomicrange <- transformDfToGr(cnvseg, seqnames = "chrName", start = "begin", end="finish")

# Keep in specific order and chromosomes
chrNames <- as.vector(unique(diffCoverage$chrName))
diffCoverage.genomicrange <- keepSeqlevels(diffCoverage.genomicrange, chrNames)
seqlengths(diffCoverage.genomicrange) <- as.numeric(seqlengths(ideoCyto$hg19)[chrNames])
chrNames <- as.vector(unique(cnvseg$chrName))
cnvseg.genomicrange <- keepSeqlevels(cnvseg.genomicrange, chrNames)
seqlengths(cnvseg.genomicrange) <- as.numeric(seqlengths(ideoCyto$hg19)[chrNames])

# Establish ploidy cutoff
diffCoverage.genomicrange$ploidycutoff <- cut(diffCoverage.genomicrange$difference,
                                              breaks = c(-Inf, 4, Inf),
                                              labels = c("<=4",">4")
)

#########################################
### PLOT GRANDLINEAR (entire chromosomes)
#########################################

# grandlinear <- autoplot(diffCoverage.genomicrange, geom = "point", coord = "genome", aes(y = difference, color=ploidycutoff)) +
#   scale_color_manual(values = c("black", "blue")) +
#   ylab("ploidy") +
#   scale_y_continuous(breaks = seq(0, 10, by = 2), limits = c(0,10)) +
#   geom_hline(aes(yintercept=2), linetype="dashed") +
#   theme(legend.position="none")

# grandlinear <- autoplot(diffCoverage.genomicrange, geom = "point", coord = "genome", aes(y = difference, color=ploidycutoff)) +
#   scale_color_manual(values = c("black", "blue")) +
#   geom_segment(cnvseg.genomicrange, aes(x=begin,xend=finish,y=ploidy,yend=ploidy), size=2, coord = "genome", stat = "identity", color="red") +
#   ylab("ploidy") +
#   scale_y_continuous(breaks = seq(0, 10, by = 2)) +
#   ylim(c(0,10)) +
#   geom_hline(aes(yintercept=2), linetype="dashed") +
#   theme(legend.position="none") +
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())


grandlinear <- plotGrandLinear(diffCoverage.genomicrange, aes(y = difference, color=ploidycutoff), spaceline = TRUE, ylim = c(0,10)) +
                            scale_color_manual(values = c("black", "blue")) +
                            geom_hline(aes(yintercept=2), linetype="dashed")

ggsave(filename= paste0(base_output, ".grandlinear.png", collapse = ""), height=8, width=8)


########################################
### COMPARE GENE BOUNDARIES WITH SEGMENTS
########################################

compare_genes_with_segments <- function(genomerange1, genomicrange2, seq, ggplot_item){
  
  if(head(start(sort(subsetByOverlaps(genomerange1, genomicrange2[seq], ignore.strand=TRUE), ignore.strand=TRUE)), n=1) == start(genomicrange2[seq])){
    
    ggplot_item <- ggplot_item + geom_vline(xintercept = head(start(sort(subsetByOverlaps(genomerange1, genomicrange2[seq], ignore.strand=TRUE), ignore.strand=TRUE)), n=1),
                                            linetype="dotted",
                                            color="green")
    
  } else{
    ggplot_item <- ggplot_item + geom_vline(xintercept = c(head(start(sort(subsetByOverlaps(genomerange1, genomicrange2[seq], ignore.strand=TRUE), ignore.strand=TRUE)), n=1),
                                                           start(genomicrange2[seq])), 
                                            linetype="dotted", 
                                            color=factor(c("blue", "red"))
    )
    
  }
  if(tail(end(sort(subsetByOverlaps(genomerange1, genomicrange2[seq], ignore.strand=TRUE), ignore.strand=TRUE)), n=1) == end(genomicrange2[seq])){
    
    ggplot_item <- ggplot_item + geom_vline(xintercept = tail(end(sort(subsetByOverlaps(genomerange1, genomicrange2[seq], ignore.strand=TRUE), ignore.strand=TRUE)), n=1),
                                            linetype="dotted",
                                            color="green")
    
  } else{
    ggplot_item <- ggplot_item + geom_vline(xintercept = c(tail(end(sort(subsetByOverlaps(genomerange1, genomicrange2[seq], ignore.strand=TRUE), ignore.strand=TRUE)), n=1),
                                                           end(genomicrange2[seq])), 
                                            linetype="dotted", 
                                            color=factor(c("blue", "red"))
    )
    
  }
  return(ggplot_item)
}



########################################
### LOOP OVER EACH SEGMENT AND PLOT IF MEETS COVERAGE AND PLOIDY THRESHOLDS
########################################

for (seq in seq_along(cnvseg.genomicrange)){
  if (cnvseg.genomicrange[seq]$ploidy >= 4 & cnvseg.genomicrange[seq]$confidence >= confidence_cutoff){
    wh.segment <- sort(subsetByOverlaps(genesymbol, cnvseg.genomicrange[seq], ignore.strand=TRUE), ignore.strand=TRUE)

    # Fix gene aliases
    wh.segment.substitutes <- subset(wh.segment$symbol, !(wh.segment$symbol %in% aliasSymbol$symbol) & wh.segment$symbol %in% aliasSymbol$alias_symbol)
    wh.segment$symbol <- unlist(replace(wh.segment$symbol, which(wh.segment$symbol %in% wh.segment.substitutes), values=lapply(wh.segment.substitutes, FUN=checkGeneAlias)))

    for (gene in unique(wh.segment$symbol)){
      if (length(target_gene_list) > 0){
        if (gene %in% target_gene_list){

          #wh.gene <- genesymbol[gene]
          # 
          # if (length(target_gene_list) > 0 & !is.na(target_gene_list)){
          #   wh <- subset(wh, wh$symbol %in% target_gene_list)
          # }
          #     ideogram <- Ideogram(genome = "hg19", xlabel=TRUE, color="red")
          #     ideogram <- ideogram + xlim(wh)

          # Query Ensembl BioMart for gene boundaries and convert to GR object
          wh.gene <- getGeneBoundaries(gene)
          wh.gene$chromosome_name <- paste("chr", wh.gene$chromosome_name, sep="")
          wh.gene <- subset(wh.gene, wh.gene$chromosome_name %in% chrNames)
          wh.gene <- transformDfToGr(wh.gene, seqnames = "chromosome_name", start = "start_position", end="end_position")

          transcripts <- autoplot(Homo.sapiens, which=wh.gene) + xlim(wh.gene)
          
          target_regions.plot <- autoplot(sort(subsetByOverlaps(target_regions.genomicrange, wh.gene, ignore.strand=TRUE), ignore.strand=TRUE),
                                          geom="segment", size=4, aes(y=0), stat = "identity", color="purple") + 
                                          scale_x_sequnit("bp") +
                                          theme(axis.text.y=element_blank(),
                                                axis.title.y=element_blank(),
                                                axis.ticks.y=element_blank())
          #     diffCoverage.genomicrange.subrange <- autoplot(sort(subsetByOverlaps(cnvseg.genomicrange, wh.gene, ignore.strand=TRUE), ignore.strand=TRUE),
          #                                     geom="segment", size=1, aes(y=ploidy), stat = "identity") + scale_x_sequnit("bp") +
          #                                     ylim(c(0,max(cnvseg.genomicrange$ploidy) + 2)) +
          #                                     geom_point(data=subset(diffCoverage, cnvseg.genomicrange[seq]$chrName == diffCoverage$chrName), aes(x=begin,y=difference)) +
          #                                     xlim(wh)
          #print(sort(subsetByOverlaps(cnvseg.genomicrange, wh.gene, ignore.strand=TRUE), ignore.strand=TRUE))
          diffCoverage$begin
          diffCoverage.genomicrange.subrange <- autoplot(sort(subsetByOverlaps(cnvseg.genomicrange, wh.gene, ignore.strand=TRUE), ignore.strand=TRUE),
                                                         geom="segment", size=1, aes(y=ploidy), stat = "identity") + scale_x_sequnit("bp") +
            ylim(c(0,max(cnvseg.genomicrange$ploidy)*1.5)) +
            geom_point(data=subset(diffCoverage, (min(wh.gene@ranges@start) <= diffCoverage$begin & max(wh.gene@ranges@start + wh.gene@ranges@width) >= diffCoverage$finish) & cnvseg.genomicrange[seq]$chrName == diffCoverage$chrName), aes(x=begin,y=difference)) +
            xlim(wh.gene)

          # diffCoverage.genomicrange.subrange <- autoplot(sort(subsetByOverlaps(cnvseg.genomicrange, wh, ignore.strand=TRUE), ignore.strand=TRUE),
          #                                                geom="segment", size=1, aes(y=ploidy), stat = "identity") + scale_x_sequnit("bp") +
          #   ylim(c(0,max(cnvseg.genomicrange$ploidy) + 2)) +
          #   geom_point(data=subset(diffCoverage, cnvseg.genomicrange[seq]$chrName == diffCoverage$chrName), aes(x=begin,y=difference)) +
          #   xlim(wh)
          
          if (length(subsetByOverlaps(target_regions.genomicrange, wh.gene, ignore.strand=TRUE)) > 0){
              tracks("Target\nRegions"=target_regions.plot, 
                     "Copy Number Segments"=diffCoverage.genomicrange.subrange, 
                     "Transcripts\n(all genes)"=transcripts, 
                     heights=c(1,3,2.5))
          } else {
              tracks("Copy Number Segments"=diffCoverage.genomicrange.subrange, 
                     "Transcripts\n(all genes)"=transcripts, 
                     heights=c(1,3,2.5))
          }
              
          
          file_name <- c(gene, ".hmm", ".png")
          file_name <- paste0(file_name, collapse="")
          ggsave(filename=file_name, height=8, width=8)
        }
      }
    }
    
    if (length(target_gene_list) > 0 & !is.na(target_gene_list)){
      wh.segment <- subset(wh.segment, wh.segment$symbol %in% target_gene_list)
    }

    #     ideogram <- Ideogram(genome = "hg19", xlabel=TRUE, color="red")
    #     ideogram <- ideogram + xlim(wh.segment)
    
    targeted_genes.plot <- ggplot(wh.segment) + geom_segment(stat="identity", size=6, aes(y=0, color=symbol)) + 
      scale_x_sequnit("bp") +
      theme(axis.text.y=element_blank(),
            axis.title.y=element_blank(),
            axis.ticks.y=element_blank(),
            legend.title=element_blank(),
            legend.position="top")
    targeted_genes.plot <- compare_genes_with_segments(wh.segment, cnvseg.genomicrange, seq, targeted_genes.plot) 
    
    transcripts <- autoplot(Homo.sapiens, which=wh.segment) + xlim(wh.segment)
    transcripts <- compare_genes_with_segments(wh.segment, cnvseg.genomicrange, seq, transcripts)
    
    target_regions.plot <- autoplot(sort(subsetByOverlaps(target_regions.genomicrange, wh.segment, ignore.strand=TRUE), ignore.strand=TRUE),
                                    geom="segment", size=4, aes(y=0), stat = "identity", color="purple") + 
      scale_x_sequnit("bp") +
      theme(axis.text.y=element_blank(),
            axis.title.y=element_blank(),
            axis.ticks.y=element_blank())
    target_regions.plot <- compare_genes_with_segments(wh.segment, cnvseg.genomicrange, seq, target_regions.plot)
    
    #     diffCoverage.genomicrange.subrange <- autoplot(sort(subsetByOverlaps(cnvseg.genomicrange, wh.segment, ignore.strand=TRUE), ignore.strand=TRUE),
    #                                     geom="segment", size=1, aes(y=ploidy), stat = "identity") + scale_x_sequnit("bp") +
    #                                     ylim(c(0,max(cnvseg.genomicrange$ploidy) + 2)) +
    #                                     geom_point(data=subset(diffCoverage, cnvseg.genomicrange[seq]$chrName == diffCoverage$chrName), aes(x=begin,y=difference)) +
    #                                     xlim(wh.segment)

    diffCoverage.genomicrange.subrange <- autoplot(sort(subsetByOverlaps(cnvseg.genomicrange, wh.segment, ignore.strand=TRUE), ignore.strand=TRUE),
                                                   geom="segment", size=1, aes(y=ploidy), stat = "identity") + scale_x_sequnit("bp") +
      ylim(c(0,max(CNA.segmented.df.pvalue.genomicrange$ploidy)*1.5)) +
      geom_point(data=subset(diffCoverage, (min(wh.segment@ranges@start) <= diffCoverage$begin & max(wh.segment@ranges@start + wh.segment@ranges@width) >= diffCoverage$finish) & cnvseg.genomicrange[seq]$chrName == diffCoverage$chrName), aes(x=begin,y=difference)) +
      xlim(wh.segment)
    
    # diffCoverage.genomicrange.subrange <- autoplot(sort(subsetByOverlaps(cnvseg.genomicrange, wh, ignore.strand=TRUE), ignore.strand=TRUE),
    #                                                geom="segment", size=1, aes(y=ploidy), stat = "identity") + scale_x_sequnit("bp") +
    #   ylim(c(0,max(cnvseg.genomicrange$ploidy) + 2)) +
    #   geom_point(data=subset(diffCoverage, cnvseg.genomicrange[seq]$chrName == diffCoverage$chrName), aes(x=begin,y=difference)) +
    #   xlim(wh)
    diffCoverage.genomicrange.subrange <- compare_genes_with_segments(wh.segment, cnvseg.genomicrange, seq, diffCoverage.genomicrange.subrange)
    
    tracks("Target\nRegions"=target_regions.plot, 
           "Targeted\nGenes"=targeted_genes.plot, 
           "Copy Number Segments"=diffCoverage.genomicrange.subrange, 
           "Transcripts\n(all genes)"=transcripts, 
           heights=c(1,2.5,3,2.5))
    
    file_name <- c("chr", cnvseg.genomicrange[seq]$chrName, "-", cnvseg.genomicrange[seq]$begin, "-", cnvseg.genomicrange[seq]$finish, ".hmm", ".png")
    file_name <- paste0(file_name, collapse="")
    #ggsave(filename=file_name)
    ggsave(filename=file_name, height=8, width=8)
  }
}


for (seq in seq_along(CNA.segmented.df.pvalue.genomicrange)){
  if (as.integer(CNA.segmented.df.pvalue.genomicrange[seq]$ploidy) >= 4){
    wh.segment <- sort(subsetByOverlaps(genesymbol, CNA.segmented.df.pvalue.genomicrange[seq], ignore.strand=TRUE), ignore.strand=TRUE)
    
    # Fix gene aliases
    wh.segment.substitutes <- subset(wh.segment$symbol, !(wh.segment$symbol %in% aliasSymbol$symbol) & wh.segment$symbol %in% aliasSymbol$alias_symbol)
    wh.segment$symbol <- unlist(replace(wh.segment$symbol, which(wh.segment$symbol %in% wh.segment.substitutes), values=lapply(wh.segment.substitutes, FUN=checkGeneAlias)))
    
    for (gene in unique(wh.segment$symbol)){
      if (length(target_gene_list) > 0){
        if (gene %in% target_gene_list){
          
          # Query Ensembl BioMart for gene boundaries and convert to GR object
          wh.gene <- getGeneBoundaries(gene)
          wh.gene$chromosome_name <- paste("chr", wh.gene$chromosome_name, sep="")
          wh.gene <- subset(wh.gene, wh.gene$chromosome_name %in% chrNames)
          wh.gene <- transformDfToGr(wh.gene, seqnames = "chromosome_name", start = "start_position", end="end_position")

          # 
          # if (length(target_gene_list) > 0 & !is.na(target_gene_list)){
          #   wh <- subset(wh, wh$symbol %in% target_gene_list)
          # }
          #     ideogram <- Ideogram(genome = "hg19", xlabel=TRUE, color="red")
          #     ideogram <- ideogram + xlim(wh)

          transcripts <- autoplot(Homo.sapiens, which=wh.gene) + xlim(wh.gene)

          target_regions.plot <- autoplot(sort(subsetByOverlaps(target_regions.genomicrange, wh.gene, ignore.strand=TRUE), ignore.strand=TRUE),
                                          geom="segment", size=4, aes(y=0), stat = "identity", color="purple") + 
            scale_x_sequnit("bp") +
            theme(axis.text.y=element_blank(),
                  axis.title.y=element_blank(),
                  axis.ticks.y=element_blank())

          #     diffCoverage.genomicrange.subrange <- autoplot(sort(subsetByOverlaps(cnvseg.genomicrange, wh.gene, ignore.strand=TRUE), ignore.strand=TRUE),
          #                                     geom="segment", size=1, aes(y=ploidy), stat = "identity") + scale_x_sequnit("bp") +
          #                                     ylim(c(0,max(cnvseg.genomicrange$ploidy) + 2)) +
          #                                     geom_point(data=subset(diffCoverage, cnvseg.genomicrange[seq]$chrName == diffCoverage$chrName), aes(x=begin,y=difference)) +
          #                                     xlim(wh)
          #print(sort(subsetByOverlaps(CNA.segmented.df.pvalue.genomicrange, wh.gene, ignore.strand=TRUE), ignore.strand=TRUE))
          diffCoverage.genomicrange.subrange <- autoplot(sort(subsetByOverlaps(CNA.segmented.df.pvalue.genomicrange, wh.gene, ignore.strand=TRUE), ignore.strand=TRUE),
                                                         geom="segment", size=1, aes(y=ploidy), stat = "identity") + scale_x_sequnit("bp") +
            ylim(c(0,max(CNA.segmented.df.pvalue.genomicrange$ploidy)*1.5)) +
            geom_point(data=subset(diffCoverage, (min(wh.gene@ranges@start) <= diffCoverage$begin & max(wh.gene@ranges@start + wh.gene@ranges@width) >= diffCoverage$finish) & CNA.segmented.df.pvalue.genomicrange[seq]$chrName == diffCoverage$chrName), aes(x=begin,y=difference)) +
            xlim(wh.gene)
          
          # diffCoverage.genomicrange.subrange <- autoplot(sort(subsetByOverlaps(cnvseg.genomicrange, wh, ignore.strand=TRUE), ignore.strand=TRUE),
          #                                                geom="segment", size=1, aes(y=ploidy), stat = "identity") + scale_x_sequnit("bp") +
          #   ylim(c(0,max(cnvseg.genomicrange$ploidy) + 2)) +
          #   geom_point(data=subset(diffCoverage, cnvseg.genomicrange[seq]$chrName == diffCoverage$chrName), aes(x=begin,y=difference)) +
          #   xlim(wh)

          tracks("Target\nRegions"=target_regions.plot, 
                 "Copy Number Segments"=diffCoverage.genomicrange.subrange, 
                 "Transcripts\n(all genes)"=transcripts, 
                 heights=c(1,3,2.5))
          
          file_name <- c(gene, ".cbs", ".png")
          file_name <- paste0(file_name, collapse="")
          ggsave(filename=file_name, height=8, width=8)
        }
      }
    }

    if (length(target_gene_list) > 0 & !is.na(target_gene_list)){
      wh.segment <- subset(wh.segment, wh.segment$symbol %in% target_gene_list)
    }
    #     ideogram <- Ideogram(genome = "hg19", xlabel=TRUE, color="red")
    #     ideogram <- ideogram + xlim(wh.segment)
    targeted_genes.plot <- ggplot(wh.segment) + geom_segment(stat="identity", size=6, aes(y=0, color=symbol)) + 
      scale_x_sequnit("bp") +
      theme(axis.text.y=element_blank(),
            axis.title.y=element_blank(),
            axis.ticks.y=element_blank(),
            legend.title=element_blank(),
            legend.position="top")
    targeted_genes.plot <- compare_genes_with_segments(wh.segment, CNA.segmented.df.pvalue.genomicrange, seq, targeted_genes.plot) 
    targeted_genes.plot
    
    transcripts <- autoplot(Homo.sapiens, which=wh.segment) + xlim(wh.segment)
    transcripts <- compare_genes_with_segments(wh.segment, CNA.segmented.df.pvalue.genomicrange, seq, transcripts)
    
    target_regions.plot <- autoplot(sort(subsetByOverlaps(target_regions.genomicrange, wh.segment, ignore.strand=TRUE), ignore.strand=TRUE),
                                    geom="segment", size=4, aes(y=0), stat = "identity", color="purple") + 
      scale_x_sequnit("bp") +
      theme(axis.text.y=element_blank(),
            axis.title.y=element_blank(),
            axis.ticks.y=element_blank())
    target_regions.plot <- compare_genes_with_segments(wh.segment, CNA.segmented.df.pvalue.genomicrange, seq, target_regions.plot)
    
    #     diffCoverage.genomicrange.subrange <- autoplot(sort(subsetByOverlaps(cnvseg.genomicrange, wh.segment, ignore.strand=TRUE), ignore.strand=TRUE),
    #                                     geom="segment", size=1, aes(y=ploidy), stat = "identity") + scale_x_sequnit("bp") +
    #                                     ylim(c(0,max(cnvseg.genomicrange$ploidy) + 2)) +
    #                                     geom_point(data=subset(diffCoverage, cnvseg.genomicrange[seq]$chrName == diffCoverage$chrName), aes(x=begin,y=difference)) +
    #                                     xlim(wh.segment)
    diffCoverage.genomicrange.subrange <- autoplot(sort(subsetByOverlaps(CNA.segmented.df.pvalue.genomicrange, wh.segment, ignore.strand=TRUE), ignore.strand=TRUE),
                                                   geom="segment", size=1, aes(y=ploidy), stat = "identity") + scale_x_sequnit("bp") +
      ylim(c(0,max(CNA.segmented.df.pvalue.genomicrange$ploidy)*1.5)) +
      geom_point(data=subset(diffCoverage, (min(wh.segment@ranges@start) <= diffCoverage$begin & max(wh.segment@ranges@start + wh.segment@ranges@width) >= diffCoverage$finish) & CNA.segmented.df.pvalue.genomicrange[seq]$chrName == diffCoverage$chrName), aes(x=begin,y=difference)) +
      xlim(wh.segment)
    
    # diffCoverage.genomicrange.subrange <- autoplot(sort(subsetByOverlaps(cnvseg.genomicrange, wh, ignore.strand=TRUE), ignore.strand=TRUE),
    #                                                geom="segment", size=1, aes(y=ploidy), stat = "identity") + scale_x_sequnit("bp") +
    #   ylim(c(0,max(cnvseg.genomicrange$ploidy) + 2)) +
    #   geom_point(data=subset(diffCoverage, cnvseg.genomicrange[seq]$chrName == diffCoverage$chrName), aes(x=begin,y=difference)) +
    #   xlim(wh)
    diffCoverage.genomicrange.subrange <- compare_genes_with_segments(wh.segment, CNA.segmented.df.pvalue.genomicrange, seq, diffCoverage.genomicrange.subrange)
    
    tracks("Target\nRegions"=target_regions.plot, 
           "Targeted\nGenes"=targeted_genes.plot, 
           "Copy Number Segments"=diffCoverage.genomicrange.subrange, 
           "Transcripts\n(all genes)"=transcripts, 
           heights=c(1,2.5,3,2.5))
    
    file_name <- c(CNA.segmented.df.pvalue.genomicrange[seq]$chrName, "-", CNA.segmented.df.pvalue.genomicrange[seq]$loc.start, "-", CNA.segmented.df.pvalue.genomicrange[seq]$loc.end, ".cbs", ".png")
    file_name <- paste0(file_name, collapse="")
    ggsave(filename=file_name, height=8, width=8)
  }
}

