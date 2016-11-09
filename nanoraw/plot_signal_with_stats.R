suppressPackageStartupMessages( library(dplyr) )
suppressPackageStartupMessages( library(stringr) )
suppressPackageStartupMessages( library(ggplot2) )
suppressPackageStartupMessages( library(cowplot) )

load('region.signal_data.RData')
load('region.base_data.RData')
statsDat <- read.csv('region.stats.csv')

compBases <- c('A'='T', 'C'='G','G'='C', 'T'='A')
revComp <- function(seq){
    str_c(compBases[rev(str_split(seq, "")[[1]])], collapse="")
}

motif <- "CCTGG"
poss <- c(14, 16)
expandBases <- 4

statsDat <- statsDat %>% mutate(Region=str_pad(Region, 3, pad="0"))
baseDat <- baseDat %>% mutate(Region=as.character(Region))
sigDat <- sigDat %>% mutate(Region=as.character(Region))
baseDat <- left_join(baseDat, statsDat, by=c('Region', 'Position'))

regStrandsTmp <- baseDat %>% group_by(Region) %>%
    summarize(Strand=first(Strand))
regStrands <- regStrandsTmp$Strand
names(regStrands) <- regStrandsTmp$Region

## find regions with motif of interest
regionSeq <- baseDat %>% group_by(Region, Strand) %>%
    summarize(seq=str_c(Base, collapse=""))
regionSeq <- regionSeq %>%
    mutate(seq=ifelse(Strand == '+', seq, revComp(seq)))
regionSeq <- regionSeq %>%
    mutate(motifPos=str_locate(seq, motif)[,'start'])
containsRegs <- (regionSeq %>% filter(motifPos %in% poss))$Region

motifBaseDat <- baseDat %>% filter(Region %in% containsRegs)
motifSigDat <- sigDat %>%
    filter(Region %in% containsRegs,
           ifelse(Strand == "Forward Strand", '+', '-') ==
           regStrands[Region])

## do just forward strand first
plusBaseDat <- motifBaseDat %>% filter(Strand == '+')
plusSigDat <- motifSigDat %>% filter(Strand == 'Forward Strand')

sPlusBaseDat <- split(plusBaseDat, plusBaseDat$Region)
sPlusSigDat <- split(plusSigDat, plusSigDat$Region)

regsBaseDat <- list()
regsSigDat <- list()
for (region in unique(plusBaseDat$Region)){
    regBaseDat <- sPlusBaseDat[[region]]
    regSigDat <- sPlusSigDat[[region]]
    regMotifStart <- (regionSeq %>% filter(Region==region))$motifPos
    regStart <- min(regBaseDat$Position)
    subRegBaseDat <- regBaseDat %>%
        filter(Position >= regStart + regMotifStart - 1 - expandBases,
               Position <= regStart + regMotifStart +
               nchar(motif) + expandBases - 2)
    subRegSigDat <- regSigDat %>%
        filter(Position >= regStart + regMotifStart - 1 - expandBases,
               Position < regStart + regMotifStart +
               nchar(motif) + expandBases - 1)
    subRegBaseDat$Position <- subRegBaseDat$Position - regStart -
        regMotifStart + expandBases + 1
    subRegSigDat$Position <- subRegSigDat$Position - regStart -
        regMotifStart + expandBases + 1
    regsBaseDat[[length(regsBaseDat) + 1]] <- subRegBaseDat
    regsSigDat[[length(regsSigDat) + 1]] <- subRegSigDat
}

regsBaseDat <- do.call(rbind.data.frame, regsBaseDat)
regsSigDat <- do.call(rbind.data.frame, regsSigDat)

save(regsBaseDat, file='regsBaseDat.RData')
save(regsSigDat, file='regsSigDat.RData')

plotRegs <- unique(regsBaseDat$Region)[1:5]

subRegsBaseDat <- regsBaseDat %>% filter(Region %in% plotRegs)
subRegsSigDat <- regsSigDat %>% filter(Region %in% plotRegs)

ps <- lapply(plotRegs, function(region){
    regBaseDat <- subRegsBaseDat %>% filter(Region==region)
    regSigDat <- subRegsSigDat %>% filter(Region==region)
    (ggplot(regSigDat) +
     geom_path(aes(x=Position, y=Signal, color=Group, group=Read),
               alpha=0.3, size=0.05, show.legend=FALSE) +
     geom_text(aes(x=Position+0.5, y=-5,
                   label=Base, color=Base),
               data=regBaseDat,
               hjust=0.5, vjust=0, size=3, show.legend=FALSE) +
     scale_color_manual(
         values=c(
             'A'='#00CC00', 'C'='#0000CC', 'G'='#FFB300',
             'T'='#CC0000', '-'='black', 'N'='black',
             'Group1'='blue', 'Group2'='red')) +
     geom_vline(
         xintercept=min(regBaseDat$Position):
         (max(regBaseDat$Position) + 1),
         size=0.01) +
     theme_bw() +
     theme(axis.text.x=element_blank(), axis.text.y=element_blank(),
           axis.title.x=element_blank()))
})

ps[[length(ps) + 1]] <- ggplot(regsBaseDat) +
    geom_boxplot(aes(x=Position+0.5, y=Score,
                     group=cut_width(Position, 0.9999),
                     closed='right')) +
    scale_x_continuous()
    xlab('Position') + theme_bw() +
    theme(axis.text.x=element_text(hjust=0))

pdf('motif_centered_signal.w_test_scores.pdf', height=5)
do.call(plot_grid, c(ps, list(ncol=1, align='v', rel_heights=c(rep(1, 5), 2))))
foo <- dev.off()



pdf('motif_centered_signal.pdf', height=8)
ggplot(subRegsSigDat) +
    geom_path(aes(x=Position, y=Signal, color=Group, group=Read),
              alpha=0.3, size=0.05, show.legend=FALSE) +
    facet_grid(Region ~ .) +
    geom_text(aes(x=Position+0.5, y=-5,
                  label=Base, color=Base),
              data=subRegsBaseDat,
              hjust=0.5, vjust=0, size=3, show.legend=FALSE) +
              scale_color_manual(
                  values=c(
                      'A'='#00CC00', 'C'='#0000CC', 'G'='#FFB300',
                      'T'='#CC0000', '-'='black', 'N'='black',
                      'Group1'='blue', 'Group2'='red')) +
    geom_vline(
        xintercept=min(subRegsBaseDat$Position):
        (max(subRegsBaseDat$Position) + 1),
        size=0.01) +
    theme_bw() + theme(axis.text.x=element_text(hjust=0))
foo <- dev.off()
