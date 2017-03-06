plotMotifStats <- function(
        PlotDat, BaseDat, StatsFDat, StatsDat, PlotType, QuantWidth){
    ylim <- 3.5
    regions <- unique(PlotDat$Region)
    midReg <- regions[(length(regions) + 1) / 2]
    ps <- lapply(regions, function(region){
        rBaseDat <- BaseDat[BaseDat$Region==region,]
        rPlotDat <- PlotDat[PlotDat$Region==region,]
        if(PlotType %in% c('Signal', 'Downsample')){
            ## randomize so all of one group isn't on top
            sRPlotDat <- split(rPlotDat, rPlotDat$Read)
            rPlotDat <- do.call(rbind.data.frame, sample(sRPlotDat))
            p <- ggplot(rPlotDat) + geom_path(
                      aes(x=Position, y=Signal, color=Group, group=Read),
                      alpha=0.5, size=0.1, show.legend=FALSE)
        } else if(PlotType == 'Quantile'){
            p <- ggplot(rPlotDat) +
                geom_rect(aes(xmin=Position, xmax=Position + QuantWidth,
                              ymin=Lower, ymax=Upper, fill=Group),
                          alpha=0.2, show.legend=FALSE) +
                scale_fill_manual(
                    values=c('Group1'='blue', 'Group2'='red')) +
                ylab('Signal')
        } else if(PlotType == 'Boxplot'){
            p <- ggplot(rPlotDat) +
                geom_boxplot(
                    aes(Position + 0.5, ymin=SigMin, lower=Sig25,
                        middle=SigMed, upper=Sig75, ymax=SigMax,
                        fill=Group), size=0.2, alpha=1,
                    stat="identity", show.legend=FALSE) +
                scale_fill_manual(
                    values=c('Group1'='blue', 'Group2'='red')) +
                ylab('Signal')
        } else if(PlotType == 'Violin'){
            p <- ggplot(rPlotDat) +
                geom_violin(aes(x=Position + 0.5, y=Signal, fill=Group,
                                group=paste0(Group, Position)),
                            size=0, show.legend=FALSE) +
                scale_fill_manual(
                    values=c('Group1'='blue', 'Group2'='red')) +
                ylab('Signal')
        }
        p <- p + geom_text(aes(x=Position+0.5, y=-ylim,
                          label=Base, color=Base),
                      data=rBaseDat,
                      hjust=0.5, vjust=0, size=3, show.legend=FALSE) +
            scale_color_manual(
                values=c(
                    'A'='#00CC00', 'C'='#0000CC', 'G'='#FFB300',
                    'T'='#CC0000', '-'='black', 'N'='black',
                    'Group1'='red', 'Group2'='black')) +
            geom_vline(
                xintercept=min(rBaseDat$Position):
                (max(rBaseDat$Position) + 1), size=0.01) +
           scale_x_continuous(expand=c(0,0)) +
           coord_cartesian(ylim=c(-ylim, ylim)) +
           theme_bw() +
           theme(axis.text.x=element_blank(),
                 axis.text.y=element_blank(),
                 axis.title.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 axis.ticks.y=element_blank(),
                 plot.margin=margin(0,0,0,0,'lines'))
        if(region != midReg){
            p <- p + theme(axis.title.y=element_blank())
        }
        return(p)
    })

    maxStat <- max(StatsDat$NegLogPValue)
    if(maxStat < 1){ tickVals <- c(0,0.2,0.4,0.6,0.8,1)
    } else if(maxStat < 10){ tickVals <- seq(0,10,by=2)
    } else { tickVals <- seq(0,100,by=5) }
    ps[[length(ps) + 1]] <- ggplot(StatsDat) +
        geom_violin(aes(
            x=Position+0.5, y=NegLogPValue,
            group=cut_width(Position, 0.9999)), size=0.1, fill='black') +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(breaks=tickVals) +
        theme_bw() +
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.x=element_blank())
    maxStatF <- max(StatsFDat$NegLogFishersPValue)
    if(maxStatF < 1){ tickVals <- c(0,0.2,0.4,0.6,0.8,1)
    } else if(maxStatF < 10){ tickVals <- seq(0,10,by=2)
    } else { tickVals <- seq(0,100,by=5) }
    ps[[length(ps) + 1]] <- ggplot(StatsFDat) +
        geom_violin(aes(
            x=Position+0.5, y=NegLogFishersPValue,
            group=cut_width(Position, 0.9999)), size=0.1, fill='black') +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(breaks=tickVals) +
        xlab('Position') + theme_bw() +
        theme(axis.text.x=element_text(hjust=0))
    print(do.call(
        plot_grid,
        c(ps, list(ncol=1, align='v',
                   rel_heights=c(rep(1, length(regions)), 2, 3)))))
}
