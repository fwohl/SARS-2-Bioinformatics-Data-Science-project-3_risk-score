plot.iris<-function(x, labels,
                   main="A UMAP visualization of the variants in Corona",
                   colors=c("#ff7f00", "#e377c2", "#17becf","#1EE629","#CE1A63","#132A67","#D10934" ),
                   pad=0.1, cex=0.6, pch=19, add=FALSE, legend.suffix="",
                   cex.main=1, cex.legend=0.85) {
          
            layout <- x
            if (is(x, "umap")) {
              layout <- x$layout
            } 
            
            xylim <- range(layout)
            xylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
            if (!add) {
              par(mar=c(0.2,0.7,1.2,0.7), ps=10)
              plot(xylim, xylim, type="n", axes=F, frame=F)
              rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd=0.25)  
            }
            points(layout[,1], layout[,2], col=colors[as.integer(labels)],
                   cex=cex, pch=pch)
            mtext(side=3, main, cex=cex.main)
          
            labels.u <- unique(labels)
            legend.pos <- "topleft"
            legend.text <- as.character(labels.u)
            if (add) {
              legend.pos <- "bottomleft"
              legend.text <- paste(as.character(labels.u), legend.suffix)
            }
          
            legend(legend.pos, legend=legend.text, inset=0.03,
                   col=colors[as.integer(labels.u)],
                   bty="n", pch=pch, cex=cex.legend)
          }
