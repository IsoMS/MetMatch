#' plotScores
#'
#' @param data input data
#' @param score.method score method
#' @param run.type run type
#'
#' @export
#'
plotScores = function(scored.data, run.type="decoy"){

   require(MASS,quietly = TRUE)
   require(RColorBrewer,quietly = TRUE)
   require(fields,quietly = TRUE)

   score.mass = scored.data$scored.data$score.mass.Q
   score.int = scored.data$scored.data$score.int.Q
   z = which(score.mass > 0 & score.int>0)

   kd <- kde2d(x=score.mass[z],
               y=score.int[z], n = 256,
               h = c(width.SJ(score.mass[z]),
                     width.SJ(score.int[z])))

   nbcol=128

   rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
   rc <- rf(nbcol-1)
   rc = c("#ffffff",rc)

   zcol  = cut(kd$z, nbcol)

   image.plot(kd, col=rc,xlab="Score Mass", ylab="Score Intensity")


   if(run.type == "decoy"){
      w = which(scored.data$dec.scored.data$hit.type == "decoy")
      dec.data = scored.data$dec.scored.data[w,]
      points(dec.data$score.mass.Q[w],dec.data$score.int.Q[w],cex=0.75,pch=3)
   }

   score=scored.data$scored.data$score.Q
   hist(score[z],n=40,freq=F,xlab="Score",main="")
   lines(density(score[z]),lwd=2,col=4)

   if(run.type=="decoy"){
      cutoffs = compute.thresholds(scored.data$dec.scored.data)
      abline(v=cutoffs, lty=2,lwd=2,col=2)
   }
}


#' plotDeltas
#'
#' @param Data input delta values in dataframe
#'
#' @export
#'

plotDeltas = function(delta.models,mass.tol){
   xm = seq(-mass.tol,mass.tol,by=2*mass.tol/200)
   xi = seq(-1,1,by=0.01)

   mod.mass = delta.models$mod.mass
   mod.int = delta.models$mod.int

   d.lap.mass = dlaplace(xm, mod.mass$mu, mod.mass$b)
   d.lap.int = dlaplace(xi, mod.int$mu, mod.int$b)

   h.mass = hist(delta.models$delta.mass,freq=F,n=100,main="",xlab = "delta mass",border="grey30")
   lines(xm, d.lap.mass*max(h.mass$density)/max(d.lap.mass),col=4, lwd=2,lty=1)
   abline(v=mod.mass$mu,lty=2, lwd=2,col=2)


   h.int = hist(delta.models$delta.int,freq=F,n=100,main="",xlab="delta Int",border="grey30",xlim=c(-0.5,0.5))
   lines(xi, d.lap.int*max(h.int$density)/max(d.lap.int),col=4, lwd=2,lty=1)
   abline(v=mod.int$mu,lty=2, lwd=2, col=2)

}


#' plotHeatMap
#'
#' @param Dataframe from align.MetMatchResults
#'
#' @return Writes a png file of the heatmap
#' @export
#'
#'
plotHeatMap = function(mergedData,which="intensity",name){

   require(gplots)


   form = mergedData$formula

   mu = apply(mergedData[,grep(".metab.csv",names(mergedData))],1,mean,na.rm=T)

   mergedData = mergedData[order(-mu),]
   sc = mergedData$score
   mergedData = as.matrix(mergedData[,grep(".metab.csv",names(mergedData))])
   if(which=="intensity"){
      mergedData = log(mergedData)
   }
   mn.mergedData = min(mergedData,na.rm=T)
   mergedData[which(mergedData==0)] = NA

   mergedData[is.infinite(mergedData)] = NA
   colnames(mergedData) = NULL
   rownames(mergedData) = form

   mx.mergedData = max(mergedData,na.rm=T)
   bks <- seq(mn.mergedData,mx.mergedData,length=63)
   cols <- colorRampPalette(c("darkblue","cornflowerblue","orange","red","pink"))(62)
   #dev.new()


   png(name,    # create PNG for the heat map
       width = 5*600,        # 5 x 600 pixels
       height = 8*600,
       res = 600,            # 600 pixels per inch
       pointsize = 8)
   if(which== "intensity"){
      xname = "Log Intensity"
   }else{
      xname = "Score"
   }
   tryCatch(heatmap.2(mergedData,trace="none",  scale="none",xlab="",  na.color="darkblue", dendrogram = "none",
                      sepcolor="white",colsep=1:dim(mergedData)[2],lhei=c(.75,5),lwid=c(1,4),key.title=NA,
                      key.xlab=xname,col=cols, breaks=bks,key=T,Rowv=F,Colv=F, margins=c(10,20),cexRow=0.25),
            error = function(e) print(""))
   dev.off()


}



