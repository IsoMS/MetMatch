#' getMetMatchResults
#'
#' Extracts a data frame of results from the group of metmatch results.
#'
#' @param wd
#' @param metFiles
#' @param cutoff
#'
#' @return Dataframe of values used for alignent by rt and formulae.
#' @export
#'
#'
getMetMatchResults = function(metFiles,cutoff=0){

   ####################################################################
   ## Get the total list of formula and feature vars across all samples
   fh = NULL
   formula = NULL
   rt = NULL
   mz = NULL
   M = NULL
   ion = NULL

   int = NULL
   score = NULL
   for(s in metFiles){
      d = read.csv(s)
      fh = c(fh,   rep(s,length(d$rt.Q)) )
      formula = c(formula,as.character(d$best.hit))
      mz = c(mz,d$exactMass.Q)
      ion = c(ion,as.character(d$ion.Q))
      int = c(int,d$intensity)
      score = c(score,d$score.Q)
      rt = c(rt,d$rt.Q)
   }

   w = order(rt)
   int = int[w]
   rt = rt[w]
   mz = mz[w]
   ion = ion[w]
   score = score[w]
   formula = formula[w]
   fh = fh[w]
   M = convertIonToNeutralMass(ion,mz)

   #ID = unlist(lapply(fh, function(x) unlist(strsplit(x,"/"))[2]))
   #ID = sub(".metab.csv","",ID)
   #ID = sub(".mzXML_cam.csv","",ID)
   ID = fh

   ##################################################################
   df = data.frame(ID ,mz,M,ion,rt,int,formula,score)
   df = df[1:length(df[,1])-1,]
   #####################################################################
   ### Apply a score cutoff based on 5% FDR (or similar)
   s.df = df[which(df$score>=cutoff),]
   s.df
}


#' groupIonsByRT
#'
#' Groups memtmatch results into data chunks organized by formula and separated by a retention time clutering method.
#'
#' @param results
#' @param maxRtDiff
#'
#' @return List of aggregated data, aggregated by formula and RT.
#' @export
#' @keywords internal
#'
groupIonsByRT = function(results,maxRtDiff = 60){

   #################################################################
   #### 1. get the list of u formula and get the data in chuncks by formula
   u.formula = unique(results$formula)
   f.chunks = lapply(u.formula,function(x) results[results$formula == x,])


   ######################################################################
   ## cluster peaks by RT, then create separate entries by class
   chunks = list()
   i=1
   for(tm in f.chunks){

      if(length(tm$rt)>1){
         d = dist(tm$rt,method="euclidean")
         h = hclust(d)
         cut = cutree(h,h=maxRtDiff)
         tm = cbind(tm,"class" = cut)
         nC = length(unique(cut))
         if(nC > 1){
            for(c in unique(cut)){
               s.tm = tm[tm$class==c,]
               s.tm$formula = paste(s.tm$formula,"_",c,sep="")
               chunks[[i]] = s.tm
               i=i+1
            }
         }else{
            chunks[[i]] = tm
            i=i+1
         }

      }else{
         chunks[[i]] = tm
         i=i+1
      }

   }
   chunks
}





#' filterDuplicateIons
#'
#' Post processing helper function.
#'
#' @param g.Ions
#'
#' @return
#' @export
#' @keywords internal
#'
#'
filterDuplicateIons = function(g.Ions){

   removeExtraPeaks = function(tm){
      df = data.frame()
      mtch = match(tm$ID,tm$ID)
      for(m in unique(mtch)){

         get = tm[which(mtch==m),]
         df = rbind(df,get[which(get$int == max(get$int)),])

      }
      df
   }

   chunks = list()
   for(i in 1:length(g.Ions)){
      tm = g.Ions[[i]]
      mtch = match(tm$ID,tm$ID)
      if( length(mtch) > length(unique(mtch)) ){
         chunks[[i]] = removeExtraPeaks(tm)
      }else{
         chunks[[i]] = tm
      }
   }
   chunks
}


#' filterIonsByCutoff
#'
#' @param g.Ions
#' @param N
#' @param cutoff
#'
#' @return
#' @export
#' @keywords internal
#'
filterIonsByCutoff = function(g.Ions,N=1,cutoff=0){
   j=1
   chunks = list()
   for(i in 1:length(g.Ions)){
      tm = g.Ions[[i]]
      if(max(tm$score)>cutoff && length(tm$rt) >= N){
         chunks[[j]] = tm
         j=j+1
      }
   }
   chunks
}





#' groupIonsByFilename
#'
#' Helper function for prostprocessing
#'
#' @param g.ions
#' @param metFiles
#'
#' @return
#' @export
#'
#' @keywords internal
#'
groupIonsByFilename = function(g.ions,metFiles){

   #####################################################################
   ## REpopulate the DATA frame with empty spots:
   ##       For each chunk, create an empty data frame then populate
   ##        it with the missing chunk data
   for(i in 1:length(g.ions)){
      tm = g.ions[[i]]
      slots =  match(metFiles,tm$ID)
      tm = cbind(metFiles,tm[slots,])
      ### Cleanup
      row.names(tm) = NULL
      tm$ID = metFiles
      tm = tm[,-1]
      g.ions[[i]] = tm
   }
   g.ions
}


#' mergeData
#'
#' @param g.ions
#' @param which
#'
#' @return
#' @export
#'
#' @keywords internal
#'
mergeData = function(g.ions,which = "intensity"){

   #########################################################################
   ###   Finally lets build a raw data table for analysis
   out = data.frame()
   form = featID = ion = mn.rt=mn.mz=mn.M= mn.score = cv.score= mn.int = cv.int = n.samples = array()

   for(i in 1:length(g.ions)){

      tm = g.ions[[i]]
      tm$score[which(tm$score==0)] = NA
      form[i] = max(as.character(tm$formula),na.rm=T)
      ion[i] = as.character(unique(tm$ion[which(!is.na(tm$ion))]))
      mn.rt[i] = round(mean(tm$rt,na.rm=T),3)
      mn.mz[i] = round(mean(tm$mz,na.rm=T),4)
      mn.M[i] =  convertIonToNeutralMass(ion[i],mn.mz[i])
      mn = round(mean(tm$int,na.rm=T),1)
      sd = sd(tm$int,na.rm=T)
      mn.int[i] = mn
      cv.int[i] = round(sd/mn * 100,2)
      n.samples[i] = length(tm$int[!is.na(tm$int)])
      mn.score[i] = round(mean(tm$score,na.rm=T),2)
      sd.score = sd(tm$score,na.rm=T)
      cv.score[i] = round(sd.score/mn.score * 100,2)

      featID[i] = paste(mn.mz[i],"@",mn.rt[i],sep="")
      if(which == "intensity"){
         out = rbind(out, t(tm$int))
      }else{
         out = rbind(out, t(tm$score))
      }
   }

   colnames(out) =tm$ID
   meta = data.frame(ID = featID,mz = as.numeric(mn.mz),ion=ion,M = as.numeric(mn.M),rt = as.numeric(mn.rt),int = as.numeric(mn.int),
                     cv.int = as.numeric(cv.int),formula = as.character(form), score = as.numeric(mn.score),
                     cv.score = as.numeric(cv.score),
                     n.samples = as.numeric(n.samples))

   w.table = cbind(meta,out)
   w.table = w.table[order(w.table$rt),]
   w.table

}

#' align.MetMatchResults
#'
#' Postprocesses groups of MetMatch results by formula and retention time.
#'
#' @param wd  The working directory containing the metmatch results
#' @param cutoff  A minimum score.  Must be >=0;
#' @param N Minimum number of samples in which an alinged FFM must be present.
#' @param maxRtDiff maximum retention time deviation (sec) in which the formulae must be present between samples
#'
#' @return
#' @export
#'
align.MetMatchResults = function(wd,cutoff=0.1,N=1,maxRtDiff=60){
   startwd = getwd()
   setwd(wd)

   metFiles = list.files(pattern = ".metab.csv$")
   results = getMetMatchResults(metFiles,cutoff=cutoff)
   g.ions = groupIonsByRT(results)
   g.ions = filterDuplicateIons(g.ions)
   g.ions = filterIonsByCutoff(g.ions,N=N,cutoff=cutoff)
   g.ions = groupIonsByFilename(g.ions,metFiles)
   s.df = mergeData(g.ions,which="score")
   i.df = mergeData(g.ions,which="intensity")
   s.df = s.df[order(-s.df$score),]
   i.df = i.df[order(-i.df$score),]
   write.csv(file = "MetMatch_groupedData_SCORES.csv", s.df,row.names=FALSE,na = "")
   write.csv(file = "MetMatch_groupedData_INTENSITIES.csv", i.df,row.names=FALSE,na = "")
   setwd(startwd)
   list("scores" = s.df, "intensities" = i.df)
}


