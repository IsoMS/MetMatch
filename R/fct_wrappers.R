

#' Annotate feature isotopes.
#'
#' Runs the CAMERA pipeline.
#'
#' @param xs The XCMS xset
#' @param ppm CAMERA mz tolerance window in ppm
#' @param mzabs CAMERA mz deviation in mz
#' @param sigma sigma for CAMERA groupFWHM function
#'
#' @return
#' @export
#'
#' @examples
getIsotopes = function(xs,fh,ppm=5,mzabs=0.02,sigma=3){
   require(CAMERA)
   an <- xsAnnotate(xs)
   an = groupFWHM(an,sigma=sigma)
   an = findIsotopes(an,maxiso=6, ppm=ppm, mzabs=mzabs)
   pl = getPeaklist(an)
   write.csv(file=fh,pl)
}


#' run.CAMERA
#'
#' Wrapper function to run the XCMS/CAMERA feature finding software.This is a wrapper for a single or group of mzXML files in a directory.
#'
#' @param wd   Working directory containing the raw mzxml files
#' @param save Saves the XCMS set as a RDA file.
#' @param ...
#'
#' @return
#' @export
#'
#'
run.CAMERA = function(wd,method="massifquant",prefilter=c(3,500),peakwidth=c(3,45),snthresh=1.5,withWave = 1,save=FALSE,...){

   require(CAMERA)
   startwd = getwd()
   setwd(wd)
   mzxmlFiles = list.files(pattern = ".mzXML$")
   camNames = sub(".mzXML",".cam.csv",mzxmlFiles)

   xs = list();

   i = 1;
   for(mzfh in mzxmlFiles){
      xs[[i]] = xcmsSet(mzxmlFiles[i],method=method,ppm=10,
              prefilter=prefilter,peakwidth=peakwidth,snthresh=snthresh,withWave=withWave)
      getIsotopes(xs[[i]],fh=camNames[i])
      i=i+1
   }
   if(save==TRUE){
      save(file = "xcmsResults.rda",xs)
   }

   setwd(startwd)

}


#' run.metmatch
#'
#' Wrapper function for searching all isotope files in a directory.
#'
#' @param wd
#' @param DB
#' @param save
#' @param nCore
#' @param score.method
#' @param run.type
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
run.metmatch = function(wd,DB,save=TRUE,score.method="best",run.type="nodecoy",score.type = "both"){

   #mclapply = getSys()
   startwd = getwd()
   setwd(wd)

   camFiles = list.files(pattern="\\.cam.csv$")
   sampleNames = sub(".cam.csv","",camFiles)

   print("Converting CAMERA files to searchable objects and running searches.")

   res = list()
   i = 1
   for (cf in camFiles) {
      Q = convert.camera.file(cf,1)
      res[[i]] = metmatch(Q,DB,score.method=score.method,run.type=run.type,score.type=score.type,out.file.name=sampleNames[i])
      i=i+1
   }

   if(save==TRUE){
      save(file="Q.rda",Q)
      save(file="metmatch.results.RDA",res)
   }

   print("Runs complete!")
   setwd(startwd)

}





