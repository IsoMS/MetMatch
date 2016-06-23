# MetMatch
MS1 Global Search Engine for Feature to Formula Matching 

MetMatch is a MS1 search engine for global mass spectrometry feature to formula matching. MetMatch takes as input a user supplied list of formula and features produced from the XCMS/CAMERA feature finding  and isotope annotation software.  Included are methods to build a target feature library (TFL), preprocessing routines, running the metmatch algorithm, and postprocessing routines. In this vignette you will find a generic workflow for analyzing the data produced from a replicate mass spectrometry experiment.


###1. Quickstart
For those wanting to dive in, the following commands should allow you to quickly analyze a single mzxml file.  MetMatch was developed with a replicate series of MS samples, and users are encouraged to read further for analyzing replicate sets of MS metabolomics data.
```
#1. Read the list of formulae to build the TFL
formulae = read.delim("formulae.txt")


#2. Build the TFL (database)
TFL = build.TFL.emass(k=6,formulae)

#3. run XCMS, enter your own mzXML file
xs = xcmsSet("data.mzXML",method="massifquant",withWave=1,ppm=10,
              prefilter=c(3,300),peakwidth=c(3,45),snthresh=1.25)



#4. Run CAMERA to annotate isotopes

#Wrapper function
getIsotopes = function(x){
  require(CAMERA)
  an = xsAnnotate(x)
  an = groupFWHM(an,sigma=3)
  an = findIsotopes(an,maxiso=6, ppm=5, mzabs=0.02)
  an = groupCorr(an,calcIso=TRUE)
  getPeaklist(an)
 }

pl = getPeakList(xs)\n
write.csv(file="isotopes.cam.csv",pl)

#5. Convert camera file to query list
Q = convert.camera.file("isotopes.cam.csv",1)

#6. Run MetMatch!
res = metmatch(Q.data,TFL,score.method="best",run.type="nodecoy",out.file.name="results")
```

###2. Target Feature Library
Target Feature Libraries (TFL) are the samplewise specific databases (herein called a search-space) and are one of the two required inputs for MetMatch.  Building a TFL can take several minutes to hours depending on the number of formulae, but is completed once for a list of unique formulae expected to be encountered in a sample.  For our example, formulae are derived from the Human metabolome database and LipidMaps, since our samples are derioved from a human cell line.  Additionally, once the TFL is built, the user has the option of saving the TFL to text for later rapid loading.  Routines can be scripted to remove or add TFL entries simply by running the build command, writing to text the new TFL entries, and inserting these entries into the searchspace text file.

#####Building the TFL
A list of formulae is required to initiate the TFL building proces. We then initiate the TFL build. k is the expected maximum number of isotopes (we use 6).  Before running this, please realize that a precompiled searchspace has been built to save the user time (see below):

```
formulae = read.delim("formulae.txt")
TFL = build.TFL.emass(k=6,formulae,nCore=64)
```



