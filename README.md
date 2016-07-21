
<img src="/images/logo.png" width="600px"/>

##MS1 Global Search Engine for Feature to Formula Matching 

MetMatch is a MS1 search engine for global mass spectrometry feature to formula matching, developed by Scott Walmsley,PhD (University of Colorado-Anschutz Medical School) and Hyungwon Choi,PhD (National University of Singapore School of Public Heatlh).  MetMatch takes as input a user supplied list of formula and annotated features produced from the XCMS/CAMERA feature finding  and isotope annotation software.  Included are methods to build a target feature library (TFL), preprocessing routines, running the MetMatch algorithm, and postprocessing routines. In this vignette you will find a generic workflow for analyzing the data produced from a replicate mass spectrometry experiment.  In replicate mass spectrometry experiments, each sample is independently and rapidly searched for correct formula matching (no prior feature alignments necessary like current tools) and are aligned after the MS1 search.   The FDR is estimated per sample using a ***novel application of a decoy database*** aiding the user in selecting a minimum cutoff score. Additionally, this protocol can be run for MS1 data acquired in MS1 only mode, or for simultaneous acquired MS1 and MS2 data (MetMatch only searched the MS1 data) as demonstrated in our co-submitted MetaboDIA software framework (Genbo et.al, manuscript pending). In this approach, the coacquired MS2 structural data were demonstrated to have improved accuracy of MS2 identification due to MetMatch's MS1 accuracy in formula matching and in dealing with retention time differentiable isobars.

##Installation 
#####(Source currently not available--7/21/2016--check back soon)
After downloading the compressed archive (tar.gz), install with R as you would any other package:

```
install.packages(path_to_file, repos = NULL, type="source")
```

##Usage
###1. Quickstart
For those wanting to dive in, the following commands should allow you to quickly analyze a single mzxml file.  MetMatch was developed with a replicate series of MS samples, and users are encouraged to read further for analyzing replicate sets of MS metabolomics data.
```
#1. Read the list of formulae to build the TFL
formulae = read.delim("formulae.txt")

#2. Read the TFL (database)
data(TFL)

#3. run XCMS/CAMERA:
wd = "<Path to your folder containing the mzxml file(s)>"
run.CAMERA(wd)

#4. Run MetMatch!
run.metmatch(wd,TFL)

#5. Open your results files!
```

###2. Target Feature Library
Target Feature Libraries (TFL) are the samplewise specific databases (herein called a search-space) and are one of the two required inputs for MetMatch.  Building a TFL can take up to several minutes depending on the number of formulae, but is completed once for a list of unique formulae expected to be encountered in a sample.  For our example, formulae are derived from the Human metabolome database and LipidMaps, since our samples are derived from a human cell line.  Additionally, once the TFL is built, the user has the option of saving the TFL to a text format file for later rapid loading.  Routines can be scripted to remove or add TFL entries simply by running the build command, writing to text the new TFL entries, and inserting these entries into the searchspace text file.

#####Building the TFL
A list of formulae is required to initiate the TFL building proces. After reading formulae, a simple command is provided to clean the formulae of unwanted entries (deuterated, containing R groups, single atom entries, etc...). We then initiate the TFL build. k is the expected maximum number of isotopes (we use 6).  Before running this, please realize that a precompiled searchspace has been built to save the user time (see below):

```
formulae = read.delim("formulae.txt")
formulae = clean.formulae(formulae)
TFL = build.DB.emass(k=6,formulae,nCore=64)
```
#####Writing TFLs and reading pre-built TFLs}
Two simple commands are provided to write and read the TFL:

```
write.DB(file="TFL.txt")
read.DB(file="TFL.txt")
```
Or you can use our precompiled library (HMDB+Lipidmaps):
```
data(TFL)
```

The format of each entry in the TFL text file is <\# formula \#peaks> on the first line followed by the isotopic information: <Formula monoisotopic mass, \#peaks, isotope exact mass, isotope relative abundance>.
```
#  C7H11N3O2 6
C7H11N3O2	169.08457	6	169.08457	0.909003
C7H11N3O2	169.08457	6	170.08722	0.083435
C7H11N3O2	169.08457	6	171.08928	0.007122
C7H11N3O2	169.08457	6	172.09161	0.000421
C7H11N3O2	169.08457	6	173.09385	1.9e-05
C7H11N3O2	169.08457	6	174.09607	1e-06
#  C3H10N2 5
C3H10N2	74.08384	5	74.08384	0.958645
C3H10N2	74.08384	5	75.08621	0.04068
C3H10N2	74.08384	5	76.08824	0.00067
C3H10N2	74.08384	5	77.08985	5e-06
C3H10N2	74.08384	5	78.0911	0
```

###3. MS Data Preprocessing
#####Peak finding with XCMS/CAMERA
At this time in the data analysis, the user will preprocess data to find isotopic clusters (features).  It is recommended to use the ProteoWizard MSConvert.exe tool to convert your data to the mzXML format. The pipeline is run using mzXML formatted files, but it is assumed that any format compatible with XCMS can be used.  Once peak finding is completed, the resultant peaks are grouped into isotope clusters using CAMERA. The user is expected to have some working knwoledge of the settings in CAMERA, but default setting are provided which should be applicable to most QTOF instruments. The peaklists are written to a csv file and can be analyzed later with the MetMatch algorithm.  These steps by far are the most time consuming steps in the MetMatch pipeline.  We provide a wrapper function to perform this task.  Adapting MetMatch to other feature finding and isotope grouping software (eg vendor software) is in the works and will be released in sequential versions of the software package.  

```
run.CAMERA(wd)
```

### 4. Performing Feature to Formula Matching (FFM) using the MetMatch search engine
After data preprocessing, the features have been written to comma delimited format (.csv) and are ready to read and convert to a format that MetMatch can search. This is a simple, two step process: 1) read and convert the CAMERA /annotated isotope  output to a list of query features containing isotopic masses and intensities, and 2) run the MetMatch algorithm.  We provide a wrapper function to run these tasks.
MetMatch requires the user to set the ppm and intensity tolerances.  The most important setting here is the mass setting.  To model delta mass and delta intensities correctly, use a wide enough mass tolerance to help model the *incorrect* matches. For example, in QTOF data, metabolites with a mass error of less than 10ppm are routinely acquired, however we use a default *50ppm* mass tolerance window.  Intensity tolerances are on a *relative intensity* scale of 0-1, and the default setting of 0.5 will work for most data.  **MetMatch** takes as input the query list of experimental isotopes (Q), and the target feature library created earlier. Please see the manual for details on specific settings in the MetMatch function.

```
wd = <relative or full path to directory of experimental files>
run.metmatch(wd,TFL)
```
### 5. Estimating the FFM error rate
In the MetMatch algorithm, the FFM error rate is computed by running the *"decoy"* option for either the run.metmatch  wrapper function, or within the MetMatch function itself.  This biggest disadvantage for this step is that the decoy option takes longer to run (~3 times longer than without the "decoy" option).   The advantage, however, is similar to a decoy peptide database in it's usage with proteomics MSMS search engines, MetMatch also utilizes a similar approach to aid determining the cutoff "FFM score" that will be computed for each sample.   This method is adaptive to the quality of the sample's features, and consistent in producing cutoffs (ideal minimum scores) between similar samples. 
#####Decoy Feature Library Contruction
The DFL is constructed automatically by learning the correct parameters for DFL construction from the samples themselves.  Specifically, models computed from the isotope delta masses and intensities are used to parameterize the DFL contruction.   More details can be found in the manuscript online methods (pending). We chose to construct 5 DFL entries per TFL entry to increase the precision in computing the threshold, especially since the TFL database is only ~13000 entries small, and this helps to accurize the detected score cutoff. 
#####Determining the cutoff score
The minimum score required or "trusted" by the experimentor is computed within the MetMatch algorithm.   Just simply use the "decoy" option and a small text file with the extension ".cutoff.txt" will be produced for each file searched in your working directory.   The results will include scores for FDR of 1 and 5%, and are plotted in the diagnostic plots showing the distribution of scores (see below).  Typically, a 5% FDR has been indicative of a good quality FFM.
To compute your FDR, simply run the decoy option:
For all files in a working directory:
```
wd = <relative or full path to files>
run.metmatch(wd,TFL,run.type="decoy")
```

Or optionally, for advanced users who want fine grain control over the metmatch algorithm:
```
?metmatch  # look at the help file for options.
metmatch = function(Q.data,TFL,run.type = "decoy",out.file.name)
```
-->Cutoff scores produced will be in the cutoffs.txt file.

### 6. Interpreting search results
MetMatch produces a table of results together with a plot of results for visualizing the delta models and the final score.  The tables are fairly comprehensive. However the most important columns to the typical experimentor will be the "query" ion, mz, and retention time, assigned formula, and the relevant score.  These tables are indicated with the extenstion ".metab.csv".   At this point, the user maps the best formula to all possible compound names in their database, and is at an advantage if they have coacquired MS2 data.
#####Tables
Below is 3 lines from the top hits of a search result contained in a "metab.csv" file.   Q stands for "query" or the experimental ion matched against the TFL.   The type of ion detected, it's mass and rt and total intensity is copied to the search output.  The FFMs are sometimes several per ion, with the highest scoring FFM presented as the *best.hit .  The FFMs DB indexes and formulae are shown, and are followed by the scores.  The independent mass and intensity scores as well as the final scores are shown. 

index.Q | ion.Q | exactMass.Q | rt.Q | intensity | index.DB | hit.DB | ... | score.mass.Q	| score.int.Q	| score.Q	| best.hit
----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | -----
1557|M+H|720.590851|166.809|1308855.385|10756;10844|C40H82NO7P;C41H86NO6P|...|12.26601007|8.441255533|20.7072656|C40H82NO7P
1622|M+H|706.5395207|163.312|3614270.539|1949;10765|C38H76NO8P;C39H80NO7P|...|9.566656057|9.167408573|18.73406463|C38H76NO8P
1844|M+H|792.5528487|89.386|15171.64471|1999;2057;2948;11015|C45H78NO8P;C46H82NO7P;C42H82NO10P;C44H74NO9P|...|11.93744504|6.768302973|18.70574801|C45H78NO8P

#####Plots
Four diagnostic plots are included. The first two (A and B) are plots of the delta isotope masses and delta isotope intensities from the entire sample.   These have a model distribution fit to them and indicate nearest the peak of the model (blue line) which FFMs are likely "good" FFMs versus the rest of which are bad. The log ratio score for each FFM is computed from this distribution over the mean of the remaining "bad" delta masses and intensities.   The third plot (C) shows a heatmap distribution of these scores, which when added together become the feature's final score. The final distribution of scores (D) for a sample are shown.  When the "decoy" option is run, the decoy FFMs (indicating a false hit) are shown as "+"'s in plot C and are used to compute the FDR and derive a cutoff score for 1 and 5% FDR (red dash lines in plot D).   
 
 
Plot A. Delta Masses  | Plot B. Delta Intensities
--------------- | ----------------
<img src="/images/A.png" width="300px"/>|<img src="/images/B.png" width="300px"/>


Plot C.  Mass Scores vs. Intensity Scores|Plot D. Distribution of FFM scores with smooth line fitted (blue)
--------------- | ----------------
<img src="/images/C.png" width="300px"/>|<img src="/images/D.png" width="225px"/>


### 7. Aligning results from multiple samples  
In our study, we searched multiple samples which were aligned together based on their FFM and retention time (using a hierarchial clustering method). Different extraction types are clustered separately. This was because all were derived from the same cell line (HEK293), but different chromatographies were used (C18 versus HILIC).  Hence the user will want to keep alignment of the data in groups determined by the extraction type and type of chromatographic separation completed in the experiment (eg. do NOT align C18 data with HILIC data).  We provide several postprocessing functions and a wrapper function to align your data.   The user will want some knowledge of the reproducibility of their chromatographic separations and subsequent MS data runs, and this is the necessary work required of any user in any good mass spectrometry experiment.  This will enable correct selection of the RT window for clustering of FFMs from the multiple samples.  The result of this workflow is a table of  aligned FFMs, and there is also a ploting function to visualiZe the distribution of scores and intensities in a heatmap.
Alinging samples is easy, with just one set of commands:
```
wd = <relative or full path to .matab.csv files>
results = align.MetMatchResults(wd,maxRtDiff=60,cutoff=1,N=1) #cutoff: minimum MetMatch score (see if lower scoring FFMs align with higher scoring ones. N: minimum number of replicates required for the observed FFM for clustering.
```
These command produce the tabel of merged scores and intensities between samples: 
##### Table of merged data:

ID|M|rt|int|cv.int|formula|score|cv.score|n.samples|r001.metab.csv|r002.metab.csv|r003.metab.csv|...
----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | -----
720.5913@166.947|720.5913|166.947|1303148.9|2.14|C40H82NO7P_1|18.75|33.73|10|1305446.52|1322333.483|1328024.739|...
690.5071@88.719|690.5071|88.719|55103.6|3.26|C37H72NO8P_1|17.94|21.59|2|53832.8515| | |...		
706.5397@163.151|706.5397|163.151|3626055.3|1.34|C38H76NO8P_1|16.91|50.01|10|3693877.078|3649969.047|3665835.568|...
577.5198@48.066|577.5198|48.066|191926.6|3.29|C37H68O4_1|16.46|33.85|10|179228.4324|189135.7611|200455.9577|...
768.5536@97.104|768.5536|97.104|31742.2|4.62|C43H78NO8P_1|16.34|81.92|9|30703.91272|33240.90292|32650.57507|...
863.5645@49.022|863.5645|49.022|41230.9|4.73|C45H83O13P|16.08|59.99|9|41450.61353|41800.97931|45552.54431|...

#### HeatMap of merged data:
Finally, you can visualize the relative distribution of intensities and scores as aligned using the previous functions.  The rest is up to the user/scientist to determine how best to analyze their data further, but MetMatch to this point provides accurate Feature to Formulae matching.  

To plot the features using the FFM information:
```
plotHeatMap(results$intensities,which="intensity", name="aligned_heatmap_INTENSITIES.png")
plotHeatMap(results$scores,which="score",name="aligned_heatmap_SCORES.png")
```

Aligned Scores|Aligned Intensities
--------------- | ----------------
<img src="/images/score_mat.png" width="400px"/>|<img src="/images/int_mat.png" width="400px"/>


#CONCLUSION
MetMatch provides feature to formula matching in a global analysis framework.   The algorithm provides rapid formulae assignments per sample.  In replicate mass spectrometry experiments, each sample is independently and rapidly searched for correct formula matching (no prior feature alignments necessary like current tools) and are aligned after the MS1 search.   The FDR is estimated per sample using a novel application of a decoy database aiding the user in selecting a minimum cutoff score.












