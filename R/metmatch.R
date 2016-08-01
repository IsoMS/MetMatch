
#' MetMatch, the molecular formula scoring algorithm with local error rate estimation
#'
#' @param Q.data  The Query data, a list of query objects
#' @param DB  The search database.  When using a decoy database, be sure that run.type option is set to "decoy"
#' @param mass.tol Mass tolerance in part per million (ppm)
#' @param int.rel.tol Relative intensity tolerance (0-1)
#' @param score.method Scoring method, options are "best" (default) or "delta"
#' @param score.type Score type, options are "both", "mass" or "intensity"
#' @param run.type Options are "nodecoy" (default) or "decoy".  Choose decoy when using a decoy database.
#' @export
#' @return List of query-results objects.
metmatch = function(Q.data,DB,mass.tol = 50, rel.int.tol = 0.5, Mono = FALSE, score.method =
                       "best", run.type = "nodecoy",score.type="both",out.file.name) {

   db.data = get.DB.dataframe(DB,6)
   q.data = get.Query.data(Q.data)
   query.isotope.data = get.query.Isotope.data(q.data)


   ##########################################################################
   ######   Compute delta masses and intensities on experimental data vs TFL
   delta.data = compute.Deltas(query.isotope.data,q.data, db.data,mass.tol,rel.int.tol,score.type)

   ### EM the mixture model
   delta.models = model.Deltas(delta.data)

   ### Compute the scores for each feature-formula-match
   scored.data = compute.model.Scores(delta.data,delta.models,score.method, score.type)


   if(run.type == "decoy"){


      ##########################################################################################################
      #### Generate the decoy database.   This is built using model parameters learned from the second best hits.
      print("Computing Decoy database, this may take a while");
      ptm <- proc.time()
      decoy.DB = generate.decoy.spectra(DB,mass.tol,delta.data$second.dI[,1],5);
      decoy.DB = c(DB,decoy.DB);
      proc.time() - ptm
      dec.db.data = get.DB.dataframe(decoy.DB,6)
      ## END decoy generation


      #############################################################################
      ######   Compute delta masses and intensities on experimental data vs TFL+DFL
      dec.delta.data = compute.Deltas(query.isotope.data,q.data, dec.db.data,mass.tol,rel.int.tol,score.type)

      ### EM the mixture model
      dec.delta.models = model.Deltas(dec.delta.data)
      dec.scored.data = compute.model.Scores(dec.delta.data,dec.delta.models,score.method,score.type)
      computed.FDR.data = compute.FDR(dec.scored.data,dec.db.data,score.method)


      scored.data = list(scored.data = scored.data$scored.data,
                         unscored.data = scored.data$unscored.data,
                         dec.scored.data = computed.FDR.data$scored.data,
                         dec.unscored.data = computed.FDR.data$unscored.data)

   }

   write.csv(file = paste(out.file.name, ".metab.csv", sep = ""),row.names=FALSE,
             rbind(scored.data$scored.data,scored.data$unscored.data))

   ###Plot the report
   pdf(paste(out.file.name,".pdf",sep=""))

      par(mfrow=c(2,2))
      par(mar=c(5,5,5,5))

      print("Plotting results")

      plotDeltas(delta.models,mass.tol)
      if(run.type=="decoy"){
          cutoffs = compute.thresholds(scored.data$dec.scored.data)
             write.table(file=paste(out.file.name, ".cutoffs.txt",
                                     sep = ""),cutoffs,sep=" ")
         plotScores(scored.data,run.type="decoy")
      }
      else{
         plotScores(scored.data)
      }

   dev.off();

}
