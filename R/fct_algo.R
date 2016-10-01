#'  Laplace Distribution
#'
#'  Given a set of x values, mean and scale parameter, returns a Laplace distribution.
#'
#'
#' @param x  Single value or array of values
#' @param m  mu
#' @param b  scale
#'
#' @return y value of Laplace distribution
#' @export
#' @keywords internal
dlaplace = function(x, m, b) {
  1 / (2 * b) * exp(-abs(x - m) / b)
}

#' log Laplace Distribution
#'
#' Given a set of x values, mean and scale parameter, returns a log Laplace distribution.
#'
#' @usage logdlaplace(x,m,b)
#'
#' @param x  Single value or array of values
#' @param m  mu
#' @param b  scale
#'
#' @return y value of Laplace distribution
#' @export
#' @keywords internal
logdlaplace = function(x, m, b) {
  -log(2 * b) - abs(x - m) / b
}


#' Mixture model for Laplace Uniform distribution
#'
#' Performs expectation-maximization to fit the Laplace-Uniform distributions given input set of X values.
#'
#' @param x  Input array of values of X
#' @param tol Tolerance for minimizaton
#' @param maxiter Maximum number of iterations
#'
#' @return List of parameters
#' @export
#' @keywords internal
laplace.unif.EM = function(x, tol = 1e-5, maxiter = 100) {
  n = length(x)
  tmp = density(x)
  m = tmp$x[tmp$y == max(tmp$y)]

  b = 1  ## initial value
  z = rep(0.15, n)
  pi.match = 0.5

  par.diff = 1e5
  iter = 0

  ### mode m is fixed, just update the variance of laplace 2b^2
  while ((iter < maxiter) & (abs(par.diff) > tol)) {
    new.v = sum(z * (x - m) ^ 2) / sum(z)
    new.b = sqrt(new.v / 2)

    p1 = pi.match * dlaplace(x, m, sqrt(new.b))
    p2 = (1 - pi.match) * dunif(x, min(x) - 0.01, max(x) + 0.01)
    z = p1 / (p1 + p2)
    new.pi.match = mean(z)

    par.diff = max(abs(c(new.pi.match - pi.match, new.b - b)))
    iter = iter + 1
    if (iter %% 10 == 0)
      cat(iter)

    b = new.b
    pi.match = new.pi.match

  }
  cat("\n")

  res = list(
    mu = m, b = b, xmin = min(x) - 0.01, xmax = max(x) + 0.01, pos = pi.match,unif = mean(p2,na.rm =
                                                                                            TRUE)
  )
  res
}


#' Compute Mass Tolerance
#'
#' Computes the mass tolerance window for given mass and desired ppm
#'
#' @param M   Input neutral mass
#' @param ppm   Desired ppm tolerance
#'
#' @return Mass tolerance
#' @export
#' @keywords internal
massTol = function(M,ppm) {
  as.numeric((ppm * M) / 1e6)
}

#' Compute Isotope Delta Masses
#'
#' Computes the delta masses for a single query isotope masses compared to all masses in the database
#'
#' @param mass.Q   query mass
#' @param mass.DB   Database of masses
#'
#' @return Array of delta masses
#' @export
#' @keywords internal
compute.delta.mass = function(mass.Q,mass.DB) {
  deltaM = NULL
  if (is.null(dim(mass.DB))) {
    deltaM = data.frame((mass.Q - mass.DB) / mass.DB * 1e6)
  }else{
    deltaM = data.frame((mass.Q - mass.DB[1,]) / mass.DB[1,] * 1e6)
    for (i in 2:length(mass.DB[,1])) {
      deltaM = rbind(deltaM, (mass.Q - mass.DB[i,]) / mass.DB[i,] * 1e6)
    }
  }
  deltaM
}


#' Compute Isotope Delta intensities
#'
#' Computes the delta intensities for a single query isotope intensities compared to all intensities in the database
#'
#' @param int.Q  The query intensity
#' @param int.DB  The database of intensities
#'
#' @return Array of delta intensities
#' @export
#' @keywords internal
compute.delta.int = function(int.Q,int.DB) {
  deltaI = NULL
  if (is.null(dim(int.DB))) {
    deltaI = data.frame(int.Q - int.DB)
  }else{
    deltaI = data.frame(int.Q - int.DB[1,])
    for (i in 2:length(int.DB[,1])) {
      deltaI = rbind(deltaI, (int.Q - int.DB[i,]))
    }
  }
  deltaI
}



#' Format dataframe of masses or Intensities
#'
#' @param X  Input of query masses or intensites
#' @param n  Number of isotope peaks allowed
#' @param type "mass" or "I" intensity
#'
#' @return Array of mass or intensity vectors
#' @export
#' @keywords internal
get.dataframe = function(X,n,type="mass"){

   X = sapply(X,function(x) strsplit(as.character(x),";"))
   X = sapply(X,function(x) as.numeric(x))

   if(type=="mass") out = data.frame(M0 =sapply(X,function(x) unlist(x[1]) ))
   else out = data.frame(I0 =sapply(X,function(x) unlist(x[1]) ))

   for(i in 2:n){
      numbers = sapply(X,function(x) unlist(x[i]))
      out = cbind(out,numbers)
      if(type=="mass") colnames(out)[i] = paste("M",i-1,sep="")
      else colnames(out)[i] =  paste("I",i-1,sep="")
   }
   out
}



#' get.Query.data
#'
#' @param Q.data
#'
#' @return
#' @export
#' @keywords internal
get.Query.data = function(Q.data){

   QueryIsoMass = sapply(Q.data,function(x)
      paste(x$IsoMass,collapse = ";"))
   QueryIsoInt = sapply(Q.data,function(x)
      paste(x$IsoInt,collapse = ";"))
   QueryMZ = sapply(Q.data,function(x)
      x$mz)
   QueryRT = sapply(Q.data,function(x)
      x$RT)
   QueryIon = sapply(Q.data,function(x)
      x$Ion)
   QueryIntensity = sapply(Q.data,function(x)
      x$sumInt)
   inputIndex = seq(1, length(QueryMZ), by = 1)

   data.frame(inputIndex,QueryIsoMass,QueryIsoInt,QueryMZ,QueryRT,QueryIon,QueryIntensity)

}


#' get.query.Isotope.data
#'
#' Organizes lists of experimental isotope data. Internal function.
#'
#' @param q.data
#'
#' @return list of experimental isotopic peak  data.
#' @export
#' @keywords internal
get.query.Isotope.data = function(q.data){
   M = get.dataframe(q.data$QueryIsoMass,n = 6)
   I = get.dataframe(q.data$QueryIsoInt,n = 6,type = "I")
   sumInt = q.data$QueryIntensity
   mz = q.data$QueryMZ
   rt = q.data$QueryRT
   ion = q.data$QueryIon
   list(inputIndex = q.data$inputIndex,M=M,I=I,sumInt=sumInt,mz=mz,rt=rt,ion=ion)
}



#' compute.Deltas
#'
#' Computes the delta scores for each isotopic peak.
#'
#'
#' @param query.isotope.data
#' @param q.data
#' @param db.data
#' @param mass.tol
#' @param rel.int.tol
#'
#' @return  lists of data frames containing the computed delta values, and unscored peaks.
#' @export
#'
#'
compute.Deltas = function(query.isotope.data,q.data,db.data,mass.tol,rel.int.tol,score.type,Mono=FALSE){
   mass.DB = as.matrix(db.data[, grep("^M", colnames(db.data))])
   int.DB = as.matrix(db.data[, grep("^I", colnames(db.data))])

   d = query.isotope.data
   inputIndex = d$inputIndex
   delta.mass = delta.int = index.DB = index.Q = exactMass.Q = first.dM = second.dM =first.dI = second.dI = sum.dI  =sum.dM = NULL

   for (i in 1:length(d$M[,1])) {
      dM = dI = dM0 = dM1 = w = NULL;
      mass.Q = d$M[i,]
      int.Q = d$I[i,]

      ### filter database for near masses
      max.diff.mass = massTol(mass.Q[1],mass.tol)
      if(score.type=="both"){
         w = which(
            db.data$exactMass > as.numeric(mass.Q$M0 - max.diff.mass)
            &
               db.data$exactMass < as.numeric(mass.Q$M0 + max.diff.mass)
            &
               db.data$Io        < as.numeric(int.Q$I0  + rel.int.tol)
            &
               db.data$Io        > as.numeric(int.Q$I0  - rel.int.tol)
         )
      }
      if(score.type == "mass"){
         w = which(
            db.data$exactMass > as.numeric(mass.Q$M0 - max.diff.mass)
            &
               db.data$exactMass < as.numeric(mass.Q$M0 + max.diff.mass)
         )
      }
      if(score.type == "intensity"){
         w = which(
            db.data$Io        < as.numeric(int.Q$I0  + rel.int.tol)
            &
               db.data$Io        > as.numeric(int.Q$I0  - rel.int.tol)
         )
      }

      sub.mass.DB = mass.DB[w,] #>
      sub.int.DB =  int.DB[w,] #>

      if (!is.na(sub.mass.DB[1])){
         index.DB   = c(index.DB,w)
         index.Q  = c(index.Q, rep(i,length(w)))

         dM = compute.delta.mass(mass.Q,sub.mass.DB);
         dI = compute.delta.int(int.Q,sub.int.DB);

         delta.mass = rbind(delta.mass,dM)
         delta.int  = rbind(delta.int,dI)
         if(Mono == FALSE){
            s.M = apply(dM,1,sum,na.rm=T)
         }
         if(Mono == TRUE){
            s.M = dM[,1]   #### TEST for monoisotopic mass
         }
         s.M = abs(s.M)
         s.M = s.M[order(s.M)]

         s.I = apply(dI,1,sum,na.rm=T)
         s.I = abs(s.I)
         s.I[order(s.I)]

         tp =  dI[order(abs(dI$I0)),]
         first.dI = rbind(first.dI,tp[1,])

         if(length(w) > 1){
            sum.dI = rbind(sum.dI,s.I[2])
            sum.dM = rbind(sum.dM,s.M[2])

            second.dI = rbind(second.dI,tp[2,])
         }
      }
   }

   delta.mass = as.matrix(delta.mass)
   delta.int = as.matrix(delta.int)

   w = which(abs(delta.mass) > mass.tol)

   delta.mass[w] = NA
   delta.int[w] = NA

   temp.df = data.frame(
      index.Q = index.Q,
      ion.Q = d$ion[index.Q],
      exactMass.Q = d$mz[index.Q],
      rt.Q=d$rt[index.Q],
      intensity = d$sumInt[index.Q],
      index.DB = index.DB,
      hit.DB = db.data$formula[index.DB],
      delta.mass,
      delta.int
   )
   list(temp.df=temp.df,query.isotope.data=d,
        delta.mass = as.data.frame(delta.mass),
        sum.dM = sum.dM[,1],
        delta.int = as.data.frame(delta.int),
        sum.dI = sum.dI[,1],
        first.dM = first.dM,
        second.dM = second.dM,
        first.dI=first.dI,
        second.dI = second.dI)

}


#' model.Deltas.
#'
#' Returns the Laplace-uniform model parameters after EM.
#'
#' @param delta.data
#'
#' @return List of model parameters
#' @export
#'
#'
model.Deltas = function(delta.data){
   v.delta.mass = as.vector(delta.data$delta.mass)
   v.delta.mass = v.delta.mass[!is.na(v.delta.mass)]

   v.delta.int = as.vector(delta.data$delta.int)
   v.delta.int = v.delta.int[!is.na(v.delta.int)]

   mod.mass = laplace.unif.EM(v.delta.mass) ## mixture model: mass
   mod.int = laplace.unif.EM(v.delta.int)   ## mixture model:intensity
   list(mod.mass=mod.mass,
       mod.int=mod.int,
       delta.mass=v.delta.mass,
       delta.int = v.delta.int)
}


#' compute.model.Scores
#'
#' Performs the log liklihood scoring
#'
#' @param delta.data
#' @param delta.models
#' @param score.method
#' @param score.type
#'
#' @return lists of scored and unscored data.
#' @export
#'
#'
compute.model.Scores = function(delta.data,delta.models,score.method, score.type){

   ##############################################################
   ####  compute the model scores and then add to the data
   temp.df = delta.data$temp.df
   query.isotope.data = delta.data$query.isotope.data
   mod.mass=delta.models$mod.mass
   mod.int=delta.models$mod.int

   u.index.Q = unique(temp.df$index.Q)

   out.df = NULL  ## the output data frame

   for (q in u.index.Q) {

      score.mass.Q = score.int.Q = score.Q = df.E = E = NULL;

      #### get the entry's multiple values
      E = temp.df[which(temp.df$index.Q == q),]


      #############################################################################################
      ####  Mass and intensity distrbutions are fitted with a Laplace over Uniform random curve.

      ## compute mass(SM scores)
      score.mass = as.matrix(logdlaplace(E[,grep("^M",colnames(E))], mod.mass$mu, mod.mass$b) + log(mod.mass$xmax - mod.mass$xmin))

      ## compute int (SI) scores
      score.int  = as.matrix(logdlaplace(E[,grep("^I",colnames(E))], mod.int$mu, mod.int$b) + log(mod.int$xmax - mod.int$xmin))

      w = which(score.mass < 0)
      w2 = which(score.int < 0)

      if(score.type == "both"){
         score.mass[w] = 0
         score.int[w2] = 0
         score.mass[w2] = 0
         score.int[w] = 0
      }
      if(score.type == "mass"){
         score.mass[w] = 0
         score.int[w] = 0
      }
      if(score.type== "intensity"){
         score.int[w2] = 0
         score.mass[w2] = 0
      }

      score.mass.Q = apply(score.mass,1,sum,na.rm = T)
      score.int.Q  = apply(score.int, 1,sum,na.rm = T)

      if(score.type == "both"){
         score.Q = score.mass.Q  +  score.int.Q
      }
      if(score.type == "mass"){
         score.Q = score.mass.Q

      }
      if(score.type == "intensity"){
         score.Q = score.int.Q
      }

      df.E = cbind(E, score.mass.Q,score.int.Q,score.Q)
      df.E = df.E[order(-df.E$score.Q),]

      rank.hit = seq(1,length(df.E[,1]),by = 1)
      df.E = cbind(df.E,rank.hit)
      out.df = rbind(out.df,df.E)
   }


   out.df = out.df[order(-out.df$score.Q),]
   out.df$hit.DB = as.character(out.df$hit.DB)

   ###################################################
   ###   Collapse multi hits to one line
   tmp.out.df = data.frame()

   u.index.Q = unique(out.df$index.Q)
   delta.score = NULL

   for (e in u.index.Q) {
      tmp.df=NULL
      tmp.df = out.df[which(out.df$index.Q == e),]

      if (length(tmp.df$index.Q) > 1) {
         tmp.df = tmp.df[order(tmp.df$rank.hit),]
         best.hit = tmp.df$hit.DB[1]
         tmp.df[1,]$hit.DB = paste(as.character(tmp.df$hit.DB),collapse = ";")
         tmp.df[1,]$index.DB = paste(as.character(tmp.df$index.DB),collapse = ";")

         delta.score = c(delta.score,tmp.df$score.Q[1] - tmp.df$score.Q[2])
         tmp.df = cbind(tmp.df,best.hit)
         tmp.out.df = rbind(tmp.out.df,tmp.df[1,])

      }else{
         delta.score = c(delta.score,tmp.df$score.Q)
         best.hit = tmp.df$hit.DB[1]
         tmp.df = cbind(tmp.df,best.hit)
         tmp.out.df = rbind(tmp.out.df,tmp.df)
      }
   }

   out.df = tmp.out.df[,-which(colnames(tmp.df) == "rank.hit")]
   out.df = cbind(out.df,delta.score)

   if (score.method == "best") {
      out.df = out.df[order(-as.numeric(out.df$score.Q)),]
   }else{
      out.df = out.df[order(-as.numeric(out.df$delta.score)),]
   }

   w = which(is.na(match(query.isotope.data$inputIndex, out.df$index.Q)))
   unscored.data = data.frame(
      index.Q = query.isotope.data$inputIndex[w],
      ion.Q = query.isotope.data$ion[w],
      exactMass.Q = query.isotope.data$mz[w],
      rt.Q = query.isotope.data$rt[w],
      intensity = query.isotope.data$sumInt[w],
      index.DB = NA,
      hit.DB = NA,
      M0 = NA, M1 = NA, M2 = NA, M3 = NA, M4 = NA, M5 = NA,
      I0 = NA, I1 = NA, I2 = NA, I3 = NA, I4 = NA, I5 = NA,
      score.mass.Q = NA,
      score.int.Q = NA,
      score.Q = NA,
      best.hit = NA,
      delta.score = NA)

   list(scored.data = out.df,unscored.data = unscored.data)
}


#' compute.FDR
#'
#' computed the FDR when using the decoy option in  metmatch.
#'
#' @param dec.scored.data
#' @param dec.db.data
#'
#' @return Lists of scored and unscored data with FDR
#' @export
#'
#'
compute.FDR = function(dec.scored.data,dec.db.data,score.method){
   scored.data = dec.scored.data$scored.data

   PSIR = FDR = hit.type = index.score = NA

   hit.type[grep("decoy",scored.data$best.hit)] = "decoy"
   hit.type[grep("decoy",scored.data$best.hit,invert = T)] = ""

   n.decoy = length(grep("decoy",dec.db.data$formula))
   k = n.decoy / (length(dec.db.data$formula) - n.decoy)

   for (i in 1:nrow(scored.data)) {
      index.score[i] = i
      PSIR[i] = sum(hit.type[1:i] == "decoy") / k
      FDR[i] = PSIR[i] / index.score[i]

   }

   scored.data = cbind(scored.data,hit.type,PSIR,FDR)
   #### Minorize the FDR estimates
   if (score.method == "best") {
      u.score = scored.data$score.Q
   }else{
      ## assuming this is monotone decreasing
      u.score = scored.data$delta.score
   }

   u.fdr = FDR
   u.fdr.m = u.fdr
   is.minor.point = rep(TRUE, length(u.score))
   for (i in 1:(length(u.score) - 1)) {
      tmpid = (i + 1):length(u.score)
      if (any(u.fdr[tmpid] <= u.fdr[i])) {
         is.minor.point[i] = FALSE
      }
   }

   is.minor.point[length(u.score)] = FALSE

   #### First minor point to the last minor point
   mp.id = which(is.minor.point)
   n.mp = length(mp.id)
   u.fdr.m[1:(mp.id[1])] = u.fdr[mp.id[1]]

   for (k in 1:(n.mp - 1)) {
      u.fdr.m[(mp.id[k] + 1):(mp.id[k + 1])] = u.fdr[mp.id[k + 1]]
   }


   u.fdr.m[(mp.id[n.mp]):(length(u.score))] = u.fdr[mp.id[n.mp]]
   scored.data$FDR = u.fdr.m

   dec.scored.data$scored.data = scored.data


   dec.scored.data$unscored.data$hit.type =
   dec.scored.data$unscored.data$PSIR =
   dec.scored.data$unscored.data$FDR = NA
   dec.scored.data

}


#' compute.thresholds
#'
#' @param dec.scored.data
#'
#' @return
#' @export
#'
#' @keywords internal
compute.thresholds = function(dec.scored.data){
   cuts =seq(0.001,0.20,by=0.001)

   cutoffs = sapply(cuts,function(x)  min(dec.scored.data$score.Q[which(dec.scored.data$FDR<x)]))
   names(cutoffs) = paste("alpha",cuts,sep="_")
   cutoffs
}



