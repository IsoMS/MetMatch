#' Load the CAMERA output file into a list of query objects
#'
#' @param fh file handle
#' @param N N samples
#'
#' @return Q class
#' @export
#'
convert.camera.file  = function(fh,N){
   d = read.csv(fh)
   e = d[grepl("M",d$isotopes)==FALSE,]
   e = e[order(e$mz),]
   z1.e = sapply(e$mz, function(x) convertIonToNeutralMass("M+H",x))
   z2.e = sapply(e$mz, function(x) convertIonToNeutralMass("M+2H",x))


   d = d[grep("M",d$isotopes),]
   d = d[order(d$isotopes),]
   #write.csv(file=sub(".cam.csv",".isocam.csv",fh),d,row.names=FALSE);
   #write.csv(file=sub(".cam.csv",".mzcam.csv",fh),cbind(e,z1=z1.e,z2=z2.e),row.names=FALSE);

   isotopes = as.character(d$isotopes)
   id= unlist( lapply( isotopes, function(x) strsplit( x,"\\[")[[1]][2] ))
   id = as.numeric(sub("\\]","",id) )
   u.id = unique(id)
   ions = unlist( lapply( isotopes, function(x) strsplit( x,"\\[")[[1]][3] ))
   ions = fix.camera.ions(ions)

   d$isotopes = ions
   d = cbind(d,id)
   Q = list()

   for(i in u.id){
      pattern = paste("^",i,"$",sep="")
      tm = d[grepl(pattern,d$id),]
      if(N>1){
         int = apply(tm[,10:(10+N-1)],1,mean)
      }else{
         int = tm$maxo
      }

      sumInt = sum(int)
      int = int / sumInt
      type = tm$isotopes[1]

      M = sapply(tm$mz , function(x) convertIonToNeutralMass(type,x))

      Q[[i]] = list(File = fh,Ion = type, Mass = M[1], RT = tm$rt[1], mz = tm$mz[1], sumInt = sumInt, Nisotope = length(int),
                    IsoMass = M,IsoInt = int )

   }

   Q

}


#' Fix CAMERA (xcms) annotations
#'
#' Helper function for loading CAMERA annotations
#'
#' @param ions the input ion type
#' @export
#' @keywords internal
fix.camera.ions = function(ions){

   ions = sub("M\\]\\+","M+H",ions)
   ions = sub("M\\+1\\]\\+","M+H+1",ions)
   ions = sub("M\\+2\\]\\+","M+H+2",ions)
   ions = sub("M\\+3\\]\\+","M+H+3",ions)
   ions = sub("M\\+4\\]\\+","M+H+4",ions)
   ions = sub("M\\+5\\]\\+","M+H+5",ions)
   ions = sub("M\\+6\\]\\+","M+H+6",ions)

   ions = sub("M\\]2\\+","M+2H",ions)
   ions = sub("M\\+1\\]2\\+","M+2H+1",ions)
   ions = sub("M\\+2\\]2\\+","M+2H+2",ions)
   ions = sub("M\\+3\\]2\\+","M+2H+3",ions)
   ions = sub("M\\+4\\]2\\+","M+2H+4",ions)
   ions = sub("M\\+5\\]2\\+","M+2H+5",ions)
   ions = sub("M\\+6\\]2\\+","M+2H+6",ions)


   ions = sub("M\\]3\\+","M+3H",ions)
   ions = sub("M\\+1\\]3\\+","M+3H+1",ions)
   ions = sub("M\\+2\\]3\\+","M+3H+2",ions)
   ions = sub("M\\+3\\]3\\+","M+3H+3",ions)
   ions = sub("M\\+4\\]3\\+","M+3H+4",ions)
   ions = sub("M\\+5\\]3\\+","M+3H+5",ions)
   ions = sub("M\\+6\\]3\\+","M+3H+6",ions)

   #negative
   ions = sub("M\\]\\-","M-H",ions)
   ions = sub("M\\+1\\]\\-","M-H+1",ions)
   ions = sub("M\\+2\\]\\-","M-H+2",ions)
   ions = sub("M\\+3\\]\\-","M-H+3",ions)
   ions = sub("M\\+4\\]\\-","M-H+4",ions)
   ions = sub("M\\+5\\]\\-","M-H+5",ions)
   ions = sub("M\\+6\\]\\-","M-H+6",ions)

   ions = sub("M\\]2\\-","M-2H",ions)
   ions = sub("M\\+1\\]2\\-","M-2H+1",ions)
   ions = sub("M\\+2\\]2\\-","M-2H+2",ions)
   ions = sub("M\\+3\\]2\\-","M-2H+3",ions)
   ions = sub("M\\+4\\]2\\-","M-2H+4",ions)
   ions = sub("M\\+5\\]2\\-","M-2H+5",ions)
   ions = sub("M\\+6\\]2\\-","M-2H+6",ions)

   ions = sub("M\\]3\\-","M-3H",ions)
   ions = sub("M\\+1\\]3\\-","M-3H+1",ions)
   ions = sub("M\\+2\\]3\\-","M-3H+2",ions)
   ions = sub("M\\+3\\]3\\-","M-3H+3",ions)
   ions = sub("M\\+4\\]3\\-","M-3H+4",ions)
   ions = sub("M\\+5\\]3\\-","M-3H+5",ions)
   ions = sub("M\\+6\\]3\\-","M-3H+6",ions)

   ions

}



#' convertIonsToNeutralMass
#'
#' Converts ion mz to M neutral charge mass.
#'
#' @param ionType the input ion type
#' @param mz the input mz
#' @export
#' @return  the neutral mass
#' @keywords internal
convertIonToNeutralMass = function(ionType,mz){

   #Constants
   PROTON_MASS = 1.007276 # Mass of Proton
   SODIUM_MASS = 22.989218 # Mass of Na+
   POTASSIUM_MASS = 38.963158 # Mass of K+
   NH4_MASS = 18.033823 # Mass of Ammonium

   # Important ions and adducts...negative mode
   STRIP_PROTON_MASS = -1.007276 # Subtraction of proton
   BROMIUM_MASS = 78.918885
   CHLORINE_MASS = 34.969402
   ACETATE_MASS = 0.00
   TRI_FLOURO_ACETATE_MASS = 0.00

   # Neutral loss
   WATER_MASS = 18.010565

   # Positive ion types


   # Helper function
   removeNmer = function(mz,constant,N){
      return ( (mz - constant) /  N )
   }


   # Main Function
   nominalMass = 0
   if(ionType == "M+H"){  # single charge ion
      nominalMass = mz - PROTON_MASS
   }
   if(ionType == "M+2H"){  # double charge
      nominalMass = mz*2 - 2*PROTON_MASS
   }
   if(ionType == "M+3H"){  # double charge
      nominalMass = mz*3 - 3*PROTON_MASS
   }
   if(ionType == "2M+H"){  # single charge ion
      nominalMass = removeNmer(mz, PROTON_MASS,2)
   }
   if(ionType == "3M+H"){  # single charge ion
      nominalMass = removeNmer(mz, PROTON_MASS,3)
   }
   if(ionType == "M+Na"){ # single charge Na adduct
      nominalMass = mz - SODIUM_MASS
   }
   if(ionType == "2M+Na"){ # Dimer with Na adduct
      nominalMass = removeNmer(mz, SODIUM_MASS,2)
   }
   if(ionType == "M+K"){ # single charge K adduct
      nominalMass = mz - POTASSIUM_MASS
   }
   if(ionType == "2M+K"){ # Dimer with K adduct
      nominalMass = removeNmer(mz, POTASSIUM_MASS,2)
   }
   if(ionType == "M+NH4"){ # Monomer with NH4 adduct
      nominalMass = removeNmer(mz,NH4_MASS,1)
   }
   if(ionType == "2M+NH4"){ # Dimer with NH4 adduct
      nominalMass = removeNmer(mz,NH4_MASS,2)
   }
   if(ionType == "M+H+[-H2O]"){  #Water loss
      nominalMass = mz + WATER_MASS - PROTON_MASS
   }

   #Negative ions
   if(ionType == "M-H"){  #Water loss
      nominalMass = mz + PROTON_MASS
   }
   if(ionType == "M-2H"){  #z=2
      nominalMass = mz*2 + 2*PROTON_MASS
   }
   if(ionType == "M-3H"){  #z=3
      nominalMass = mz*3 + 3*PROTON_MASS
   }
   if(ionType == "2M-H"){  #Dimer
      nominalMass = removeNmer(mz,PROTON_MASS,2)
   }
   if(ionType == "3M-H"){  #Trimer
      nominalMass = removeNmer(mz,PROTON_MASS,3)
   }
   if(ionType == "2M+Cl"){  #Chlorine + Dimer
      nominalMass = removeNmer(mz,CHLORINE_MASS,2)
   }
   if(ionType == "2M+HCOO"){  #Water loss
      nominalMass = removeNmer(mz,ACETATE_MASS,2)
   }
   if(ionType == "2M+CF3COO"){  #Water loss
      nominalMass = removeNmer(mz,TRI_FLOURO_ACETATE_MASS,2)
   }
   return (nominalMass)
}


