.LOG <- function(..., file, print = TRUE) {
  msg <- paste0(..., "\n")
  if (print) cat(msg)
  cat(msg, file = file, append = TRUE)
}

.get_renamed_colnames <- function(hold_names, userprovided, checkforsingle=c(), filename, N_provided, log.file) {
  interpreted_names <- list(
    SNP=c("SNP","SNPID","RSID","RS_NUMBER","RS_NUMBERS", "MARKERNAME", "ID","PREDICTOR","SNP_ID"),
    A1=c("A1", "ALLELE1","EFFECT_ALLELE","INC_ALLELE","REFERENCE_ALLELE","EA","REF"),
    A2=c("A2","ALLELE2","ALLELE0","OTHER_ALLELE","REF","NON_EFFECT_ALLELE","DEC_ALLELE","OA","NEA", "ALT"),
    effect=c("OR","B","BETA","LOG_ODDS","EFFECTS","EFFECT","SIGNED_SUMSTAT", "Z","ZSCORE","EST","ZSTAT","ZSTATISTIC", "BETA1", "LOGOR"),
    INFO=c("INFO", "IMPINFO"),
    P=c("P","PVALUE","PVAL","P_VALUE","P-VALUE","P.VALUE","P_VAL","GC_PVALUE","WALD_P"),
    N=c("N","WEIGHT","NCOMPLETESAMPLES", "TOTALSAMPLESIZE", "TOTALN", "TOTAL_N","N_COMPLETE_SAMPLES", "SAMPLESIZE", "NEFF", "N_EFF", "N_EFFECTIVE"),
    N_CAS=c("NCASE","N_CASE","N_CASES","N_CAS", "NCAS", "NCA"),
    N_CON=c("NCONTROL","N_CONTROL","N_CONTROLS","N_CON","CONTROLS_N", "NCON", "NCO"),
    MAF=c("MAF", "CEUAF", "FREQ1", "EAF", "FREQ1.HAPMAP", "FREQALLELE1HAPMAPCEU", "FREQ.ALLELE1.HAPMAPCEU", "EFFECT_ALLELE_FREQ", "FREQ.A1"),
    Z=c("ZSCORE", "Z-SCORE", "ZSTATISTIC", "Z-STATISTIC")
  )
  full_names <- list(
    P="P-value",
    A1="effect allele",
    A2="other allele",
    effect="beta or effect",
    SNP="rs-id"
  )
  if (N_provided) {
    interpreted_names[["N"]] <- NULL
  }
  for (col in names(interpreted_names)) {
    if (col %in% names(userprovided)) {
      .LOG("Interpreting the ",userprovided[[col]]," column as the ",col, " column, as requested",file=log.file)
      hold_names[ hold_names == toupper(userprovided[[col]]) ] <- col
    } else if (col %in% hold_names) {
      .LOG("Interpreting the ",col," column as the ",col, " column.",file=log.file)
    } else if (any(interpreted_names[[col]] %in% hold_names)) {
      .LOG("Interpreting the ", hold_names[ hold_names %in% interpreted_names[[col]] ], " column as the ",col," column.",file=log.file)
      hold_names[ hold_names %in% interpreted_names[[col]] ] <- col
    } else {
      # Print a message for missing columns
      .LOG('Cannot find ', col, ' column, try renaming it to ', col, ' in the summary statistics file for:',filename,file=log.file)
    }
  }
  # Print log and throw warning messages if multiple or no columns were found for those specified in checkforsingle
  if (length(checkforsingle) > 0) {
    for (col in checkforsingle) {
      if(sum(hold_names == col) == 0) {
        .LOG('Cannot find ',full_names[[col]],' column, try renaming it ', col, ' in the summary statistics file for:',filename,file=log.file)
        warning(paste0('Cannot find ',full_names[[col]],' column, try renaming it ', col, ' in the summary statistics file for:', filename))
      }
      if(sum(hold_names == col) > 1) {
        .LOG('Multiple columns are being interpreted as the ',full_names[[col]],' column, try renaming the column you dont want interpreted to ', col, '2 in the summary statistics file for:',filename,file=log.file)
        warning(paste0('Multiple columns are being interpreted as the ',full_names[[col]],' column, try renaming the column you dont want interpreted to ', col, '2 in the summary statistics file for:', filename))
      }
    }
  }
  return(hold_names)
}
#function to rearrange the sampling covariance matrix from original order to lavaan's order:
#'k' is the number of variables in the model
#'fit' is the fit function of the regression model
#'names' is a vector of variable names in the order you used
.rearrange <- function (k, fit, names) {
    order1 <- names
    order2 <- rownames(inspect(fit)[[1]]) #order of variables
    kst <- k*(k+1)/2
    covA <- matrix(NA, k, k)
    covA[lower.tri(covA, diag = TRUE)] <- 1:kst
    covA <- t(covA)
    covA[lower.tri(covA, diag = TRUE)] <- 1:kst
    colnames(covA) <- rownames(covA) <- order1 #give A actual variable order from lavaan output
    #reorder A by order2
    covA <- covA[order2, order2] #rearrange rows/columns
    vec2 <- lav_matrix_vech(covA) #grab new vectorized order
    return(vec2)
}

##modification of trycatch that allows the results of a failed run to still be saved
.tryCatch.W.E <- function(expr) {
    W <- NULL
    w.handler <- function(w){ # warning handler
      W <<- w
      invokeRestart("muffleWarning")
    }
    list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                     warning = w.handler), warning = W)
}

.get_V_full <- function(k, V_LD, varSNPSE2, V_SNP) {
    ##create shell of full sampling covariance matrix
    V_Full<-diag(((k+1)*(k+2))/2)

    ##input the ld-score regression region of sampling covariance from ld-score regression SEs
    V_Full[(k+2):nrow(V_Full),(k+2):nrow(V_Full)]<-V_LD

    ##add in SE of SNP variance as first observation in sampling covariance matrix
    V_Full[1,1]<-varSNPSE2

    ##add in SNP region of sampling covariance matrix
    V_Full[2:(k+1),2:(k+1)]<-V_SNP
    return(V_Full)
}

.get_V_SNP <- function(SE_SNP, I_LD, varSNP, GC, coords, k, i) {
     V_SNP<-diag(k)
    #loop to add in the GWAS SEs, correct them for univariate and bivariate intercepts, and multiply by SNP variance from reference panel
    if(GC == "conserv"){
      for (p in 1:nrow(coords)) {
        x<-coords[p,1]
        y<-coords[p,2]
        if (x != y) {
          V_SNP[x,y]<-(SE_SNP[i,y]*SE_SNP[i,x]*I_LD[x,y]*I_LD[x,x]*I_LD[y,y]*varSNP[i]^2)}
        if (x == y) {
          V_SNP[x,x]<-(SE_SNP[i,x]*I_LD[x,x]*varSNP[i])^2
        }
      }
    }

    if(GC == "standard"){
      for (p in 1:nrow(coords)) {
        x<-coords[p,1]
        y<-coords[p,2]
        if (x != y) {
          V_SNP[x,y]<-(SE_SNP[i,y]*SE_SNP[i,x]*I_LD[x,y]*sqrt(I_LD[x,x])*sqrt(I_LD[y,y])*varSNP[i]^2)}
        if (x == y) {
          V_SNP[x,x]<-(SE_SNP[i,x]*sqrt(I_LD[x,x])*varSNP[i])^2
        }
      }
    }

    if(GC == "none"){
      for (p in 1:nrow(coords)) {
        x<-coords[p,1]
        y<-coords[p,2]
        if (x != y) {
          V_SNP[x,y]<-(SE_SNP[i,y]*SE_SNP[i,x]*I_LD[x,y]*varSNP[i]^2)}
        if (x == y) {
          V_SNP[x,x]<-(SE_SNP[i,x]*varSNP[i])^2
        }
      }
    }
    return(V_SNP)
}

.get_Z_pre <- function(i, beta_SNP, SE_SNP, I_LD, GC) {
    if(GC == "conserv"){
        Z_pre<-beta_SNP[i,]/(SE_SNP[i,]*diag(I_LD))
    }
    if(GC=="standard"){
        Z_pre<-beta_SNP[i,]/(SE_SNP[i,]*sqrt(diag(I_LD)))
    }
    if(GC=="none"){
        Z_pre<-beta_SNP[i,]/SE_SNP[i,]
    }
    return(Z_pre)
}

