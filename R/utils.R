.LOG <- function(..., file, print = TRUE) {
  msg <- paste0(..., "\n")
  if (print) cat(msg)
  cat(msg, file = file, append = TRUE)
}

.get_renamed_colnames <- function(hold_names, userprovided, checkforsingle=c(), filename, N_provided, log.file,
                                  warnz=FALSE, warn_for_missing=c(), stop_on_missing=c(), utilfuncs=NULL) {
  interpreted_names <- list(
    SNP=c("SNP","SNPID","RSID","RS_NUMBER","RS_NUMBERS", "MARKERNAME", "ID","PREDICTOR","SNP_ID"),
    A1=c("A1", "ALLELE1","EFFECT_ALLELE","INC_ALLELE","REFERENCE_ALLELE","EA","REF"),
    A2=c("A2","ALLELE2","ALLELE0","OTHER_ALLELE","REF","NON_EFFECT_ALLELE","DEC_ALLELE","OA","NEA", "ALT"),
    effect=c("OR","B","BETA","LOG_ODDS","EFFECTS","EFFECT","SIGNED_SUMSTAT","EST", "BETA1", "LOGOR"),
    INFO=c("INFO", "IMPINFO"),
    P=c("P","PVALUE","PVAL","P_VALUE","P-VALUE","P.VALUE","P_VAL","GC_PVALUE","WALD_P"),
    N=c("N","WEIGHT","NCOMPLETESAMPLES", "TOTALSAMPLESIZE", "TOTALN", "TOTAL_N","N_COMPLETE_SAMPLES", "SAMPLESIZE", "NEFF", "N_EFF", "N_EFFECTIVE"),
    MAF=c("MAF", "CEUAF", "FREQ1", "EAF", "FREQ1.HAPMAP", "FREQALLELE1HAPMAPCEU", "FREQ.ALLELE1.HAPMAPCEU", "EFFECT_ALLELE_FREQ", "FREQ.A1"),
    Z=c("Z", "ZSCORE", "Z-SCORE", "ZSTATISTIC", "ZSTAT", "Z-STATISTIC"),
    SE=c("STDERR", "SE", "STDERRLOGOR")
  )
  full_names <- list(
    P="P-value",
    A1="effect allele",
    A2="other allele",
    effect="beta or effect",
    SNP="rs-id",
    SE="standard error"
  )
  if (!is.null(utilfuncs)) {
    for (j in names(utilfuncs)) {
        assign(j, utilfuncs[[j]], envir=environment())
    }
  }
  if (N_provided) {
    interpreted_names[["N"]] <- NULL
  } else {
    if ("NEFF" %in% hold_names) {
      .LOG("Found an NEFF column for sample size. \n
Please note that this is likely effective sample size and should only be used for liability h^2 conversion for binary traits and that it should reflect the sum of effective sample sizes across cohorts.\n
Be aware that some NEFF columns reflect half of the effective sample size, in which case sample size values should be doubled prior to running munge.", file=log.file)
    }
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
    } else if ((col == "effect")){
      if (any(interpreted_names[["Z"]] %in% hold_names)) {
        if (!warnz) {
          .LOG("Interpreting the ", hold_names[hold_names %in% interpreted_names[["Z"]] ] , " column as the ",col," column.",file=log.file)
          hold_names[hold_names %in% interpreted_names[["Z"]] ] <- col
        } else {
          .LOG("There appears to be a Z-statistic column in the summary statistic file ", filename, ". Please set linprob to TRUE for binary traits or OLS to true for continuous traits in order to back out the betas or if betas are already available remove this column.", print=FALSE, file=log.file)
          warning(paste0("There appears to be a Z-statistic column in the summary statistic file ", filename, ". Please set linprob to TRUE for binary traits or OLS to true for continuous traits in order to back out the betas or if betas are already available remove this column."))
        }
      }
    } else {
      if (col %in% warn_for_missing) {
        .LOG('Cannot find ', col, ' column, try renaming it to ', col, ' in the summary statistics file for:',filename,file=log.file)
      } else if (col %in% stop_on_missing) {
        stop(paste0('Cannot find ', col, ' column, try renaming it to ', col, ' in the summary statistics file for:',filename))
      }
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

