sumstats <- function(filenames,reference,trait.names=NULL,se.logit=NULL,model=NULL,prev=NA,N=NULL,info.filter=.6,maf.filter=0.01,keep.indel=FALSE,parallel=TRUE,cores=NULL,log.prefix='sumstats') {
## Pre-process summary statistics for models including individual variants

  begin.time <- Sys.time()

  model.types <- c( 'LIN', 'LOG', 'LPM', 'OLS' )

  n.traits <- length(filenames)

  if (is.null(trait.names)) {
    trait.names <- paste0("trait.",1:n.traits)
  } else {
    if ( length(trait.names) < n.traits ) {
      cat("Please specify 'trait.names' for all traits.\n")
      return(NULL)
    }
  }

  if (is.null(model))
    model=rep('LIN', n.traits)
  if (sum(model %in% model.types) < n.traits) {
    cat("Please specify one of", paste( model.types, collapse=', ' ), "for each trait.\n")
    return(NULL)
  }

  if (any(model=='LOG')) {
    if (is.null(se.logit) || any(!is.finite(se.logit[model=='LOG']))) {
      cat("Please specify one 'se.logit' for each 'LOG' model.\n")
      return(NULL)
    }
  }

  if (any(model=='LPM')) {
    if (is.null(prev) || any(!is.finite(prev[model=='LPM']))) {
      cat("Please specify one numeric 'prev' for each LPM trait (set the others, e.g., to NA).\n")
      return(NULL)
    }
  }

  if (any(model %in% c('LPM','OLS'))) {
    if (is.null(N) || any(!is.finite(N[model %in% c('LPM','OLS')]))) {
      cat("Please specify a finite 'N' for each LPM or OLS trait.\n")
      return(NULL)
    }
  }

  # set number of cores
  if ( is.null(cores) ) {
    ## if no default provided, use 1 less than the total number of cores available so the computer
    ## will still function
    cores <- detectCores() - 1
    if ( cores > n.traits ) {
      cores <- n.traits
    }
  }

  if ( parallel ) {
    split.log = FALSE  ## no standard output if processing all in parallel
    cat( 'note: standard output is partially suppressed for parallel processing.\n' )
  } else {
    cores = 1  ## this defaults mclapply() to lapply()
    split.log = TRUE  ## preserve standard output if processing all sequentially
  }


  .sumstats.process <- function(X) {
  ## Work horse for the mclapply() call

    i = X

    log.file <- paste0( log.prefix, '_', make.names( trait.names[i] ), '.log' )
    sink(log.file, append=FALSE, split=split.log)

    names.beta <- paste0("beta.",trait.names)
    names.se <- paste0("se.",trait.names)

    writeLines( strwrap( paste0(
      "Reading summary statistics for '", filenames[i], "'. ",
      "Please note that this step usually takes a few minutes due to the size of the summary ",
      "statistics files."
    ) ) )
    ##note that fread is not used here as we have observed different formatting for headers causing mismatched columns
    tmpfile <- read.table( filenames[i], header=T, quote="\"", fill=T, na.string=c(".",NA,"NA","") )

    writeLines( strwrap( paste0(
      "Pre-processing summary statistics for '", filenames[i], "'."
    ) ) )

    ss.tags <- GenomicSEM:::.sumstats.tags()
    names(tmpfile) <- GenomicSEM:::.sumstats.parse.header( tmpfile, maf.override=T )

    # Print a warning message when multiple columns are interpreted as p-values, effects, rsIDs or alleles columns
    if(sum(names(tmpfile) %in% "P") > 1) warning(paste0(
      'Multiple columns are being interpreted as the P-value column. Try renaming the column you dont want interpreted"
      "as P, e.g., to P2 for:', filenames[i]))
    if(sum(names(tmpfile) %in% "effect") > 1) warning(paste0(
      'Multiple columns are being interpreted as the effect column. Try renaming the column you dont want interpreted"
      "as effect, e.g., to effect2 for:', filenames[i]))
    if(sum(names(tmpfile) %in% "SNP") > 1) warning(paste0(
      'Multiple columns are being interpreted as the variant ID column. Try renaming the column you dont want"
      "interpreted as variant ID, e.g., to SNP2 for:', filenames[i]))
    if(sum(names(tmpfile) %in% "A1") > 1) warning(paste0(
      'Multiple columns are being interpreted as the effect allele column. Try renaming the column you dont want"
      "interpreted as effect allele, e.g., to A1_2 for:', filenames[i]))
    if(sum(names(tmpfile) %in% "A2") > 1) warning(paste0(
      'Multiple columns are being interpreted as the non-effect allele column. Try renaming the column you dont want"
      "interpreted as non-effect allele, e.g., to A2_2 for:', filenames[i]))

    if ( model[i] == 'LOG' ) {
      if ( length( intersect( names(tmpfile), ss.tags$zscore ) ) > 0 ) {
        warning( paste0(
          "There appears to be a Z-statistic column in the summary statistic file for ", trait.names[i], ".",
          "Transformations for case/control traits require an OR or a logistic beta column. Please remove/replace the",
          "Z-statistic column"
        ) )
        writeLines( strwrap( paste0(
          "There appears to be a Z-statistic column in the summary statistic file for ", trait.names[i], ".",
          "Transformations for case/control traits require an OR or a logistic beta column. Please remove/replace the",
          "Z-statistic column"
        ) ) )
      }
    }

    # Compute N as N cases plus N controls if reported:
    if ( "N_CAS" %in% colnames(tmpfile) && "N_CON" %in% colnames(tmpfile) ) {
      tmpfile$N <- tmpfile$N_CAS + tmpfile$N_CON
      cat(
        "As the file includes both N_CAS and N_CON columns, the sum of these two columns will be",
        "used as the total sample size.",
        sep="\n"
      )
    }

    if ( !is.null(N) && !is.na(N[i]) ) {
      if ( "N" %in% colnames(tmpfile) ) {
        writeLines( strwrap( paste(
          "As the summary statistics file includes a sample size column, this is being used in",
          "place of the user-provided sample size. If you wish to still use the provided sample",
          "size, instead of the sample size reported in the summary statistics file, please change",
          "the sample size column header to N2, for example. However, note that the sample size",
          "is only used for LPM and OLS conversions."
        ) ) )
      } else {
        tmpfile$N<-N[i]
      }
    }

    if ( keep.indel ) {
      tmpfile$A1 <- factor(toupper(tmpfile$A1))
      tmpfile$A2 <- factor(toupper(tmpfile$A2))
      cat(
        "Keeping variants other than SNPs: this may cause problems when aligning alleles across",
        "trait summary statistics and reference file.",
        sep="\n"
      )
    } else {
      ##make sure all alleles are upper case for matching
      tmpfile$A1 <- factor(toupper(tmpfile$A1), c("A", "C", "G", "T"))
      tmpfile$A2 <- factor(toupper(tmpfile$A2), c("A", "C", "G", "T"))
#       cat(
#         "Removing variants other than SNPs, e.g., indels. To change this behavior set",
#         "keep.indel=TRUE in the function call.",sep="\n"
#       )
    }

    ##merge with reference file
    writeLines( strwrap( paste("Merging file:", filenames[i], "with the reference file:", reference) ) )
    b<-nrow(tmpfile)
    writeLines( strwrap( paste(b, "rows present in the full", filenames[i], "summary statistics file.") ) )
    tmpfile <- suppressWarnings(inner_join(ref,tmpfile,by="SNP",all.x=F,all.y=F))
    writeLines( strwrap( paste(
      (b-nrow(tmpfile)), "rows were removed from the", filenames[i], "summary statistics file as",
      "the rsIDs for these SNPs were not present in the reference file."
    ) ) )

    ##remove any rows with missing p-values
    b<-nrow(tmpfile)
    if ("P" %in% colnames(tmpfile)) {
      tmpfile<-subset(tmpfile, !(is.na(tmpfile$P)))
    }
    if (b-nrow(tmpfile) > 0) writeLines( strwrap( paste(
      b-nrow(tmpfile), "rows were removed from the", filenames[i], "summary statistics file due to",
      "missing values in the P-value column."
    ) ) )

    ##check that p-value column does not contain an excess of 1s/0s
    if ((sum(tmpfile$P > 1) + sum(tmpfile$P < 0)) > 100) writeLines( strwrap(
      "More than 100 SNPs have P-value above 1 or below 0. The P column may be mislabled!"
    ) )

    ##remove any rows with missing effects
    b<-nrow(tmpfile)
    if ("EFFECT" %in% colnames(tmpfile)) {
      tmpfile<-subset(tmpfile, !(is.na(tmpfile$EFFECT)))
    }
    if (b-nrow(tmpfile) > 0) writeLines( strwrap( paste(
      b-nrow(tmpfile), "rows were removed from the", filenames[i], "summary statistics file due to",
      "missing values in the effect column."
    ) ) )

    ##remove SNPs that don't match A1 or A2 in the reference file.
    b<-nrow(tmpfile)
    tmpfile<-subset(tmpfile, !(tmpfile$A1.x != (tmpfile$A1.y) & tmpfile$A1.x != (tmpfile$A2.y)))
    if (b-nrow(tmpfile) > 0) writeLines( strwrap( paste(
      b-nrow(tmpfile), "row(s) were removed from the", filenames[i], "summary statistics file due",
      "to the effect allele (A1) column not matching A1 or A2 in the reference file."
    ) ) )
    b<-nrow(tmpfile)
    tmpfile<-subset(tmpfile, !(tmpfile$A2.x != (tmpfile$A2.y) & tmpfile$A2.x != (tmpfile$A1.y)))
    if (b-nrow(tmpfile) > 0) writeLines( strwrap( paste(
      b-nrow(tmpfile), "row(s) were removed from the", filenames[i], "summary statistics file due",
      "to the other allele (A2) column not matching A1 or A2 in the reference file."
    ) ) )

    ##check the presence of SE, needed for some of GenomicSEM's functionality
    if (sum(colnames(tmpfile) %in% "SE") == 0) {
      writeLines( strwrap( paste( 'Cannot find a SE column, try renaming it SE in', filenames[i] ) ) )
      warning( paste( 'Cannot find a SE column, try renaming it SE in', filenames[i] ) )
      tmpfile$SE = NA
    }

    ##check the presence of INFO
    if ("INFO" %in% colnames(tmpfile)) {
      b<-nrow(tmpfile)
      tmpfile <- tmpfile[tmpfile$INFO >= info.filter,]
      writeLines( strwrap( paste0(
        b-nrow(tmpfile), "rows were removed from the '", filenames[i], "' summary statistics due ",
        "to INFO values below the designated threshold of ", info.filter, "."
      ) ) )
    } else {
      writeLines( strwrap( "No INFO column, cannot filter on INFO, which may influence results." ) )
    }

    ##determine whether it is OR or logistic/continuous effect based on median effect size
    EFFECT.untransf<-tmpfile$EFFECT[[1]]
    tmpfile$EFFECT<-ifelse(
      rep(round(median(tmpfile$EFFECT,na.rm=T)) == 1,nrow(tmpfile)),
      log(tmpfile$EFFECT),
      tmpfile$EFFECT
    )
    EFFECT.transf<-tmpfile$EFFECT[[1]]
    if (EFFECT.untransf == EFFECT.transf) {
      writeLines( strwrap( paste(
        "The values in the effect column were NOT identified as odds ratios (OR) for the", filenames[i],
        "summary statistics file based on the median of the effect column being close to 0.",
        "Please, make sure that this was the correct interpretation."
      ) ) )
    } else {
      writeLines( strwrap( paste(
        "The values in the effect column were identified as odds ratios (OR) for the", filenames[i],
        "summary statistics file based on the median of the effect column being close to 1.",
        "Please, make sure that this was the correct interpretation."
      ) ) )
    }

    ##flip effect to match ordering in reference file
    tmpfile$EFFECT <- ifelse(
      tmpfile$A1.x != (tmpfile$A1.y) & tmpfile$A1.x == (tmpfile$A2.y),
      -tmpfile$EFFECT,
      tmpfile$EFFECT
    )

    ##remove any rows printed as exactly 0
    b<-nrow(tmpfile)
    if("EFFECT" %in% colnames(tmpfile)) {
      tmpfile<-subset(tmpfile, tmpfile$EFFECT != 0)
    }
    if(b-nrow(tmpfile) > 0) writeLines( strwrap( paste(
      b-nrow(tmpfile), "rows were removed from the", filenames[i], "summary statistics file due to",
      "effect values estimated as exactly 0, which cause problems for matrix inversion in later",
      "Genomic SEM analyses.") ) )

    ##compute the variants' variances
    varSNP <- 2 * tmpfile$MAF * (1-tmpfile$MAF)

    if (model[i] == 'LOG') {
      writeLines( strwrap( paste( "LOG transformation is being used for file:", filenames[i] ) ) )
      ##correction factor for logistic regression coefficients
      logCF <- 1 / ((tmpfile$EFFECT^2) * varSNP + (pi^2)/3)^.5
      tmpfile$SE <- tmpfile$SE * logCF
      if (se.logit[i]) {
        writeLines( strwrap( paste(
          "Performing transformation under the assumption that the effect column is an odds ratio",
          "or a logarithm of one (please see output above to determine how it was interpreted) and",
          "the SE column is the standard error of the latter (i.e., in the logistic scale) for:",
          filenames[i]
        ) ) )
      } else {
        writeLines( strwrap( paste(
          "Performing transformation under the assumption that the effect column is an odds ratio",
          "or a logarithm of one (please see output above to determine how it was interpreted) and",
          "the SE column is the standard error of the former (i.e. NOT in the logistic scale) for:",
          filenames[i]
        ) ) )
        tmpfile$SE <- tmpfile$SE / exp(tmpfile$EFFECT)
      }
      tmpfile$EFFECT <- tmpfile$EFFECT * logCF
    }

    if (model[i] == 'LPM') {
      writeLines( strwrap( paste( "LPM transformation is being used for file:", filenames[i] ) ) )
      tmpfile$Z <- sign(tmpfile$EFFECT) * sqrt(qchisq(tmpfile$P,1,lower=F))
      if ("N" %in% colnames(tmpfile) && is.finite(prev[i])) {
        tmpfile$SE <- 1 / sqrt( prev[i]*(1-prev[i]) * tmpfile$N * varSNP )
        tmpfile$EFFECT <- tmpfile$Z * tmpfile$SE
        ##correction factor for logistic regression coefficients (uses the re-estimated effect)
        logCF <- 1 / ((tmpfile$EFFECT^2) * varSNP + (pi^2)/3)^.5
        tmpfile$SE <- tmpfile$SE * logCF
        tmpfile$EFFECT <- tmpfile$EFFECT * logCF
      } else {
        writeLines( strwrap( paste(
          "ERROR: sample size (N) and prevalence (prev) are needed for LPM Standardization.",
          "Please, provide prev as argument to the function and either do the same for N or try",
          "changing the name of the sample size column to N."
        ) ) )
      }
    }

    if (model[i] == 'OLS') {
      writeLines( strwrap( paste( "OLS transformation is being used for file:", filenames[i] ) ) )
      tmpfile$Z <- sign(tmpfile$EFFECT) * sqrt(qchisq(tmpfile$P,1,lower=F))
      if ("N" %in% colnames(tmpfile)) {
        tmpfile$SE <- 1 / sqrt( tmpfile$N * varSNP )
        tmpfile$EFFECT <- tmpfile$Z * tmpfile$SE
      } else {
        writeLines( strwrap( paste(
          "ERROR: the sample size (N) is needed for OLS Standardization. Please, either provide",
          "the sample size N as argument to the function or try changing the name of the sample",
          "size column to N."
        ) ) )
      }
    }

    output <- cbind.data.frame( tmpfile$SNP,
                                tmpfile$EFFECT,
                                tmpfile$SE )
    output <- na.omit(output)
    output <- output[apply(output!=0, 1, all),]
    colnames(output) <- c("SNP",names.beta[i],names.se[i])

    writeLines( strwrap( paste(
      nrow(output), "SNPs are left from the summary statistics file", filenames[i], "after QC and",
      "merging with the reference file."
    ) ) )

    if ( mean( abs( output[,names.beta[i]] / output[,names.se[i]] ) ) > 5 ) {
      writeLines( strwrap( paste0(
        "WARNING: The average value of estimate over standard error (i.e., Z) is > 5 for ", trait.names[i],
        ". This suggests a column was misinterpreted or arguments were misspecified. Please, post on the google group",
        "if you are unable to figure this issue out yourself."
      ) ) )
      warning( paste0( 
        "The average value of estimate over standard error (i.e., Z) is > 5 for ", trait.names[i],
        ". This suggests a column was misinterpreted or arguments were misspecified. Please, post on the google group",
        "if you are unable to figure this issue out yourself."
      ) )
    }

    cat( '\n\n' )

    sink()

    output

  }

  writeLines( strwrap( paste0(
    "Please note that the files should be in the same order in which they were listed in the ldsc function call."
  ) ) )

  cat("Reading in reference file.\n")
  ref <- fread(file=reference,header=T,data.table=F)

  ##filter reference file on user provided maf.filter
  cat("Applying MAF filer of", maf.filter, "to the reference file.\n")
  ref<-subset(ref, ref$MAF >= maf.filter)

  data.frame.out <- ref

  data.frame.tmp <- mclapply(1:n.traits, .sumstats.process, mc.cores=cores)

  for ( i in 1:length(data.frame.tmp) ) data.frame.out <-
    suppressWarnings(inner_join(data.frame.out,data.frame.tmp[[i]],by="SNP",all.x=F,all.y=F))

  b<-nrow(data.frame.out)
  #data.frame.out<-data.frame.out[!duplicated(data.frame.out$BP),]

  end.time <- Sys.time()

  total.time <- difftime(time1=end.time,time2=begin.time,units="sec")
  mins <- floor(floor(total.time)/60)
  secs <- total.time-mins*60

  writeLines( strwrap( paste0( 
    "After merging all summary statistics using listwise deletion, performing QC, and merging ",
    "with the reference, ", nrow(data.frame.out), " SNPs are left in the final multivariate ",
    "summary statistics file."
  ) ) )
  writeLines( strwrap( paste0( 
    "Finished running at ", end.time, "."
  ) ) )
  writeLines( strwrap( paste0( 
    "Running sumstats for all files took ", mins," minutes and ", secs," seconds."
  ) ) )
  writeLines( strwrap( paste0( 
    "Please check the log file(s) '", log.prefix, "_*.log' to ensure that all columns were ",
    "interpreted correctly and no warnings were issued for any of the summary statistics files."
  ) ) )

  data.frame.out

}

