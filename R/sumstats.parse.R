## Summary statistics parsing utilities

.sumstats.tags <- function() {
  ## Pseudo-class of default tags expected in summary statistics files
  ## It's just a bunch of string arrays for now but allows in theory to define custom methods, too

  m <- list()

  m$variant <- c(
    "snp",
    "SNP",
    "snpid",
    "rsID",
    "SNPID",
    "rsid",
    "RSID",
    "RS_NUMBER",
    "rs_number",
    "RS_NUMBERS",
    "rs_numbers",
    "MarkerName",
    "Markername",
    "markername",
    "MARKERNAME",
    "ID"
  )
  
  m$effect.allele <- c(
    "a1",
    "A1",
    "allele1",
    "Allele1",
    "ALLELE1",
    "alt",
    "ALT",
    "EA",
    "Effect_allele",
    "Effect_Allele",
    "EFFECT_ALLELE",
    "INC_ALLELE"
  )
  
  m$noneffect.allele <- c(
    "a2",
    "A2",
    "allele2",
    "Allele2",
    "ALLELE2",
    "ALLELE0",
    "DEC_ALLELE",
    "NEA",
    "OA",
    "Other_allele",
    "Other_Allele",
    "OTHER_ALLELE",
    "Non_Effect_allele",
    "Non_Effect_Allele",
    "NON_EFFECT_ALLELE",
    "ref",
    "REF"
  )
  
  m$zscore <-c(
    "Z",
    "Zscore",
    "Z.score",
    "Zstatistic",
    "Z.statistic"
  )

  m$effect <- c(
    "b",
    "B",
    "beta",
    "Beta",
    "BETA",
    "est",
    "Effect",
    "Effect_Beta",
    "Effect.Beta",
    "EFFECT",
    "EFFECTS",
    "LOG_ODDS",
    "LOG.ODDS",
    "or",
    "OR",
    "SIGNED_SUMSTAT",
    "SIGNED.SUMSTAT",
    m$zscore
  )

  m$se <- c(
    "se",
    "sebeta",
    "se_beta",
    "se.beta",
    "StdErr",
    "stdErr",
    "stderr",
    "Standard_Err",
    "Standard.Err",
    "Standard_Error",
    "Standard.Error",
    "Standard_err",
    "Standard.err",
    "Standard_error",
    "Standard.error",
    "standard_err",
    "standard.err",
    "standard_error",
    "standard.error",
    "SE",
    "SEBETA",
    "SE_BETA",
    "SE.BETA",
    "SE_Beta",
    "SE.Beta",
    "SE_beta",
    "SE.beta"
  )

  m$info <- c(
    "info",
    "Info",
    "INFO",
    "r2",
    "r_2",
    "R2",
    "R_2",
    "rsq",
    "r_sq",
    "Rsq",
    "R_sq"
  )

  pval <- c(
    "p",
    "P",
    "pval",
    "Pval",
    "PVal",
    "PVAL",
    "p_val",
    "P_val",
    "P_Val",
    "P_VAL",
    "pvalue",
    "Pvalue",
    "PValue",
    "PVALUE",
    "p_value",
    "P_value",
    "P_Value",
    "P_VALUE",
    "p.value",
    "P.value",
    "P.Value",
    "P.VALUE"
  )

  m$pval <- c(
    paste0( 'gc_', pval ),
    paste0( 'GC_', pval ),
    pval
  )

  m$n <- c(
    "N",
    "n",
    "WEIGHT",
    "Weight",
    "weight",
    "nCompleteSamples",
    "n_complete_samples",
    "n.complete.samples",
    "total_sample_size",
    "total.sample.size",
    "TotalSampleSize",
    "TotalN",
    "totalN",
    "Total_N",
    "Total_N",
    "total.N",
    "total.N",
    "Total_n",
    "Total_n",
    "total.n",
    "total.n"
  )

  m$n.cases <- c(
    "NCAS",
    "NCASE",
    "N_CAS",
    "N_CASE",
    "N_CASES",
    "Nca",
    "Ncas",
    "Ncase",
    "Ncases",
    "N_cas",
    "N_case",
    "N_cases",
    "N.cas",
    "N.case",
    "N.cases",
    "n_cas",
    "n_case",
    "n_cases",
    "n.cas",
    "n.case",
    "n.cases"
  )

  m$n.controls <- c(
    "NCON",
    "NCONTROL",
    "N_CON",
    "N_CONTROL",
    "N_CONTROLS",
    "Nco",
    "Ncon",
    "Ncontrol",
    "Ncontrols",
    "N_con",
    "N_control",
    "N_controls",
    "N.con",
    "N.control",
    "N.controls",
    "n_con",
    "n_control",
    "n_controls",
    "n.con",
    "n.control",
    "n.controls"
  )

  m$maf <- c(
    ##Note that this is *supposed* to match any allele frequency field, since
    ##MAF is only used to filter variants or compute allelic/genotypic variance
    "MAF",
    "maf",
    "CEUaf",
    "Freq1",
    "EAF",
    "Freq1.Hapmap",
    "FreqAllele1HapMapCEU",
    "Freq.Allele1.HapMapCEU",
    "EFFECT_ALLELE_FREQ",
    "Freq.A1"
  )

  m <- list2env(m)
  class(m) <- "sumstatTags"
  return(m)

}


.sumstats.parse.header <- function( tmpfile, maf.override=FALSE ) {
  ## Summary statistics file header parser

  ss.tags <- .sumstats.tags()

  hold_names <- names(tmpfile)

  names1<-hold_names
  if("SNP" %in% hold_names) cat( "Interpreting the SNP column as the SNP column.\n" )
  hold_names[hold_names %in% ss.tags$variant] <- "SNP"
  if(length(base::setdiff(names1,hold_names)) > 0)
    cat( paste("Interpreting the", setdiff(names1, hold_names), "column as the SNP column.\n") )

  names1<-hold_names
  if("A1" %in% hold_names) cat( "Interpreting the A1 column as the A1 column.\n" )
  hold_names[hold_names %in% ss.tags$effect.allele] <- "A1"
  if(length(base::setdiff(names1,hold_names)) > 0)
    cat( paste("Interpreting the", setdiff(names1, hold_names), "column as the A1 column.\n") )

  names1<-hold_names
  if("A2" %in% hold_names) cat( "Interpreting the A2 column as the A2 column.\n" )
  hold_names[hold_names %in% ss.tags$noneffect.allele] <- "A2"
  if(length(base::setdiff(names1,hold_names)) > 0)
    cat( paste("Interpreting the", setdiff(names1, hold_names), "column as the A2 column.\n") )

  names1<-hold_names
  if("effect" %in% hold_names) cat( "Interpreting the effect column as the effect column.\n" )
  hold_names[hold_names %in% ss.tags$effect] <- "EFFECT"
  if(length(base::setdiff(names1,hold_names)) > 0)
    cat( paste("Interpreting the", setdiff(names1, hold_names), "column as the effect column.\n") )

  names1<-hold_names
  if("INFO" %in% hold_names) cat( "Interpreting the INFO column as the INFO column.\n" )
  hold_names[hold_names %in% ss.tags$info] <- "INFO"
  if(length(base::setdiff(names1,hold_names)) > 0)
    cat( paste("Interpreting the", setdiff(names1, hold_names), "column as the INFO column.\n") )

  names1<-hold_names
  if("SE" %in% hold_names) cat( "Interpreting the SE column as the SE column.\n" )
  hold_names[hold_names %in% ss.tags$se] <- "SE"
  if(length(base::setdiff(names1,hold_names)) > 0)
    cat( paste("Interpreting the", setdiff(names1, hold_names), "column as the SE column.\n") )

  names1<-hold_names
  if("P" %in% hold_names) cat( "Interpreting the P column as the P column.\n" )
  hold_names[hold_names %in% ss.tags$pval] <- "P"
  if(length(base::setdiff(names1,hold_names)) > 0)
    cat( paste("Interpreting the", setdiff(names1, hold_names), "column as the P column.\n") )

  names1<-hold_names
  if("N" %in% hold_names) cat( "Interpreting the N column as the N (sample size) column.\n" )
  hold_names[hold_names %in% ss.tags$n] <- "N"
  if(length(base::setdiff(names1,hold_names)) > 0)
    cat( paste("Interpreting the", setdiff(names1, hold_names), "column as the N (sample size) column.\n") )

  names1<-hold_names
  if("N_CAS" %in% hold_names) cat( "Interpreting the N_CAS column as the N_CAS (sample size for cases) column.\n" )
  hold_names[hold_names %in% ss.tags$n.cases] <- "N_CAS"
  if(length(base::setdiff(names1,hold_names)) > 0)
    cat( paste("Interpreting the", setdiff(names1, hold_names), "column as the N_CAS (sample size for cases) column.\n") )

  names1<-hold_names
  if("N_CON" %in% hold_names) cat( "Interpreting the N_CON column as the N_CON (sample size for controls) column.\n" )
  hold_names[hold_names %in% ss.tags$n.controls] <- "N_CON"
  if(length(base::setdiff(names1,hold_names)) > 0)
    cat( paste("Interpreting the", setdiff(names1, hold_names), "column as the N_CON (sample size for controls) column.\n") )

  names1<-hold_names
  if("MAF" %in% hold_names) cat( "Interpreting the MAF column as the MAF (minor allele frequency) column.\n" )
  hold_names[hold_names %in% ss.tags$maf] <- "MAF"
  if(length(setdiff(names1,hold_names)) > 0)
    cat( paste("Interpreting the", setdiff(names1, hold_names), "column as the MAF (minor allele frequency) column.\n") )
  
  if ( maf.override ) {
    ##override MAF with reference [NOTE: this assumes the reference has a 'MAF' field]
    ##rename common MAF labels to MAF_Other so the MAF from the reference file is used across traits for conversions
    hold_names[hold_names %in% ss.tags$maf] <- "MAF_Other"
  }

  return( hold_names )

}

