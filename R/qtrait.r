#' Run the QTrait function to test trait-specific heterogeneity in genetic correlation models
#'
#' \code{QTrait} tests whether the genetic association between a latent factor (defined by a set of indicator traits) and an external trait can be fully explained by the common factor, or whether specific indicator traits show heterogeneity in this relationship.
#'
#' @param LDSCoutput An object from the \code{ldsc()} function containing multivariate LDSC output.
#' @param indicators A character vector specifying the names of indicator traits that define the latent factor. Must match trait names in \code{LDSCoutput}.
#' @param traits A character vector of external traits for which QTrait statistics and genetic correlations will be estimated. Must match trait names in \code{LDSCoutput}.
#' @param mresid A proportion (default = 0.25) of the RMS genetic correlation used to flag outlier indicator traits.
#' @param mresidthreshold An absolute threshold (default = 0.10) for residual genetic correlation to define meaningful outliers.
#' @param lsrmr Proportion (default = 0.25) used to define meaningful lSRMR based on the RMS genetic correlation.
#' @param lsrmrthreshold Absolute threshold (default = 0.10) to define meaningful lSRMR.
#' @param save.plots Logical. If TRUE, saves scatterplots of trait–indicator associations against loadings. Default is TRUE.
#' @param stdout Logical. If TRUE, plots use standardized output (correlations vs. standardized loadings). If FALSE, uses covariances vs. unstandardized loadings. Default is TRUE.
#'
#' @return A data frame summarizing genetic correlations between traits and the latent factor, QTrait statistics for both common and follow-up models, heterogeneity flags, lSRMR values, and identified outlier indicators.
#'
#' @details
#' \strong{Common Pathway Model Output} includes:
#' \itemize{
#'   \item \code{rGF1Trait_CPM}, \code{SErGF1Trait_CPM}, \code{pvalrGF1Trait_CPM}: Genetic correlation, standard error, and p-value.
#'   \item \code{QTrait_CPM}, \code{df_CPM}, \code{p_value_CPM}, \code{Qsignificant_CPM}: Q statistic, degrees of freedom, and p-value.
#'   \item \code{lSRMR_CPM}, \code{lSRMR_above_threshold_CPM}, \code{heterogeneity_CPM}: Local SRMR and heterogeneity status.
#' }
#'
#' \strong{Follow-Up Model Output} includes analogous statistics when outlier indicators are freed:
#' \itemize{
#'   \item \code{rGF1Trait_FUM}, \code{SErGF1Trait_FUM}, \code{pvalrGF1Trait_FUM}
#'   \item \code{QTrait_FUM}, \code{df_FUM}, \code{p_value_FUM}, \code{Qsignificant_FUM}
#'   \item \code{lSRMR_FUM}, \code{lSRMR_above_threshold_FUM}, \code{heterogeneity_FUM}
#'   \item \code{Unconstrained_paths}: Names of outlier indicator traits.
#' }
#'
#' @examples
#' \dontrun{
#' # Load example LDSC output
#' load("LDSC_G_factor_QTRAIT_tutorial.RData")
#'
#' # Define indicators of genetic g
#' indicators <- c("Matrix", "Memory", "RT", "Symbol_Digit", "TMTB", "Tower", "VNR")
#'
#' # Define external correlate
#' traits <- "EA"
#'
#' # Run QTrait
#' qtrait_out <- QTrait(
#'   LDSCoutput = LDSC_G_factor_QTRAIT_tutorial,
#'   indicators = indicators,
#'   traits = traits,
#'   mresid = 0.25,
#'   mresidthreshold = 0.10,
#'   lsrmr = 0.25,
#'   lsrmrthreshold = 0.10,
#'   save.plots = TRUE,
#'   stdout = TRUE
#' )
#'
#' print(qtrait_out)
#' }
#'
#' @author Javier de la Fuente
#' @references \url{https://rpubs.com/JaFuente/QTrait}
#' @export
QTrait <- function(LDSCoutput,indicators,traits,
                   mresid=.25,mresidthreshold=.10,
                   lsrmr=.25,lsrmrthreshold = .10,
                   save.plots=TRUE,stdout = TRUE){
  suppressWarnings({
  #Load required packages  
  list.of.packages <- c("data.table","GenomicSEM","ggplot2","ggrepel","corrplot","dplyr")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) stop("Missing package(s) ", paste0(new.packages, collapse=" and "))
  lapply(list.of.packages, library,character.only = TRUE)
  # Load LDSC output
  # check if variable is a file path
  if(is.character(LDSCoutput)) {
    if(file.exists(LDSCoutput)) {
      # load serialized LDSC output from the file and 
      # capture the loaded object's name
      object_name <- load(LDSCoutput)
      # assign loaded object by name
      LDSCoutput <- get(object_name)
    } else {
      stop(paste("LDSC file", LDSCoutput, "does not exist."))
    }
  }
  
  # check if LDSCoutput is a list with the expected elements
  if(!is.list(LDSCoutput) | any(!c("V", "S", "I", "N", "m") %in% names(LDSCoutput))) {
    LDSCoutput_name <- deparse(substitute(LDSCoutput))
    stop(paste(LDSCoutput_name, "is not an ldsc() object"))
  }
  
  # check if LDSCoutput has genetic correlation matrices
  if(!all(c("V_Stand", "S_Stand") %in% names(LDSCoutput))) {
    stop(paste("QTrait() requires a genetic correlation and sampling correlation matrix. Rebuild the covariance structure with ldsc(..., stand = TRUE)"))
  }
  
  # Write Common Factor Model syntax 
  CF_model_syntax <- paste(
    paste("F1 =~ ", paste(indicators, collapse = " + ",sep=""), collapse = " "),
    paste(indicators, "~~ rv", 1:length(indicators), "*", indicators, collapse = "\n",sep=""),
    paste("rv", 1:length(indicators), ">0.001", collapse = "\n",sep=""),
    sep = "\n") 
  # Fit Common Factor Model without external correlates
  sink("NUL")
  CFM_fit <- usermodel(covstruc = LDSCoutput, estimation = "DWLS", model = CF_model_syntax,#Run CFM
                         CFIcalc = F, std.lv = T, imp_cov = F)
  sink()
  # Extract unstandardized lambdas from common factor model without external correlates
  unstd_lambdas_CFM = CFM_fit$results[CFM_fit$results$lhs=="F1" & 
                                      CFM_fit$results$op=="=~" 
                                      ,c("rhs","Unstand_Est")]
  # Write Common Factor Model syntax with fixed measurement model
  CF_model_syntax <- paste(
    paste("F1 =~", paste(unstd_lambdas_CFM[,"Unstand_Est"], "*", unstd_lambdas_CFM$rhs, collapse = " + ")),
    paste0(indicators, "~~ rv", 1:length(indicators), "*", indicators, collapse = "\n"),
    paste0("rv", 1:length(indicators), ">0.001", collapse = "\n"),
    #paste(paste0("l", 1:length(unstd_lambdas_CFM[,"Unstand_Est"]), 
    #            " := ", unstd_lambdas_CFM[,"Unstand_Est"], 
    #            "/sqrt(", unstd_lambdas_CFM[,"Unstand_Est"], "^2 + rv", 1:length(unstd_lambdas_CFM[,"Unstand_Est"]), ")"), 
    #      collapse = "\n"),
    sep = "\n"
  )

  #Create matrix to store QTrait statistics 
  Q_mat <- data.frame(
      rGF1Trait_CPM = numeric(),
      SErGF1Trait_CPM = numeric(),
      pvalrGF1Trait_CPM = numeric(),
      rGF1Trait_significat_CPM = numeric(),
      QTrait_CPM = logical(), 
      df_CPM = numeric(),
      p_value_CPM = numeric(),
      Qsignificant_CPM = numeric(),
      lSRMR_CPM = numeric(),
      lSRMR_above_threshold_CPM = numeric(),
      heterogeneity_CPM = character(),
      rGF1Trait_FUM = numeric(),
      SErGF1Trait_FUM = numeric(),
      pvalrGF1Trait_FUM = numeric(),
      rGF1Trait_significat_FUM = numeric(),
      QTrait_FUM = numeric(),  
      df_FUM = numeric(),
      p_value_FUM = numeric(),
      Qsignificant_FUM = numeric(),
      lSRMR_FUM = numeric(),
      lSRMR_above_threshold_FUM = numeric(),
      reduction_lSRMR = numeric(),
      heterogeneity_FUM = character(),
      Unconstrained_paths = character(),
      stringsAsFactors = FALSE  # To prevent strings from being converted to factors
    )

  plot_list <- list()
  rGs_plots <- list()
  betaF1Trait <- list()
  SEF1Trait <- list()
  pvalbetaF1Trait <- list()
  rGF1Trait <- list()
  SErGF1Trait <- list()
  pvalrGF1Trait <- list()

  if(save.plots & !file.exists("Plots")){
    dir.create("Plots")
  }
  
  for (i in traits){
    cat('\n',"----------------------------","\n") 
    cat('\n',"Fitting common pathways model for external trait", i,"\n","\n")
    #Fit Common Pathways model
    l1 <- paste('\n',i,"F"," =~ NA*",i,sep = "")
    l2 <- paste('\n',i," ~~ 0*",i," ",sep = "")
    l3 <- paste('\n',"F1"," ~ ",i,"F"," ",sep = "")
    l4 <- paste('\n',i,"F"," ~~ 1*",i,"F"," ",sep = "")
    full.mod <- paste(CF_model_syntax,l1,l2,l3,l4," \n ", sep = "") #Store full model syntax
    sink("NUL")
    CPM_fit <- usermodel(covstruc = LDSCoutput, estimation = "DWLS", model = full.mod,#Run CPM
                         CFIcalc = TRUE, std.lv = F, imp_cov = T)
    sink()
    CPchi <- CPM_fit$modelfit$chisq #Store chisq of CPM
    CPdf <-  CPM_fit$modelfit$df #Store df of CPM   
   
   # Extract the residual variances (lhs == rhs and op == "~~")
   residuals <- CPM_fit$results[CPM_fit$results$op == "~~" & CPM_fit$results$lhs == CPM_fit$results$rhs, c("lhs", "Unstand_Est")]
   # Create a vector to store the standardized lambdas
   lambdas <- data.frame(rhs = unstd_lambdas_CFM$rhs, Unstand_Est = unstd_lambdas_CFM$Unstand_Est, STD_Lambda = NA)
   # Loop through lambdas and calculate standardized lambdas
   for(r in 1:nrow(lambdas)) {
     trait <- lambdas$rhs[r]  # Get the trait name
     residual_variance <- residuals[residuals$lhs == trait, "Unstand_Est"]  # Get the corresponding residual variance
     # Calculate the standardized lambda
     std_lambda <- lambdas$Unstand_Est[r] / sqrt(lambdas$Unstand_Est[r]^2 + residual_variance)
     lambdas$STD_Lambda[r] <- std_lambda
   }

    #Extract beta coefficient relating the external correlate with the factor
    betaF1Trait[[i]][["CPM"]] <- CPM_fit$results[CPM_fit$results$lhs == "F1" 
                                    & CPM_fit$results$op == "~" 
                                    & CPM_fit$results$rhs == paste(i,"F",sep = ""),
                                    c("Unstand_Est","STD_All")]

    SEF1Trait[[i]][["CPM"]] <- as.numeric(CPM_fit$results[CPM_fit$results$lhs == "F1" 
                                    & CPM_fit$results$op == "~" 
                                    & CPM_fit$results$rhs == paste(i,"F",sep = ""),
                                    "Unstand_SE"])
     

    pvalbetaF1Trait[[i]][["CPM"]] <- ifelse(CPM_fit$results[CPM_fit$results$lhs == "F1" 
                                              & CPM_fit$results$op == "~" 
                                              & CPM_fit$results$rhs == paste(i,"F",sep = ""),
                                              "p_value"]=="< 5e-300",5e-300,as.numeric(CPM_fit$results[CPM_fit$results$lhs == "F1" 
                                              & CPM_fit$results$op == "~" 
                                              & CPM_fit$results$rhs == paste(i,"F",sep = ""),
                                              "p_value"]))
    
    #Fit Independent Pathways model
    cat('\n',"Fitting independent pathways model for external trait", i,"\n","\n")
    l1 <- paste('\n',i,"F"," =~ NA*",i," ",sep = "")
    l2 <- paste('\n',i," ~~ 0*",i," ",sep = "")
    l3 <- paste('\n',paste(indicators, collapse = " + ",sep="")," ~ ",i,"F"," ",sep = "")
    l4 <- paste('\n',i,"F"," ~~ 1*",i,"F"," ",sep = "")
    l5 <- paste('\n',"F1~~0*" ,i,"F"," ",sep = "")
    l6 <- paste('\n',"F1~~rvl*F1"," ",sep = "")
    l7 <- paste('\n',"rvl>0.001"," ",sep = "")
    full.mod <- paste(CF_model_syntax,l1,l2,l3,l4,l5,l6,l7," \n ", sep = "") #Store full model syntax
    sink("NUL")
    IPM_fit <- usermodel(covstruc = LDSCoutput, estimation = "DWLS", model = full.mod,
                         CFIcalc = TRUE, std.lv = F, imp_cov = F)
    sink() 
    IPchi <- IPM_fit$modelfit$chisq
    IPdf <- IPM_fit$modelfit$df
    # Check if the chisq from the IPM is NA and assign 0 if it is
    if (is.na(IPchi)) {
      IPchi <- 0
    }
    # Nested model comparison
    nested_chi <- CPchi - IPchi
    nested_df <- CPdf - IPdf
  
    #LDSC genetic correlations vs factor loadings
    #Store observed genetic correlations
    S_Stand <- LDSCoutput$S_Stand
    rownames(S_Stand) <- colnames(S_Stand)
    S_Stand <- S_Stand[i,indicators]
    #Store model-implied genetic correlations
    MImp_S_Stand <- cov2cor(CPM_fit$resid_cov$`Model Implied Covariance Matrix`)
    MImp_S_Stand <- MImp_S_Stand[i,indicators]
    #Calculate residual genetic correlation matrix
    S_Stand_Resid <- S_Stand - MImp_S_Stand
    #Extract genetic covariances
    S <- LDSCoutput$S
    rownames(S) <- colnames(S)
    S <- S[i,indicators]
    #Store matrix containing observed, model-implied, and residual genetic correlations
    rGs <- as.data.frame(cbind(S_Stand,MImp_S_Stand,S_Stand_Resid,S))    
    #Local SRMR
    # Square, mean, sqrt residual rGs to get lsrmr
    diffsq <- rGs$S_Stand_Resid^2
    lsrmr_cpm <- sqrt(mean(diffsq))
    # Bonferroni p-value threshold
    bfpvalt <- .05/length(traits)
    # Store whether BetaF1Trait is significant based on Bonferroni p-value threshold
    BetaF1Trait_significat <- ifelse(pvalbetaF1Trait[[i]][["CPM"]] < bfpvalt, "*", "NS")
    BetaF1Trait_significat_CPM <- BetaF1Trait_significat
    # Store whether QTrait is significant based on Bonferroni p-value threshold
    Qsignificant <- ifelse(pchisq(nested_chi, nested_df, lower.tail = FALSE) < bfpvalt, "*", "NS")
    # Store whether lSRMR is above threshold
    lSRMR_above_threshold <- ifelse(lsrmr_cpm > lsrmrthreshold & lsrmr_cpm > lsrmr*sqrt(mean(rGs$S_Stand^2)), "Yes", "No")

    #Use three-criteria test of heterogeneity to determine significant heterogeneity
    SigHet <- ifelse(BetaF1Trait_significat == "*" & Qsignificant == "*" & lSRMR_above_threshold == "Yes", "Yes", "No")

    #If there is significant heterogeneity based on the two-criteria
    #test (Qtrait statistic + local SRMR) proceed to identify outliers
    #and run follow-up models
    if(SigHet=="Yes"){
    # 1. Identify outliers with residuals that are more than  half the mean 
    # of the observed rG and surpassing the mresidthreshold
    mrG <- sqrt(mean(rGs$S_Stand^2)) #root mean square genetic correlation
    inside <- ifelse(abs(rGs$S_Stand_Resid) < mresidthreshold |
                     abs(rGs$S_Stand_Resid) < mresid*mrG, "", indicators)  
    #Store outliers
    `%nin%` = Negate(`%in%`)
    if(identical(inside, rep("", length(indicators)))){
      Outliers <- "None"
    } else {
      Outliers <- paste(inside[inside%nin%""],collapse = ",")
    }
   
   #2. Run follow-up models
   # If there are outliers
   if(Outliers!="None"){
        #1. Identify the most extreme outlier
        outlier_fum <- rownames(rGs)[which.max(abs(rGs$S_Stand_Resid))]
        #2. Fit follow-up model for most extreme outlier
        cat('\n',"Fitting follow-up model for external trait", i,"\n","\n")
        cat('\n',"Most egregious outlier:", outlier_fum,"\n","\n")
        l1 <- paste('\n',i,"F"," =~ NA*",i,sep = "")
        l2 <- paste('\n',i," ~~ 0*",i," ",sep = "")
        l3 <- paste('\n',"F1"," ~ ",i,"F"," ",sep = "")
        l4 <- paste('\n',i,"F"," ~~ 1*",i,"F"," ",sep = "")
        l5 <- paste('\n', outlier_fum," ~ ",i,"F"," ",sep = "")#This is the uncontrained direct pathway
        l7 <- paste('\n',"F1~~rvl*F1"," ",sep = "")
        l8 <- paste('\n',"rvl>0.001"," ",sep = "")
        followup.mod <- paste(CF_model_syntax,l1,l2,l3,l4,l5,l7,l8," \n ", sep = "") #Store full model syntax
        sink("NUL")
        FUM_fit <- usermodel(covstruc = LDSCoutput, estimation = "DWLS", model = followup.mod,#Run follow-up model
                             CFIcalc = TRUE, std.lv = F, imp_cov = T)
        sink()
        FUMchi <- FUM_fit$modelfit$chisq #Store chisq of follow-up model
        FUMdf <-  FUM_fit$modelfit$df #Store df of follow-up model
        cat('\n',"----------------------------","\n") 
        # Check if the chisq from the IPM is NA and assign 0 if it is
        if (is.na(FUMchi)) {
          FUMchi <- 0
        }
        # Nested model comparison
        nested_chi_FUM <- FUMchi - IPchi
        nested_df_FUM <- FUMdf - IPdf

        #Extract beta coefficient relating the external correlate with the factor
        betaF1Trait[[i]][[outlier_fum]] <- FUM_fit$results[FUM_fit$results$lhs == "F1" 
                                    & FUM_fit$results$op == "~" 
                                    & FUM_fit$results$rhs == paste(i,"F",sep = ""),
                                    c("Unstand_Est","STD_All")]

        SEF1Trait[[i]][[outlier_fum]] <- as.numeric(FUM_fit$results[FUM_fit$results$lhs == "F1" 
                                    & FUM_fit$results$op == "~" 
                                    & FUM_fit$results$rhs == paste(i,"F",sep = ""),
                                    "Unstand_SE"])                          
    
        pvalbetaF1Trait[[i]][[outlier_fum]] <- as.numeric(FUM_fit$results[FUM_fit$results$lhs == "F1" 
                                              & FUM_fit$results$op == "~" 
                                              & FUM_fit$results$rhs == paste(i,"F",sep = ""),
                                              "p_value"])

        #Compare observed vs model-implied correlations from the follow-up model
        # Store model-implied genetic correlations
        MImp_S_Stand_FUM <- cov2cor(FUM_fit$resid_cov$`Model Implied Covariance Matrix`)[i, indicators]
        # Calculate residual genetic correlation matrix
        S_Stand_Resid_FUM <- S_Stand - MImp_S_Stand_FUM
        # Store observed, model-implied, and residual genetic correlations
        rGs_FUM <- data.frame(S_Stand, MImp_S_Stand_FUM, S_Stand_Resid_FUM,S)
        # Local SRMR of follow-up model
        lsrmr_FUM <- sqrt(mean((rGs_FUM$S_Stand_Resid) ^ 2))        
        # Store whether QTrait is significant based on Bonferroni p-value threshold
        Qsignificant_FUM <- ifelse(pchisq(nested_chi_FUM, nested_df_FUM, lower.tail = FALSE) < bfpvalt, "*", "NS")
        # Store whether lSRMR is significant
        lSRMR_above_threshold_FUM <- ifelse(lsrmr_FUM > lsrmrthreshold & lsrmr_FUM > lsrmr*sqrt(mean(rGs$S_Stand^2)), "Yes", "No")
        # Store whether BetaF1Trait is significant based on Bonferroni p-value threshold
        BetaF1Trait_significat <- ifelse(pvalbetaF1Trait[[i]][[outlier_fum]] < bfpvalt, "*", "NS")
        #Use two-criteria test of heterogeneity to determine significant heterogeneity in follow-up model
        SigHet_FUM <- ifelse(BetaF1Trait_significat == "*" & Qsignificant_FUM == "*" & lSRMR_above_threshold_FUM == "Yes", "Yes", "No")
        
        #Check if there is significant heterogenity in the follow-up model
        k <- 1
        outliers_fum <- outlier_fum
        while(SigHet_FUM=="Yes"&& 
              nested_df_FUM!=0 &&
              !all(abs(rGs_FUM$S_Stand_Resid) < mresid*mrG | abs(rGs_FUM$S_Stand_Resid) < mresidthreshold)){
          #Identify outliers with residuals that are more than  half the mean of the observed rG
          inside <- ifelse(abs(rGs_FUM$S_Stand_Resid) < mresid*mrG | abs(rGs_FUM$S_Stand_Resid) < mresidthreshold, "", indicators)
          ordered_names <- rownames(rGs_FUM)[order(abs(rGs$S_Stand_Resid), decreasing = TRUE)]
          filtered_names <- ordered_names[ordered_names %nin% unlist(strsplit(outliers_fum, ","))]
          #Store outliers
          if(identical(inside, rep("", length(indicators)))){
            outliers_fum <- "None"
            } else {
           outliers_fum <- paste(outlier_fum,filtered_names[1],sep=",")
          }
          # Split the string into individual names
          outliers_fum <- unlist(strsplit(outliers_fum, ","))
          # Remove duplicates
          outliers_fum <- unique(outliers_fum)
          # Concatenate the unique names back into a string
          outliers_fum <- paste(outliers_fum, collapse = ",")     
          
          #Fit follow-up model
          cat('\n',"Fitting follow-up model for external trait", i,"\n","\n")
          cat('\n',"Most egregious outliers:", outliers_fum,"\n","\n")
          l1 <- paste('\n',i,"F"," =~ NA*",i,sep = "")
          l2 <- paste('\n',i," ~~ 0*",i," ",sep = "")
          l3 <- paste('\n',"F1"," ~ ",i,"F"," ",sep = "")
          l4 <- paste('\n',i,"F"," ~~ 1*",i,"F"," ",sep = "")
          l5 <- paste('\n', paste(paste(strsplit(outliers_fum, ",")[[1]], collapse = " + ",sep=""))," ~ ",i,"F"," ",sep = "")#This is the uncontrained direct pathway
          l7 <- paste('\n',"F1~~rvl*F1"," ",sep = "")
          l8 <- paste('\n',"rvl>0.001"," ",sep = "")
          followup.mod <- paste(CF_model_syntax,l1,l2,l3,l4,l5,l7,l8," \n ", sep = "") #Store full model syntax
          sink("NUL")
          FUM_fit <- usermodel(covstruc = LDSCoutput, estimation = "DWLS", model = followup.mod,#Run follow-up model
                               CFIcalc = TRUE, std.lv = F, imp_cov = T)
          sink()
          FUMchi <- FUM_fit$modelfit$chisq #Store chisq of follow-up model
          FUMdf <-  FUM_fit$modelfit$df #Store df of follow-up model
          cat('\n',"----------------------------","\n") 
                               
          if (FUMdf==0){
            nested_chi_FUM <- 0
            nested_df_FUM <- 0
          }

          #Extract beta coefficient relating the external correlate with the factor
          betaF1Trait[[i]][[paste(outliers_fum, collapse = ",")]] <- FUM_fit$results[FUM_fit$results$lhs == "F1" 
                                    & FUM_fit$results$op == "~" 
                                    & FUM_fit$results$rhs == paste(i,"F",sep = ""),
                                    c("Unstand_Est","STD_All")]

          SEF1Trait[[i]][[paste(outliers_fum, collapse = ",")]] <- as.numeric(FUM_fit$results[FUM_fit$results$lhs == "F1" 
                                    & FUM_fit$results$op == "~" 
                                    & FUM_fit$results$rhs == paste(i,"F",sep = ""),
                                    "Unstand_SE"])

          #Extract p value for beta coefficient relating the external correlate with the factor in the follow-up model
          pvalbetaF1Trait[[i]][[paste(outliers_fum, collapse = ",")]] <- as.numeric(FUM_fit$results[FUM_fit$results$lhs == "F1" 
                                              & FUM_fit$results$op == "~" 
                                              & FUM_fit$results$rhs == paste(i,"F",sep = ""),
                                              "p_value"])
          # Nested model comparison
          nested_chi_FUM <- FUMchi - IPchi
          nested_df_FUM <- FUMdf - IPdf
          # Model-implied genetic correlations
          MImp_S_Stand_FUM <- cov2cor(FUM_fit$resid_cov$`Model Implied Covariance Matrix`)[i, indicators]
          # Residual genetic correlation matrix
          S_Stand_Resid_FUM <- S_Stand - MImp_S_Stand_FUM
          # Observed, model-implied, and residual genetic correlations
          rGs_FUM <- data.frame(S_Stand, MImp_S_Stand_FUM, S_Stand_Resid_FUM,S)
          # Local SRMR
          diffsq_FUM <- rGs_FUM$S_Stand_Resid^2
          lsrmr_FUM <- sqrt(mean(diffsq_FUM))   
          # Store whether QTrait is significant based on Bonferroni p-value threshold
          Qsignificant_FUM <- ifelse(pchisq(nested_chi_FUM, nested_df_FUM, lower.tail = FALSE) < bfpvalt, "*", "NS")
          #Store whether lSRMR is significant
          lSRMR_above_threshold_FUM <- ifelse(lsrmr_FUM > lsrmrthreshold & lsrmr_FUM > lsrmr*sqrt(mean(rGs$S_Stand^2)), "Yes", "No")
          
          #Store whether BetaF1Trait is significant based on Bonferroni p-value threshold
          BetaF1Trait_significat <- ifelse(pvalbetaF1Trait[[i]][[paste(outliers_fum, collapse = ",")]] < bfpvalt, "*", "NS")
          SigHet_FUM <- ifelse(BetaF1Trait_significat == "*" & Qsignificant_FUM == "*" & lSRMR_above_threshold_FUM == "Yes", "Yes", "No")
          Unconstrained_paths <- paste(outliers_fum, collapse = ",")
          outlier_fum <- Unconstrained_paths 
          k <- length(strsplit(outliers_fum, ",")[[1]])    
        }                  
        outliers_fum = outlier_fum  
        Unconstrained_paths <- paste(outliers_fum, collapse = ",") 

    pct_reduction_lSRMR_FUM <- paste0(round(((lsrmr_cpm - lsrmr_FUM) /  lsrmr_cpm) * 100,2),"%")

    #Store QTrait results
    Q_mat[i,] <- c(betaF1Trait[[i]][["CPM"]][[1]],
                   SEF1Trait[[i]][["CPM"]][[1]],
                   pvalbetaF1Trait[[i]][["CPM"]],
                   BetaF1Trait_significat_CPM,
                   nested_chi,nested_df,pchisq(nested_chi, nested_df, lower.tail = FALSE),Qsignificant,
                   lsrmr_cpm,lSRMR_above_threshold,SigHet,
                   betaF1Trait[[i]][[length(betaF1Trait[[i]])]][[1]],
                   SEF1Trait [[i]][[length(betaF1Trait[[i]])]][[1]],
                   pvalbetaF1Trait[[i]][[length(betaF1Trait[[i]])]],
                   BetaF1Trait_significat,
                   nested_chi_FUM,nested_df_FUM,pchisq(nested_chi_FUM,nested_df_FUM,lower.tail = F),
                   Qsignificant_FUM,lsrmr_FUM,lSRMR_above_threshold_FUM,pct_reduction_lSRMR_FUM,SigHet_FUM,Unconstrained_paths)

   } else {
    #Store QTrait results
    Q_mat[i,] <- c(betaF1Trait[[i]][["CPM"]][[1]],
                   SEF1Trait[[i]][["CPM"]][[1]],
                   pvalbetaF1Trait[[i]][["CPM"]],
                   BetaF1Trait_significat_CPM,
                   nested_chi,nested_df,pchisq(nested_chi, nested_df, lower.tail = FALSE),Qsignificant,
                   lsrmr_cpm,lSRMR_above_threshold,SigHet,
                   "-","-","-","-",
                   "-","-","-","-",
                   "-","-","-","None")
   }  
} 

if(SigHet=="No"){
#Store QTrait results
    Q_mat[i,] <- c(betaF1Trait[[i]][["CPM"]][[1]],
                   SEF1Trait[[i]][["CPM"]][[1]],
                   pvalbetaF1Trait[[i]][["CPM"]],
                   BetaF1Trait_significat_CPM,
                   nested_chi,nested_df,pchisq(nested_chi, nested_df, lower.tail = FALSE),Qsignificant,
                   lsrmr_cpm,lSRMR_above_threshold,SigHet,
                   "-","-","-","-",
                   "-","-","-","-",
                   "-","-","-","None")
    }
  
    if(SigHet == "No"){
      Outliers  <- "None"
      Unconstrained_paths <- "None"
      outliers_fum <- Unconstrained_paths
    }

    #3. Plots: rGs and lambdas
    #Correlated vectors plot          
    #Extract unstandardized betas from external correlate to individual traits
    betas_ind <- IPM_fit$results[IPM_fit$results$lhs%in%indicators & 
                      IPM_fit$results$op=="~" 
                      ,c("lhs","Unstand_Est")]
    #Extract unstandardized SEs for the betas from external correlate to individual traits
    SE_betas_ind <- IPM_fit$results[IPM_fit$results$lhs%in%indicators & 
                      IPM_fit$results$op=="~" 
                      ,c("lhs","Unstand_SE")]
    #Extract SE for rGs
    SE_rGs <- matrix(0, nrow(LDSCoutput$S), nrow(LDSCoutput$S))
    SE_rGs[lower.tri(SE_rGs, diag = TRUE)] <- sqrt(diag(LDSCoutput$V_Stand))
    rownames(SE_rGs) <- colnames(LDSCoutput$S);colnames(SE_rGs) <- colnames(LDSCoutput$S)
    SE_rGs <- SE_rGs[i,indicators]
    SE_rGs <- as.data.frame(cbind(indicators,SE_rGs))
    
    # Convert row names to a column
    rGs$Factor = rownames(rGs)
    # Merge the three data frames
    # Merge unstd lambdas with rGs
    rGs <- merge(rGs, lambdas, by.x = "Factor", by.y = "rhs", all.x = TRUE)
    colnames(rGs)[c(6,7)] = c("Lambda_unstd","Lambda_std")
    # Merge betas_ind with rGs
    rGs <- merge(rGs, betas_ind, by.x = "Factor", by.y = "lhs", all.x = TRUE)
    colnames(rGs)[8] = "Beta"
    # Merge betas_ind with rGs
    rGs <- merge(rGs, SE_betas_ind, by.x = "Factor", by.y = "lhs", all.x = TRUE)
    colnames(rGs)[9] = "Beta_SE"
    # Merge SE rGs
    rGs <- merge(rGs, SE_rGs, by.x = "Factor", by.y = "indicators", all.x = TRUE)
    colnames(rGs)[9] = "Beta_SE"
    # Reorder the rows based on the original row names of rGs
    rownames(rGs) <- rGs$Factor

    #Unstandardized weights based on 1/SE^2
    rGs$w_unstd <- 1/(as.numeric(rGs$Beta_SE)^2)
    #Standardized weights
    rGs$w_std <- 1/(as.numeric(rGs$SE_rGs)^2)

#Scatterplot
#Create data.frame to store betas across models
betas <- data.frame(model = character(), slope = numeric(), std = character(), stringsAsFactors = FALSE)

  # Loop through each model
  for (model in names(betaF1Trait[[i]])) {
    # Extract the data for each model
    model_data <- betaF1Trait[[i]][[model]]
    # Replace "CPM" with "Common pathway model" and add "Unconstrained paths: " for follow-up models
    if (model == "CPM") {
      model_name <- "Common pathway model"
    } else {
      model_name <- paste("Unconstrained paths:", model)
    }   
    # Add Unstand_Est (Unstandardized) and STD_All (Standardized)
    betas <- rbind(betas,
                   data.frame(model = model_name, slope = model_data$Unstand_Est[1], std = "Unstandardized"),
                   data.frame(model = model_name, slope = model_data$STD_All[1], std = "Standardized"))
  }


# Get unique models and generate a red color palette
unique_models <- unique(betas$model)
color_palette <- colorRampPalette(c("black", "#BEBEBE"))(length(unique_models))

# Create a named vector mapping models to colors
color_mapping <- setNames(color_palette, unique_models)
colnames(rGs) <- c("Indicator","Observed","Model implied","Residual",
                  "Covariance","Lambda_unstd","Lambda_std","Beta","Beta SE","SE_rG",
                  "Weight_unstd","Weight_std")

# Highlight outliers
rGs$highlight <- ifelse(rGs$Indicator %in% strsplit(outliers_fum, ",")[[1]], "yes", "no")

# Create the plot with the ablines
if (stdout == FALSE) {
  # Determine the range of Lambda and Beta values, handling potential NA values
  x_range <- range(rGs$Lambda_unstd, na.rm = TRUE)
  y_range <- range(rGs$Beta, na.rm = TRUE)

  # Add a buffer (10% of the range)
  buffer <- 0.10  
  x_buffer <- buffer * diff(x_range)
  y_buffer <- buffer * diff(y_range)

  # Define x and y axis limits with buffer
  x_limits <- c(x_range[1] - x_buffer, x_range[2] + x_buffer)
  y_limits <- c(y_range[1] - y_buffer, y_range[2] + y_buffer)

} else {
  # Default axis limits when stdout is TRUE
  x_limits <- c(0, 1)
  
  # Select the appropriate dataset and column names for standardized values
  betas_selected <- betas[betas$std == "Unstandardized",]
  x_var <- "Lambda_std"
  y_var <- "Observed"
  x_label <- "Lambda (Standardized)"
  y_label <- "Genetic correlation (rG)"
  weight <- "Weight_std"

  # Set x-axis breaks
  x_breaks <- seq(0, 1, 0.2)

  # Check the sign of rGs$Observed for dynamic y-axis limits and breaks
  if (all(rGs$Observed < 0)) {
    y_breaks <- seq(-1, 0, 0.2)
    y_limits <- c(-1, 0)
  } else if (all(rGs$Observed > 0)) {
    y_breaks <- seq(0, 1.2, 0.2)
    y_limits <- c(0, 1)
  } else {
    y_breaks <- seq(-1, 1, 0.2)
    y_limits <- c(-1, 1)
  }
}

# Select the appropriate dataset for unstandardized values
if (stdout == FALSE) {
  betas_selected <- betas[betas$std == "Unstandardized",]
  x_var <- "Lambda_unstd"
  y_var <- "Beta"
  x_label <- "Lambda (Unstandardized)"
  y_label <- "Beta (Unstandardized)"
  weight <- "Weight_unstd"
  
  # Dynamically set x and y axis limits based on range
  x_limits <- range(rGs[[x_var]], na.rm = TRUE)
  y_limits <- range(rGs[[y_var]], na.rm = TRUE)
  
  # Define breaks dynamically
  x_breaks <- pretty(x_limits, n = 5)  # Creates about 5 evenly spaced breaks
  y_breaks <- pretty(y_limits, n = 5)
}


betas_selected$model_complexity <- sapply(betas_selected$model, function(x) {
  if (grepl("Unconstrained paths", x)) {
    # Remove redundant "Unconstrained paths:" and split by commas to count paths
    unconstrained_paths <- sub("Unconstrained paths: ", "", x)  # Remove first "Unconstrained paths:"
    unconstrained_paths <- strsplit(unconstrained_paths, ",")[[1]]
    return(length(unconstrained_paths))  # Count number of unconstrained paths
  } else {
    return(0)  # For "Common pathway model"
  }
})

# Reorder the rows by 'model_complexity'
betas_selected <- betas_selected[order(betas_selected$model_complexity), ]

# Start with the basic plot
my_graph <- ggplot(rGs, aes(x = .data[[x_var]], y = .data[[y_var]], size = .data[[weight]], label = Indicator)) +
  scale_fill_manual(values = c("yes" = "red", "no" = "black"), guide = "none") +  
  labs(color = "Model beta") +  
  scale_color_manual(values = color_palette, 
                     breaks = unique(betas_selected$model), 
                     labels = paste(unique(betas_selected$model), round(unique(betas_selected$slope), 3))) + 
  scale_size(guide = "none") +
  scale_x_continuous(breaks = x_breaks, limits = x_limits) +  
  scale_y_continuous(breaks = y_breaks, limits = y_limits) +  
  theme_minimal(base_size = 10, base_line_size = 10/22, base_rect_size = 10/22) +
  labs(title = " ",
       x = x_label,
       y = y_label) + 
  theme_classic(base_family = "serif") +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),     
    axis.title = element_text(size = 14, face = "bold"),   
    axis.ticks = element_line(size = 0.5)    
  )

# Add geom based on stdout value
if (stdout == TRUE) {
  my_graph <- my_graph + 
    geom_segment(data = betas_selected, 
                 aes(x = 0, y = 0, xend = 1, yend = slope * 1, color = model), 
                 inherit.aes = FALSE, linetype = "solid", size = 1) +
    geom_point(aes(fill = highlight), shape = 21) + 
    geom_text(vjust = -1, size = 5)  # Add text after points
  
} else {
  my_graph <- my_graph + 
    geom_abline(data = betas_selected, 
                aes(intercept = 0, slope = slope, color = model), 
                linetype = "solid", size = 1) + 
    geom_point(aes(fill = highlight), shape = 21) + 
    geom_text(vjust = -1, size = 5)  # Add text after points
}
 
# Scatterplot
plot_list[[i]] <- my_graph

    # Convert data to long format for ggplot
     data_long <- tidyr::pivot_longer(rGs[,c(2:4,1)], -Indicator, names_to = "Variable", values_to = "Value")
      # Create the heatmap with custom labels and colors
      rGs_plots[[i]] <- ggplot(data_long, aes(x = factor(Variable, levels = c("Observed", "Model implied", "Residual")),
                                              y = Indicator, fill = Value, label = round(Value, 2))) +
        geom_tile() +
        geom_text(color = "black", size = 5) +
        scale_fill_gradient2(low = "#D73027", mid = "white", high = "#4575B4", midpoint = 0,
                             guide = "legend", limits = c(-1, 1)) +
        theme_minimal(base_family = "serif") + 
        labs(title = paste("Observed, model-implied, and residual genetic correlations across indicators and", i),
             x = "Genetic correlations", y = "Indicator") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"), # Adjust size for x-axis text
              axis.text.y = element_text(size = 12, face = "bold"),  # Adjust size for y-axis text
              axis.title = element_text(size = 14, face = "bold"),  # Adjust size for axis labels
              text = element_text(size = 12),  # Adjust size for all other text
              legend.text = element_text(size = 10),  # Adjust size for legend text
              legend.title = element_text(size = 14),
              plot.title = element_text(size = 14, face = "bold")) +  # Adjust size for legend title
        guides(fill = guide_colorbar(title = "Correlation")) +
        scale_x_discrete(labels = c("Observed", "Model-implied", "Residual"))
      
     vectors_plot <- suppressWarnings(ggpubr::ggarrange(rGs_plots[[i]],plot_list[[i]],ncol = 1,nrow = 2))
     
     if(save.plots){
     # Save as a PDF file
     ggsave(paste("QTrait_Statistics",paste(indicators,collapse = "_"),"_",i,"",".pdf",sep = ""), vectors_plot, 
              path = "Plots",width = 12, height = 13, units = "in")
     }  
    }  
  
  cat("----------------------------","\n") 

  
  # List of non-numeric columns to exclude from conversion
  non_numeric_columns <- c("rGF1Trait_significat_CPM", "Qsignificant_CPM", "lSRMR_above_threshold_CPM","heterogeneity_CPM",
                         "rGF1Trait_significat_FUM", "Qsignificant_FUM", "heterogeneity_FUM", "lSRMR_above_threshold_FUM","reduction_lSRMR",
                         "Unconstrained_paths")
  Q_mat <- Q_mat %>%
  mutate(across(-all_of(non_numeric_columns), ~ as.numeric(.)))
  Q_mat <<- Q_mat
  return(Q_mat)
  cat('\n',"Bonferonni p-value threshold based on the number of tests:",bfpvalt,"\n")
  cat('\n',"----------------------------","\n") 
  })
}  
