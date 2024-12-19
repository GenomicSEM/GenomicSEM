summaryGLSbands <- function(OBJECT = NULL, 
                            Y = NULL, V_Y = NULL, 
                            PREDICTORS, INTERVALS = 20,
                            CONTROLVARS = NULL, INTERCEPT = T,
                            QUAD = F, BANDS = T, BAND_SIZE = 1, 
                            XLAB = "", YLAB = "",
                            XCOORDS = c(-2,2),
                            YCOORDS = c(0,1)){
  
  
  if (is.vector(PREDICTORS)) {
    PREDICTORS <- matrix(PREDICTORS, ncol = 1)
  }
  
  num_predictors <- nrow(PREDICTORS)
  confintvalpredictors <- PREDICTORS  
  
  if (QUAD) {
    # Add a second column which is the square of the first column
    PREDICTORS <- cbind(PREDICTORS, PREDICTORS^2)
  }
  
  X <- PREDICTORS
  
  # Check if CONTROLVARS is not NULL
  if (!is.null(CONTROLVARS)) {
    # Ensure CONTROLVARS is a matrix
    if (is.vector(CONTROLVARS)) {
      CONTROLVARS <- matrix(CONTROLVARS, ncol = 1)
    }
    
    # Combine PREDICTORS with CONTROLVARS
    X <- cbind(X, CONTROLVARS)
  }
  
  
  colnames(PREDICTORS) <- paste("b", 1:ncol(PREDICTORS), sep = "")
  colnames(X) <- paste("b", 1:ncol(X), sep = "")
  
  ###Creating object for outcome variables for GLS analysis
  ###Creating matrix of predictors, defaulting to including intercept
  if (INTERCEPT) {
    X <- cbind(rep(1, num_predictors), X)
    colnames(X)[1] <- "b0"
  } else {
    X <- X
  }
  
  
  ###Creating a matrix of subsetted V matrix
  if (is.null(Y)){
    Ohm <- matrix(unlist(OBJECT[1]), nrow = num_predictors, ncol = num_predictors)
  }  
  
  if (is.null(OBJECT)){
    Ohm <- V_Y
  }  
  
  ###Creating a vector of genetic correlations
  if (is.null(Y)){
    y   <- unlist(OBJECT[2])
  }  
  
  if (is.null(OBJECT)){
    y <- Y
  }  
  
  ###Estimating betas with GLS equation
  BETAS <- solve(t(X) %*% solve(Ohm) %*% as.matrix(X)) %*% t(X) %*% solve(Ohm) %*% y #beta      
  
  ###Grabbing standard errors of 
  SE <- sqrt(diag(solve(t(X) %*% solve(Ohm) %*% as.matrix(X)))) #SE of beta
  Z <- BETAS/SE
  Pvals <- 2*pnorm(abs(BETAS)/SE,lower.tail=FALSE) #p value 
  
  outputGLS <- cbind(BETAS[,1], Pvals, SE, (BETAS/SE))
  colnames(outputGLS) <- c("betas", "pvals", "SE", "Z")
  
  
  ###Creating sets of predictors to estimate bands along the length of predictors by shifting
  if (QUAD) {
    confintval <- matrix(NA, nrow = num_predictors, ncol = INTERVALS*2)
    for (i in 1:INTERVALS){
      confintval[,2*i-1] <- confintvalpredictors -  
        (i)*(max(range(confintvalpredictors))-min(range(confintvalpredictors)))/INTERVALS
      confintval[,(2*i)] <- (confintvalpredictors -    
                               (i)*(max(range(confintvalpredictors))-min(range(confintvalpredictors)))/INTERVALS)^2
    }
  } else { 
    confintval <- matrix(NA, nrow = num_predictors, ncol = INTERVALS)
    for (i in 1:INTERVALS){
      confintval[,i] <- confintvalpredictors - 
        (i)*(max(range(confintvalpredictors))-min(range(confintvalpredictors)))/INTERVALS
    }
  }
  
  ###### DO I NEED TO CHANGE THIS
  #Adding an intercept column to the object of all shifted predictors
  confintval <- cbind(rep(1,num_predictors),confintval)
  
  #Looping through GLS for all predictors and popping out vector of standard errors
  standerrors <- matrix(NA, nrow = INTERVALS, ncol = 1)
  
  if (QUAD) {
    for (i in 1:INTERVALS){
      XX <- cbind(confintval[,1], confintval[,2*i], confintval[,((2*i)+1)])
      if (exists("CONTROLVARS") && length(CONTROLVARS) != 0) {
        # If controlvariables is a vector, convert it to a matrix (if necessary)
        if (is.vector(CONTROLVARS)) {
          CONTROLVARS <- matrix(CONTROLVARS, ncol = 1)
        }
        # Concatenate CONTROLVARS onto the matrix XX
        XX <- cbind(XX, CONTROLVARS)
      }
      
      ###Creating a matrix of subsetted V matrix
      if (is.null(Y)){
        Ohm <- matrix(unlist(OBJECT[1]), nrow = num_predictors, ncol = num_predictors)
      }  
      
      if (is.null(OBJECT)){
        Ohm <- V_Y
      }  
      
      ###Creating a vector of genetic correlations
      if (is.null(Y)){
        y   <- unlist(OBJECT[2])
      }  
      
      if (is.null(OBJECT)){
        Ohm <- Y
      }  
      
      ###Estimating betas with GLS equation
      BETA <- solve(t(XX) %*% solve(Ohm) %*% XX) %*% t(XX) %*% solve(Ohm) %*% y #beta        
      ###Grabbing standard errors of 
      SE <- sqrt(diag(solve(t(XX) %*% solve(Ohm) %*% XX))) #SE of beta
      Pvals <- 2*pnorm(abs(BETA)/SE,lower.tail=FALSE) #p value 
      standerrors[i,1] <- SE[1]
    }        
    
  } else { 
    for (i in 1:INTERVALS){
      XX <- cbind(confintval[,1], confintval[, i + 1])
      
      if (exists("CONTROLVARS") && length(CONTROLVARS) != 0) {
        # If controlvariables is a vector, convert it to a matrix (if necessary)
        if (is.vector(CONTROLVARS)) {
          CONTROLVARS <- matrix(CONTROLVARS, ncol = 1)
        }
        # Concatenate CONTROLVARS onto the matrix XX
        XX <- cbind(XX, CONTROLVARS)
      }
      ###Creating a matrix of subsetted V matrix
      if (is.null(Y)){
        Ohm <- matrix(unlist(OBJECT[1]), nrow = num_predictors, ncol = num_predictors)
      }  
      
      if (is.null(OBJECT)){
        Ohm <- V_Y
      }  
      
      ###Creating a vector of genetic correlations
      if (is.null(Y)){
        y   <- unlist(OBJECT[2])
      }  
      
      if (is.null(OBJECT)){
        Ohm <- Y
      }  
      
      ###Estimating betas with GLS equation
      BETA <- solve(t(XX) %*% solve(Ohm) %*% XX) %*% t(XX) %*% solve(Ohm) %*% y #beta        
      ###Grabbing standard errors of 
      SE <- sqrt(diag(solve(t(XX) %*% solve(Ohm) %*% XX))) #SE of beta
      Pvals <- 2*pnorm(abs(BETA)/SE,lower.tail=FALSE) #p value 
      standerrors[i,1] <- SE[1]
    }      
  }
  
  evenly_spaced_values <- seq(from = min(PREDICTORS), to = max(PREDICTORS), length.out = INTERVALS)
  
  ###getting error bar values 
  if (QUAD) {
    ##creating standard error bar object for plotting
    PREDICTORSSES <- seq(from = min(PREDICTORS[,1]), to = max(PREDICTORS[,1]), length.out = INTERVALS)
    PREDICTORSSES <- c(PREDICTORSSES)
    eqval <- BETAS[1] + BETAS[2] * PREDICTORSSES + BETAS[3] * PREDICTORSSES^2
    
    plotdatases <- data.frame(PREDICTORSSES, eqval + BAND_SIZE*standerrors, eqval - BAND_SIZE*standerrors)
    names(plotdatases) <- c("Predictors", "SEss", "SEs")
    
    #creating object with rgs for plotting
    plotdata <- data.frame(X[,2], y)
    names(plotdata) <- c("Predictors","rGs")
    
  } else { 
    PREDICTORSSES <- seq(from = min(PREDICTORS[,1]), to = max(PREDICTORS[,1]), length.out = INTERVALS)
    eqval <- BETAS[1] + BETAS[2] * PREDICTORSSES
    plotdatases <- data.frame(PREDICTORSSES, eqval + BAND_SIZE*standerrors, eqval - BAND_SIZE*standerrors)
    names(plotdatases) <- c("Predictors", "SEss", "SEs")
    plotdata <- data.frame(X[,2], y)
    names(plotdata) <- c("Predictors","rGs")
  }
  
  
  if (QUAD) {
    if (BANDS) {
      plot <- ggplot(plotdata, aes(x = Predictors, y = rGs)) +
        #geom_line(aes(x = Predictors, y = (rGs + SEs)), color = "#FFC300") +
        geom_point() +
        stat_function(fun = function(Predictors) BETAS[1] + BETAS[2] * Predictors + BETAS[3] * Predictors^2, 
                      color = "#7A9083", linewidth = 0.8, linetype = 2) +
        geom_point(size = 3, color = "black") +
        geom_point(size = 1.2, color = "#DFD0B7") +
        # scale_x_continuous(breaks = seq(-2, 2, by = 0.5)) +
        coord_cartesian(ylim = YCOORDS,
                        xlim = XCOORDS) +
        theme_classic() +
        labs(x = XLAB, y = YLAB) +
        geom_line(data = plotdatases, aes(x = Predictors, y = SEss), linewidth = 0.4, linetype = 3) +
        geom_line(data = plotdatases, aes(x = Predictors, y = SEs), linewidth = 0.4, linetype = 3) +
        geom_ribbon(data = plotdatases, inherit.aes = FALSE, aes(x = Predictors, ymin = SEss, ymax = SEs),                fill = "#7A9083", alpha = 0.1) 
    } else {
      plot <- ggplot(plotdata, aes(x = Predictors, y = rGs)) +
        #geom_line(aes(x = Predictors, y = (rGs + SEs)), color = "#FFC300") +
        geom_point() +
        stat_function(fun = function(Predictors) BETAS[1] + BETAS[2] * Predictors + BETAS[3] * Predictors^2, 
                      color = "#7A9083", linewidth = 0.8, linetype = 2) +
        geom_point(size = 3, color = "black") +
        geom_point(size = 1.2, color = "#DFD0B7") +
        # scale_x_continuous(breaks = seq(-2, 2, by = 0.5)) +
        coord_cartesian(ylim = YCOORDS,
                        xlim = XCOORDS) +
        theme_classic() +
        labs(x = XLAB, y = YLAB)}
  } else { 
    
    if (BANDS) {
      plot <- ggplot(plotdata, aes(x = Predictors, y = rGs)) +
        #geom_line(aes(x = Predictors, y = (rGs + SEs)), color = "#FFC300") +
        geom_point() +
        stat_function(fun = function(Predictors) BETAS[1] + BETAS[2] * Predictors, 
                      color = "#7A9083", linewidth = 0.8, linetype = 2) +
        geom_point(size = 3, color = "black") +
        geom_point(size = 1.2, color = "#DFD0B7") +
        # scale_x_continuous(breaks = seq(-2, 2, by = 0.5)) +
        coord_cartesian(ylim = YCOORDS,
                        xlim = XCOORDS) +
        theme_classic() +
        labs(x = XLAB, y = YLAB) +
        geom_line(data = plotdatases, aes(x = Predictors, y = SEss), linewidth = 0.4, linetype = 3) +
        geom_line(data = plotdatases, aes(x = Predictors, y = SEs), linewidth = 0.4, linetype = 3) +
        geom_ribbon(data = plotdatases, inherit.aes = FALSE, aes(x = Predictors, ymin = SEss, ymax = SEs),                fill = "#7A9083", alpha = 0.1) 
    } else {
      plot <- ggplot(plotdata, aes(x = Predictors, y = rGs)) +
        #geom_line(aes(x = Predictors, y = (rGs + SEs)), color = "#FFC300") +
        geom_point() +
        stat_function(fun = function(Predictors) BETAS[1] + BETAS[2] * Predictors, 
                      color = "#7A9083", linewidth = 0.8, linetype = 2) +
        geom_point(size = 3, color = "black") +
        geom_point(size = 1.2, color = "#DFD0B7") +
        # scale_x_continuous(breaks = seq(-2, 2, by = 0.5)) +
        coord_cartesian(ylim = YCOORDS,
                        xlim = XCOORDS) +
        theme_classic() +
        labs(x = XLAB, y = YLAB)}
  }
  print(outputGLS)
  plot
}
