install.packages(c("testit", "memoise", "pbapply", "R.cache", "glmmTMB", 
                   "truncnorm", "insight", "performance", "AICcmodavg", 
                   "ggplot2", "gridExtra"))
library(testit)
library(memoise)
suppressPackageStartupMessages(library(R.cache))
library(pbapply)
library(glmmTMB)
library(truncnorm)
library(insight)
library(performance)
library(AICcmodavg)
library(ggplot2)
library(gridExtra)


####################### NEGATIVE BINOMIAL ######################################

###############  DATA GENERATION  ###############

.getNBReads <- function(reps, alpha, beta, sigma2, phi, FUN ,...) {
  
  testit::assert("reps must be a vector", is.vector(reps))
  testit::assert("sigma2 must be non-negative", sigma2 >= 0)
  testit::assert("phi must be positive", phi > 0)
  
  cov_ <<- R.cache::evalWithMemoization({do.call(FUN, list(...))}, key =list(reps, FUN, ...))
  
  num <- length(cov_)
  num.strains <- length(reps)
  intercept <- alpha
  coveffect <- beta*cov_
  straineffect <- rnorm(num.strains, sd = sqrt(sigma2))
  
  means <- exp(intercept + rep.int(straineffect, reps) + coveffect)
  
  NBcounts <- lapply(1:num, function(x){
    counts.x <- MASS::rnegbin(n = 1,
                              mu = means[x],
                              theta = 1/phi)
    return(counts.x)
  })
  NBcounts <- matrix(do.call(c, NBcounts), nrow = 1)
  
  sample.names <- lapply(1:num.strains, function(x){
    sample.x <- paste0("S", x, "_", 1:(reps[x]))
    return(sample.x)
  })
  sample.names <- as.vector(do.call(c, sample.names))
  colnames(NBcounts) <- sample.names
  return(NBcounts)
  
}

getReadMatrix.NB <- function(reps, alphas, betas, sigma2s, phis, FUN, ...){
  num.probes <- length(alphas)
  CountMatrix <- lapply(1:num.probes, function(x){
    probe.x <- .getNBReads(reps, alphas[x], betas[x], sigma2s[x], phis[x], FUN, ...)
  })
  CountMatrix <- do.call(rbind, CountMatrix)
  rownames(CountMatrix) <- paste0("Gene ", 1:num.probes)
  return(CountMatrix)
}

###############  MODEL FITTING  ###############


fit.NB <- function(CountMatrix, Strains, test = FALSE){
  
  if(is.null(dim(CountMatrix))){
    print('Fitting a single feature.')
    CountMatrix <- matrix(CountMatrix, nrow = 1)
  }
  
  GeneIDs <- rownames(CountMatrix)
  paras <- t(pbapply::pbsapply(1:nrow(CountMatrix), function(x){
    CountVector <- CountMatrix[x, ]
    GeneID <- GeneIDs[x]
    dat_sub <- data.frame(expr = as.numeric(CountVector), covariate= cov_, strain = Strains)
    model_sub <- try({glmmTMB::glmmTMB(formula = expr ~ 1 + covariate + (1 | strain),
                                       data = dat_sub, family='nbinom2', REML = TRUE)})
    
    
    # print(AICcmodavg::checkConv(model_sub))
    # print(!performance::check_singularity(model_sub, tolerance = 1e-05))
    # print(class(model_sub) != "try-error")
    if ((!performance::check_singularity(model_sub, tolerance = 1e-05)) && (class(model_sub) != "try-error") && (AICcmodavg::checkConv(model_sub)$converged)){
      fixed <- fixef(model_sub)$cond
      as <- fixed[1]
      bs <- fixed[2]
      sigma_a2 <- get_variance_random(model_sub, tolerance = 1e-05)
      phi <- sigma(model_sub)
      para_sub <- c(as, bs, sigma_a2, phi)
    }
    else{
      print(paste("Fitting problem for feature", x,"returning NA"))
      para_sub <- rep(NA, 4)
    }
    
    if (test){
      model_sub_red <- try({glmmTMB::glmmTMB(formula = expr ~ 1 + (1 | strain),
                                             data = dat_sub, family='nbinom', link='log')}, silent = T)
      if (class(model_sub) != "try-error" &
          class(model_sub_red) != "try-error"){
        pval <- anova(model_sub_red, model_sub)$"Pr(>Chi)"[2]
      }else{
        print(paste("Cannot do test for feature", x,"fitting problem."))
        pval <- NA
      }
      para_sub <- c(para_sub, pval)
    }
    
    return(para_sub)
    
  }))
  paras1 <- matrix(paras[ , 1:4], ncol = 4)
  rownames(paras1) <- GeneIDs
  colnames(paras1) <- c("Intercept", "covariate", "strain", "phi_g")
  if (test){
    pvals <- matrix(paras[,5])
    rownames(pvals) <- GeneIDs
    colnames(pvals) <- "pvalue"
    return(list(paras = paras1, pvals = pvals))
  }else{
    return(list(paras = paras1, pvals = NULL))
  }
  
}

################ COMPUTE VPC ###############

compute1NBVPC <- function(alpha_g, beta, sigma2_g, phi_g){
  if(is.na(sigma2_g)){
    vpc <- NA
  }else{
    if(sigma2_g < 0){
      stop("Random effect variance needs to be non-negative.")
    }
    if(phi_g <= 0){
      stop("Invalid dispersion value.")
    }
    denom <- (exp(sigma2_g) - 1) + (phi_g * exp(sigma2_g)) + (exp(-(alpha_g + beta) - sigma2_g / 2))
    #denom <- exp(sigma2_g) - 1 + exp(sigma2_g)*phi_g + exp(-alpha_g - sigma2_g / 2)
    if(denom == 0){
      vpc <- NA
    }else{
      vpc <- (exp(sigma2_g) - 1) / denom
    }
  }
  return(vpc)
}

computeVPC.NB <- function(para){
  
  if(is.null(dim(para))){
    vpcs <- compute1NBVPC(unlist(para[1]), unlist(para[2]), unlist(para[3]), unlist(para[4]))
  }else{
    vpcs <- apply(para, 1, function(x){
      vpc <- compute1NBVPC(unlist(x[1]), unlist(x[2]), unlist(x[3]), unlist(x[4]))
      return(vpc)
    })
  }
  
  vpcs <- matrix(vpcs, ncol = 1)
  rownames(vpcs) <- rownames(para)
  colnames(vpcs) <- 'NB-fit'
  
  return(vpcs)
}



sim <- function(n){
  vec.num <- c(3, 5, 2, 3, 4, 2)
  as <- rnorm(n, 5)#c(-1, 1, 2, 5, 10)
  sig2s <- rtruncnorm(n, 10)#c(10, 0.2, 0.1, 0.03, 0.01)
  ps <- runif(n, 1.5, 2.5)
  phis <- rtruncnorm(n, 10)#c(1.5, 1, 0.5, 0.1, 0.1)
  betas <- rnorm(n)
  true_param = matrix(c(as, betas, sig2s, phis), ncol=4)
  strains <- factor(c(rep("S1",3),rep("S2",5),rep("S3",2),rep("S4",3),rep("S5",4), rep("S6",2)))
  model_Matrix <- getReadMatrix.NB(vec.num, as, betas, sig2s, phis, "rbinom", n=sum(vec.num), size=1, prob=0.3)
  fit <- fit.NB(model_Matrix, strains, FALSE)
  #print(fit$paras)
  #print(as)
  #print(cor(as,fit$paras[,1], use = "complete.obs"))
  #print(cor(as,fit$paras[,1], use = "complete.obs"))
  #print(as)
  #print(fit$paras[,1])
  #fit_test <- fit.NB(model_Matrix, strains, TRUE)
  true_h2 <- computeVPC.NB(true_param)
  est_h2 <- computeVPC.NB(fit$paras)
  print(paste("Intercept: ", cor(as,fit$paras[,1], use = "complete.obs")))
  print(paste("covariate: ", cor(betas,fit$paras[,2], use = "complete.obs")))
  print(paste("strain: ", cor(sig2s,fit$paras[,3], use = "complete.obs")))
  print(paste("phi: ", cor(phis,fit$paras[,4], use = "complete.obs")))
  print(paste("h2: ", cor(true_h2,est_h2, use = "complete.obs")))
  dat <- data.frame(true_h2, est_h2, as, fit$paras[,1], betas, fit$paras[,2],
                    sig2s, fit$paras[,3], phis, fit$paras[,4])
  colnames(dat) <- c("true_h2", "est_h2", "true_as", "est_as", "true_bs", "est_bs",
                     "true_sigs", "est_sigs", "true_phis", "est_phis")
  g1 = ggplot(dat, aes(x = true_as, y = est_as)) + geom_point()+ggtitle(paste("Intercept: ", round(cor(as,fit$paras[,1], use = "complete.obs"),4)))
  g2 = ggplot(dat, aes(x = true_bs, y = est_bs)) + geom_point()+ggtitle(paste("covariate: ", round(cor(betas,fit$paras[,2], use = "complete.obs"),4)))
  g3 = ggplot(dat, aes(x = true_sigs, y = est_sigs)) + geom_point()+ggtitle(paste("Random Effect Variance: ", round(cor(sig2s,fit$paras[,3], use = "complete.obs"),4)))
  g4 = ggplot(dat, aes(x = true_phis, y = est_phis)) + geom_point()+ggtitle(paste("Dispersion parameter: ", round(cor(phis,fit$paras[,4], use = "complete.obs"),4)))
  g5 = ggplot(dat, aes(x = true_h2, y = est_h2)) + geom_point()+ggtitle(paste("h2: ", round(cor(true_h2,est_h2, use = "complete.obs"),4)))
  grid.arrange(arrangeGrob(g1, g2, ncol=2), arrangeGrob(g3, g4, ncol=2), g5, nrow=3)
  
}


sim(1000)


