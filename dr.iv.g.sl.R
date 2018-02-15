### Doubly-Robust g-estimator ###
### with data-adaptive fits for the instrument propensity score and the expsoure model
### 1st version (Feb 2018) uses parametric models for treatment effect curve and treatment-free outcome model
### assumes continuous outcome and binary exposure and instrument
### based on Vansteelandt Didilez 2015 eprint arXiv:1510.01770
### Programmed by K DiazOrdaz
### if using parametric models for nuisance parameters, needs to bootstrap
## require(SuperLearner)
## 
dr.iv.g.sl <- function(data,
                      eff.mod = T, ## change to F if no effect modification needs to be estimated
                      oneway_adh = F, ## only important when not using SL and for randomised trials where those randomised to control never have access to active
                      Exposure = c("A"), #name of endogeneous exposure
                      eff.mod.cov =c("V"), #name of effect modifier (currently needs to be only one variable)
                      instrument = c("Z"), #name of instrument (only one instrument)
                      my_form  = paste("Y ~  V + W1 + A  ", collapse=""), #specify the variables to be included in the outcome model, if parametric specify non-linear terms
                      ma_form  = paste("A ~ Z + V + W1 ", collapse=""), #specify the variables to be included in the exposure model, if parametric specify non-linear terms
                      mz_form = paste(" Z ~  V + W1 ", collapse=""), #instrument propensity score
                      sl= T, #T for SL, F for parametric. If parametric, only logit models currently supported
                      library = my.sl.library, #specify the machine learning algorithms to be used
                      method = "method.NNLS" #argument to be passed to the Super Learner
                      )
{
  n <- dim(data)[1]
  Zcol <- (1:ncol(data))[colnames(data) %in% instrument]
  Z<- data[,Zcol]
  Acol <- (1:ncol(data))[colnames(data) %in% Exposure]
  A<- data[,Acol]
  exog.vars<-setdiff(attr(terms(as.formula(my_form)), "term.labels"),  "A")
  Exo.covcols <- (1:ncol(data))[colnames(data) %in% exog.vars]
  Exo.matrix <- data.frame(data[, Exo.covcols])
  
  Z_covnames <- attr(terms(as.formula(mz_form)), "term.labels")
  A_covnames <- attr(terms(as.formula(ma_form)), "term.labels")
  Y_covnames <- attr(terms(as.formula(my_form)), "term.labels")
  
  Zcovcols <- (1:ncol(data))[colnames(data) %in% Z_covnames]
  Acovcols <- (1:ncol(data))[colnames(data) %in% A_covnames]
  Ycovcols <- (1:ncol(data))[colnames(data) %in% Y_covnames]
  
  Z.matrix <- data.frame(data[,Zcol])
  Z.cov.matrix <- data.frame(data[,Zcovcols])
  A.cov.matrix <- data.frame(data[,c(Acovcols)])
  Y.cov.matrix <- data.frame(data[,c(Ycovcols)])
  
  if (sl){
    mz <- SuperLearner(data$Z, Z.cov.matrix, newX = Z.cov.matrix, SL.library = library,verbose = FALSE,
                       family =binomial(),cvControl = list(shuffle = F, V=10),method=method)
    Zhat <- mz$SL.predict
    ma <- SuperLearner(data$A, A.cov.matrix, newX = A.cov.matrix, SL.library = library, verbose = FALSE,
                       family =binomial(),cvControl = list(shuffle = F, V=10),method=method)
    Ahat <- ma$SL.predict
    ####
    newdata1 <- A.cov.matrix
    newdata2 <- A.cov.matrix
    newdata1$Z <- rep(0,n)
    newdata2$Z <- rep(1,n)
    m1 <- SuperLearner(data$A, A.cov.matrix, newX = newdata1, SL.library = library, verbose = FALSE,
                       family =binomial(),cvControl = list(shuffle = F, V=10),method=method)
    pi_preds2 <- m1$SL.predict
    m2 <- SuperLearner(data$A, A.cov.matrix, newX = newdata2, SL.library = library, verbose = FALSE,
                       family =binomial(),cvControl = list(shuffle = F, V=10),method=method)
    pi_preds3 <- m2$SL.predict
    A_Z0W <-pi_preds2
    A_Z1W <-pi_preds3
    EAhat_gVW  <- Zhat*A_Z1W + (1-Zhat)*A_Z0W
    Ah_m_EAhgVW <- Ahat - EAhat_gVW
    } else {
      d0_glm =glm(mz_form, data=data, family="binomial")
      Zhat <- predict(d0_glm, newdata=data, type="response")
    
      if (oneway_adh){
      Acov<-setdiff(attr(terms(as.formula(ma_form)), "term.labels"),  instrument)
      ma_form <- paste(Exposure,"~", paste(Acov,collapse=" + "), collapse="")
      pi_glm =(glm(ma_form , data=subset(data,data$Z==1), family="binomial"))
      Ahat <- ifelse(data$Z==1,predict(pi_glm, newdata=data, type="response"),0)
      EAhat_gVW  <- predict(pi_glm, newdata=data, type="response") * Zhat
      Ah_m_EAhgVW = Ahat - EAhat_gVW
      } else {
        newdata1 <- A.cov.matrix
        newdata2 <- A.cov.matrix
        newdata1$Z <- rep(0,n)
        newdata2$Z <- rep(1,n)
        
        pi_glm =(glm(ma_form, data=data, family="binomial"))
        Ahat <- predict(pi_glm, newdata=data, type="response")
        Ahat_Z0  <- predict(pi_glm, newdata=newdata1, type="response") 
        Ahat_Z1  <- predict(pi_glm, newdata=newdata2, type="response") 
        EAhat_gVW  <- Zhat*Ahat_Z1 + (1-Zhat)*Ahat_Z0
        Ah_m_EAhgVW = Ahat - EAhat_gVW
      }
    }
    
  Y<-as.matrix(data$Y)
  ones<-as.matrix(rep(1,dim(data)[1]))
  
  if (eff.mod){
    EMcol <- (1:ncol(data))[colnames(data) %in% eff.mod.cov]
    AV <- data$AV <- data[,Acol]*data[,EMcol]
    ZV <- data[,Zcol]*data[,EMcol]
    V <- data[,EMcol]
    instruments<-c(exog.vars,instrument, "ZV" )
    AhV_m_EAhgVW_V = Ah_m_EAhgVW * V
    Exo.covcols <- (1:ncol(data))[colnames(data) %in% exog.vars]
    Exo.matrix <- data.frame(data[, Exo.covcols])
    X<-as.matrix(cbind(ones,Exo.matrix,A, AV), ncol=dim(Y.cov.matrix)[2]+1+length(EMcol)) 
    Xtilde<-as.matrix(cbind(ones, Exo.matrix, Ah_m_EAhgVW, AhV_m_EAhgVW_V), ncol=dim(Y.cov.matrix)[2]+1+length(EMcol))  
    psihat.P5<-solve(t(Xtilde) %*% X)%*%t(Xtilde) %*% Y
  }else{
    instruments<-c(exog.vars,instrument)
    X<-as.matrix(cbind(ones,Exo.matrix,A), ncol=dim(Y.cov.matrix)[2]+1 )
    Xtilde<-as.matrix(cbind(ones, Exo.matrix, Ah_m_EAhgVW), ncol=dim(Y.cov.matrix)[2]+1 )
    psihat.P5<-solve(t(Xtilde) %*% X)%*%t(Xtilde) %*% Y
  }
  
  psihat.0 <- psihat.P5[row.names(psihat.P5) %in% Exposure]
  psihat.1 <- psihat.P5[row.names(psihat.P5) %in% "AV"]
  psi.hat.dr<-c( psihat.0, psihat.1)
  
  theta.vars <- data.frame(cbind(ones, data[, colnames(Exo.matrix)]))
  m.vars.col<-setdiff(colnames(X), colnames(theta.vars))
  m.vars <- data.frame(data[, m.vars.col])
  theta.hat <- as.matrix(theta.vars)  %*% psihat.P5[row.names(psihat.P5) %in%   colnames(theta.vars)]
  m.hat <- as.matrix(m.vars)  %*% psihat.P5[row.names(psihat.P5) %in% m.vars.col ]
  
  sigma2 <- var(Y-m.hat)
  
  #    sigma_2 = p$Pr_Z1_W*Pi_minus_exp[,3]^2 + (1-p$Pr_Z1_W)*Pi_minus_exp[,2]^2

  kappa <- 1/(as.numeric(sigma2))*(Ah_m_EAhgVW)
  
  
  f0 <- as.matrix(kappa,nrow=n)  
  f1 <- as.matrix(kappa * V ,nrow=n) 
 
  f<-lapply(1:n, function(x){f=matrix(cbind(f0[x,],f1[x,]),ncol=1);return(f)})
  AV.T <- lapply(1:n, function(x){av=matrix(m.vars[x,],ncol=2);return(av)})
  
  Mn <- lapply(1:n, function(x){b=matrix(as.numeric(f[[x]]),ncol=1) %*% matrix(as.numeric(AV.T[[x]]),ncol=2);return(b)})
  
  
  M0.hat<- ginv(Reduce('+', Mn)/n)
  
  
  
  
  
  IF<-lapply(1:n, function(x){IF=t(((M0.hat %*% matrix(as.numeric(f[[x]]),ncol=1)) * (Y[x]-theta.hat[x]) )-  psi.hat.dr)
  return(IF)})
  
  IC<-matrix(unlist(IF),ncol=2, byrow=T)
  
  IC_stats <- function(IC) {
    IC_mean= apply(IC, 2, mean) #mean by columns
    IC_var = apply(IC, 2, var)	
    IC_sd = apply(IC, 2, sd)
    return(list(IC_mean=IC_mean, IC_var=IC_var, IC_sd=IC_sd))
  }

  SE <-  sqrt(IC_stats(IC)$IC_var/n)  
  
  #names(psi.hat.dr)<-names(SE)<-c("psi0", "psi1")
    
  return(list(psi=psi.hat.dr, SE=SE))

}
