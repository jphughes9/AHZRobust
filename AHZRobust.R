library(lme4)
library(Matrix)
library(numDeriv)
library("sandwich")
library("clubSandwich")
library("MASS")
library("expm")

##########################
# Compute AHZ degrees of freedom for an lmer/glmer object using robust 
# variance structure
AHZ.glmerMod = function(obj,L,cluster,type="classic",kadjust=FALSE,Fadjust=FALSE){
  #
  ######################
  # Helper functions (from clubSandwich)
  ######################
  get_outer_group <- function(obj) {
    group_n <- lme4::getME(obj, "l_i")
    group_facs <- lme4::getME(obj, "flist")
    group_facs[[which.min(group_n)]]
  }
  check_nested <- function(inner_grp, outer_grp) {
    n_outer <- tapply(outer_grp, inner_grp, function(x) length(unique(x)))
    all(n_outer == 1)
  }
  is_nested_lmerMod <- function(obj, cluster = get_outer_group(obj)) {
    group_facs <- lme4::getME(obj, "flist")
    nested <- vapply(group_facs, check_nested, outer_grp=cluster, FUN.VALUE=TRUE)
    all(nested)
  }
  #################
  # Other helper functions
  #################
  mtx_DA <- function(D,A) {
    matrix(rep(Matrix::diag(D),ncol(A))*as.numeric(A), ncol=ncol(A))
  }
  #
  mtx_AD <- function(A,D) {
    matrix(rep(Matrix::diag(D), each=nrow(A))*as.numeric(A), ncol=ncol(A))
  }
  #
  sqrtm_inv2 = function(A){
    solve(sqrtm(A))
  }
  sqrtm_inv1 = function(A){
    R = chol(A)
    Rinv = solve(R)
    U = polar(Rinv)
    Rinv%*%t(U)
  }
  sqrtm_inv = function(A){
    tryCatch(sqrtm_inv1(A),error=function(e){sqrtm_inv2(A)})
  }
  #
  decode_type <- function(type){
    ## decode type options
    type1 = type
    r = exact = DF = d = D = NA
    ropt = substr(type,1,2)
    if (ropt=="FG") {    
      type1="FG"
      if (nchar(type)==2) {
        r = 0.75
      } else {
        r = as.numeric(substr(type,regexpr("\\(",type)[[1]]+1,regexpr("\\)",type)[[1]]-1))
      }
    } 
    if (ropt=="KC") {    
      type1="KC"
      if (nchar(type)==2) {
        exact = FALSE
      } else {
        kcstr = tolower(substr(type,regexpr("\\(",type)[[1]]+1,regexpr("\\)",type)[[1]]-1))
        if (!(kcstr=="exact" | kcstr=="e")) stop("Invalid option for KC")
        exact = TRUE
      }
    } 
    if (ropt=="MB" & substr(type,3,3)=="N") {
      type1="MBN"
      # defaults
      DF = TRUE
      d = 2
      r = 1
      D = NA
      text = substr(type,regexpr("\\(",type)[[1]]+1,regexpr("\\)",type)[[1]]-1)
      if (text!=""){
        mbnargs = lapply(strsplit(text,","),strsplit,"=")
        numargs = length(mbnargs[[1]])
        for (i in 1:numargs){
          if (!(mbnargs[[1]][[i]][1] %in% c("DF","d","D","r"))) {
            stop("Allowable arguments for MBN are 'DF','d','D','r'")
          }
          if (mbnargs[[1]][[i]][1]=="DF") {
            assign(mbnargs[[1]][[i]][1],as.logical(mbnargs[[1]][[i]][2]))
          } else {
            assign(mbnargs[[1]][[i]][1],as.numeric(mbnargs[[1]][[i]][2]))
          }
        }
      }
    } 
    # Check if type is one of the specific allowed values
    allowed_types <- c("classic", "DF", "KC", "MD", "FG", "MBN", "MBB")
    if (!(type1 %in% allowed_types)) {
      stop("The 'type' must be one of the following: 'classic', 'DF', 'KC', 'MD', 'FG', 'MBN', 'MBB'.")
    }
    return(list(type1=type1,exact=exact,r=r,DF=DF,d=d,D=D))
  }
  #
  extract_info <- function(obj,cluster,typeOpts){
    # extract info and bundle it into a list
    n = stats::nobs(obj)
    clusternames = unique(cluster)
    m = length(clusternames)
    X = stats::model.matrix(obj,type="fixed")
    np=dim(X)[2]
    Z = stats::model.matrix(obj,type="random")
    nq=dim(Z)[2]
    Y = obj@resp$y
    if (lme4::isLMM(obj)) nden=rep(1,length(Y)) else nden = obj@resp$n  #allows processing of binomial data
    eps = 1e-14
    info = list(n=n,cluster=cluster,clusternames=clusternames,m=m,X=X,np=np,Z=Z,nq=nq,Y=Y,
                type1=typeOpts$type1,exact=typeOpts$exact,DF=typeOpts$DF,r=typeOpts$r,d=typeOpts$d,D=typeOpts$D,
                nden=nden,link=stats::family(obj)$link,linkinv=stats::family(obj)$linkinv,
                variance=stats::family(obj)$variance,eps=eps)
    return(info)
  }
  #
  CalcRobust <- function(info,WB_B,XtVX,beta,sigma2,eta,ginv_eta) {
    # main robust variance calculations
    if (Matrix::isDiagonal(WB_B)) diagB=TRUE else diagB=FALSE
    WB_C1 = Matrix::solve(XtVX)
    sum=matrix(0,info$np,info$np)
    # start loop over clusters
    W = list(); A = list(); V = list()
    i=0
    for (g in info$clusternames){
      i=i+1
      grp = (info$cluster == g & info$nden>0)
      ng = sum(grp)
      if (info$link == "identity") {
        delta = Matrix::diag(ng)  
        deltainv = delta 
      } else if (info$link == "logit") {
        term = ginv_eta[grp]*(1-ginv_eta[grp])
        delta = Matrix::diag(term,ng,ng)
        deltainv = Matrix::diag(1/term,ng,ng)
      } else if (info$link == "log") {
        term = ginv_eta[grp]
        delta = Matrix::diag(term,ng,ng)
        deltainv = Matrix::diag(1/term,ng,ng)
      } else {
        stop("Link ",link," not supported")
      }
      #
      P = deltainv%*%(info$Y[grp]-ginv_eta[grp]) + eta[grp]
      e = matrix(P - info$X[grp,]%*%beta,ncol=1)
      ete = e[,,drop=FALSE]%*%Matrix::t(e[,,drop=FALSE])
      #
      Sigma = Matrix::diag(sigma2*info$variance(ginv_eta[grp])/info$nden[grp],ng,ng)
      
      if (info$link=="identity") {
        WB_A <- Matrix::diag(1/Matrix::diag(Sigma),ng,ng)
      } else {
        WB_A <- Matrix::diag(1/Matrix::diag(mtx_AD(mtx_DA(deltainv,Sigma),deltainv)),ng,ng)
      }
      
      WB_U <- info$Z[grp,,drop=FALSE] 
      WB_Ut <- Matrix::t(info$Z[grp,,drop=FALSE])
      # Compute the inverse
      if (diagB) {
        WB_AUB <- mtx_AD(mtx_DA(WB_A,WB_U),WB_B)
      } else {
        WB_AUB <- mtx_DA(WB_A,WB_U)%*%WB_B
      }
      WB_UtA <- mtx_AD(WB_Ut,WB_A)
      V[[i]] = WB_U%*%WB_B%*%WB_Ut + mtx_AD(mtx_DA(deltainv,Sigma),deltainv)
      W[[i]] <- WB_A - WB_AUB%*%Matrix::solve(Matrix::diag(info$nq) + WB_Ut%*%WB_AUB)%*%WB_UtA
      #    
      if (info$type1=="MD") {
        WB_U = -Matrix::t(W[[i]])%*%info$X[grp,,drop=FALSE]
        WB_V = Matrix::t(info$X[grp,,drop=FALSE])
        # Since WB_A is identity, the following expressions are simplified from general Woodbury
        O = WB_C1 + WB_V%*%WB_U
        A[[i]] = -(WB_U%*%MASS::ginv(matrix(as.numeric(O),dim(O)))%*%WB_V) 
        Matrix::diag(A[[i]]) = Matrix::diag(A[[i]]) + 1
        AVX = A[[i]]%*%W[[i]]%*%info$X[grp,,drop=FALSE]
        sum = sum + Matrix::t(AVX)%*%ete%*%AVX
      } else {
        if (info$type1=="KC") {
          if (info$exact) {
            H = info$X[grp,,drop=FALSE]%*%XtVX%*%Matrix::t(info$X[grp,,drop=FALSE])%*%W[[i]]  
            A[[i]] = sqrtm_inv(Matrix::diag(ng) - Matrix::t(H))
            if (is.complex(A[[i]])) stop("(I-H_g)^(-1/2) is complex")
            AVX = A[[i]]%*%W[[i]]%*%info$X[grp,,drop=FALSE]
            sum = sum + Matrix::t(AVX)%*%ete%*%AVX
          } else {
            WB_U = -Matrix::t(W[[i]])%*%info$X[grp,,drop=FALSE]
            WB_V = Matrix::t(info$X[grp,,drop=FALSE])
            # Since WB_A is identity, the following expressions are simplified from general Woodbury
            O = WB_C1 + WB_V%*%WB_U
            A[[i]] = -(WB_U%*%MASS::ginv(matrix(as.numeric(O),dim(O)))%*%WB_V) 
            Matrix::diag(A[[i]]) = Matrix::diag(A[[i]]) + 1
            VX = W[[i]]%*%info$X[grp,,drop=FALSE]
            AVX = A[[i]]%*%VX
            term = Matrix::t(AVX)%*%ete%*%VX
            sum = sum + (term + Matrix::t(term))/2
          }
        } else {
          if (info$type1=="FG") {
            A[[i]] = diag(ng)
            Q = Matrix::t(info$X[grp,,drop=FALSE])%*%W[[i]]%*%info$X[grp,,drop=FALSE]%*%XtVX
            XAA = mtx_AD(info$X[grp,],Matrix::diag(1/sqrt(1-pmin(info$r,Matrix::diag(Q))),info$np,info$np))
            sum = sum + Matrix::t(XAA)%*%W[[i]]%*%ete%*%W[[i]]%*%XAA
          } else {
            if (info$type1=="MBB"){
              sqrtVinv = sqrtm_inv(W[[i]])
              term = Matrix::t(sqrtVinv)%*%(V[[i]] - info$X[grp,,drop=FALSE]%*%XtVX%*%Matrix::t(info$X[grp,,drop=FALSE]))%*%sqrtVinv
              A[[i]] = Matrix::t(sqrtVinv)%*%sqrtm_inv(term)%*%sqrtVinv
              AVX = A[[i]]%*%W[[i]]%*%info$X[grp,,drop=FALSE]
              sum = sum + Matrix::t(AVX)%*%ete%*%AVX
            } else {
              # classic, DF and MBN
              A[[i]]=diag(ng)
              VX = W[[i]]%*%info$X[grp,,drop=FALSE]
              sum = sum + Matrix::t(VX)%*%ete%*%VX 
            }}}}
    }
    # end loop over clusters
    c = 1
    deltam = 0
    phi = 0
    if (info$type1=="DF") {
      if (info$m-info$np>0) c = info$m/(info$m-info$np) else cat("DF not valid because m-p <= 0; defaulting to classic")
    }
    if (info$type1=="MBN") {
      f = sum(info$nden)
      if (info$DF) {c = (f-1)/(f-info$np) * (info$m/(info$m-1))}
      if (is.na(info$D)) {
        # standard option for d as described in MBN
        if (info$m > (info$d+1)*info$np) {deltam = info$np/(info$m-info$np)} else {deltam = 1/info$d}
      } else {
        # force d to a particular value regardless of m
        deltam = 1/info$D
      }
      omega = XtVX %*% sum
      evals = Re(eigen(omega,only.values=TRUE)$values)
      if (info$m > info$np) {pstar = info$np} else {pstar = sum(evals>info$eps)}
      phi =  max(info$r,sum(evals)/pstar)
    }
    robustVar = c*XtVX%*%sum%*%XtVX + deltam*phi*XtVX
    return(list(A=A,W=W,V=V,var=robustVar))
  }
  #
  AHZ1df <- function(L,info,M,robust,Fadjust){
    B = list()
    bigW = .bdiag(robust$W)
    bigV = .bdiag(robust$V)
    term = (diag(info$n) - info$X%*%M%*%t(info$X)%*%bigW)
    i=0
    for (g in info$clusternames){
      i=i+1
      grp = (info$cluster == g & info$nden>0)
      ng = sum(grp)
      if (Fadjust) {
        B[[i]] = (1/sqrt(L%*%M%*%t(L)))*L%*%M%*%t(info$X[grp,,drop=FALSE])%*%
        robust$W[[i]]%*%robust$A[[i]]%*%term[grp,,drop=FALSE]
      } else {
        B[[i]] = (1/sqrt(L%*%M%*%t(L)))*L%*%M%*%t(info$X[grp,,drop=FALSE])%*%
          robust$W[[i]]%*%term[grp,,drop=FALSE]
      }
    }
    sum=0
    for (i in 1:info$m){
      for (j in 1:info$m){
        term = B[[i]] %*% bigV %*% t(B[[j]])
        sum = sum + 2*term^2
      }}
    as.numeric(2/sum)
  }
  #
  AHZqdf <- function(L,info,M,robust,Fadjust){
    q = nrow(L)
    B = list()
    bigW = .bdiag(robust$W)
    bigV = .bdiag(robust$V)
    term1 = sqrtm_inv(L%*%M%*%t(L))%*%L%*%M
    term2 = (diag(info$n) - info$X%*%M%*%t(info$X)%*%bigW)
    i=0
    for (g in info$clusternames){
      i=i+1
      grp = (info$cluster == g & info$nden>0)
      ng = sum(grp)
      if (Fadjust) {
        B[[i]] = term1%*%t(info$X[grp,,drop=FALSE])%*%
        robust$W[[i]]%*%robust$A[[i]]%*%term2[grp,,drop=FALSE]
      } else {
        B[[i]] = term1%*%t(info$X[grp,,drop=FALSE])%*%
          robust$W[[i]]%*%term2[grp,,drop=FALSE]
      }
    }
    sum=0
    for (i in 1:info$m){
      for (j in 1:info$m){
        term = B[[i]] %*% bigV %*% t(B[[j]])
        dsum = 0
        for (s in 1:q){
          for (t in 1:q) { 
            dsum = dsum + term[s,t]*term[t,s] + term[s,s]*term[t,t]
          }}
        sum = sum + dsum
      }}
    as.numeric(q*(q+1)/sum)
  }
    #######################
  ## Function starts here
  #######################
  ####  Checks ####
  if ("merMod" %in% class(obj)) {
    stop("The 'obj' should be an object fitted using lmer or glmer.")
  }
  if (!is.null(obj@call$weights)) 
    stop("Models with prior weights are not currently supported.")
  if (!missing(cluster) && !is.factor(cluster)) {
    stop("If 'cluster' is manually input, it must be of class 'factor'.")
  }
  if (missing(cluster)) 
    cluster <- get_outer_group(obj)
  if (!is_nested_lmerMod(obj, cluster)) 
    stop("Non-nested random effects detected. Method is not available for such models.")
  beta = fixef(obj)
  np = length(beta)
  if (!is.matrix(L)) {
    if (length(L) != np) {
      stop("length of L must be equal to number of fixed effects")
    } else {
      L = matrix(L,1,np) 
    }
  } else {
    if (ncol(L)!=np) stop("number of columns of lvec must be equal to number of fixed effects")
  }
  q = nrow(L)
  if (q != rankMatrix(L)) stop("contrasts in L are not independent")
  # Check if type is one of the specific allowed values
  typeOpts <- decode_type(type)
  allowed_types <- c("classic", "DF", "KC", "MD", "FG", "MBN","MBB")
  if (!(typeOpts$type1 %in% allowed_types)) {
    stop("The 'type' must be one of the following: 'classic', 'DF', 'KC', 'MD', 'FG', 'MBN', 'MBB'.")
  }
#  if (typeOpts$type1 %in% c("FG","MBN")) warning("Adjustment matrix A set to identity for df calculation")
  ####### Set things up ########
  info = extract_info(obj,cluster,typeOpts)
  beta=matrix(lme4::fixef(obj),ncol=1)
  eta = stats::predict(obj,type="link")
  ginv_eta = stats::predict(obj,type="response")
  sigma2 = stats::sigma(obj)^2
  lambdat = as.matrix(lme4::getME(obj,"Lambdat"))
  R <- t(lambdat)%*%lambdat*sigma2
  XtVX = stats::vcov(obj)
  robust = CalcRobust(info,R,XtVX,beta,sigma2,eta,ginv_eta)
  #### Compute df
  if (q==1){
    # 1 df case
    estimate <- sum(L %*% beta)
    se.estimate <- sqrt(drop(L %*% robust$var %*% t(L)))
    tstat <- estimate/se.estimate
    ddf <-  max(1,AHZ1df(L,info,XtVX,robust,Fadjust))# denominator DF
    pvalue <- 2 * pt(abs(tstat), df = ddf, lower.tail = FALSE)
    list(estimate=estimate, se=se.estimate, 't stat'= tstat, df=ddf, pvalue = pvalue, robustVar = robust$var)
  } else {
    # multi-df case
    var_Lbeta <- L %*% robust$var %*% t(L)
    z = L%*%beta
    Q = t(z)%*%solve(var_Lbeta)%*%z
    eta <- max(q,AHZqdf(L,info,XtVX,robust,Fadjust))
    Fvalue <- ifelse(kadjust,as.numeric(Q*(eta-q+1)/(eta*q)),as.numeric(Q))
    pvalue <- pf(q=Fvalue, df1=q, df2=eta-q+1, lower.tail=FALSE)
    list('F value'=Fvalue, ndf=q, ddf=eta-q+1, pvalue=pvalue,robustVar = robust$var)
  }
}
