extractFrames_new <- function (formula, data) {
  
  Terms <- terms(formula)
  term_labels <- attr(Terms, "term.labels")
  which_RE <- grep("|", term_labels, fixed = TRUE)
  namesVars <- all.vars(formula)
  respVar <- as.character(formula)[2L]
  
  # Fixed Effects
  
  if(which_RE==1){
    formYx <- paste(1, collapse = " + ")
    formYx <- as.formula(paste(respVar, "~", formYx))
    TermsX <- terms(formYx, data = data)
    mfX <- model.frame(TermsX, data = data)
    TermsX <- terms(mfX)
    X <- model.matrix(TermsX, data)
    
  }else{
    formYx <- paste(term_labels[-which_RE], collapse = " + ")
    formYx <- as.formula(paste(respVar, "~", formYx))
    TermsX <- terms(formYx, data = data)
    mfX <- model.frame(TermsX, data)
    TermsX <- terms(mfX)
    X <- model.matrix(TermsX, data)
  }
  # Random Effects
  spl <- unlist(strsplit(term_labels[which_RE], " | ", fixed = TRUE))
  idVar <- spl[2L]
  data <- data[complete.cases(data[namesVars]), ]
  id <- data[[idVar]]
  id <- match(id, unique(id))
  
  #**********************************************************************
  offset <- as.vector(c(1, 1 + cumsum(tapply(id, id, length))))
  offset_sft = offset - 1
  ZrowsStart <- offset[1:(length(offset)-1)]
  ZrowsEnd = offset_sft[2:length(offset_sft)]
  Zrows = as.matrix(data.frame(ZrowsStart, ZrowsEnd))
  #**********************************************************************
  formYz <- paste(spl[1], collapse = " + ")
  formYz <- as.formula(paste(respVar, "~", formYz))
  TermsZ <- terms(formYz, data = data)
  mfZ <- model.frame(TermsZ, data = data)
  TermsZ <- terms(mfZ)
  Z <- model.matrix(TermsZ, data)
  #**********************************************************************
  Zt = t(Z)
  Zinv = ginv(Z)
  Ztinv = ginv(Zt)
  # individual Z inverse
  Zv = matrix(0, nrow(Zinv), ncol(Zinv))
  for (i in 1:length(ZrowsStart)) {
    index1 = ZrowsStart[i]
    index2 = ZrowsEnd[i]
    Z_sub = Z[index1:index2, ]
    Z_sub_inv = ginv(Z_sub)
    Zv[, index1:index2] = Z_sub_inv
  }
  #**********************************************************************
  # response variable
  y <- model.response(mfX)
  if (is.factor(y))
    y <- as.vector(unclass(y) - 1)
  # hierarchical centering
  find_positions <- function (nams1, nams2) {
    nams1 <- gsub("^", "\\^", nams1, fixed = TRUE)
    vals <- c(glob2rx(nams1), glob2rx(paste0(nams1, ":*")),
              glob2rx(paste0("*:", nams1)))
    out <- sort(unique(unlist(lapply(vals, grep, x = nams2))))
    out
  }
  check_td <- function (x, id) {
    !all(sapply(split(x, id), function (z) all(z - z[1L] < .Machine$double.eps^0.5)))
  }
  has_interceptX <- attr(TermsX, "intercept")
  has_interceptZ <- attr(TermsZ, "intercept")
  performHC <- has_interceptX && (has_interceptX == has_interceptZ)
  if (performHC) {
    terms.labs_X <- attr(TermsX, "term.labels")
    terms.labs_Z <- attr(TermsZ, "term.labels")
    # check for time-varying covariates
    timeTerms <- if (length(terms.labs_Z))
      unlist(lapply(terms.labs_Z, FUN = function(x) grep(x, colnames(X), fixed = TRUE)))
    
    timeTerms <- if (length(terms.labs_Z))
      which(colnames(X)%in%colnames(Z)[-1L])

    
    
    
       
    which_td <- unname(which(apply(X, 2, check_td, id = id)))
    
    ## we define which time columns will be in the Xt matrix 
    
    common_td<-which(which_td%in%timeTerms)
    
    numxt<- which_td[-common_td]
    
    
    all_TDterms <- unique(c(timeTerms, which_td))
    
    
    if(length(which_td)==0){
      baseline <- seq_len(ncol(X))[-1]
      
    }else{
      baseline <- seq_len(ncol(X))[-c(which_td,1)]
    }
    
   
    
    #ind_colmns <- c(list(baseline), lapply(colnames(Z)[-1L], find_positions, 
     #                                      nams2 = colnames(X)))
    #ind_colmns2 <- seq_len(ncol(X))
    #ind_colmns2 <- ind_colmns2[!ind_colmns2 %in% unlist(ind_colmns)]
    
    
    
    #========================================================================================================================================    

    

    
 # we will create an error message. When the Z will contain a baseline variable then the code will stop and will give a message
    

    
      restr<-baseline%in%timeTerms
    
    
    
    if(sum(restr)>=1){
      stop("The package does not support yet models which contain baseline covariates in the design matrix Z of the random effects")
    }
    
  
    
      
    data.id <- data[!duplicated(id), ]
    Xhc <- if (length(terms.labs_Z)) {
      mfHC <- model.frame(TermsX, data = data.id)
      which.timevar <- unique(unlist(lapply(terms.labs_Z, 
                                            FUN = function (x) grep(x, names(mfHC), fixed = TRUE))))
      mfHC[which.timevar] <- lapply(mfHC[which.timevar], 
                                    function (x) { 
                                      if(is.numeric(x)){
                                        
                                        x[]<- 1
                                      }else{
                                        x[]<- x
                                      }; x  })
      
      model.matrix(formYx, mfHC)
    } else {
      model.matrix(formYx, model.frame(TermsX, data = data.id))
    }
  }
  
  
  
  
  
  
  
  environment(TermsX) <- environment(TermsZ) <- NULL
  #**********************************************************************
  if(dim(X)[2]>1){
  Xc = scale(X[, -1], center = TRUE, scale = FALSE) # except intercept
  Xc = cbind(intercept=X[,1], Xc)
  Xs = scale(X[, -1], center = TRUE, scale = TRUE)
  Xs = cbind(intercept=X[,1], Xs)
  }else{
   Xc=X
   Xs=X
  }

  # changed code 
  #==============================================================================================
  
  
  
  if(length(which_td)==0){
    cond <- seq_len(ncol(X))[-1]
    
  }else{
    cond <- seq_len(ncol(X))[-c(all_TDterms,1)]
  }
  
  
  if(length(cond)!=0L){
    XhcC<-Xhc
    XhcC1= scale(XhcC[,cond], center = TRUE, scale = FALSE) # covariates except time
    XhcC[,cond] = XhcC1
    XhcS<-Xhc
    XhcS1 = scale(XhcS[, cond], center = TRUE, scale = TRUE)
    XhcS[, cond] = XhcS1
  }
  
  if(length(cond)==0L){
    XhcC=Xhc
    XhcS=Xhc
  } 
  #=========================================================================================================    
  
  
  if(length(timeTerms)!=0){
    Zc = scale(Z[, -1], center = TRUE, scale = FALSE) # Zc/Zs only involves time
    Zc = cbind(intercept=Z[,1], Zc)
    Zc_inv = ginv(Zc)
    
    Zs = scale(Z[, -1], center = TRUE, scale = TRUE)
    Zs = cbind(intercept=Z[,1], Zs)
    Zs_inv = ginv(Zs)
  }
  
  if(length(timeTerms)==0){
    Zc = Z
    Zc_inv = ginv(Zc)
    Zs = Z
    Zs_inv = ginv(Zs)
  }
  
  
  
  
  
  #changed code 
  
  #**********************************************************************
  if (ncol(X) >2) {
    means_X = apply(X[, -1], 2, mean)
    SDs_X = apply(X[, -1], 2, sd)
    mean_sd_X = means_X/SDs_X
  } 
  
  if (ncol(X) ==2){
    means_X = mean(X[, -1])
    SDs_X = sd(X[, -1])
    mean_sd_X = means_X/SDs_X
  }
  
  
  if (ncol(Z) >2) {
    means_Z = apply(Z[, -1], 2, mean)
    SDs_Z = apply(Z[, -1], 2, sd)
    mean_sd_Z = means_Z/SDs_Z
  } 
  
  if (ncol(Z) ==2){
    means_Z = mean(Z[, -1])
    SDs_Z = sd(Z[, -1])
    mean_sd_Z = means_Z/SDs_Z
  }
  
  
  
  
  if (((length(cond)!=0L)&(length(cond)>1))) {
    means_Xhc = apply(Xhc[, cond], 2, mean)
    SDs_Xhc = apply(Xhc[, cond], 2, sd)
    mean_sd_Xhc = means_Xhc/SDs_Xhc
  }
  
  
  
  if (((length(cond)!=0L)&(length(cond)==1))) {
    means_Xhc = mean(Xhc[, cond])
    SDs_Xhc = sd(Xhc[, cond])
    mean_sd_Xhc = means_Xhc/SDs_Xhc
  }
  
  
  
  condition_td1<-(length(which_td)!=0)&(length(numxt)!=0) 
  if( condition_td1){
    Xhc<-Xhc[,-numxt]
    XhcC<-XhcC[,-numxt]
    XhcS<-XhcS[,-numxt]
    Xt<-as.matrix(X[,numxt])
    colnames(Xt)<-dimnames(X)[[2]][numxt]
    
  }
  
  condition_td2<-(length(which_td)==0)&(length(numxt)==0) 
  if( condition_td2){
    Xhc<-Xhc
    XhcC<-XhcC
    XhcS<-XhcS
  }
  
  condition_td3<-((length(which_td)!=0)&(length(numxt)==0)&(length(common_td)==0))
  if( condition_td3){
    Xhc<-Xhc[,-which_td]
    XhcC<-XhcC[,-which_td]
    XhcS<-XhcS[,-which_td]
    Xt<-as.matrix(X[,which_td])
    colnames(Xt)<-dimnames(X)[[2]][which_td]
  }
  
  
  condition_td4<-((length(which_td)!=0)&(length(numxt)==0)&(length(common_td)!=0))
  if( condition_td4){
    Xhc<-Xhc
    XhcC<-XhcC
    XhcS<-XhcS

  }
  
  
if(is.matrix(Xhc)==FALSE){
  Xhc<-as.matrix(Xhc)
  colnames(Xhc)<-"intercept"
  XhcC<-as.matrix(XhcC)
  colnames(XhcC)<-"intercept"
  XhcS<-as.matrix(XhcS)
  colnames(XhcS)<-"intercept"
  
}
  
  
  
  
  
  if(((length(numxt)!=0L)|condition_td3)){
    XtC = scale(Xt, center = TRUE, scale = FALSE) # covariates except time
    
    XtS = scale(Xt, center = TRUE, scale = TRUE)
    
  }

  
 
  
  
  if (((length(numxt)!=0L)|condition_td3)) {
    means_Xt = apply(Xt, 2, mean)
    SDs_Xt = apply(Xt, 2, sd)
    mean_sd_Xt = means_Xt/SDs_Xt
  }
  
  if ((((length(numxt)!=0L)&(length(numxt)==1))|(condition_td3&(length(which_td)==1)))) {
    means_Xt = mean(Xt)
    SDs_Xt = sd(Xt)
    mean_sd_Xt = means_Xt/SDs_Xt
  }
  
  
  
  

  
  
  
  
  
  # extract results
  
  if(((length(cond)!= 0L)&(length(timeTerms)!= 0L)&(length(numxt)!=0L))){
    
    l<-list(N = nrow(Z), n = length(unique(id)), idVar = idVar, respVar = respVar,
            id = id, ZrowsStart = ZrowsStart, ZrowsEnd = ZrowsEnd,
            y = y, X = X, Xc = Xc, Xs = Xs, XhcC = XhcC, XhcS = XhcS,XtC = XtC, XtS = XtS,Xt=Xt,
            Z_ = Z, Zinv = Zinv, Ztinv = Ztinv, Zv = Zv, Zc = Zc,  Zs = Zs,
            means_X = means_X, SDs_X = SDs_X, mean_sd_X = mean_sd_X,means_Xt = means_Xt, SDs_Xt = SDs_Xt, mean_sd_Xt = mean_sd_Xt,
            means_Z = means_Z, SDs_Z = SDs_Z, mean_sd_Z = mean_sd_Z,
            means_Xhc = means_Xhc, SDs_Xhc = SDs_Xhc, mean_sd_Xhc = mean_sd_Xhc,
            TermsX = TermsX, TermsZ = delete.response(TermsZ), xlev = .getXlevels(TermsX, mfX),
            Xhc = Xhc, 
            ncx = ncol(X), ncz = ncol(Z),ncxt=ncol(Xt),ncxh=ncol(Xhc))
  }
  
  
  if(((length(cond)!= 0L)&(length(timeTerms)== 0L)&(condition_td3==FALSE))){
    
    l<-list(N = nrow(Z), n = length(unique(id)), idVar = idVar, respVar = respVar,
            id = id, ZrowsStart = ZrowsStart, ZrowsEnd = ZrowsEnd,
            y = y, X = X, Xc = Xc, Xs = Xs, XhcC = XhcC, XhcS = XhcS,
            Z_ = Z, Zinv = Zinv, Ztinv = Ztinv,
            means_X = means_X, SDs_X = SDs_X, mean_sd_X = mean_sd_X,
            means_Xhc = means_Xhc, SDs_Xhc = SDs_Xhc, mean_sd_Xhc = mean_sd_Xhc,
            TermsX = TermsX, TermsZ = delete.response(TermsZ), xlev = .getXlevels(TermsX, mfX),
            Xhc = Xhc, 
            ncx = ncol(X), ncz = ncol(Z),ncxh=ncol(Xhc))
  }
  
  
  if(((length(cond)!= 0L)&(length(timeTerms)== 0L)&(condition_td3==TRUE))){
    
    l<-list(N = nrow(Z), n = length(unique(id)), idVar = idVar, respVar = respVar,
            id = id, ZrowsStart = ZrowsStart, ZrowsEnd = ZrowsEnd,
            y = y, X = X, Xc = Xc, Xs = Xs, XhcC = XhcC, XhcS = XhcS,XtC = XtC, XtS = XtS,Xt=Xt,
            Z_ = Z, Zinv = Zinv, Ztinv = Ztinv,
            means_X = means_X, SDs_X = SDs_X, mean_sd_X = mean_sd_X,mean_sd_X = mean_sd_X,means_Xt = means_Xt, SDs_Xt = SDs_Xt, mean_sd_Xt = mean_sd_Xt,
            means_Xhc = means_Xhc, SDs_Xhc = SDs_Xhc, mean_sd_Xhc = mean_sd_Xhc,
            TermsX = TermsX, TermsZ = delete.response(TermsZ), xlev = .getXlevels(TermsX, mfX),
            Xhc = Xhc, 
            ncx = ncol(X), ncz = ncol(Z),ncxh=ncol(Xhc),ncxt=ncol(Xt))
  }
  
  
  
  
  
  
  
  #if(((length(cond)!= 0L)&(length(timeTerms)== 0L)&((length(numxt)!=0L)&(dim(Xt)[2]>=1)))){
    
   # l<-list(N = nrow(Z), n = length(unique(id)), idVar = idVar, respVar = respVar,
         #   id = id, ZrowsStart = ZrowsStart, ZrowsEnd = ZrowsEnd,
          #  y = y, X = X, Xc = Xc, Xs = Xs, XhcC = XhcC, XhcS = XhcS,XtC = XtC, XtS = XtS,Xt=Xt,
          #  Z_ = Z, Zinv = Zinv, Ztinv = Ztinv,
          #  means_X = means_X, SDs_X = SDs_X, mean_sd_X = mean_sd_X,
          #  means_Xhc = means_Xhc, SDs_Xhc = SDs_Xhc, mean_sd_Xhc = mean_sd_Xhc,means_Xt = means_Xt, SDs_Xt = SDs_Xt, mean_sd_Xt = mean_sd_Xt,
           # TermsX = TermsX, TermsZ = delete.response(TermsZ), xlev = .getXlevels(TermsX, mfX),
           # Xhc = Xhc, colmns_HC = ind_colmns, colmns_nHC = ind_colmns2,
            #ncx = ncol(X), ncz = ncol(Z),ncxt=ncol(Xt))
 # }
  
  if(((length(cond)!= 0L)&(length(timeTerms)!= 0L)&(length(numxt)==0L))){
    
    l<-list(N = nrow(Z), n = length(unique(id)), idVar = idVar, respVar = respVar,
            id = id, ZrowsStart = ZrowsStart, ZrowsEnd = ZrowsEnd,
            y = y, X = X, Xc = Xc, Xs = Xs, XhcC = XhcC, XhcS = XhcS,
            Z_ = Z, Zinv = Zinv, Ztinv = Ztinv, Zv = Zv, Zc = Zc,  Zs = Zs,
            means_X = means_X, SDs_X = SDs_X, mean_sd_X = mean_sd_X,
            means_Z = means_Z, SDs_Z = SDs_Z, mean_sd_Z = mean_sd_Z,
            means_Xhc = means_Xhc, SDs_Xhc = SDs_Xhc, mean_sd_Xhc = mean_sd_Xhc,
            TermsX = TermsX, TermsZ = delete.response(TermsZ), xlev = .getXlevels(TermsX, mfX),
            Xhc = Xhc, 
            ncx = ncol(X), ncz = ncol(Z),ncxh=ncol(Xhc))
    
  }
  
  
  
  
  if(((length(cond)== 0L)&(length(timeTerms)== 0L)&(condition_td3==FALSE))){
    
    l<-list(N = nrow(Z), n = length(unique(id)), idVar = idVar, respVar = respVar,
            id = id, ZrowsStart = ZrowsStart, ZrowsEnd = ZrowsEnd,
            y = y, X = X, Xc = Xc, Xs = Xs, XhcC = XhcC, XhcS = XhcS,
            Z_ = Z, Zinv = Zinv, Ztinv = Ztinv,
            
            
            TermsX = TermsX, TermsZ = delete.response(TermsZ), xlev = .getXlevels(TermsX, mfX),
            Xhc = Xhc, 
            ncx = ncol(X), ncz = ncol(Z),ncxh=ncol(Xhc))
  }
  
  
  if(((length(cond)== 0L)&(length(timeTerms)== 0L)&(condition_td3==TRUE))){
    
    l<-list(N = nrow(Z), n = length(unique(id)), idVar = idVar, respVar = respVar,
            id = id, ZrowsStart = ZrowsStart, ZrowsEnd = ZrowsEnd,
            y = y, X = X, Xc = Xc, Xs = Xs, XhcC = XhcC, XhcS = XhcS,XtC = XtC, XtS = XtS,Xt=Xt,
            Z_ = Z, Zinv = Zinv, Ztinv = Ztinv,means_Xt = means_Xt, SDs_Xt = SDs_Xt, mean_sd_Xt = mean_sd_Xt,
            
            
            TermsX = TermsX, TermsZ = delete.response(TermsZ), xlev = .getXlevels(TermsX, mfX),
            Xhc = Xhc, 
            ncx = ncol(X), ncz = ncol(Z),ncxh=ncol(Xhc),ncxt=ncol(Xt))
  }
  
  
  
  return(l)
  
  
  #**********************************************************************
  
} 
