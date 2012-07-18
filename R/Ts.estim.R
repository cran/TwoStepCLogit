#' Two-step estimator
#' 
#' Function that computes the two-step estimator proposed in Craiu et al.
#' (2011) and its \code{print} method.
#' 
#' Calls \code{\link{coxph}} from the package \pkg{survival}.
#' 
#' @param formula A formula object, with the response on the left of a \code{~} operator, 
#'                and the covariates on the right. The right hand side of the model must also 
#'                include two special terms: a \code{strata} and a \code{cluster} term 
#'                (ex. \code{formula = Y ~ X1 + X2 + X3 + strata(var_strata) + cluster(var_cluster)}).
#'                The \code{strata} and \code{cluster} functions (from the package \pkg{survival}) are
#'                used to identify the stratification and the cluster variables, respectively.
#' @param data A data frame (or object coercible by as.data.frame to a data frame) containing the 
#'             variables in the model.
#' @param random  A formula object, with a blank on the left of a \code{~} operator, 
#'                and, on the right, the covariates with random coefficients among the 
#'                covariate listed in the model \code{formula} (ex. \code{random = ~ X1 + X3}).
#'                The default is to add random coefficients for every covariates listed in 
#'                the model \code{formula}.
#' @param all.m.1 \code{TRUE} if sum of Y's in all strata is 1, \code{FALSE} otherwise 
#'                (the default). When in doubt use \code{FALSE} (always works, but slower than
#'                necessary if all stratum sums are 1).
#' @param D The form of the between-cluster variance-covariance matrix of the regression 
#'        coefficients (matrix D) : either \code{"UN"} for unstructured matrix D or \code{"UN(1)"} 
#'        (the default) for diagonal matrix D.
#' @param itermax maximal number of EM iterations (default = 2000)
#' @param tole maximal distance between successive EM iterations tolerated
#'             before declaring convergence (default = 0.000001) 
#' 
#' @return \item{beta}{ A vector: the regression coefficients. } 
#' @return \item{se}{ A vector: the regression coefficients' standard errors. } 
#' @return \item{D}{ A matrix: estimate of the between-cluster variance-covariance matrix of the regression coefficients (matrix D).}
#' @return \item{r.effect}{ The random effect estimates. }
#' @return \item{coxph.warn}{ A list of character string vectors. If the \code{\link{coxph}} 
#'                            function generates one or more warnings when fitting the Cox
#'                            model to a cluster, a copy of these warnings are stored in 
#'                            \code{coxph.warn$Cluster_name} where \code{Cluster_name} is 
#'                            the identification value for the cluster in the data set. 
#'                            A \code{NULL} list element means that \code{\link{coxph}} did not 
#'                            produce any warnings for that cluster. }
#' @return \item{Call}{ The function call.}
#' 
#' @author Radu V. Craiu, Thierry Duchesne, Daniel Fortin and Sophie
#' Baillargeon
#' @references Craiu, R.V., Duchesne, T., Fortin, D. and Baillargeon, S.
#' (2011), Conditional Logistic Regression with Longitudinal Follow-up and
#' Individual-Level Random Coefficients: A Stable and Efficient Two-Step
#' Estimation Method, \emph{Journal of Computational and Graphical Statistics}. 
#' \bold{20}(3), 767-784.
#' @seealso \code{\link{ddim}}
#' @keywords models
#' @export
#' @importFrom survival coxph Surv strata cluster
#' @examples
#' data(bison)
#' 
#' # Two ways for specifying the same model
#' # Model: covariates forest, biomass and pmeadow
#' # Random effects in front of forest and biomass
#' # Main diagonal covariance structure for D (the default)
#' way1 <- Ts.estim(formula = Y ~ forest  + biomass + pmeadow + 
#'         strata(Strata) + cluster(Cluster), data = bison, 
#'         random = ~ forest + biomass)
#' way1
#' way2 <- Ts.estim(formula = bison[,3] ~ as.matrix(bison[,c(6,8:9)]) + 
#'         strata(bison[,2]) + cluster(bison[,1]), data = bison, 
#'         random = ~ as.matrix(bison[,c(6,8)]))
#' way2
#' 
#' # Unstructured covariance for D
#' Fit <- Ts.estim(formula = Y ~ forest  + biomass + pmeadow + 
#'         strata(Strata) + cluster(Cluster), data = bison, 
#'         random = ~ forest + biomass, D="UN")
#' Fit
Ts.estim <- function(formula, data, random, all.m.1=FALSE, D="UN(1)", itermax=2000, tole=0.000001){
  
  call <- scall <- match.call()
  
  simplify <- function(label, data.name){
  # Fonction pour simplifier les noms de variables :
    # Sous-fonctions utile pour conserver un seul élément d'un vecteur
    keeplast <- function(x) x[length(x)]
    keepsec <- function(x) if (length(x)>1) x[2] else x
    # Pour enlever "data.name$" si présent
    label.split <- strsplit(label, split="$", fixed = TRUE)
    label <- unlist(lapply(label.split, keeplast))
    # Si as.matrix est présent, les noms auront l'allure as.matrix(xxx)... : conserver seulement la fin ...
    label.split <- ifelse(grepl("as.matrix(", label, fixed=TRUE), strsplit(label, split=")", fixed = TRUE), label)
    label <- unlist(lapply(label.split, keeplast))
    # Pour conserver seulement ... dans "data.name[, \"...\"]"
    label.split <- strsplit(label, split="\"", fixed = TRUE)
    label <- unlist(lapply(label.split, keepsec))
    # S'il y a encore "[", c'est qu'il y a un numéro de colonne,
    # je veux le remplacer par le nom de la colonne, ci celui-ci existe
    test <- paste("!(is.null(colnames(", data.name, ")) || all(colnames(", data.name, ") == \"\"))")
    # pour tester si les noms de colonne sont NULL ou tous égaux à une chaîne de caractères vide
    if (eval(parse(text = test))) {
      label_m <- sub(data.name, paste("colnames(", data.name, ")", sep=""), label, fixed=TRUE)
      label_m <- sub(",", "", label_m, fixed=TRUE)
      submit <- function(x) if (grepl("colnames", x, fixed=TRUE)) eval(parse(text=x)) else x
      label <- unlist(lapply(label_m, submit))
    }
    return(label)
  }  
  
  ########################################################################
  # Validation des arguments et création de variables
  
  # Validation des arguments formula et data
  data.name <- info.cluster <- info.strata <- y <- mm <- var.cluster <- NULL
  sargs <- match(c("formula", "data"), names(scall), 0L)
  scall <- scall[c(1L, sargs)]
  scall[[1L]] <- as.name("shared")
  out.shared <- eval(scall)
  for(i in 1:length(out.shared)) assign(names(out.shared)[i],out.shared[[i]])  
  
  # Création des variables covar.labels et p
  if (ncol(mm)==0) stop("at least one covariate must be included in the model")  
  covar.labels <- dimnames(mm)[[2]]
  covar.labels <- simplify(covar.labels, data.name)
  p <- length(covar.labels)  # number of beta coefficients
  
  # Créations de variables relatives aux clusters
  Clusters <- unique(var.cluster)  # cluster identification
  k <- length(Clusters)  # number of clusters
    
  # Validation de l'argument 'random' et création des variables random.labels, rpos et q
  if (is.null(call$random)) {
    random.labels <- covar.labels
    rpos <- 1:length(covar.labels)
  } else {
    if (class(random)!="formula") stop("the 'random' argument must be a formula")
    random_noi <- update.formula(old=random, new = ~ . - 1)
    mfr <- model.frame(random_noi, data=data, na.action=NULL)
    mmr <- model.matrix(random_noi, data=mfr)
    random.labels <- dimnames(mmr)[[2]]
    if (length(random.labels) == 0) stop("at least one covariate must have a random coefficient")
    random.labels <- simplify(random.labels, data.name)
    if(!all(random.labels %in% covar.labels))
      stop("variables in 'random' must also be in 'formula'")
    rpos <- match(random.labels, covar.labels)
  }
  q <- length(random.labels) # number of random coefficients
  
  ########################################################################

  
  ### Pour ajuster le modèle de cox séparément par cluster
  
  # Étape 1 - Préparation pour les appels à coxph() :
  method <- if (all.m.1) "efron" else "exact"
  # Enlever de la formule le terme cluster et modifier la variable réponse
  appel.formula <- paste("update.formula(old=formula, new=Surv(2-.,.) ~ . -",info.cluster$vars,")") 
  formula_coxph <- eval(parse(text=appel.formula))
  
  # Étape 2 - appeler coxph avec la formule ci-dessus pour chaque cluster
  betas <- vector(length=k*p)  # vector of beta_{ij}, i=1,...,k, j=1,...,p
  R <- matrix(0, ncol=k*p, nrow=k*p)  # variance matrix of betas
  R.sim <- R
  coxph.warn <- vector(length=k, mode="list")
  names(coxph.warn) <- Clusters
  for(i in 1:k){
    assign(data.name, data[var.cluster==Clusters[i], ])  ## ça remplace le jeu de données data complet par le
                                                         ## sous jeu de données seulement le cluster traité
    appel.coxph <- paste("coxph(",paste(deparse(formula_coxph),collapse=""),", data=",data.name,", method=method)",sep="")
    try.model.fit <- tryCatch.W.E(eval(parse(text=appel.coxph)))
    if (identical(class(try.model.fit$value), "erreur"))
      stop("\nError message generated by the function coxph for cluster ",Clusters[i],":\n  ",
           gettext(try.model.fit$value$message),
           "\nTo solve the problem, you should remove the problematic cluster or variables from the data set.", call.=FALSE)  
    model.fit <- try.model.fit$value
    if (!is.null(try.model.fit$warnings)) coxph.warn[[i]] <- try.model.fit$warnings  
      # car affecter NULL à un élément d'une liste efface cet élément, je ne veux pas ça   
    pos <- ((i-1)*p+1):(i*p)
    betas[pos] <- model.fit$coefficients
    R[pos,pos] <- model.fit$var 
    R.sim[pos,pos] <- if (p==1) model.fit$var else diag(diag(model.fit$var))
  }  
  
  ### code for REML estimation
  W1 <- kronecker(rep(1,k),diag(rep(1,q))) # fixed effects design matrix
  M <- diag(rep(1,k*q))-W1%*%t(W1)/k # matrix to "kill" the fixed effects
  M <- M[,(1:(q*(k-1)))]
  gammA <- t(M)%*%betas[rep(p*(0:(k-1)),each=q)+rep(rpos,k)] # response with fixed effects removed
  if (q==p) { RR <- R.sim } else { 
    RR <- diag(diag(R.sim)[rep(p*(0:(k-1)),each=q)+rep(rpos,k)]) # subset of R.sim corresponding to the random regression coefficients
  }
  
  # EM algo for either diagonal (UN(1)) or 
  # unstructured (UN) form for D
  if(D=="UN(1)"){
    # EM-algorithm 
    # Initial values
    eta.0 <- rep(1,q)
    eta.1 <- eta.0
    iterations <- 0
    abs.error <- 99999
    # E & M steps ... stops when itermax reach or tole reached
    while((iterations < itermax)*(abs.error>tole)){
      D.0 <- diag(rep(eta.0,k))
      MMRM.inv <- M%*%solve(t(M)%*%RR%*%M)
      Sigma <- solve(MMRM.inv%*%t(M)+solve(D.0))
      Mu <- Sigma%*%MMRM.inv%*%gammA
      for(j in 1:q){
        A <- diag(Sigma+Mu%*%t(Mu))
        eta.1[j] <- mean(A[((1:k)-1)*q+j])
      }
      abs.error <- max(abs(eta.0-eta.1))
      iterations <- iterations+1
      eta.0 <- eta.1
    }
    DD <- diag(rep(eta.0,k))
    DD.block <- if(length(eta.0)==1) matrix(eta.0,1,1) else diag(eta.0)
  }
  if(D=="UN"){
    # EM-algorithm 
    # Initial values
    D0.block <- diag(rep(1,q))
    D1.block <- D0.block
    D.0 <- kronecker(diag(rep(1,k)),D0.block)
    iterations <- 0
    abs.error <- 99999
    # E & M steps ... stops when itermax reach or tole reached
    while((iterations < itermax)*(abs.error>tole)){
      MMRM.inv <- M%*%solve(t(M)%*%RR%*%M)
      Sigma <- solve(MMRM.inv%*%t(M)+solve(D.0))
      Mu <- Sigma%*%MMRM.inv%*%gammA
      A <- Sigma + Mu%*%t(Mu)
      SUM <- 0
      for(cc in 1:k){
        Vcc <- A[(((cc-1)*q+1):(cc*q)),(((cc-1)*q+1):(cc*q))]
        SUM <- SUM + Vcc
      }
      D1.block <- SUM/k
      abs.error <- max(abs(D0.block-D1.block))
      iterations <- iterations+1
      D.0 <- kronecker(diag(rep(1,k)),D1.block)    
      #print(cat(paste(c(iterations,abs.error))))    
    }
    DD <- D.0
    DD.block <- D1.block
  }
  # Now estimation of the regression coefficients
  if (q==p) { D.sim <- DD; D.sim.block <- DD.block } else {
    D.sim.block <- matrix(0,p,p)
    D.sim.block[rpos,rpos] <- DD.block
    D.sim <- matrix(0,k*p,k*p)
    for(cc in 1:k){
      D.sim[rpos+(cc-1)*p,rpos+(cc-1)*p] <- DD[(q*(cc-1)+1):(cc*q),(q*(cc-1)+1):(cc*q)]
    }
  }
  V.sim <- R.sim + D.sim # Normally R + ZDZ', but here Z=identity ...
  Q <- kronecker(rep(1,k),diag(rep(1,p)))
  Var.Beta.sim <- solve(t(Q)%*%solve(V.sim)%*%Q)
  BetaWLS.simREML <- Var.Beta.sim%*%t(Q)%*%solve(V.sim)%*%betas
  se <- sqrt(diag(Var.Beta.sim))
  
  names(BetaWLS.simREML) <- names(se) <- rownames(D.sim.block) <- colnames(D.sim.block) <- covar.labels
  r.effect <- matrix(Mu,byrow=TRUE,ncol=q)
  rownames(r.effect) <- 1:k
  colnames(r.effect) <- random.labels
  outp <- list(beta=c(BetaWLS.simREML),se=se,D=D.sim.block,r.effect=r.effect,coxph.warn=coxph.warn,call=call)
  class(outp) <- "Ts.estim"
  return(outp)  
}


#' @rdname Ts.estim
#' @method print Ts.estim
#' @param x An object, produced by the \code{\link{Ts.estim}} function, to print.
#' @param \dots Further arguments to be passed to \code{print.default}. 
#' @export
"print.Ts.estim" <- function(x, ...) {
  cat("Call:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
    
  cat("beta coefficients:\n")
  tableau2 <- cbind(estimate=round(x$beta,6),se=round(x$se,6))
  rownames(tableau2) <- names(x$beta)
  print.default(tableau2, print.gap = 2, quote = FALSE, right=TRUE, ...)
  
  cat("\nD = estimate of the between-cluster variance-covariance matrix D,\n    for the random coefficients only:\n")  
  keep <- colnames(x$D)[rowSums(x$D)!=0]
  D <- x$D[keep, keep, drop=FALSE]
  print.default(D, print.gap = 2, quote = FALSE, right=TRUE, ...)
    
  invisible(x)
}

