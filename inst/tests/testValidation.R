context("Arguments validation")

test_that("the 'data' argument is validated correctly", {      
      # the data argument can be a matrix
#      bison.m <- as.matrix(bison)
#      Tsout <- Ts.estim(formula = Y ~ biomass + I(biomass^2) + strata(Strata) + cluster(Cluster), data = bison.m)
#      expect_that(round(Tsout$beta,2), equals(c("biomass"=8.59, "I(biomass^2)"=-7.79)))
#      ddout <- ddim(formula = Y ~ strata(Strata) + cluster(Cluster), data = bison.m)
#      expect_that(length(ddout$Sc), equals(20))
# probl?me avec testthat que je ne comprends pas...
      
      # the data argument cannot be missing
      expect_that(Ts.estim(formula = Y ~ biomass + strata(Strata) + cluster(Cluster)),
          throws_error("a 'data' argument is required"))
      expect_that(ddim(formula = Y ~ strata(Strata) + cluster(Cluster)),
          throws_error("a 'data' argument is required"))
    
      # the data argument must be an object, not a function call
      expect_that(Ts.estim(formula = Y ~ biomass + strata(Strata) + cluster(Cluster), data = as.matrix(bison)),
          throws_error("the 'data' argument must be an object, not a function call"))
      expect_that(ddim(formula = Y ~ strata(Strata) + cluster(Cluster), data = as.matrix(bison)),
          throws_error("the 'data' argument must be an object, not a function call"))
    })


test_that("the 'formula' argument is validated correctly by Ts.estim", {      
      # the formula argument cannot be missing
      expect_that(Ts.estim(data=bison),
          throws_error("a 'formula' argument is required"))
      # the formula argument must be a formula
      expect_that(Ts.estim(formula="test", data=bison),
          throws_error("the 'formula' argument must be a formula"))
      # missing cluster term
      expect_that(Ts.estim(formula = Y ~ biomass + I(biomass^2) + strata(Strata), data = bison),
          throws_error("the formula must include a cluster term"))
      # missing strata term
      expect_that(Ts.estim(formula = Y ~ biomass + I(biomass^2) + cluster(Cluster), data = bison),
          throws_error("the formula must include a strata term"))
      # no covariates given
      expect_that(Ts.estim(formula = Y ~ strata(Strata) + cluster(Cluster), data = bison),
          throws_error("at least one covariate must be included in the model"))
      expect_that(Ts.estim(formula = Y ~ 1 + strata(Strata) + cluster(Cluster), data = bison),
          throws_error("at least one covariate must be included in the model"))
      # no response variable given
      expect_that(Ts.estim(formula = ~ biomass + strata(Strata) + cluster(Cluster), data = bison),
          throws_error("a response variable must be given"))
      # a non binary response variable is given
      expect_that(Ts.estim(formula = biomass ~ meadow + strata(Strata) + cluster(Cluster), data = bison),
          throws_error("the response variable can only take the values 0 and 1"))
      # wrong variable name given
      expect_that(Ts.estim(formula = Y ~ biomass + strata(Strato) + cluster(Cluster), data = bison),
          throws_error("variable names figuring in 'formula' should be found among column names of 'data'"))
      # more than one variable given for cluster identification
      expect_that(Ts.estim(formula = Y ~ biomass + strata(Strata) + cluster(Cluster + meadow), data = bison),
          throws_error("in 'formula', the cluster identifier must be a single variable"))
      expect_that(Ts.estim(formula = Y ~ biomass + strata(Strata) + cluster(bison[,c("Cluster","meadow")]), data = bison),
          throws_error("in 'formula', the cluster identifier must be a single variable"))
      
      Ts.estim(formula = Y ~ biomass + strata(Strata + meadow) + cluster(Cluster), data = bison)
      test <- ddim(formula = Y ~ strata(Strata + meadow) + cluster(Cluster), data = bison)
      })


test_that("the 'formula' argument is validated correctly by ddim", {      
      # ddim works even when a covariate is given
      ddout1 <- ddim(formula = Y ~ biomass + strata(Strata) + cluster(Cluster), data = bison)
      expect_that(ncol(ddout1$Ystat), equals(4))
      # ddim works even when the response variable is not standard
      ddout2 <- ddim(formula = meadow ~ strata(Strata) + cluster(Cluster), data = bison)
      expect_that(ncol(ddout2$Ystat), equals(4))
      # ddim works the addition of two variables is given to strata
      ddout3 <- ddim(formula = Y ~ strata(Strata + meadow) + cluster(Cluster), data = bison)
      expect_that(ncol(ddout3$Ystat), equals(4))
      
      # the formula argument cannot be missing
      expect_that(ddim(data=bison),
          throws_error("a 'formula' argument is required"))
      # the formula argument must be a formula
      expect_that(ddim(formula="test", data=bison),
          throws_error("the 'formula' argument must be a formula"))
      # missing cluster term
      expect_that(ddim(formula = Y ~ strata(Strata), data = bison),
          throws_error("the formula must include a cluster term"))
      # missing strata term
      expect_that(ddim(formula = Y ~ cluster(Cluster), data = bison),
          throws_error("the formula must include a strata term"))
      # no response variable given
      expect_that(ddim(formula = ~ biomass + strata(Strata) + cluster(Cluster), data = bison),
          throws_error("a response variable must be given"))
      # a non binary response variable is given
      expect_that(ddim(formula = biomass ~ strata(Strata) + cluster(Cluster), data = bison),
          throws_error("the response variable can only take the values 0 and 1"))
      # wrong variable name given
      expect_that(ddim(formula = Y ~ strata(Strato) + cluster(Cluster), data = bison),
          throws_error("variable names figuring in 'formula' should be found among column names of 'data'"))
      # more than one variable given for cluster identification
      expect_that(ddim(formula = Y ~ strata(Strata) + cluster(Cluster + meadow), data = bison),
          throws_error("in 'formula', the cluster identifier must be a single variable"))
      expect_that(ddim(formula = Y ~ strata(Strata) + cluster(bison[,c("Cluster","meadow")]), data = bison),
          throws_error("in 'formula', the cluster identifier must be a single variable"))
      # more than one variable given for strata identification
      expect_that(ddim(formula = Y ~ strata(bison[,c("Strata","meadow")]) + cluster(Cluster), data = bison),
          throws_error("in 'formula', the stratum identifier must be a single variable"))
})
    

test_that("the 'random' argument is validated correctly", {      
      # the random argument must be a formula
      expect_that(Ts.estim(formula = Y ~ biomass + strata(Strata) + cluster(Cluster), data = bison, random = "test"),
          throws_error("the 'random' argument must be a formula"))
      # no random term given
      expect_that(Ts.estim(formula = Y ~ biomass + strata(Strata) + cluster(Cluster), data = bison, random = ~ 1),
          throws_error("at least one covariate must have a random coefficient"))
      expect_that(Ts.estim(formula = Y ~ biomass + strata(Strata) + cluster(Cluster), data = bison, random = ~ 0),
          throws_error("at least one covariate must have a random coefficient"))
      # unaccepted variable name in random 
      expect_that(Ts.estim(formula = Y ~ biomass + strata(Strata) + cluster(Cluster), data = bison, random = ~ meadow),
          throws_error("variables in 'random' must also be in 'formula'"))
    })


test_that("the Ts.estim function works correctly", {    
      # only one covariate
      one <- Ts.estim(formula = Y ~ biomass + strata(Strata) + cluster(Cluster), data = bison)
      expect_that(dim(one$D), equals(c(1,1)))
      
      # data wihtout colnames
#      bison.nn <- bison
#      colnames(bison.nn) <- NULL
#      nn <- Ts.estim(formula = bison.nn[,3] ~ bison.nn[,8] + strata(bison.nn[,2]) + cluster(bison.nn[,1]), data = bison.nn,
#                     random = ~ bison.nn[,8])
#      # des warnings sont g?n?r?s, mais la r?ponse est bonne.
#      expect_that(nn$D, is_equivalent_to(one$D))
# probl?me avec testthat que je ne comprends pas...
      
      # Different ways of submitting a model with one covariate give the same result
#      logic8 <- c(rep(FALSE,7), TRUE, FALSE)
#      logic2 <- c(FALSE, TRUE, rep(FALSE,7))
#      logic1 <- c(TRUE, rep(FALSE,8))
      way1 <- Ts.estim(formula = Y ~ biomass + strata(Strata) + cluster(Cluster), data = bison, random = ~ bison$biomass)
      way2 <- Ts.estim(formula = Y ~ biomass + strata(Strata) + cluster(Cluster), data = bison, random = ~ bison[,"biomass"])
      way3 <- Ts.estim(formula = Y ~ biomass + strata(Strata) + cluster(Cluster), data = bison, random = ~ bison[,8])
#      way4 <- Ts.estim(formula = Y ~ biomass + strata(Strata) + cluster(Cluster), data = bison, random = ~ bison[,logic8])
# probl?me avec testthat que je ne comprends pas...
      way5 <- Ts.estim(formula = Y ~ bison$biomass + strata(Strata) + cluster(Cluster), data = bison, random = ~ biomass)
      way6 <- Ts.estim(formula = Y ~ bison[,"biomass"] + strata(Strata) + cluster(Cluster), data = bison, random = ~ biomass)
      way7 <- Ts.estim(formula = Y ~ bison[,8] + strata(Strata) + cluster(Cluster), data = bison, random = ~ biomass)
#      way8 <- Ts.estim(formula = Y ~ bison[,logic8] + strata(Strata) + cluster(Cluster), data = bison, random = ~ biomass)
# ne fonctionne pas, mais fonctionne si c(rep(FALSE,7), TRUE, FALSE) est donn? directement plut?t que de passer par logic8
      way9 <- Ts.estim(formula = Y ~ biomass + strata(bison$Strata) + cluster(Cluster), data = bison, random = ~ biomass)
      way10 <- Ts.estim(formula = Y ~ biomass + strata(bison[,"Strata"]) + cluster(Cluster), data = bison, random = ~ biomass)
      way11 <- Ts.estim(formula = Y ~ biomass + strata(bison[,2]) + cluster(Cluster), data = bison, random = ~ biomass)
#      way12 <- Ts.estim(formula = Y ~ biomass + strata(bison[,logic2]) + cluster(Cluster), data = bison, random = ~ biomass)
# ne fonctionne pas, mais ce n'est pas grave car il y a plein d'autres fa?on de sp?cifier les variables
# et une erreur relativement informative est g?n?r?e
      way13 <- Ts.estim(formula = Y ~ biomass + strata(Strata) + cluster(bison$Cluster), data = bison, random = ~ biomass)
      way14 <- Ts.estim(formula = Y ~ biomass + strata(Strata) + cluster(bison[,"Cluster"]), data = bison, random = ~ biomass)
      way15 <- Ts.estim(formula = Y ~ biomass + strata(Strata) + cluster(bison[,1]), data = bison, random = ~ biomass)
#      way16 <- Ts.estim(formula = Y ~ biomass + strata(Strata) + cluster(bison[,logic1]), data = bison, random = ~ biomass)
# ne fonctionne pas, mais ce n'est pas grave car il y a plein d'autres fa?on de sp?cifier les variables
# et une erreur relativement informative est g?n?r?e
      way17 <- Ts.estim(formula = bison$Y ~ biomass + strata(Strata) + cluster(Cluster), data = bison, random = ~ biomass)
      way18 <- Ts.estim(formula = bison[,"Y"] ~ biomass + strata(Strata) + cluster(Cluster), data = bison, random = ~ biomass)
      way19 <- Ts.estim(formula = bison[,3] ~ biomass + strata(Strata) + cluster(Cluster), data = bison, random = ~ biomass)
      res <- c(way1$D, way2$D, way3$D, way5$D, way6$D, way7$D, way9$D, way10$D, way11$D, way13$D, way14$D, way15$D,
          way17$D, way18$D, way19$D)
      expect_that(res, equals(rep(one$D, 15)))
      resname <- c(colnames(way1$D), colnames(way2$D), colnames(way3$D), colnames(way5$D), colnames(way6$D), colnames(way7$D))
      expect_that(resname, equals(rep(colnames(one$D), 6)))      

      
      # Different ways of submitting a model with more thanone covariate give the same result
      many <- Ts.estim(formula = Y ~ forest + biomass + pmeadow + strata(Strata) + cluster(Cluster), data = bison)
      
      expect_that(Ts.estim(formula = Y ~ bison[,c(6,8,9)] + strata(Strata) + cluster(Cluster), data = bison),
          throws_error())
      expect_that(Ts.estim(formula = Y ~ forest + biomass + pmeadow + strata(Strata) + cluster(Cluster), data = bison,
              random = ~ bison[,c(6,8,9)]), throws_error())
      expect_that(Ts.estim(formula = Y ~ bison[,c("forest","biomass","pmeadow")] + strata(Strata) + cluster(Cluster), data = bison),
          throws_error())
      expect_that(Ts.estim(formula = Y ~ forest + biomass + pmeadow + strata(Strata) + cluster(Cluster), data = bison,
              random = ~ bison[,c("forest","biomass","pmeadow")]), throws_error())
      
#      mat <- as.matrix(bison[,c(6,8,9)])
#      way <- Ts.estim(formula = Y ~ mat + strata(Strata) + cluster(Cluster), data = bison)
#      indic <- c(6,8,9)
#      way <- Ts.estim(formula = Y ~ as.matrix(bison[,indic]) + strata(Strata) + cluster(Cluster), data = bison)
# ces fa?ons d'entrer les arguments ne fonctionnent pas : ok car une erreur relativement informative est g?n?r?e

      way1 <- Ts.estim(formula = Y ~ forest + biomass + pmeadow + strata(Strata) + cluster(Cluster), data = bison,
          random = ~ forest + biomass + pmeadow)
      way2 <- Ts.estim(formula = Y ~ as.matrix(bison[,c(6,8,9)]) + strata(Strata) + cluster(Cluster), data = bison,
          random = ~ forest + biomass + pmeadow)
      way3 <- Ts.estim(formula = Y ~ as.matrix(bison[,c("forest","biomass","pmeadow")]) + strata(Strata) + cluster(Cluster), 
          data = bison, random = ~ forest + biomass + pmeadow)
      way4 <- Ts.estim(formula = Y ~ as.matrix(bison[,c(rep(FALSE,5), TRUE, FALSE, TRUE, TRUE)]) + strata(Strata) + cluster(Cluster), 
          data = bison, random = ~ forest + biomass + pmeadow)
      way5 <- Ts.estim(formula = Y ~ forest + biomass + pmeadow + strata(Strata) + cluster(Cluster), data = bison,
          random = ~ forest + biomass + pmeadow)
      way6 <- Ts.estim(formula = Y ~ as.matrix(bison[,c(6,8,9)]) + strata(Strata) + cluster(Cluster), data = bison,
          random = ~ as.matrix(bison[,c(6,8,9)]))
      way7 <- Ts.estim(formula = Y ~ as.matrix(bison[,c("forest","biomass","pmeadow")]) + strata(Strata) + cluster(Cluster), 
          data = bison, random = ~ as.matrix(bison[,c(6,8,9)]))
      way8 <- Ts.estim(formula = Y ~ as.matrix(bison[,c(rep(FALSE,5), TRUE, FALSE, TRUE, TRUE)]) + strata(Strata) + cluster(Cluster), 
          data = bison, random = ~ as.matrix(bison[,c(6,8,9)]))
      way9 <- Ts.estim(formula = Y ~ forest + biomass + pmeadow + strata(Strata) + cluster(Cluster), data = bison,
          random = ~ as.matrix(bison[,c("forest","biomass","pmeadow")]))
      way10 <- Ts.estim(formula = Y ~ as.matrix(bison[,c(6,8,9)]) + strata(Strata) + cluster(Cluster), data = bison,
          random = ~ as.matrix(bison[,c("forest","biomass","pmeadow")]))
      way11 <- Ts.estim(formula = Y ~ as.matrix(bison[,c("forest","biomass","pmeadow")]) + strata(Strata) + cluster(Cluster), 
          data = bison, random = ~ as.matrix(bison[,c("forest","biomass","pmeadow")]))
      way12 <- Ts.estim(formula = Y ~ as.matrix(bison[,c(rep(FALSE,5), TRUE, FALSE, TRUE, TRUE)]) + strata(Strata) + cluster(Cluster), 
          data = bison, random = ~ as.matrix(bison[,c("forest","biomass","pmeadow")]))
      way13 <- Ts.estim(formula = Y ~ forest + biomass + pmeadow + strata(Strata) + cluster(Cluster), data = bison,
          random = ~ as.matrix(bison[,c(rep(FALSE,5), TRUE, FALSE, TRUE, TRUE)]))
      way14 <- Ts.estim(formula = Y ~ as.matrix(bison[,c(6,8,9)]) + strata(Strata) + cluster(Cluster), data = bison,
          random = ~ as.matrix(bison[,c(rep(FALSE,5), TRUE, FALSE, TRUE, TRUE)]))
      way15 <- Ts.estim(formula = Y ~ as.matrix(bison[,c("forest","biomass","pmeadow")]) + strata(Strata) + cluster(Cluster), 
          data = bison, random = ~ as.matrix(bison[,c(rep(FALSE,5), TRUE, FALSE, TRUE, TRUE)]))
      way16 <- Ts.estim(formula = Y ~ as.matrix(bison[,c(rep(FALSE,5), TRUE, FALSE, TRUE, TRUE)]) + strata(Strata) + cluster(Cluster), 
          data = bison, random = ~ as.matrix(bison[,c(rep(FALSE,5), TRUE, FALSE, TRUE, TRUE)]))
      
      res <- c(way1$D, way2$D, way3$D, way4$D, way5$D, way6$D, way7$D, way8$D, 
               way9$D, way10$D, way11$D, way12$D, way13$D, way14$D, way15$D, way16$D)
      expect_that(res, equals(rep(many$D, 16)))
      resname <- c(colnames(way1$D), colnames(way2$D), colnames(way3$D), colnames(way4$D), colnames(way5$D), colnames(way6$D),
          colnames(way7$D), colnames(way8$D), colnames(way9$D), colnames(way10$D), colnames(way11$D), colnames(way12$D),
          colnames(way13$D), colnames(way14$D), colnames(way15$D), colnames(way16$D))
      expect_that(resname, equals(rep(colnames(many$D), 16)))      
      
      
      # If not all covariates have random coefficient
      nar <- Ts.estim(formula = Y ~ forest + biomass + pmeadow + strata(Strata) + cluster(Cluster), data = bison,
          random = ~ forest + biomass)
      nar2 <- Ts.estim(formula = Y ~ as.matrix(bison[,c(6,8,9)]) + strata(Strata) + cluster(Cluster), data = bison,
          random = ~ as.matrix(bison[,c(6,8)]))
      expect_that(nar$D, equals(nar2$D))
      expect_that(nar$D[3,], is_equivalent_to(rep(0,3)))
      expect_that(nar$D[,3], is_equivalent_to(rep(0,3)))
      
      # type of covariance structure
      UN1d <- Ts.estim(formula = Y ~ biomass + I(biomass^2) + strata(Strata) + cluster(Cluster), 
          data = bison, random = ~ biomass + I(biomass^2))
      UN1 <- Ts.estim(formula = Y ~ biomass + I(biomass^2) + strata(Strata) + cluster(Cluster), 
          data = bison, random = ~ biomass + I(biomass^2), D="UN(1)")
      UN <- Ts.estim(formula = Y ~ biomass + I(biomass^2) + strata(Strata) + cluster(Cluster), 
          data = bison, random = ~ biomass + I(biomass^2), D="UN")
      expect_that(UN1d$D, equals(UN1$D))      
      expect_that(UN1$D[2,1]==0, is_true())      
      expect_that(UN$D[2,1]!=0, is_true())
      
      # categorical covariates
      bison_c <- cbind(bison, catego=ifelse(bison$biomass>0&bison$biomass<0.5,"medium",ifelse(bison$biomass>=0.5,"high","zero")))      
      Fit_c <- Ts.estim(formula = Y ~ catego + strata(Strata) + cluster(Cluster), data = bison_c, all.m.1=FALSE)
      expect_that(length(Fit_c$beta), equals(2))
      
      })



