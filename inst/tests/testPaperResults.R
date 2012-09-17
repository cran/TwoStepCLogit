context("Paper results")

test_that("Results from table 2 of the paper are reproduced", {
      Tsout <- Ts.estim(formula = Y ~ biomass + I(biomass^2) + strata(Strata) + cluster(Cluster), 
               data = bison, all.m.1=FALSE, D="UN(1)")      
      expect_that(round(Tsout$beta,2), equals(c("biomass"=8.59, "I(biomass^2)"=-7.79)))
      expect_that(round(Tsout$se,2), equals(c("biomass"=1.13, "I(biomass^2)"=2.07)))
      expect_that(round(diag(Tsout$D),1), equals(c("biomass"=13.8, "I(biomass^2)"=48.4)))
      expect_that(round(colSums(Tsout$r.effect),5), equals(c("biomass"=0, "I(biomass^2)"=0)))

    
      ddout <- ddim(formula = Y ~ strata(Strata) + cluster(Cluster), data = bison)
      expect_that(length(ddout$Sc), equals(20))
      expect_that(sum(ddout$Sc), equals(1410))      
    })


