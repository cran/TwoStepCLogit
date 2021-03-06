\name{NEWS}
\title{News for Package 'TwoStepCLogit'}
\section{Changes in TwoStepCLogit version 1.2.5}{
  \itemize{
    \item Use of a different approximation in the calculation of the asymptotic 
          variance matrix in vcov that now gives non-zero asymptotic covariances. 
          }
}
\section{Changes in TwoStepCLogit version 1.2.4}{
  \itemize{
    \item Ts.estim() now returns the variance-covariance matrix of the regression 
          coefficients.
          }
}
\section{Changes in TwoStepCLogit version 1.2.3}{
  \itemize{
    \item Bug fix: Ts.estim() now works correctly when called inside another function.
          }
}
\section{Changes in TwoStepCLogit version 1.2.2}{
  \itemize{
    \item The rownames of the output matrix r.effects from Ts.estim() now identify
          the clusters.
    \item Test files and some files related to RStudio were removed from the source files.
    \item Bug fix: the name of the data file can now be the name of a function present
          in the search path.
    \item A meaningful error message is now returned if the value of a covariate in the
          formula remains constant within all strata of at least one cluster.
          }
}
\section{Changes in TwoStepCLogit version 1.2.1}{
  \itemize{
    \item Bug fix: Ts.estim() now works correctly with categorical covariates.
          }
}
\section{Changes in TwoStepCLogit version 1.2}{
  \itemize{
    \item The handling of errors from coxph() has been changed in order to identify 
          the cluster from which an error is generated.
    \item The warnings generated by coxph() are now stored in a new output value 
          'coxph.warn'.
    \item A file TwoStepCLogit-Ex.Rout.save has been added to the package source in 
          the subdirectory 'tests/Examples'. Now R CMD check compares the output from 
          the package examples to this file and reports any differences.
    \item The function ddim() has been modified. It nows takes as input arguments 
          'formula' and 'data'. The argument 'formula' must identify the response, 
          cluster and strata variables. It has a print method. ddim() returns 'Ystat'
          and 'Sc', which used to be returned by Ts.estim(). These elements have been
          removed from Ts.estim() output values.
    \item The print method for objects of class 'Ts.estim' now removes from the D 
          matrix the rows and columns containing only zeros before printing it.
    \item Tests have been added in the package subdirectory 'inst/tests'.
    \item The extraction of the covariate names for the output has been improved.  
    \item The package documentation has been updated according to the changes 
          mentioned here.
  }
}
\section{Changes in TwoStepCLogit version 1.1}{
  \itemize{
    \item Many validation error messages for input arguments have been added to 
          Ts.estim().    
    \item A function named ddim() has been added. It prints dimension statistics for 
          a data set processed with Ts.estim().
    \item The argument 'All.m.1' has been renamed 'all.m.1'.
    \item The output values from Ts.estim() now include 'call', the function call.
    \item The output values from Ts.estim() now include 'Ystat', the numbers of 
          observations per stratum and the sum of the responses per stratum.
    \item The output values from Ts.estim() now include 'Sc', the number of strata 
          for each cluster.
    \item The output values 'beta', 'se' and 'D' from Ts.estim() now have a 'names' 
          or 'dimnames' attribute in order to identify the covariates associated 
          with the coefficients.  
    \item A print method has been added for output of the Ts.estim() function.
    \item The input argument 'random' of Ts.estim() is now a formula.
    \item The input arguments 'Donnees' and 'Proposed.model' of Ts.estim() have,
          respectively, become 'data' and 'formula', as for many functions fitting 
          linear models in R, such as lm().
    \item The package documentation has been updated according to the changes 
          mentioned here.
  }
}
\section{Changes in TwoStepCLogit version 1.0-2}{
  \itemize{
    \item Tests written with functions from the package 'testthat' have been added in 
          the package subdirectory 'inst/tests'.
    \item The output values from Ts.estim() now include 'r.effect', the random effect 
          estimates.
    \item 'roxygen2' is now used to generate the package documentation.
    \item The examples have been changed in order to reproduce the paper results and 
          to avoid warnings.
  }
}
\section{Changes in TwoStepCLogit version 1.0-1}{
  \itemize{
    \item Bug fix: Ts.estim() works even if the cluster id numbers are not 
          consecutive numbers.
  }
}
\section{Changes in TwoStepCLogit version 1.0}{
  \itemize{
    \item First CRAN release of the TwoStepCLogit package.
  }
}
