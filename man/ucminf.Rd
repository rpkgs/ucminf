\name{ucminf}
\alias{ucminf}

\title{UCMINF - An algorithm for unconstrained non-linear optimization}

\description{
  The algorithm is of quasi-Newton type with BFGS updating of the inverse
  Hessian and soft line search with a trust region type monitoring of the
  input to the line search algorithm.
}

\usage{
ucminf(par, fn, gr = NULL, \dots, control = list(), hessian=0)
}


\arguments{    
  \item{par}{
    Initial estimate of minimum for \sQuote{fn}.
  }
  \item{fn}{
    Objective function to be minimized.
  }
  \item{gr}{
    Gradient of objective function. If NULL a finite difference
    approximation is used.
  }
  \item{\dots}{
    Optional arguments passed to the objective and gradient functions.
  }        
  \item{control}{
    A list of control parameters. See \sQuote{Details}.
  }
  \item{hessian}{
    Integer value:
    \describe{
      \item{0}{No hessian approximation is returned.}
      \item{1}{Returns a numerical approximation of the Hessian
	using \sQuote{hessian} in the package \sQuote{numDeriv}. }
      \item{2}{Returns final approximation of the inverse Hessian based on the series of BFGS
	updates during optimization. }
      \item{3}{Same at 2, but will also return the Hessian (the inverse of 2). }
    }
    If a \code{TRUE} or \code{FALSE} value is given it will switch
  between option 1 or 0.
  }
}

\details{
  The implementation of UCMINF is based on the original Fortran source
  code. The algorithm is documented in (Nielsen, 2000) (see References
  below) with a description of the algorithm and a comparison to the
  Fortran subroutine \code{MINF} and the Matlab function
  \code{fminunc}.

  The interface in R is designed so that it is very easy to switch
  between using \code{ucminf} and \code{optim}. The 
  arguments \sQuote{par}, \sQuote{fn}, \sQuote{gr}, and \sQuote{hessian}
  are all the same (with a few extra options for \sQuote{hessian} in
  \code{ucminf}). The difference is that there is no \sQuote{method} 
  argument in \code{ucminf} and that some of the components in the
  \sQuote{control} argument are different due to differences in the
  algorithms. 

  The algorithm can be given an initial estimate of the Hessian for the
  optimization and it is possible to get the final approximation of the
  Hessian based on the series of BFGS updates. This extra functionality
  may be useful for optimization in a series of related problems.

  The functions \sQuote{fn} and \sQuote{gr} can return \code{Inf} or \code{NaN}
  if the functions cannot be evaluated at the supplied value, but the
  functions must be computable at the initial value. Note, the functions
  are not allowed to return \code{NA}.
  
  The \sQuote{control} argument is a list that can supply any of the
  following components:
  
  \describe{
    \item{trace}{
      If trace is positive then tracing information is printed for each iteration.}
    \item{grtol}{
      The algorithm stops when
      \eqn{||F'(x)||_\infty \leq }{||F'(x)||_inf <=} grtol, that
      is when the largest absolute value of the gradient is less than
      grtol. Default value is \code{grtol = 1e-6}. }
    \item{xtol}{
      The algorithm stops when \eqn{||x-x_p||_2 \leq
	\textrm{xtol}\cdot(\textrm{xtol} + 
      ||x||_2)}{||x-x_p||_2 <=
	xtol*(xtol + ||x||_2)}, where \eqn{x_p} and \eqn{x} are the
      previous and current estimate of the minimizer. Thus the algorithm
      stops when the last relative step length is
      sufficiently small. Default value is \code{xtol = 1e-12}.}
    \item{stepmax}{
      Initial maximal allowed step length (radius of
      trust-region). The value is updated during the
      optimization. Default value is \code{stepmax = 1}.}
    \item{maxeval}{
      The maximum number of function evaluations. A function evaluation
      is counted as one evaluation of the objective function and of the
      gradient function. Default value is \code{maxeval = 500}.  } 
    \item{grad}{
      Either \sQuote{forward} or \sQuote{central}. Controls
      the type of finite difference approximation to be used for the
      gradient if no gradient function is given in the input argument
      \sQuote{gr}. Default value is \code{grad = 'forward'}.}
    \item{gradstep}{
      Vector of length 2. The step length in finite
      difference approximation for the gradient. Step length is
      \eqn{|x_i|\cdot\textrm{gradstep[1]+gradstep[2]}}{|x_i|*gradstep[1]+gradstep[2]}.
      Default value is \code{gradstep = c( 1e-6, 1e-8 )}.  } 
    \item{invhessian.lt}{
      A vector with an initial approximation to the lower triangle of the
      inverse Hessian. If not given, the inverse Hessian is initialized
      as the identity matrix. If \code{H0} is the initial hessian matrix then
      the lower triangle of the inverse of \code{H0} can be found as
      \code{invhessian.lt = solve(H0)[lower.tri(H0,diag=TRUE)]}.}
    \item{method.args}{
      Used when \code{hessian = 1}. Arguments for
      hessian approximation using the \sQuote{numDeriv} package, see
      \code{\link[numDeriv]{hessian}} for more details.}
  }
}

\value{
  \item{par}{Computed minimizer.}
  \item{value}{Objective function value at computed minimizer.}
  \item{convergence}{Flag for reason of termination:
    \describe{
      \item{1}{Stopped by small gradient (grtol).}
      \item{2}{Stopped by small step (xtol).}
      \item{3}{Stopped by function evaluation limit (maxeval).}
      \item{4}{Stopped by zero step from line search}
      \item{-2}{Computation did not start: length(par) = 0.}
      \item{-4}{Computation did not start: stepmax is too small.}
      \item{-5}{Computation did not start: grtol or xtol <= 0.}
      \item{-6}{Computation did not start: maxeval <= 0.}
      \item{-7}{Computation did not start: given Hessian not pos. definite.}
    }
  }
  \item{message}{String with reason of termination.}
  \item{hessian, invhessian}{Estimate of (inv.) Hessian at computed minimizer. The type of
    estimate is given by the input argument \sQuote{hessian}.}
  \item{invhessian.lt}{The lower triangle of the final approximation to the
    inverse Hessian based on the series of BFGS	updates during optimization.}
  \item{info}{Information about the search:
    \describe{
      \item{maxgradient}{\eqn{||F'(x)||_\infty}{||F'(x)||_inf}, the
	largest element in the absolute 
	value of the gradient at the computed minimizer.}
      \item{laststep}{Length of last step.}
      \item{stepmax}{Final maximal allowed step length.}
      \item{neval}{Number of calls to both objective and gradient function.}
    }
  }
}

\author{

  \sQuote{UCMINF} algorithm design and Fortran code by Hans Bruun Nielsen.

  Implementation in \R by Stig B. Mortensen, \email{sbm@imm.dtu.dk}.
    
}

\seealso{\code{\link[stats]{optim}}, \code{\link[stats]{nlminb}}, \code{\link[stats]{nlm}}.}

\references{

  Nielsen, H. B. (2000) \sQuote{UCMINF - An Algorithm For Unconstrained,
  Nonlinear Optimization}, Report IMM-REP-2000-18,
  Department of Mathematical Modelling, Technical University of
  Denmark. \url{http://www2.imm.dtu.dk/~hbn/publ/TR0019.ps}

  The original Fortran
  source code can be found at
  \url{http://www2.imm.dtu.dk/~hbn/Software/ucminf.f}. The code has been
  slightly modified in this package to be suitable for use with R.

  The general structure of the implementation in \R is based on the
  package \sQuote{FortranCallsR} by Diethelm Wuertz. 
  
}




\examples{
## Rosenbrock Banana function
fR = function(x) { (1-x[1])^2+100*(x[2]-x[1]^2)^2 }
gR = function(x) { c( -400*x[1]*(x[2]-x[1]*x[1]) - 2*(1-x[1]),
                       200*(x[2]-x[1]*x[1])) }
#  Find minimum
ucminf(par = c(2,.5), fn = fR, gr=gR)

#  Compare hessian approximations
ucminf(par = c(2,.5), fn = fR, gr=gR, hessian=1)$hessian
ucminf(par = c(2,.5), fn = fR, gr=gR, hessian=3)$hessian

#  Compare run times with optim's BFGS method
#  (chosen convergence criteria result in similar accuracy)
system.time( for(i in 1:500)
               ucminf(par = c(2,0.5), fn = fR, gr=gR)
           )
system.time( for(i in 1:500)
               optim(par = c(2,0.5), fn = fR, gr=gR,method='BFGS')
           )


## Quadratic function
fQ = function(x) { sum((4*x-1)^2) }
gQ = function(x) { 32*x -8}

#  Find minimum with too small stepmax and print trace
ucminf(par = c(20.5,20.0), fn = fQ, gr = gQ, 
       control=list(stepmax=1,trace=TRUE))

#  The same again with a larger stepmax
ucminf(par = c(20.5,20.0), fn = fQ, gr = gQ,
       control=list(stepmax=100,trace=TRUE))
}

\keyword{optimize}

