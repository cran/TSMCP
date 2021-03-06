% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cpvnts.R
\name{cpvnts}
\alias{cpvnts}
\title{Multiple change point detection and variable selection for nonstationary
time series models.}
\usage{
cpvnts(Y, method = c("mcp", "scad"), p, lam.d)
}
\arguments{
\item{Y}{an autoregressive time series.}

\item{method}{the method to be used by mcp or scad.
See \code{\link[plus]{plus}} in R packages \pkg{plus} for details.}

\item{p}{an upper bound of autoregressive orders.}

\item{lam.d}{an parameter of the refining procedure. Suggest lam.d=0.05 if
the length of Y is smaller than 5000 and lam.d=0.01 otherwise.}
}
\value{
pluscpv returns an object of class "cpvnts".
An object of class "cpvnts" is a list containing the following components:

\item{variable.selection}{a matrix of variable selection.
Element of the	matrix is 1 if the variable is selected and 0	otherwise.}

\item{change.points}{estimators of change points.}
}
\description{
This function provides a novel and fast methodology for simultaneous multiple
structural break estimation and variable selection for nonstationary time series
models by using ASCAD or AMCP to estimate the chang points proposed by Jin, Shi and Wu (2011).
}
\examples{

### Example 1: No change point with sample size 1000###

n <- 1000
Y <- rnorm(n)
cp1 <- c(3, 1000)
a0 <- c(0.7)
a1 <- c(1.3)
a2 <- c(-0.8)

for (j in 2:2) {
  for (i in (cp1[j-1] + 1):cp1[j]) {
    Y[i] <- a0[j-1] + a1[j-1] * Y[i-1] + a2[j-1] * Y[i-2] + rnorm(1)
  }
}

ts.plot(Y)

# Y(t)=0.7+1.3Y(t-1)-0.8Y(t-2)+e(t)


mcp.cp <- cpvnts(Y, "mcp", 4, 0.05)
mcp.cp  #result of AMCP

scad.cp <- cpvnts(Y, "scad", 4, 0.05)
scad.cp  #result of ASCAD


### Example 2(Davis et al (2006)):Two change points with sample size
### 1024##
n <- 1024
Y <- rnorm(n)
cp1 <- c(3, 512, 769, 1024)
a1 <- c(0.9, 1.69, 1.32)
a2 <- c(0, -0.81, -0.81)

for (j in 2:4) {
  for (i in (cp1[j - 1] + 1):cp1[j]) {
    Y[i] <- a1[j - 1] * Y[i - 1] + a2[j - 1] * Y[i - 2] + rnorm(1)
  }
}

ts.plot(Y)

# Y(t)=0.9Y(t-1)+e(t),if t<512 Y(t)=1.69Y(t-1)-0.81Y(t-2)+e(t), if
# 512<t<769 Y(t)=1.32Y(t-1)-0.81Y(t-2)+e(t), if 769<t<1024


mcp.cp <- cpvnts(Y, "mcp", 4, 0.05)
mcp.cp  #result of AMCP

scad.cp <- cpvnts(Y, "scad", 4, 0.05)
scad.cp  #result of ASCAD


### Example 3: Six change points with sample size 10000###
n <- 10000
Y <- rnorm(n)
cp1 <- c(3, 1427, 3084, 4394, 5913, 7422, 8804, 10000)

a1 <- c(1.58, 1.12, 1.61, 1.24, 1.53, 1.32, 1.69)
a2 <- c(-0.79, -0.68, -0.75, -0.66, -0.64, -0.81, -0.81)

for (j in 2:length(cp1)) {
  for (i in (cp1[j - 1] + 1):cp1[j]) {
    Y[i] <- a1[j - 1] * Y[i - 1] + a2[j - 1] * Y[i - 2] + rnorm(1)
  }
}



ts.plot(Y)


mcp.cp <- cpvnts(Y, "mcp", 4, 0.01)
mcp.cp  #result of AMCP

scad.cp <- cpvnts(Y, "scad", 4, 0.01)
scad.cp  #result of ASCAD


}
\references{
{Jin, B., Shi, X., and Wu, Y. (2013). A novel and fast methodology for
 simultaneous multiple structural break estimation and variable selection for
 nonstationary time series models. \emph{Statistics and Computing}, 23(2), 221-231.}
}
\seealso{
plus
}
