\name{network.test}
\alias{network.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
  High-Dimensional Differential Network Inference
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Testing differential connectivity in high-dimensional networks
}
\usage{
network.test(X, est.method = "MB", test.method = "GraceI", sample.split = FALSE, 
alpha = 0.05, test.level = "node", rule = "OR")
}
\arguments{
  \item{X}{a list of standardized design matrices with the same number of columns.}
  \item{est.method}{method for network estimation: neighborhood selection ("MB").}
  \item{test.method}{method for hypothesis testing: GraceI ("GraceI"), LDPE ("LDPE"), ridge ("ridge") or SKAT ("SKAT").}
  \item{sample.split}{whether samples need to be randomly splitted for estimation and testing.}
  \item{alpha}{alpha level of type-I error rate.}
  \item{test.level}{node-wise ("node") or edge-wise ("edge").}
  \item{rule}{"AND" or "OR" rule for the estimation.}
}
\details{
This function tests whether multiple networks have the same set of edges, controlling type-I error rate.}
\value{
  An R `list' with elements:
  \item{diffnet}{a matrix or vector of differential connections.}
}
\author{
  Sen Zhao
}