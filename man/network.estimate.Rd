\name{network.estimate}
\alias{network.estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
  High-Dimensional Network Estimation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Estimate multiple conditional dependency networks and their intersections.
}
\usage{
  network.estimate(X, lambda = NULL, method = "MB", rule = "OR")
}
\arguments{
  \item{X}{
	  a list of standardized design matrices with the same number of columns.
  }
  \item{method}{
  	method for network estimation: neighborhood selection ("MB").
  }
  \item{rule}{
	  "AND" or "OR" rule for the estimation.
  }
  \item{lambda}{
	  a vector of tuning parameters.
  }
}
\details{
  Estimate multiple conditional dependency networks and their intersections
}
\value{
  An R `list' with elements:
  \item{est}{a list of estimated adjacency matrices.}
  \item{common}{adjacency matrix of the common network.}
}
\author{
  Sen Zhao
}
