
# MIRT4FC

<!-- badges: start -->
<!-- badges: end -->

The goal of MIRT4FC is to fit Two-Parameter Logistic Item Response Theory(2PL-IRT) models in R by iStEM alogrithm. Functionality for extracting results, making predictions, and simulating data is provided as well. 

## Installation

You can install the development version of MIRT4FC like so:

``` r
library(MIRT4FC)
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(MIRT4FC)
D <- 6
nitem.per.dim <- 10
nblock <- D * nitem.per.dim / 3
set.seed(123456)
item.par <- data.frame(a=seq_len(D*nitem.per.dim))
item.par <- within(item.par,{
 a <- runif(D*nitem.per.dim,0.7,3)
 b <- rnorm(D*nitem.per.dim)
 d <- a*b
})
BID <- data.frame(Block=rep(1:nblock,each=3),
                 Item=rep(1:3,nblock),
                 Dim=c(combn(D,3)[,sample(choose(D,3),nblock,replace = TRUE)]))
item.par$d <- c(t(aggregate(item.par$d,by=list(BID$Block),function(x)x-mean(x))[,-1]))
N <- 1000
v <- matrix(0.5,D,D)
diag(v) <- 1
theta <- mvnfast::rmvn(N,seq(-1,1,length.out = D),sigma = v)
Y <- data.sim(item.par,theta,BID,blocksize=3,res='rank')
x <- StEM(Y,BID,maxitr = 100,blocksize=3,res='rank',fix.sigma = TRUE,cores=1)
```

