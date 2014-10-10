\name{CorrectByPca}
\alias{CorrectByPca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CorrectByPca(exp, pca)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{exp}{
%%     ~~Describe \code{exp} here~~
}
  \item{pca}{
%%     ~~Describe \code{pca} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (exp, pca) 
{
    stopifnot(require("pcaMethods"))
    exp.nrm <- matrix(nrow = nrow(exp), ncol = ncol(exp), dimnames = dimnames(exp))
    for (i in 1:nrow(exp)) {
        non.na <- array(!is.na(exp[i, ]))
        if (length(which(non.na)) < 16) {
            exp.nrm[i, ] <- NA
            next
        }
        mod <- lm(as.numeric(exp[i, ]) ~ scores(pca)[, 1:24])
        fit <- summary(mod)
        exp.nrm[i, non.na] <- fit$residuals
    }
    exp.nrm <- data.frame(exp.nrm)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
