\name{RemoveBatchEffects}
\alias{RemoveBatchEffects}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
RemoveBatchEffects(exp, process.info, sample.info)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{exp}{
%%     ~~Describe \code{exp} here~~
}
  \item{process.info}{
%%     ~~Describe \code{process.info} here~~
}
  \item{sample.info}{
%%     ~~Describe \code{sample.info} here~~
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
function (exp, process.info, sample.info) 
{
    ex.date <- as.factor(process.info$RNA_EXTRACT_DATE)
    sen.id <- as.factor(sample.info$SENTRIX_ID)
    sen.pos <- as.factor(sample.info$SENTRIX_POSITION)
    exp.nrm <- matrix(nrow = nrow(exp), ncol = ncol(exp), dimnames = dimnames(exp))
    coef <- matrix(nrow = nrow(exp), ncol = 3)
    for (i in 1:nrow(exp)) {
        non.na <- array(!is.na(exp[i, ]))
        if (length(which(non.na)) < 16) {
            exp.nrm[i, ] <- NA
            next
        }
        mod.1 <- lm(as.numeric(exp[i, ]) ~ sen.id + sen.pos + 
            ex.date)
        fit.1 <- summary(mod.1)
        exp.nrm[i, non.na] <- fit.1$residuals
        mod.2 <- lm(exp.nrm[i, ] ~ as.numeric(exp[i, ]))
        fit.2 <- summary(mod.2)
        coef[i, ] <- c(fit.2$coefficients[2], fit.2$r.squared, 
            fit.1$r.squared)
    }
    exp.nrm <- data.frame(exp.nrm)
    coef <- data.frame(FIT_EST = coef[, 1], FIT_R2 = coef[, 2], 
        CORRECT_R2 = coef[, 3])
    out <- list(exp.nrm, coef)
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
