\name{LogTransform}
\alias{LogTransform}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Log Transform
}
\description{
Perform $\log_2$ transformation on numerical data, ommitting NA values.
}
\usage{
LogTransform(x, neg.rm = TRUE, row.names = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
data to transform.
}
  \item{neg.rm}{
logical, should non-positive values be set to NA instead of an arbitrary small value.
}
  \item{row.names}{
logical, does x contain a column of character values.
}
}
\details{
If x contains a column of character values (e.g. probe IDs), the row.names arg 
should be used to assign them to the data.frame prior to stripping them for
processing.
}
\value{
data.frame of $\log_2$ transformed numeric values.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Alex Holloway
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
function (x, neg.rm = TRUE, row.names = FALSE) 
{
    if (row.names) {
        rownames(x) <- x[, 1]
        x <- x[, -1]
    }
    if (neg.rm) {
        x[x <= 0] <- NA
    }
    else {
        x[x <= 0] <- 10^(-10)
    }
    x[!is.na(x)] <- log2(x[!is.na(x)])
    x[x < 0] <- 0
    return(data.frame(x))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
