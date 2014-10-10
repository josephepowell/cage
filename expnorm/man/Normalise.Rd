\name{Normalise}
\alias{Normalise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Normalise Gene Expression data}
\description{
Perform end-to-end normalisation of gene expression data, through log2 
transformation, quantile normalisation, inverse normal transformation, 
and principal components analysis. By default, all operations are 
performed.
}
\usage{
Normalise(x, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
matrix or data.frame of expression data to be normalised.
}
  \item{log}{
logical, should log2 transformation be performed.
}
  \item{quantile}{
logical, should quantile normalisation be performed.
}
  \item{pca}{
logical, should the data be transformed have range in 0 to 1.
}
  \item{inv.norm}{
logical, should batch effects be removed by principal components analysis.
}
  \item{n.pcs}{
numeric, quantity of principal components to fit.
}
  \item{neg.rm}{
logical, should negative values be removed during log2 transformation.
}
  \item{row.names}{
logical, does the matrix contain a character column of probe IDs.
}
}
\details{
Default usage performs all transformation and normalisation operations in order.

Inverse normal transformation is rank-based, and performed according to 
Blom (1958):

    qnorm((r - c) / (n - 2c + 1))

where \code{qnorm()} is the standard normal quartile function, r is the variable's
rank, c is the Blom constant (3 / 8), and n is sample size.
}
\value{
A \code{data.frame} of gene expression data, normalised and transformed 
according to supplied arguments.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Alex Holloway}
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
function (x, log = TRUE, quantile = TRUE, pca = TRUE, inv.norm = TRUE, 
    n.pcs = 25, neg.rm = TRUE, row.names = FALSE) 
{
    if (log) {
        x <- LogNormalise(x, neg.rm, row.names)
    }
    if (quantile) {
        x <- QuantileNormalise(x, row.names)
    }
    if (pca) {
        x <- x[apply(x, 1, function(x) all(is.na(x))), ]
        pc <- ppca(t(x), n.pcs)
        x <- CorrectByPca(x, pc)
    }
    if (inv.norm) {
        x <- t(apply(x, 1, InverseNormal))
    }
    return(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
