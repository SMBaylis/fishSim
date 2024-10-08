% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fishSim_dev.R
\name{PoNG}
\alias{PoNG}
\title{Find a Point of No Growth (by tweaking first-year mortality)}
\usage{
PoNG(
  mateType = "flat",
  mortType = "flat",
  batchSize,
  firstBreed = 1,
  maxClutch = Inf,
  osr = c(0.5, 0.5),
  maturityCurve,
  femaleCurve,
  maxAge = Inf,
  mortRate,
  ageMort,
  stockMort,
  ageStockMort
)
}
\arguments{
\item{mateType}{the value of \code{type} used in the altMate() call. Must be one
of \code{flat}, \code{age}, or \code{ageSex}. If \code{flat}, \code{batchSize} must be provided. If
\code{age}, \code{maturityCurve} and \code{batchSize} must be provided. If \code{ageSex},
\code{femaleCurve} and \code{batchSize} must be provided. Defaults to \code{flat}.}

\item{mortType}{the value of \code{type} used in the \code{\link[=mort]{mort()}} call. Must be one of
\code{flat}, \code{age}, \code{stock}, or \code{ageStock}. If \code{flat}, \code{mortRate} must be
provided. If \code{age}, \code{ageMort} must be provided. If \code{stock}, \code{stockMort} must
be provided. If \code{ageStock}, \code{ageStockMort} must be provided. Defaults to
\code{flat}.}

\item{batchSize}{the value of \code{batchSize} used in the \code{\link[=altMate]{altMate()}} call.
Cannot be blank.}

\item{firstBreed}{the value of \code{firstBreed} used in the \code{\link[=altMate]{altMate()}} call.
Defaults to 1.}

\item{maxClutch}{the value of \code{maxClutch} used in the \code{\link[=altMate]{altMate()}} call.
Defaults to \code{Inf}. If non-\code{Inf}, \emph{effective} \code{batchSize} is estimated as the
mean of 1000000 draws from the distribution of \code{batchSize}, subsetted to
those <= \code{maxAge}.}

\item{osr}{the value of \code{osr} used in the \code{\link[=altMate]{altMate()}} call. Female
proportion is used as a multiplier on the fecundities. Defaults to \code{c(0.5, 0.5)}.}

\item{maturityCurve}{the value of \code{maturityCurve} used in the \code{\link[=altMate]{altMate()}}
call. \code{\link[=check_growthrates]{check_growthrates()}} only uses female fecundities in its estimates,
so \code{femaleCurve} is equivalent to \code{maturityCurve} in \code{\link[=check_growthrates]{check_growthrates()}},
but \code{maturityCurve} is used when mateType is \code{age}. If both mortality and
maturity are specified as vectors, they can be of different lengths. If the
maturity vector is shorter, it is "padded" to the same length as the
mortality vector by repeating the last value in the vector.}

\item{femaleCurve}{the value of \code{femaleCurve} used in the \code{\link[=altMate]{altMate()}} call.
\code{\link[=check_growthrates]{check_growthrates()}} only uses female fecundities in its estimates, so
\code{femaleCurve} is equivalent to \code{maturityCurve} in \code{\link[=check_growthrates]{check_growthrates()}}, but
\code{femaleCurve} is used when \code{mateType} is \code{ageSex}. If both mortality and
maturity are specified as vectors, they can be of different lengths. If the
maturity vector is shorter, it is "padded" to the same length as the
mortality vector by repeating the last value in the vector.}

\item{maxAge}{the value of \code{maxAge} used in the \code{mort()} call. Defaults to
\code{Inf}.}

\item{mortRate}{the value of \code{mortRate} used in the \code{mort()} call}

\item{ageMort}{the value of \code{ageMort} used in the \code{mort()} call. If both
mortality and maturity are specified as vectors, they can be of different
lengths. If the mortality vector is shorter, it is "padded" to the same
length as the maturity vector by repeating the last value in the vector.}

\item{stockMort}{the value of \code{stockMort} used in the \code{mort()} call}

\item{ageStockMort}{the value of \code{ageStockMort} used in the \code{mort()} call. If
both mortality and maturity are specified as vectors, they can be of
different lengths. If the mortality vector is shorter, it is "padded" to the
same length as the maturity vector by repeating the last value in the
vector.}
}
\description{
This is essentially a wrapper for \code{\link[=check_growthrate]{check_growthrate()}}, taking all the same
inputs. It returns a plot of projected growth-rates across all possible
first-year survival rates and the numeric first-year survival rate at which
zero population growth is expected (via \link{uniroot}, or something).  The
reasoning behind this utility is that often in biological systems, adult
survival rates and fecundities may be quite well-characterised, and
population growth rates may also be well characterised, but first-year
survival may be nearly impossible to assess. This utility allows a value of
first-year survival to be chosen such that the population size does not
change, while leaving all well-characterised adult survival parameters
unchanged. It is also useful for answering questions of the form: 'how high
would our first-year survival have to be, in order for this population to
\emph{not} be in decline?', which will surely be familiar in applied management
situations. Note that \code{\link[=PoNG]{PoNG()}} uses some brute-force methods on the back
end, so it's not terribly efficient. It takes around 5 - 10 seconds per
stock on a fairly-modern (vintage 2018) laptop.
}
\examples{
batchSize = 0.8
firstBreed = 1
mortRate = 0.2
PoNG(batchSize = batchSize, firstBreed = firstBreed, mortRate = mortRate)

mortType = "stock"
stockMort = c(0.2, 0.3, 0.5)
firstBreed = 1
batchSize = 0.9
PoNG(mortType = "stock", batchSize = batchSize, firstBreed = firstBreed, stockMort = stockMort)
 ## note that only two of the stocks return a valid PoNG - with 0.5 mortality, stock 3 cannot
 ## reach null growth with any first-year survival rate between 0 and 1.
}
\seealso{
\code{\link[=check_growthrate]{check_growthrate()}}
}
