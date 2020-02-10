\name{found.spill.CIPHE}
\alias{found.spill.CIPHE}
\title{Found the spill keyword inside flowFrome}
\usage{
found.spill.CIPHE(flowFrame)
}
\description{
Found the keywords name of a square matrixe with colnames of this matrice are containes in flowFrame colnames expression matrice.
The function retourne the name of this keywords and suppose to be the spill matrice. The function return also the number of dimension in square matrice.
}
\examples{
fcs <- read.FCS.CIPHE("path.fcs")
spill <- found.spill.CIPHE(fcs)
fcs.com <- compensate.CIPHE(fcs, spill)
}
