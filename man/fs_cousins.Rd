\name{fs_cousins}
\alias{fs_cousins}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
fs_cousins
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
Data generated using ped-sim from either offspring of first cousin or second cousin (500 for each).
Poisson recombination model was used.
}
\usage{data("fs_cousins")}
\format{
  A data frame with 9208 observations generated using ped-sim https://github.com/williamslab/ped-sim.
  \describe{
    \item{\code{subject}}{The ID of a specific individual}
    \item{\code{chr}}{The chromosome from which the segment came from}
    \item{\code{start}}{The physical starting position of the segment}
    \item{\code{end}}{The physical end position of the segment}
    \item{\code{start_cm}}{The genetic starting position of the segment (in centimorgan)}
    \item{\code{end_cm}}{The genetic end position of the segment (in centimorgan)}
    \item{\code{length_cm}}{The genetic length of the segment (in centimorgan)}
    \item{\code{length}}{The physical length of the segment}
    \item{\code{Nm}}{Number of male ancestors}
    \item{\code{Nf}}{Number of female ancestors}
    \item{\code{male_meiosis}}{Is there an extra male meiosis? either 1 if there was an extra male meiosis, or -1 if there was an extra female meiosis}
    \item{\code{mutation}}{The number of mutations in segment}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
  https://github.com/williamslab/ped-sim
}
\examples{
data(fs_cousins)
## maybe str(fs_cousins) ; plot(fs_cousins) ...
}
\keyword{datasets}
