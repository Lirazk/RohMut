\name{fs_cousins}
\alias{fs_cousins}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
fs_cousins
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
Data generated using ped-sim from either offspring of first cousin or second cousin (500 for first cousin and 478 for second).
Poisson recombination model was used.
}
\usage{data("fs_cousins")}
\format{
  A data frame with 9208 observations generated using ped-sim \url{https://github.com/williamslab/ped-sim}.
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
  \url{https://github.com/williamslab/ped-sim}
  
  Caballero, M., Seidman, D. N., Qiao, Y., Sannerud, J., Dyer, T. D., Lehman, D. M., Curran, J. E., Duggirala, R., Blangero, J., Carmi, S., &amp; Williams, A. L. (2019). Crossover interference and sex-specific genetic maps shape identical by descent sharing in close relatives. PLOS Genetics, 15(12). \cr\url{https://doi.org/10.1371/journal.pgen.1007979}
}
\examples{
data(fs_cousins)
## maybe str(fs_cousins) ; plot(fs_cousins) ...
}
\keyword{datasets}
