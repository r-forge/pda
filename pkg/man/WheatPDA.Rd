\name{WheatPDA}
\alias{WheatPDA}
\title{Wheat virus in plant sap}
\usage{
data(WheatPDA)
data(BactRoom)
}
\format{
    Wheat data frame with 192 observations on 8 variables.
    \tabular{rlll}{
	[,1] \tab id \tab  factor \tab  bird identifier \cr
	[,2] \tab bact \tab  factor \tab  wheat or control \cr
	[,3] \tab temp \tab  factor \tab  room temperature \cr
	[,4] \tab inoc \tab  factor \tab when innoculated\cr
	[,5] \tab bill \tab  numeric \tab bill length\cr
	[,6] \tab leg \tab   numeric \tab leg length\cr
	[,7] \tab code \tab  factor \tab plot code\cr
	[,8] \tab trt \tab   trt \tab bact*trt combination
    }
    BactRoom data frame with 2685 observations on 8 variables.
    \tabular{rlll}{
	[,1] \tab ploidy \tab ordered \tab ploidy number \cr
	[,2] \tab species \tab  factor \tab species name \cr
	[,3] \tab cross \tab  factor \tab cross identifier\cr
	[,4] \tab accession \tab factor \tab accession number\cr
	[,5] \tab plant \tab factor \tab plant number\cr
	[,6] \tab tr1 \tab  numeric \tab trial 3 yield\cr
	[,7] \tab tr2 \tab  numeric \tab trial 3 yield\cr
	[,8] \tab tr3 \tab  numeric \tab trial 3 yield
    }
}
\description{An agronomist (Katsiotis 1993) studied the effect
of \code{ploidy} (number of copies of DNA chromosomes -- humans have 2)
on the yield in several species of alfalfa.  The \code{ploidy} ranged
from 2 to 8, with one to five \code{species} per \code{ploidy}.  Within
each species several matings (known as \code{cross}es) were arranged.
Each mating produced about five \code{accession}s.  Roughly 15
offspring from each accession were planted and later harvested to
measure yield.  This completely nested design has a mix of fixed
effects (\code{ploidy, species}) and random effects (\code{cross}es,
\code{accession}s, \code{offspring}).}
\source{A Katsiotis,
  UW-Madison Agronomy Department
  (\url{http://agronomy.wisc.edu})
}
\references{Katsiotis A (1993)
	`The determination and use of pollen grain size in four
	ploidy levels and the cytogenetics of tetraploid-octoploid
	hybrids in search of \emph{2n} gametes of \emph{Avena}',
	PhD Dissertation, Department of Agronomy, University of Wisconsin-Madison.}
\examples{
data( WheatPDA )

# H:22.3 Wheat BLUEs and BLUPs
class( WheatPDA$ploidy ) = factor( WheatPDA$ploidy, ordered = FALSE )
WheatPDA$sp = interaction( WheatPDA$ploidy, WheatPDA$species)
WheatPDA$cr = interaction( WheatPDA$sp, WheatPDA$cross)
WheatPDA$ac = interaction( WheatPDA$cr, WheatPDA$accession)

\dontrun{
## this does not work properly: singular X'X matrix
library( lme4 )
Wheat.lme = lmer( tr1 ~ ploidy + sp + (1|cr) + (1|ac),
  data = WheatPDA, na.action = na.omit )
Wheat.lme = lmer( tr1 ~ ploidy + sp + (1|cross),
  data = WheatPDA, na.action = na.omit )
}
}
% Converted by  version .
\keyword{datasets}
