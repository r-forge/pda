\name{uncollide}
\alias{uncollide}
\title{spread out values to avoid collision in plotting}
\usage{
uncollide( x, small )
}
\description{}
\arguments{
  \item{x}{numeric vector}
  \item{small}{minimum separation of values}
}
\value{vector of x's with slight adjustment to prevent overlap but stay
    within range}
  \examples{
\dontrun{
text( x, uncollide(y), state )
}
}
% Converted by  version .
\keyword{manip}
