\name{ellipse}
\alias{circle}
\alias{ellipse}
\alias{norm.ellipse}
\alias{panel.ellipse}
\title{draw circles and ellipses}
\usage{
circle(x=0, y=0, radius=1, sx, sy, \dots)
ellipse(x=0, y=0, sx=1, sy=1, rho=0, fine=0.005, center=0.1,
  show=TRUE, add=TRUE, xlab, ylab, bty, xaxt, yaxt, type, lty, pty, \dots)
norm.ellipse(x, y, level=0.05,
  radius=sqrt(qchisq(1-level,2)), \dots)
panel.ellipse(x=0, y=0, sx=1, sy=1, rho=0, fine=0.005, center=0.1,
  xlab, ylab, bty, xaxt, yaxt, type, lty, pty, \dots)
}
\description{}
\arguments{
  \item{x,y}{center of circle or ellipse before
    orienting}\item{radius}{length of radius (positive)}
  \item{sx,sy,rho}{standardizing units for x,y (positive)and
    correlation between x and y (between -1 and 1)}
  \item{fine}{fineness of arcs for drawing (between 0 and 1)}
  \item{center}{if positive, fraction of radius for main axes}
  \item{show}{plot ellipse or circle if TRUE}
  \item{add}{add to existing plot if TRUE}
  \item{level}{probability level (default 0.05)}
  \item{xlab,ylab,xaxt,yaxt,type,lty,pty,bty}{plot parameters}
  \item{\dots}{other plot parameters (many are set if add=F)}
}
\value{list of dataframes containing ellipse or circle arc(s)}
\examples{
# circle with unit radius
circle()
# return coordinates of arc for ellipse
ellipse(sx=2,sy=3,rho=.5,show=FALSE)
# ellipse rotated 45 degrees, centered at (1,2)
ellipse(xo=1,yo=2,sx=.25,sy=.5,orient=pi/2)
# ellipse at (1,2), but then rotated 45 degrees
ellipse(x=1,y=2,sx=.25,sy=.5,orient=pi/2)
\dontrun{
# 95\% contour for bivariate normal from data (x,y)
norm.ellipse(x,y)
}
}
\keyword{hplot}
% Converted by  version .
