\name{Running}
\alias{Running}
\title{Running and Oxygen Consumption}
\usage{
data(Running)
}
\format{
    Running data frame with 24 observations on 13 variables.
    \tabular{rlll}{
	[,1] \tab subject \tab factor  \tab subject identifier \cr
	[,2] \tab order \tab ordered \tab order of treatment \cr
	[,3] \tab speed \tab factor \tab running speed  \cr
	[,4] \tab stride  \tab numeric \tab  stride frequency \cr
	[,5] \tab mass  \tab numeric \tab attached mass in g \cr
	[,6] \tab vo2 \tab numeric \tab oxygen consumption  \cr
 	[,7] \tab tm  \tab numeric \tab time \cr
 	[,8] \tab etot  \tab numeric \tab etot \cr
 	[,9] \tab cn \tab numeric \tab cn \cr
 	[,10] \tab ar  \tab numeric \tab ar \cr
 	[,11] \tab ecn  \tab numeric \tab ecn \cr
 	[,12] \tab ahm  \tab numeric \tab ahm \cr
 	[,13] \tab vhm  \tab numeric \tab vhm
    }
}
\description{This study examined how the metabolic cost of
locomotion varied with speed, stride frequency and body mass.  Cost
was determined by measuring oxygen consumption (\code{vo2}),
analyzing the oxygen content in air inhaled and exhaled by the
subjects through a mask.  The rate of oxgen consumption was measured
for three \code{subject}s locomoting at all combinations of low and
high levels of the three factors --- running \code{speed}, \code{stride}
frequency, and \code{mass} distribution in the leg.  The first factor
was set using a treadmill, the second by synchronizing the subjects'
pace with a metronome, and the third by varying the positions of
weights strapped onto the legs or waist.  In addition to these three
design factors, other variables were measured by filming each trial
with a high-speed motion camera.  Ignore these other measurements for
this problem.

The \code{order} of the eight combinations was randomly assigned for
each test session.  For each test combination, \code{subjects}
continued running until their rate of oxygen consumption levelled off,
between eight and 15 minutes into the run, signalling that
steady state had been reached.  Subjects were given daily sessions
for an initial training period, during which the rate of oxygen
consumption at each test combination decreased to a relatively
constant level.  \code{Subject}s were considered trained when their
oxygen consumption for all test combinations were consistent and
repeatable between test sessions.  All subjects achieved the trained
state within two weeks.  The data used in this study were collected in
single sessions from trained \code{subject}s.  However, there still
might be some concern about incomplete training causing systematic
differences in consumption during the session.
}
\source{MJ Myers}
\references{Myers MJ, Steudel K and White SC (1993)
	`Uncoupling the correlates of locomotor costs: a factorial
	approach',
	\emph{J Experimental Zoology 265}, 211-223.}
\examples{
data( Running )

Running$code <- interaction( Running$subject, floor( Running$speed ),
  factor( Running$stride ), factor( Running$mass ))

Running.fit <- lm( vo2 ~ subject+speed+stride+mass, Running )
Running.intfit <- lm( vo2 ~ speed*stride*mass, Running )

print( xyplot( vo2 ~ order, Running, groups = subject,
  type = "p", pch = levels( Running$code ), cex = 2,
  xlab = "order of treatment", ylab = "VO2",
  main = "Running Problem" ),
  more = TRUE, split = c(1,1,2,1) )

Running$pred = predict( Running.fit )
print( xyplot( pred ~ order, Running, groups = subject,
  type = "p", pch = levels( Running$code ), cex = 2,
  xlab = "order of treatment", ylab = "VO2",
  main = "(sub.speed.stride.mass)" ),
  split = c(2,1,2,1) )
} 
% Converted by  version .
\keyword{datasets}
