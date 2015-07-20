\name{Nuclei}
\alias{Nuclei}
\docType{data}
\title{
	Brain activity in 14 brain regions for 71 individuals
}
\description{
The activity in 14 brain nuclei were measured in rats that were in one of four groups: 
1) Non-stressed, Control 
2) Stressed, Control
3) Non-stressed, Vinclozolin
4) Stressed, Vinclozolin
  
}
\usage{
data(Nuclei)
}

\format{
  A numeric matrix with 71 individuals as rows and the activity of 14 brain nuclei as columns.  NAs indicate missing data.
}

\details{
Two different cohorts of male rats of the F3 generation of Vinclozolin (Vinclozolin-Lineage) and Vehicle Control (Control-Lineage) Lineages produced at Washington State University are shipped to the University of Texas on the day after weaning. Rats are randomly pair-housed (one Control-Lineage and one Vinclozolin-Lineage animal) and remain in these dyads throughout the duration of the study. Half of the dyads are randomly chosen to receive chronic restraint stress (CRS) treatment for 6 hours daily for 21 consecutive days commencing 1 hr after lights off.  Activity in 14 brain nuclei were measured at the end of the study.
}

\source{
The data are provided courtesy of David Crews at the University of Texas at Austin.
}

\references{
Crews, D, R Gillette, SV Scarpino, M Manikkam, MI Savenkova, MK Skinner. 2012.
Epigenetic Transgenerational Alterations to Stress Response in Brain Gene Networks and 
Behavior. Proc. Natl. Acad. Sci. USA. 109 (23). 9143 - 9148.
}

\examples{
data(Nuclei)
}
