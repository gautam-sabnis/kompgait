% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/komp_lmer_dispersion.R
\name{komp_lmer_dispersion}
\alias{komp_lmer_dispersion}
\title{Fit a LMM to linear gait traits}
\usage{
komp_lmer_dispersion(data, CtrlStrain = "C57BL/6NJ", model)
}
\arguments{
\item{data}{data_per_animal dataframe}

\item{CtrlStrain}{specify the control strain}

\item{model}{specify the model}
}
\value{
a list containing the results of the model
}
\description{
Fit a linear mixed model to identify significant differences in
trait's dispersion value (sd) between mutant and control strain
}
\examples{
komp_lmer_mean(data = data_per_animal_disp_line, CtrlStrain = 
"C57BL/6NJ", model = "M3")
}
