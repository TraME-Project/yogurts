\name{ipfp}

\alias{ipfp}

\title{IPFP}

\description{
    IPFP.
}

\usage{
ipfp(obs_choice_prob,eps_mat,ipfp_temper=1.0,err_tol=1E-06,max_iter=1000)
}

\arguments{
	\item{obs_choice_prob}{
	Observed choice probabilities.}
	\item{eps_mat}{
	Matrix of random terms.}
    \item{ipfp_temper}{
	Temperature parameter.}
    \item{err_tol}{
    Error tolerance.}
    \item{max_iter}{
    Maximum number of iterations.}
}

\details{
    See the paper.
}

\value{
Returns delta parameters.
}

\author{Odran Bonnet, Alfred Galichon, Keith O'Hara, and Matt Shum}

\examples{
\dontrun{
ipfp(s_j,eps_ij,1.0)
}
}

\keyword{ipfp}
\concept{ipfp}