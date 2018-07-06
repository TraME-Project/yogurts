\name{msa}

\alias{msa}

\title{MSA}

\description{
    Auction method.
}

\usage{
msa(obs_choice_prob,eps_mat,tol_auction)
}

\arguments{
	\item{obs_choice_prob}{
	Observed choice probabilities.}
	\item{eps_mat}{
	Matrix of random terms.}
    \item{tol_auction}{
    Auction tolerance.
    }
}

\details{
    See the paper.
}

\value{
Returns upper and lower bounds on the delta parameters.
}

\author{Odran Bonnet, Alfred Galichon, Keith O'Hara, and Matt Shum}

\examples{
\dontrun{
msa(s_j,eps_ij,tol_auction)
}
}

\keyword{msa}
\concept{msa}