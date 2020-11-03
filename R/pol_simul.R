#' Simulate data using the quadratic-normal voting model
#'
#' Generate simulated data using the quadratic-normal voting model (Poole, 2005) with ideological homogeneity among same party members and ideological distance between different party members.
#'
#' @param party \eqn{N\times 1} vector of voter party affiliation.
#' @param M Number of voting items.
#' @param partyMean \eqn{P\times 1} vector of party mean ideology for \eqn{P} parties.
#' @param partySD \eqn{P\times 1} vector for standard deviation of intra-party ideology.
#'
#' @return
#' \itemize{
#'   \item \code{votes} \eqn{N\times M} binary vote matrix with yes: 1 and no: 0 generated from the quadratic normal voting model (Poole, 2005) with Gaussian noise drawn from \eqn{\mathcal{N}(0,0.2^2)}.
#'   \item \code{ideology} Ideology of legislators used in the simulation.
#'   \item \code{a} \eqn{M\times 1} vector of ideological locations of bill proposals drawn from \eqn{\mathcal{N}(0,1)}.
#'   \item \code{q} \eqn{M\times 1} vector of ideological locations of corresponding status quos drawn from \eqn{\mathcal{N}(0,1)}.
#'  }
#' @export pol_simul
#'
#' @examples
#' ## Generate a vote matrix with 1000 votes in a two party legislature consisting of 100 legislators
#' party <- rbind(matrix(1,50,1),matrix(2,50,1))
#' vote_data <- pol_simul(party = party, M = 1000, partyMean = c(-1,1), partySD = c(.5,.5))
#'
#' @references Poole, Keith T. \emph{Spatial models of parliamentary voting}. Cambridge University Press, 2005.
#'
#' @importFrom stats rnorm

pol_simul <- function(party,M,partyMean,partySD) {

  N <- max(dim(party)[1],dim(party)[2])
  plist <- as.matrix(as.numeric(names(table(party))))
  P <- max(dim(plist)[1],dim(plist)[2])

  ideology=matrix(0,N,1)
  for (i in 1:P){
    ideology[party==plist[i],1]=rnorm(sum(party==plist[i]),mean=partyMean[i],sd=partySD[i])
  }
  a <- matrix(rnorm(M,mean=0,sd=1), 1, M)
  q <- matrix(rnorm(M,mean=0,sd=1), 1, M)

  votes <- ((-(matrix(rep(ideology,M),ncol=M)-t(matrix(rep(a,N),ncol=N)))^2 + (matrix(rep(ideology,M),ncol=M)-t(matrix(rep(q,N),ncol=N)))^2 + matrix( rnorm(N*M,mean=0,sd=1/5), N, M))>0)

  return(votes = list(votes = as.matrix(votes), ideology = ideology, a = a, q = q))
}
