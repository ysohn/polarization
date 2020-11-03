#' Vote profile correlation matrix converter
#'
#' Convert an \eqn{N\times M} vote matrix into an \eqn{N\times N} vote profile correlation matrix of \eqn{N} voters following Newman (2001).
#'
#' @param V \eqn{N\times M} vote matrix with yes: 1 and no: 0.
#' @return \eqn{N\times N} adjacency matrix of voter-to-voter vote profile correlation for yes votes with normalization following Newman (2001).
#' @export net_newman
#' @examples
#' ## Generate a vote matrix with 1000 votes in a two party legislature consisting of 100 members
#' party <- rbind(matrix(1,50,1),matrix(2,50,1))
#' vote_data <- pol_simul(party = party, M = 1000, partyMean = c(-1,1), partySD = c(.5,.5))
#' V <- vote_data$votes
#'
#' ## Sum up vote profile correlation matrices for yes and no votes respectively
#' A <- net_newman(V) + net_newman(!V)
#'
#' @references Newman, Mark EJ. "Scientific collaboration networks. II. Shortest paths, weighted networks, and centrality." \emph{Physical review E} 64.1 (2001): 016132.
#'
#' Porter, Mason A., et al. "A network analysis of committees in the US House of Representatives." \emph{Proceedings of the National Academy of Sciences} 102.20 (2005): 7057-7062.
#'

net_newman <- function(V){
  B <- t(as.matrix(colSums(V)))
  V <- V[,B>1]
  N <- dim(V)[1]
  A <- (V/(matrix(1,N,1)%*%(t(as.matrix(colSums(V)))-1)))%*%t(V)
  diag(A) <- 0
  return(A)
}
