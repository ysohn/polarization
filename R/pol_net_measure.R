#' Network-based polarization measure
#'
#' Polarization measure calculation from an \eqn{N\times N} vote profile correlation matrix.
#'
#' @param A \eqn{N\times N} symmetric weighted adjacency matrix of vote profile correlation.
#' @param party \eqn{N\times 1} vector of voter party affiliation.
#' @param method Polarization measure.
#' \itemize{
#'   \item \code{QNewman} Modularity value calculated following Newman (2006).
#'   \item \code{NCut} \eqn{-1\times} Normalized cut calculated following Shi and Malik (2000).
#' }
#' @return Polarization measure.
#' @export pol_net_measure
#' @examples
#' ## Generate a vote matrix with 1000 votes in a two party legislature
#' party <- rbind(matrix(1,50,1),matrix(2,50,1))
#' vote_data <- pol_simul(party = party, M = 1000, partyMean = c(-1,1), partySD = c(.5,.5))
#' V <- vote_data$votes
#'
#' ## Sum up vote profile correlation matrices for yes and no votes respectively
#' A <- net_newman(V)+net_newman(!V)
#'
#' ## Calculate Modularity and - Normalized cut for matrix A
#' pol_value <- matrix(0,2,1)
#' for (i in 1:2){
#'   pol_value[i,1] <- pol_net_measure(A, party, i)
#' }
#'
#' @references Newman, Mark EJ. "Modularity and community structure in networks." \emph{Proceedings of the national academy of sciences} 103.23 (2006): 8577-8582.
#'
#' Shi, Jianbo, and Jitendra Malik. "Normalized cuts and image segmentation." \emph{IEEE Transactions on pattern analysis and machine intelligence} 22.8 (2000): 888-905.

pol_net_measure <- function(A, party, method) {
  switch(method,
         QNewman={
           print('QNewman')
           N <- dim(A)[1]
           K <- as.matrix(rowSums(A))
           m <- sum(A)
           B <- A-(K%*%t(K))/m
           s <- matrix(rep(party,N),ncol=N)
           Q <- (!(s-t(s)))*B/m
           Q=sum(Q)
           return(Q)
         },
         NCut={
           print('-NCut')
           N <- dim(A)[1]
           plist <- as.matrix(as.numeric(names(table(party))))
           P <- max(dim(plist)[1],dim(plist)[2])
           d <- l <- n <- 0
           for (i in 1:P){
             pp <- party==plist[i]
             d[i]=sum(sum(A[pp,]))
             l[i]=sum(A[pp,pp])
             n[i]=sum(party==plist[i])
           }
           c=-sum((d-l)/d)
           return(c)
         },
         warning("no match\n")
  )
}
