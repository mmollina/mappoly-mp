#' HMM estimation function
#'
#' @param void internal function to be documented
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
hmm_map_reconstruction <- function(ploidy, n.mrk, n.ind, haplo, emit = NULL,
                                   rf_vec, verbose = TRUE, use_H0 = FALSE,
                                   highprec = FALSE, tol = 10e-4) {
  ## Checking capabilities
  if (verbose && !capabilities("long.double") && highprec){
    cat("This function uses high precision calculations, but your system's architecture doesn't support long double allocation ('capabilities('long.double') = FALSE'). Running in low precision mode.\n")
    highprec = FALSE
  }
  ## In case no genotypic probabilities distributions are provided
  if(is.null(emit)){
    emit <- vector("list", length(haplo))
    for(i in  1:length(haplo)){
      tempemit <- vector("list", length(haplo[[i]]))
      for(j in 1:length(haplo[[i]])){
        tempemit[[j]] <- rep(1, nrow(haplo[[i]][[j]]))
      }
      emit[[i]] <- tempemit
    }
  }
  if(highprec){
    res.temp <-
      .Call("est_haplotype_map_highprec",
            ploidy,
            n.mrk,
            n.ind,
            haplo,
            emit,
            rf_vec,
            verbose,
            tol,
            use_H0,
            PACKAGE = "mappoly2")
    return(res.temp)

  } else {
    res.temp <-
      .Call("est_haplotype_map",
            ploidy,
            n.mrk,
            n.ind,
            haplo,
            emit,
            rf_vec,
            verbose,
            tol,
            use_H0,
            PACKAGE = "mappoly2")
    return(res.temp)
  }
}
