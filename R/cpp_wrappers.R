#' HMM estimation function
#'
#' @param void internal function to be documented
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
hmm_map_reconstruction <- function(input.obj,
                                   verbose = TRUE,
                                   use_H0 = FALSE,
                                   tol = 10e-4) {

  ## Checking capabilities
  if (verbose && !capabilities("long.double")){
    cat("This function uses high precision calculations, but your system's architecture doesn't support long double allocation ('capabilities('long.double') = FALSE'). Running in low precision mode.\n")
    highprec = FALSE
  }
  ploidy1_id <- input.obj$ploidy$pl1/2-1
  ploidy2_id <- input.obj$ploidy$pl2/2-1
  res.temp <-
    .Call("est_hmm_map",
          ploidy1_id,
          ploidy2_id,
          input.obj$n.mrk,
          input.obj$n.ind,
          input.obj$states,
          input.obj$emit,
          rep(0.01, input.obj$n.mrk-1),
          verbose,
          tol,
          use_H0,
          PACKAGE = "mappolymp")
  structure(list(n.mrk = input.obj$n.mrk,
                 n.ind = input.obj$n.ind,
                 rf = res.temp[[2]],
                 loglike = res.temp[[1]],
                 states = input.obj$states,
                 emit = input.obj$emit,
                 phases = input.obj$merged.phases,
                 pedigree = input.obj$pedigree), class = "multi.mappoly.map")
}

#' @export
#' @importFrom dplyr  mutate select
calc_genoprob_mutli_fam <- function(consensus.map,
                                    verbose = TRUE,
                                    step = 0) {
  ploidy1_id <- consensus.map$pedigree$pl1/2-1
  ploidy2_id <- consensus.map$pedigree$pl2/2-1
  ##FIXME: change this for multiple ploidy leves
  pl <- consensus.map$pedigree$pl1[1]
  if(round(step, 1)  ==  0){
    x <-
      .Call("calc_genoprob_multi_fam",
            ploidy1_id,
            ploidy2_id,
            consensus.map$n.mrk,
            consensus.map$n.ind,
            consensus.map$states,
            consensus.map$emit,
            consensus.map$rf,
            verbose,
            PACKAGE = "mappolymp")
    prob.out <- as.data.frame(x)
    colnames(prob.out) <- c("ind", "snp", "par", "hom", "prob")
    pos <- cumsum(mappoly::imf_h(c(0, consensus.map$rf)))
    names(pos) <- names(consensus.map$states)
    prob.out$ind <- names(consensus.map$states[[1]])[prob.out$ind]
    prob.out$snp <- names(consensus.map$states)[prob.out$snp]
    prob.out$pos <- pos[prob.out$snp]
    head(prob.out)
    prob.out <- prob.out %>%
      mutate(homolog = paste0("P", par, "_h", hom)) %>%
      select(-one_of(c("par", "hom")))
    return(structure(list(genoprob = prob.out, pedigree = consensus.map$pedigree),
              class = "multi.mappoly.genoprob"))
  }
  else if(round(step, 1)  >  0){
    map <- create_map_multi_fam(consensus.map, step)
    R <- generate_biallelic_indices(pl)
    Rtemp <- lapply(1:consensus.map$n.ind, function(x) R[[1]][[1]][[1]])
    Etemp <- lapply(1:consensus.map$n.ind, function(x) matrix(rep(1/nrow(Rtemp[[1]]), nrow(Rtemp[[1]])), ncol = 1))
    names(Rtemp) <- names(Etemp) <- names(consensus.map$states[[1]])
    St <- Em <- vector("list", length(map))
    names(St) <- names(Em) <- names(map)
    for(i in names(map)){
      if(i %in% names(consensus.map$states)){
        St[[i]] <- consensus.map$states[[i]]
        Em[[i]] <- consensus.map$emit[[i]]
      } else {
        ###FIXME
        ###for different ploidies, include an extra for loop here
        St[[i]] <- Rtemp
        Em[[i]] <- Etemp
      }
    }
    x <-
      .Call("calc_genoprob_multi_fam",
            ploidy1_id,
            ploidy2_id,
            length(St),
            consensus.map$n.ind,
            St,
            Em,
            mf_h(diff(map)),
            verbose,
            PACKAGE = "mappolymp")
    prob.out <- as.data.frame(x)
    colnames(prob.out) <- c("ind", "snp", "par", "hom", "prob")
    prob.out$ind <- names(St[[1]])[prob.out$ind]
    prob.out$snp <- names(St)[prob.out$snp]
    prob.out$pos <- map[prob.out$snp]
    head(prob.out)
    prob.out <- prob.out %>%
      mutate(homolog = paste0("P", par, "_h", hom)) %>%
      select(-one_of(c("par", "hom")))
    return(structure(list(genoprob = prob.out, pedigree = consensus.map$pedigree),
              class = "multi.mappoly.genoprob"))
  } else stop("step must be positive")
}
