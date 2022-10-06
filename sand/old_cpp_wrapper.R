#' HMM estimation function
#'
#' @param void internal function to be documented
#' @examples
#' require(tidyverse)
#' p <- c(2, 2)
#' names(p) <- c("P1", "P2")
#' cm <- matrix(c("P1","P2"),
#'              ncol = 2, byrow = T)
#' D <- simulate_multiple_crosses(ploidy = p, #four parents
#'                                cross.mat = cm,
#'                                n.ind = c(20), # per cross
#'                                n.mrk= c(3, 3), # per parent
#'                                alleles = list(c(0:1), c(0:1)),
#'                                map.length = 5)
#' D$dat
#' H <- states_to_visit(D, err = 0.0, is.log = TRUE)
#' names(H)
#' head(H$hmm.info)
#' mrk.sample <- names(D$map)
#' i1 <- H$hmm.info %>%
#'       filter(ind == "Ind_P1xP2_1", mrk %in% mrk.sample) %>%
#'       #filter(mrk %in% mrk.sample[2]) %>%
#'       mutate(st.all = paste(st.p1, st.p2, sep = "-"))
#' i1$mrk <- as_factor(i1$mrk)
#' p1 <- ggplot(i1, aes(x=st.all, y = emit, fill = emit)) +
#'      scale_fill_distiller(palette = "RdYlBu") +
#'      geom_bar(stat="identity") + facet_grid(.~mrk)
#' p1 + coord_flip()
#' x<-est_map_R(H)
#' round(x,2)
#' @export
# est_map_R <- function(states.hmm, rf.vec = NULL, tol = 10e-5, verbose = FALSE) {
#   ploidy.mat <- states.hmm$ploidy.cross.id
#   ##FIXME: check is there is change of order when dealing with multiple families
#   crosstype <- states.hmm$hmm.info %>% group_by(ind) %>% summarise(crosstype = unique(pl.id))
#   n.ind <- length(crosstype$crosstype)
#   o.mrk <-unique(states.hmm$hmm.info$mrk)
#   n.mrk <- length(o.mrk)
#   state.lengths <- states.hmm$hmm.info %>%
#     group_by(ind, mrk) %>%
#     summarise(n = n(), .groups = 'drop')
#   if(is.null(rf.vec))
#     rf.vec <- rep(0.01, n.mrk - 1)
#   res.temp <- .Call("est_log_hmm_map",
#                     ploidy.mat,
#                     n.mrk,
#                     n.ind,
#                     as.double(rf.vec),
#                     as.integer(crosstype$crosstype-1),
#                     as.integer(state.lengths$n),
#                     as.integer(cumsum(c(0, state.lengths$n))),
#                     as.integer(states.hmm$hmm.info$st.all),
#                     as.double(states.hmm$hmm.info$emit),
#                     tol = tol,
#                     verbose = verbose,
#                     PACKAGE = "mappoly2")
#   res.temp
# }

