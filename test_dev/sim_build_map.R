myfunc <- function(ploidy.vec,
                   cross.mat,
                   n.ind,
                   n.mrk,
                   alleles,
                   map.length)
{
  sim.cross <- simulate_multiple_crosses(ploidy.vec,
                                         cross.mat,
                                         n.ind,
                                         n.mrk,
                                         alleles,
                                         map.length)
  states <- states_to_visit(sim.cross)
  restemp <- hmm_map_reconstruction(ploidy1 = states$ploidy$pl1,
                                    ploidy2 = states$ploidy$pl2,
                                    n.mrk = states$n.mrk,
                                    n.ind = states$n.ind,
                                    haplo = states$states,
                                    emit = states$emit,
                                    rf_vec = rep(0.01, states$n.mrk-1),
                                    verbose = FALSE,
                                    use_H0 = FALSE,
                                    tol = 1e-3)
  mpm <- cumsum(mappoly::imf_h(c(0, restemp[[2]])))
  names(mpm) <- names(states$states)
  return(data.frame(sim.cross$map, est.map = mpm))
}
