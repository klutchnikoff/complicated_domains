library(tidyverse)
library(doParallel)
library(doFuture)
library(future)
library(doRNG)
library(future.batchtools)
library(rngtools)
library(bernoulliR1)

source('./polynomial_sector.R')

##
## PARAMETERS
##

N.NODE <- 42
N.CPU <- 48
walltime <- "01:45:00"

NN <- c(200, 500, 1000, 2000)
HH <- seq(0.01, 0.6, by = 0.001)
KK <- c(1, 2.1)
MM <- 0:5

##
## PARALLEL PLAN
##

if (Sys.getenv("USER") == "klutchnikoff") {
  plan(list(multisession, multisession))
} else {
  myslurm <- tweak(batchtools_slurm, resources = list("ncpus" = N.CPU, "walltime" = walltime))
  plan(list(myslurm, multisession))
}

##
## LOOPS
##

rng <- RNGseq(N.CPU * N.NODE, 1234)

risk_poly <-
  foreach(
    r_node = 1:N.NODE,
    .combine = "rbind",
    .options.future = list(seed = TRUE)
  ) %dofuture% {
    foreach(
      r_cpu = 1:N.CPU,
      r = rng[(r_node - 1) * N.CPU + 1:N.CPU],
      .combine = "rbind",
      .options.future = list(seed = TRUE)
    ) %dofuture% {
      setRNG(r)
      res_tmp <-
        tibble(
          rep_node = numeric(),
          rep_cpu = numeric(),
          method = character(),
          k = numeric(),
          n = numeric(),
          h = numeric(),
          value = numeric()
        )
      for (n in NN) {
        for (k in KK) {
	  # closest points that are not NA
          if (k == 1) {idx <- c(1, 1)} else {idx <- c(1, 10)}
	  domain <- polynomial_sector(k)
          f <- function(x, y) {
            f_poly(x, y, k)
          }
          f00 <- f(0, 0)
          
	  # Strangely (0,0) does not belong to domain for spatstat ! Bug?
          eps <- 0.001
          zero <- spatstat.geom::ppp(eps, eps ^ k / 2, domain)
          data <- spatstat.random::rpoint(n, f, win = domain)
	
          for (h in HH) {
            ##
            ## Risk of sparr method
            ##

            f_sparr <- sparr::bivariate.density(data, h)$z
	    f_sparr <- as.vector(f_sparr[idx[1], idx[2]])
            ##
            ## res_tmp
            ##
            res_tmp <- rbind(
              res_tmp,
              tibble(
                rep_node = r_node,
                rep_cpu = r_cpu,
                h = h,
                n = n,
                k = k,
                method = "SPARR",
                value = (f_sparr - f00) ** 2
              )
            )
            for (m in MM) {
              mylp <- paste0("LP", m)
              assign(
                mylp,
                density_estimation(
                  data,
                  domain,
                  at_points = zero,
                  bandwidth = h,
                  degree = m
                )
              )
              res_tmp <- rbind(
                res_tmp,
                tibble(
                  rep_node = r_node,
                  rep_cpu = r_cpu,
                  h = h,
                  n = n,
                  k = k,
                  method = mylp,
                  value = (get(mylp) - f00) ** 2
                )
              )
            }
          }
        }
      }
      res_tmp
    }
  }

save(risk_poly, file = "../data/risk_poly_2000.rda")
