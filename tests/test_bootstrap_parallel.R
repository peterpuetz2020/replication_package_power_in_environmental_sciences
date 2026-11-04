library(doParallel)
source("scripts/functions.R")
dat <- data.frame(
  yi = c(.1, -.2, .3, -.4, .5, -.6), vi = 1,
  tau2 = rep(c(0, .1, .3), each = 2), GE = rep(c(.2, .8, 1.5), each = 2),
  cID = rep(c("a", "b", "c"), each = 2),
  metaID = c("one", "one", "two", "two", "two", "two")
)
dat <- split(dat, dat$cID)
clusters <- counterfactual_meta_clusters(dat)
registerDoSEQ()
serial <- cf.ci.cluster(dat, counterfactual_p_grid(), 40L, clusters, .5)
cl <- makeCluster(2L)
tryCatch({
  registerDoParallel(cl)
  parallel <- cf.ci.cluster(dat, counterfactual_p_grid(), 40L, clusters, .5)
  batches <- lapply(list(1:17, 18:40), function(ids) {
    cf.ci.cluster(dat, counterfactual_p_grid(), length(ids), clusters, .5,
                  iteration_ids = ids)
  })
  batched <- lapply(1:3, function(i) do.call(rbind, lapply(batches, `[[`, i)))
  stopifnot(isTRUE(all.equal(serial, parallel, tolerance = 0)),
            isTRUE(all.equal(serial, batched, tolerance = 0)))
}, finally = {
  stopCluster(cl)
  registerDoSEQ()
})
cat("PASS: article bootstrap is identical across serial, parallel and batched execution\n")
