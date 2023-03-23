##  This version of tkmeans is a special case of tclust(), i.e.
##  equal.weights=TRUE, restr="eigen" (default) and rest.fact=1 (maximal restriction),
##  however, it is slower and cannot handle data with p>n, for which purpose was
##  developed the new tkmeans() function.
.tkmeans.old <-
function (x, k = 3, alpha = 0.05, nstart = 50, iter.max = 20, center = 0, scale = 1,
          store.x = TRUE, drop.empty.clust = TRUE, trace = 0, warnings = 2,
          zero.tol = 1e-16)
{
  ret <- tclust(x = x, k = k, restr.fact = 1, restr = "eigen", alpha = alpha,
                  nstart = nstart, iter.max = iter.max, equal.weights = TRUE,
                  store.x = store.x, center = center, scale = scale,
                  drop.empty.clust = drop.empty.clust, trace = trace,
                  zero.tol = zero.tol, warnings = min (warnings, 2))

  ret$obj = ret$cov[1,1,1]
  ret
}
