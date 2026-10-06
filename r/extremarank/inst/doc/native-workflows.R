## -----------------------------------------------------------------------------
library(extremarank)
effects <- matrix(c(5, 4, 6, 5, 1, 2, 1, 2, -2, -1, -3, -2), nrow = 4,
  dimnames = list(paste0("donor", 1:4), c("A", "B", "C")))
result <- extremarank_effects(effects, k = 1, budget = 1,
  search_order = "influence", max_diagnostics = 4)
summary(result)
as.data.frame(result)
plot(result)
plot(result, type = "influence")
curve <- robustness_curve(effects, budgets = 0:2, k = 1, diagnostics = FALSE)
curve$data
plot(curve)

## -----------------------------------------------------------------------------
set.seed(5006)
x <- matrix(rnorm(30 * 8), 30, dimnames = list(paste0("g", 1:30), paste0("s", 1:8)))
meta <- data.frame(sample_id = colnames(x), donor_id = colnames(x),
  group = rep(c("treated", "control"), each = 4), site = rep(1:4, 2))
preflight_extremarank(x, meta, "treated", "control", design = "welch", budget = 2)

## ----eval=requireNamespace("limma", quietly=TRUE)-----------------------------
fits <- refit_extremarank(x, meta, "treated", "control", model = "limma",
  design = "welch", budget = 1, k = 5, all_features = TRUE)
summary(fits)
head(as.data.frame(fits))
zero <- refit_extremarank(x, meta, "treated", "control", model = "limma",
  design = "welch", max_refits = 0)
zero$status

