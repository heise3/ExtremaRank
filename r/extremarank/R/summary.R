as.data.frame.extremarank_result <- function(x, row.names = NULL, optional = FALSE, ...) {
    if (!is.null(x$features)) return(x$features)
    x$fit_status
}

summary.extremarank_result <- function(object, ...) {
    status <- if (is.null(object$features)) NULL else table(object$features$membership_status, useNA = "ifany")
    structure(list(status = object$status, method = object$method, candidates = status,
        direction_statuses = if (is.null(object$features)) NULL else table(object$features$direction_status),
        coverage = object$coverage, plan = object$plan, timing = object$timing,
        elapsed_seconds = object$elapsed_seconds,
        minimum_change = if (!is.null(object$audit)) object$audit$topk_minimum_change else object$minimum_observed_change,
        diagnostic_coverage = if (is.null(object$diagnostics)) NULL else
            c(evaluated = object$diagnostics$n_scenarios, feasible = object$diagnostics$feasible_scenarios),
        scope = object$scope), class = "summary.extremarank_result")
}

print.summary.extremarank_result <- function(x, ...) {
    cat(x$method, ": ", x$status, "\n", sep = "")
    if (!is.null(x$candidates)) { cat("Candidate membership:\n"); print(x$candidates) }
    if (!is.null(x$direction_statuses)) { cat("Candidate effect direction:\n"); print(x$direction_statuses) }
    if (!is.null(x$coverage)) {
        cat("Deletion refits: ", x$coverage$evaluated, " / ", x$coverage$planned,
            "; valid ", x$coverage$valid, "; failed ", x$coverage$failed,
            "; not run ", x$coverage$not_run, "\n", sep = "")
    }
    if (!is.null(x$diagnostic_coverage)) print(x$diagnostic_coverage)
    cat("Elapsed seconds: ", x$elapsed_seconds, "\n", x$scope, "\n", sep = "")
    invisible(x)
}

plot.extremarank_result <- function(x, type = c("status", "influence"), ...) {
    type <- match.arg(type)
    if (type == "status") {
        if (is.null(x$features)) stop("No evaluable candidate results to plot")
        counts <- table(x$features$membership_status)
        colors <- c(CERTIFIED = "#1B9E77", REFUTED = "#D95F02", UNRESOLVED = "#7570B3",
            OBSERVED_STABLE = "#66A61E", OBSERVED_CHANGED = "#E7298A",
            PARTIALLY_EVALUATED = "#7570B3", NOT_EVALUATED = "#666666")
        graphics::barplot(counts, col = colors[names(counts)], ylab = "Queried features",
            main = paste("Candidate membership:", x$status), ...)
    } else {
        scenarios <- if (!is.null(x$diagnostics)) x$diagnostics$scenarios else x$scenarios[-1L]
        valid <- Filter(function(s) !is.null(s$jaccard) && is.finite(s$jaccard), scenarios)
        if (!length(valid)) stop("No evaluated influence diagnostics; enable diagnostics or valid refits")
        values <- vapply(valid, function(s) 1 - s$jaccard, numeric(1))
        labels <- vapply(valid, function(s) if (!is.null(s$deleted_index) && !is.null(x$audit)) x$audit$unit_ids[s$deleted_index] else paste(s$deleted_units, collapse = ";"), character(1))
        graphics::barplot(values, names.arg = labels, las = 2, col = "#377EB8", ylab = "1 - Top-K Jaccard",
            main = paste(length(valid), "evaluated deletion scenarios"), ...)
    }
    invisible(x)
}

robustness_curve <- function(x, budgets = 0:2, ...) {
    if (!is.numeric(budgets) || any(!is.finite(budgets)) || any(budgets != floor(budgets))) stop("budgets must be nonnegative integers")
    budgets <- sort(unique(as.integer(budgets)))
    if (!length(budgets) || anyNA(budgets) || any(budgets < 0)) stop("budgets must be nonnegative integers")
    dots <- list(...)
    if ("budget" %in% names(dots)) stop("use budgets rather than budget")
    results <- lapply(budgets, function(b) do.call(extremarank_effects, c(list(effects = x, budget = b), dots)))
    data <- do.call(rbind, lapply(seq_along(results), function(i) {
        z <- results[[i]]; s <- z$features$membership_status
        data.frame(budget = budgets[i], certified = sum(s == "CERTIFIED"), refuted = sum(s == "REFUTED"),
            unresolved = sum(s == "UNRESOLVED"), elapsed_seconds = z$elapsed_seconds)
    }))
    structure(list(data = data, results = results), class = "extremarank_curve")
}

plot.extremarank_curve <- function(x, ...) {
    graphics::matplot(x$data$budget, x$data[, c("certified", "refuted", "unresolved")], type = "b",
        pch = c(16, 17, 15), lty = 1, col = c("#1B9E77", "#D95F02", "#7570B3"),
        xlab = "Maximum deleted donors", ylab = "Queried features", ...)
    graphics::legend("right", c("Certified", "Refuted", "Unresolved"), col = c("#1B9E77", "#D95F02", "#7570B3"), lty = 1, pch = c(16, 17, 15))
    invisible(x)
}
