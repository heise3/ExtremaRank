.er_capture_fit <- function(fun, keep, seed = 51006L) {
    old_kind <- RNGkind(); had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({ do.call(RNGkind, as.list(old_kind)); if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv) }, add = TRUE)
    set.seed(seed, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
    warnings <- messages <- character()
    result <- tryCatch(withCallingHandlers(fun(keep),
        warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") },
        message = function(m) { messages <<- c(messages, conditionMessage(m)); invokeRestart("muffleMessage") }),
        error = function(e) e)
    list(result = result, warnings = unique(warnings), messages = unique(messages))
}

.er_refit_engine <- function(x, meta, unit_ids, sample_units, plan, weights, paired,
    model, direction, k, budget, max_refits, fit, preparation, versions,
    features, all_features, store_scores, score_path, checkpoint, resume, progress, workers, seed = 51006L) {
    store_scores <- match.arg(store_scores, c("memory", "disk", "none"))
    workers <- .er_int(workers, "workers", 1, 4)
    seed <- .er_int(seed, "seed", 0, .Machine$integer.max - 1000001L)
    if (all_features && !is.null(features)) stop("use features or all_features, not both")
    if (!is.null(features)) {
        features <- .er_ids(features, "features")
        if (!all(features %in% rownames(x))) stop("queried features are outside the frozen universe")
    }
    signature <- digest::digest(list(x, meta, unit_ids, sample_units, weights, paired, model,
        direction, k, budget, preparation, versions, features, all_features, seed,
        package = as.character(utils::packageVersion("extremarank"))), algo = "sha256")
    if (resume && is.null(checkpoint)) stop("resume requires checkpoint")
    if (!is.null(checkpoint)) {
        if (length(checkpoint) != 1L || !nzchar(checkpoint)) stop("declare one checkpoint directory")
        marker <- file.path(checkpoint, "contract.rds")
        if (resume) {
            if (!file.exists(marker) || !identical(readRDS(marker), signature)) stop("checkpoint input, parameters or package/model versions changed")
        } else {
            if (dir.exists(checkpoint) && length(list.files(checkpoint, all.files = TRUE, no.. = TRUE))) stop("checkpoint must be new or empty unless resume=TRUE")
            dir.create(checkpoint, recursive = TRUE, showWarnings = FALSE)
            saveRDS(signature, marker)
        }
        checkpoint <- normalizePath(checkpoint)
    }
    if (store_scores == "disk") {
        if (is.null(score_path)) stop("store_scores='disk' requires score_path")
        marker <- file.path(score_path, "contract.rds")
        if (file.exists(marker)) {
            if (!resume || !identical(readRDS(marker), signature)) stop("score_path contract differs or resume was not requested")
        } else {
            if (dir.exists(score_path) && length(list.files(score_path, all.files = TRUE, no.. = TRUE))) stop("score_path must be new or empty")
            dir.create(score_path, recursive = TRUE, showWarnings = FALSE)
            saveRDS(signature, marker)
        }
        score_path <- normalizePath(score_path)
    }
    masks <- list(integer()); capped <- FALSE
    totals <- colSums(weights)
    visit <- function(prefix, start, need) {
        if (capped) return(invisible(NULL))
        if (!need) {
            removed <- if (length(prefix)) colSums(weights[prefix, , drop = FALSE]) else c(0, 0)
            if (totals[1L] - removed[1L] < 2L || (!paired && totals[2L] - removed[2L] < 2L)) return(invisible(NULL))
            if (length(masks) - 1L >= max_refits) { capped <<- TRUE; return(invisible(NULL)) }
            masks[[length(masks) + 1L]] <<- prefix
        } else if (start <= length(unit_ids) - need + 1L) for (i in seq.int(start, length(unit_ids) - need + 1L)) {
            visit(c(prefix, i), i + 1L, need - 1L)
            if (capped) break
        }
        invisible(NULL)
    }
    if (budget) for (r in seq_len(budget)) { visit(integer(), 1L, r); if (capped) break }
    keep_sets <- lapply(masks, function(d) which(!sample_units %in% unit_ids[d]))
    sids <- c("baseline", if (length(masks) > 1L) paste0("delete_", seq_len(length(masks) - 1L)))
    cached <- vector("list", length(masks))
    if (!is.null(checkpoint) && resume) for (i in seq_along(masks)) {
        path <- file.path(checkpoint, paste0(sids[i], ".rds"))
        if (file.exists(path)) cached[[i]] <- path
    }
    fit_status <- scores <- score_files <- scenarios <- failed_features <- list()
    baseline_topk <- baseline_signs <- baseline_member <- baseline_rank <- NULL
    selected <- if (all_features) rownames(x) else features
    member_changed <- sign_changed <- member_min <- sign_min <- rank_best <- rank_worst <- selection_count <- NULL
    member_witness <- sign_witness <- list()
    started <- proc.time()[[3L]]
    cluster <- NULL
    if (workers > 1L && any(vapply(cached, is.null, logical(1)))) {
        cluster <- parallel::makePSOCKcluster(workers)
        on.exit(parallel::stopCluster(cluster), add = TRUE)
        parallel::clusterCall(cluster, function(paths) .libPaths(paths), .libPaths())
    }
    for (batch_start in seq.int(1L, length(masks), by = workers)) {
        batch <- seq.int(batch_start, min(length(masks), batch_start + workers - 1L))
        missing <- batch[vapply(cached[batch], is.null, logical(1))]
        if (length(missing)) {
            runner <- .er_fit_task
            tasks <- lapply(missing, function(i) list(keep = keep_sets[[i]], seed = as.integer((as.double(seed) + i) %% .Machine$integer.max)))
            fresh <- if (is.null(cluster)) lapply(tasks, runner, fit = fit) else
                parallel::parLapply(cluster, tasks, runner, fit = fit)
            for (j in seq_along(missing)) {
                i <- missing[j]; cached[[i]] <- fresh[[j]]
                if (!is.null(checkpoint)) {
                    tmp <- tempfile(tmpdir = checkpoint); saveRDS(fresh[[j]], tmp)
                    if (!file.rename(tmp, file.path(checkpoint, paste0(sids[i], ".rds")))) stop("checkpoint write failed")
                }
            }
        }
        for (i in batch) {
            record <- cached[[i]]; if (is.character(record)) record <- readRDS(record)
            result <- record$result; sid <- sids[i]
            detail <- "Full fixed-universe ranking fitted"; invalid <- character()
            if (!inherits(result, "error")) {
                valid <- is.data.frame(result) && all(c("feature_id", "score", "effect") %in% names(result)) &&
                    !anyDuplicated(result$feature_id) && setequal(result$feature_id, rownames(x)) &&
                    is.numeric(result$score) && is.numeric(result$effect)
                if (!valid) result <- simpleError("fit must return each frozen feature exactly once, with numeric score and effect")
                else { result$feature_id <- as.character(result$feature_id)
                    result <- result[match(rownames(x), result$feature_id), , drop = FALSE]
                    invalid <- result$feature_id[!is.finite(result$score) | !is.finite(result$effect)] }
            }
            if (inherits(result, "error")) { status <- "NOT_EVALUABLE"; detail <- conditionMessage(result); result <- NULL }
            else if (length(invalid)) { status <- "NOT_EVALUABLE"; detail <- paste(length(invalid), "nonfinite scores/effects in the frozen universe") }
            else status <- "OK"
            fit_status[[i]] <- data.frame(scenario_id = sid, status = status, detail = detail,
                warnings = paste(record$warnings, collapse = " | "), messages = paste(record$messages, collapse = " | "), stringsAsFactors = FALSE)
            failed_features[sid] <- list(invalid)
            if (store_scores == "memory") scores[sid] <- list(result)
            if (store_scores == "disk" && !is.null(result)) {
                path <- file.path(score_path, paste0(sid, "_scores.rds")); saveRDS(result, path)
                score_files[sid] <- list(path)
            }
            topk <- NULL; changed <- NA; sign_changes <- NA_integer_
            deleted <- unit_ids[masks[[i]]]
            if (status == "OK") {
                ranking <- switch(direction, up = result$score, down = -result$score, absolute = abs(result$score))
                order <- order(-ranking, enc2utf8(result$feature_id), method = "radix")
                topk <- result$feature_id[order[seq_len(k)]]
                signs <- stats::setNames(sign(result$effect), result$feature_id)
                if (i == 1L) {
                    baseline_topk <- topk; baseline_signs <- signs
                    if (is.null(selected)) selected <- topk
                    baseline_member <- selected %in% topk
                    baseline_rank <- match(selected, result$feature_id[order])
                    rank_best <- rank_worst <- baseline_rank; selection_count <- integer(length(selected))
                    member_changed <- sign_changed <- integer(length(selected))
                    member_min <- sign_min <- rep(NA_integer_, length(selected))
                    member_witness <- sign_witness <- vector("list", length(selected))
                } else if (!is.null(baseline_topk)) {
                    changed <- !setequal(topk, baseline_topk); sign_changes <- sum(signs != baseline_signs)
                    ranks <- match(selected, result$feature_id[order]); rank_best <- pmin(rank_best, ranks); rank_worst <- pmax(rank_worst, ranks)
                    selection_count <- selection_count + (selected %in% topk)
                    a <- (selected %in% topk) != baseline_member; b <- signs[selected] != baseline_signs[selected]
                    member_changed <- member_changed + a; sign_changed <- sign_changed + b
                    for (j in which(a & is.na(member_min))) { member_min[j] <- length(deleted); member_witness[[j]] <- deleted }
                    for (j in which(b & is.na(sign_min))) { sign_min[j] <- length(deleted); sign_witness[[j]] <- deleted }
                }
            }
            scenarios[[i]] <- list(scenario_id = sid, deleted_units = deleted, deleted_indices = masks[[i]], deleted_samples = meta$sample_id[sample_units %in% deleted], status = status, topk = topk,
                topk_changed = changed, effect_sign_changes = sign_changes,
                jaccard = if (!is.null(topk) && !is.null(baseline_topk)) length(intersect(topk, baseline_topk)) / length(union(topk, baseline_topk)) else NA_real_)
            if (progress) message("ExtremaRank ", i, "/", length(masks), ": ", sid, " ", status)
            cached[i] <- list(NULL)
        }
    }
    table <- do.call(rbind, fit_status)
    evaluated <- nrow(table) - 1L; valid <- sum(table$status[-1L] == "OK"); failed <- evaluated - valid
    remaining <- cpp_remaining(plan$feasible_deletions, evaluated)
    changed <- any(vapply(scenarios[-1L], function(s) isTRUE(s$topk_changed), logical(1)))
    overall <- if (table$status[1L] != "OK") "NOT_EVALUABLE" else if (!evaluated) "NOT_EVALUATED" else
        if (changed) "OBSERVED_CHANGED" else if (failed || capped) "PARTIALLY_EVALUATED" else "OBSERVED_STABLE"
    candidate_status <- function(n) ifelse(n > 0L, "OBSERVED_CHANGED",
        if (!evaluated) "NOT_EVALUATED" else if (failed || capped) "PARTIALLY_EVALUATED" else "OBSERVED_STABLE")
    candidates <- NULL
    if (!is.null(baseline_topk)) {
        candidates <- data.frame(feature_id = selected, baseline_rank = baseline_rank, baseline_topk = baseline_member,
            membership_status = candidate_status(member_changed), direction_status = candidate_status(sign_changed),
            membership_observed_changes = member_changed, direction_observed_changes = sign_changed,
            best_observed_rank = rank_best, worst_observed_rank = rank_worst,
            observed_selection_fraction = if (valid) selection_count / valid else rep(NA_real_, length(selected)),
            membership_minimum_observed_change = member_min, direction_minimum_observed_change = sign_min,
            valid_deletion_refits = valid, failed_deletion_refits = failed, stringsAsFactors = FALSE)
        candidates$membership_witness <- I(member_witness); candidates$direction_witness <- I(sign_witness)
    }
    witnesses <- Filter(function(s) isTRUE(s$topk_changed), scenarios)
    minimum <- if (length(witnesses)) min(vapply(witnesses, function(s) length(s$deleted_units), integer(1))) else NA_integer_
    structure(list(status = overall, method = paste("native R", model, "refits"), model = model, features = candidates,
        direction = direction, top_k = k, budget = budget, fit_status = table, failed_features = failed_features,
        scores = scores, score_files = score_files, store_scores = store_scores, scenarios = scenarios,
        baseline_topk = baseline_topk, minimum_observed_change = minimum, max_refits = max_refits,
        enumeration_capped = capped, coverage = list(planned = plan$feasible_deletions, evaluated = evaluated,
            valid = valid, failed = failed, not_run = remaining, baseline_status = table$status[1L]),
        plan = plan, deletion_units = data.frame(sample_id = meta$sample_id, deletion_unit = sample_units),
        unit_ids = unit_ids, workers = workers, seed = seed, checkpoint = checkpoint, contract_sha256 = signature,
        preparation = preparation, versions = versions, elapsed_seconds = proc.time()[[3L]] - started,
        version = as.character(utils::packageVersion("extremarank")),
        scope = "Whole Top-K membership status; separate per-candidate observed membership/effect-sign results. No certificate for unseen refits, FDR or biological validity"),
        class = "extremarank_result")
}

refit_extremarank_custom <- function(x, metadata, fit_function, unit_column = "donor_id",
    k = NULL, budget = 1L, direction = c("up", "down", "absolute"), features = NULL,
    all_features = FALSE, max_refits = 1000L, store_scores = c("memory", "disk", "none"),
    score_path = NULL, checkpoint = NULL, resume = FALSE, progress = FALSE, workers = 1L,
    max_dense_bytes = 512 * 1024^2, block_rows = 1000L, seed = 51006L) {
    if (!is.function(fit_function)) stop("fit_function must be a function(x, metadata)")
    z <- .er_input(x, metadata, NULL); meta <- z$metadata
    if (is.null(meta) || !unit_column %in% names(meta)) stop("declare an existing deletion unit column")
    members <- .er_ids(meta[[unit_column]], unit_column, FALSE); units <- unique(members)
    x <- .er_realize(z$x, max_dense_bytes, block_rows)
    .er_memory(nrow(x), ncol(x) * workers, max_dense_bytes)
    if (any(!is.finite(x))) stop("custom model inputs must be finite")
    budget <- .er_int(budget, "budget", 0, length(units) - 2L)
    max_refits <- .er_int(max_refits, "max_refits", 0)
    k <- if (is.null(k)) min(20L, nrow(x)) else .er_int(k, "k", 1, nrow(x))
    plan <- cpp_plan(integer(length(units)), TRUE, budget)
    prep <- list(unit_column = unit_column, metadata = meta, filter_frozen = TRUE,
        fit_function_sha256 = digest::digest(fit_function, algo = "sha256"),
        contract = "User function refits each retained sample matrix and returns frozen feature_id, score and effect; external data and package versions must be controlled by the user")
    fit <- function(keep) fit_function(x[, keep, drop = FALSE], meta[keep, , drop = FALSE])
    out <- .er_refit_engine(x, meta, units, members, plan, cbind(rep(1L, length(units)), 0L), TRUE,
        "custom", match.arg(direction), k, budget, max_refits, fit, prep,
        c(R = as.character(getRversion())), features, all_features, store_scores,
        score_path, checkpoint, resume, progress, workers, seed)
    out$design <- "declared custom model"; out
}

.er_fit_task <- function(task, fit) .er_capture_fit(fit, task$keep, task$seed)
