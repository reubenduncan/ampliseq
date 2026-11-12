process DADA2_DENOISING {
    tag "${meta.run}"
    label 'process_high'
    label 'process_long'
    label 'error_retry'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/81/81153df5d53322e6d91b2c4c9bc4da50774fb1d101ead002fe75bb75fc6f036c/data'
        : 'community.wave.seqera.io/library/bioconductor-dada2_r-base_r-digest_tbb:38acac09bac46f36'}"

    input:
    tuple val(meta), path("filtered/*"), path(errormodel)

    output:
    tuple val(meta), path("*.dada.rds"), emit: denoised
    tuple val(meta), path("*.seqtab.rds"), emit: seqtab
    tuple val(meta), path("*.mergers.rds"), emit: mergers
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml", emit: versions_dada2_denoising, topic: versions
    path "*.args.txt", emit: args

    script:
    def prefix = task.ext.prefix ?: "prefix"
    def quality_type = task.ext.quality_type ?: "Auto"
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def quantile = task.ext.quantile ?: 0.001
    def pool = params.sample_inference == 'pooled' ? 'TRUE' : 'FALSE'
    def pool_label = params.sample_inference == 'pseudo' ? 'pseudo' : pool
    def parallel_samples = task.ext.parallel_samples ?: task.cpus
    def pseudo_denoise_function = '''
    dada_pseudo_worker <- function(input, err, dada_options, priors, thread_count, quality_type) {
        if (quality_type != "Auto") {
            input <- dada2::derepFastq(input, qualityType = quality_type)
        }
        call_args <- c(
            list(
                derep = input,
                err = err,
                pool = FALSE,
                priors = priors,
                multithread = thread_count,
                verbose = FALSE
            ),
            dada_options
        )
        do.call(dada2::dada, call_args)
    }

    accumulate_dada_transitions <- function(results) {
        trans <- lapply(results, function(result) result$trans)
        max_columns <- max(vapply(trans, ncol, integer(1)))
        counts <- matrix(0, nrow = 16, ncol = max_columns)
        rownames(counts) <- c(
            "A2A", "A2C", "A2G", "A2T", "C2A", "C2C", "C2G", "C2T",
            "G2A", "G2C", "G2G", "G2T", "T2A", "T2C", "T2G", "T2T"
        )
        colnames(counts) <- seq(0, max_columns - 1)
        for (sample_trans in trans) {
            columns <- seq_len(ncol(sample_trans))
            counts[, columns] <- counts[, columns] + sample_trans
        }
        counts
    }

    estimate_shared_error <- function(results) {
        tryCatch(
            suppressWarnings(dada2::loessErrfun(accumulate_dada_transitions(results))),
            error = function(condition) NULL
        )
    }

    dada_pseudo_parallel <- function(files, err, dada_options, ncpus, max_workers, log_prefix, quality_type) {
        files <- sort(files)
        sample_names <- basename(files)
        ncpus <- max(1L, as.integer(ncpus))
        if (length(files) == 1L) {
            return(dada_pseudo_worker(
                files[[1]], err, dada_options, character(0), ncpus, quality_type
            ))
        }

        max_workers <- max(1L, as.integer(max_workers))
        workers <- min(ncpus, max_workers, length(files))
        threads_per_worker <- max(1L, floor(ncpus / workers))
        cat(
            "DADA2 pseudo-pooling for", length(files), "samples with", workers,
            "workers and", threads_per_worker, "DADA2 thread(s) per worker.\\n"
        )
        cluster <- parallel::makePSOCKcluster(
            workers,
            outfile = paste0(log_prefix, ".workers.log")
        )
        on.exit(parallel::stopCluster(cluster), add = TRUE)
        parallel::clusterEvalQ(cluster, suppressPackageStartupMessages(library(dada2)))

        first_pass <- parallel::parLapplyLB(
            cluster,
            files,
            dada_pseudo_worker,
            err = err,
            dada_options = dada_options,
            priors = character(0),
            thread_count = threads_per_worker,
            quality_type = quality_type
        )
        names(first_pass) <- sample_names

        transition_counts <- accumulate_dada_transitions(first_pass)
        shared_err <- tryCatch(
            suppressWarnings(dada2::loessErrfun(transition_counts)),
            error = function(condition) NULL
        )
        if (is.null(shared_err)) {
            stop("DADA2 could not estimate the shared error model after pseudo-pooling pass 1.")
        }

        first_pass_table <- dada2::makeSequenceTable(first_pass)
        pseudo_prevalence <- if ("PSEUDO_PREVALENCE" %in% names(dada_options)) {
            dada_options[["PSEUDO_PREVALENCE"]]
        } else {
            2
        }
        pseudo_abundance <- if ("PSEUDO_ABUNDANCE" %in% names(dada_options)) {
            dada_options[["PSEUDO_ABUNDANCE"]]
        } else {
            Inf
        }
        pseudo_priors <- colnames(first_pass_table)[
            colSums(first_pass_table > 0) >= pseudo_prevalence |
                colSums(first_pass_table) >= pseudo_abundance
        ]
        rm(first_pass, first_pass_table, transition_counts)

        second_pass <- parallel::parLapplyLB(
            cluster,
            files,
            dada_pseudo_worker,
            err = shared_err,
            dada_options = dada_options,
            priors = pseudo_priors,
            thread_count = threads_per_worker,
            quality_type = quality_type
        )
        names(second_pass) <- sample_names
        shared_err_out <- estimate_shared_error(second_pass)
        for (i in seq_along(second_pass)) {
            second_pass[[i]][["err_out"]] <- shared_err_out
        }
        second_pass
    }
    '''
    if (!meta.single_end) {
        """
        #!/usr/bin/env Rscript
        suppressPackageStartupMessages(library(dada2))

        ${pseudo_denoise_function}

        errF <- readRDS("${errormodel[0]}")
        errR <- readRDS("${errormodel[1]}")

        filtFs <- list.files("./filtered/", pattern = "_1.filt.fastq.gz", full.names = TRUE)
        filtRs <- sub("_1.filt.fastq.gz", "_2.filt.fastq.gz", filtFs)

        #denoising
        sink(file = "${prefix}.dada.log")
        if ("${params.sample_inference}" == "pseudo") {
            dada_options <- list(${args})
            dadaFs <- dada_pseudo_parallel(filtFs, errF, dada_options, ${task.cpus}, ${parallel_samples}, "${prefix}.forward", "${quality_type}")
            dadaRs <- dada_pseudo_parallel(filtRs, errR, dada_options, ${task.cpus}, ${parallel_samples}, "${prefix}.reverse", "${quality_type}")
        } else if ("${quality_type}" == "Auto") {
            # Avoid using memory-inefficient derepFastq() if not necessary
                dadaFs <- dada(filtFs, err = errF, pool = ${pool}, ${args}, multithread = ${task.cpus})
                dadaRs <- dada(filtRs, err = errR, pool = ${pool}, ${args}, multithread = ${task.cpus})
        } else {
            derepFs <- derepFastq(filtFs, qualityType="${quality_type}")
            dadaFs <- dada(derepFs, err = errF, pool = ${pool}, ${args}, multithread = ${task.cpus})
            derepRs <- derepFastq(filtRs, qualityType="${quality_type}")
            dadaRs <- dada(derepRs, err = errR, pool = ${pool}, ${args}, multithread = ${task.cpus})
        }
        saveRDS(dadaFs, "${prefix}_1.dada.rds")
        saveRDS(dadaRs, "${prefix}_2.dada.rds")
        sink(file = NULL)

        # merge
        if ("${params.mergepairs_strategy}" == "consensus") {
            mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, ${args2}, justConcatenate = FALSE, verbose=TRUE)
            concats <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, ${args2}, justConcatenate = TRUE, verbose=TRUE)

            # in case there is only one sample in the entire run
            if (is.data.frame(mergers)) {
                mergers <- list(sample = mergers)
                concats <- list(sample = concats)
            }

            # define the overlap threshold to decide if concatenation or not
            min_overlap_obs <- lapply(mergers, function(X) {
                mergers_accepted <- X[["accept"]]
                if (sum(mergers_accepted) > 0) {
                    min_overlap_obs <- X[["nmatch"]][mergers_accepted] + X[["nmismatch"]][mergers_accepted]
                    rep(min_overlap_obs, X[["abundance"]][mergers_accepted])
                } else {
                    NA
                }
            })

            min_overlap_obs <- Reduce(c, min_overlap_obs)
            min_overlap_obs <- min_overlap_obs[!is.na(min_overlap_obs)]
            min_overlap_obs <- quantile(min_overlap_obs, ${quantile})

            for (x in names(mergers)) {
                to_concat <- !mergers[[x]][["accept"]] & (mergers[[x]][["nmismatch"]] + mergers[[x]][["nmatch"]]) < min_overlap_obs

                if (sum(to_concat) > 0) {
                    mergers[[x]][to_concat, ] <- concats[[x]][to_concat, ]
                    # filter out unaccepted non concatenated sequences
                    mergers[[x]] <- mergers[[x]][mergers[[x]][["accept"]], ]
                }

            }

            # if one sample, need to convert back to df for next steps

            if(length(mergers) == 1) {
                mergers <- mergers[[1]]
            }

        } else {
            mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, ${args2}, verbose=TRUE)
        }

        saveRDS(mergers, "${prefix}.mergers.rds")

        # make table
        seqtab <- makeSequenceTable(mergers)
        saveRDS(seqtab, "${prefix}.seqtab.rds")

        write.table('dada\t${args}, pool = ${pool_label}', file = "dada.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
        write.table('mergePairs\t${args2}', file = "mergePairs.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
        writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
        """
    }
    else {
        """
        #!/usr/bin/env Rscript
        suppressPackageStartupMessages(library(dada2))

        ${pseudo_denoise_function}

        errF = readRDS("${errormodel}")

        filtFs <- sort(list.files("./filtered/", pattern = ".fastq.gz", full.names = TRUE))

        #denoising
        sink(file = "${prefix}.dada.log")
        if ("${params.sample_inference}" == "pseudo") {
            dada_options <- list(${args})
            dadaFs <- dada_pseudo_parallel(filtFs, errF, dada_options, ${task.cpus}, ${parallel_samples}, "${prefix}.single", "${quality_type}")
        } else if ("${quality_type}" == "Auto") {
            # Avoid using memory-inefficient derepFastq() if not necessary
                dadaFs <- dada(filtFs, err = errF, pool = ${pool}, ${args}, multithread = ${task.cpus})
        } else {
            derepFs <- derepFastq(filtFs, qualityType="${quality_type}")
            dadaFs <- dada(derepFs, err = errF, pool = ${pool}, ${args}, multithread = ${task.cpus})
        }
        saveRDS(dadaFs, "${prefix}.dada.rds")
        sink(file = NULL)

        #make table
        seqtab <- makeSequenceTable(dadaFs)
        saveRDS(seqtab, "${prefix}.seqtab.rds")

        #dummy file to fulfill output rules
        saveRDS("dummy", "dummy_${prefix}.mergers.rds")

        write.table('dada\t${args}, pool = ${pool_label}', file = "dada.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
        writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
        """
    }
}
