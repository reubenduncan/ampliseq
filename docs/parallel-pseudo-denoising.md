# Parallel DADA2 pseudo-pooling

The pseudo-pooling path in `DADA2_DENOISING` runs independent samples concurrently inside the existing Nextflow task. It uses local PSOCK workers, so all workers remain inside the same task allocation and on the same compute node.

For each read direction, the task follows the two inference passes used by DADA2 pseudo-pooling:

- Pass 1 denoises each sample with the supplied error model and no pseudo-priors.
- The task sums transition counts across samples and estimates one shared error model with `loessErrfun`.
- The task builds a sequence table from pass 1 and selects pseudo-priors using DADA2's `PSEUDO_PREVALENCE` and `PSEUDO_ABUNDANCE` options.
- Pass 2 denoises each sample with the shared error model and shared pseudo-priors.
- The existing merge-pairs and sequence-table steps then run as before.

The forward and reverse reads are processed in separate two-pass groups, matching the pipeline's previous DADA2 calls. Samples remain grouped by sequencing run and region, as they are in the standard workflow.

The worker count defaults to the smaller of the DADA2 task's `task.cpus` value and the number of samples in that group. Each worker receives at least one DADA2 thread, with the allocated threads divided among workers. A one-sample group uses the ordinary single-sample DADA2 call. Set `ext.parallel_samples` to cap the worker count independently of the CPU allocation.

For a single-node run with 100 CPUs and 256 GB of memory, configure the denoising task in the run's Nextflow config:

```groovy
process {
    withName: DADA2_DENOISING {
        cpus = 100
        memory = 200.GB
        ext.parallel_samples = 100
    }
}
```

The task memory request should account for the combined memory of active R workers. If memory use is too high, lower `ext.parallel_samples` while keeping the CPU allocation. For example, 100 CPUs with 32 workers gives each active worker about three DADA2 threads. This approach has not been benchmarked on the Thames dataset. Compare a small multi-sample run against the unmodified v3.1.0 pipeline before processing the full cohort.

The PBS head job must also reserve at least the CPU and memory amounts assigned to `DADA2_DENOISING`. With the `local` executor, Nextflow does not create separate PBS jobs for these workers.
