# Running EM-seq Pipeline with Slurm

This guide explains how to run the EM-seq pipeline using Slurm.

## Slurm Configuration

The EM-seq pipeline is configured to use Slurm as an executor through the `slurm` profile in `nextflow.config`. The default configuration is set to use the `normal` partition, which is the default partition in the Docker-based Slurm cluster provided with this repository.

### Default Slurm Settings

```groovy
slurm {
    process.executor = 'slurm'
    process.queue = 'normal'  // Default queue/partition
    process.clusterOptions = '--account=your_account'  // Replace with your SLURM account
    
    // SLURM-specific resource parameters
    process.memory = { ... }
    process.time = { ... }
    process.cpus = { ... }
    
    // Limit job submission rate
    executor.queueSize = 100
    
    // Adjust resource limits for SLURM environment
    params.max_memory = 128.GB
    params.max_cpus = 32
    params.max_time = 72.h
}
```

## Customizing for Your Environment

If you're running the pipeline on a different Slurm cluster, you may need to customize the configuration:

1. **Partition/Queue**: Change `process.queue` to match a valid partition on your Slurm cluster
2. **Account**: Update `process.clusterOptions` with your Slurm account information
3. **Resource Limits**: Adjust `params.max_memory`, `params.max_cpus`, and `params.max_time` based on your cluster's limits

You can check available partitions on your Slurm cluster with:

```bash
sinfo
```

## Running with Slurm

To run the pipeline with Slurm, use the provided script:

```bash
./run_nextflow_slurm.sh
```

Or run Nextflow directly with the Slurm profile:

```bash
nextflow run main.nf -profile slurm [other parameters]
```

## Troubleshooting

If you encounter issues:

1. Check if the specified partition exists on your Slurm cluster
2. Verify that your account has access to the specified partition
3. Check Slurm job status with `squeue` and `sacct`
4. Review Slurm job logs in the output directory
5. Ensure that the Slurm daemons are running properly

For Docker-based Slurm, see the documentation in the `slurm-docker-cluster` directory.