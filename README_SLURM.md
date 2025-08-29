<p><a target="_blank" href="https://app.eraser.io/workspace/jyA6Y2LLi1pxZeg0U25k" id="edit-in-eraser-github-link"><img alt="Edit in Eraser" src="https://firebasestorage.googleapis.com/v0/b/second-petal-295822.appspot.com/o/images%2Fgithub%2FOpen%20in%20Eraser.svg?alt=media&amp;token=968381c8-a7e7-472a-8ed6-4a6626da5501"></a></p>

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
    // Memory limits are set low for Docker-based Slurm compatibility
    params.max_memory = 900.MB
    params.max_cpus = 2
    params.max_time = 72.h
}
```
## Customizing for Your Environment
If you're running the pipeline on a different Slurm cluster, you may need to customize the configuration:

1. **Partition/Queue**: Change `process.queue`  to match a valid partition on your Slurm cluster
2. **Account**: Update `process.clusterOptions`  with your Slurm account information
3. **Resource Limits**: Adjust `params.max_memory` , `params.max_cpus` , and `params.max_time`  based on your cluster's limits
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
## Memory Limitations in Docker-based Slurm
The Docker-based Slurm cluster provided with this repository has very limited resources:

- Each node has only 1000 MB (1 GB) of memory (`RealMemory=1000`  in slurm.conf)
- Each node has only 1 CPU available (not explicitly defined in slurm.conf, so defaults to 1)
- The default memory per CPU is 500 MB (`DefMemPerCPU=500`  in slurm.conf)
- There are only 2 compute nodes available (`c[1-2]`  in slurm.conf)
Due to these limitations, the memory settings in the Slurm profile have been significantly reduced:

- Maximum memory per job: 900 MB (reduced from 128 GB)
- Maximum CPUs per job: 2 (reduced from 32)
- Default process memory: 500 MB per task attempt (reduced from 6 GB)
If you're running on a production Slurm cluster with more resources, you may want to increase these limits in `nextflow.config`.

## Troubleshooting
If you encounter issues:

1. Check if the specified partition exists on your Slurm cluster
2. Verify that your account has access to the specified partition
3. Check Slurm job status with `squeue`  and `sacct` 
4. Review Slurm job logs in the output directory
5. Ensure that the Slurm daemons are running properly
6. If you see "Memory specification can not be satisfied" errors, check the memory limits in your Slurm cluster and adjust the memory settings in `nextflow.config`  accordingly
7. If you see "CPU count per node can not be satisfied" errors, check the CPU limits in your Slurm cluster and adjust the CPU settings in `nextflow.config`  accordingly. For the Docker-based Slurm cluster, each node has only 1 CPU available.
For Docker-based Slurm, see the documentation in the `slurm-docker-cluster` directory.


<!-- eraser-additional-content -->
## Diagrams
<!-- eraser-additional-files -->
<a href="/README_SLURM-Methylation Sequencing Repository Flowchart-1.eraserdiagram" data-element-id="VDYutNPCpy0OVF7UcFOXf"><img src="/.eraser/jyA6Y2LLi1pxZeg0U25k___Y8ljVFLrIvgxRnuhvSMgSjkCJqt1___---diagram----29d0891a0647b62d24ce3ffcb66965d4-Methylation-Sequencing-Repository-Flowchart.png" alt="" data-element-id="VDYutNPCpy0OVF7UcFOXf" /></a>
<!-- end-eraser-additional-files -->
<!-- end-eraser-additional-content -->
<!--- Eraser file: https://app.eraser.io/workspace/jyA6Y2LLi1pxZeg0U25k --->