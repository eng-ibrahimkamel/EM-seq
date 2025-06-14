This repository contains tools and data related to [Enzymatic Methylation Sequencing](https://www.neb.com/products/e7120-nebnext-enzymatic-methyl-seq-kit) and [Enzymatic 5hmC-seq (E5hmC-seq)](https://www.neb.com/en-us/products/e3350nebnext-enzymatic-methyl-seq-5hmc-kit)

There are 3 [nextflow](https://www.nextflow.io/) scripts:
 - em-seq.nf (to align reads, filter and call methylation)
 - bins.nf (to calculate binned coverage around the TSS
 - cov_vs_meth.nf (to generate the "coverage by feature type" figure from the [EM-seq paper](https://genome.cshlp.org/content/31/7/1280))

Reference genomes containing spike-in methylation controls are available via an amazon s3 bucket: s3://neb-em-seq-sra/
 - GRCh38: https://neb-em-seq-sra.s3.amazonaws.com/grch38_core%2Bbs_controls.fa
 - T2T chm13 (hs1): https://neb-em-seq-sra.s3.amazonaws.com/T2T_chm13v2.0%2Bbs_controls.fa

## Execution Environments

The pipeline supports both local and SLURM execution environments through Nextflow profiles:

### Local Execution (Default)

```bash
nextflow run main.nf --input_glob '*_R1.fastq*' --path_to_genome_fasta /path/to/genome.fa --email your.email@example.com
```

### SLURM Execution

To run on a SLURM cluster, use the `slurm` profile:

```bash
nextflow run main.nf --input_glob '*_R1.fastq*' --path_to_genome_fasta /path/to/genome.fa --email your.email@example.com -profile slurm
```

You may need to customize the SLURM settings in `nextflow.config` to match your cluster configuration:

```groovy
// In nextflow.config
profiles {
    slurm {
        process.executor = 'slurm'
        process.queue = 'your_queue'  // Change to your SLURM queue/partition
        process.clusterOptions = '--account=your_account'  // Change to your SLURM account
    }
}
```

To use the Nextflow v1 scripts in this repository you need an older version of nextflow. 
```
NXF_VER=22.10.4 nextflow run em-seq.nf --genome em-seq_ref_files/T2T_chm13v2.0+bs_controls.fa --flowcell AAC27FDF --fastq_glob '*_R{1,2}.fastq*' -resume
```
We hope to upgrade to the NextFlow v2 syntax in the future.

You may also be interested in the [nf-core methylseq project](https://nf-co.re/methylseq/2.5.0)
