## JLU
conda activate phylofiller
snakemake -p  --cluster-config snakemake/cluster_jlu.json --cluster "qsub -A {cluster.account} -q {cluster.queue} -l mem={cluster.mem} -l walltime={cluster.time},nodes={cluster.nodes}:ppn={cluster.ppn}" -j 100 --latency-wait 900 --use-conda --cluster-status phylofiller/qstat_jlu.py --max-status-checks-per-second 1 --keep-going -n

## HHU
conda activate phylofiller
snakemake -p -s snakemake/Snakefile --cluster-config snakemake/cluster_hhu.json --cluster "qsub -A {cluster.account} -q {cluster.queue} -l mem={cluster.mem} -l walltime={cluster.time},nodes={cluster.nodes}:ppn={cluster.ppn}" -j 100 --latency-wait 900 --use-conda --cluster-status phylofiller/qstat_hhu.py --max-status-checks-per-second 1 --keep-going /gpfs/project/jansses/PhyloFiller/fungi.complete -n
