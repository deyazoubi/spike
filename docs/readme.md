# Spike documentation
### For all users
 1. [Access pipeline data through SFTP](access.md): Most data of the pipeline is made available to you via the web page [Snupy](https://snupy.hpc.rz.uni-duesseldorf.de/). However, you have the possibility to get *all* data via SFTP, but you need to sign up for this first. Useful to get the "bam" files or CNV plots.

### For admins
 2. [The Samplesheet](samplesheet.md): this table is the most valuable piece of information, since it stored all you need to know about all samples ever processed with the pipeline.
   1. First, you should get a little bit familiar with a Pandas DataFrame. Here is a very good and [short introduction](https://pandas.pydata.org/pandas-docs/stable/getting_started/10min.html).
 3. [Snakemake](hhu.md): Spike is implemented via the Snakemake framework. You should read our documentation about the used command line flags *before* you execute Spike.
 4. [Process flowcells @UKD](walkthrough_hhu.md): gives detailed instruction and best practises on how to routinely process new flowcell data at the UKD.
 5. Create a new [project](new_project.md).
