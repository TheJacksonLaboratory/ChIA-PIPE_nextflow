# ChIA-PIPE: A Fully Automated Nextflow Pipeline for Comprehensive ChIA-PET Data Processing and Analysis

This [Nextflow](https://www.nextflow.io/) pipeline is available on JAX Github:
https://github.com/TheJacksonLaboratory/ChIA-PIPE_nextflow

Run [ChIA-PIPE_nextflow](https://github.com/TheJacksonLaboratory/ChIA-PIPE_nextflow) on paired-end ChIA-PET data.

This is adapted to Sumner cluster (Slurm, Singularity) and uses the jaxreg.jax.org containers repository. Follow the User Guide below to add it to your library to enable you running the pipeline directly using the containers from jaxreg.

User Guide:
([https://jacksonlaboratory.sharepoint.com/sites/ResearchIT/SitePages/JAX-Singularity-Container-Registry-User-Guide.aspx](https://jacksonlaboratory.sharepoint.com/sites/ResearchIT/SitePages/JAX-Singularity-Container-Registry-User-Guide.aspx))



## To test the pipeline: install Nextflow, ChIA-PIPE_nextflow and run


### [Install Nextflow](https://www.nextflow.io/index.html#GetStarted)

Install Nextflow with the following command:

```bash
wget -qO- https://get.nextflow.io | bash
```

This will create the nextflow main executable file in the current directory. Make the binary executable on your system by running `chmod +x nextflow` and move the nextflow file to a directory accessible by your $PATH. For example: `/home/$USER/bin/nextflow`


### Install ChIA-PIPE_nextflow

```
git clone https://github.com/TheJacksonLaboratory/ChIA-PIPE_nextflow.git
```

### Run ChIA-PIPE_nextflow

You must create a parameter file before running the pipeline. An example of parameter file [(can be found here)](https://github.com/TheJacksonLaboratory/ChIA-PIPE_nextflow/blob/main/params.config).

params.config file contains at least the parameter entries:

```
params {
   fastqInputs = [path with comma separated R1 and R2]
   genome      = [name] 
   bwaIndex    = [path]
   run         = [name]
   runType     = [name]
   ipFactor    = [name]
   cellType    = [name]
   inputCtrl   = [path]
   blackList   = [path]
   peakCaller  = [name]
   outputDir   = [path]
}
```

Plus any other optional parameters, which can be shown by passing the --help argument.

```   
nextflow run /path/to/installed/ChIA-PIPE_nextflow --help
```

DO NOT use the full input fastq file names like: 
LHG0146_GT21-16861_CTCTCTAC-AAGGAGTA_S1_R1_001.fastq.gz
LHG0146_GT21-16861_CTCTCTAC-AAGGAGTA_S1_R2_001.fastq.gz


Make sure the input files are in this format: 
LibraryID_*R{1,2}*.fastq.gz

E.g.
LHG0146_R1_001.fastq.gz
LHG0146_R2_001.fastq.gz


Make the symbolic links for the full input fastq file names with the following commands. Here LHG0146 is a LibraryID.

ln -s LHG0146_GT21-16861_CTCTCTAC-AAGGAGTA_S1_R1_001.fastq.gz LHG0146_R1_001.fastq.gz
ln -s LHG0146_GT21-16861_CTCTCTAC-AAGGAGTA_S1_R2_001.fastq.gz LHG0146_R2_001.fastq.gz


This is an example of fastqInputs (Comma delimited).
fastqInputs = "${_fqPath}/LHG0146_R1_001.fastq.gz,${_fqPath}/LHG0146_R2_001.fastq.gz"


You can use an example parameter file [(can be found here)](https://github.com/TheJacksonLaboratory/ChIA-PIPE_nextflow/blob/main/params.config) for testing.



#### Launch Pipeline with SLURM Wrapper Script::

The pipeline can be launched using a wrapper script [(provided here)](https://github.com/TheJacksonLaboratory/ChIA-PIPE_nextflow/blob/main/submit_chiapipe.sh).

Copy params.config and wrapper script to your run directory, modify/edit them and execute the following command.

```bash
sbatch submit_chiapipe.sh
```


### Output Files

The pipeline output files will be located in the outputDir within your run directory.
For example: `ChIA-PIPE_nxf_output/*`



**NOTE 1: You must specify a location for the pipeline in the wrapper script.**
`export projectDir=/path/to/installed_ChIA-PIPE_nextflow_pipeline`

**NOTE 2: You must make the symbolic links for the full input fastq file names as shown above.**

**NOTE 3: You must set up your access to [The Jackson Laboratory Container Registry](https://jaxreg.jax.org/) (jaxreg) as per above User Guide.**
