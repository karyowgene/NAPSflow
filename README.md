# NAPSflow: a snakemake pipeline for neoantigen prediction

A Snakemake pipeline to predict breast cancer neoantigens from WES and RNA data of matched tumor and normal tissue samples \
The workflow can be run as a single pipeline, but each of the data processing steps, are also stand along snakemake pipelines \
This allows those who just have WES data to process it upto vcf annotation and do some other downstream analysis, those in need to \
RNA seq analysis to determine expression levels from their data. \
The pipeline is under active development, please feel free to contribute. 

To run the pipeline, you need to edit the `config.yaml`\
Then: \
`git clone https://github.com/karyowgene/NAPSflow.git` \
`cd NAPSflow` \
`conda create -n napsflow-env env.yaml` \
`conda activate napsflow-env` \
`cd workflow` \
`snakemake --cores 8` 

The pipeline was used in the analysis of data published in this paper: (https://pmc.ncbi.nlm.nih.gov/articles/PMC11668681/)