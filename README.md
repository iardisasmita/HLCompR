# HLC_CompR
A computational pipeline to compare RNA sequencing (RNA-seq) data of hepatocyte *in vitro* models across studies. This pipeline involves read count normalization using DESeq2 variance-stabilizing transformation and quantile normalization.

HLC_CompR can be performed with or without an addition of new samples. HLC_CompR can also be performed either through a [web application]() or locally.

## Table of content

- [Addition of new sample](#addition-of-new-samples)
    - [Read mapping](#Read-mapping)
    - [CountData formatting](#CountData-formatting)
    - [ColData formatting](#ColData-formatting)
- [Running HLC_CompR web app](#Running-HLC_CompR-web-app)
- [Running HLC_CompR locally](#Running-HLC_CompR-locally)

## Addition of new samples
In order to compare RNA-seq data between studies, the data should be processed uniformly. Please use the following guideline to process additional samples to be included in the comparison analysis.

**NOTE**: It is recommended to take along PHH sample(s) to serve as a control for the presence of study-specific batch effects.

### Read mapping
#### Data upload to Galaxy
After quality controls have been performed, upload the FASTQ file(s) to the [Galaxy web platform](https://usegalaxy.eu/). For tutorials on how to upload data to Galaxy, see [here](https://galaxyproject.org/support/loading-data/) or [here](https://training.galaxyproject.org/training-material/topics/galaxy-interface/tutorials/upload-rules/tutorial.html).

Additionally, upload the reference genome and anotation that will be used for mapping. To ensure uniform processing, please use the files below:

| Item | Description | Link |
| --- | --- | ---|
| Reference Genome | Gencode_human_release_33_GRCh38.p13.genome.fa.gz | [Download](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/GRCh38.p13.genome.fa.gz) |
| Annotation | Gencode_human_release_33_gencode.v33.annotation.gtf.gz | [Download](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.annotation.gtf.gz)

**NOTE** : The top five rows of the annotation file need to be removed before it can be used in the RNA STAR tool.

#### Mapping with RNA STAR tool
Perform mapping of the FASTQ file(s) using the RNA STAR tool (Galaxy Version 2.7.2b). See figures below for the setting of RNA STAR tool.

![Figure_1](https://github.com/iardisasmita/HLC_CompR/blob/main/Images/fig1.jpg)

The mapping will take some time and after it is done, four files will be created for each sample.

![Figure_2](https://github.com/iardisasmita/HLC_CompR/blob/main/Images/fig2.jpg)

Download the read count file(s) by clicking the download button.

![Figure_3](https://github.com/iardisasmita/HLC_CompR/blob/main/Images/fig3.jpg)

### CountData formatting
Read counts obtained from RNA STAR tool will have one Ensembl gene ID column and three read count columns as shown below.

![Figure_4](https://github.com/iardisasmita/HLC_CompR/blob/main/Images/fig4.jpg)

Choose the correct column depending on the library preparation protocol (strandedness) and format them as shown below.

![Figure_5](https://github.com/iardisasmita/HLC_CompR/blob/main/Images/fig5.jpg)

This will be the CountData for the HLC_CompR web app.

### ColData formatting
The ColData contains descriptions of the sample(s) in the CountData. The format of the ColData is shown below.

![Figure_6](https://github.com/iardisasmita/HLC_CompR/blob/main/Images/fig6.jpg)

## Running HLC_CompR web app 

## Running HLC_CompR locally
