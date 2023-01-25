# Synechocystis RNAlocSeq Analysis

This repository contains a pipeline to analyze RNAloc-Seq Data from Synechocystis.
The Pipeline can be run automatically via installing the conda environment from `env.yml`

```shell
conda env create -f env.yml
```

If you have the `Data` folder you can run the pipeline via

> **Warning**
> The Data is not publically available yet


```shell
snakemake -s deseqpipeline.smk --cores 1 -k
```

This produces the `PipelineData` directory where you can find all important files.

The files within this directory contain wildcards, which map the to the corresponding
analysis pipeline. The following Table describes the abbreviations used in file names:

| Abbreviation | Values                           | Description                                                                                                                                                                                                                                                        |
|--------------|----------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| set          | Light, LightDark, LightDarkLight | The Dataset used.  Light contains everything but the Puromycin columns from the Light Dataset. LightDark represents the Light part of the LightDark Dataset LightDarkLight is the merged Dataset containing all Light columns from the Light and LightDark Dataset |
| design       | rf                               | A mapping mapping the abbreviation to the design formula used in DESeq2. rf: "~ replicate + fraction"                                                                                                                                                              |
| cond         | M, TC, C                         | Which fraction to compare: M: Membrane C: Cytosol TC: Total Cell                                                                                                                                                                                                   |
| base         | M, TC, C                         | The base to compare cond to. See cond for description                                                                                                                                                                                                              |
| ud           | upregulated, downregulated       | Only present in Enrichment files. Shows whether the file shows up or downregulated enriched set.                                                                                                                                                                   |

