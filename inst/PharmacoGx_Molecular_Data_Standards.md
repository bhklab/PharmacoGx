---
output:
  pdf_document: default
  html_document: default
---
# Stardardized Molecular Data names and annotations for PSets

This document will describe two things: First, we document the currently required standards for naming, annotating, and including molecular data into psets.\
Second, we will describe the desired future state.

This is designed to be a living document, which can be updated as things change, progress is made on the "todo" list, and new ideas are added.

## Current necessary considerations in molecular profiles for PSets

There are a few things that are necessary for proper functioning of PSets in PharmacoGx. Here, we focus on the requirements from the molecularData slot only.

### Molecular Data Names

This is what is returned when you run `mDataNames` on a PSet. It is also the names of the elements of the list making up the `@molecularData` slot.

Currently, there is a mixed bag of molecular data names used across psets, as no fixed vocabulary is required here. Nothing is enforced, and only the "show" function in pharmacogx assumes certain names, and it handles missing entries (for example, missing profile types, or even no molecular profiles) gracefully.

Currently, we meet the following names in PSets downloaded from orcestra:

| Name                            | Data type                                                    |
|---------------------------------|--------------------------------------------------------------|
| rnaseq                          | rnaseq (gene tpm??)                                          |
| rna                             | microarray rna                                               |
| Kallisto_0.46.1.rnaseq          | rnaseq gene level tpm                                        |
| Kallisto_0.46.1.rnaseq.counts   | rnaseq gene level counts                                     |
| Kallisto_0.46.1.isoforms        | rnaseq isoform level tpm                                     |
| Kallisto_0.46.1.isoforms.counts | rnaseq isoform level counts                                  |
| cnv                             | snp array derived copy number                                |
| mutation                        | mutation data, at protein level                              |
| fusion                          | presence of fusion data, as binary                           |
| mutationall                     | mutation data, at protein level                              |
| mutationchosen                  | mutation data, at protein level                              |
| mutation_exome                  | mutation data, at protein level, specifically from exome seq |


Of these, `show` is aware of:

-   rna
-   rnaseq
-   dna
-   snp
-   cnv

Clearly, this is outdated and should be fixed.

### Annotation of SummarizedExperiments

PharmacoSets require a specific entry in the metadata of a SummarizedExperiment to be set to a valid value to work with the functions in the package.

Specifically, we require `metadata(SE)$annotation` to be set to one of the valid values described below.

Functions within the package use the value of this field to determine how to treat the data within the SummarizedExperiment (is this gene expression, copy number, mutation, etc)

The currently allowed vocabulary is:


| Name                            | Data type                                                                                                                                                 |
|---------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------|
| rna                             | gene expression values                                                                                                                                    |
| cnv                             | logR copy number                                                                                                                                          |
| mutation                        | gene level binary mutation, encoded as "wt" for 0, or the protein change for a 1. Multiple mutations in same gene are separated with '///' .              |
| fusion                          | gene partner level binary fusion status, encoded as "wt" for 0, or the fusion partners for a 1. Multiple fusions in same genes are separated with '///' . |


### Annotation of the right version for the PSet

Finally, if in the creation of your PSet, you include SummarizedExperiments and not the older ExpressionSet objects, you must set an annotation on the PSet to signify this. The original PSet used ExpressionSets.

To do this, set `annotation(PSet)$version <- 2` or larger.

## Future Changes

This section is TBD. However, some considerations:

We should have a more standardized set of annotations for each SE. Proposed additions:

- pipeline: describing the pipeline
- transformation: any transformation applied to the raw output of the pipeline


We should support methylation data soon.
