## Data
Raw sequencing data and metadata not included in this repository. To run the pipeline, provide the following resources locally:

### Required Input
- **Raw sequencing reads:** Place paired end FASTQ files in data/raw_data/. Expected input includes paired forward and reverse reads for each sample.
- **Sample metadata:** Provide a metadata file containing sample identifiers and the variables required for downstream analyses. The sample identifiers must match the sequencing sample names used in the phyloseq object.
- **SILVA reference files:** Download the SILVA reference files required for DADA2 taxonomic assignment and place them in data/reference/. This pipeline uses:
  - SILVA training set for taxonomic assignment
  - SILVA species assignment reference database

 ## Data Availability
This repository does not redistribute restricted or potentially sensitive study data. Users wishing to reproduce the workflow should use appropriately authorized data or publicly available 16S rRNA sequencing datasets with a compatible metadata structure.
