## Dataset and Data Availability
This project uses publicly available 16S rRNA gene sequencing data and associated sample metadata from the NCBI BioProject PRJNA622999, titled Gut Microbiota of women with Polycystic Ovary Syndrome. The original study investigated gut microbial composition in women with polycystic ovary syndrome (PCOS) and healthy controls.

### Dataset Details:
- Repository: NCBI BioProject / Sequence Read Archive (SRA)
- BioProject accession: PRJNA622999
- Study population: 20 women with PCOS and 20 matched healthy controls
- Sample type: Fecal samples
- Sequencing method: 16S rRNA gene sequencing
- Target region: V4 region
- Sequencing data: 40 SRA experiments associated with 40 BioSamples

The original dataset remains available through NCBI. This repository contains the analysis code and selected outputs generated from the publicly available data; it does not claim ownership of the original sequencing data or metadata.
Please refer to the original BioProject record and associated publication for details of the study design, sample collection, and data generation.

To run the pipeline with your data, provide the following resources locally:

### Required Input
- **Raw sequencing reads:** Place paired end FASTQ files in data/raw_data/. Expected input includes paired forward and reverse reads for each sample. 
- **Sample metadata:** Provide a metadata file containing sample identifiers and the variables required for downstream analyses. The sample identifiers must match the sequencing sample names used in the phyloseq object.
- **SILVA reference files:** This pipeline uses the SILVA 16S rRNA reference database (version 138.1).
  Download the following files and place them in data/reference/:
  - silva_nr99_v138.1_train_set.fa.gz
  - silva_species_assignment_v138.1.fa.gz

Reference database: https://www.arb-silva.de/

 ## Data Availability
This repository does not redistribute restricted or potentially sensitive study data. Users wishing to reproduce the workflow should use appropriately authorized data or publicly available 16S rRNA sequencing datasets with a compatible metadata structure.
