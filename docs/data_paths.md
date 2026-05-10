# Data Paths and File Structure

## HPC Location
All data resides on UQ Bunya HPC under:
Edward
## Core Analysis Files

| File | Path | Description |
|------|------|-------------|
| Sample metadata | `sra_ww_mobilization/results/sample_map_complete.tsv` | 185 samples, country/year/category |
| RGI + MAG taxonomy | `sra_ww_mobilization/results/arg_analysis/rgi_with_mag_taxonomy.csv` | Core integrated dataset |
| Bracken genus abundance | `sra_ww_mobilization/results/bracken_genus_abundance.csv` | Community profiles |
| ARGs-OAP abundance | `sra_ww_mobilization/results/arg_analysis/argoap_type_abundance.csv` | Normalised ARG copies/16S |

## Per-Sample Data Structure
/QRISdata/Q6636/data/{category}/{sample_id}/
├── clean_fastp/       fastp QC reads
├── assemblies/        MEGAHIT contigs
├── kraken2/           Kraken2 classification
├── bracken/           Bracken abundance
├── rgi/               RGI ARG annotation
├── argoap/            ARGs-OAP quantification
├── genomad/           geNomad MGE detection
├── binning/           MetaBAT2 MAGs
├── checkm2/           CheckM2 quality reports
└── gtdbtk/            GTDB-Tk taxonomy
## Sample Categories
- `ww_influent_municipal` (n=108)
- `ww_effluent_municipal` (n=61)
- `ww_influent_hospital` (n=3)
- `ww_effluent_hospital` (n=2)
- `ww_sludge` (n=9)
- `ww_AU_municipal` (n=2)

## Key Database Paths
| Database | Path |
|----------|------|
| CARD (RGI) | Built into RGI conda env |
| geNomad | `/scratch/user/uqnzhai/genomad_db/genomad_db` |
| CheckM2 | `/QRISdata/Q6636/db/checkm2_db_fresh/CheckM2_database/` |
| GTDB-Tk r220 | `/scratch/opendata/genomics/GTDB/releases/release220/` |
| Kraken2 | `/QRISdata/Q6636/metawrap/kraken2_db` |
