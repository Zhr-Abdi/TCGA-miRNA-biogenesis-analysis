# cBioPortal alteration frequency raw data

This folder contains raw alteration frequency files downloaded from cBioPortal.

Each file corresponds to one miRNA biogenesis-related gene and contains alteration frequency data across selected TCGA Firehose Legacy cancer studies.

Columns in the raw files:

- Cancer Study
- Alteration Frequency
- Alteration Type
- Alteration Count

Alteration types include:

- amp: amplification
- homdel: deep deletion
- mutated: mutation
- multiple: multiple alterations

Rows with zero frequency may be absent from the original cBioPortal downloads and will be added during data processing.
