# RSNA_kaggle_CT_wrangle

Download the "Stage 2" data from 
https://www.kaggle.com/c/rsna-intracranial-hemorrhage-detection/data

## DCM -> NifTI

1) Uncompress the data (`.dcm` slices) into a folder, e.g. `data_dir = '/mnt/wwn-0x5000c500b6b1dabd/rsna-intracranial-hemorrhage-detection'` in `wrangle_RSNA_Kaggle.R`
2) Run `wrangle_RSNA_Kaggle(data_dir)` with the proper `data_dir`

### `dcm2niix`?

Will fail with default settings as the _z_-distances are not constant per volume, so manual slice sorting has been done and NifTIs are saved with basic metadata, i.e. voxel dimensions

## HU Range Check

Most of the `.dcm`s had air set at 0 HU (and the noise above it) and they are corrected to "correct" value of -1024 (plus the noise above it). You could re-check this, but most of the air values should be at and above -1000. This util reads in the HU intensity histograms saved during first pass.

## `kaggle_check_headers.R`

You can after the processing check if the headers were corrected actually straight from the intermediate `.RData`s which is faster disk I/O-wise then re-reading all the .dcms

