# Hypertension, Diabetes and HRQoL in Korea

Standalone SAS programs examining associations of hypertension, diabetes and treatment factors with the five EQ-5D-3L dimensions in two Korean population surveys.

## Programs

- **KCHS_ALL.sas**: Korea Community Health Survey analysis.
- **KNHANES_ALL.sas**: Korea National Health and Nutrition Examination Survey analysis.

Both programs include data preparation, baseline characteristics, regional and annual summaries, primary associations, five additional-covariate sensitivity analyses, fixed 90%/70%/50% subsamples and exploratory treatment analyses. Regression models account for survey weights, strata and primary sampling units.

## Requirements

- SAS 9.4 with SAS/STAT; use a UTF-8 session.
- Study-specific merged datasets: `KCHS_Data.sas7bdat` and `KNHANES_Data.sas7bdat`.
- The 13 annual KCHS SPSS files used to reconstruct sampling identifiers.
- Saved survey-specific subsample IDs referenced in each program.

These data inputs are not bundled with the scripts.

## Usage

1. Update the project, input and output paths at the beginning of the selected program. Match the input filenames and directory layout to your local files; encoded path expressions preserve the original Unicode directory names.
2. Open the program in SAS and submit the entire file. Each program contains its own statistical procedures and requires no external SAS include files.
3. Review the generated log, HTML results, SAS datasets and CSV summaries. Each run writes to a separate output directory.
