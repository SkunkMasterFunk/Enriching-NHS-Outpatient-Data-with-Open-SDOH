# Beyond Deprivation Indices: Enriching NHS Outpatient Data with Open Social Determinants of Health to Predict Missed Appointments

Analysis code for a study of missed outpatient appointments ("did not attend", or
DNA) across specialties at Cambridge University Hospitals, linking hospital
referral records to publicly available small-area measures of deprivation,
transport access, and other social determinants of health.

The analysis is split across two files by language. `DNA_analysis_pipeline.R`
carries the cleaning, linkage, descriptive statistics, regression modelling, and
mapping; `DNA_analysis_pipeline.py` carries the geographic linkage, the distance
calculations, and the two machine learning analyses. The two interleave, so the
running order below is essential.

---

## 1. Repository layout

Both scripts assume they are run from the project root, and every path in them is
relative to it. The expected structure is:

```
project_root/
├── DNA_analysis_pipeline.R          analysis pipeline (R)
├── DNA_analysis_pipeline.py         analysis pipeline (Python)
├── README.md
│
├── data/                            raw specialty extracts (.xlsx, not shared)
├── open_data/                       publicly available source files
│   ├── ahah_v5.csv
│   ├── iod_2025.csv
│   ├── NAPTAN_2022.csv
│   ├── vehicle_availability_2021.csv
│   ├── childcare_accessibility_2023.csv
│   ├── digital_propensity_index_2021.csv
│   ├── marital_status_2021.csv
│   ├── ethnicity_LSOA_21.csv
│   ├── lsoa_pop.csv
│   ├── LSOA shapefile.geojson
│   └── multi_csv/                   ONS postcode directory, one file per area
│
├── clean_data/                      created by the pipeline
└── outputs/                         created by the pipeline
```

`clean_data/` and `outputs/` are created automatically on first run and need not
exist in advance.

---

## 2. Running order

The pipeline runs in seven stages. Stages 3 and 4 are Python and sit in the
middle of the R script, which is why the R file is not simply sourced end to end.

| Stage | What happens | Run |
|---|---|---|
| 1 | Combine all referral data into one aggregate file | R, Section 1 |
| 2 | Clean the individual public datasets | R, Section 2 |
| 3 | Link LSOA and lat/long to postcode | `python DNA_analysis_pipeline.py postcode-link` |
| 4 | Calculate haversine distances | `python DNA_analysis_pipeline.py distances` |
| 5 | Link public data to the aggregate dataset | R, Section 3 |
| 6 | Derived variables, missingness, descriptive statistics | R, Sections 4 to 6 |
| 7 | Statistical analysis and heatmaps | R, Sections 7 and 8 |
| 7b | Clustering | `python DNA_analysis_pipeline.py clustering` |

In practice this means running the R script down to the marked stop after Section
2, switching to Python for the two geography steps, then returning to R for
everything from Section 3 onwards, and finally running the clustering
against the complete case file that R Section 3 produces. The R script carries a
conspicuous banner at the point where it expects the Python steps to have run, so
the handover is difficult to miss.

The Python steps can be run individually, as above, or chained with
`python DNA_analysis_pipeline.py all`, though the latter is only correct if the R
script has already been taken through Section 3, since the clustering reads the
complete case file.

There is one further optional step,
`python DNA_analysis_pipeline.py verify-postcodes`, which re-runs the postcode
lookup independently and writes a list of rows whose stored geography disagrees
with the reference file. It is slow and is not part of `all`; it is worth running
once after any change to the lookup source, and otherwise not at all.

---

## 3. Dependencies

### R

```r
install.packages(c(
  "dplyr", "tidyverse", "readxl", "ggplot2", "scales",
  "splines", "lmtest", "car",
  "grid", "gridExtra",
  "sf", "ggspatial", "cowplot"
))
```

`splines` and `grid` ship with R and need no installation. `sf` requires the
system libraries GDAL, GEOS, and PROJ. These may be installed natively.
If not, Debian or Ubuntu these come from
`libgdal-dev libgeos-dev libproj-dev`, and on macOS from `brew install gdal geos proj`.

### Python

```bash
pip install pandas numpy scikit-learn scipy statsmodels matplotlib seaborn hdbscan
```

`hdbscan` is needed only for the clustering, which degrades gracefully if it is
absent: the step
prints a notice and continues with the remaining four algorithms, and recent
versions of scikit-learn ship their own HDBSCAN implementation which the script
will use in preference if it is available.

---

## 4. Data availability

The hospital referral extracts in `data/` are patient level records and are **not
included in this repository**. They cannot be shared publicly. Everything in
`open_data/` is publicly available and can be downloaded from the sources listed
in the manuscript; the file names above are the ones the scripts expect, so
downloaded files may need renaming.

The ONS postcode directory in `open_data/multi_csv/` is likewise not included, as
it is large and separately licensed. The scripts expect it split by postcode
area, one CSV per area, each carrying `pcd`, `lat`, `long`, and `lsoa21` columns.
File names are matched on their trailing area code, so `ONSPD_CB.csv` is where
Cambridge postcodes are looked up; the naming otherwise does not matter. Where a
directory uses some other convention, the indexer falls back to the leading
letters of the file name. If the directory is missing altogether, the
`postcode-link` step raises an explicit error pointing back to this section
rather than failing obscurely further down.

**How the lookup works.** Each postcode is parsed for its area prefix (`CB1 2AB`
gives `CB`), and only that area's file is opened to find the matching row and its
coordinates. Because the same area recurs across thousands of patients, the
cohort's postcodes are grouped by area first, so each file is opened exactly once
and files for areas nobody lives in are never opened at all. Scanning stops early
once every postcode wanted from a file has been found. This keeps memory
proportional to the size of the cohort rather than to the ONS directory, which is
several gigabytes in full.

The step reports its matching per area, so a systematically unmatched area is
visible immediately:

```
  CB ONSPD_CB.csv       12,431 / 12,433 matched
  PE ONSPD_PE.csv        3,092 /  3,092 matched
  no lookup file for areas: QQ
```

Postcodes that are blank, malformed, or belong to an area with no lookup file are
retained with NA geography rather than dropped, so that the count of rows lost to
linkage failure stays visible and is reported explicitly by R Section 3.

---

## 5. Pipeline outputs

### Intermediate data (`clean_data/`)

The analysis dataset passes through five named stages, and both scripts define
these names as constants at the top so that the two languages cannot disagree
about them:

| Constant | File | Created by |
|---|---|---|
| `F_AGG` | `aggregate_specialty_data_24.csv` | R Section 1 |
| `F_AGG_PCD` | `..._pcd_link.csv` | Python `postcode-link` |
| `F_AGG_NAPTAN` | `..._pcd_NAPTAN_link.csv` | Python `distances` |
| `F_AGG_LINKED` | `..._linked.csv` | R Section 3 |
| `F_AGG_CCA` | `..._linked_CCA.csv` | R Section 3 |

`F_AGG_LINKED` retains every row, including those whose linkage failed, and is
the file used for the missingness analysis. `F_AGG_CCA` is the complete case
dataset and is what every analysis from Section 4 onwards reads.

Alongside these, Section 2 writes one cleaned file per public dataset, and the
Python distance step writes `add_dist.csv`, which is keyed on postcode rather
than on LSOA.

### Results

Descriptive tables, model summaries, odds ratio tables, predicted probabilities,
likelihood ratio tests, and the model comparison are all **printed to the
terminal**. To keep a copy of a run, redirect it:

```bash
Rscript DNA_analysis_pipeline.R > outputs/analysis_log.txt 2>&1
```

The only files written to `outputs/` are those that cannot be read as text:

```
outputs/
├── map_dna_aggregate.png                DNA rate by LSOA, England
├── map_patients_aggregate.png           patient distribution by LSOA, England
├── map_dna_by_specialty*.png            DNA rate, Cambridgeshire, faceted
├── map_patients_by_specialty*.png       patient distribution, faceted
├── map_bivariate_health.png             deprivation against vehicle ownership
├── map_bivariate_employment.png
├── map_bivariate_crime.png
├── distance_decay.png
└── clustering/                          comparison tables, profiles, OR tables, figures
```

The specialty maps are split into pages of twelve where there are more
specialties than fit legibly on one grid, hence the `*` in those names.

---

## 6. Analytical notes

A few decisions are embedded in the code and are worth stating plainly, since
they affect how the results should be read.

**Three analysis samples, used consistently.** Aggregate analyses use one row per
patient, taking their first referral, so that tests of independence hold.
Specialty-stratified analyses use one row per patient per specialty, so that a
patient referred to two services contributes to both without contributing twice
to either. Repeat non-attendance is examined on the full appointment level
sample. Section 4 of the R script builds all three, and everything downstream
draws from them rather than re-deriving its own.

**Missed appointments include patient cancellations.** An appointment cancelled
by the patient, or one they could not attend, is counted as a DNA; cancellations
by the service are excluded entirely, on the grounds that they are not patient
behaviour. Only completed appointments and DNAs are retained.

**Recorded and linked ethnicity are different variables.** `ethnicity2` is the
individual's recorded ethnicity, with missing kept as its own level rather than
dropped, since ethnicity recording is itself patterned. `non_white_pct` is the
LSOA level percentage from the census. They measure different things and are
reported separately throughout.

**Deprivation is modelled with natural splines.** A linear term would impose a
constant effect per decile, which is not plausible, while a ten level factor
spends nine degrees of freedom to say what a smooth curve says in two or three.
The degrees of freedom are chosen by AIC rather than fixed in advance, and
whole-term significance is judged by likelihood ratio test. The individual spline
basis coefficients in the OR tables are **not** interpretable as the effect of a
one unit change; read those terms from the partial effect plots instead. The two
genuinely continuous covariates, vehicle ownership and childcare accessibility,
are tested for non-linearity rather than assumed to need a spline, and the
resulting choice is carried through every model that uses them.

**Distances are straight line, not routed.** Both the distance to Addenbrooke's
and the distance to the nearest transport node are great circle distances. This
understates real travel, and understates it more for rural patients than urban
ones. It is worth stating explicitly wherever these estimates are reported.

**Five clustering algorithms, not one.** Each makes different assumptions about
cluster shape, so structure that survives all five is more credible than
structure only one method finds. Agreement is quantified by the Adjusted Rand
Index, and the substantive question is settled by whether cluster membership
improves cross-validated AUC over the covariates alone. A good silhouette score
describes geometry, not clinical meaning, and a delta AUC below roughly 0.005
indicates the covariates had already captured that structure.

---

## 7. Changes made during consolidation

This code was consolidated from thirteen Jupyter notebooks. The reorganisation
was mostly mechanical, but three substantive corrections were made, and they are
recorded here rather than buried, since each changes results.

PDF generation was likewise dropped in favour of terminal output, the point of
this code being reproducibility rather than typesetting.

**The interpreter variable was always missing.** The cleaning function contained
a copy-paste error, testing `sex == "Y"` rather than `interpreter == "Y"`, so the
recode produced NA for every row. It is now correct, which means interpreter
requirement becomes a usable variable for the first time.

**The heatmap notebook referenced objects before defining them.** `lsoa_shapes`
and `bivar_lsoa` were both used above their definitions and only worked because
of out-of-order cell execution. The consolidated script defines everything before
use, and this has been checked mechanically across the whole file.

Three notebooks were fully superseded and are not represented separately:
`data_clean_single` was the single-specialty prototype of `data_clean_function`,
`heatmaps_clean` was superseded by `heatmaps_new` (which also reads the correct
complete case file), and the first clustering cell of `KMeans` was superseded by
the multi-algorithm cell that follows it. Beyond that, three near-identical
bivariate map blocks were reduced to one constructor, and roughly fifteen
repeated DNA cross-tabulation blocks to three helpers, which accounts for most of
the reduction in length. Neither refactor changes any output.
