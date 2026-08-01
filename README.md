# Housing Data Lab

**A modular, weighted, decision-support platform for identifying and ranking
Local Layer Super Output Areas (LSOAs) for land acquisition and builder
partnerships — built for the Homes England *Designing a Housing Data Lab*
brief (Cambridge Policy Hackathon 2026).**

The tool turns fragmented Indices of Deprivation (IoD) sub-domains, land-use
data, an affordability-derived "market failure" proxy, and (optionally)
flood-risk and IoD-2025 official scores into a single composite *Opportunity
Score* per LSOA — with every variable's weight tunable from a sidebar of
sliders and toggles.

```
  IoD sub-domains  ─┐
  Land-use pct/num  ├──►  Standardise (min-max → deciles/quintiles)
  Flood risk        │            │
  IoD-2025 main    ─┘            ▼
                            Wide LSOA frame
                                 │
                          ┌──────┴──────┐
                          ▼             ▼
                    GeoJSON join    Weighted-index engine
                          │             │
                          └──────┬──────┘
                                 ▼
                     Heatmap  +  Top-N rank table
                     (live, in browser, Dash UI)
```

---

## Quick start

```bash
# 1. Create + activate a virtual environment (recommended)
python -m venv .venv
source .venv/bin/activate          # on Windows: .venv\Scripts\activate

# 2. Install dependencies
pip install -r requirements.txt

# 3. Drop the source CSVs in data/
#    (see "Data" section below for the full schema)
cp /path/to/filtered_*_domain.csv data/
cp /path/to/filtered_land_use_*.csv data/

# 4. Optional: drop the LSOA boundary GeoJSON in too
#    (without it the app falls back to a ranked-bar visualisation)
cp /path/to/lsoa_boundaries.geojson data/

# 5. Smoke-test the data pipeline (recommended first run)
python tests/test_pipeline.py

# 6. Launch the dashboard
python main.py                     # http://127.0.0.1:8050
```

---

## Directory layout

```
housing_data_lab/
├── README.md                       # this file
├── requirements.txt
├── main.py                         # entry point — `python main.py`
│
├── config/
│   └── settings.py                 # paths, default weights, schema map
│                                   # — edit this to add a new covariate
│
├── data/                           # ← drop the source CSVs here
│   ├── filtered_income_domain.csv          (required)
│   ├── filtered_employment_domain.csv      (required)
│   ├── filtered_education_domain.csv       (required)
│   ├── filtered_health_domain.csv          (required)
│   ├── filtered_crime_domain.csv           (required)
│   ├── filtered_barriers_domain.csv        (required)
│   ├── filtered_living_env_domain.csv      (required)
│   ├── filtered_land_use_pct.csv           (recommended)
│   ├── filtered_land_use_number.csv        (recommended)
│   ├── iod2025_main_scores.csv             (optional)
│   ├── flood_risk_lsoa.csv                 (optional)
│   └── lsoa_boundaries.geojson             (optional — enables choropleth)
│
├── src/
│   ├── pipeline.py                 # ◀── orchestrator (one call → everything)
│   ├── data_processing/
│   │   ├── loader.py               # CSV ingest + LSOA-ID normalisation
│   │   ├── standardiser.py         # min-max + decile/quintile binning
│   │   ├── market_failure.py       # affordability-derived proxy
│   │   └── land_use.py             # standardisation + threshold filter
│   ├── spatial/
│   │   └── joiner.py               # GeoJSON join (graceful fallback)
│   ├── analysis/
│   │   └── weighted_index.py       # composite scoring engine
│   └── ui/
│       ├── app.py                  # Dash factory
│       ├── components.py           # layout / sliders / toggles
│       └── callbacks.py            # interactivity (re-score on slider move)
│
└── tests/
    └── test_pipeline.py            # end-to-end smoke test
```

The strict layering means each module has one job:

| Module | Inputs | Output | Used by |
|---|---|---|---|
| `loader.py`        | filesystem CSVs                            | `RawData` bundle of DataFrames | `pipeline.py` |
| `standardiser.py`  | raw sub-domain frames                       | wide frame with `*_scaled`, `*_decile` cols | `pipeline.py` |
| `market_failure.py`| Barriers frame                              | wide frame with `market_failure_*` cols | `pipeline.py` |
| `land_use.py`      | pct + number frames                         | wide frame with `*_pct__scaled` + filter helper | `pipeline.py`, `callbacks.py` |
| `joiner.py`        | wide tabular frame, GeoJSON                 | GeoDataFrame (or unchanged tabular) | `pipeline.py` |
| `weighted_index.py`| wide frame, weight dict, active-toggle dict| `ScoringResult` (scored frame + contribution matrix) | `callbacks.py` |
| `components.py`    | nothing (static layout)                     | Dash `html.Div` tree                | `app.py` |
| `callbacks.py`     | the pipeline product                        | dynamic `figure` + `DataTable`      | `app.py` |

---

## Data sources & schema contracts

### Required: IoD sub-domain CSVs

Seven files, one per Indices of Deprivation sub-domain. Each must have
`lsoa21` (or `LSOA21CD` / `LSOA Code`) as its first column and arbitrary
numeric sub-domain measures in subsequent columns. The standardiser auto-
detects every non-ID column and produces matching `*_scaled` (continuous
[0,1]) and `*_decile` (integer 1..10) columns.

If you add a new sub-domain (say, "transport accessibility"), append one
line to `SUBDOMAIN_FILES` in `config/settings.py` — no code change needed.

### Required: Land-use files (pct + number)

`filtered_land_use_pct.csv` is the one the standardiser ingests for the
weighted index; counts are kept for the rank-table tooltip. The loader
renames the pandas-deduplicated `undeveleoped_pct.1` column to
`nondeveloped_pct` per the challenge prompt wording.

### Optional: IoD-2025 main scores

Drop the MHCLG official summary file as `iod2025_main_scores.csv` to get
the headline domain scores alongside the standardised sub-domain rollups
in the rank table.

### Optional: Flood risk CSV

Two columns: LSOA ID + a numeric risk score. The pipeline auto-detects
the score column (first non-ID numeric col).

### Optional: LSOA boundary GeoJSON

Any FeatureCollection with `LSOA21CD` (or one of the candidate IDs)
in each feature's properties. The choropleth uses Plotly's carto-positron
tile-set; ensure your file is in EPSG:4326 (lon/lat).

---

## The math, briefly

For each LSOA *i* and active covariate *j*:

```
                  ┌  +1  if covariate j is an OPPORTUNITY signal
   sign  s_j  =  ┤            (e.g. vacant land, market failure)
                  └  -1  if it is a DEPRIVATION signal
                              (e.g. income, health, crime)

   raw_i  =  Σ_j  s_j · w_j · x_ij        where  x_ij ∈ [0, 1]

   opportunity_score_i  =  100 ·  (raw_i − min(raw)) / (max(raw) − min(raw))
```

Why the sign-flip trick? It lets one engine handle both kinds of signal
in a single composite: an analyst can simultaneously demand "high
opportunity" (vacant land, market failure) AND "low deprivation" (avoid
deeply income-poor areas) by leaving the deprivation sliders at moderate
values. Without the explicit `s_j` it would be impossible to mix the two
without writing per-variable boilerplate.

Polarity per covariate is declared once in `config/settings.WEIGHT_POLARITY`.

---

## Extending the platform

| Want to ... | Edit |
|---|---|
| Add a new sub-domain CSV | `config/settings.SUBDOMAIN_FILES` |
| Add a new derived metric (e.g. transport accessibility) | `src/data_processing/<new_module>.py`, then add a column registration in `WeightedIndexEngine.COVARIATE_COLUMNS` |
| Change a default weight | `config/settings.DEFAULT_WEIGHTS` |
| Change decile/quintile bin counts | `config/settings.N_DECILES` / `N_QUINTILES` |
| Switch from quantile to uniform binning | `config/settings.BINNING_STRATEGY` |
| Replace the UI framework | rebuild `src/ui/*`; the rest of the pipeline is UI-agnostic |

---

## Ethical considerations (per the challenge brief)

* **Transparency.** The rank-table tooltip shows the top-5 covariate
  contributions to every score. Analysts can audit *why* an LSOA was
  prioritised, not just *that* it was.
* **Aggregation risk.** Domain rollups average sub-domain deciles, which
  loses local nuance. The rank table preserves the underlying decile per
  domain so deep-dive is one click away.
* **Equity.** Deprivation covariates default to *negative* weights in the
  composite (high deprivation pushes the score *down*), reflecting the
  policy view that housing-led intervention alone is insufficient in
  deeply deprived places — but this is a defaulted assumption, not a
  hard-coded one. Reverse the polarity in `settings.WEIGHT_POLARITY` if
  your fund explicitly targets the most-deprived areas.
* **Data lineage.** Every transformation is traceable: raw → scaled →
  decile columns coexist in the wide frame, so any score is reproducible
  from the source CSV without losing information along the way.

---

## License

Built for the Cambridge Policy Hackathon 2026. Released under MIT for the
educational use of Homes England, partner institutions, and other public-
sector teams that find the architecture useful.
