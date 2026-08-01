"""
DNA_analysis_pipeline.py

Missed outpatient appointments (DNAs) at Cambridge University Hospitals:
the geographic linkage, the distance calculations, and the clustering analysis.

This file consolidates the Python half of the project. It is organised as a set
of named steps run from the command line, since two of them belong in the middle
of the R pipeline and the rest belong at the end:

    python DNA_analysis_pipeline.py postcode-link     # step 3, after R Section 2
    python DNA_analysis_pipeline.py verify-postcodes  # optional linkage audit
    python DNA_analysis_pipeline.py distances         # step 4, before R Section 3
    python DNA_analysis_pipeline.py clustering        # after R Section 3
    python DNA_analysis_pipeline.py all               # every step in order

The R half lives in DNA_analysis_pipeline.R; see README.md for the full ordering
and for the note on the postcode lookup reference data, which is not included in
this repository.
"""

import argparse
import csv
import os
import re
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

# Threads are pinned before numpy or sklearn are imported anywhere. The KMeans
# implementation in older scikit-learn builds leaks memory on Windows when
# OMP_NUM_THREADS is unset, and pinning it also makes the clustering
# reproducible run to run.
for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

import numpy as np
import pandas as pd


###############################################################################
# CONFIGURATION
#
# All paths are relative to the project root, matching the R script. It is
# important to run this file from the project root, since every read and write
# assumes that location.
###############################################################################

OPEN_DIR  = Path("open_data")
CLEAN_DIR = Path("clean_data")
OUT_DIR   = Path("outputs")

# The same pipeline stage names used in the R script
F_AGG        = CLEAN_DIR / "aggregate_specialty_data_24.csv"
F_AGG_PCD    = CLEAN_DIR / "aggregate_specialty_data_24_pcd_link.csv"
F_AGG_NAPTAN = CLEAN_DIR / "aggregate_specialty_data_24_pcd_NAPTAN_link.csv"
F_AGG_CCA    = CLEAN_DIR / "aggregate_specialty_data_24_linked_CCA.csv"

# The ONS postcode directory, split into one CSV per postcode area. Each file
# carries pcd, lat, long, and lsoa21 columns.
PCD_LOOKUP_DIR = OPEN_DIR / "multi_csv"

F_NAPTAN   = CLEAN_DIR / "NAPTAN_24_clean.csv"
F_ADD_DIST = CLEAN_DIR / "add_dist.csv"

# Addenbrooke's Hospital, Cambridge
ADD_LAT, ADD_LON = 52.176, 0.140

EARTH_RADIUS_KM = 6371
RANDOM_STATE = 42

for _d in (CLEAN_DIR, OUT_DIR):
    _d.mkdir(exist_ok=True)


###############################################################################
# SHARED HELPERS
###############################################################################

def haversine(lat1, lon1, lat2, lon2):
    """Great circle distance in km between two points given in degrees."""
    lat1, lon1, lat2, lon2 = map(np.radians, (lat1, lon1, lat2, lon2))
    dlat, dlon = lat2 - lat1, lon2 - lon1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    return EARTH_RADIUS_KM * 2 * np.arcsin(np.sqrt(a))


_POSTCODE_CLEAN_RE = re.compile(r"\s+")
_AREA_RE = re.compile(r"^([A-Z]{1,2})(?=\d)")


def normalize_uk_postcode(pc):
    """
    Strip a UK postcode down to a comparable form: upper case, no internal
    spacing, with the various spellings of missing returned as None. Postcodes
    are recorded inconsistently across the specialty extracts, so normalising
    before any comparison is what makes the linkage rates believable.
    """
    # Missing postcodes arrive as None, as a float nan, or as pandas NA
    # depending on how the column was read, so all three are caught here rather
    # than relying on the caller to have cleaned the column first.
    if pc is None or pc is pd.NA or (isinstance(pc, float) and np.isnan(pc)):
        return None

    pc = str(pc).strip().upper()
    if pc in ("NA", "N/A", "NAN", "NONE", ""):
        return None
    pc = _POSTCODE_CLEAN_RE.sub("", pc)
    if not any(ch.isdigit() for ch in pc):
        return None
    return pc


def postcode_area_prefix(pc_norm):
    """Return the leading alphabetic area code, which selects the lookup file."""
    m = _AREA_RE.match(pc_norm)
    if m:
        return m.group(1)

    letters = []
    for ch in pc_norm:
        if ch.isalpha():
            letters.append(ch)
        else:
            break
    if not letters:
        raise ValueError(f"Cannot parse postcode area from: {pc_norm!r}")
    return "".join(letters)[:2]


###############################################################################
# STEP 3 — POSTCODE TO LSOA AND COORDINATE LINKAGE
#
# The referral data carry a postcode; everything geographic downstream needs a
# latitude, a longitude, and a 2021 LSOA code. The ONS directory is large enough
# that it is held as one file per postcode area, so the lookup is built as a
# dictionary keyed on the normalised postcode and the referral data mapped
# through it in one pass rather than searched file by file per row.
###############################################################################

def index_lookup_files(path=PCD_LOOKUP_DIR, file_key_regex=r"_([A-Z]{1,2})\.csv$"):
    """
    Map each postcode area to its file in the ONS directory. The directory is
    split by area, so CB postcodes live in a file whose name ends in _CB.csv,
    and only that file needs opening to resolve a Cambridge postcode.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"Postcode lookup directory not found: {path}. "
            "See README.md for how to obtain the ONS postcode directory."
        )

    if path.is_file():
        return {None: [path]}

    pattern = re.compile(file_key_regex, re.IGNORECASE)
    index = {}

    for fp in path.glob("*.csv"):
        m = pattern.search(fp.name)
        if m:
            index.setdefault(m.group(1).upper(), []).append(fp)

    # Fall back to matching on the leading letters of the file name, in case the
    # directory uses a naming convention the regex above does not anticipate.
    if not index:
        for fp in path.glob("*.csv"):
            letters = "".join(c for c in fp.stem if c.isalpha())[:2].upper()
            if letters:
                index.setdefault(letters, []).append(fp)

    return index


def resolve_area_file(prefix, index):
    """
    Pick the file for a postcode area. Two-letter areas (CB) are tried first,
    then the single-letter form (C), since the directory may be split either
    way. Where several files match, the most recently modified one wins.
    """
    candidates = index.get(prefix) or index.get(prefix[0]) or []
    if not candidates:
        return None
    return max(candidates, key=lambda fp: fp.stat().st_mtime)


def read_area_file(csv_path, wanted_postcodes, fields, postcode_col="pcd",
                   encoding="utf-8-sig"):
    """
    Scan one area file and return the rows for the postcodes asked for. The file
    is read once and only the requested postcodes are retained, so memory stays
    proportional to the cohort rather than to the ONS directory.
    """
    found = {}

    with csv_path.open("r", encoding=encoding, newline="") as f:
        reader = csv.DictReader(f)
        if postcode_col not in (reader.fieldnames or []):
            print(f"  {csv_path.name}: no '{postcode_col}' column, skipped")
            return found

        for row in reader:
            pcd_norm = normalize_uk_postcode(row.get(postcode_col, ""))
            if pcd_norm in wanted_postcodes and pcd_norm not in found:
                found[pcd_norm] = {k: row.get(k) for k in fields}
                if len(found) == len(wanted_postcodes):
                    break   # every postcode in this area is accounted for

    return found


def lookup_postcodes(postcodes, path=PCD_LOOKUP_DIR,
                     fields=("lat", "long", "lsoa21")):
    """
    Resolve a set of postcodes against the ONS directory.

    The postcodes are grouped by area first, so each area file is opened exactly
    once and only the files actually needed are opened at all. Looking a
    postcode up one at a time would re-read the same file for every patient in
    the same town, which for a cohort of this size is the difference between
    seconds and hours.
    """
    index = index_lookup_files(path)
    print(f"  {len(index)} postcode area files indexed")

    # Group the cohort's postcodes by their area prefix
    by_area = {}
    unparseable = 0
    for pc in postcodes:
        try:
            by_area.setdefault(postcode_area_prefix(pc), set()).add(pc)
        except ValueError:
            unparseable += 1

    if unparseable:
        print(f"  {unparseable} postcodes could not be parsed to an area")

    lookup = {}
    missing_areas = []

    for prefix in sorted(by_area):
        wanted = by_area[prefix]
        csv_path = resolve_area_file(prefix, index)

        if csv_path is None:
            missing_areas.append(prefix)
            continue

        found = read_area_file(csv_path, wanted, fields)
        lookup.update(found)
        print(f"  {prefix:<2s} {csv_path.name:<28s} "
              f"{len(found):>6,} / {len(wanted):>6,} matched")

    if missing_areas:
        print(f"  no lookup file for areas: {', '.join(missing_areas)}")

    return lookup


def step_postcode_link():
    """Attach lat, long, and lsoa21 to the aggregate referral data."""
    print("=" * 70)
    print("STEP 3 - POSTCODE LINKAGE")
    print("=" * 70)

    df = pd.read_csv(F_AGG)
    print(f"Loaded {len(df):,} appointment rows")

    # Normalise once and keep the result, since the same postcode recurs across
    # a patient's appointments and there is no reason to parse it repeatedly.
    # The original pcd column is left untouched, being the join key for the
    # distance file later on.
    pcd_norm = df["pcd"].astype("string").map(normalize_uk_postcode)

    unique_pcds = {p for p in pcd_norm.dropna().unique()}
    print(f"Unique postcodes to resolve: {len(unique_pcds):,}")

    lookup = lookup_postcodes(unique_pcds)

    df["lat"]    = pcd_norm.map(lambda p: (lookup.get(p) or {}).get("lat"))
    df["long"]   = pcd_norm.map(lambda p: (lookup.get(p) or {}).get("long"))
    df["lsoa21"] = pcd_norm.map(lambda p: (lookup.get(p) or {}).get("lsoa21"))

    df["lat"]  = pd.to_numeric(df["lat"], errors="coerce")
    df["long"] = pd.to_numeric(df["long"], errors="coerce")

    matched = int(df["lsoa21"].notna().sum())
    print(f"\nResolved {len(lookup):,} of {len(unique_pcds):,} unique postcodes")
    print(f"Matched:   {matched:,} rows ({matched / len(df):.1%})")
    print(f"Unmatched: {len(df) - matched:,} rows (retained with NA geography)")

    df.to_csv(F_AGG_PCD, index=False)
    print(f"File created at: {F_AGG_PCD}")


def step_verify_postcodes():
    """
    Optional audit of the linkage. Re-runs the lookup independently and flags
    rows whose stored geography disagrees with the reference file, or whose
    postcode is absent from it. Worth running once after any change to the
    lookup source; the output is a short exception list, not a full copy.
    """
    print("=" * 70)
    print("POSTCODE LINKAGE VERIFICATION")
    print("=" * 70)

    output_file = CLEAN_DIR / "postcode_verification_errors.csv"
    fields_to_check = ["lat", "long", "lsoa21"]

    df = pd.read_csv(F_AGG_PCD, dtype=str)

    pcd_norm = df["pcd"].astype("string").map(normalize_uk_postcode)
    lookup = lookup_postcodes({p for p in pcd_norm.dropna().unique()},
                              fields=tuple(fields_to_check))

    print(f"Checking {len(df):,} rows...")
    flags = []

    for _, row in df.iterrows():
        pc_norm = normalize_uk_postcode(row.get("pcd"))
        info = lookup.get(pc_norm) if pc_norm else None

        if info is None:
            flags.append({"ptid": row.get("ptid"), "flag": "MISSING"})
            continue

        for field in fields_to_check:
            recorded = "NA" if pd.isna(row.get(field)) else str(row.get(field))
            expected = "NA" if info.get(field) is None else str(info.get(field))
            if recorded != expected:
                flags.append({"ptid": row.get("ptid"), "flag": "MISMATCH"})
                break

    pd.DataFrame(flags, columns=["ptid", "flag"]).to_csv(output_file, index=False)

    n_missing  = sum(f["flag"] == "MISSING" for f in flags)
    n_mismatch = sum(f["flag"] == "MISMATCH" for f in flags)
    print(f"Missing postcodes:  {n_missing:,}")
    print(f"Mismatched values:  {n_mismatch:,}")
    print(f"Results written to: {output_file}")


###############################################################################
# STEP 4 — DISTANCE CALCULATIONS
#
# Two distances are derived for every patient: the straight line distance to
# Addenbrooke's, which is the appointment location, and the distance to the
# nearest public transport access node, which stands in for the practical
# difficulty of getting there without a car. Both are great circle rather than
# routed distances, which understates travel for rural patients; that limitation
# is worth stating explicitly when the estimates are reported.
###############################################################################

def step_distances():
    """Nearest NAPTAN node distance and distance to Addenbrooke's."""
    print("=" * 70)
    print("STEP 4 - DISTANCE CALCULATIONS")
    print("=" * 70)

    from sklearn.neighbors import BallTree

    df     = pd.read_csv(F_AGG_PCD)
    access = pd.read_csv(F_NAPTAN)

    df["sex"]       = df["sex"].astype("Int64")
    df["ethnicity"] = df["ethnicity"].astype("Int64")

    # Rows without coordinates cannot have a distance and are set aside rather
    # than dropped, so that the output keeps one row per input row.
    valid_access = access.dropna(subset=["Latitude", "Longitude"]).copy()
    valid_mask   = df[["lat", "long"]].notna().all(axis=1)

    df_valid   = df[valid_mask].copy()
    df_invalid = df[~valid_mask].copy()

    print(f"Valid rows:   {len(df_valid):,}")
    print(f"Invalid rows: {len(df_invalid):,}  (these will stay NA)")

    # ---- Nearest transport node ------------------------------------------
    # A BallTree on the haversine metric finds the nearest of roughly 400,000
    # nodes for each patient in log time; the pairwise alternative would be
    # tens of billions of comparisons.
    print("Building BallTree...")
    points = np.radians(df_valid[["lat", "long"]].values)
    nodes  = np.radians(valid_access[["Latitude", "Longitude"]].values)
    tree   = BallTree(nodes, metric="haversine")

    print("Calculating nearest-node distances...")
    distances, indices = tree.query(points, k=1)

    nearest_idx = indices[:, 0]
    df_valid["NAPTAN_distance_km"]        = distances[:, 0] * EARTH_RADIUS_KM
    df_valid["nearest_NAPTAN_lat"]        = valid_access.iloc[nearest_idx]["Latitude"].values
    df_valid["nearest_NAPTAN_lon"]        = valid_access.iloc[nearest_idx]["Longitude"].values
    df_valid["nearest_NAPTAN_commonname"] = valid_access.iloc[nearest_idx]["CommonName"].values

    for col in ("NAPTAN_distance_km", "nearest_NAPTAN_lat",
                "nearest_NAPTAN_lon", "nearest_NAPTAN_commonname"):
        df_invalid[col] = np.nan

    # Recombine in the original row order so that nothing downstream depends on
    # the order in which the valid and invalid rows were processed.
    df_final = pd.concat([df_valid, df_invalid]).sort_index()
    df_final.to_csv(F_AGG_NAPTAN, index=False)
    print(f"File created at: {F_AGG_NAPTAN}")

    # ---- Distance to Addenbrooke's ---------------------------------------
    # Written as its own postcode keyed file, since it is joined on pcd rather
    # than on LSOA in the R linkage step.
    print("Calculating distances to Addenbrooke's...")
    df_valid["add_dist"] = haversine(df_valid["lat"].values,
                                     df_valid["long"].values,
                                     ADD_LAT, ADD_LON)
    df_invalid["add_dist"] = np.nan

    add_out = pd.concat([df_valid[["pcd", "add_dist"]],
                         df_invalid[["pcd", "add_dist"]]]).sort_index()
    add_out.to_csv(F_ADD_DIST, index=False)
    print(f"File created at: {F_ADD_DIST}")


###############################################################################
# CLUSTERING
#
# The regression models in the R script assume that each covariate acts on its
# own. Clustering asks the complementary question: are there recognisable types
# of patient in this deprivation and access space, and does knowing which type a
# patient belongs to tell us anything the individual covariates do not?
#
# Five algorithms are run rather than one, because each makes different
# assumptions about cluster shape, and structure that survives all five is more
# credible than structure that only one method finds. Agreement between them is
# quantified by the Adjusted Rand Index, and the substantive question is settled
# by whether cluster membership improves out of sample AUC over the covariates
# alone. It is important to be cautious with the result either way: a good
# silhouette score describes geometry, not clinical meaning.
###############################################################################

CLUSTER_FEATURES = [
    "IMD_decile", "income_decile", "employment_decile", "education_decile",
    "health_decile", "crime_decile", "vehicle_pct", "childcare_overall",
    "NAPTAN_distance_km", "age",
]

K_RANGE = range(2, 11)
N_INIT = 20
HDBSCAN_MIN_CLUSTER_SIZES = [50, 100, 200, 500]
OPTICS_MIN_SAMPLES = 50
OPTICS_XI = 0.05


def _load_first_referral():
    """
    One row per patient, taking their earliest referral. This mirrors agg_data in
    the R script, so that the clustering and the regressions describe the same
    sample rather than differing by how repeat appointments were handled.
    """
    df_raw = pd.read_csv(F_AGG_CCA)
    df = (df_raw.sort_values("referral_date")
                .groupby("ptid", sort=False)
                .first()
                .reset_index())
    print(f"Patients:       {df.shape[0]:,}")
    print(f"DNA prevalence: {df['DNA'].mean():.3f}  ({int(df['DNA'].sum()):,} events)")
    return df


def _prepare_cluster_matrix(df):
    """Median impute and z-score the clustering features."""
    from sklearn.impute import SimpleImputer
    from sklearn.preprocessing import StandardScaler

    df_cluster = df[CLUSTER_FEATURES + ["DNA"]].copy()
    df_cluster[CLUSTER_FEATURES] = df_cluster[CLUSTER_FEATURES].apply(
        pd.to_numeric, errors="coerce")

    missing = df_cluster[CLUSTER_FEATURES].isnull().mean().mul(100).round(1)
    print("\nMissingness (%):")
    print(missing[missing > 0].to_string() if missing.any() else "  None")

    # Median rather than mean imputation, since the deciles are ordinal and a
    # mean would place patients at values the scale does not contain.
    X_imp = SimpleImputer(strategy="median").fit_transform(df_cluster[CLUSTER_FEATURES])

    # Scaling is essential here: distances in kilometres and deciles on a one to
    # ten scale would otherwise contribute to the distance metric in proportion
    # to their units rather than their importance.
    X = StandardScaler().fit_transform(X_imp)
    y = df_cluster["DNA"].fillna(0).astype(int).values

    print(f"\nFeature matrix: {X.shape[0]:,} x {X.shape[1]}")
    return X, X_imp, y


def _silhouette_safe(X, labels, sample_size=10_000):
    """Silhouette on non-noise points only; None if fewer than two clusters."""
    from sklearn.metrics import silhouette_score

    mask = labels >= 0
    if mask.sum() < 2 or len(set(labels[mask])) < 2:
        return None
    return silhouette_score(X[mask], labels[mask],
                            sample_size=min(sample_size, int(mask.sum())),
                            random_state=RANDOM_STATE)


def _cluster_profiles(labels, X_original, y, feature_names):
    """Mean of each feature per cluster, on the original scale, plus DNA rate."""
    d = pd.DataFrame(X_original, columns=feature_names)
    d["cluster"] = labels
    d["DNA"] = y
    profile = d.groupby("cluster")[feature_names + ["DNA"]].mean().round(3)
    sizes = d.groupby("cluster").size().rename("n_patients")
    return pd.concat([sizes, profile], axis=1)


def _extract_or_table(fitted_model, feature_names):
    """Tidy odds ratio table from a statsmodels Logit fit, matching the R output."""
    params, conf, pvalues = fitted_model.params, fitted_model.conf_int(), fitted_model.pvalues

    def sig(p):
        if p < 0.001:
            return "***"
        if p < 0.01:
            return "**"
        if p < 0.05:
            return "*"
        if p < 0.10:
            return "."
        return ""

    return pd.DataFrame({
        "Term":     feature_names,
        "OR":       np.exp(params).round(3),
        "CI_lower": np.exp(conf[0]).round(3),
        "CI_upper": np.exp(conf[1]).round(3),
        "p_value":  pvalues.round(4),
        "sig":      [sig(p) for p in pvalues],
    })


def _logistic_validation(labels, y, X_covariates, method_name):
    """
    Test whether cluster membership carries information the covariates do not.
    Three models are compared by cross validated AUC: covariates alone, clusters
    alone, and both together. The delta is the number that matters; a cluster
    solution that raises AUC by less than about 0.005 is telling us the
    covariates already captured that structure.
    """
    import statsmodels.api as sm
    from scipy.stats import chi2_contingency
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold, cross_val_score

    mask = labels >= 0
    if mask.sum() < 50 or len(set(labels[mask])) < 2:
        return None

    lab, yy = labels[mask], y[mask]
    Xc = X_covariates[mask]

    dummies = pd.get_dummies(lab, prefix="cluster", drop_first=True).astype(int)

    chi2, p_chi2, dof, _ = chi2_contingency(pd.crosstab(lab, yy))

    logit_A = sm.Logit(yy, sm.add_constant(dummies.values)).fit(disp=False)
    X_B_df  = pd.concat([pd.DataFrame(Xc, columns=CLUSTER_FEATURES).reset_index(drop=True),
                         dummies.reset_index(drop=True)], axis=1)
    logit_B = sm.Logit(yy, sm.add_constant(X_B_df.values)).fit(disp=False)

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_STATE)

    def cv_auc(X_arr):
        clf = LogisticRegression(max_iter=1000, solver="lbfgs", class_weight="balanced")
        return cross_val_score(clf, X_arr, yy, cv=cv, scoring="roc_auc").mean()

    auc_covars  = cv_auc(Xc)
    auc_cluster = cv_auc(dummies.values)
    auc_both    = cv_auc(X_B_df.values)

    return {
        "method": method_name,
        "chi2": chi2, "chi2_p": p_chi2, "chi2_df": dof,
        "auc_covariates": auc_covars,
        "auc_cluster": auc_cluster,
        "auc_both": auc_both,
        "delta_auc": auc_both - auc_covars,
        "or_A": _extract_or_table(logit_A, ["Intercept"] + list(dummies.columns)),
        "or_B": _extract_or_table(logit_B, ["Intercept"] + CLUSTER_FEATURES + list(dummies.columns)),
    }


def step_clustering():
    """Run all five clustering algorithms and compare them."""
    import matplotlib
    matplotlib.use("Agg")   # headless backend, safe on a server
    import matplotlib.pyplot as plt
    import seaborn as sns
    from itertools import combinations

    from sklearn.cluster import KMeans, AgglomerativeClustering, OPTICS
    from sklearn.mixture import GaussianMixture
    from sklearn.metrics import (adjusted_rand_score, davies_bouldin_score,
                                 calinski_harabasz_score)

    print("=" * 70)
    print("CLUSTERING")
    print("=" * 70)

    out_dir = OUT_DIR / "clustering"
    out_dir.mkdir(exist_ok=True)

    df = _load_first_referral()
    X, X_imp, y = _prepare_cluster_matrix(df)
    X_raw_df = pd.DataFrame(X_imp, columns=CLUSTER_FEATURES)

    labels_by_method = {}
    sweeps = {}

    def metrics(labels):
        mask = labels >= 0
        n_clusters = len(set(labels[mask]))
        return {
            "MSC": _silhouette_safe(X, labels),
            "DB":  davies_bouldin_score(X[mask], labels[mask]) if n_clusters > 1 else None,
            "CH":  calinski_harabasz_score(X[mask], labels[mask]) if n_clusters > 1 else None,
            "n_clusters": n_clusters,
            "noise_frac": float((labels < 0).mean()),
        }

    # ---- 1. K-Means ------------------------------------------------------
    # Partitions into k spherical, roughly equal variance clusters. k is chosen
    # by the mean silhouette coefficient rather than by the elbow, since the
    # elbow is a judgement call and the silhouette is not.
    print("\n" + "=" * 70)
    print("ALGORITHM 1: K-MEANS")
    print("=" * 70)

    km_rows = []
    for k in K_RANGE:
        km = KMeans(n_clusters=k, n_init=N_INIT, random_state=RANDOM_STATE)
        lab = km.fit_predict(X)
        m = metrics(lab)
        km_rows.append({"k": k, **m, "Inertia": km.inertia_, "_labels": lab})
        print(f"  k={k:2d}  MSC={m['MSC']:.4f}  Inertia={km.inertia_:,.0f}")

    km_df = pd.DataFrame([{a: b for a, b in r.items() if a != "_labels"} for r in km_rows])
    best_km_row = max(km_rows, key=lambda r: r["MSC"] or -1)
    best_k = best_km_row["k"]
    labels_by_method["kmeans"] = best_km_row["_labels"]
    sweeps["kmeans"] = km_df
    print(f"\n  -> Optimal k = {best_k}  (MSC = {best_km_row['MSC']:.4f})")

    # ---- 2. HDBSCAN ------------------------------------------------------
    # Density based, so it finds arbitrarily shaped clusters and, unlike the
    # partitional methods, is allowed to label a point as noise rather than
    # forcing every patient into a group. The noise fraction is reported
    # alongside the silhouette, since a high silhouette on a small non-noise
    # core is not the same as a good clustering of the cohort.
    print("\n" + "=" * 70)
    print("ALGORITHM 2: HDBSCAN")
    print("=" * 70)

    hdb_df = None
    try:
        try:
            from sklearn.cluster import HDBSCAN as _HDBSCAN

            def _fit_hdbscan(mcs):
                return _HDBSCAN(min_cluster_size=mcs).fit_predict(X)
        except ImportError:
            import hdbscan as _hdbscan_pkg

            def _fit_hdbscan(mcs):
                return _hdbscan_pkg.HDBSCAN(min_cluster_size=mcs).fit_predict(X)

        hdb_rows = []
        for mcs in HDBSCAN_MIN_CLUSTER_SIZES:
            lab = _fit_hdbscan(mcs)
            m = metrics(lab)
            hdb_rows.append({"min_cluster_size": mcs, **m, "_labels": lab})
            print(f"  min_cluster_size={mcs:4d}  clusters={m['n_clusters']:3d}  "
                  f"MSC={m['MSC'] if m['MSC'] is None else round(m['MSC'], 4)}  "
                  f"noise={m['noise_frac']:.2%}")

        hdb_df = pd.DataFrame([{a: b for a, b in r.items() if a != "_labels"} for r in hdb_rows])
        valid = [r for r in hdb_rows if r["MSC"] is not None]
        if valid:
            best_hdb = max(valid, key=lambda r: r["MSC"])
            labels_by_method["hdbscan"] = best_hdb["_labels"]
            sweeps["hdbscan"] = hdb_df
    except ImportError:
        print("  HDBSCAN unavailable; skipping (pip install hdbscan)")

    # ---- 3. Gaussian mixture ---------------------------------------------
    # Probabilistic and elliptical, so it can represent clusters of unequal
    # spread that K-Means would split or merge.
    print("\n" + "=" * 70)
    print("ALGORITHM 3: GAUSSIAN MIXTURE MODEL")
    print("=" * 70)

    gmm_rows = []
    for k in K_RANGE:
        gmm = GaussianMixture(n_components=k, random_state=RANDOM_STATE)
        lab = gmm.fit_predict(X)
        m = metrics(lab)
        gmm_rows.append({"k": k, **m, "BIC": gmm.bic(X), "_labels": lab})
        print(f"  k={k:2d}  MSC={m['MSC']:.4f}  BIC={gmm.bic(X):,.0f}")

    gmm_df = pd.DataFrame([{a: b for a, b in r.items() if a != "_labels"} for r in gmm_rows])
    best_gmm = max(gmm_rows, key=lambda r: r["MSC"] or -1)
    labels_by_method["gmm"] = best_gmm["_labels"]
    sweeps["gmm"] = gmm_df
    print(f"\n  -> Optimal k = {best_gmm['k']}  (MSC = {best_gmm['MSC']:.4f})")

    # ---- 4. Agglomerative -------------------------------------------------
    # Hierarchical with Ward linkage, making no distributional assumption at all,
    # which is a useful contrast to the two model based methods above.
    print("\n" + "=" * 70)
    print("ALGORITHM 4: AGGLOMERATIVE (WARD)")
    print("=" * 70)

    agg_rows = []
    for k in K_RANGE:
        lab = AgglomerativeClustering(n_clusters=k, linkage="ward").fit_predict(X)
        m = metrics(lab)
        agg_rows.append({"k": k, **m, "_labels": lab})
        print(f"  k={k:2d}  MSC={m['MSC']:.4f}")

    agg_df = pd.DataFrame([{a: b for a, b in r.items() if a != "_labels"} for r in agg_rows])
    best_agg = max(agg_rows, key=lambda r: r["MSC"] or -1)
    labels_by_method["agglomerative"] = best_agg["_labels"]
    sweeps["agglomerative"] = agg_df
    print(f"\n  -> Optimal k = {best_agg['k']}  (MSC = {best_agg['MSC']:.4f})")

    # ---- 5. OPTICS --------------------------------------------------------
    # Density based like HDBSCAN, but tolerant of clusters whose densities
    # differ, which matters where dense urban LSOAs sit alongside sparse rural
    # ones in the same feature space.
    print("\n" + "=" * 70)
    print("ALGORITHM 5: OPTICS")
    print("=" * 70)

    optics = OPTICS(min_samples=OPTICS_MIN_SAMPLES, xi=OPTICS_XI)
    optics_labels = optics.fit_predict(X)
    m_opt = metrics(optics_labels)
    print(f"  clusters={m_opt['n_clusters']}  noise={m_opt['noise_frac']:.2%}  MSC={m_opt['MSC']}")
    if m_opt["n_clusters"] > 1:
        labels_by_method["optics"] = optics_labels

    optics_reachability = pd.DataFrame({
        "ordering": optics.ordering_,
        "reachability": optics.reachability_[optics.ordering_],
    })

    # ---- Cross-method stability ------------------------------------------
    # The Adjusted Rand Index between every pair of solutions. High agreement
    # means the structure is a property of the data rather than of one
    # algorithm's assumptions; low agreement across the board is itself a
    # finding, and argues against reporting any single solution as definitive.
    print("\n" + "=" * 70)
    print("CROSS-METHOD STABILITY (Adjusted Rand Index)")
    print("=" * 70)

    names = list(labels_by_method)
    ari_matrix = pd.DataFrame(np.eye(len(names)), index=names, columns=names)
    for a, b in combinations(names, 2):
        la, lb = labels_by_method[a], labels_by_method[b]
        mask = (la >= 0) & (lb >= 0)
        ari = adjusted_rand_score(la[mask], lb[mask]) if mask.sum() > 1 else np.nan
        ari_matrix.loc[a, b] = ari_matrix.loc[b, a] = round(ari, 4)
    print(ari_matrix.to_string())

    # ---- Validation and profiles per method -------------------------------
    print("\n" + "=" * 70)
    print("MASTER COMPARISON TABLE")
    print("=" * 70)

    comparison_rows, lr_results, profiles = [], {}, {}

    for name, lab in labels_by_method.items():
        m = metrics(lab)
        prof = _cluster_profiles(lab, X_imp, y, CLUSTER_FEATURES)
        profiles[name] = prof
        prof.to_csv(out_dir / f"profiles_{name}.csv")

        lr = _logistic_validation(lab, y, X_imp, name)
        lr_results[name] = lr

        row = {"Method": name, "Clusters": m["n_clusters"],
               "MSC": None if m["MSC"] is None else round(m["MSC"], 4),
               "Noise": round(m["noise_frac"], 4)}
        if lr:
            row.update({"Chi2_p":    f"{lr['chi2_p']:.2e}",
                        "AUC_covar": round(lr["auc_covariates"], 4),
                        "AUC_clust": round(lr["auc_cluster"], 4),
                        "AUC_both":  round(lr["auc_both"], 4),
                        "Delta_AUC": round(lr["delta_auc"], 4)})
            lr["or_A"].to_csv(out_dir / f"or_clusteronly_{name}.csv", index=False)
            lr["or_B"].to_csv(out_dir / f"or_covariates_{name}.csv", index=False)
        comparison_rows.append(row)

    comparison = pd.DataFrame(comparison_rows)
    print(comparison.to_string(index=False))

    # ---- Figures ----------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
    for ax, (nm, sw) in zip(axes, [("K-Means", km_df), ("GMM", gmm_df),
                                   ("Agglomerative", agg_df)]):
        ax.plot(sw["k"], sw["MSC"], marker="o", linewidth=2, color="#00BDB6")
        ax.set_xlabel("Number of clusters (k)")
        ax.set_ylabel("Mean silhouette coefficient")
        ax.set_title(nm, fontweight="bold")
        ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_dir / "fig1_parameter_sweeps.png", dpi=150, bbox_inches="tight")
    plt.close()

    fig, ax = plt.subplots(figsize=(6, 5))
    sns.heatmap(ari_matrix.astype(float), annot=True, fmt=".2f", cmap="RdBu_r",
                center=0, vmin=-1, vmax=1, ax=ax)
    ax.set_title("Cross-method agreement (ARI)", fontweight="bold")
    plt.tight_layout()
    plt.savefig(out_dir / "fig2_ari_heatmap.png", dpi=150, bbox_inches="tight")
    plt.close()

    n_meth = len(profiles)
    fig, axes = plt.subplots(n_meth, 1, figsize=(12, max(3, 3 * n_meth)))
    if n_meth == 1:
        axes = [axes]
    for ax, (name, prof) in zip(axes, profiles.items()):
        z = prof[CLUSTER_FEATURES].apply(lambda c: (c - c.mean()) / c.std(), axis=0)
        sns.heatmap(z, annot=prof[CLUSTER_FEATURES].round(1), fmt=".1f",
                    cmap="RdBu_r", center=0, linewidths=0.5, ax=ax,
                    cbar_kws={"label": "Z-score"})
        ax.set_title(f"Cluster profiles: {name}", fontweight="bold")
    plt.tight_layout()
    plt.savefig(out_dir / "fig3_cluster_profiles.png", dpi=150, bbox_inches="tight")
    plt.close()

    if any(lr for lr in lr_results.values()):
        fig, ax = plt.subplots(figsize=(7, 4))
        d = [(n, lr["delta_auc"]) for n, lr in lr_results.items() if lr]
        ax.barh([x[0] for x in d], [x[1] for x in d], color="#00BDB6")
        ax.axvline(0, color="grey", linestyle="--")
        ax.set_xlabel("Delta AUC from adding cluster membership")
        ax.set_title("Does cluster membership add predictive value?", fontweight="bold")
        plt.tight_layout()
        plt.savefig(out_dir / "fig4_delta_auc.png", dpi=150, bbox_inches="tight")
        plt.close()

    # ---- Exports ----------------------------------------------------------
    comparison.to_csv(out_dir / "clustering_comparison.csv", index=False)
    ari_matrix.to_csv(out_dir / "ari_matrix.csv")
    optics_reachability.to_csv(out_dir / "optics_reachability.csv", index=False)

    for name, sw in sweeps.items():
        sw.to_csv(out_dir / f"{name}_sweep.csv", index=False)

    assignment_table = pd.DataFrame({"ptid": df["ptid"].values})
    for name, lab in labels_by_method.items():
        assignment_table[f"cluster_{name}"] = lab
    assignment_table.to_csv(out_dir / "patient_cluster_assignments.csv", index=False)

    print(f"\nOutputs written to {out_dir}/")


###############################################################################
# COMMAND LINE ENTRY POINT
###############################################################################

STEPS = {
    "postcode-link":    step_postcode_link,
    "verify-postcodes": step_verify_postcodes,
    "distances":        step_distances,
    "clustering":       step_clustering,
}

# "all" deliberately omits verify-postcodes, which is a slow optional audit
# rather than part of the pipeline proper.
ALL_STEPS = ["postcode-link", "distances", "clustering"]


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("step", choices=list(STEPS) + ["all"],
                        help="which pipeline step to run")
    args = parser.parse_args()

    steps = ALL_STEPS if args.step == "all" else [args.step]

    for name in steps:
        STEPS[name]()
        print()


if __name__ == "__main__":
    main()
