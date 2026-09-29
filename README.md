# Farcast Biosciences - Sample Activity Labeler

An interactive Streamlit application for labeling sample-level biological activity across multimodal assays (Histopathology, Cytokines, Flow Cytometry, and NanoString Gene Expression Signatures).

## Overview

- **Assay Integration**: Combines 4 assay readouts (H&E / IHC, Cytokine profiling, Flow cytometry, NanoString GES).
- **Automated GES Scoring**: Computes mean log2 expression across gene signatures.
- **Log2 Fold-Change Analysis**: Computes Treatment vs. Control arm changes.
- **Responsive Feature Selection**: Identifies responding features based on user-defined fold-change cutoffs and sample thresholds.
- **Direction & Signature Activity Calling**: Classifies feature directionality (Increased / Decreased / No meaningful change) and scores biological activity.
- **Publication-Ready Export**: Generates categorized Excel matrices (`Y` = Active) and flat binary tables.

## Running Locally

1. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

2. **Run Streamlit**:
   ```bash
   streamlit run rxb_activity_app.py
   ```

3. Open your browser at `http://localhost:8501`.

## Deployment to Railway

1. Push this repository to GitHub.
2. Log into [Railway.app](https://railway.app/).
3. Click **New Project** → **Deploy from GitHub repo**.
4. Select `sample_activity_labeler`.
5. Railway will automatically detect Python and use `railway.toml` / `Procfile` to build and deploy.
6. In **Settings** → **Networking**, click **Generate Domain** to get your public URL.
