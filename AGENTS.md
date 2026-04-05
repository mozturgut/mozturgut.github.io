# AGENTS.md

## Cursor Cloud specific instructions

### Repository overview

This repository contains two independent projects sharing a single Git repo:

1. **Portfolio** (`Portfolio/`) — A static HTML/CSS/JS personal portfolio website (Bootstrap 5.3.3, iPortfolio template). All vendor dependencies are pre-bundled; no package manager or build step is needed.
2. **Cell-Cell Communication** (`Cell-Cell Communication/`) — R-based bioinformatics pipeline (~5 scripts) for analyzing cell-cell communication in spinal cord injury datasets using CellChat, CellCall, LIANA, NicheNet, etc.

### Portfolio site

Serve with any static HTTP server. Quickest option:

```
cd Portfolio && python3 -m http.server 8080
```

Then open `http://localhost:8080/`. No build, lint, or test tooling exists for this project.

### Cell-Cell Communication pipeline

- **R ≥ 4.x** is required but is **not installed** in the Cloud Agent VM by default.
- Scripts require external data files (`LeeDat.rds`, `WangDat.rds`, `QinDat.rds`) that are **not included** in the repository.
- Scripts expect high memory (originally targeted at a 256 GB RAM workstation).
- Working directory is controlled by the `CELLCOMM_PROJECT_ROOT` environment variable; if unset, scripts use `getwd()`.
- Run order: `CellChat_ARG1_Final.R` → `CellChat_ARG1_FinalQin.R` → `CELLCALL_NEUTROPHIL_6_ANALYSES Linear.R` → `CELLCALL_QIN_ONLY_ANALYSES.R` → `MASTER_POST_COMMUNICATION_PIPELINE.r`.
- R packages come from CRAN, Bioconductor, and GitHub — there is no `renv.lock` or other lockfile.

### Notes

- No linter, test framework, or CI/CD pipeline is configured for either project.
- No `package.json`, `requirements.txt`, `Makefile`, `Dockerfile`, or `docker-compose.yml` exists.
- The update script is intentionally a no-op (`echo "No dependencies to install"`) since all Portfolio assets are vendored and the R pipeline has no automated dependency management.
