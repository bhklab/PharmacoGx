# Reporting contract

Every executed analysis should produce enough information to reproduce the
result and audit its interpretation.

## Minimum report

- execution route: MCP, direct R, or planning only
- host and model identifiers when available
- R and PharmacoGx versions, package repository, Git reference, and commit SHA
- CoreGx version, installed package path or explicit source checkout, dirty
  status, dependency lock/session information and complete Skill resource hashes
- exact execution command, start/end times, exit status and host transcript or
  replay location; record operator interventions and route changes
- PSet names, versions, source records, filenames, expected and observed file
  sizes and checksums, persistent locations, and download approvals
- treatment, sample, endpoint, molecular-profile, and feature identifiers
- all candidate molecular profiles, their usable coverage, the pre-specified
  selection rule, and the reason for the selected profile
- metric definition and sensitivity direction
- cohort field, stored disease or histology values, label used in outputs,
  exclusions, initial matches, missingness, and final complete-case counts
- primary endpoint and method; secondary or sensitivity analyses; modeling,
  inference, summary, standardization, transformation, and covariate settings;
  multiplicity family and adjustment decision
- per-PSet effect size, direction, machine-readable uncertainty, raw p-value,
  and adjusted p-value when applicable
- warnings, errors, unavailable capabilities, and user approvals
- distinction among executed evidence, attached metadata, and literature

## Artifacts

Use machine-readable formats for tables and manifests. Prefer:

```text
capability_audit.csv
download_manifest.csv
profile_selection.csv
analysis_manifest.json
included_samples.csv
per_dataset_results.csv
analysis.R
README.md
session-info.txt
```

Add PNG plots when they help interpretation. Do not make a figure the only
record of a numerical result.

For multiple PSets, keep dataset-level results in separate rows and preserve
the identifiers used in each source. If a combined summary is scientifically
justified, retain the individual inputs and document the integration model.
