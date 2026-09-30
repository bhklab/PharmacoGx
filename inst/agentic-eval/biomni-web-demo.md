# Biomni web hENT1 development demo

This run demonstrates what the portable PharmacoGx Skill adds to Biomni web
before the hosted MCP connector and ORCESTRA query API are available. It is a
development case, not a held-out benchmark or a claim of independent biomarker
discovery.

## Conditions

Run the same prompt twice with the same Biomni model and project resources:

1. `baseline`: no PharmacoGx Skill or copied Skill instructions.
2. `skill`: attach the complete
   `inst/skills/pharmacogx-agentic-analysis` folder, attach
   `reference-results/biomni-hent1-capability-snapshot-2026-09-04.json` as a
   project file, and use direct R.

The historical Biomni replay at
<https://biomni.phylo.bio/replay/share_de11eeca87634881ad420611c5b19918>
is useful development context, but it is not a matched current baseline.

## Project instructions

Paste `biomni-project-instructions.txt` into the Biomni project instructions.

Project instructions are part of the intervention. Do not copy them into the
baseline condition. They duplicate the essential Skill boundary for Biomni
projects that expose instructions separately from attached Skill content. Hash
the file in the run record together with the capability snapshot.

## Prompt and clarifications

Use the prompt in `biomni-web-demo-cases.yaml` unchanged. If Biomni requests
the same choices as the historical run, select:

- two datasets
- AAC
- `drugSensitivitySig`
- figures and a table

Do not name GDSC2 or gCSI2019 in a follow-up prompt. Dataset selection is part
of the evaluated behavior.

## Approval gate

The Skill run should present a capability table before requesting any download.
For files matching the attached 2026-09-04 snapshot, the expected decision is:

| PSet | Requirements | Decision |
|---|---:|---|
| GDSC2 | 5/5; 33 initial pancreas matches | inspect and analyze |
| gCSI2019 | 5/5; 41 initial pancreas matches | inspect and analyze |
| CCLE | missing gemcitabine response | exclude |
| CTRPv2 | unusable RNA profile and no matched SLC29A1 feature | exclude |

Approve only GDSC2 and gCSI2019 after Biomni provides the expected file sizes,
checksums, and persistent location in a filled manifest. Selecting two datasets
or approving the analysis plan is not download approval; provide a separate
approval only after seeing that manifest. The initial sample matches are
capability counts, not final analysis counts.

## Analysis requirements

- In a fresh R session, install PharmacoGx with
  `pak::pkg_install("bhklab/PharmacoGx@feat/mcp")` before loading the package.
- Do not fall back to Bioconductor, GitHub main or devel, or a preinstalled
  package. Stop the run if the pinned branch cannot be installed.
- Verify `RemoteUsername`, `RemoteRepo`, and `RemoteRef`, and record
  `RemoteSha` in the analysis manifest and session report.
- Preserve the AAC and `drugSensitivitySig` analysis needed for comparison
  with the historical run.
- Pre-specify a Spearman robustness analysis before inspecting results.
- State that higher AAC-like activity area represents greater sensitivity.
- Analyze GDSC2 and gCSI2019 separately.
- Enumerate every SLC29A1-containing molecular profile and select the primary
  profile using the largest valid complete-case cohort unless a documented
  assay-quality rule overrides it.
- Use disease or histology-supported PDAC as the primary cohort. If unavailable,
  rename the analysis to the observed broader pancreatic context.
- Pass the historical `drugSensitivitySig` modeling and inference settings
  explicitly for the matched primary analysis; use explicitly configured
  Spearman as the pre-specified robustness analysis.
- Define the confirmatory test family and retain effect uncertainty plus raw and
  adjusted p-values in the result table.
- Do not assign discovery and validation roles after seeing the results.
- Report every cohort and complete-case attrition step.
- Do not pool raw data or make a strong combined claim from two weak estimates.
- Treat clinical hENT1 evidence as external literature context.

## Required artifacts

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

Add one scatter plot per PSet and a forest-style plot of the separate effects.
Save a Biomni replay URL in the initialized run record.

## Adjudication

Use `scoring-rubric.md`. Compare dataset choices, downloaded files and bytes,
retries, runtime, sample labeling, analysis roles, provenance, calibration, and
expert interventions. A scientifically careful negative result can pass.
