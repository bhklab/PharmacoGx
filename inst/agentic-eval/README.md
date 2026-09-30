# PharmacoGx Plugin Evaluation

This directory defines a host-neutral evaluation contract for the PharmacoGx
Skill and MCP plugin. Claude, Biomni, and future compatible hosts receive the
same prompts in matched comparisons and write results using the same run-record
schema.

The development suite is intentionally separate from any locked manuscript
benchmark. Development cases may reveal implementation gaps and may be used to
improve the plugin. A case cannot remain held out after its prompt, reference
answer, or output has been used for tuning.

## Files

- `development-cases.yaml`: initial debugging and scientific-QC questions.
- `biomni-web-demo-cases.yaml`: matched Skill-only Biomni development demo.
- `biomni-web-demo.md`: portal setup, execution gates, and artifact checklist.
- `biomni-project-instructions.txt`: exact project-level intervention text for
  the Biomni Skill condition.
- `benchmark-case.schema.json`: machine-readable case contract.
- `run-record.schema.json`: output contract for host adapters.
- `scoring-rubric.md`: scoring dimensions and automatic failures.
- `reference-results/`: instructions for independently generated references.
- `scripts/validate-cases.R`: structural and leakage checks for case manifests.
- `scripts/init-run.R`: initializes an identical run record for each host and
  condition.
- `scripts/export-requirements-snapshot.R`: exports a dated, reviewable local
  PSet capability check for Skill-only hosts that lack the MCP preflight.

## Conditions

The contract supports four conditions:

1. `baseline`: host without the PharmacoGx plugin.
2. `skill`: host with the shared Skill and a declared direct-execution route,
   but no PharmacoGx MCP tools.
3. `mcp`: host with MCP tools but without the PharmacoGx Skill.
4. `plugin`: host with both the shared Skill and MCP tools.

The immediate Biomni web evaluation uses `skill`. Use `plugin` only after the
hosted PharmacoGx MCP connector is active. Do not compare conditions that used
different prompts, PSet versions, model settings, or unrecorded assistance as
though they were matched.

Record the exact host, model, model version, plugin commit, PSet versions,
temperature, resource limits, raw response, tool trace, artifacts, and timing.
Do not manually rescue a failed run; record a separate assisted run if needed.

## Usage

```sh
Rscript inst/agentic-eval/scripts/validate-cases.R

Rscript inst/agentic-eval/scripts/init-run.R \
  --host biomni-web \
  --condition skill \
  --model MODEL_ID \
  --plugin-commit GIT_COMMIT \
  --output /path/to/biomni-skill-run.json
```

The initialized JSON record contains every prompt with a `pending` result.
Host adapters should fill those records without changing the source manifest.

For the hENT1 development demo, validate and initialize the dedicated manifest:

```sh
Rscript inst/agentic-eval/scripts/validate-cases.R \
  inst/agentic-eval/biomni-web-demo-cases.yaml

Rscript inst/agentic-eval/scripts/init-run.R \
  --host biomni-web \
  --condition skill \
  --model MODEL_ID \
  --plugin-commit GIT_COMMIT \
  --manifest inst/agentic-eval/biomni-web-demo-cases.yaml \
  --resource-paths \
  /path/to/project-instructions.txt,/path/to/capability-snapshot.json \
  --output /path/to/biomni-skill-run.json
```

This YAML is the canonical PharmacoGx scientific evaluation contract, not a
universal LLM test-suite format. Adapters may translate it to another runner,
but should preserve prompts, limits, expected observations, and hard failures.

Generate a fresh capability snapshot from locally available PSets before a
Skill-only portal run:

```sh
Rscript inst/agentic-eval/scripts/export-requirements-snapshot.R \
  --datasets gdsc2,gcsi2019,ctrpv2,ccle_pharmacoset \
  --drug gemcitabine \
  --sample-query pancreas \
  --metric aac_recomputed \
  --mdata-type rna \
  --features SLC29A1 \
  --local-pset-dirs /absolute/path/to/PSets \
  --output /path/to/capability-snapshot.json
```

The snapshot is an evaluation input, not a reference answer. Review it, remove
machine-specific paths, and verify source file checksums before attaching it to
a portal task.
