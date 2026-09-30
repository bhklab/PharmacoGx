# Scoring Rubric

Score each dimension from 0 to 2:

- `0`: failed, unsupported, or materially misleading.
- `1`: partially correct but incomplete or insufficiently supported.
- `2`: correct, complete for the case, and supported by the execution trace.

## Dimensions

1. **Dataset selection** — selects compatible PSets and distinguishes missing,
   unloadable, unchecked, and incompatible data.
2. **Execution** — completes the requested workflow through the declared route.
   MCP and Skill+MCP runs must use the typed tools; Skill-only direct-R runs
   must preserve their complete R script and must not claim MCP execution.
3. **Result correctness** — agrees with the independently generated reference
   within its stated tolerance.
4. **Scientific validity** — uses appropriate endpoints, signed directions,
   matching rules, covariates, and multiple-testing treatment.
5. **Provenance** — reports dataset and package versions, identifiers, sample
   counts, response fields, and generated artifacts.
6. **Calibration** — separates supported results, exploratory evidence,
   missing capabilities, and external biological context.
7. **Operational efficiency** — avoids repeated discovery, unnecessary or
   repeated downloads, repeated cold-loading of large PSets, and unrecorded
   recovery attempts; remains within the case limits.

Critical dimensions are declared per case. A case passes only when every
critical dimension scores 2 and no automatic failure occurs.

## Automatic Failures

- Fabricates dataset coverage, measurements, identifiers, or results.
- Treats a load or worker failure as biological absence.
- Pools raw observations across PSets without an approved integration model.
- Reverses the direction of a response or synergy metric.
- Ignores required multiple-testing adjustment or sample support.
- Claims validation from weak, directionally conflicting, or single-PSet
  evidence.
- Downloads data without explicit approval.
- Downloads a PSet again when the same verified file is available in persistent
  storage.
- Silently substitutes toy data for a real-data question.
- Claims MCP or ORCESTRA execution in a Skill-only direct-R run.
- Uses general knowledge as though it were executed PSet evidence.
- Uses manuscript benchmark cases to tune the plugin after the benchmark is
  locked.
