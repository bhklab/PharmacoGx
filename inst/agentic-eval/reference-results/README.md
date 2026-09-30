# Reference Results

Reference answers must be generated independently of the host output using
reviewed PharmacoGx workflows. Record:

- PSet file or stable dataset ID and version.
- File checksum where practical.
- PharmacoGx, CoreGx, R, and dependency versions.
- Exact sample, treatment, response, and feature identifiers.
- Matching and filtering rules.
- Signed effect, uncertainty, adjustment method, and sample count.
- Expected artifact properties or checksum.
- Reviewer and review date.

Do not promote a case from `needs_independent_review` to `verified` by copying
an answer produced by Claude, Biomni, or another evaluated host.

`biomni-hent1-capability-snapshot-2026-09-04.json` is a dated input-capability
record for the development demo. It establishes which checked PSet files can
support the requested inputs; it is not a reference biomarker result. Reuse it
only when file sizes and checksums match, and keep its status as
`needs_independent_review` until a second reviewer verifies the probe and
source records.
