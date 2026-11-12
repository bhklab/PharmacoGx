# Discovery backend contract

Verified 2026-10-27 with bounded unauthenticated requests:

- `POST https://pharmacodb.ca/graphql`: `{ datasets { id name } }` returned 10
  datasets. This is a dated observation of the legacy collection, not a claim
  of complete current PSet coverage.
- `GET https://orcestra.ca/api/pset/available`: returned 37 records. Observed
  fields include `name`, `doi`, `downloadLink`, `dateCreated`, nested `dataset`
  (including `name` and `versionInfo`) and `availableDatatypes`. File size and
  checksum were not exposed as top-level fields; consult the source record if
  needed. Missing metadata does not authorize guessed provenance.
- Guide revision: `mattbocc/pharmacodb-graphql` commit
  `193e18c4a93d46e1906bb1195ec0d3bac4758cb2`. Consult the repository for extended
  queries, but validate them against the live server. The helper's fixed query
  selections deliberately exclude precomputed association results.

## R use

Dependencies: `curl` and `jsonlite`; install only in the intended analysis
environment when needed. Sourcing this file performs no network request:

```r
source("/path/to/orcestra-discovery/scripts/discover.R")
record <- pgx_discover("datasets")
stopifnot(identical(record$status, "ok"))
record$data
# Additional bounded choices:
# pgx_discover("dataset_counts", "<exact observed dataset name>")
# pgx_discover("compound", "<requested compound name>")
# pgx_discover("gene", "<requested gene symbol>")
# pgx_discover("psets")
# Save record with jsonlite::write_json(..., auto_unbox=TRUE, pretty=TRUE).
```

`datasets` has no pagination; it is currently a small metadata list. Dataset
counts and single-entity lookups avoid nested collection expansion. The helper
uses variables, fixed query documents, HTTPS, a 30-second timeout and a 2 MiB
response limit. It performs no automatic retries. GraphQL errors suppress the
usable `data` field even if partial results exist; the raw response is retained
for diagnosis. Operation-specific shape checks reject unexpected payloads or
missing required record fields. Empty lists or null entities return `status = "empty"`, with no usable
`data`; this means the query returned no match, not that biology is absent.
Missing expected root fields return `missing_field`; malformed records return
`invalid_shape`. Set a new endpoint only after verifying its schema/permissions.

`ok` means the response has the expected basic shape, not that every nullable
field is informative. Null counts are unknown, not zero; missing compound UIDs
or gene annotations do not establish cross-source identity. A name-only PSet
catalog record is a candidate, not evidence of a downloadable release: inspect
its actual download link, release metadata and provenance before acquisition.

## Source-specific cautions

- GraphQL dataset IDs are not ORCESTRA PSet IDs. A single source dataset can
  correspond to multiple releases/curations; retain an explicit join ledger.
- A gene lookup proves entity identity, not expression coverage in a dataset.
  Dataset marginal counts likewise do not prove joint drug/gene/cohort support.
- Database association rows, when explicitly allowed, are precomputed results
  with their own model and endpoint definitions. They are neither raw molecular
  coverage nor independent reproduction of an analysis.
- Nested annotation resolvers can fail when annotations are missing. Avoid
  unnecessary nested fields; inspect returned errors. Do not bulk-fetch all
  compounds merely to work around a missing annotation on one candidate.

## Replaceable boundary

Keep backend-specific query documents and decoding in `scripts/discover.R`.
The returned record carries `backend`, `endpoint`, `operation`, `query`,
`variables`, `retrieved_at`, `http_status`, `status`, `errors`, `data` and
`raw_response`; `data` retains the source-native payload. Consumers must not
assume a unified registry or infer local compatibility from it.
The envelope also records the adapter schema and guide revision; backend release
is explicitly unknown unless the service supplies it. PSet release information
is retained inside the original catalog record.

A future ORCESTRA measurement API requires a separately verified adapter and
field/identity mapping. Preserve the provenance envelope and the three evidence
states, and test error/partial/empty behavior before switching. Do not change
only the URL and assume GraphQL semantics or IDs remain equivalent.
