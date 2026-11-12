---
name: orcestra-discovery
description: Discover candidate pharmacogenomic measurements and downloadable PharmacoSets using PharmacoDB GraphQL and the ORCESTRA PSet catalog, with explicit source/version mapping and local PSet verification. Use before selecting or acquiring data for PharmacoGx analyses.
---

# ORCESTRA ecosystem discovery

Use direct R for discovery. MCP is optional. Keep two distinct sources:

| Source | Establishes | Does not establish |
| --- | --- | --- |
| PharmacoDB GraphQL | Candidate datasets, entities and reported measurement coverage | Downloadable PSet identity, release equivalence, DOI or file checksum |
| ORCESTRA PSet catalog | Listed PSet names, dataset metadata, DOI and download link | Joint usable treatment/sample/endpoint/molecular coverage in the local object |

The current measurement backend is `https://pharmacodb.ca/graphql`, guided by
Matthew's [pharmacodb-graphql](https://github.com/mattbocc/pharmacodb-graphql).
It is a temporary discovery adapter, not a true ORCESTRA measurement API.
Read [backend contract](references/backend-contract.md) for verified query
shapes, limitations and the replacement boundary.

For disease-specific discovery, consider relevant single-disease PSets even
when tissue labels are missing or use equivalent terminology. The companion
analysis skill's [dataset handling](../pharmacogx-agentic-analysis/references/dataset-handling.md)
lists release-specific model types and tissue contexts; verify these against
the actual source and object. A catalog label alone does not establish sample
eligibility. Honor local-only task scope and reuse the exact verified local
files; skip remote discovery and acquisition when they are outside scope.

## Workflow

1. Record the requested drug, sample/disease context, endpoint and molecular
   feature. Resolve identifiers without conflating database IDs, public UIDs,
   PubChem identifiers, symbols and Ensembl IDs.
2. Query the smallest useful metadata selection. The supplied
   [R helper](scripts/discover.R) supports dataset enumeration/counts, one
   compound/gene lookup, and the PSet catalog. It never downloads PSets or asks
   for precomputed biomarker answers. Read its status before using its data.
3. For further queries, verify the live schema and use explicit page/row limits
   where supported. Avoid `all: true`, unbounded nested entity lists and bulk
   association queries. Some resolvers ignore pagination; narrow the query or
   stop rather than claiming a bounded complete scan.
4. Treat HTTP failure, invalid JSON, GraphQL `errors`, partial data, empty
   matches and missing fields as different outcomes. HTTP 200 does not prove
   success. Keep error evidence, but do not interpret partial/error responses
   as biological absence. Stop after a diagnostic query; retry only a justified
   transient failure within the task budget.
5. If acquisition is needed, consult `https://orcestra.ca/api/pset/available`
   (or `PharmacoGx::availablePSets(canonical = FALSE)`). Match candidate datasets
   to exact PSet records through explicit identity/release evidence. Preserve
   ambiguous/missing mappings; names alone do not establish equivalence.
6. Present the exact candidate file, source record, release, expected size and
   checksum if available, destination, compatibility evidence and remaining
   unknowns before any download. Honor existing authorization; obtain approval
   for downloads outside it. Missing checksums stay unknown, never invented.
7. Reuse verified persistent files. Preserve original bytes before migration:
   current `downloadPSet()` updates and rewrites the downloaded file. Its exact
   PSet name is not necessarily the dataset name, and its default destination
   is temporary. Use a persistent explicit destination and retain provenance
   for both source and transformed artifacts.
8. Before analysis, inspect the actual local PSet using documented accessors,
   resolve IDs and measure complete-case coverage. Catalog/schema compatibility
   is insufficient. Keep compatible, incompatible, unchecked and load-failed
   states distinct.

Save endpoint/backend identity, retrieval time, exact query and variables,
response/status/errors, source record, guide revision and mapping decisions.
For scientific execution, also retain package commits and PSet checksums.
Report these three conclusions separately: candidate measurement coverage,
listed downloadable PSet, and verified local analysis compatibility.

Do not upload data, change services, or promise future endpoints. For benchmark
tasks, obey the allowed-source policy; remote association statistics may expose
answers and must not enter agent inputs when the task evaluates their derivation.
