# PharmacoGx

==========

[![R build status](https://github.com/bhklab/PharmacoGx/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/bhklab/PharmacoGx/actions)

**Bioc-Release**: ![Bioconductor RELEASE](http://bioconductor.org/shields/build/release/bioc/PharmacoGx.svg)
**Bioc-Devel**: ![Bioconductor DEVEL](http://bioconductor.org/shields/build/devel/bioc/PharmacoGx.svg)

R package to analyze large-scale pharmacogenomic datasets.

## Installation

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("PharmacoGx")
```

### Development version

```R
if (!requireNamespace("pak", quietly = TRUE))
    install.packages("pak")
pak::pkg_install("bhklab/PharmacoGx")
```

## Agentic Integration Demo

Optional MCP and agent-skill demo assets live in `inst/mcp/` and
`inst/skills/`. They expose a small PharmacoGx tool surface for bundled
preclinical pharmacogenomic examples, including sensitivity summaries,
responder ranking, dose-response metrics, and simple synergy references.

See `inst/mcp/README.md` for setup with `mcptools`, Biomni, ToolUniverse, and
agent `SKILL.md` workflows.
