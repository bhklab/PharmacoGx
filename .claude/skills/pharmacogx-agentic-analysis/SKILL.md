---
name: pharmacogx-agentic-analysis
description: Use PharmacoGx for reproducible preclinical pharmacogenomic analyses, including PharmacoSet feasibility checks, drug response, and biomarker questions.
---

# PharmacoGx workflow

Read and follow the canonical Skill at
`inst/skills/pharmacogx-agentic-analysis/SKILL.md` in this repository.
Resolve its reference links relative to that canonical Skill's directory.

Use the project's registered `pharmacogx` MCP server when available. Confirm
execution with `pgx_session_info`; if tools are unavailable, report the missing
connection instead of implying MCP execution. Respect the user's chosen route
and permissions, including download and output-file approvals.
