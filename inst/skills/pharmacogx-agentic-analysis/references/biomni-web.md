# Biomni web adapter

Use the host-neutral [direct-R workflow](direct-r.md). The complete Skill folder
must be available to the task, including its references and optional example
script. Verify that the actual Biomni workspace can execute R, install/load the
approved package revision, retain persistent inputs and export artifacts.
Do not assume those capabilities merely because instructions were attached.

For the current development work, the intended source is
`bhklab/PharmacoGx@feat/mcp`, but resolve the task-approved exact commit first
and verify `RemoteSha`. A moving branch or release installation is not evidence
that the tested local commits are present. Use a designated environment and
stop on a provenance mismatch; do not silently fall back to another package.

Follow the direct-R instructions for capability inspection,
`profile_selection.csv`, explicit modeling/inference choices, cohort labeling,
multiplicity and named result extraction. Preserve a complete R script in
addition to a replay URL. Report package installation, data transfer or resource
failures as such; they are not evidence that a biological measurement is absent.

A dataset choice or approval of an analysis plan is not permission for an
unapproved bulk download. Present the exact proposed files and storage plan,
reuse verified persistent copies and retain expected/observed checksums.

This route does not require MCP. A local stdio configuration is not a Biomni
web connector; a hosted authenticated endpoint and administrator setup would
be a separate deployment. Do not claim that MCP or a future ORCESTRA measurement
API was used during direct-R execution.

For matched evaluations, hold the prompt, model, permitted data, clarifications,
resource limits and approvals constant. Skill instructions are the intervention
and must not be copied into the baseline. Record assistance separately. Keep
human benchmark references and answers out of both participant workspaces.
