test_that("agentic evaluation supports the Skill-only condition", {
  skip_if_not_installed("yaml")

  eval_dir <- system.file("agentic-eval", package = "PharmacoGx")
  schema_path <- file.path(eval_dir, "run-record.schema.json")
  demo_path <- file.path(eval_dir, "biomni-web-demo-cases.yaml")
  snapshot_path <- file.path(
    eval_dir,
    "reference-results",
    "biomni-hent1-capability-snapshot-2026-09-04.json"
  )
  exporter_path <- file.path(
    eval_dir,
    "scripts",
    "export-requirements-snapshot.R"
  )
  instructions_path <- file.path(eval_dir, "biomni-project-instructions.txt")

  expect_true(file.exists(schema_path))
  expect_true(file.exists(demo_path))
  expect_true(file.exists(snapshot_path))
  expect_true(file.exists(exporter_path))
  expect_true(file.exists(instructions_path))

  instructions <- paste(
    readLines(instructions_path, warn = FALSE),
    collapse = "\n"
  )
  instructions_flat <- gsub("[[:space:]]+", " ", instructions)
  expect_match(
    instructions,
    "bhklab/PharmacoGx@feat/mcp",
    fixed = TRUE
  )
  expect_match(instructions, "RemoteSha", fixed = TRUE)
  expect_match(
    instructions_flat,
    "a dataset selection or plan approval is not download approval",
    ignore.case = TRUE
  )
  expect_match(instructions, "every molecular profile", fixed = TRUE)

  schema <- jsonlite::read_json(schema_path, simplifyVector = TRUE)
  conditions <- schema$properties$environment$properties$condition$enum
  expect_setequal(conditions, c("baseline", "skill", "mcp", "plugin"))

  demo <- yaml::read_yaml(demo_path)
  expect_identical(demo$case_set, "development")
  expect_length(demo$cases, 1)
  expect_true(demo$cases[[1]]$data_scope$allow_download)
  expect_identical(
    demo$cases[[1]]$reference$status,
    "needs_independent_review"
  )

  snapshot <- jsonlite::read_json(snapshot_path, simplifyVector = TRUE)
  expect_identical(snapshot$reference_status, "needs_independent_review")
  expect_setequal(
    snapshot$scope$compatible_datasets,
    c("gdsc2", "gcsi2019")
  )
  snapshot_text <- paste(
    readLines(snapshot_path, warn = FALSE),
    collapse = "\n"
  )
  expect_false(grepl("/Users/|/home/", snapshot_text))
})

test_that("portable Skill routes execution without development-case leakage", {
  skill_dir <- system.file(
    "skills",
    "pharmacogx-agentic-analysis",
    package = "PharmacoGx"
  )
  skill_path <- file.path(skill_dir, "SKILL.md")
  mcp_reference <- file.path(skill_dir, "references", "mcp-tools.md")
  biomni_reference <- file.path(skill_dir, "references", "biomni-web.md")

  expect_true(file.exists(skill_path))
  expect_true(file.exists(mcp_reference))
  expect_true(file.exists(biomni_reference))

  skill <- paste(readLines(skill_path, warn = FALSE), collapse = "\n")
  biomni <- paste(
    readLines(biomni_reference, warn = FALSE),
    collapse = "\n"
  )
  expect_match(skill, "Choose the execution route", fixed = TRUE)
  expect_match(skill, "references/mcp-tools.md", fixed = TRUE)
  expect_match(skill, "references/biomni-web.md", fixed = TRUE)
  expect_false(grepl("hENT1|SLC29A1|gemcitabine", skill, ignore.case = TRUE))
  expect_match(
    biomni,
    'pak::pkg_install("bhklab/PharmacoGx@feat/mcp")',
    fixed = TRUE
  )
  expect_match(biomni, "RemoteSha", fixed = TRUE)
  expect_match(
    biomni,
    "Earlier choices or approval of an analysis plan do not satisfy",
    fixed = TRUE
  )
  expect_match(biomni, "profile_selection.csv", fixed = TRUE)
  expect_match(biomni, 'modeling.method = "spearman"', fixed = TRUE)
  expect_match(biomni, "internal `.Data` slot", fixed = TRUE)
})
