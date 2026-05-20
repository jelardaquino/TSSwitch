# Contributor Covenant Code of Conduct and Contributing Guidelines — TSSwitch

Thank you for your interest in contributing to TSSwitch! This document describes how to report issues, suggest improvements, and submit code changes.

---

## How to Contribute

### Reporting Bugs or Issues

1. Search [existing issues](https://github.com/jelardaquino/TSSwitch/issues) first to avoid duplicates.
2. Open a new issue and include:
   - A clear, descriptive title
   - Steps to reproduce the problem
   - Expected vs. actual behavior
   - Your R version (`R.version`), OS, and relevant package versions (`sessionInfo()`)
   - Any error messages (full text, not just a screenshot)

### Suggesting Enhancements

Open an issue with the label `enhancement`. Describe:
- The analysis or feature you'd like to add
- Why it would benefit the pipeline
- Any relevant tools, packages, or methods you have in mind

### Submitting Code Changes

1. **Fork** the repository and create a new branch from `main`:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Follow the coding conventions** used in this project:
   - Descriptive variable names (`isoform_counts`, not `df2`)
   - Function header comments describing purpose, parameters, and return values
   - Section headers (e.g., `# ── Data Loading ──`) for readability
   - `dir.create(..., showWarnings = FALSE, recursive = TRUE)` for output directories
   - `stringsAsFactors = FALSE` in all `read.table()`/`read.delim()` calls
   - `set.seed()` for any stochastic steps

3. **Do not commit raw data** — never commit `.rds`, `.bam`, `.fastq`, `.fastq.gz`, `.gtf`, or large count matrix `.txt` files. See `.gitignore`.

4. **Do not hard-code absolute paths.** Use a configurable variable at the top of each script:
   ```r
   BASE_DIR <- "/path/to/your/data"
   ```

5. **Test your changes** on at least one of the two datasets (`PRJNA1206164` or `syn52047893`) before submitting.

6. **Open a Pull Request** against `main` with:
   - A summary of what was changed and why
   - Any new R package dependencies added
   - Example output (figure or table) if applicable

---

## Code of Conduct

This project follows the [Contributor Covenant](https://www.contributor-covenant.org/version/2/1/code_of_conduct/) Code of Conduct. Be respectful and constructive in all interactions.

---

## Questions?

Open a [GitHub Issue](https://github.com/jelardaquino/TSSwitch/issues) or start a [Discussion](https://github.com/jelardaquino/TSSwitch/discussions).
