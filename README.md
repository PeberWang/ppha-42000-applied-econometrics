# PPHA 42000 · Applied Econometrics I — Open Collaboration

Shared repository for **PPHA 42000 Applied Econometrics I (PhD)**, Autumn 2026
(MACRM, Harris School of Public Policy, University of Chicago). Instructor: Steven N. Durlauf.

A place to build up, over the quarter, a systematic body of work: rigorous proofs,
problem-set solutions, self-written notes, and data cleaning / analysis code.

> **Scope & copyright.** This repository contains only original, publicly shareable
> material. Instructor lecture notes, slides, the syllabus, and textbook PDFs are **not**
> included here — they live outside this repo and may be copyrighted. Please do not add
> them. Reading-list papers are likewise not redistributed here.

## Layout

| Path | What goes here |
|---|---|
| `proofs/` | Derivations and proofs (identification results, estimator properties, asymptotics) |
| `problem-sets/` | Solutions, organized as `ps01/`, `ps02/`, … |
| `notes/` | Self-written notes (LaTeX / Markdown) |
| `data/` | Datasets you are allowed to share, plus data dictionaries |
| `code/` | Cleaning and analysis code (Python / R) |

## How to contribute

1. **Branch** — `git checkout -b ps03/alice` (or `proof/ols-finite-sample`, `code/cleaning-pset02`).
2. **Commit** — small, focused commits. Suggested messages: `ps03: solve exercise 2`,
   `proof: add Gauss-Markov derivation`.
3. **Pull request** — open a PR into `main` and tag a maintainer. Keep each PR scoped to
   one problem set, one proof, or one code change.
4. **Math** — write math in LaTeX. Use separate `.tex` files for long derivations,
   Markdown for short write-ups.

### Conventions

- One result per file in `proofs/`; state assumptions before the result.
- Always state your sources (lecture-note number, Hayashi / Greene section) so others can verify.
- **Never** commit copyrighted course material — problem-set PDFs, slides, textbook
  scans, or reading-list papers. Paraphrase problems and write solutions in your own words.

## Maintainers

- Mingpei Wang ([@PeberWang](https://github.com/PeberWang)) — *add co-maintainers here as needed*

## License

- **Code** (`code/`, scripts, notebooks): MIT — see [`LICENSE`](LICENSE).
- **Prose & notes** (proofs, solutions, notes, this README): CC BY 4.0 — see
  [`LICENSE-CC-BY-4.0.md`](LICENSE-CC-BY-4.0.md).
