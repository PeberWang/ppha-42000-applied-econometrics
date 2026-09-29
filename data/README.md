# Data

Datasets used for exercises, plus their documentation.

## Rules

- **Do not commit** restricted, private, or unlicensed data — e.g. Canvas-internal
  datasets or micro data containing personal information. Ask before adding anything
  whose license you are unsure of.
- Keep large files out of Git: use Git LFS or an external link documented in the
  relevant README.
- Every dataset gets a short README covering: **source**, **license**, **how to obtain
  it**, and a **variable dictionary**.

## Suggested structure

```
data/
  raw/         # original files (only if redistribution is allowed)
  processed/   # cleaned / derived data
  README.md
```
