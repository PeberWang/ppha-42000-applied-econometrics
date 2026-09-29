# Code

Cleaning and analysis code for the course (Python / R).

## Conventions

- Declare dependencies: include `requirements.txt` (Python) or `renv.lock` (R).
- Start each script with a header comment: purpose, inputs, outputs, author.
- Notebooks: name them `NN-topic.ipynb`; keep them reproducible (fixed seeds) and strip
  noisy outputs before committing.
- Prefer functions over copy-paste; keep a script runnable top-to-bottom.
- Paths should be relative to the repo root — no absolute local paths.
