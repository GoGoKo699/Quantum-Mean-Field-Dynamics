# Hosted mathematical review

The read-only collector checks mathematics in actual commit-pinned GitHub pages. It creates a reference MathML representation from the exact Git source, collects the rendered expressions in an anonymous Chromium context, and compares their symbols and structure at viewport widths 1440 and 760.

Install the runtime in a temporary directory outside this repository:

```bash
npm install --prefix /absolute/temporary/dependencies playwright@1.56.1 mathjax-full@3.2.1
/absolute/temporary/dependencies/node_modules/.bin/playwright install chromium
NODE_PATH=/absolute/temporary/dependencies/node_modules node verification/hosted/scripts/capture_hosted_math.cjs \
  --root /absolute/repository/checkout \
  --repository GoGoKo699/Quantum-Mean-Field-Dynamics \
  --revision FULL_40_CHARACTER_COMMIT \
  --output /absolute/new-empty-evidence-directory
```

The default scope is every tracked Markdown page at that commit. Repeated `--page PATH` arguments deliberately restrict a diagnostic; a representative-page pass does not establish complete coverage. Use a new output directory for each attempt. The output records the source tree, source hashes, actual observation times, runtime, HTML, screenshots and comparison results.

The GitHub Actions job uses a read-only repository token for checkout, disables persistent credentials and opens an anonymous browser for the public page captures. It requests no repository secrets or write permissions. Its uploaded evidence is temporary; the job's success is distinct from a human comparison of the screenshots with the source.

A pass checks the declared source and expression scope. It does not prove the mathematics, assess priority, establish legibility without visual inspection or guarantee future renderer behavior. Numerical verification and the scientific discussion remain in [the verification overview](../../README.md).

For an interrupted transport-only capture, `--resume-from /absolute/previous-evidence-directory` can reuse matching, sufficiently fresh passing captures into a new empty output directory. Policy or access failures require resolution of the actual condition and cannot be retried through a different route.
