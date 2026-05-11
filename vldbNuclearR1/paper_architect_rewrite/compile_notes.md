# Compile Notes

Command run from repository root:

```sh
latexmk -pdf paper_architect_rewrite/main_rewrite.tex
```

Result:
Compilation succeeded after allowing LaTeX to write generated artifacts in the Dropbox-backed paper directory.
The output PDF is `main_rewrite.pdf` at the repository root.
The final output has 13 pages.

PDF integrity note:
Direct output into the Dropbox/Overleaf-backed directory was observed once as a truncated 18KB PDF even though the LaTeX log reported a 5.9MB output.
To avoid this Dropbox File Provider issue, the final verified PDF was rebuilt with:

```sh
latexmk -g -pdf -outdir=/private/tmp/vldbNuclearR1-build paper_architect_rewrite/main_rewrite.tex
```

The verified PDF has size 5,907,320 bytes, starts with `%PDF-1.7`, ends with `%%EOF`, and `pdfinfo` reports 13 pages.
It was copied to `main_rewrite.pdf` and `paper_architect_rewrite/main_rewrite_compiled.pdf`.

Hard-error scan:
No fatal LaTeX errors, undefined-control-sequence errors, unresolved-reference errors, or unresolved-citation errors remained in `main_rewrite.log`.

Path scan:
`main_rewrite.tex` inputs the rewritten sections under `paper_architect_rewrite/sections/`.
The figure paths referenced by the rewritten sections were found during the successful LaTeX compile.

Warnings to review:
The hard-warning scan found no fatal errors, no unresolved references, no unresolved citations, no hyperref PDF-string warnings, and no overfull boxes after the caption rewrite.
The log contains a non-fatal `balance` warning, acmart `\vspace` warnings, relsize warnings from tiny algorithm fonts, underfull box warnings from float placement, and BibTeX metadata warnings inherited from `references.bib`.
These do not block compilation, but they should be addressed or accepted before submission.

Caption audit:
All 15 figure captions now state what the figure shows, what the reader should notice, and which claim or RQ the figure supports.
For broad results whose raw source is not fully present locally, the captions state only the qualitative plotted claim and preserve the raw-source TODOs.

TODO scan:
The rewrite intentionally leaves TODO comments where experimental evidence is not fully verified against local source files.
The final grep found 12 evidence TODO comments in compiled rewrite TeX files, plus one example TODO in `analysis/sentence_level_rules.md`.
