#!/usr/bin/env python3
"""Stitch the per-test analysis figures into one downloadable PDF report.

The regression analysis saves PNGs via `--save <rundir>/ci`, producing
`<rundir>/ci_*.png` under test_cases/<test>/<rundir>/. This walks that tree and
writes one PDF page per figure, titled `<test> -- <figure>`, so a CI run yields
a single PDF artifact (e.g. changa-regression.pdf). Uses matplotlib only (no
LaTeX / extra dependency); if no figures are found it still writes a one-page
note so the artifact always exists.

    python3 setupfiles/make_pdf_report.py --root test_cases --out report.pdf
"""
import argparse
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def test_name(path, root):
    """test_cases/<test>/<rundir>/ci_x.png -> '<test>' (first dir below root)."""
    rel = os.path.relpath(path, root)
    parts = rel.split(os.sep)
    return parts[0] if len(parts) > 1 else "report"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", default="test_cases",
                    help="directory tree to search for figures")
    ap.add_argument("--glob", default="ci_*.png",
                    help="figure filename pattern saved by the analysis")
    ap.add_argument("--out", default="regression-report.pdf",
                    help="output PDF path")
    args = ap.parse_args()

    pngs = sorted(glob.glob(os.path.join(args.root, "**", args.glob),
                            recursive=True))

    with PdfPages(args.out) as pdf:
        if not pngs:
            fig = plt.figure(figsize=(8.3, 11.7))  # A4 portrait
            fig.text(0.5, 0.5, f"No figures found under {args.root}/**/{args.glob}",
                     ha="center", va="center")
            pdf.savefig(fig)
            plt.close(fig)
            print(f"no figures matched; wrote placeholder -> {args.out}")
            return

        for png in pngs:
            img = mpimg.imread(png)
            h, w = img.shape[0], img.shape[1]
            dpi = 100.0
            fig = plt.figure(figsize=(w / dpi, h / dpi + 0.4), dpi=dpi)
            ax = fig.add_axes([0, 0, 1, h / (h + 0.4 * dpi)])
            ax.imshow(img)
            ax.axis("off")
            fig.suptitle(f"{test_name(png, args.root)} -- {os.path.basename(png)}",
                         fontsize=10, y=0.995)
            pdf.savefig(fig, dpi=dpi)
            plt.close(fig)
            print(f"added page: {png}")

    print(f"wrote {len(pngs)} figure(s) -> {args.out}")


if __name__ == "__main__":
    main()
