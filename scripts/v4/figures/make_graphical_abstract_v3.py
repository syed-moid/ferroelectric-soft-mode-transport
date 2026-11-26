#!/usr/bin/env python3
"""Graphical abstract: graphical_abstract_v3.png / .pdf and graphical_abstract_v3_AB.png / .pdf.

Successor of make_graphical_abstract_v2.py: panel 2 shows two segments per bar, "At least one soft partner" (one + two soft partners, #eb6834) and
"No soft partner" (gray), printed with one decimal and the percent sign, the soft-partner value in bold; no hatch and
no callout. New headings and subtitles for panels 2 and 3, "polar ≈ 0.5%", "neither soft family" in the key, the
model names A and B once, and a new bottom line. Every plotted and printed number is read from the processed outputs
that feed the manuscript tables (FIT_FILES and FIELDS below: file and field for each number); panel 3's callout and
the bottom line give the range of the plotted ratios below 80 cm⁻¹. --without-archived draws models A and B only and
reads nothing that needs the files of the 2015 study; the figure with the archived fit needs the outputs of those
steps (--held; not in this release).

Three questions at 300 K (16^3 transport mesh, primary family rule c_min 0.5, TRISYM = 1), for model A, model B and
the archived 2015 fit under the settings of model A:
  1. share of kappa_P carried by each family                          (main Tables 5 and 8; SI Table S37)
  2. carrier-weighted linewidth of the main carriers with at least one / no soft partner (main Table 5, Sec. 5.1,
     Table 8; SI Table S36)
  3. linewidth share per available mode, rotational / other partners  (main Tables 6, 7 and 8; SI Table S34)

Usage (repository root):
  python scripts/v4/figures/make_graphical_abstract_v3.py --out DIR [--without-archived | --held HELD]
      [--a DIR] [--b DIR] [--without-names] [--size W H] [--dpi DPI]
--a / --b: processed outputs of parts a and b (default data/processed/v4/a and data/processed/v4/b); --held: the
processed outputs of the steps that need the files of the 2015 study (t4/t4_P_summary.json, t5_P.json, t5b_P.json). The layout checks at the end stop the script with an error if any text lies outside
the figure or two texts overlap; text is never rescaled for other sizes.
"""
import argparse
import importlib.util
import json
import sys
from pathlib import Path

import matplotlib
import matplotlib.text

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import FancyArrowPatch, Polygon, Rectangle  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent

# ---- colours: families as in Fig. 4 ---------------------------------------------------------------------------------
_spec = importlib.util.spec_from_file_location("make_t4_figures", HERE / "make_t4_figures.py")
_t4 = importlib.util.module_from_spec(_spec)
with matplotlib.rc_context():                      # its module-level rcParams (grid on, font size …) must not leak
    _spec.loader.exec_module(_t4)                  # module has a __main__ guard: importing draws nothing
ROT, POL = _t4.FAMCOL["rot"], _t4.FAMCOL["pol"]
INK, INK2, GRID, SURF = _t4.INK, _t4.INK2, _t4.GRID, _t4.SURF
INK3 = "#8a8984"                                   # legend marker fill (neutral; models differ by shape only)
# "neither" bars: lighter than Fig. 4's #7a7974 because the values are printed on these bars (as in v1).
NEI = "#b9b8b2"
SOFT = "#eb6834"   # "At least one soft partner" (separated from the blue and aqua families also for CVD)

# ---- sources --------------------------------------------------------------------------------------------------------
FIT_FILES = {   # fit label -> (directory option, (T4 summary, T5 bands, T5b bands) relative to that directory)
    "model A": ("a", ("t4/t4_A_summary.json", "t5/t5_A.json", "t5b/t5b_A.json")),
    "model B": ("b", ("t4/t4_B_summary.json", "t5/t5_B.json", "t5b/t5b_B.json")),
    "archived 2015 fit": ("held", ("t4/t4_P_summary.json", "t5_P.json", "t5b_P.json")),
}
DEF = "definitions/primary_c0.5"
# quantity -> (which file of FIT_FILES, field path, factor to percent or 1)
FIELDS = {
    ("KP", "rot"): (0, f"{DEF}/kappa_P_share_by_family/rot", 100),
    ("KP", "pol"): (0, f"{DEF}/kappa_P_share_by_family/pol", 100),
    ("KP", "neither"): (0, f"{DEF}/kappa_P_share_by_family/neither", 100),
    ("KP", "mixed"): (0, f"{DEF}/kappa_P_share_by_family/mixed", 100),
    ("LW", "one soft partner"): (0, f"{DEF}/kappa_weighted_nsoft_frac/1", 100),
    ("LW", "two soft partners"): (0, f"{DEF}/kappa_weighted_nsoft_frac/2", 100),
    ("LW", "no soft partner"): (0, f"{DEF}/kappa_weighted_nsoft_frac/0", 100),
    ("RATIO", "< 80"): (2, "bands/<80/ratio_rot_over_nonrot/total", 1),
    ("RATIO", "80–95"): (2, "bands/80-95/ratio_rot_over_nonrot/total", 1),
    ("RATIO", "95–105"): (2, "bands/95-105/ratio_rot_over_nonrot/total", 1),
    ("RATIO", "105–145"): (1, "bands/105-145/ratio_rot_over_nonrot/total", 1),
    ("RATIO", "145–160"): (1, "bands/145-160/ratio_rot_over_nonrot/total", 1),
}
BANDS = ["< 80", "80–95", "95–105", "105–145", "145–160"]


def _field(d, path):
    """Field path with '/' separators; pair labels like 'neither/rot' are joined back when needed."""
    parts, cur, i = path.split("/"), d, 0
    while i < len(parts):
        for j in range(len(parts), i, -1):
            key = "/".join(parts[i:j])
            if isinstance(cur, dict) and key in cur:
                cur, i = cur[key], j
                break
        else:
            raise KeyError(path)
    return cur


def load(dirs, fits):
    """Return KP, LW (soft / none), RATIO and the source table [(quantity, fit, file, field, value)]; dirs: {"a", "b",
    "held"}."""
    path = {f: [Path(dirs[FIT_FILES[f][0]]) / p for p in FIT_FILES[f][1]] for f in fits}
    js = {f: [json.loads(p.read_text()) for p in path[f]] for f in fits}
    val, src = {}, []
    for (grp, name), (i, fpath, fac) in FIELDS.items():
        for f in fits:
            v = float(_field(js[f][i], fpath)) * fac
            val.setdefault((grp, name), []).append(v)
            src.append((f"{grp}:{name}", f, str(path[f][i]), fpath, v))
    for f in fits:                                 # the files belong to the fit they are read for
        want = {"model A": "A", "model B": "B", "archived 2015 fit": "P"}[f]
        assert all(j.get("model") in (want, f"model {want}") for j in js[f]), (f, [j.get("model") for j in js[f]])
    KP = {k: val[("KP", k)] for k in ("rot", "pol", "neither")}
    # "At least one soft partner" = one + two soft partners; the zero-, one- and two-soft-partner shares partition
    # each carrier's linewidth, so the two segments close to 100 %
    LW = {"At least one soft partner": [a + b for a, b in zip(val[("LW", "one soft partner")],
                                                               val[("LW", "two soft partners")])],
          "No soft partner": val[("LW", "no soft partner")]}
    RATIO = {f: [val[("RATIO", b)][n] for b in BANDS] for n, f in enumerate(fits)}
    for n in range(len(fits)):
        assert val[("KP", "mixed")][n] == 0, fits[n]
    return KP, LW, RATIO, src


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--without-archived", action="store_true", help="models A and B only")
    ap.add_argument("--without-names", action="store_true", help="omit the line naming models A and B")
    ap.add_argument("--size", nargs=2, type=float, metavar=("W", "H"), default=(7.2, 3.9), help="inches")
    ap.add_argument("--dpi", type=float, default=300)
    ap.add_argument("--a", default="data/processed/v4/a")
    ap.add_argument("--b", default="data/processed/v4/b")
    ap.add_argument("--held")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    if not a.without_archived and not a.held:
        ap.error("--held is required for the figure with the archived fit (or use --without-archived)")

    FITS = ["model A", "model B"] + ([] if a.without_archived else ["archived 2015 fit"])
    NF = len(FITS)
    KP, LW, RATIO, src = load({"a": a.a, "b": a.b, "held": a.held}, FITS)
    print("source of every number (quantity, fit, file, field, value):")
    for q, f, p, fld, v in src:
        print(f"  {q:34s} {f:18s} {p:48s} {fld:58s} {v:.10g}")

    # ---- closing checks: each bar closes to 100 % within rounding -------------------------------------------------
    for i in range(NF):
        assert abs(sum(v[i] for v in KP.values()) - 100) < 0.02, i
        assert abs(sum(v[i] for v in LW.values()) - 100) < 1e-9, i
    # the published panel 2 values of models A and B (one decimal); stop if the sources give anything else
    EXPECT = {"model A": ("23.9", "76.1"), "model B": ("26.3", "73.7")}
    got = {f: (f"{LW['At least one soft partner'][i]:.1f}", f"{LW['No soft partner'][i]:.1f}") for i, f in enumerate(FITS)}
    if any(got[f] != EXPECT[f] for f in EXPECT):
        sys.exit(f"panel 2 values differ from the published ones: {got}")
    # the fixed words still describe the plotted values
    assert all(abs(x - 0.5) < 0.1 for x in KP["pol"]), KP["pol"]       # "polar ≈ 0.5%"
    assert all(20 <= x <= 30 for x in LW["At least one soft partner"])  # bottom line: "About a quarter"
    lo80 = [RATIO[f][0] for f in FITS]
    RANGE = f"{min(lo80):.1f}–{max(lo80):.1f}"
    RANGE80 = f"{RANGE} ×"

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5, "axes.edgecolor": INK2,
                         "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
                         "axes.spines.right": False, "hatch.linewidth": 0.9, "pdf.fonttype": 42})
    W, H = a.size
    fig = plt.figure(figsize=(W, H), dpi=a.dpi, facecolor=SURF)

    def head(x, y, title, sub):
        fig.text(x, y, title, fontsize=9, fontweight="bold", color=INK, va="top")
        fig.text(x, y - 0.047, sub, fontsize=7, color=INK2, va="top")

    # ---- wording (fixed); "\n" marks where a line is wrapped to fit -----------------------------------------------
    P2_TITLE = "2  How much scattering of the main heat carriers\ninvolves a soft partner?"
    P2_SUBTITLE = ("Main carriers together carry half the particle heat current.\n"
                   "Linewidth shares are weighted by the heat each carrier carries.")
    P3_TITLE = "3  How do rotational and\nnon-rotational partners compare?"
    P3_SUBTITLE = "Linewidth share per available mode:\nrotational ÷ non-rotational in the same band."
    BOTTOM = ("About a quarter of the main carriers’ linewidth involves a soft partner. Below 80 cm⁻¹, the linewidth "
              "share per available\nrotational mode is {RANGE} times that for non-rotational modes.")
    # ---- layout (figure fractions) ------------------------------------------------------------------------------------
    BAR_H = 0.115                                  # height of the two bar panels
    P1_HEAD, P1_AX = 0.775, 0.520                  # panel 1 heading and axes bottom (panel 3 heading at P1_HEAD too)
    P2_HEAD, P2_SUB, P2_AX = 0.475, 0.385, 0.175   # panel 2 heading, subtitle and axes bottom (= panel 3 bottom)
    P3_SUB, P3_TOP = 0.685, 0.585                  # panel 3 subtitle and axes top
    NAMES_X = 0.785                                # model names, right of the key

    # ---- headline and family key -------------------------------------------------------------------------------------
    fig.text(0.015, 0.975, "Soft phonons in cubic SrTiO₃ at 300 K: carrying heat and scattering its carriers",
             fontsize=10.5, fontweight="bold", color=INK, va="top")

    def icon_rot(ax):
        """Four corner-sharing octahedra seen along the rotation axis, neighbours turned in opposite senses."""
        ax.set_xlim(-1.15, 1.15), ax.set_ylim(-1.15, 1.15), ax.set_aspect("equal"), ax.axis("off")
        base = np.array([[0.5, 0], [0, 0.5], [-0.5, 0], [0, -0.5]])
        for cx, cy, sgn in ((-0.5, 0.5, 1), (0.5, 0.5, -1), (-0.5, -0.5, -1), (0.5, -0.5, 1)):
            ang = np.deg2rad(13 * sgn)
            rot = np.array([[np.cos(ang), -np.sin(ang)], [np.sin(ang), np.cos(ang)]])
            ax.add_patch(Polygon(base @ rot.T + [cx, cy], closed=True, fc=ROT, ec=SURF, lw=0.8))

    def icon_pol(ax):
        """One octahedron with its central cation displaced against the oxygen cage."""
        ax.set_xlim(-1.15, 1.15), ax.set_ylim(-1.15, 1.15), ax.set_aspect("equal"), ax.axis("off")
        ax.add_patch(Polygon(np.array([[0.95, 0], [0, 0.95], [-0.95, 0], [0, -0.95]]), closed=True, fc="none",
                             ec=POL, lw=1.3))
        ax.plot([0], [-0.30], "o", ms=3.4, color=POL)
        ax.add_patch(FancyArrowPatch((0, -0.30), (0, 0.62), arrowstyle="-|>", mutation_scale=6.5, color=POL, lw=1.2,
                                     shrinkA=0, shrinkB=0))

    def icon_box(ax, fc, hatch=None):
        ax.set_xlim(0, 1), ax.set_ylim(0, 1), ax.set_aspect("equal"), ax.axis("off")
        ax.add_patch(Rectangle((0.12, 0.12), 0.76, 0.76, fc=fc, ec=SURF, lw=0, hatch=hatch))

    KEY_Y, KEY_H = 0.825, 0.078
    for x, draw, t1, t2 in ((0.015, icon_rot, "rotational", "octahedra rotate (R–M edges)"),
                            (0.262, icon_pol, "polar", "cation–oxygen displacement (near Γ)"),
                            (0.560, lambda ax: icon_box(ax, NEI), "neither soft family", None)):
        draw(fig.add_axes([x, KEY_Y, KEY_H * H / W, KEY_H]))
        if t2 is None:
            fig.text(x + 0.050, KEY_Y + 0.037, t1, fontsize=8, fontweight="bold", color=INK, va="center")
        else:
            fig.text(x + 0.050, KEY_Y + 0.056, t1, fontsize=8, fontweight="bold", color=INK, va="center")
            fig.text(x + 0.050, KEY_Y + 0.018, t2, fontsize=6.6, color=INK2, va="center")
    if not a.without_names:   # the model names, once
        fig.text(NAMES_X, KEY_Y + 0.056, "A: ALAMODE tutorial fit", fontsize=6.6, color=INK2, va="center")
        fig.text(NAMES_X, KEY_Y + 0.018, "B: fit to independent force data", fontsize=6.6, color=INK2, va="center")
    fig.add_artist(Line2D([0.015, 0.985], [0.805, 0.805], color=GRID, lw=0.8))

    # ---- left: two stacked-bar panels --------------------------------------------------------------------------------
    LX, LW_AX = 0.150, 0.415
    Y = list(range(NF - 1, -1, -1))                # top bar = model A
    TOP = NF - 1

    def bars(ax, parts, colors):
        left = np.zeros(NF)
        for (name, vals), col in zip(parts.items(), colors):
            thin = max(vals) < 1                       # a sliver keeps its full width: no surface gap drawn over it
            ax.barh(Y, vals, left=left, height=0.56, color=col, edgecolor=SURF, linewidth=0 if thin else 1.1,
                    zorder=3 if thin else 2)
            left += np.array(vals)
        ax.set_xlim(0, 100), ax.set_ylim(-0.48, TOP + 0.48)
        ax.set_yticks(Y), ax.set_yticklabels(FITS, fontsize=7, color=INK)
        ax.tick_params(axis="y", length=0, pad=4)
        ax.set_xticks([0, 25, 50, 75, 100])
        ax.tick_params(axis="x", labelsize=6.5, length=2)
        ax.spines["left"].set_visible(False)

    head(0.015, P1_HEAD, "1  Which modes carry the heat?", "share of the particle heat current (%)")
    ax1 = fig.add_axes([LX, P1_AX, LW_AX, BAR_H])
    bars(ax1, {"rot": KP["rot"], "pol": KP["pol"], "neither": KP["neither"]}, [ROT, POL, NEI])
    ax1.set_xticklabels([])
    for y, r, n in zip(Y, KP["rot"], KP["neither"]):
        ax1.text(r / 2, y, f"{r:.1f}", ha="center", va="center", fontsize=6.8, color=INK)
        ax1.text(r + KP["pol"][0] + n / 2, y, f"{n:.1f}", ha="center", va="center", fontsize=6.8, color=INK)
    ax1.annotate("polar ≈ 0.5%", xy=(KP["rot"][0] + 0.25, TOP + 0.30), xytext=(21.5, TOP + 0.95), fontsize=6.6,
                 color=INK2, va="center", annotation_clip=False,
                 arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6, shrinkA=1, shrinkB=0))
    ax1.text(6.8, TOP + 0.95, "rotational", fontsize=6.6, color=INK2, ha="center", va="center", clip_on=False)
    ax1.text(57, TOP + 0.95, "neither family", fontsize=6.6, color=INK2, ha="center", va="center", clip_on=False)

    fig.text(0.015, P2_HEAD, P2_TITLE, fontsize=9, fontweight="bold", color=INK, va="top", linespacing=1.05)
    fig.text(0.015, P2_SUB, P2_SUBTITLE, fontsize=7, color=INK2, va="top", linespacing=1.1)
    ax2 = fig.add_axes([LX, P2_AX, LW_AX, BAR_H])
    soft, none_ = LW["At least one soft partner"], LW["No soft partner"]
    bars(ax2, LW, [SOFT, NEI])
    for y, sv, nv in zip(Y, soft, none_):
        ax2.text(sv / 2, y, f"{sv:.1f}%", ha="center", va="center", fontsize=6.8, color=INK, fontweight="bold")
        ax2.text(100 - nv / 2, y, f"{nv:.1f}%", ha="center", va="center", fontsize=6.8, color=INK)
    ax2.text(0, TOP + 0.95, "At least one soft partner", fontsize=6.6, color=INK2, ha="left", va="center",
             clip_on=False)
    ax2.text(100, TOP + 0.95, "No soft partner", fontsize=6.6, color=INK2, ha="right", va="center", clip_on=False)
    ax2.set_xlabel("percent", fontsize=6.8, color=INK2, labelpad=1)

    # ---- right: ratio per available mode -------------------------------------------------------------------------------
    RX = 0.655
    fig.text(RX - 0.045, P1_HEAD, P3_TITLE, fontsize=9, fontweight="bold", color=INK, va="top", linespacing=1.05)
    fig.text(RX - 0.045, P3_SUB, P3_SUBTITLE, fontsize=7, color=INK2, va="top", linespacing=1.1)
    ax3 = fig.add_axes([RX, P2_AX, 0.325, P3_TOP - P2_AX])
    ax3.axhline(1, color=INK2, lw=0.8)
    ax3.text(-0.45, 1.10, "equal per mode", fontsize=6.4, color=INK2, ha="left", va="bottom")
    MARK = {"model A": ("o", 5.6), "model B": ("s", 5.2), "archived 2015 fit": ("D", 4.8)}
    OFF = {3: (-0.20, 0.0, 0.20), 2: (-0.10, 0.10)}[NF]
    for (fit, vals), off in zip(RATIO.items(), OFF):
        mk, ms = MARK[fit]
        ax3.plot(np.arange(5) + off, vals, ls="none", marker=mk, ms=ms, mfc=ROT, mec=SURF, mew=0.9, zorder=3)
    ax3.set_xlim(-0.55, 4.55), ax3.set_ylim(0, 7.2)
    ax3.set_xticks(range(5)), ax3.set_xticklabels(BANDS, fontsize=6.6)
    ax3.set_yticks([0, 1, 2, 4, 6])
    ax3.tick_params(axis="both", labelsize=6.6, length=2)
    ax3.yaxis.grid(True, color=GRID, lw=0.6), ax3.set_axisbelow(True)
    ax3.set_xlabel("partner frequency band (cm⁻¹)", fontsize=6.8, color=INK2, labelpad=2)
    ax3.set_ylabel("ratio", fontsize=6.8, color=INK2, labelpad=2)
    ax3.text(0.48, 6.0, RANGE80, fontsize=9, fontweight="bold", color=INK, va="center")
    ax3.text(0.48, 5.25, "below 80 cm⁻¹", fontsize=6.6, color=INK2, va="center")
    ax3.legend([Line2D([], [], ls="none", marker=MARK[f][0], ms=4.6, mfc=INK3, mec=SURF, mew=0.6) for f in FITS],
               FITS, loc="upper right", bbox_to_anchor=(1.02, 0.80), fontsize=6.4, frameon=False, handletextpad=0.2,
               labelspacing=0.25, borderaxespad=0)

    fig.text(0.015, 0.018, BOTTOM.format(RANGE=RANGE), fontsize=7.2, color=INK, va="bottom", linespacing=1.15)

    # ---- layout checks (as v1): every text inside the figure, no two texts overlapping -------------------------------
    fig.canvas.draw()
    PW, PH = fig.canvas.get_width_height()
    texts = [t for t in fig.findobj(matplotlib.text.Text) if t.get_text().strip() and t.get_visible()]
    bad = []
    for t in texts:                                # whole artist, arrows of annotations included (as v1)
        b = t.get_window_extent()
        if not (b.x0 >= 0 and b.x1 <= PW and b.y0 >= 0 and b.y1 <= PH):
            bad.append(f"outside the figure: {t.get_text()!r}")
    # overlap of the text boxes themselves (as v1: Text.get_window_extent, without an annotation's arrow)
    boxes = [(t.get_text(), matplotlib.text.Text.get_window_extent(t)) for t in texts]
    for i, (ta, p) in enumerate(boxes):
        for tb, q in boxes[i + 1:]:
            if p.x0 < q.x1 - 2 and q.x0 < p.x1 - 2 and p.y0 < q.y1 - 2 and q.y0 < p.y1 - 2:
                bad.append(f"overlap: {ta!r} / {tb!r}")
    if bad:
        sys.exit(f"layout checks failed at {W} x {H} in, {a.dpi} dpi ({len(bad)}):\n  " + "\n  ".join(bad))
    print(f"layout checks passed: {len(boxes)} texts, {W} x {H} in, {a.dpi} dpi, {PW} x {PH} px")

    stem = "graphical_abstract_v3" + ("_AB" if a.without_archived else "")
    a.out.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(a.out / f"{stem}.{ext}", facecolor=SURF)
    print("written:", a.out / f"{stem}.png", a.out / f"{stem}.pdf")


if __name__ == "__main__":
    main()
