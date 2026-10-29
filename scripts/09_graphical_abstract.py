"""Graphical abstract (main figure, option 4): study schematic.

Draws the four-step schematic used as the graphical abstract:
(1) data sources, (2) diagnostic groups, (3) case-case odds ratio estimator,
(4) age-adjusted Mantel-Haenszel ORs for Chile and Brazil.

Inputs
    - results/tables/A20.csv: observed counts by country, age group and sex
      (adults >=20 y, 2015-2025), written by scripts/07_tables.py.
    - data/external/ne_110m_chile_brazil.geojson: country outlines extracted from
      Natural Earth 1:110m Admin 0 Countries, natural-earth-vector v5.1.2
      (public domain). Optional: the figure is drawn without them if missing.

Outputs
    - results/figures/graphical_abstract.png (300 dpi)
    - results/figures/graphical_abstract.pdf

Every number shown in the figure is computed here from Table A20 or checked
against it, and the statements in the figure are checked against the data;
the script stops with an error if a check fails.
No network access is required.
"""

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
A20_FILE = ROOT / "results" / "tables" / "A20.csv"
SHAPES_FILE = ROOT / "data" / "external" / "ne_110m_chile_brazil.geojson"
OUT_DIR = ROOT / "results" / "figures"
OUT_NAME = "graphical_abstract"

plt.rcParams.update({"font.size": 10, "font.family": "DejaVu Sans",
                     "axes.spines.top": False, "axes.spines.right": False})

# ---------------------------------------------------------------------------
# Data: Table A20, observed counts
# ---------------------------------------------------------------------------
AGES = ["20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "≥80"]
COUNTRIES = ["Chile", "Brazil"]
SEXES = {"Women": "W", "Men": "M"}
ARTERIAL = ["Coronary", "Cerebrovascular", "Other arterial"]
# A20.csv column -> group used here (Varices_I83 is not part of the main contrast)
A20_COLUMNS = {"PE": "PE", "I80_I82": "I80_82", "Coronary": "Coronary",
               "Cerebrovascular": "Cerebrovascular",
               "Other_arterial": "Other arterial"}


def load_a20(path):
    """Read Table A20 into {(country, sex): {group: counts ordered by AGES}}."""
    if not path.exists():
        raise SystemExit(f"{path} not found; run scripts/07_tables.py first.")
    with path.open(encoding="utf-8", newline="") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        rows = list(csv.DictReader(fh, dialect=csv.Sniffer().sniff(sample, ",;\t")))
    table = {}
    for r in rows:
        age = r["age"].strip().replace("\u2013", "-")  # en dash in the CSV
        key = (r["country"].strip(), SEXES[r["sex"].strip()])
        for col, group in A20_COLUMNS.items():
            table.setdefault(key, {}).setdefault(group, {})[age] = int(r[col])
    raw = {}
    for key, groups in table.items():
        raw[key] = {}
        for group, by_age in groups.items():
            missing = set(AGES) - set(by_age)
            assert not missing, (key, group, missing)
            raw[key][group] = np.array([by_age[a] for a in AGES], dtype=float)
    assert set(raw) == {(c, s) for c in COUNTRIES for s in SEXES.values()}, set(raw)
    return raw


RAW = load_a20(A20_FILE)

COMPARATORS = {"Pooled arterial": ARTERIAL, "Coronary": ["Coronary"],
               "Cerebrovascular": ["Cerebrovascular"],
               "Other arterial": ["Other arterial"]}

COLORS = {"Venous": "#0072B2", "Coronary": "#D55E00",
          "Cerebrovascular": "#E69F00", "Other arterial": "#009E73",
          "Pooled arterial": "#555555"}

# Chile, ages 20-29: female-to-male rate ratios from the Chile rate analysis
# (scripts/03_rates.py), shown rounded in the explanatory box. Their ratio is
# checked against the case-case OR computed from Table A20.
CHILE_2029_VENOUS_RR = 1.53
CHILE_2029_CORONARY_RR = 0.23


# ---------------------------------------------------------------------------
# Estimates
# ---------------------------------------------------------------------------
def counts(country, sex):
    """Counts by age for the main contrast (venous = I26 + I80-I82; I83 excluded)."""
    col = RAW[(country, sex)]
    return {"Venous": col["PE"] + col["I80_82"],
            "Coronary": col["Coronary"],
            "Cerebrovascular": col["Cerebrovascular"],
            "Other arterial": col["Other arterial"]}


def cells(country, comparator):
    """2x2 cells by age: a, b = women venous, arterial; c, d = men venous, arterial."""
    w, m = counts(country, "W"), counts(country, "M")
    cols = COMPARATORS[comparator]
    return (w["Venous"], sum(w[k] for k in cols),
            m["Venous"], sum(m[k] for k in cols))


def stratum_or(country, comparator):
    """Female-to-male case-case OR in each age stratum."""
    a, b, c, d = cells(country, comparator)
    return (a * d) / (b * c)


def mh_or(country, comparator):
    """Mantel-Haenszel age-adjusted case-case OR."""
    a, b, c, d = cells(country, comparator)
    n = a + b + c + d
    return (a * d / n).sum() / (b * c / n).sum()


N_TOTAL = {c: int(sum(v.sum() for s in ("W", "M") for v in counts(c, s).values()))
           for c in COUNTRIES}
OR = {c: {k: stratum_or(c, k) for k in COMPARATORS} for c in COUNTRIES}
MH = {c: {k: mh_or(c, k) for k in COMPARATORS} for c in COUNTRIES}
OR_CHILE_2029_CORONARY = OR["Chile"]["Coronary"][AGES.index("20-29")]


def check_against_manuscript():
    """Stop if the transcribed counts or the figure's statements do not hold."""
    assert N_TOTAL == {"Chile": 834_780, "Brazil": 7_409_594}, N_TOTAL

    i2029, i4049, i80 = AGES.index("20-29"), AGES.index("40-49"), AGES.index("≥80")
    assert round(OR["Chile"]["Coronary"][i2029], 2) == 6.69
    assert round(OR["Brazil"]["Coronary"][i2029], 2) == 6.62
    assert round(OR["Chile"]["Pooled arterial"][i4049], 2) == 2.26
    assert round(OR["Brazil"]["Pooled arterial"][i80], 2) == 1.61

    for country, comp, expected in [("Chile", "Pooled arterial", 1.80),
                                    ("Brazil", "Pooled arterial", 2.20),
                                    ("Chile", "Coronary", 2.54),
                                    ("Brazil", "Coronary", 2.61)]:
        assert abs(MH[country][comp] - expected) < 0.006, (country, comp,
                                                           MH[country][comp])

    # "Larger venous share in women in every age group, in both countries"
    for country in COUNTRIES:
        for comp in COMPARATORS:
            assert (OR[country][comp] > 1).all(), (country, comp)

    # "Largest vs coronary, smallest vs cerebrovascular" (age-adjusted ORs)
    single = ARTERIAL
    for country in COUNTRIES:
        ranked = sorted(single, key=lambda k: MH[country][k])
        assert ranked[0] == "Cerebrovascular" and ranked[-1] == "Coronary", (
            country, {k: MH[country][k] for k in single})

    # Rounded rate ratios in the explanatory box are consistent with the OR
    lo = (CHILE_2029_VENOUS_RR - 0.005) / (CHILE_2029_CORONARY_RR + 0.005)
    hi = (CHILE_2029_VENOUS_RR + 0.005) / (CHILE_2029_CORONARY_RR - 0.005)
    assert lo <= OR_CHILE_2029_CORONARY <= hi, (lo, OR_CHILE_2029_CORONARY, hi)


# ---------------------------------------------------------------------------
# Country outlines (optional)
# ---------------------------------------------------------------------------
def load_shapes(path, names):
    if not path.exists():
        print(f"Country outlines not found ({path}); drawing without them.")
        return {}
    gj = json.loads(path.read_text(encoding="utf-8"))
    shapes = {}
    for f in gj.get("features", []):
        name = f.get("properties", {}).get("ADMIN")
        if name in names:
            g = f["geometry"]
            polys = (g["coordinates"] if g["type"] == "MultiPolygon"
                     else [g["coordinates"]])
            shapes[name] = [np.array(poly[0])[:, :2] for poly in polys]
    return shapes


# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------
def box(ax, x, y, w, h, fc="white", ec="0.6", lw=1.0):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle="round,pad=0,rounding_size=1.0",
                                fc=fc, ec=ec, lw=lw))


def arrow(ax, x1, x2, y):
    ax.add_patch(FancyArrowPatch((x1, y), (x2, y), arrowstyle="-|>",
                                 mutation_scale=18, color="0.4", lw=2))


def silhouette(ax, shapes, name, bounds, color="0.45"):
    if name not in shapes:
        return
    ins = ax.inset_axes(bounds, transform=ax.transData)
    for ring in shapes[name]:
        ins.fill(ring[:, 0], ring[:, 1], color=color, lw=0)
    ins.set_aspect("equal")
    ins.axis("off")


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
def draw(shapes):
    fig, ax = plt.subplots(figsize=(14, 7.84))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 56)
    ax.set_aspect("equal")
    ax.axis("off")

    for x, t in [(1, "1  Data sources"), (26, "2  Diagnostic groups"),
                 (50, "3  Estimator"), (78, "4  Shared pattern")]:
        ax.text(x, 52.5, t, fontsize=12, fontweight="bold", va="center")

    for x1, x2 in [(22.6, 25.4), (46.6, 49.4), (74.6, 77.4)]:
        arrow(ax, x1, x2, 28)

    # 1. Data sources
    sources = {
        "Chile": dict(y=29, sil=[2, 30.5, 5, 17],
                      lines=["DEIS, 2015-2025", "Individual records",
                             "Public + private"], unit="discharges"),
        "Brazil": dict(y=6, sil=[1.5, 9, 8, 14],
                       lines=["SIH/SUS, 2015-2025", "Aggregate counts",
                              "SUS-financed only"], unit="hospitalizations"),
    }
    for country, s in sources.items():
        y0 = s["y"]
        box(ax, 1, y0, 21, 20)
        silhouette(ax, shapes, country, s["sil"])
        ax.text(10, y0 + 17, country, fontsize=12, fontweight="bold", va="center")
        for i, line in enumerate(s["lines"]):
            ax.text(10, y0 + 14 - 2.5 * i, line, fontsize=8.5, va="center")
        ax.text(10, y0 + 5.5, f"{N_TOTAL[country]:,}", fontsize=10.5,
                fontweight="bold", va="center")
        ax.text(10, y0 + 3, s["unit"], fontsize=8.5, va="center")

    # 2. Diagnostic groups
    box(ax, 26, 38, 20, 11, fc=COLORS["Venous"], ec="none")
    ax.text(36, 46, "Selected venous", fontsize=10.5, fontweight="bold",
            color="white", ha="center", va="center")
    ax.text(36, 43, "I26, I80-I82", fontsize=9, color="white",
            ha="center", va="center")
    ax.text(36, 40.3, "PE, phlebitis, thrombosis", fontsize=8, color="white",
            ha="center", va="center")
    ax.text(36, 35.5, "versus", fontsize=9, style="italic", color="0.3",
            ha="center", va="center")

    for name, codes, y0, tc in [("Coronary", "I20-I25", 25, "white"),
                                ("Cerebrovascular", "I60-I69", 16, "black"),
                                ("Other arterial", "I70-I79", 7, "white")]:
        box(ax, 26, y0, 20, 8, fc=COLORS[name], ec="none")
        ax.text(36, y0 + 5.2, name, fontsize=10.5, fontweight="bold", color=tc,
                ha="center", va="center")
        ax.text(36, y0 + 2.4, codes, fontsize=9, color=tc,
                ha="center", va="center")

    # 3. Estimator
    box(ax, 50, 31, 24, 18)
    ax.text(62, 46, "Case-case odds ratio", fontsize=10.5, fontweight="bold",
            ha="center", va="center")
    ax.text(62, 40,
            r"$\mathrm{OR}=\frac{(\mathrm{venous}\,/\,\mathrm{arterial})_{\mathrm{women}}}"
            r"{(\mathrm{venous}\,/\,\mathrm{arterial})_{\mathrm{men}}}$",
            fontsize=15, ha="center", va="center")
    ax.text(62, 34, "Within each of seven age strata;\n"
            "Mantel-Haenszel age-adjusted summary",
            fontsize=8, ha="center", va="center", color="0.25")

    box(ax, 50, 6, 24, 21, fc="#FFF8E1", ec="#C9A227")
    ax.text(62, 24, "Reading the OR", fontsize=10.5, fontweight="bold",
            ha="center", va="center")
    ax.text(62, 20.3, "OR = venous rate ratio\n÷ arterial rate ratio",
            fontsize=8.5, ha="center", va="center", linespacing=1.3)
    ax.text(62, 16.4, "Chile, ages 20-29, venous vs coronary:", fontsize=8.5,
            ha="center", va="center")
    ax.text(62, 13.2,
            f"{CHILE_2029_VENOUS_RR:.2f} ÷ {CHILE_2029_CORONARY_RR:.2f} "
            f"≈ {OR_CHILE_2029_CORONARY:.2f}",
            fontsize=13, fontweight="bold", ha="center", va="center")
    ax.text(62, 9.0, "Moderately higher venous rate and\n"
            "much lower coronary rate; not a sixfold\n"
            "venous excess. Unrounded rates used.",
            fontsize=7.5, ha="center", va="center", color="0.25")

    # 4. Shared pattern
    box(ax, 78, 6, 21, 43)
    ins = ax.inset_axes([86.5, 24, 11.5, 19], transform=ax.transData)
    rows = ["Coronary", "Other arterial", "Cerebrovascular", "Pooled arterial"]
    row_labels = ["Coronary", "Other arterial", "Cerebrovasc.", "Pooled arterial"]
    dodge = 0.17  # Chile above, Brazil below the row, so close values stay visible
    for i, comp in enumerate(rows):
        ch, br = MH["Chile"][comp], MH["Brazil"][comp]
        ins.plot([ch, br], [i - dodge, i + dodge], "-", color=COLORS[comp],
                 lw=1.5, alpha=0.5)
        ins.scatter(ch, i - dodge, s=45, marker="o", color=COLORS[comp],
                    zorder=3)
        ins.scatter(br, i + dodge, s=45, marker="s", facecolor="white",
                    edgecolor=COLORS[comp], linewidth=1.6, zorder=3)
        ins.text(max(ch, br) + 0.14, i, f"{ch:.2f} / {br:.2f}", fontsize=7,
                 va="center")
    ins.scatter([], [], s=45, marker="o", color="0.3", label="Chile")
    ins.scatter([], [], s=45, marker="s", facecolor="white", edgecolor="0.3",
                linewidth=1.8, label="Brazil")
    ins.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2,
               frameon=False, fontsize=8, handletextpad=0.2, columnspacing=1.0)
    ins.axvline(1, color="0.5", ls="--", lw=1)
    ins.set_xlim(0.9, 3.9)
    ins.set_xticks([1, 2, 3])
    ins.set_ylim(len(rows) - 0.5, -0.5)
    ins.set_yticks(np.arange(len(rows)))
    ins.set_yticklabels(row_labels, fontsize=8)
    ins.tick_params(axis="x", labelsize=8)
    ins.tick_params(axis="y", length=0)
    ins.set_xlabel("Age-adjusted OR", fontsize=8)
    ins.spines["left"].set_visible(False)

    ax.text(88.5, 14.5, "Larger venous share in women\nin every age group,\n"
            "in both countries", fontsize=9, fontweight="bold",
            ha="center", va="center")
    ax.text(88.5, 9, "Largest vs coronary,\nsmallest vs cerebrovascular",
            fontsize=8, ha="center", va="center", color="0.25")

    # Footer
    ax.text(50, 2.5, "Units are discharge events, not individuals. Odds ratios "
            "describe diagnostic composition among selected vascular discharges, "
            "not population risk.", fontsize=8, ha="center", va="center",
            color="0.3")
    return fig


def main():
    check_against_manuscript()
    fig = draw(load_shapes(SHAPES_FILE, COUNTRIES))
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_DIR / f"{OUT_NAME}.png", dpi=300, bbox_inches="tight")
    # No creation date in the PDF, so re-running does not produce a spurious diff
    fig.savefig(OUT_DIR / f"{OUT_NAME}.pdf", bbox_inches="tight",
                metadata={"CreationDate": None})
    plt.close(fig)
    print(f"Graphical abstract written to {OUT_DIR.relative_to(ROOT)}/"
          f"{OUT_NAME}.png and .pdf")


if __name__ == "__main__":
    main()
