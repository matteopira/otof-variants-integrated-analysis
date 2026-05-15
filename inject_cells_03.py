"""
Inject new analysis cells into notebook 03-domain-mapping.ipynb
"""
import json
import uuid

NB_PATH = "notebooks/03-domain-mapping.ipynb"

with open(NB_PATH) as f:
    nb = json.load(f)

def md_cell(source: str) -> dict:
    return {
        "cell_type": "markdown",
        "id": str(uuid.uuid4())[:36],
        "metadata": {},
        "source": source,
    }

def code_cell(source: str) -> dict:
    return {
        "cell_type": "code",
        "execution_count": None,
        "id": str(uuid.uuid4())[:36],
        "metadata": {},
        "outputs": [],
        "source": source,
    }

# ===========================================================================
# TASK 1 — Fisher's exact test per domain
# ===========================================================================

MD_TASK1 = (
    "## Fisher's Exact Test per Domain vs Linker\n\n"
    "The permutation test used in the preceding section draws randomly from the pool of all "
    "named-domain residues, so the reference distribution is the average density across the "
    "entire set of C2 and TM domains. This makes the null hypothesis 'domain X is like all "
    "other domains combined.'\n\n"
    "Fisher's exact test instead uses the linker region as a specific, well-defined background. "
    "For each domain we ask: is the ratio of P/LP to non-P/LP variants in this domain different "
    "from the same ratio in the linker? This gives an odds ratio directly interpretable as "
    "enrichment (OR > 1) or depletion (OR < 1) relative to the linker, which is the largest "
    "continuous region with a well-sampled pathogenic density. Bonferroni correction is applied "
    "for the 7 domain comparisons."
)

CODE_TASK1 = """\
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path

results_dir = Path("../results")
results_dir.mkdir(exist_ok=True)

# Domain lengths (UniProt Q9HC10, canonical isoform)
domain_lengths = {
    "C2A": 122, "C2B": 121, "C2C": 116,
    "C2D": 115, "C2E": 116, "C2F": 117, "TM": 32,
}

# Observed P/LP counts per domain (from domain mapping analysis above)
plp_counts = {
    "C2A": 8, "C2B": 20, "C2C": 18,
    "C2D": 18, "C2E": 15, "C2F": 15, "TM": 3,
}

# Linker: 158 P/LP, length = PROTEIN_LENGTH - sum(domain lengths)
PROTEIN_LENGTH = 1997
linker_length = PROTEIN_LENGTH - sum(domain_lengths.values())
linker_plp    = 158
linker_non    = linker_length - linker_plp

print(f"Linker length: {linker_length} aa | P/LP: {linker_plp} | non-P/LP: {linker_non}")
print()

# Run Fisher's exact test for each domain vs linker
N_DOMAINS = len(domain_lengths)
rows = []

for domain in domain_lengths:
    d_plp  = plp_counts[domain]
    d_len  = domain_lengths[domain]
    d_non  = d_len - d_plp

    table = [[d_plp, d_non],
             [linker_plp, linker_non]]

    odds_ratio, p_raw = stats.fisher_exact(table, alternative="two-sided")
    p_bonf = min(p_raw * N_DOMAINS, 1.0)
    sig = ("***" if p_bonf < 0.001 else
           "**"  if p_bonf < 0.01  else
           "*"   if p_bonf < 0.05  else "n.s.")

    rows.append({
        "Domain": domain,
        "P/LP": d_plp,
        "Length_aa": d_len,
        "Odds_Ratio": round(odds_ratio, 3),
        "p_raw": round(p_raw, 5),
        "p_Bonferroni": round(p_bonf, 5),
        "Significance": sig,
    })

fisher_df = pd.DataFrame(rows)

header = f"{'Domain':6s} {'P/LP':>5s} {'Length':>7s} {'Odds Ratio':>12s} {'p (raw)':>10s} {'p (Bonf)':>11s} {'Sig':>6s}"
print(header)
print("-" * 60)
for _, r in fisher_df.iterrows():
    line = (f"{r['Domain']:6s} {int(r['P/LP']):>5d} {int(r['Length_aa']):>7d} "
            f"{r['Odds_Ratio']:>12.3f} {r['p_raw']:>10.5f} {r['p_Bonferroni']:>11.5f} {r['Significance']:>6s}")
    print(line)

out_csv = results_dir / "domain_fisher_exact.csv"
fisher_df.to_csv(out_csv, index=False)
print()
print(f"Saved: {out_csv}")

print()
print("Interpretation:")
for _, r in fisher_df.iterrows():
    direction = "enriched" if r["Odds_Ratio"] > 1 else "depleted"
    print(f"  {r['Domain']}: OR={r['Odds_Ratio']:.3f} ({direction}) | "
          f"Bonferroni p={r['p_Bonferroni']:.4f} {r['Significance']}")

linker_density = linker_plp / linker_length
print()
print(f"Linker reference density: {linker_density:.4f} P/LP per aa")
print("Note: OR < 1 indicates depletion relative to linker; OR > 1 indicates enrichment.")
"""

# ===========================================================================
# TASK 2 — Fraction pathogenic per domain
# ===========================================================================

MD_TASK2 = (
    "## Fraction Pathogenic per Domain: Ascertainment-Corrected Metric\n\n"
    "Absolute pathogenic variant density (P/LP variants per amino acid) is confounded by "
    "ascertainment bias: ClinVar coverage is uneven across domains because submitters and "
    "researchers preferentially study variants in functionally characterized regions. A domain "
    "with many submitted variants of any class will appear denser in all categories.\n\n"
    "Fraction pathogenic -- defined as P/LP / (P/LP + B/LB + VUS + Conflicting) -- conditions "
    "on the total number of submitted variants per domain, controlling for ascertainment. It asks: "
    "of all variants that have been classified in this domain, what proportion are pathogenic? "
    "Wilson confidence intervals account for the binomial uncertainty in each domain's estimate."
)

CODE_TASK2A = """\
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
from pathlib import Path
from scipy import stats

results_dir = Path("../results")
data_dir    = Path("../data")

cols_drop = [
    "Somatic clinical impact",
    "Somatic clinical impact date last evaluated",
    "Somatic clinical impact review status",
    "Oncogenicity classification",
    "Oncogenicity date last evaluated",
    "Oncogenicity review status",
    "Unnamed: 24",
]
clinvar = pd.read_csv(data_dir / "clinvar_result.txt", sep="\\t")
clinvar = clinvar.drop(columns=[c for c in cols_drop if c in clinvar.columns])

def extract_aa_pos(name):
    if pd.isna(name):
        return np.nan
    text = str(name)
    embedded = re.search(r"\\(p\\.([^)]+)\\)", text)
    if embedded:
        text = "p." + embedded.group(1)
    m = re.search(r"p\\.(?:[A-Za-z]{1,3})(\\d+)", text)
    if m:
        try:
            return int(m.group(1))
        except ValueError:
            pass
    return np.nan

clinvar["aa_position"] = clinvar["Name"].apply(extract_aa_pos)

domains_list = [
    {"name": "C2A", "start": 1,    "end": 122,  "color": "#4C72B0"},
    {"name": "C2B", "start": 360,  "end": 480,  "color": "#55A868"},
    {"name": "C2C", "start": 481,  "end": 596,  "color": "#C44E52"},
    {"name": "C2D", "start": 940,  "end": 1054, "color": "#8172B2"},
    {"name": "C2E", "start": 1158, "end": 1273, "color": "#CCB974"},
    {"name": "C2F", "start": 1481, "end": 1597, "color": "#64B5CD"},
    {"name": "TM",  "start": 1942, "end": 1973, "color": "#777777"},
]
domain_color_map = {d["name"]: d["color"] for d in domains_list}
domain_color_map["Linker"] = "#999999"

def assign_domain(pos):
    if pd.isna(pos):
        return None
    for d in domains_list:
        if d["start"] <= pos <= d["end"]:
            return d["name"]
    return "Linker"

clinvar["domain"] = clinvar["aa_position"].apply(assign_domain)

pathogenic_classes = ["Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"]
benign_classes     = ["Benign", "Likely benign", "Benign/Likely benign"]
vus_classes        = ["Uncertain significance"]
conflict_classes   = ["Conflicting classifications of pathogenicity"]

def classify(cls):
    if cls in pathogenic_classes: return "PLP"
    if cls in benign_classes:     return "BLB"
    if cls in vus_classes:        return "VUS"
    if cls in conflict_classes:   return "Conflict"
    return "Other"

clinvar["cat"] = clinvar["Germline classification"].apply(classify)

domain_order = ["C2A", "C2B", "C2C", "C2D", "C2E", "C2F", "TM", "Linker"]
rows = []
for dom in domain_order:
    sub     = clinvar[clinvar["domain"] == dom]
    n_total = len(sub)
    n_plp   = (sub["cat"] == "PLP").sum()
    n_blb   = (sub["cat"] == "BLB").sum()
    n_vus   = (sub["cat"] == "VUS").sum()
    n_conf  = (sub["cat"] == "Conflict").sum()
    n_other = (sub["cat"] == "Other").sum()
    denom   = n_plp + n_blb + n_vus + n_conf
    frac    = n_plp / denom if denom > 0 else np.nan
    rows.append({
        "Domain": dom,
        "Total_variants": n_total,
        "PLP": n_plp, "BLB": n_blb, "VUS": n_vus,
        "Conflicting": n_conf, "Other": n_other,
        "Denominator": denom,
        "Fraction_pathogenic": frac,
    })

frac_df = pd.DataFrame(rows)

def wilson_ci(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p  = k / n
    dn = 1 + z**2 / n
    c  = (p + z**2 / (2*n)) / dn
    margin = z * np.sqrt(p*(1-p)/n + z**2/(4*n**2)) / dn
    return (max(0.0, c - margin), min(1.0, c + margin))

frac_df["CI_lo"], frac_df["CI_hi"] = zip(
    *frac_df.apply(lambda r: wilson_ci(r["PLP"], r["Denominator"]), axis=1)
)
frac_df["CI_err_lo"] = frac_df["Fraction_pathogenic"] - frac_df["CI_lo"]
frac_df["CI_err_hi"] = frac_df["CI_hi"]  - frac_df["Fraction_pathogenic"]

print("Fraction pathogenic per domain:")
print(frac_df[["Domain", "PLP", "Denominator", "Fraction_pathogenic", "CI_lo", "CI_hi"]].to_string(index=False))

out_csv = results_dir / "domain_fraction_pathogenic.csv"
frac_df.to_csv(out_csv, index=False)
print()
print(f"Saved: {out_csv}")
"""

CODE_TASK2B = """\
import matplotlib.pyplot as plt
import numpy as np

plot_df = frac_df.dropna(subset=["Fraction_pathogenic"]).sort_values("Fraction_pathogenic")

fig, ax = plt.subplots(figsize=(9, 6))

colors  = [domain_color_map.get(d, "#aaaaaa") for d in plot_df["Domain"]]
xerr_lo = plot_df["CI_err_lo"].values
xerr_hi = plot_df["CI_err_hi"].values

ax.barh(
    plot_df["Domain"],
    plot_df["Fraction_pathogenic"],
    xerr=[xerr_lo, xerr_hi],
    color=colors,
    edgecolor="black",
    linewidth=0.6,
    capsize=4,
    error_kw={"elinewidth": 1.2, "capthick": 1.2, "ecolor": "#333333"},
    height=0.65,
)

linker_frac = frac_df.loc[frac_df["Domain"] == "Linker", "Fraction_pathogenic"].values[0]
ax.axvline(linker_frac, color="#333333", linestyle="--", linewidth=1.5,
           label=f"Linker reference ({linker_frac:.3f})")

for i, (_, row) in enumerate(plot_df.iterrows()):
    ax.text(
        row["Fraction_pathogenic"] + row["CI_err_hi"] + 0.003,
        i,
        f"{int(row['PLP'])}/{int(row['Denominator'])}",
        va="center", ha="left", fontsize=9, color="#333333",
    )

ax.set_xlabel("Fraction pathogenic (P/LP / classified variants)", fontsize=11)
ax.set_ylabel("Domain", fontsize=11)
ax.set_title(
    "OTOF: fraction pathogenic per domain\\n"
    "Error bars = 95% Wilson CI  |  Dashed line = linker reference",
    fontsize=12, fontweight="bold",
)
ax.legend(fontsize=9, loc="lower right")
ax.set_xlim(0, (plot_df["Fraction_pathogenic"] + plot_df["CI_err_hi"]).max() * 1.35)
ax.spines[["top", "right"]].set_visible(False)

plt.tight_layout()
out_png = results_dir / "domain_fraction_pathogenic.png"
plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.show()
print(f"Saved: {out_png}")
"""

# ===========================================================================
# TASK 3 — Not-provided distribution per domain
# ===========================================================================

MD_TASK3 = (
    "## Not-Provided Condition Annotation: Domain-Level Bias Check\n\n"
    "A prior analysis (notebook 01) showed that 63.6% of OTOF ClinVar records lack a specific "
    "condition annotation ('not provided'). If this missingness were concentrated in particular "
    "domains, domain-level analyses could be confounded: a domain with many uncharacterized "
    "submissions might appear artificially benign or pathogenic depending on the classification "
    "of those unannotated records.\n\n"
    "The chi-square test of independence below tests whether the probability of lacking a "
    "condition annotation is uniform across domains. If not significant, the not-provided pattern "
    "is homogeneous and domain-density conclusions are robust to this potential confound."
)

CODE_TASK3A = """\
import pandas as pd
import numpy as np
import re
from scipy import stats
from pathlib import Path

results_dir = Path("../results")
data_dir    = Path("../data")

cols_drop = [
    "Somatic clinical impact",
    "Somatic clinical impact date last evaluated",
    "Somatic clinical impact review status",
    "Oncogenicity classification",
    "Oncogenicity date last evaluated",
    "Oncogenicity review status",
    "Unnamed: 24",
]
clinvar_np = pd.read_csv(data_dir / "clinvar_result.txt", sep="\\t")
clinvar_np = clinvar_np.drop(columns=[c for c in cols_drop if c in clinvar_np.columns])

def extract_aa_pos(name):
    if pd.isna(name):
        return np.nan
    text = str(name)
    embedded = re.search(r"\\(p\\.([^)]+)\\)", text)
    if embedded:
        text = "p." + embedded.group(1)
    m = re.search(r"p\\.(?:[A-Za-z]{1,3})(\\d+)", text)
    if m:
        try:
            return int(m.group(1))
        except ValueError:
            pass
    return np.nan

clinvar_np["aa_position"] = clinvar_np["Name"].apply(extract_aa_pos)

domains_np = [
    {"name": "C2A", "start": 1,    "end": 122},
    {"name": "C2B", "start": 360,  "end": 480},
    {"name": "C2C", "start": 481,  "end": 596},
    {"name": "C2D", "start": 940,  "end": 1054},
    {"name": "C2E", "start": 1158, "end": 1273},
    {"name": "C2F", "start": 1481, "end": 1597},
    {"name": "TM",  "start": 1942, "end": 1973},
]

def assign_domain_np(pos):
    if pd.isna(pos):
        return None
    for d in domains_np:
        if d["start"] <= pos <= d["end"]:
            return d["name"]
    return "Linker"

clinvar_np["domain"] = clinvar_np["aa_position"].apply(assign_domain_np)
cv_mapped_np = clinvar_np[clinvar_np["domain"].notna()].copy()

NOT_PROVIDED_TERMS = {"not provided", "not specified", "phenotype not yet classified"}

def is_not_provided(cond):
    if pd.isna(cond):
        return True
    return str(cond).strip().lower() in NOT_PROVIDED_TERMS

cv_mapped_np["not_provided"] = cv_mapped_np["Condition(s)"].apply(is_not_provided)

domain_order = ["C2A", "C2B", "C2C", "C2D", "C2E", "C2F", "TM", "Linker"]
rows_np = []
for dom in domain_order:
    sub       = cv_mapped_np[cv_mapped_np["domain"] == dom]
    n_total   = len(sub)
    n_prov    = (~sub["not_provided"]).sum()
    n_np      = sub["not_provided"].sum()
    frac_np   = n_np / n_total if n_total > 0 else np.nan
    rows_np.append({
        "Domain": dom,
        "Total": n_total,
        "Provided": n_prov,
        "Not_provided": n_np,
        "Frac_not_provided": frac_np,
    })

bias_df = pd.DataFrame(rows_np)
print("Not-provided fraction per domain:")
print(bias_df.to_string(index=False))

contingency = bias_df[["Provided", "Not_provided"]].values
chi2_np, pval_np, dof_np, _ = stats.chi2_contingency(contingency)
print()
print(f"Chi-square test of independence: chi2={chi2_np:.3f}, df={dof_np}, p={pval_np:.4f}")

out_csv = results_dir / "domain_not_provided_bias.csv"
bias_df.to_csv(out_csv, index=False)
print(f"Saved: {out_csv}")
"""

CODE_TASK3B = """\
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(figsize=(10, 5))

x = np.arange(len(bias_df))
prov_frac  = bias_df["Provided"]     / bias_df["Total"]
notpr_frac = bias_df["Not_provided"] / bias_df["Total"]

ax.bar(x, prov_frac,  label="Condition provided", color="#55A868", edgecolor="black", linewidth=0.5)
ax.bar(x, notpr_frac, bottom=prov_frac, label="Not provided",
       color="#CCCCCC", edgecolor="black", linewidth=0.5)

ax.set_xticks(x)
ax.set_xticklabels(bias_df["Domain"], fontsize=10)
ax.set_ylabel("Fraction of variants", fontsize=11)
ax.set_xlabel("Domain", fontsize=11)
ax.set_ylim(0, 1.10)

sig_str = f"p = {pval_np:.4f}" if pval_np >= 0.0001 else "p < 0.0001"
ax.set_title(
    f"Not-provided condition annotation by domain\\n"
    f"Chi-square: chi2={chi2_np:.2f}, df={dof_np}, {sig_str}",
    fontsize=11, fontweight="bold",
)
ax.legend(fontsize=9, loc="upper right")
ax.spines[["top", "right"]].set_visible(False)

for i, (_, row) in enumerate(bias_df.iterrows()):
    ax.text(i, 1.02, f"n={int(row['Total'])}", ha="center", va="bottom",
            fontsize=8, color="#333333")

plt.tight_layout()
out_png = results_dir / "domain_not_provided_bias.png"
plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.show()
print(f"Saved: {out_png}")

print()
print("Interpretation:")
if pval_np >= 0.05:
    print("  Chi-square is NOT significant (p >= 0.05).")
    print("  The not-provided fraction is uniform across domains.")
    print("  Domain-density conclusions are robust to this confound.")
else:
    print(f"  Chi-square IS significant (p = {pval_np:.4f}).")
    print("  The not-provided fraction differs by domain -- interpret domain densities with caution.")
    max_dom = bias_df.loc[bias_df['Frac_not_provided'].idxmax(), 'Domain']
    print(f"  Domain with highest not-provided fraction: {max_dom}")
"""

# ===========================================================================
# TASK 5 — HGVS non-parseable sensitivity analysis
# ===========================================================================

MD_TASK5 = (
    "## HGVS Non-Parseable Sensitivity Analysis\n\n"
    "The robust HGVS parser defined above maps 65.8% of ClinVar records (1,600 of 2,432) to "
    "an amino acid position. The remaining 34.2% (832 records) could not be mapped and are "
    "therefore excluded from all domain-level analyses. This exclusion is only acceptable if "
    "the non-parseable variants are not systematically biased by functional category or clinical "
    "classification.\n\n"
    "The expected reason for non-parseability is that many ClinVar entries describe intronic or "
    "non-coding variants (e.g., splice site changes named with c-dot notation) that have no "
    "protein consequence notation. If so, non-parseable variants should be enriched in "
    "intronic/splicing consequences but should NOT be biased toward any particular germline "
    "classification -- otherwise domain results could overestimate or underestimate pathogenic "
    "density in a non-random way."
)

CODE_TASK5A = """\
import pandas as pd
import numpy as np
import re
from pathlib import Path

results_dir = Path("../results")
data_dir    = Path("../data")

cols_drop = [
    "Somatic clinical impact",
    "Somatic clinical impact date last evaluated",
    "Somatic clinical impact review status",
    "Oncogenicity classification",
    "Oncogenicity date last evaluated",
    "Oncogenicity review status",
    "Unnamed: 24",
]
clinvar_s = pd.read_csv(data_dir / "clinvar_result.txt", sep="\\t")
clinvar_s = clinvar_s.drop(columns=[c for c in cols_drop if c in clinvar_s.columns])

def extract_aa_pos(name):
    if pd.isna(name):
        return np.nan
    text = str(name)
    embedded = re.search(r"\\(p\\.([^)]+)\\)", text)
    if embedded:
        text = "p." + embedded.group(1)
    m = re.search(r"p\\.(?:[A-Za-z]{1,3})(\\d+)", text)
    if m:
        try:
            return int(m.group(1))
        except ValueError:
            pass
    return np.nan

clinvar_s["aa_position"] = clinvar_s["Name"].apply(extract_aa_pos)
clinvar_s["parseable"]   = clinvar_s["aa_position"].notna()

n_total   = len(clinvar_s)
n_parse   = clinvar_s["parseable"].sum()
n_nonpar  = n_total - n_parse
print(f"Total ClinVar variants: {n_total}")
print(f"Parseable:              {n_parse} ({100*n_parse/n_total:.1f}%)")
print(f"Non-parseable:          {n_nonpar} ({100*n_nonpar/n_total:.1f}%)")

mc_col = "Molecular consequence"
if mc_col not in clinvar_s.columns:
    mc_col = next((c for c in clinvar_s.columns if "consequence" in c.lower()), None)
print()
print(f"Consequence column: {mc_col}")
if mc_col:
    print(clinvar_s[mc_col].value_counts().head(10))
"""

CODE_TASK5B = """\
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
from scipy import stats

def simplify_consequence(cons):
    if pd.isna(cons):
        return "Unknown"
    c = str(cons).lower()
    if "missense" in c:
        return "Missense"
    if "stop" in c or "nonsense" in c:
        return "Nonsense/stop"
    if "frameshift" in c or "frame" in c:
        return "Frameshift"
    if "splice" in c or "intron" in c:
        return "Splice/intronic"
    if "synonymous" in c or "silent" in c:
        return "Synonymous"
    if "utr" in c or "5'" in c or "3'" in c:
        return "UTR"
    if "deletion" in c or "duplication" in c or "copy" in c:
        return "CNV/indel"
    return "Other"

mc_col = "Molecular consequence"
if mc_col not in clinvar_s.columns:
    mc_col = next((c for c in clinvar_s.columns if "consequence" in c.lower()), None)

if mc_col:
    clinvar_s["cons_simple"] = clinvar_s[mc_col].apply(simplify_consequence)
else:
    def infer_from_name(name):
        if pd.isna(name):
            return "Unknown"
        n = str(name).lower()
        if "p." in n and "ter" not in n and "*" not in n:
            return "Missense"
        if "ter" in n or "*" in n:
            return "Nonsense/stop"
        if "c." in n and "p." not in n:
            return "Splice/intronic"
        return "Other"
    clinvar_s["cons_simple"] = clinvar_s["Name"].apply(infer_from_name)

def simplify_cls(cls):
    if pd.isna(cls): return "Other"
    if cls in ["Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"]:
        return "Pathogenic"
    if cls in ["Benign", "Likely benign", "Benign/Likely benign"]:
        return "Benign"
    if cls == "Uncertain significance":
        return "VUS"
    if cls == "Conflicting classifications of pathogenicity":
        return "Conflicting"
    return "Other"

clinvar_s["cls_simple"] = clinvar_s["Germline classification"].apply(simplify_cls)

cons_ct = pd.crosstab(clinvar_s["cons_simple"], clinvar_s["parseable"])
cls_ct  = pd.crosstab(clinvar_s["cls_simple"],  clinvar_s["parseable"])

chi2_cons, p_cons = stats.chi2_contingency(cons_ct)[:2]
chi2_cls,  p_cls  = stats.chi2_contingency(cls_ct)[:2]

print(f"Consequence vs parseability: chi2={chi2_cons:.2f}, p={p_cons:.4e}")
print(f"Classification vs parseability: chi2={chi2_cls:.2f}, p={p_cls:.4e}")

cons_frac = cons_ct.div(cons_ct.sum(axis=1), axis=0)
cls_frac  = cls_ct.div(cls_ct.sum(axis=1), axis=0)
cons_frac.columns = ["Non-parseable" if not c else "Parseable" for c in cons_frac.columns]
cls_frac.columns  = ["Non-parseable" if not c else "Parseable" for c in cls_frac.columns]
"""

CODE_TASK5C = """\
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

results_dir = Path("../results")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
w = 0.35

# Left: consequence
ax = axes[0]
cs = cons_frac.sort_values("Parseable", ascending=False)
x  = np.arange(len(cs))
ax.bar(x - w/2, cs["Non-parseable"], width=w, label="Non-parseable",
       color="#CC6677", edgecolor="black", linewidth=0.5)
ax.bar(x + w/2, cs["Parseable"],     width=w, label="Parseable",
       color="#4477AA", edgecolor="black", linewidth=0.5)
ax.set_xticks(x)
ax.set_xticklabels(cs.index, rotation=35, ha="right", fontsize=9)
ax.set_ylabel("Fraction of variants in group", fontsize=10)
ax.set_title(f"Molecular consequence\\nchi2={chi2_cons:.1f}, p={p_cons:.2e}",
             fontsize=11, fontweight="bold")
ax.legend(fontsize=9)
ax.spines[["top", "right"]].set_visible(False)

# Right: classification
ax = axes[1]
cl = cls_frac.sort_values("Parseable", ascending=False)
x2 = np.arange(len(cl))
ax.bar(x2 - w/2, cl["Non-parseable"], width=w, label="Non-parseable",
       color="#CC6677", edgecolor="black", linewidth=0.5)
ax.bar(x2 + w/2, cl["Parseable"],     width=w, label="Parseable",
       color="#4477AA", edgecolor="black", linewidth=0.5)
ax.set_xticks(x2)
ax.set_xticklabels(cl.index, rotation=35, ha="right", fontsize=9)
ax.set_ylabel("Fraction of variants in group", fontsize=10)
ax.set_title(f"Germline classification\\nchi2={chi2_cls:.1f}, p={p_cls:.2e}",
             fontsize=11, fontweight="bold")
ax.legend(fontsize=9)
ax.spines[["top", "right"]].set_visible(False)

fig.suptitle(
    "HGVS non-parseable variants: sensitivity analysis\\n"
    "Parseable vs non-parseable variants by molecular consequence and classification",
    fontsize=12, fontweight="bold",
)
plt.tight_layout()
out_png = results_dir / "hgvs_nonparseable_sensitivity.png"
plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.show()
print(f"Saved: {out_png}")
"""

CODE_TASK5D = """\
import pandas as pd
import numpy as np
from pathlib import Path

results_dir = Path("../results")

missense_mask  = clinvar_s["cons_simple"] == "Missense"
n_missense     = missense_mask.sum()
n_miss_parse   = (missense_mask & clinvar_s["parseable"]).sum()
miss_frac_p    = n_miss_parse / n_missense if n_missense > 0 else np.nan

print(f"Missense variants total:     {n_missense}")
print(f"Parseable among missense:    {n_miss_parse} ({100*miss_frac_p:.1f}%)")
print()

print("Interpretation:")
if p_cons < 0.05:
    print(f"  Consequence vs parseability: SIGNIFICANT (p={p_cons:.4e})")
    print("  Non-parseable variants are enriched in certain consequence categories")
    print("  (expected: intronic/splice variants lack protein notation).")
else:
    print(f"  Consequence vs parseability: not significant (p={p_cons:.4e})")

if p_cls < 0.05:
    print(f"  Classification vs parseability: SIGNIFICANT (p={p_cls:.4e})")
    print("  Non-parseable variants differ in classification -- domain results may be biased.")
else:
    print(f"  Classification vs parseability: not significant (p={p_cls:.4e})")
    print("  Classification is homogeneous between parseable and non-parseable variants.")
    print("  Domain-density results are robust to HGVS non-parseability.")

summary_rows = []
for cls_val in clinvar_s["cls_simple"].unique():
    mask = clinvar_s["cls_simple"] == cls_val
    n    = mask.sum()
    np_  = (mask & ~clinvar_s["parseable"]).sum()
    pp   = (mask & clinvar_s["parseable"]).sum()
    summary_rows.append({
        "Classification": cls_val,
        "Total": n,
        "Parseable": pp,
        "NonParseable": np_,
        "FracParseable": pp / n if n > 0 else np.nan,
    })

summary_df = pd.DataFrame(summary_rows).sort_values("Total", ascending=False)
out_csv    = results_dir / "hgvs_nonparseable_summary.csv"
summary_df.to_csv(out_csv, index=False)
print()
print(f"Saved: {out_csv}")
print(summary_df.to_string(index=False))
"""

# ===========================================================================
# TASK 4 — p.Arg1172Gln deep dive
# ===========================================================================

MD_TASK4 = (
    "## p.Arg1172Gln: Deep-Dive Analysis\n\n"
    "p.Arg1172Gln (NM_194248.2:c.3515G>A) is the only OTOF variant with homozygotes in gnomAD "
    "(1 observed homozygote, African/African-American ancestry). The variant sits at position 1172 "
    "(within the C2E domain) and has a complex classification history in ClinVar, with at least one "
    "submission citing 'Likely pathogenic' while gnomAD population frequency is non-trivially high "
    "for a recessive deafness gene. This deep dive integrates (1) all ClinVar submissions for this "
    "variant, (2) gnomAD allele frequencies broken down by ancestry, and (3) ConSurf evolutionary "
    "conservation score for residue 1172 — all from data already present in the repository."
)

CODE_TASK4A = """\
import pandas as pd
import numpy as np
import re
from pathlib import Path

results_dir = Path("../results")
data_dir    = Path("../data")

TARGET_VARIANT = "p.Arg1172Gln"
TARGET_POS     = 1172

# ------------------------------------------------------------------
# 1. ClinVar submissions for this variant
# ------------------------------------------------------------------
cols_drop = [
    "Somatic clinical impact",
    "Somatic clinical impact date last evaluated",
    "Somatic clinical impact review status",
    "Oncogenicity classification",
    "Oncogenicity date last evaluated",
    "Oncogenicity review status",
    "Unnamed: 24",
]
clinvar_dv = pd.read_csv(data_dir / "clinvar_result.txt", sep="\\t")
clinvar_dv = clinvar_dv.drop(columns=[c for c in cols_drop if c in clinvar_dv.columns])

def extract_aa_pos(name):
    if pd.isna(name):
        return np.nan
    text = str(name)
    embedded = re.search(r"\\(p\\.([^)]+)\\)", text)
    if embedded:
        text = "p." + embedded.group(1)
    m = re.search(r"p\\.(?:[A-Za-z]{1,3})(\\d+)", text)
    if m:
        try:
            return int(m.group(1))
        except ValueError:
            pass
    return np.nan

clinvar_dv["aa_position"] = clinvar_dv["Name"].apply(extract_aa_pos)

# Match by position and by name substring
mask_pos  = clinvar_dv["aa_position"] == TARGET_POS
mask_name = clinvar_dv["Name"].str.contains("1172Gln|1172Q|R1172Q", na=False, case=False)
var_cv    = clinvar_dv[mask_pos | mask_name].copy()

print(f"ClinVar entries for {TARGET_VARIANT}: {len(var_cv)}")
if len(var_cv) == 0:
    print("  No exact match found. Checking position 1172 broadly:")
    var_cv = clinvar_dv[mask_pos]
    print(f"  Variants at position 1172: {len(var_cv)}")

if len(var_cv) > 0:
    show_cols = [c for c in [
        "Name", "Germline classification", "Review status",
        "Condition(s)", "Submitter", "Last evaluated",
    ] if c in var_cv.columns]
    print()
    print(var_cv[show_cols].to_string(index=False))
"""

CODE_TASK4B = """\
import pandas as pd
import numpy as np
from pathlib import Path

results_dir = Path("../results")
data_dir    = Path("../data")

TARGET_POS = 1172

# ------------------------------------------------------------------
# 2. gnomAD allele frequencies by ancestry
# ------------------------------------------------------------------
gnomad = pd.read_csv(data_dir / "gnomad_otof_variants.csv")
print(f"gnomAD columns: {list(gnomad.columns)}")
print(f"gnomAD shape: {gnomad.shape}")
"""

CODE_TASK4C = """\
import pandas as pd
import numpy as np
from pathlib import Path

results_dir = Path("../results")
data_dir    = Path("../data")

TARGET_POS = 1172

gnomad = pd.read_csv(data_dir / "gnomad_otof_variants.csv")

# Try to find the variant in gnomAD
# Check which column holds the position info
pos_col = None
for c in gnomad.columns:
    if "pos" in c.lower() or "position" in c.lower():
        pos_col = c
        break

name_col = None
for c in gnomad.columns:
    if "name" in c.lower() or "hgvs" in c.lower() or "variant" in c.lower() or "consequence" in c.lower():
        name_col = c
        break

print(f"Position column: {pos_col}")
print(f"Name/HGVS column: {name_col}")

# Match by position or by name substring
match_mask = pd.Series([False] * len(gnomad), index=gnomad.index)

if pos_col:
    try:
        match_mask |= (gnomad[pos_col].astype(str).str.extract(r'(\\d+)')[0].astype(float) == TARGET_POS)
    except Exception:
        pass

for c in gnomad.columns:
    try:
        match_mask |= gnomad[c].astype(str).str.contains("1172", na=False)
    except Exception:
        pass

var_gnomad = gnomad[match_mask]
print(f"\\ngnomAD rows matching position 1172: {len(var_gnomad)}")
if len(var_gnomad) > 0:
    print(var_gnomad.to_string(index=False))
else:
    print("Not found in gnomAD extract. Showing columns for reference:")
    print(gnomad.dtypes)
    print(gnomad.head(3).to_string())

# Check for ancestry-specific AF columns
anc_cols = [c for c in gnomad.columns if any(
    pop in c.lower() for pop in ["afr", "amr", "eas", "fin", "nfe", "sas", "mid", "asj", "oth", "ancestry", "pop"]
)]
print(f"\\nAncestry-specific columns: {anc_cols}")
"""

CODE_TASK4D = """\
import pandas as pd
import numpy as np
from pathlib import Path

results_dir = Path("../results")
data_dir    = Path("../data")

TARGET_POS = 1172

# ------------------------------------------------------------------
# 3. ConSurf conservation score for residue 1172
# ------------------------------------------------------------------
consurf = pd.read_csv(data_dir / "consurf_Q9HC10_grades.csv")
print(f"ConSurf columns: {list(consurf.columns)}")
print(f"Shape: {consurf.shape}")

# Find the position column
pos_col_cs = None
for c in consurf.columns:
    if "pos" in c.lower() or "resid" in c.lower() or "seq" in c.lower() or "aa" in c.lower():
        pos_col_cs = c
        print(f"  Candidate position col: {c}")

# Show a few rows to understand structure
print(consurf.head(3).to_string())
"""

CODE_TASK4E = """\
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

results_dir = Path("../results")
data_dir    = Path("../data")

TARGET_POS  = 1172
VARIANT_LBL = "p.Arg1172Gln"

consurf = pd.read_csv(data_dir / "consurf_Q9HC10_grades.csv")

# Detect position and score columns
pos_col_cs = consurf.columns[0]  # usually first column
for c in consurf.columns:
    try:
        vals = pd.to_numeric(consurf[c], errors="coerce")
        if vals.between(1, 9997).sum() > 100:
            pos_col_cs = c
            break
    except Exception:
        pass

score_col = None
for c in consurf.columns:
    if "grade" in c.lower() or "score" in c.lower() or "conserv" in c.lower():
        score_col = c
        break
if score_col is None:
    # pick numeric column that looks like 1-9 conservation scores
    for c in consurf.columns:
        try:
            vals = pd.to_numeric(consurf[c], errors="coerce").dropna()
            if vals.between(1, 9).mean() > 0.7:
                score_col = c
                break
        except Exception:
            pass

print(f"Position column: {pos_col_cs}")
print(f"Score column: {score_col}")

consurf["_pos"]   = pd.to_numeric(consurf[pos_col_cs], errors="coerce")
consurf["_score"] = pd.to_numeric(consurf[score_col], errors="coerce") if score_col else np.nan

row_1172 = consurf[consurf["_pos"] == TARGET_POS]
if len(row_1172) == 0:
    print(f"\\nResidue {TARGET_POS} not found in ConSurf data.")
    score_1172 = np.nan
else:
    score_1172 = row_1172["_score"].values[0]
    print(f"\\nConSurf grade for residue {TARGET_POS}: {score_1172}")
    print(f"  (1 = most variable, 9 = most conserved)")

all_scores = consurf["_score"].dropna()
pct = (all_scores <= score_1172).mean() * 100 if not np.isnan(score_1172) else np.nan
print(f"  Percentile among all OTOF residues: {pct:.1f}% ({int(score_1172)}/9 conservation grade)")

# Visualise: distribution of all ConSurf grades + position of 1172
fig, ax = plt.subplots(figsize=(8, 4))
bins = np.arange(0.5, 10.5, 1)
ax.hist(all_scores, bins=bins, color="#4477AA", edgecolor="black", linewidth=0.5,
        alpha=0.8, label="All OTOF residues")
if not np.isnan(score_1172):
    ax.axvline(score_1172, color="#CC3311", linewidth=2.5,
               label=f"{VARIANT_LBL} (grade {int(score_1172)})")
ax.set_xlabel("ConSurf conservation grade (1=variable, 9=conserved)", fontsize=11)
ax.set_ylabel("Number of residues", fontsize=11)
ax.set_title(f"ConSurf conservation landscape\\n{VARIANT_LBL} highlighted", fontsize=12, fontweight="bold")
ax.legend(fontsize=10)
ax.spines[["top", "right"]].set_visible(False)
ax.set_xticks(range(1, 10))
plt.tight_layout()
out_png = results_dir / "arg1172gln_consurf.png"
plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.show()
print(f"Saved: {out_png}")

# Summary table
summary_1172 = pd.DataFrame([{
    "Variant": VARIANT_LBL,
    "Position": TARGET_POS,
    "Domain": "C2E (aa 1158-1273)",
    "ConSurf_grade": int(score_1172) if not np.isnan(score_1172) else "N/A",
    "Percentile_conserved": f"{pct:.1f}%" if not np.isnan(pct) else "N/A",
}])
out_csv = results_dir / "arg1172gln_summary.csv"
summary_1172.to_csv(out_csv, index=False)
print(f"Saved: {out_csv}")
print(summary_1172.to_string(index=False))
"""

# ===========================================================================
# Assemble and save
# ===========================================================================

new_cells = [
    md_cell(MD_TASK1),
    code_cell(CODE_TASK1),
    md_cell(MD_TASK2),
    code_cell(CODE_TASK2A),
    code_cell(CODE_TASK2B),
    md_cell(MD_TASK3),
    code_cell(CODE_TASK3A),
    code_cell(CODE_TASK3B),
    md_cell(MD_TASK4),
    code_cell(CODE_TASK4A),
    code_cell(CODE_TASK4B),
    code_cell(CODE_TASK4C),
    code_cell(CODE_TASK4D),
    code_cell(CODE_TASK4E),
    md_cell(MD_TASK5),
    code_cell(CODE_TASK5A),
    code_cell(CODE_TASK5B),
    code_cell(CODE_TASK5C),
    code_cell(CODE_TASK5D),
]

nb["cells"].extend(new_cells)

with open(NB_PATH, "w") as f:
    json.dump(nb, f, indent=1)

print(f"Appended {len(new_cells)} cells. Total cells: {len(nb['cells'])}")
