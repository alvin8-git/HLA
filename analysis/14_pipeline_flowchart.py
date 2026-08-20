"""
14_pipeline_flowchart.py
Generates the methods pipeline flowchart for the HLA registry size report.
Output: figures/pipeline_flowchart.png
"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

HERE    = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(HERE, 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

# --- Layout constants ---------------------------------------------------------
BOX_W  = 7.0    # box width in data units
BOX_H  = 0.70   # box height in data units
GAP    = 0.40   # gap between boxes
STEP   = BOX_H + GAP   # y distance per step
BOX_X  = 0.5   # left edge of boxes
CX     = BOX_X + BOX_W / 2   # horizontal centre

STEPS = [
    ('Raw HLA typing data\n'
     'BMDP 53,297 + SCBB 5,889 donors + HSA 564 patient-donor pairs\n'
     'Five loci typed: HLA-A, -B, -C, -DRB1, -DQB1',
     '#D6EAF8', None),
    ('EM haplotype phasing (in-house Python algorithm)\n'
     'Input: unphased diplotypes, ALL five-locus-typed individuals (no cap)\n'
     'Retain haplotypes at f ≥ 1e-6, below every group 1/(2n)  →  100% of mass',
     '#D5F5E3',
     'Validated vs GENE[RATE]\n(1,227–3,011 shared\nhaplotypes, r = 0.75–0.95)'),
    ('Hardy–Weinberg expansion to diplotype frequencies\n'
     'F(hᵢ, hᵢ) = f̂ᵢ²     [homozygous]\n'
     'F(hᵢ, hⱼ) = 2·f̂ᵢ·f̂ⱼ   [heterozygous, i ≠ j]',
     '#D5F5E3', None),
    ('Coverage function  C(N) = Σₖ Fₖ · [1 − (1 − Fₖ)ᴺ]\n'
     'Expected fraction of patients finding ≥1 matched donor in N draws\n'
     'Computed for same-ethnicity and cross-ethnic donor pools',
     '#D5F5E3',
     'Same-ethnicity &\ncross-ethnic variants'),
    ('Binary search for N*  (50 iterations, range 10³–10¹⁰ donors)\n'
     'Minimum N such that C(N) ≥ θ\n'
     'θ ∈ {75%, 85%, 90%, 95%} coverage targets',
     '#FCF3CF', None),
    ('Dirichlet parametric bootstrap  (B = 1,000 resamples)\n'
     'Resample {f̂ₖ} with concentration nₑff × f̂ₖ per ethnicity\n'
     'Report: bootstrap median N* ± 95% CI (2.5th–97.5th percentile)',
     '#FADBD8',
     'Bias-corrected via\nbootstrap median'),
]

N = len(STEPS)
fig_h = N * STEP + GAP + 0.6   # total height
fig, ax = plt.subplots(figsize=(9, fig_h), facecolor='white')
ax.set_facecolor('white')
ax.set_xlim(0, BOX_X + BOX_W + 2.8)
ax.set_ylim(-GAP / 2, N * STEP + 0.5)
ax.axis('off')

for i, (label, color, note) in enumerate(STEPS):
    # y coordinate of box bottom (step 0 at top)
    y_bot = (N - 1 - i) * STEP

    # Draw box
    box = FancyBboxPatch((BOX_X, y_bot), BOX_W, BOX_H,
                         boxstyle='round,pad=0.06',
                         facecolor=color, edgecolor='#444444', linewidth=1.2)
    ax.add_patch(box)

    # Label inside box
    ax.text(CX, y_bot + BOX_H / 2, label,
            ha='center', va='center', fontsize=8,
            multialignment='center', linespacing=1.45)

    # Arrow FROM bottom of this box TO top of next box
    if i < N - 1:
        y_arrow_top    = y_bot                         # bottom of current box
        y_arrow_bottom = y_bot - GAP                   # top of next box = y_bot - GAP
        ax.annotate(
            '',
            xy=(CX, y_arrow_bottom),           # arrowhead at top of next box
            xytext=(CX, y_arrow_top),           # tail at bottom of current box
            arrowprops=dict(
                arrowstyle='->', color='#222222',
                lw=1.6, mutation_scale=14,
            ),
        )

    # Side note
    if note:
        note_x  = BOX_X + BOX_W + 0.15
        note_y  = y_bot + BOX_H / 2
        ax.annotate(
            note,
            xy=(BOX_X + BOX_W, note_y),
            xytext=(note_x + 0.08, note_y),
            fontsize=7.5, color='#444444', va='center',
            arrowprops=dict(arrowstyle='-', color='#888888', lw=0.9),
        )

ax.set_title('Figure 1. Analysis pipeline — from raw HLA typing data to '
             'registry size estimates with bootstrap CIs',
             fontsize=10, fontweight='bold', pad=8, loc='left')

plt.tight_layout(pad=0.5)
out = os.path.join(FIG_DIR, 'pipeline_flowchart.png')
fig.savefig(out, dpi=150, bbox_inches='tight')
plt.close(fig)
print(f'Saved: {out}')
