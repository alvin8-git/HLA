"""
14_pipeline_flowchart.py
Generates a methods pipeline flowchart for the HLA registry size report.
Output: figures/pipeline_flowchart.png
"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

HERE    = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(HERE, 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

STEPS = [
    ('Raw HLA typing data\n(BMDP 53,297 + SCBB 5,889 + HSA 564\ndonors & cord blood units)',
     '#D6EAF8'),
    ('EM haplotype phasing\n(5-locus: A·B·C·DRB1·DQB1)\nOutput: haplotype frequencies {f̂ₖ}',
     '#D5F5E3'),
    ('Hardy–Weinberg expansion\nDiplotype frequencies:\nF(hᵢ,hᵢ) = fᵢ²;  F(hᵢ,hⱼ) = 2fᵢfⱼ',
     '#D5F5E3'),
    ('Coverage function C(N)\nC(N) = Σₖ Fₖ·[1−(1−Fₖ)ᴺ]\n(Beatty et al. 1995 framework)',
     '#D5F5E3'),
    ('Binary search for N*\nMinimum N satisfying C(N) ≥ θ\n(θ = 75%, 85%, 90%, 95%)',
     '#FCF3CF'),
    ('Dirichlet bootstrap (B=1,000)\n95% CI on N*\nConcentration: nₑff × f̂ₖ',
     '#FADBD8'),
]

SIDE_NOTES = {
    1: 'GENE[RATE] validation\n(max discrepancy 0.27%)',
    3: 'Same-ethnicity &\ncross-ethnic variants',
    5: 'Bias-corrected via\nbootstrap median',
}

fig, ax = plt.subplots(figsize=(9, 7), facecolor='white')
ax.set_xlim(0, 10)
ax.set_ylim(0, len(STEPS) + 0.5)
ax.axis('off')

box_w, box_h = 6.5, 0.72
box_x = 1.75

for i, (label, color) in enumerate(STEPS):
    y = len(STEPS) - 1 - i
    box = FancyBboxPatch((box_x, y + 0.05), box_w, box_h,
                         boxstyle='round,pad=0.04',
                         facecolor=color, edgecolor='#555555', linewidth=1.0)
    ax.add_patch(box)
    ax.text(box_x + box_w / 2, y + box_h / 2 + 0.05, label,
            ha='center', va='center', fontsize=8.5, multialignment='center',
            fontfamily='monospace')

    # Arrow between steps
    if i < len(STEPS) - 1:
        ax.annotate('', xy=(box_x + box_w / 2, y + 0.05),
                    xytext=(box_x + box_w / 2, y + 0.9),
                    arrowprops=dict(arrowstyle='->', color='#333333', lw=1.5))

    # Side notes
    if i in SIDE_NOTES:
        note_x = box_x + box_w + 0.25
        note_y = y + box_h / 2 + 0.05
        ax.annotate(SIDE_NOTES[i],
                    xy=(box_x + box_w, note_y),
                    xytext=(note_x + 0.1, note_y),
                    fontsize=7.5, color='#555555', va='center',
                    arrowprops=dict(arrowstyle='-', color='#999999', lw=0.8))

ax.set_title('Analysis pipeline', fontsize=11, fontweight='bold', pad=6)

plt.tight_layout()
out = os.path.join(FIG_DIR, 'pipeline_flowchart.png')
fig.savefig(out, dpi=150, bbox_inches='tight')
plt.close(fig)
print(f'Saved: {out}')


if __name__ == '__main__':
    pass
