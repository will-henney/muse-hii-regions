TAB=[["", "I([C I] 8727)", "I(Max DRL)", "I(Sum DRL)", "Number DRL", "Max DRL"], ["Orion South", 0.17, 0.014, 0.026, 3, 9114], ["Orion Bar", 0.06, 0.009, 0.017, 3, 9114], ["M 17 A", 0.27, 0.019, 0.16, 15, 8459], ["M 17 B", 0.73, 0.08, 0.63, 16, 8894], ["30 Dor Globule", 0.12, 0.022, 0.16, 18, 9029], ["30 Dor Clump", 0.06, 0.086, 0.71, 20, 9029], ["NGC 346 MYSO", 0.24, 0.2, 5.93, 80, 9297], ["NGC 346 Zone 0", 0.21, 0.399, 14.83, 114, 9029]]
import numpy as np
from astropy.table import Table
from matplotlib import pyplot as plt
import matplotlib
import seaborn as sns
sns.set_context("talk")
sns.set_palette("icefire", n_colors=6)
tab = Table(rows=TAB[1:], names=TAB[0])
figfile = "drl-regions-comparison.pdf"
fig, ax = plt.subplots(figsize=(6, 6))
xxx = tab["I([C I] 8727)"]
yyy = tab["I(Max DRL)"]
zzz = tab["Number DRL"]
xx = np.vstack([xxx[0::2], xxx[1::2]])
yy = np.vstack([yyy[0::2], yyy[1::2]])
zz = np.vstack([zzz[0::2], zzz[1::2]])
icols = [0, 1, -2, -1]
for i, x, y, z in zip(icols, xx.T, yy.T, zz.T):
    color = sns.color_palette()[i]
    line, = ax.plot(x, y, linewidth=3, zorder=-100, color=color)
    ax.scatter(x[0], y[0], s=20 + z[0]*10, color=color)
    ax.scatter(x[1], y[1], s=20 + z[1]*10, color=color, facecolor="white", linewidths=3)
ax.set(
    xlabel=r"100 $\times$  [C I] λ8727 / Hβ",
    ylabel=r"100 $\times$  max( DRL ) / Hβ",
    xlim=[0.03, 2.0],
    ylim=[0.003, 0.9],
    xscale="log",
    yscale="log",
)
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_xticks([0.1, 0.3, 1.0])
ax.set_yticks([0.01, 0.03, 0.1, 0.3])
sns.despine()
fig.tight_layout()
fig.savefig(figfile, bbox_inches="tight")
print(figfile, end="")
