TAB=[["Line", "MP", "eMP", "BS", "eBS", "BG", "eBG"], ["""[S II]
6731""", 0.078, 0.037, 0, 0, 0.93, 0.104], ["""[S III]
9069""", 0.429, 0.045, 0.201, 0.058, 0.29, 0.038], ["""[Ar III]
7136""", 0.209, 0.034, 0.308, 0.036, 0.21, 0.025], ["""H I
4861""", 0.135, 0.039, 0.363, 0.026, 0.217, 0.027], ["""[O III]
5007""", 0.07, 0.023, 0.36, 0.016, 0.117, 0.017]]
from astropy.table import Table
from matplotlib import pyplot as plt
import seaborn as sns
sns.set_context("talk")
tab = Table(rows=TAB[1:], names=TAB[0])
figfile = "emline-excess-plot.pdf"
fig, ax = plt.subplots()
line, = ax.plot(tab["Line"], tab["MP"])
ax.errorbar(tab["Line"], tab["MP"], yerr=tab["eMP"], fmt="o", color=line.get_color())
line, = ax.plot(tab["Line"], tab["BS"])
ax.errorbar(tab["Line"], tab["BS"], yerr=tab["eBS"], fmt="D", color=line.get_color())
line, = ax.plot(tab["Line"], tab["BG"])
ax.errorbar(tab["Line"], tab["BG"], yerr=tab["eBG"], fmt="s", color=line.get_color())
ax.set(ylabel="Brightness / Background")
sns.despine()
fig.tight_layout()
fig.savefig(figfile)
print(figfile, end="")
