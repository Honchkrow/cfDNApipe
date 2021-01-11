from cfDNApipe import *
import glob
import pysam
import numpy
import pickle
import matplotlib.pyplot as plt

pipeConfigure(
    threads=20,
    genome="hg19",
    refdir=r"/home/wzhang/genome/hg19",
    outdir=r"/data/wzhang/pcs_final/pcs_fraglen",
    data="WGS",
    type="paired",
    JavaMem="10G",
    build=True,
)

verbose = False

case_bed = glob.glob("/data/wzhang/pcs_final/HCC/*.bed")
ctrl_bed = glob.glob("/data/wzhang/pcs_final/Healthy/*.bed")

res_fraglenplot_comp = fraglenplot_comp(
    casebedInput=case_bed, ctrlbedInput=ctrl_bed, caseupstream=True, verbose=verbose
)

caseInput = res_fraglenplot_comp.getOutput("casepickleOutput")
ctrlInput = res_fraglenplot_comp.getOutput("ctrlpickleOutput")

labelInput = ["case", "control"]

caseprop = []
ctrlprop = []
fig = plt.figure(figsize=(10, 8))
for i in range(len(caseInput)):
    with open(caseInput[i], "rb") as f:
        dataInput = pickle.load(f)
        dataInput = dict(sorted(dataInput.items()))
        keys = np.fromiter(dataInput.keys(), dtype=int)
        vals = np.fromiter(dataInput.values(), dtype=int)
        idx = np.where(keys >= 250)
        vals = vals / np.sum(vals)
        caseprop.append(np.sum(vals[idx]))
    (p1,) = plt.plot(keys, vals, c="r", linewidth=0.5)
for i in range(len(ctrlInput)):
    with open(ctrlInput[i], "rb") as f:
        dataInput = pickle.load(f)
        dataInput = dict(sorted(dataInput.items()))
        keys = np.fromiter(dataInput.keys(), dtype=int)
        vals = np.fromiter(dataInput.values(), dtype=int)
        idx = np.where(keys >= 250)
        vals = vals / np.sum(vals)
        ctrlprop.append(np.sum(vals[idx]))
    (p2,) = plt.plot(keys, vals, c="b", linewidth=0.5)

plt.tick_params(labelsize=15)
font = {
    "family": "Times New Roman",
    "weight": "normal",
    "size": 20,
}
plt.xlabel("DNA Fragment Size (base pair)", font)
plt.ylabel("Density", font)

font_legend = {
    "family": "Times New Roman",
    "weight": "normal",
    "size": 15,
}
plt.legend([p1, p2], labelInput, loc="best", prop=font_legend)
plt.savefig("111.png")
plt.close(fig)

fig = plt.figure(figsize=(10, 8))
casegory = [labelInput[0] for i in range(len(caseprop))]
ctrlgory = [labelInput[1] for i in range(len(ctrlprop))]
propdf = pd.DataFrame({"category": casegory + ctrlgory, "proportion": caseprop + ctrlprop})
bp = sns.violinplot(x="category", y="proportion", data=propdf)
bp.set_xlabel("", fontsize=20)
bp.set_ylabel("Proportion of fragments below 150bp", fontsize=20)
bp.tick_params(labelsize=15)
y, h = propdf["proportion"].max() + 0.1, 0.02
t, p = stats.ttest_ind(caseprop, ctrlprop, equal_var=False)
if p >= 0.05:
    text = "p = " + str(p)
elif p >= 0.01:
    text = "p < 0.05"
elif p >= 0.001:
    text = "p < 0.01"
elif p >= 0.0001:
    text = "p < 0.001"
elif p >= 0:
    text = "p < 0.0001"
plt.plot([0, 0, 1, 1], [y, y + h, y + h, y], lw=1, c="k")
plt.text(0.5, y + h, text, ha="center", va="bottom", color="k", fontdict=font_legend)
plt.savefig("222.png")
# plt.savefig(os.path.splitext(plotOutput[1])[0] + ".pdf")
plt.close(fig)
