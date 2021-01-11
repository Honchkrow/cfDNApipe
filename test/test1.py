from cfDNApipe import *
import glob
import pysam
import numpy
import pickle
import matplotlib.pyplot as plt


caseInput = glob.glob("/opt/tsinghua/bettyzhong/novo1101-5G/HBV_HCC/intermediate_result/step_08_fraglenplot/*.pickle")


caseprop = []


for i in range(len(caseInput)):
    with open(caseInput[i], "rb") as f:
        dataInput = pickle.load(f)
        dataInput = dict(sorted(dataInput.items()))
        keys = np.fromiter(dataInput.keys(), dtype=int)
        vals = np.fromiter(dataInput.values(), dtype=int)
        idx = np.where(keys >= 250)
        vals = vals / np.sum(vals)
        print(caseInput[i])
        print(np.sum(vals[idx]))