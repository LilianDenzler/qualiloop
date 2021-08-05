#!/usr/bin/env python
import sys
import os
import pandas as pd
import numpy as np
import seaborn as sns
from yellowbrick.features import JointPlotVisualizer
import matplotlib.pyplot as plt

#################################################################
#GRAPH: Sequence Similarity/Identity vs RMSD
#DO FITTING!!!!!!!!!!!!!!!!
#################################################################
def plot_graph(data, feature):
	data = pd.read_csv(data)
	data.dropna(inplace=True)
	feature_list=[col for col in data.columns if "bin" not in col and "nom" not in col and "local" not in col and "global" not in col]
	plot=JointPlotVisualizer(columns=feature, size=(900, 700))

	X=data[feature_list]
	y=data["local_CA"].astype(float)
	plot.fit_transform(X, y)
	plot.finalize()
	plt.ylabel('RMSD LocalCA (Å)')
	plt.xlabel(feature)
	plt.title("Correlation of "+feature+" and RMSD LocalCA (Å)")


	plt.show()
	#fig.savefig("seqId_similarity_RMSD.png") 