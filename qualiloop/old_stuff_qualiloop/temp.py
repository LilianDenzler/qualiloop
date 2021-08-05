#!/usr/bin/env python
import sys
import os
import numpy as np
import csv
import pandas as pd

df=pd.read_csv("~/sync_project/WWW/CDRH3loop/CDRH3loop/qualiloop/full_dataset_uniform_08_04.csv")
cols=[col for col in df.columns if "nom" not in col and "bin" not in col]
df=df[cols]
df.to_csv(os.path.join("~/sync_project/WWW/CDRH3loop/CDRH3loop/feature_dir","full_features.csv"))