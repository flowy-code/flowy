import numpy as np
import pandas as pd

dataframe = pd.read_csv(
    "/home/moritz/lavaflow/flowtastic/examples/ETNA_LFS1/tinit_33.asc",
    skiprows=6,
    delimiter=" ",
    header=None
)

dataframe.replace(float("nan"), -9999)
data = dataframe.to_numpy()

np.savetxt("tinit_33_n.asc", data)