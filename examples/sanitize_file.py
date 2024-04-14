import numpy as np
import pandas as pd
from linecache import getline
import argparse


def sanitize(path_to_file, path_to_sanitized_file):
    header = [getline(str(path_to_file), i) for i in range(1, 7)]
    values = [float(h.split(" ")[-1].strip()) for h in header]
    (
        cols,
        rows,
        lx,
        ly,
        cell,
        nd,
    ) = values

    cols = int(cols)
    rows = int(rows)

    dataframe = pd.read_csv(path_to_file, skiprows=6, delimiter=" ", header=None)

    dataframe.replace(float("nan"), -9999, inplace=True)
    data = dataframe.to_numpy()

    # Recompute number of rows and columns
    rows_data = data.shape[0]
    cols_data = data.shape[1]

    if cols != cols_data:
        print(f"Warning: `cols` from header is {cols}, but data has {cols_data} cols")

    if rows != rows_data:
        print(f"Warning: `rows` from header is {rows}, but data has {rows_data} rows")

    with open(path_to_sanitized_file, "w") as f:
        f.write(f"ncols {cols_data}\n")
        f.write(f"nrows {rows_data}\n")
        f.write(f"xllcorner {lx}\n")
        f.write(f"yllcorner {ly}\n")
        f.write(f"cellsize {cell}\n")
        f.write(f"NODATA_value {nd}\n")
        np.savetxt(f, data, fmt="%0.4f")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        prog="sanitize_asc_file",
        description="Sanitizes an .asc file",
    )

    parser.add_argument("src", help="the input .asc file (the file to be sanitized)")
    parser.add_argument("dest", help="the output .asc file (the sanitized file)")

    args = parser.parse_args()

    sanitize(args.src, args.dest)
