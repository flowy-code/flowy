import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from linecache import getline
import pyvista as pv
from dataclasses import dataclass


@dataclass
class AscFile:
    cols: int = 0
    rows: int = 0
    lx: float = 0
    ly: float = 0
    cell: float = 0
    nd: float = -9999
    height_data = np.array([], dtype=float)

    def __init__(self, path_to_file: Path):
        header = [getline(path_to_file, i) for i in range(1, 7)]
        values = [float(h.split(" ")[-1].strip()) for h in header]
        (
            self.cols,
            self.rows,
            self.lx,
            self.ly,
            self.cell,
            self.nd,
        ) = values
        self.height_data = np.loadtxt(path_to_file, skiprows=6, delimiter=" ")

    def x_data(self):
        return np.arange(
            self.lx,
            self.lx + self.height_data.shape[0] * self.cell,
            self.cell,
        )

    def y_data(self):
        return np.arange(
            self.lx,
            self.lx + self.height_data.shape[0] * self.cell,
            self.cell,
        )

    def min_height(self):
        return np.min(self.height_data[self.height_data > self.nd + 1])

    def filter_height_data(self, replace=float("nan")):
        mask = self.height_data > self.nd + 1
        self.height_data[~mask] = replace


def plot_height_data_pyvista(
    path_asc_file_initial: Path,
    path_asc_file_final: Path,
    warp_factor=1.0,
    threshold=0.01,
):
    asc_file_initial = AscFile(path_asc_file_initial)
    asc_file_final = AscFile(path_asc_file_final)

    asc_file_initial.filter_height_data()
    asc_file_final.filter_height_data()

    x_data = asc_file_initial.x_data()
    y_data = asc_file_initial.y_data()
    min_height = asc_file_initial.min_height()

    grid = pv.ImageData()
    grid.dimensions = [len(x_data), len(y_data), 1]

    grid.point_data["height_initial"] = (
        asc_file_initial.height_data.flatten(order="F") - min_height
    )
    grid.point_data["height_final"] = (
        asc_file_final.height_data.flatten(order="F") - min_height
    )
    grid.point_data["flow"] = (
        grid.point_data["height_final"] - grid.point_data["height_initial"]
    )

    grid = grid.threshold(asc_file_final.nd + 1, scalars="height_initial")

    grid_initial = grid.warp_by_scalar("height_initial", warp_factor)
    grid_final = grid.warp_by_scalar("height_final", warp_factor).threshold(
        threshold, scalars="flow"
    )

    p = pv.Plotter()
    p.add_mesh(grid_initial, scalars="height_initial", color=True, smooth_shading=True)
    p.add_mesh(grid_final, color="red", opacity=0.8, smooth_shading=True)
    p.show()


def plot_height_data_pyplot(path_asc_file_initial: Path, path_asc_file_final: Path):
    asc_file_initial = AscFile(path_asc_file_initial)
    asc_file_final = AscFile(path_asc_file_final)

    asc_file_initial.filter_height_data()
    asc_file_final.filter_height_data()

    x_data = asc_file_initial.x_data()
    y_data = asc_file_initial.y_data()
    min_height = asc_file_initial.min_height()

    plt.contourf(x_data, y_data, asc_file_initial.height_data - min_height)
    plt.contour(
        x_data, y_data, asc_file_final.height_data - asc_file_initial.height_data
    )

    plt.show()


if __name__ == "__main__":
    path_asc_file_initial = "/home/moritz/lavaflow/flowy/output_hawai/initial.asc"
    path_asc_file_final = "/home/moritz/lavaflow/flowy/output_hawai/output.asc"
    plot_height_data_pyvista(path_asc_file_initial, path_asc_file_final)
    plot_height_data_pyplot(path_asc_file_initial, path_asc_file_final)
