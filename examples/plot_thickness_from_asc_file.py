from typing import Optional
import typer
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from linecache import getline
import pyvista as pv
from dataclasses import dataclass

app = typer.Typer()


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
        header = [getline(str(path_to_file), i) for i in range(1, 7)]
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

    def max_height(self):
        return np.max(self.height_data)

    def min_height(self):
        return np.min(self.height_data[self.height_data > self.nd + 1])

    def filter_height_data(self, replace=float("nan")):
        mask = self.height_data > self.nd + 1
        self.height_data[~mask] = replace


@app.command()
def plot_pyvista(
    path_asc_file_initial: Path = typer.Argument(help="Path to initial ASC file"),
    path_asc_file_final: Path = typer.Argument(
        help="Path to final ASC file after lava emplacement"
    ),
    warp_factor: float = typer.Option(
        1.0,
        "--warp",
        min=0.0,
        max=1.0,
        help="Factor by which to warp the terrain height.",
    ),
    threshold: float = typer.Option(
        0.01, "--thresh", help="Any elevation smaller than the threshold is not plotted"
    ),
    image_outfile: Path = typer.Option(
        None,
        "-o",
        help="Path to the output image file. If not set, an image called image.png will be created.",
    ),
    interactive: bool = typer.Option(False, "-i", help="Opens an interactive window"),
):
    asc_file_initial = AscFile(path_asc_file_initial)
    asc_file_final = AscFile(path_asc_file_final)

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
    grid_initial = grid.warp_by_scalar("height_initial", warp_factor)
    # For the lava flow
    grid_final = grid.warp_by_scalar("height_final", warp_factor).threshold(
        threshold, scalars="flow"
    )

    # For the terrain
    z_cells = np.linspace(0, 10, 10)
    xx = np.repeat(grid_initial.x, len(z_cells), axis=-1)
    yy = np.repeat(grid_initial.y, len(z_cells), axis=-1)
    zz = np.repeat(grid_initial.z, len(z_cells), axis=-1) - np.cumsum(z_cells).reshape(
        (1, 1, -1)
    )
    mesh = pv.StructuredGrid(xx, yy, zz)
    mesh["Elevation"] = zz.ravel(order="F")
    mesh = mesh.threshold(asc_file_final.nd + 1, scalars="Elevation")
    # opacity_terrain = [0, 0.2, 0.9, 0.6, 0.3]
    # opacity_terrain = "sigmoid_r"
    opacity_terrain = None

    if interactive:
        p = pv.Plotter(off_screen=False)
    else:
        p = pv.Plotter(off_screen=True)
    p.add_mesh(
        mesh,
        scalars="Elevation",
        color=True,
        smooth_shading=True,
        opacity=opacity_terrain,
        cmap="gray",
    )
    p.add_mesh(grid_final, color="red", opacity=0.8, smooth_shading=True)

    if interactive:
        p.show()
    if image_outfile is not None:
        p.screenshot(image_outfile)
    else:
        p.screenshot("./image.png")


@app.command()
def plot_pyplot(
    path_asc_file_initial: Path = typer.Argument(help="Path to initial ASC file"),
    path_asc_file_final: Path = typer.Argument(
        help="Path to final ASC file after lava emplacement"
    ),
    image_outfile: Path = typer.Option(
        None,
        "-o",
        help="Path to the output image file. If not set, an image called image.png will be created.",
    ),
    interactive: bool = typer.Option(False, "-i", help="Opens an interactive window"),
):
    asc_file_initial = AscFile(path_asc_file_initial)
    asc_file_final = AscFile(path_asc_file_final)

    asc_file_initial.filter_height_data()
    asc_file_final.filter_height_data()

    x_data = asc_file_initial.x_data()
    y_data = asc_file_initial.y_data()
    min_height = asc_file_initial.min_height()

    plt.contourf(
        x_data,
        y_data,
        asc_file_initial.height_data - min_height,
        levels=20,
        cmap="gray",
    )

    lava_thickness = asc_file_final.height_data - asc_file_initial.height_data
    lava_thickness[lava_thickness < 1e-10] = float("nan")
    plt.contourf(x_data, y_data, lava_thickness, cmap="magma")
    if interactive:
        plt.show()
    if image_outfile is not None:
        plt.savefig(image_outfile, dpi=300)
    else:
        plt.savefig("./image.png")


if __name__ == "__main__":
    app()
