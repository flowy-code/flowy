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
        self.height_data = np.transpose(np.flipud(self.height_data))

    def x_data(self):
        return np.linspace(
            self.lx,
            self.lx + self.height_data.shape[0] * self.cell,
            self.height_data.shape[0],
        )

    def y_data(self):
        return np.linspace(
            self.ly,
            self.ly + self.height_data.shape[1] * self.cell,
            self.height_data.shape[1],
        )

    def max_height(self):
        return np.max(self.height_data)

    def min_height(self):
        return np.min(self.height_data[self.height_data > self.nd + 1])

    def filter_height_data(self, replace=float("nan")):
        mask = self.height_data > self.nd
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
        help="Factor by which to warp the terrain height.",
    ),
    threshold: float = typer.Option(
        1e-12,
        "--thresh",
        help="Any flow thickness smaller than the threshold is not plotted",
    ),
    image_outfile: Path = typer.Option(
        None,
        "-o",
        help="Path to the output image file. If not set, an image called image.png will be created.",
    ),
    interactive: bool = typer.Option(False, "-i", help="Opens an interactive window"),
    contour_spacing: int = typer.Option(
        -1, "-c", help="Contour spacing in meters. Negative to disable contouring"
    ),
    contour_color: str = typer.Option("white", "--ccolor", help="Contour color"),
):
    # We use this offset to plot the final grid slightly above the initial grid
    # This gets rid of artifacts due to rounding errors
    height_offset = 0.1

    asc_file_initial = AscFile(path_asc_file_initial)
    asc_file_final = AscFile(path_asc_file_final)

    x_data = asc_file_initial.x_data()
    y_data = asc_file_initial.y_data()

    # ====================================================================
    # Build mesh for the terrain
    # ====================================================================
    thickness_terrain = 500  # gives the terrain a finite thickness [unit is meter]
    # Warp the heights
    warped_initial_height_data = (
        warp_factor * (asc_file_initial.height_data - asc_file_initial.min_height())
        + asc_file_initial.min_height()
    )

    # Here we add a second layer, to give the plot a thickness in the z direction
    x, y, z = np.meshgrid(x_data, y_data, np.ones(2), indexing="ij")
    z[:, :, 0] = warped_initial_height_data
    z[:, :, 1] = warped_initial_height_data - thickness_terrain
    # The shape of x,y,z is [len(x), leny(y), 2]

    # We compute the elevation without the thickness in z,
    # becaus we want to avoid distorting the colormap.
    # Elevation needs to have the same shape as x,y,z
    # therefore, we copy the height data twice.
    # Also note that we are not using the warped height data here
    elevation = np.zeros_like(z)
    elevation[:, :, 0] = asc_file_initial.height_data
    elevation[:, :, 1] = asc_file_initial.height_data

    grid_initial = pv.StructuredGrid(x, y, z)
    grid_initial["elevation"] = elevation.flatten(order="F")

    # ====================================================================
    # Build mesh for the flow
    # ====================================================================
    # When computing the flow thickness, we dont use the warped heights
    flow_thickness = (asc_file_final.height_data).flatten(order="F")

    warped_final_height_data = (
        warp_factor * (asc_file_final.height_data) + warped_initial_height_data
    )

    # This time we dont need a finite thickness
    # Therefore, we just use x_data and y_data
    x, y = np.meshgrid(x_data, y_data, indexing="ij")
    z = warped_final_height_data + height_offset
    # The shape of x,y,z is [len(x_data), len(y_data)]

    grid_surface = pv.StructuredGrid(x, y, z)
    grid_surface.point_data["flow_thickness"] = flow_thickness
    grid_surface.point_data["elevation"] = asc_file_initial.height_data.flatten(
        order="F"
    )

    grid_final = grid_surface.threshold(threshold, scalars="flow_thickness")

    if interactive:
        p = pv.Plotter(off_screen=False)
    else:
        p = pv.Plotter(off_screen=True, window_size=[1024, 768], image_scale=4)

    if contour_spacing > 0:
        min_height = asc_file_initial.min_height()
        max_height = asc_file_initial.max_height()
        height_contours = grid_surface.contour(
            np.arange(min_height, max_height, contour_spacing), scalars="elevation"
        )

        p.add_mesh(height_contours, color=contour_color, line_width=2)

    light_args = dict(
        line_width=3,
        show_edges=False,
        metallic=False,
        specular=0.05,
        ambient=0.3,
        pickable=False,
    )

    p.add_mesh(
        grid_initial,
        scalars="elevation",
        cmap="gray",
        opacity=1.0,
        clim=[asc_file_initial.min_height() * 1.1, asc_file_initial.max_height() * 1.1],
        scalar_bar_args={"title": r"Elevation [m]"},
        show_scalar_bar=True,
        smooth_shading=True,
        **light_args
    )

    p.add_mesh(
        grid_final,
        scalars="flow_thickness",
        cmap="magma",
        opacity=1.0,
        scalar_bar_args={"title": r"Flow thickness [m]"},
        show_scalar_bar=True,
        smooth_shading=True,
        **light_args
    )

    # Callback for getting the last camera position
    def take_screenshot(x, path):
        if path is not None:
            x.screenshot(path)

    initial_cpos = None

    cpos = p.show(
        before_close_callback=lambda x: take_screenshot(x, image_outfile),
        auto_close=True,
        return_cpos=True,
        cpos=initial_cpos,
    )

    print("Camera position")
    print(cpos)


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
    contour_spacing: int = typer.Option(
        -1, "-c", help="Contour spacing in meters. Negative to disable contouring"
    ),
    contour_color: str = typer.Option("white", "--ccolor", help="Contour color"),
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
        asc_file_initial.height_data.T,
        levels=50,
        cmap="gray",
    )

    if contour_spacing > 0:
        min_height = asc_file_final.min_height()
        max_height = asc_file_final.max_height()
        plt.contour(
            x_data,
            y_data,
            asc_file_initial.height_data.T,
            levels=np.arange(min_height, max_height, contour_spacing),
            colors=contour_color,
        )

    plt.xlabel("x [m]")
    plt.ylabel("y [m]")

    lava_thickness = asc_file_final.height_data
    lava_thickness[lava_thickness < 1e-10] = float("nan")
    plt.contourf(x_data, y_data, lava_thickness.T, cmap="magma")

    if interactive:
        plt.show()
    if image_outfile is not None:
        plt.savefig(image_outfile, dpi=300)
    else:
        plt.savefig("./image.png")


if __name__ == "__main__":
    app()
