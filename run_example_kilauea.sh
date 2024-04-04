flowy_bin=./build/flowy
input_file=./examples/KILAUEA2014-2015/input.toml
asc_file=./examples/KILAUEA2014-2015/test20m.asc
output_folder=./output_KILAUEA
$flowy_bin $input_file -a $asc_file -o $output_folder

python3 examples/plot_thickness_from_asc_file.py plot-pyvista $output_folder/initial.asc $output_folder/output.asc -i