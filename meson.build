project('flowy', 'cpp',
  version : '0.1',
  default_options : ['warning_level=3', 'cpp_std=c++20', 'optimization=3'])

add_global_arguments(['-Wno-unused-local-typedefs', '-Wno-array-bounds'], language : 'cpp')

incdir = include_directories('include')

flowy_sources = [
  'src/asc_file.cpp',
  'src/simulation.cpp',
  'src/topography.cpp',
  'src/config_parser.cpp'
]

# Library dependencies
pdflib_dep = dependency('pdf_cpplib', fallback : ['pdf_cpplib', 'pdflib_dep'])

flowy_deps = [
  pdflib_dep,
  dependency('xtensor'), 
  dependency('xtensor-blas'), 
  dependency('fmt'), 
  dependency('tomlplusplus')
]

flowy_executable = executable('flowy', flowy_sources + 'src/main.cpp',
  install : true,
  dependencies : flowy_deps,
  include_directories : incdir
)

flowy_shared_lib = library('flowy', 
  flowy_sources, 
  install:true, 
  dependencies : flowy_deps, 
  include_directories : incdir
)

tests = [
  ['Test_ASC', 'test/test_asc.cpp'],
  ['Test_Slope', 'test/test_slope.cpp'],
  ['Test_Simulation', 'test/test_simulation.cpp'],
  ['Test_Topography', 'test/test_topography.cpp'],
  ['Test_Lobe', 'test/test_lobe.cpp'],
]

Catch2 = dependency('Catch2', method : 'cmake', modules : ['Catch2::Catch2WithMain', 'Catch2::Catch2'])

foreach t : tests
  exe = executable(t.get(0), flowy_sources + t.get(1),
    dependencies : flowy_deps + Catch2,
    include_directories : incdir
  )
  test(t.get(0), exe, workdir : meson.project_source_root())
endforeach