project('flowy', 'cpp',
  version : '0.1',
  default_options : ['warning_level=3', 'cpp_std=c++20', 'optimization=3'])

cpp_args = ['-Wno-unused-local-typedefs', '-Wno-array-bounds', '-ffast-math']

incdir = include_directories('include')

_sources = [
  'src/asc_file.cpp',
  'src/simulation.cpp',
  'src/topography.cpp',
  'src/config_parser.cpp'
]

# Library dependencies
pdflib_dep = dependency('pdf_cpplib', fallback : ['pdf_cpplib', 'pdflib_dep'])

_deps = [
  pdflib_dep,
  dependency('xtensor'), 
  dependency('xtensor-blas'), 
  dependency('fmt'), 
  dependency('tomlplusplus')
]

flowyexe = executable('flowy', 
  _sources + 'src/main.cpp',
  install : true,
  dependencies : _deps,
  include_directories : incdir,
  cpp_args : cpp_args
)

flowylib = library('flowylib', 
  _sources, 
  install:true, 
  dependencies : _deps, 
  include_directories : incdir,
  cpp_args : cpp_args
)

flowylib_dep = declare_dependency(include_directories : incdir,
  link_with : flowylib, dependencies: _deps)

tests = [
  ['Test_ASC', 'test/test_asc.cpp'],
  ['Test_Slope', 'test/test_slope.cpp'],
  ['Test_Simulation', 'test/test_simulation.cpp'],
  ['Test_Topography', 'test/test_topography.cpp'],
  ['Test_Lobe', 'test/test_lobe.cpp'],
]

Catch2 = dependency('Catch2', method : 'cmake', modules : ['Catch2::Catch2WithMain', 'Catch2::Catch2'])

foreach t : tests
  exe = executable(t.get(0), _sources + t.get(1),
    dependencies : _deps + Catch2,
    include_directories : incdir,
    cpp_args : cpp_args
  )
  test(t.get(0), exe, workdir : meson.project_source_root())
endforeach