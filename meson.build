project('flowtastic', 'cpp',
  version : '0.1',
  default_options : ['warning_level=3', 'cpp_std=c++20'])

add_global_arguments(['-Wno-unused-local-typedefs', '-Wno-array-bounds'], language : 'cpp')

incdir = include_directories('include')

sources = [
  'src/flowtastic.cpp',
  'src/asc_file.cpp',
  'src/simulation.cpp',
  'src/topography.cpp'
]

# Library dependencies
pdflib_dep = dependency('pdf_cpplib', fallback : ['pdf_cpplib', 'pdflib_dep'])

exe = executable('flowtastic', sources + 'src/main.cpp',
  install : true,
  dependencies : [pdflib_dep, dependency('xtensor'), dependency('xtensor-blas'), dependency('fmt'), dependency('tomlplusplus')],
  include_directories : incdir
  )

tests = [
  ['Test_ASC', 'test/test_asc.cpp'],
  ['Test_Slope', 'test/test_slope.cpp'],
]

Catch2 = dependency('Catch2', method : 'cmake', modules : ['Catch2::Catch2WithMain', 'Catch2::Catch2'])

foreach t : tests
  exe = executable(t.get(0), sources + t.get(1),
    dependencies : [pdflib_dep, dependency('xtensor'), dependency('xtensor-blas'),dependency('fmt'), dependency('tomlplusplus'), Catch2],
    include_directories : incdir
  )
  test(t.get(0), exe, workdir : meson.project_source_root())
endforeach