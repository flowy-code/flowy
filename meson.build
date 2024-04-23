project('flowy', 'cpp',
  version : '0.1',
  default_options : ['warning_level=3', 'cpp_std=c++20', 'optimization=0'])

cpp_args = []

cpp = meson.get_compiler('cpp')

# clang and gcc specific flags
if cpp.get_id() == 'gcc' or cpp.get_id() == 'clang'
  cpp_args += ['-Wno-unused-local-typedefs', '-Wno-array-bounds', '-ffast-math']
endif

if cpp_args.length() > 0
  message('Adding compiler flags', cpp_args)
endif

incdir = include_directories('.')

_sources = [
  'src/asc_file.cpp',
  'src/simulation.cpp',
  'src/topography.cpp',
  'src/config_parser.cpp',
  'src/netcdf_file.cpp'
]

# Library dependencies
pdflib_dep = dependency('pdf_cpplib', fallback : ['pdf_cpplib', 'pdflib_dep'])

_deps = [
  pdflib_dep,
  dependency('xtensor'),
  dependency('xtensor-blas'),
  dependency('fmt'),
  dependency('tomlplusplus'),
  dependency('netcdf', language : 'cpp'),
]

# Static libraries for the executable and tests, shared for the bindings
# Generates (linux) libflowy.so && libflowy.a
flowylib = both_libraries('flowy',
  _sources,
  install:true,
  dependencies : _deps,
  include_directories : incdir,
  cpp_args : cpp_args,
)
# Generating both dependencies here is not a heavy ask
flowylib_static_dep = declare_dependency(include_directories : incdir,
  link_with : flowylib.get_static_lib(), dependencies: _deps)
flowylib_shared_dep = declare_dependency(include_directories : incdir,
  link_with : flowylib.get_shared_lib(), dependencies: _deps)

# Declare the flowy executable
if get_option('build_exe')
  flowyexe = executable('flowy',
    'src/main.cpp',
    install : true,
    dependencies : [flowylib_static_dep],
    cpp_args : cpp_args
  )
endif

if get_option('build_tests')
  # Tests
  tests = [
    ['Test_ASC', 'test/test_asc.cpp'],
    ['Test_Slope', 'test/test_slope.cpp'],
    ['Test_Simulation', 'test/test_simulation.cpp'],
    ['Test_Topography', 'test/test_topography.cpp'],
    ['Test_Lobe', 'test/test_lobe.cpp'],
    ['Test_Trunc_Normal', 'test/test_truncated_normal.cpp'],
    ['Test_Raster', 'test/test_rasterization.cpp'],
    ['Test_FileIO', 'test/test_file_io.cpp']
  ]

  Catch2 = dependency('Catch2', method : 'cmake', modules : ['Catch2::Catch2WithMain', 'Catch2::Catch2'])

  foreach t : tests
    exe = executable(t.get(0), t.get(1),
      dependencies : [flowylib_static_dep, Catch2],
      cpp_args : cpp_args
    )
    test(t.get(0), exe, workdir : meson.project_source_root())
  endforeach
endif
