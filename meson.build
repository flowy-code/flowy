project('flowy', 'cpp',
  version : '0.1',
  default_options : ['warning_level=3',
                     'cpp_std=c++20',
                     # 'buildtype=debugoptimize', # debug=true, optimize=2
                     'debug=true',
                     'optimization=3'])

cpp_args = []

cppc = meson.get_compiler('cpp')

cpp_args += cppc.get_supported_arguments([
  '-Wno-unused-local-typedefs',  # Ignore unused local typedefs warnings
  '-Wno-array-bounds',           # Suppress out-of-bounds array access warnings
  '-ffast-math',                 # Enable faster, non-IEEE math calculations
  '-fno-finite-math-only',       # Allow Inf and NaN
  # These are based on recommendations from
  # https://youtu.be/_enXuIxuNV4?si=LvtMqPwJ6jYbDY66
  '-fno-semantic-interposition', # Assume no interposition for module functions
  '-fno-plt',                    # Avoid PLT for function calls within shared libs
  '-Bsymbolic',                  # Resolve symbols to local definitions
])

if cpp_args.length() > 0
  message('Adding compiler flags', cpp_args)
endif

incdir = include_directories(['.', 'thirdparty/tsl'])

_sources = [
  'src/asc_file.cpp',
  'src/simulation.cpp',
  'src/topography.cpp',
  'src/config_parser.cpp',
]

# Library dependencies
_deps = []
pdflib_dep = dependency('pdf_cpplib', fallback : ['pdf_cpplib', 'pdflib_dep'])

## Netcdf
if get_option('with_netcdf')
  # On windows the conda-forge netcdf pkg-config file tells us to link against a mysterious debug.lib
  # Therefore, we use cmake as the preferred method, which does not make problems
  netcdf_dep = dependency('netcdf',  language : 'cpp', method: 'cmake', required: false)
  # On the linux CI, however, cmake cant find netcdf. So we still use pkg-config as a fallback.
  if not netcdf_dep.found()
    netcdf_dep = dependency('netcdf',  language : 'cpp', method: 'pkg-config')
  endif
  _deps += netcdf_dep
  cpp_args += ['-DWITH_NETCDF']
  _sources += ['src/netcdf_file.cpp']
endif

_deps += [
  pdflib_dep,
  dependency('xtensor'),
  dependency('xtensor-blas'),
  dependency('fmt'),
]

symbol_visibility = 'default'
if get_option('default_library') == 'static'
  # Ensure that if we link a static library into a shared library,
  # private symbols don't get re-exported.
  # Kanged from https://github.com/ERGO-Code/HiGHS/pull/1737/files
  symbol_visibility = 'inlineshidden'
endif

# Static libraries for the executable and tests, shared for the bindings
# Generates (linux) libflowy.so && libflowy.a
flowylib = both_libraries('flowy',
  _sources,
  install:true,
  dependencies : _deps,
  include_directories : incdir,
  gnu_symbol_visibility: symbol_visibility,
  pic: true,
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
  ]

  if get_option('with_netcdf')
    tests += [['Test_FileIO', 'test/test_file_io.cpp']]
  endif

  Catch2 = dependency('Catch2', method : 'cmake', modules : ['Catch2::Catch2WithMain', 'Catch2::Catch2'])

  foreach t : tests
    exe = executable(t.get(0), t.get(1),
      dependencies : [flowylib_static_dep, Catch2],
      cpp_args : cpp_args
    )
    test(t.get(0), exe, workdir : meson.project_source_root())
  endforeach
endif
