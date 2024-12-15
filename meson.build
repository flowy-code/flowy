project('flowy', 'cpp',
  version : '0.1',
  default_options : ['warning_level=3',
                     'cpp_std=c++20',
                     # 'buildtype=debugoptimize', # debug=true, optimize=2
                     'debug=true',
                     'optimization=3'])

cpp_args = []

cppc = meson.get_compiler('cpp')

cppc_id = cppc.get_id()
message('Got compiler', cppc_id)

if cppc.get_id() != 'intel-llvm'
  cpp_args += cppc.get_supported_arguments([
    # GCC flags
    '-ffast-math',                 # Enable faster, non-IEEE math calculations
    '-fno-finite-math-only',       # Allow Inf and NaN
    # These are based on recommendations from
    # https://youtu.be/_enXuIxuNV4?si=LvtMqPwJ6jYbDY66
    '-fno-semantic-interposition', # Assume no interposition for module functions
    '-fno-plt',                    # Avoid PLT for function calls within shared libs
    '-Bsymbolic',                  # Resolve symbols to local definitions
  ])
endif

if cppc.get_id() == 'intel-llvm'
  # Here: https://www.intel.com/content/www/us/en/docs/cpp-compiler/developer-guide-reference/2021-10/alphabetical-option-list.html
  cpp_args += cppc.get_supported_arguments([
    # Intel flags
    '-xHost',
    '-ipo',
    '-fp-model=fast',
    '-fma',
    '-fiopenmp',
  ])
endif

# Conditional arguments
cpp_args += cppc.get_supported_arguments([
  '-Wno-unused-local-typedefs',  # Ignore unused local typedefs warnings
  '-Wno-array-bounds',           # Suppress out-of-bounds array access warnings
  '-Wno-unused-variable',
  '-Wno-unused-but-set-variable',
  '-Wno-unused-label',
  # From xtensor
  '-Wno-unused-parameter',
  '-Wno-tautological-constant-compare',
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
_libs = []
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

# mkl_dep = dependency('mkl', required: true)
# mkl_libdirs = '/home/rgoswami/intel/oneapi/mkl/latest' #get_option('mkl_libdirs')
# mkl_lib = cppc.find_library('mkl_rt', dirs: mkl_libdirs, required: true)
# _libs += mkl_lib
# static/dynamic or iomp/seq
mkl_dep = dependency('mkl-dynamic-lp64-iomp', required: true)
_deps += mkl_dep

# TODO(rg): Use OpenBlas as a fallback
cpp_args += ['-DHAVE_CBLAS=1', '-DHAVE_MKL']

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
  link_with: _libs,
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
