lib_sources = [  # All source files (excluding headers)
    'hello.cpp',
    'verlet.cpp',
    'lj_direct_summation.cpp',
    'xyz.cpp',
    'thermostat.cpp',
    'lattice.cpp',
    'neighbors.cpp',
    'lj.cpp',
    'ducastelle.cpp'
]

lib = static_library(
    'lib',  # Unique name
    lib_sources,
    dependencies : [eigen, mpi]
)

lib_incdirs = include_directories('.')
