google_test = dependency('gtest')

test_sources = [
  'test_hello_world.cpp',
  'test_verlet.cpp',
  'gtest_mpi_main.cpp',
  'test_lj_direct_summation.cpp',
  'test_thermostat.cpp',
  'test_neighbors.cpp',
  'test_ducastelle.cpp'
]

tests = executable(
	'tests',
	test_sources,
	include_directories : [lib_incdirs],
	link_with : [lib],
	dependencies : [eigen, google_test, mpi]
)

test('Tests', tests)
