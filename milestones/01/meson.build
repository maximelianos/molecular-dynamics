fs = import('fs')
copy = fs.copyfile('simulation_test_input.txt')

executable(
	'milestone01',
	'main.cpp',
	include_directories : [lib_incdirs],
	link_with : [lib],
	dependencies : [eigen, mpi]
)