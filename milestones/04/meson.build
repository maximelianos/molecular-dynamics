fs = import('fs')
fs.copyfile('lj54.xyz')

executable(
	'milestone04',
	'main.cpp',
	include_directories : [lib_incdirs],
	link_with : [lib],
	dependencies : [eigen, mpi]
)