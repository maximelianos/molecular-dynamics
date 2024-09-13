fs = import('fs')
copy = fs.copyfile('cluster_923.xyz')
copy = fs.copyfile('cluster_3871.xyz')

executable(
	'milestone07',
	'main.cpp',
	include_directories : [lib_incdirs],
	link_with : [lib],
	dependencies : [eigen, mpi]
)