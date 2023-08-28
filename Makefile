solution: solution.c
	# Modify this using any compiler flags you choose, or clang.
	# Note that this will not currently use OpenMP.
	gcc solution.c -o solution

PQ-Dijkstra: PQ-Dijkstra.c
	gcc -O3 PQ-Dijkstra.c -o PQ-Dijkstra
	@echo "**Note** This is only used to build the skelton."
	@echo "To compile your own code, use 'make solution'."
