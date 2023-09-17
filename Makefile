solution: solution.c
	gcc -Wall -O3 -fopenmp solution.c -o solution

PQ-Dijkstra: PQ-Dijkstra.c
	gcc -O3 PQ-Dijkstra.c -o PQ-Dijkstra
	@echo "**Note** This is only used to build the skelton."
	@echo "To compile your own code, use 'make solution'."
