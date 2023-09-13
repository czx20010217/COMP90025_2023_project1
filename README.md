# COMP90025_2023_project1

gcc -O3 -fopenmp -o test1 test1.c
.\test1 < input1.txt 
.\test1 < grid10x10faster.txt
./test1 < input1.txt

gcc -O3 -fopenmp -o delta delta.c
.\delta < input1.txt



make PQ-Dijkstra
.\PQ-Dijkstra < input1.txt

the performance of code not using O3 is 19-8 10-4