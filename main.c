#include <stdio.h>
#include <time.h>
#include "common.h"
#include "greedy.h"



int main (int argc, const char * argv[]) {
	int * matrix, * colors;
	int nbrnodes, nbrcolors;
	int i;
	clock_t ticks1, ticks2;
	ticks1 = clock();

	read_graph_to_adjacency_matrix(stdin, &matrix, &nbrnodes);
	color_greedy(matrix, nbrnodes, &colors, &nbrcolors);
	printf("%d\n%d\n", nbrnodes, nbrcolors);
	for (i = 0; i < nbrnodes; i++)
		printf("%d\n", colors[i]);

	ticks2	= clock();
    printf("Time Elapsed = %lf\n", 1.0 * (ticks2 - ticks1) / CLOCKS_PER_SEC);
    return 0;
}
