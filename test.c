/*
 *  test.c
 *  P2
 *
 *  Created by Filipe Araujo on 4/28/10.
 *  Copyright 2010 Universidade de Coimbra. All rights reserved.
 *
 */

#include <stdio.h>


#define NODES 10

int colors_OK(int * expectedcolors, int * realcolors, int expnumcolors, int realnumcolors, int nodes) {
	int i = 0;
	if (expnumcolors == realnumcolors)
		for (; i < nodes && expectedcolors[i] == realcolors[i]; i++);
	return i == nodes;
}


int main(int argc, char *argv[]) {
	/* Petersen graph */
	int nodes = NODES;
	int adjmatrix[][NODES] = { 
		/*outer polygon   inner star */
		{0, 1, 0, 0, 1,  1, 0, 0, 0, 0},
		{1, 0, 1, 0, 0,  0, 1, 0, 0, 0},
		{0, 1, 0, 1, 0,  0, 0, 1, 0, 0},
		{0, 0, 1, 0, 1,  0, 0, 0, 1, 0},
		{1, 0, 0, 1, 0,  0, 0, 0, 0, 1},
		{1, 0, 0, 0, 0,  0, 0, 1, 1, 0},
		{0, 1, 0, 0, 0,  0, 0, 0, 1, 1},
		{0 ,0, 1, 0, 0,  1, 0, 0, 0, 1},
		{0, 0, 0, 1, 0,  1, 1, 0, 0, 0},
		{0, 0, 0, 0, 1,  0, 1, 1, 0, 0}
	};
	int expectedcolors[] = {1, 2, 1, 2, 3,  2, 1, 3, 3, 2};
	int *realcolors, numcolors;
	color_greedy(adjmatrix, nodes, &realcolors, &numcolors);
	

	free(realcolors);
	
	if (!colors_OK(expectedcolors, realcolors, 3, numcolors, nodes))
		fprintf(stderr, "Test failed: Petersen with greedy coloring\n");
	
	
	/**** Second test **********/
	color_greedy(NULL, 0, &realcolors, &numcolors);
	if (!colors_OK(NULL, NULL, 0, numcolors, 0))
		fprintf(stderr, "Test failed: empty graph coloring\n");

	/**** Third test **********/
	int adjmatrix2[][2] = { {0, 0}, {0, 0} };
	int expectedcolors2[] = {1, 1};
	color_greedy(adjmatrix2, 2, &realcolors, &numcolors);
	if (!colors_OK(expectedcolors2, realcolors, 1, numcolors, 2))
		fprintf(stderr, "Test failed: Petersen with greedy coloring\n");

	fprintf(stdout, "No news is good news...\n");
	
	return 0;
}