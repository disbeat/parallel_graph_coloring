/*
 *  common.c
 *  P2
 *
 *  Created by Filipe Araujo on 4/28/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include "common.h"

void error(char * msg) {
	fprintf(stderr, "%s\n", msg);
	exit(1);
}


/**
 *** We assume that there are no arcs (i.e., edges are undirected). This makes no difference, even if graphs have arcs.
 **/
void read_graph_to_adjacency_matrix(FILE * f, int ** pmatrix, int * pnodes) {
	int edges, i;
	fscanf(f, " %d", pnodes);
	fscanf(f, " %d", &edges);
	*pmatrix = (int *) calloc(*pnodes * *pnodes, sizeof(int));
	if (*pmatrix == NULL) {
		error("Out of memory in read_graph_to_adjacency_matrix()\n");
		MPI_Finalize();
		exit(-1);
	}
	
	/* read edges one by one */
	for (i = 0; i < edges; i++) {
		int node1, node2;
		int n = fscanf(f, " %d %d", &node1, &node2);
		if (n != 2 || node1 < 0 || node1 >= *pnodes || node2 < 0 || node2 >= *pnodes) {
			fprintf(stderr, "Wrong edge %d in the graph definition file\n", i);
			exit(1);
		}
		(*pmatrix)[POS(node1, node2, *pnodes)] = 1;
		(*pmatrix)[POS(node2, node1, *pnodes)] = 1;
	}
}
