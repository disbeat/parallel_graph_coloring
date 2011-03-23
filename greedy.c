/*
 *  greedy.c
 *  P2
 *
 *  Created by Filipe Araujo on 4/28/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <pbl.h>

#include "common.h"

/**
 *** Assigns colors 1 to n to nodes. Stores result in a matrix called *pcollors
 */
void color_greedy(int * madjacency, int nodes, int ** pcolors, int * numcolors) {
	int node, highestcolor = 1;
	
	if (madjacency == NULL || nodes == 0) {
		*pcolors = NULL;
		*numcolors = 0;
		return;
	}
		
	*pcolors = (int *) calloc(nodes, sizeof(int));
	if (*pcolors == NULL)
		error("Out of memory in color_greedy");
	
	for (node = 0; node < nodes; node++) {
		int color, othernode;
		
		PblSet * psetcolors = (PblSet *) pblSetNewHashSet();
		if (psetcolors == NULL) 
			error("Failed to allocate new set in color_greedy()\n");
		
		for (othernode = 0; othernode < nodes; othernode++) {
			if (othernode != node) {
				if (madjacency[POS(othernode, node, nodes)] == 1) {
					int othercolor = (*pcolors)[othernode];
					if (othercolor != 0)
					/* the neighbor already has a color */
						pblSetAdd(psetcolors, (void *) othercolor);
				}
				
				/* look for the first available color */
				for (color = 1; pblSetContains(psetcolors, (void *) color) && color < nodes; color++);
				(*pcolors)[node] = color;
				if (color > highestcolor)
					highestcolor = color;
			}
		}
		
		pblSetFree(psetcolors);
	}
	*numcolors = highestcolor;
}
