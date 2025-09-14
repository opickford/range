#ifndef COMMON_H
#define COMMON_H

// TODO: NOt sure what to call this. common.h seems wrong as it's not in common!

#include <stdlib.h>

// Returns a random float from 0-1. TODO: Rename?
inline float random_float()
{
	return (float)rand() / (float)RAND_MAX;
}

#define SWAP(T, a, b) do { T temp = a; a = b; b = temp; } while (0)

#endif