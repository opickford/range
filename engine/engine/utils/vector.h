#ifndef VECTOR_H
#define VECTOR_H

#include <stdint.h>


/*

TODO: Topfile comment.

This vector implementation wraps a raw data array and a capacity.

*/


// TODO: Init and Destroy functions?

// Generates an anonymous struct, 
#define Vector(T) struct { T* v; size_t capacity; }

// Reserves bytes for the given capacity.
#define Vector_reserve(v, capacity) reserve(&(v).data, &(v).capacity, capacity, sizeof(*(v).data)

// Internal helper function to aid debugging.
// Updates the given data array and old capacity.
static void reserve(
    void** data, 
    size_t* old_capacity,
    size_t new_capacity, 
    size_t type_size
);

#endif
