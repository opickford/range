#ifndef VECTOR_H
#define VECTOR_H

#include <stdint.h>
#include <stdlib.h>


/*

TODO: Topfile comment.

This vector implementation wraps a raw data array and a capacity.

*/


// TODO: Init and Destroy functions?

// Generates an anonymous struct, 
#define Vector(T) struct { T* data; size_t capacity; }

// Reserves bytes for the given capacity.
#define Vector_reserve(v, new_capacity) reserve(&(v).data, &(v).capacity, new_capacity, sizeof(*(v).data))

// Internal helper function to aid debugging.
// Updates the given data array and old capacity.
static void reserve(
    void** data,
    size_t* old_capacity,
    size_t new_capacity,
    size_t type_size
)
{
    // No more capacity needed.
    if (*old_capacity >= new_capacity) return;

    void* temp = realloc(*data, new_capacity * type_size);
    if (!temp)
    {
        // If the allocation failed, just return without updating the 
        // old capacity, this is the only way I can think of returning an
        // error, the user could then assert that the capacity is correct.

        // TODO: Log debug message?
        return;
    }

    *data = temp;
    *old_capacity = new_capacity;
 }

#endif
