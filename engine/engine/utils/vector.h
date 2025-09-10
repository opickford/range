#ifndef VECTOR_H
#define VECTOR_H

#include <stdint.h>
#include <stdlib.h>


/*

TODO: Topfile comment.

This vector implementation wraps a raw data array and a capacity.

*/

// TODO: Init function?
#define Vector_destroy(v) free((v).data)

// Define a vector representation of a given type.
// Generates an anonymous struct to represent a vector for a type.
#define Vector(T) struct { T* data; size_t capacity; }

// Reserves bytes for the given capacity.
#define Vector_reserve(v, new_capacity) reserve(&(v).data, &(v).capacity, new_capacity, sizeof(*(v).data))

// Resizes the given vector.
#define Vector_resize(v, new_capacity) resize(&(v).data, &(v).capacity, new_capacity, sizeof(*(v).data))

// Internal helper functions to aid debugging.

// TODO: COMMENTS!!!!!!
static void set_size(void** data,
    size_t* old_capacity,
    size_t new_capacity,
    size_t type_size)
{
    void* temp = realloc(*data, new_capacity * type_size);
    if (!temp)
    {
        // If the allocation failed, just return without updating the 
        // old capacity, this is the only way I can think of returning an
        // error, the user could then assert that the capacity is correct.

        // TODO: Log debug message?
        return;
    }

    // Update the given vector struct componnents.
    *data = temp;
    *old_capacity = new_capacity;
}

// Updates the given data array and old capacity.
static void reserve(
    void** data,
    size_t* old_capacity,
    size_t new_capacity,
    size_t type_size
)
{
    // No more capacity needed.
    if (*data && *old_capacity >= new_capacity) return;

    set_size(data, old_capacity, new_capacity, type_size);
 }

static void resize(void** data,
    size_t* old_capacity,
    size_t new_capacity,
    size_t type_size
)
{
    // No more capacity needed.
    if (*data && *old_capacity == new_capacity) return;

    set_size(data, old_capacity, new_capacity, type_size);
}



#endif
