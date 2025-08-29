#include "vector.h"

#include "utils/logger.h"

#include <Windows.h>

/*
void Vector_init(Vector* a, uint32_t element_size)
{
    memset(a, 0, sizeof(Vector));
    a->element_size = element_size;
}





void Vector_destroy(Vector* a)
{
    if (a->size > 0 && a->data)
    {
        free(a->data);
    }
    
}*/

// Static functions
static void reserve(void** data, size_t* old_capacity, size_t new_capacity, size_t type_size)
{
    // No more capacity needed.
    if (*old_capacity <= new_capacity) return;

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