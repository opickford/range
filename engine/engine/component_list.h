#ifndef COMPONENT_LIST_H
#define COMPONENT_LIST_H

#include "utils/memory_utils.h"

/*

could be called a ComponentList?

then we use defines to generate the functions?


*/

#include "mesh_instance.h"

#include <string.h>
#include <stdlib.h>


/*

TODO: To make these functions non-inline, i would have to define a macro to 
      generate the functions as well right???

TODO: The annoying naming convention here could be a case for snake_case types
      also it's already inconsistent because I dont use pascal case for the 
      function names


or instead of defining a specific ComponentList for each, define a
ComponentList? use jda's as inspiration?

the issue is the void* if we use a generic ComponentList, we would need to 
worry about alignment, or define functions and cast the void* to the component.

then we define functions to add and remove specifically, and get? keeping the 
alignment easy (and fast?)




TODO: Im a bit worried about shared components, for example a position, 
physics and renderer needs this, but if theyre not aligned then iterating
through would be annoying.

*/

typedef struct
{
    int count;
    int capacity;

    void* elements;

    int* id_to_index;
    int* index_to_id;

    int* free_ids;
    int free_ids_capacity;
    int free_ids_count;

} ComponentList;

void component_list_init(ComponentList* list);
void component_list_destroy(ComponentList* list);

#define DECLARE_COMPONENT(T)            \
typedef int T##ID;                      \
T##ID T##_add(ComponentList* list);

// TODO: Implement removing as well, not tested add yet actually.
// TODO: We can do typedef ComponentList MeshInstances? 
// TODO: TEST with point lights first as not implemented yet!!!!!!!

#define DEFINE_COMPONENT(T) \
T##ID T##_add(ComponentList* list)                                          \
{                                                                           \
    T##ID id;                                                               \
                                                                            \
    /* Grow arrays if at capacity. */                                       \
    if ((list)->count == (list)->capacity)                                  \
    {                                                                       \
        ++(list)->capacity;                                                 \
        T* new_elements = (T*)realloc(                                      \
            (list)->elements,                                               \
            (list)->capacity * sizeof(T));                                  \
                                                                            \
        if (!new_elements)                                                  \
        {                                                                   \
            log_error("Failed to allocate new elements");                   \
            return -1;                                                      \
        }                                                                   \
                                                                            \
        (list)->elements = new_elements;                                    \
                                                                            \
        resize_int_array(&(list)->id_to_index,                              \
            (list)->capacity);                                              \
        resize_int_array(&(list)->index_to_id,                              \
            (list)->capacity);                                              \
                                                                            \
        id = (list)->capacity - 1;                                          \
    }                                                                       \
    else                                                                    \
    {                                                                       \
        id = (list)->free_ids[--(list)->free_ids_count];                    \
    }                                                                       \
                                                                            \
    T* element = &(((T*)(list)->elements)[(list)->count]);                  \
    (list)->id_to_index[id] = (list)->count;                                \
    (list)->index_to_id[(list)->count] = id;                                \
    ++(list)->count;                                                        \
    return id;                                                              \
}

#endif