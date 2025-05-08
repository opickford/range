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

#define DECLARE_COMPONENT(T)                          \
typedef int T##ID;                                    \
                                                      \
T##ID T##s_add(ComponentList* list);                  \
void T##s_remove(ComponentList* list, T##ID id);      \
void T##_destroy(T* component);                       \
                                                      \
inline T* T##s_get(ComponentList* list)               \
{                                                     \
    return (T*)((list)->elements);                    \
}                                                     



// TODO: I'm realy not sure about this 2 levels of indirection with the ids, but is there
//       a way around it??? Unless I actually update the ids in the 'entity'
//       or wherever we're storing it?


// TODO: Implement removing as well, not tested add yet actually.
// TODO: We can do typedef ComponentList MeshInstances? 
// TODO: TEST with point lights first as not implemented yet!!!!!!!

#define DEFINE_COMPONENT(T) \
T##ID T##s_add(ComponentList* list)                                         \
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
}                                                                           \
                                                                            \
void T##s_remove(ComponentList* list, T##ID id)                             \
{                                                                           \
    /* Ensure id is valid. */                                               \
    if ((id) >= (list)->capacity || (id) == -1) return;                     \
                                                                            \
    /* Read the actual index of the component. */                           \
    const int index_to_remove = (list)->id_to_index[(id)];                  \
                                                                            \
    /* Check if component has already been removed. */                      \
    if (index_to_remove == -1) return;                                      \
                                                                            \
    /* Free up the id to be re-used. */                                     \
    if ((list)->free_ids_count == (list)->free_ids_capacity)                \
    {                                                                       \
        ++(list)->free_ids_capacity;                                        \
        resize_int_array(&(list)->free_ids, (list)->free_ids_capacity);     \
    }                                                                       \
    (list)->free_ids[(list)->free_ids_count++] = (id);                      \
                                                                            \
    /* Clear the maps of the component. */                                  \
    (list)->id_to_index[(id)] = -1;                                         \
    (list)->index_to_id[index_to_remove] = -1;                              \
                                                                            \
    /* Destroy the component. */                                            \
    T* elements = (T*)(list)->elements;                                     \
    T##_destroy(&elements[index_to_remove]);                                \
                                                                            \
    const int last_component_index = (list)->count - 1;                     \
                                                                            \
    --(list)->count;                                                        \
                                                                            \
    /* If we removed the last one, no need to do any swapping. */           \
    if (index_to_remove == last_component_index) return;                    \
                                                                            \
    /* Copy the last instance into the place of the one we just removed, */ \
    /* this ensures tight packing of the array. */                          \
    T##ID last_id = (list)->index_to_id[last_component_index];              \
                                                                            \
    /* This should never fail. */                                           \
    if (last_id != -1)                                                      \
    {                                                                       \
        (list)->id_to_index[last_id] = index_to_remove;                     \
        (list)->index_to_id[index_to_remove] = last_id;                     \
    }                                                                       \
                                                                            \
    /* Copy the data over. */                                               \
    memcpy(&elements[index_to_remove],                                      \
        &elements[last_component_index],                                    \
        sizeof(T));                                                         \
}

#endif