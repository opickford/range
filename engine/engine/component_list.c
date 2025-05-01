#include "component_list.h"

#include <string.h>
#include <stdlib.h>

void component_list_init(ComponentList* list)
{
    memset(list, 0, sizeof(ComponentList));
}

void component_list_destroy(ComponentList* list)
{
    free(list->elements);
    free(list->free_ids);
    free(list->id_to_index);
    free(list->index_to_id);
}