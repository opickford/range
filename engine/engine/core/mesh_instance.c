#include "mesh_instance.h"

#include "strides.h"

#include "utils/logger.h"

#include <stdlib.h>
#include <string.h>

Status mesh_instance_init(MeshInstance* mi)
{
	memset(mi, 0, sizeof(MeshInstance));

	mi->scale = (V3){ 1,1,1 };
	mi->has_scale_changed = 1; // TODO: Rename to recalc bounding sphere?

	return STATUS_OK;
}

Status mesh_instance_set_base(MeshInstance* mi, const MeshBase* mb)
{
	// Grow the vertex albedos buffer, there should be one albedo
	// per vertex.
	Status status = resize_float_array(&mi->vertex_alebdos, mb->num_faces * STRIDE_FACE_VERTICES * STRIDE_COLOUR);

	if (STATUS_OK == status)
	{
		mi->mb_id = mb->id;

		// Default albedo to white.
		mesh_instance_set_albedo(mi, mb, (V3) { 1.f, 1.f, 1.f });
	}

	return status;
}

void mesh_instance_set_albedo(MeshInstance* mi, const MeshBase* mb, V3 albedo)
{
	for (int i = 0; i < mb->num_faces * STRIDE_FACE_VERTICES * STRIDE_COLOUR; i += STRIDE_COLOUR)
	{
		mi->vertex_alebdos[i] = albedo.x;
		mi->vertex_alebdos[i + 1] = albedo.y;
		mi->vertex_alebdos[i + 2] = albedo.z;
	}
}

void mesh_instance_destroy(MeshInstance* mi)
{
	free(mi->vertex_alebdos);
}

Status mesh_instances_init(MeshInstances* mis)
{
	memset(mis, 0, sizeof(MeshInstances));

	return STATUS_OK;
}

void mesh_instances_destroy(MeshInstances* mis)
{
	free(mis->instances);
}

MeshInstance* mesh_instances_get(MeshInstances* mis, MeshInstanceID mi_id)
{
    return &mis->instances[mis->id_to_index[mi_id]];
}

MeshInstanceID mesh_instances_add(MeshInstances* mis)
{
    MeshInstanceID mi_id;

	// Grow the arrays if we're at maximum capacity already.
    if (mis->count == mis->capacity)
    {
        // Resize the array of mesh instances.
        ++mis->capacity;
        MeshInstance* new_instances = realloc(
            mis->instances, mis->capacity * sizeof(MeshInstance));

        if (!new_instances)
        {
            log_error("Failed to allocate for new mesh instance.");
            return 0;
        }

        mis->instances = new_instances;

        // Resize the ids map.
        resize_int_array(&mis->id_to_index, mis->capacity);
        resize_int_array(&mis->index_to_id, mis->capacity);

        // Assign a new id.
        mi_id = mis->capacity - 1;
    }
    else
    {
        // Spare capacity means there must be a free id, so reuse it.
        mi_id = mis->free_ids[--mis->free_ids_count];
    }

    // Create a new instance.
    MeshInstance* mi = &mis->instances[mis->count];
    mesh_instance_init(mi);

    mis->id_to_index[mi_id] = mis->count;
    mis->index_to_id[mis->count] = mi_id;

	++mis->count;

	return mi_id;
}

void mesh_instances_remove(MeshInstances* mis, MeshInstanceID mi_id)
{
    // Ensure id is valid.
    if (mi_id >= mis->capacity || mi_id == -1)
    {
        return;
    }

    // Read the actual index of the instance.
    const int index_to_remove = mis->id_to_index[mi_id];

    if (index_to_remove == -1)
    {
        // Instance already removed.
        return;
    }

    // Free up the instance id to be re-used.
    if (mis->free_ids_count == mis->free_ids_capacity)
    {
        ++mis->free_ids_capacity;
        resize_int_array(&mis->free_ids, mis->free_ids_capacity);
    }
    mis->free_ids[mis->free_ids_count++] = mi_id;

    // Clear the maps of the mesh instance.
    mis->id_to_index[mi_id] = -1;
    mis->index_to_id[index_to_remove] = -1;

    // Destroy the instance.
    mesh_instance_destroy(&mis->instances[index_to_remove]);

    const int last_mi_index = mis->count - 1;

    --mis->count;

    // If we removed the last one, no need to do any swapping.
    if (index_to_remove == last_mi_index)
    {   
        return;
    }

    // Copy the last instance into the place of the one we just removed, this
    // ensures tight packing of the array.
    MeshInstanceID last_id = mis->index_to_id[last_mi_index];

    // This should never fail.
    if (last_id != -1)
    {
        mis->id_to_index[last_id] = index_to_remove;
        mis->index_to_id[index_to_remove] = last_id;
    }

    // Copy the data over.
    memcpy(&mis->instances[index_to_remove],
        &mis->instances[last_mi_index],
        sizeof(MeshInstance));
}
