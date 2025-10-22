#ifndef SCENE_H
#define SCENE_H

#include "mesh_instance.h"
#include "lights.h"

#include "maths/vector3.h"

#include "common/status.h"

// TODO: This should more be storing metadata of the scene, the ecs contains lights and instances etc.
typedef struct scene
{
	mesh_bases_t mesh_bases; // TODO: Should these simply be global? 
                          // TODO: Make chds_vec

	v3_t ambient_light;
    int bg_colour;

} scene_t;

// TODO: refactor to remove the 

status_t scene_init(scene_t* scene);
status_t scene_destroy(scene_t* scene);

// mesh_base_t API Wrappers
mesh_base_id_t scene_add_mesh_base(scene_t* scene);

#endif