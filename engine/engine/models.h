#ifndef MODELS_H
#define MODELS_H

#include "renderer/render_buffers.h"

#include "common/status.h"
#include "maths/vector3.h"

#include <stdio.h>
#include <stdlib.h>

/*

ModelBase
- Defines object space data: positions, normals, uvs.
- Read from .obj
- Used to create instances.

ModelInstance
- Created from a ModelBase.
- Defines a transform: Position, Direction, Scale.
- Defines visual stuff like colours and textures (although uvs from base).

There are also intermediate buffers to store indexed results for less overall
computations. Also, these can be used by other engine components, like the 
physics system, to avoid recomputation of world space positions etc.

*/

/*

TODO: 

One day, I reckon I will need models to be made up of meshes. Where each mesh can be 
moved seperately. Maybe this is something that would be abstracted from the render functions
I'm not sure. I will need to think about it properly. But let's say we had a basic character
and wanted to do a running animation, it would be nice to be able to animate it's bodyparts,
for example, we could store transforms for parts of the animation and lerp each mesh in the 
model for those. Or reloading etc. This will be quite important.

*/

/*

TODO: Gotta think about the buffers again. Do I really need this many? Surely I can reuse the clipping one because that's the biggest?
	  And I can treat this like a memory arena, this could be a great way to reduce memory. Although remember I need at least two arenas for swapping stuff.

	  Need to think about this more though. Not really sure tbf.

	  Or at least this doesn't need to be one big struct of them right? IDRK.

TODO: I also hate how all of this works really, I'm sure there is some stuff that could be structs where we access it together constantly.
	  Need to learn more about memory. 'What every programmer should know about memory'


*/


typedef struct
{	
	// TODO: Comments and naming. Currently mb = model base, mi = model instance

	// Buffer sizes.
	int mbs_count;
	int mis_count;
	int max_mb_faces;					// The highest number of faces in a mesh, out of all the models. Used for the temporary clipping buffers.

	// ModelBase data
	int mbs_total_faces;				// The total number of faces defined by all mbs.
	int mbs_total_positions;
	int mbs_total_normals;
	int mbs_total_uvs;

	// TODO: Refactor? When iterating through mis, it would be nice to have an offset to the positions of each mi?
	int* mbs_faces_offsets;
	int* mbs_positions_offsets;
	int* mbs_normals_offsets;
	int* mbs_uvs_offsets;

	int* mbs_positions_counts;			// The model matrix transform is model specific, therefore, we must store how many positions the mesh has.
	int* mbs_normals_counts;			// The model normal matrix transform is model specific, therefore, we must store how many normals the mesh has.
	int* mbs_faces_counts;				// Number of faces in the model.
	int* mbs_uvs_counts;

	// TODO: Get rid of face prefix??
	int* mbs_face_position_indices;		// The indices to positions that make up the faces, used for indexed rendering.
	int* mbs_face_normal_indices;		// The indices to normals that make up the faces, used for indexed rendering.
	int* mbs_face_uvs_indices;

	float* mbs_object_space_positions;	// Original vertex positions without any transforms applied.
	float* mbs_object_space_normals;
	float* mbs_uvs;
	float* mbs_lightmap_uvs;
	float* mbs_object_space_centres;

	// Instance data
	// TODO: Naming.
	int mis_total_faces;				// Total number of faces from all mis, keeps track of the size of the buffers.
	int mis_total_positions;
	int mis_total_normals;

	int* mis_base_ids;						// The id of the model base.
	int* mis_texture_ids;					// The id of the texture.
	int* mis_dirty_bounding_sphere_flags;	// If a mi's scale has changed, the bounding sphere centre needs to be recalculated.
	int* mis_intersected_planes;			// For each mi, the number of planes intersected, then the indices of the planes.
	int* mis_passed_broad_phase_flags;		// Whether the mi is visible after broad phase culling. TODO: Name.

	float* mis_vertex_colours;			// Per vertex colours for the instances.
	float* mis_transforms;				// The instance world space transforms: [ Position, Direction, Scale ]
	float* mis_bounding_spheres;		// The bounding sphere for each instance in world space.

	// Transform results buffers.
	// TODO: These are specific to mis, should prefix. - or move to RenderBuffers.
	float* view_space_positions;
	float* view_space_normals;	
	
	
} Models;

// Initialises the models struct.
void models_init(Models* models);

// Parses the obj file for the number of each component.
void parse_obj_counts(FILE* file, int* num_vertices, int* num_uvs, int* num_normals, int* num_faces);

Status mb_from_obj(Models* models, RenderBuffers* rbs, const char* filename);

// TODO: It would be nice to be able to create different model
//		 instances without memory allocating each time. I think
//		 allocating a bigger pool of memory would be nice, then
//		 we can allocate when we need more capacity. So we would
//		 have like a capacity for each buffer size.

// TODO: For this, I can look into memory arenas: https://www.rfleury.com/p/untangling-lifetimes-the-arena-allocator
//		 But we don't need this for now. Essentially just allocate big blocks, store the used and total capacity etc.

// Allocates memory for n instances of the ModelBase at mb_index.
void mi_create(Models* models, RenderBuffers* rbs, int mb_index, int n);

void free_models(Models* models);

// Helpers
void mi_set_transform(Models* models, int mi_index, V3 position, V3 eulers, V3 scale);



#endif