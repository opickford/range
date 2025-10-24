#ifndef TRANSFORM_H
#define TRANSFORM_H

#include "maths/vector3.h"

// TODO: Should this just go into components?

typedef struct transform
{
    v3_t position;
    v3_t scale;
    v3_t rotation;

    // Store previous state for lerp between state updated by physics.
    v3_t previous_position;
    v3_t previous_rotation;
    //v3_t previous_scale;

} transform_t;


inline void transform_init(transform_t* t)
{
    t->position = (v3_t){ 0.f, 0.f, 0.f };
    t->scale = (v3_t){ 1.f, 1.f, 1.f };
    t->rotation = (v3_t){ 0.f, 0.f, 0.f };

    t->previous_position = t->position;
    t->previous_rotation = t->rotation;
    
    /*

TODO: Do we want to lerp scale????

      The issue is that if we do not set the previous_scale to equal scale,
      and we lerp scale (maybe from 1,1,1 to 10,10,10 for example),
      then we do not calculate a new bounding sphere as scale_has_changed
      hasn't been updated, this can cause clipping issues which causes
      out of buffer writes!!!!!!

      To fix this we could just not lerp scale as it doesn't seem necessary.


*/
    //t->previous_scale = t->scale;
}



#endif