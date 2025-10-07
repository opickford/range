#ifndef TIMER_H
#define TIMER_H

#include <time.h>

// Extremely barebones timer functionality.
// TODO: Would be nice to be able to 
typedef struct 
{
	clock_t start;

} timer_t;

inline timer_t timer_start()
{
	timer_t timer =
	{
		.start = clock()
	};
	
	return timer;
}

inline int timer_get_elapsed(timer_t* timer)
{

	return (int)((clock() - timer->start) * 1000 / CLOCKS_PER_SEC);
}

inline void timer_restart(timer_t* timer)
{
	timer->start = clock();
}

#endif

