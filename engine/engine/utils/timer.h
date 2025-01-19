#ifndef TIMER_H
#define TIMER_H

#include <time.h>

// Extremely barebones timer functionality.
// TODO: Would be nice to be able to 
typedef struct 
{
	clock_t start;

} Timer;

inline Timer timer_start()
{
	Timer timer =
	{
		.start = clock()
	};
	
	return timer;
}

inline clock_t timer_get_elapsed(Timer* timer)
{

	return (clock() - timer->start) * 1000 / CLOCKS_PER_SEC;
}

inline void timer_restart(Timer* timer)
{
	timer->start = clock();
}

#endif

