#ifndef STATUS_H
#define STATUS_H

#include <assert.h>

// TODO: In the future it would be nice to have some prefix for engine/game code to differentiate it 
//       from Windows code for example there are lots of STATUS_ from <Windows.h>.

// TODO: I think this should be in a 'core' folder?

typedef enum
{
	STATUS_OK,
	STATUS_FAILURE,
	STATUS_INVALID_ARGUMENT,
	STATUS_WIN32_FAILURE,
	STATUS_ALLOC_FAILURE,
    STATUS_FILE_FAILURE

} Status;

inline const char* status_to_str(Status status)
{
    switch (status) 
    {
        case STATUS_OK:                 return "Success";
        case STATUS_FAILURE:            return "Failure";
        case STATUS_INVALID_ARGUMENT:   return "Invalid Argument";
        case STATUS_WIN32_FAILURE:      return "Win32 Failure";
        case STATUS_ALLOC_FAILURE:      return "Alloc Failure";
        case STATUS_FILE_FAILURE:       return "File Failure";
        default:                        return "Unknown Status";
    }
}

#ifndef NDEBUG
inline void Assert(Status status) { assert(status == STATUS_OK); };
#else
inline void Assert(Status status) {};
#endif



#endif