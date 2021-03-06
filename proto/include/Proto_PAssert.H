#pragma once

#include <cassert>
#include <cstdio>
#include <cstdlib>

// Custom Assert macro with printf-style reporting
#ifdef __CUDA_ARCH__

#define PROTO_ASSERT(stmt,args...)
#define PR_assert(stmt) 
#define PR_error(stmt) 
#define PR_warning(stmt) 

#elif NDEBUG

#define PROTO_ASSERT(stmt,args...)
#define PR_assert(stmt)
#define PR_error(stmt)
#define PR_warning(stmt)


#else

#define PROTO_ASSERT(stmt,args...) \
  assert( (stmt) ||                                         \
     (fprintf(stderr,"%s:%d: error: ",__FILE__,__LINE__) && \
     fprintf(stderr,"Assertion `" #stmt "` failed.\n\t") && \
     fprintf(stderr,args) && \
     fprintf(stderr,"\n"))); \
  if (!(stmt)){std::abort();}


#define PR_assert(stmt) \
  assert( (stmt) ||                                         \
     (fprintf(stderr,"%s:%d: error: ",__FILE__,__LINE__) && \
      fprintf(stderr,"Assertion `" #stmt "` failed.\n\t")));    \
  if (!(stmt)){std::abort();}

#define PR_error(stmt) \
  fprintf(stderr,"error thrown =  `" #stmt "`!!!\n\t"); \
  std::abort();

#define PR_warning(stmt) \
  fprintf(stderr,"warning thrown =  `" #stmt "` !!!\n\t");


#endif
