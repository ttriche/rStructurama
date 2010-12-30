#include "machineInfo.h"

#ifdef WIN32
#include <windows.h>
#else
#include <sys/sysctl.h>
#include <libiberty.h>
#include <cstdio>
#ifdef OMP
#include <omp.h>
#endif
#endif

#define ARRAY_SIZE(a) (sizeof (a) / sizeof ((a)[0]))

int StMachineInfo::numProcs(void) {

#ifdef WIN32
	int numCPU;
	SYSTEM_INFO sysinfo;
	GetSystemInfo( &sysinfo );
	numCPU = sysinfo.dwNumberOfProcessors;
	return numCPU;
#else
#ifdef OMP
  int numCPU = omp_get_num_procs();
#else 
	int numCPU = 0;
  numCPU = sysconf( _SC_NPROCESSORS_ONLN );
  if( numCPU < 1 ) {
    nt mib[4];
    size_t len; 
    mib[0] = CTL_HW;
    mib[1] = HW_AVAILCPU;  // alternatively, try HW_NCPU;
    sysctl(mib, 2, &numCPU, &len, NULL, 0);
    if( numCPU < 1 ) {
      mib[1] = HW_NCPU;
      sysctl( mib, 2, &numCPU, &len, NULL, 0 );
      if( numCPU < 1 ) {
        numCPU = 1;
      }
    }
  }
	if( numCPU < 1 ) numCPU = 1;
	return numCPU;
#endif
#endif
}
