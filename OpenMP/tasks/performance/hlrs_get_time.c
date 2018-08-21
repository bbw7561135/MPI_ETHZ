#include <sys/times.h> 
#include <sys/time.h>
#include <unistd.h>        
#include <stdio.h>
#include <errno.h>

#define BODY \
	struct timeval tv;\
	struct tms buf;\
\
	if(gettimeofday(&tv,NULL)!=0) {\
		perror("get_time: ");\
		exit(-1);\
	}\
	*real = (double)tv.tv_sec + ((double)tv.tv_usec/1000000.0);\
\
	times(&buf);\
	*user = (double)buf.tms_utime/(double)CLK_TCK \
	      + (double)buf.tms_stime/(double)CLK_TCK;


void hlrs_get_time(double *user, double *real)
{
	BODY
}

void hlrs_get_time_(double *user, double *real)
{
	BODY
}

void HLRS_GET_TIME_(double *user, double *real)
{
	BODY
}

void HLRS_GET_TIME(double *user, double *real)
{
	BODY
}

void hlrs_get_time__(double *user, double *real)
{
	BODY
}

