/*
 * memstuff.cpp
 *
 *  Created on: Oct 15, 2013
 *      Author: mon
 */
#include <stdio.h>
#include <stdlib.h>

void check_pointer(void *ptr,char *info)
{
	if(!ptr)
	{
		fprintf(stderr,"Memory allocation failed!\n");
		fprintf(stderr,"%s\n",info);
		fprintf(stderr,"\nForcing exit from check_pointer()\n");
		exit(1);
	}
}






