/*
 * time.h
 *
 *  Created on: Feb 10, 2017
 *      Author: mon
 */

#ifndef TIME_H_
#define TIME_H_


#pragma once

#ifdef _WIN32
// Windows version of timer
#include <windows.h>
class timer
{
private:
	LARGE_INTEGER start;
    LARGE_INTEGER stop;
	LARGE_INTEGER frequency;
public:
	timer(void);
	void startTimer();
	void stopTimer();
	double getElapsedTime(); // Elapsed time in seconds
};
#else
// Linux version of timer
#include <time.h>
class timerReal
{
private:
	struct timespec start;
	struct timespec stop;
public:
	void startTimer();
	void stopTimer();
	double getElapsedTime(); // Elapsed time in seconds
};
#endif


#ifdef _WIN32
// Windows version of the timer
void timer::startTimer()
{
    QueryPerformanceCounter(&start);
}

void timer::stopTimer()
{
    QueryPerformanceCounter(&stop);
}

double timer::getElapsedTime()
{
	LARGE_INTEGER time;
	time.QuadPart = stop.QuadPart - start.QuadPart;
    return ((double)time.QuadPart /(double)frequency.QuadPart);
}

timer::timer(void)
{
	QueryPerformanceFrequency(&frequency);
}

#else
// Linux version of the timer
void timerReal::startTimer()
{
    clock_gettime(CLOCK_REALTIME,&start);
}

void timerReal::stopTimer()
{
    clock_gettime(CLOCK_REALTIME,&stop);
}

double timerReal::getElapsedTime()
{
	double elapsed;
	elapsed = (double)(stop.tv_sec - start.tv_sec) + (double)(stop.tv_nsec - start.tv_nsec) / 1000000000.0;
	return elapsed;
}
#endif



#endif /* TIME_H_ */
