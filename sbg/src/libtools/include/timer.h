/***************************************************************************
                          timer.h  -  description
                             -------------------
    begin                : Mon May 17 2004
    copyright            : (C) 2004 by Jose Ignacio Garzon
    email                :
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef TIMER_H
#define TIMER_H

#include <ctime>
#include <limits>
#include <stdio.h>
/**
  *@author Jose Ignacio Garzon
  *
  * Class that allows to measure the run time
  */

//
// not accurated timer....good for >> 57 minutes
//
class timer {
public:
  ///Constructor
  timer()
  {
    time(&_start_time);
  }
  ///Destructor
  ~timer(){};
  ///Resets the initial time for measure
  void restart()
  {
    time(&_start_time);
  }

  ///Returns the seconds elapsed from the initial time to the call of this method
  double elapsed() const
  {
    time_t _finish_time;
    time(&_finish_time);
    return difftime(_finish_time,_start_time);
  }

  ///Returns a string with a representation of seconds, minutes and hours of the
  ///time elapsed from the initial time to the call of this method
  char *print_time() {

  static char buf[50];
  tm time_elapsed;
  double seconds;
  time_elapsed.tm_sec=0;
  time_elapsed.tm_min=0;
  time_elapsed.tm_hour=0;
  time_elapsed.tm_mday=0;
  time_elapsed.tm_mon=0;
  time_elapsed.tm_year=0;
  time_elapsed.tm_wday=0;
  time_elapsed.tm_yday=0;
  time_elapsed.tm_isdst=0;
  time_elapsed.tm_sec=0;

  time_t _finish_time;

  time(&_finish_time);

  seconds=(int)(difftime(_finish_time,_start_time));
  if(seconds>=3600)
  {
    time_elapsed.tm_hour=((int)(seconds))/3600;
    seconds=(int)(seconds)%3600;
  }
  if(seconds>=60)
  {
    time_elapsed.tm_min=((int)(seconds))/60;
    seconds=(int)(seconds)%60;
  }
  time_elapsed.tm_sec=(int)seconds;

  strftime(buf,50," %Hh. %M' %S''",&time_elapsed);
  return buf;
}

///Returns a string with a representation of seconds of the
///time elapsed from the initial time to the call of this method
char *print_time_sec() {
/* timing function modified from FFTW */
double x=elapsed();
static char buf[128];
sprintf(buf, "%f sec", x);
return buf;
}

private:
  ///Initial time
  //std::clock_t _start_time;
  time_t _start_time;
};

//
// high accuracy timer but only works for times less igual to 57 minutes
//
class Htimer {
public:
  ///Constructor
  Htimer()
  {
    _start_time= std::clock();
  }
  ///Destructor
  ~Htimer(){};
  ///Resets the initial time for measure
  void restart()
  {
    _start_time= std::clock();
  }

  ///Returns the seconds elapsed from the initial time to the call of this method
  double elapsed() const
  {
    return double(std::clock() - _start_time)/ CLOCKS_PER_SEC;

  }

  ///Returns the maximum number of seconds measurable by this class
  double elapsed_max() const
  {
    return(double(std::numeric_limits<std::clock_t>::max()) - double(_start_time)) / double(CLOCKS_PER_SEC);
  }

  ///Returns the minimum measurable unit
  double elapsed_min() const
  {
    return double(1)/ double(CLOCKS_PER_SEC);
  }

  ///Returns a string with a representation of seconds, minutes and hours of the
  ///time elapsed from the initial time to the call of this method
  char *print_time() {
    /* timing function modified from FFTW */
    double x=elapsed();

    static char buf[128];
    if (x < 1.0E-6)
      sprintf(buf, "%f ns", x * 1.0E9);
    else if (x < 1.0E-3)
      sprintf(buf, "%f us", x * 1.0E6);
    else if (x < 1.0)
      sprintf(buf, "%f ms", x * 1.0E3);
    else if (x < 60.0)
      sprintf(buf, "%f s", x);

    return buf;

  }

//Returns a string with a representation of seconds of the
///time elapsed from the initial time to the call of this method
  char *print_time_sec() {
    /* timing function modified from FFTW */
    double x=elapsed();
    static char buf[128];

     sprintf(buf, "%.2f sec", x);


    return buf;
  }

  // Mon made (18/9/2008)
  // Returns the number of clocks elapsed since "restart()" or begining.
    time_t clocks() {
    	return time_t(std::clock() - _start_time);
    }

  // Mon made (18/9/2008)
  // Returns the number of seconds elapsed since "restart()" or begining.
    double secs() {
    	return double( (std::clock() - _start_time)/(double)CLOCKS_PER_SEC );
    }


  private:
    ///Initial time
    std::clock_t _start_time;
};


#endif
