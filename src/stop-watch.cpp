/*
Copyright (c) 2010-2013 Tommaso Urli

Tommaso Urli    tommaso.urli@uniud.it   University of Udine

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/

#include "hpp/centroidal-dynamics/Stdafx.hh"

#ifndef WIN32
#include <sys/time.h>
#else
#include <Windows.h>

#include <iomanip>
#endif

#include <iomanip>  // std::setprecision

#include "hpp/centroidal-dynamics/stop-watch.hh"

using std::map;
using std::ostringstream;
using std::string;

Stopwatch& getProfiler() {
  static Stopwatch s(REAL_TIME);  // alternatives are CPU_TIME and REAL_TIME
  return s;
}

Stopwatch::Stopwatch(StopwatchMode _mode) : active(true), mode(_mode) {
  records_of = new map<string, PerformanceData>();
}

Stopwatch::~Stopwatch() { delete records_of; }

void Stopwatch::set_mode(StopwatchMode new_mode) { mode = new_mode; }

bool Stopwatch::performance_exists(string perf_name) {
  return (records_of->find(perf_name) != records_of->end());
}

long double Stopwatch::take_time() {
  if (mode == CPU_TIME) {
    // Use ctime
    return clock();

  } else if (mode == REAL_TIME) {
    // Query operating system

#ifdef WIN32
    /*	In case of usage under Windows */
    FILETIME ft;
    LARGE_INTEGER intervals;

    // Get the amount of 100 nanoseconds intervals elapsed since January 1, 1601
    // (UTC)
    GetSystemTimeAsFileTime(&ft);
    intervals.LowPart = ft.dwLowDateTime;
    intervals.HighPart = ft.dwHighDateTime;

    long double measure = intervals.QuadPart;
    measure -= 116444736000000000.0;  // Convert to UNIX epoch time
    measure /= 10000000.0;            // Convert to seconds

    return measure;
#else
    /* Linux, MacOS, ... */
    struct timeval tv;
    gettimeofday(&tv, NULL);

    long double measure = tv.tv_usec;
    measure /= 1000000.0;  // Convert to seconds
    measure += tv.tv_sec;  // Add seconds part

    return measure;
#endif

  } else {
    // If mode == NONE, clock has not been initialized, then throw exception
    throw StopwatchException("Clock not initialized to a time taking mode!");
  }
}

void Stopwatch::start(string perf_name) {
  if (!active) return;

  // Just works if not already present
  records_of->insert(make_pair(perf_name, PerformanceData()));

  PerformanceData& perf_info = records_of->find(perf_name)->second;

  // Take ctime
  perf_info.clock_start = take_time();

  // If this is a new start (i.e. not a restart)
  //  if (!perf_info.paused)
  //    perf_info.last_time = 0;

  perf_info.paused = false;
}

void Stopwatch::stop(string perf_name) {
  if (!active) return;

  long double clock_end = take_time();

  // Try to recover performance data
  if (!performance_exists(perf_name))
    throw StopwatchException("Performance not initialized.");

  PerformanceData& perf_info = records_of->find(perf_name)->second;

  perf_info.stops++;
  long double lapse = clock_end - perf_info.clock_start;

  if (mode == CPU_TIME) lapse /= (double)CLOCKS_PER_SEC;

  // Update last time
  perf_info.last_time = lapse;

  // Update min/max time
  if (lapse >= perf_info.max_time) perf_info.max_time = lapse;
  if (lapse <= perf_info.min_time || perf_info.min_time == 0)
    perf_info.min_time = lapse;

  // Update total time
  perf_info.total_time += lapse;
}

void Stopwatch::pause(string perf_name) {
  if (!active) return;

  long double clock_end = clock();

  // Try to recover performance data
  if (!performance_exists(perf_name))
    throw StopwatchException("Performance not initialized.");

  PerformanceData& perf_info = records_of->find(perf_name)->second;

  long double lapse = clock_end - perf_info.clock_start;

  // Update total time
  perf_info.last_time += lapse;
  perf_info.total_time += lapse;
}

void Stopwatch::reset_all() {
  if (!active) return;

  map<string, PerformanceData>::iterator it;

  for (it = records_of->begin(); it != records_of->end(); ++it) {
    reset(it->first);
  }
}

void Stopwatch::report_all(int precision, std::ostream& output) {
  if (!active) return;

  output << "\n*** PROFILING RESULTS [ms] (min - avg - max - lastTime - "
            "nSamples) ***\n";
  map<string, PerformanceData>::iterator it;
  for (it = records_of->begin(); it != records_of->end(); ++it) {
    report(it->first, precision, output);
  }
}

void Stopwatch::reset(string perf_name) {
  if (!active) return;

  // Try to recover performance data
  if (!performance_exists(perf_name))
    throw StopwatchException("Performance not initialized.");

  PerformanceData& perf_info = records_of->find(perf_name)->second;

  perf_info.clock_start = 0;
  perf_info.total_time = 0;
  perf_info.min_time = 0;
  perf_info.max_time = 0;
  perf_info.last_time = 0;
  perf_info.paused = false;
  perf_info.stops = 0;
}

void Stopwatch::turn_on() {
  std::cout << "Stopwatch active." << std::endl;
  active = true;
}

void Stopwatch::turn_off() {
  std::cout << "Stopwatch inactive." << std::endl;
  active = false;
}

void Stopwatch::report(string perf_name, int precision, std::ostream& output) {
  if (!active) return;

  // Try to recover performance data
  if (!performance_exists(perf_name))
    throw StopwatchException("Performance not initialized.");

  PerformanceData& perf_info = records_of->find(perf_name)->second;

  const long unsigned int MAX_NAME_LENGTH = 60;
  string pad = "";
  for (long unsigned int i = perf_name.length(); i < MAX_NAME_LENGTH; i++)
    pad.append(" ");

  output << perf_name << pad;
  output << std::fixed << std::setprecision(precision)
         << (perf_info.min_time * 1e3) << "\t";
  output << std::fixed << std::setprecision(precision)
         << (perf_info.total_time * 1e3 / (long double)perf_info.stops) << "\t";
  output << std::fixed << std::setprecision(precision)
         << (perf_info.max_time * 1e3) << "\t";
  output << std::fixed << std::setprecision(precision)
         << (perf_info.last_time * 1e3) << "\t";
  output << std::fixed << std::setprecision(precision) << perf_info.stops
         << std::endl;

  //	ostringstream stops;
  //	stops << perf_info.stops;
  //	output << "  *  Stops " << stops.str() << std::endl;
  //	output << std::endl;
}

long double Stopwatch::get_time_so_far(string perf_name) {
  // Try to recover performance data
  if (!performance_exists(perf_name))
    throw StopwatchException("Performance not initialized.");

  long double lapse =
      (take_time() - (records_of->find(perf_name)->second).clock_start);

  if (mode == CPU_TIME) lapse /= (double)CLOCKS_PER_SEC;

  return lapse;
}

long double Stopwatch::get_total_time(string perf_name) {
  // Try to recover performance data
  if (!performance_exists(perf_name))
    throw StopwatchException("Performance not initialized.");

  PerformanceData& perf_info = records_of->find(perf_name)->second;

  return perf_info.total_time;
}

long double Stopwatch::get_average_time(string perf_name) {
  // Try to recover performance data
  if (!performance_exists(perf_name))
    throw StopwatchException("Performance not initialized.");

  PerformanceData& perf_info = records_of->find(perf_name)->second;

  return (perf_info.total_time / (long double)perf_info.stops);
}

long double Stopwatch::get_min_time(string perf_name) {
  // Try to recover performance data
  if (!performance_exists(perf_name))
    throw StopwatchException("Performance not initialized.");

  PerformanceData& perf_info = records_of->find(perf_name)->second;

  return perf_info.min_time;
}

long double Stopwatch::get_max_time(string perf_name) {
  // Try to recover performance data
  if (!performance_exists(perf_name))
    throw StopwatchException("Performance not initialized.");

  PerformanceData& perf_info = records_of->find(perf_name)->second;

  return perf_info.max_time;
}

long double Stopwatch::get_last_time(string perf_name) {
  // Try to recover performance data
  if (!performance_exists(perf_name))
    throw StopwatchException("Performance not initialized.");

  PerformanceData& perf_info = records_of->find(perf_name)->second;

  return perf_info.last_time;
}
