/*
 * Copyright 2015, LAAS-CNRS
 * Author: Andrea Del Prete
 */

#ifndef WIN32
#include <sys/time.h>
#else
#include <Windows.h>
#endif

#include <stdio.h>
#include <iostream>
#include <iomanip>      // std::setprecision
#include <boost/algorithm/string.hpp>
#include <centroidal-dynamics-lib/logger.hh>

namespace robust_equilibrium
{
  using namespace std;

  Logger& getLogger()
  {
    static Logger l(0.001, 1.0);
    return l;
  }

  Logger::Logger(double timeSample, double streamPrintPeriod)
    : m_timeSample(timeSample),
      m_streamPrintPeriod(streamPrintPeriod),
      m_printCountdown(0.0)
  {
#ifdef LOGGER_VERBOSITY_ERROR
    m_lv = VERBOSITY_ERROR;
#endif
#ifdef LOGGER_VERBOSITY_WARNING_ERROR
    m_lv = VERBOSITY_WARNING_ERROR;
#endif
#ifdef LOGGER_VERBOSITY_INFO_WARNING_ERROR
    m_lv = VERBOSITY_INFO_WARNING_ERROR;
#endif
#ifdef LOGGER_VERBOSITY_ALL
    m_lv = VERBOSITY_ALL;
#endif
  }

  void Logger::countdown()
  {
    if(m_printCountdown<0.0)
      m_printCountdown = m_streamPrintPeriod;
    m_printCountdown -= m_timeSample;
  }

  void Logger::sendMsg(string msg, MsgType type, const char* file, int line)
  {
//    if(m_lv==VERBOSITY_NONE ||
//       (m_lv==VERBOSITY_ERROR && !isErrorMsg(type)) ||
//       (m_lv==VERBOSITY_WARNING_ERROR && !(isWarningMsg(type) || isErrorMsg(type))) ||
//       (m_lv==VERBOSITY_INFO_WARNING_ERROR && isDebugMsg(type)))
//      return;

    // if print is allowed by current verbosity level
    if(isStreamMsg(type))
    {
      // check whether counter already exists
      string id = file+toString(line);
      map<string,double>::iterator it = m_stream_msg_counters.find(id);
      if(it == m_stream_msg_counters.end())
      {
        // if counter doesn't exist then add one
        m_stream_msg_counters.insert(make_pair(id, 0.0));
        it = m_stream_msg_counters.find(id);
      }

      // if counter is greater than 0 then decrement it and do not print
      if(it->second>0.0)
      {
        it->second -= m_timeSample;
        return;
      }
      else  // otherwise reset counter and print
        it->second = m_streamPrintPeriod;
    }

    vector<string> fields;
    boost::split(fields, file, boost::is_any_of("/"));
    const char* file_name = fields[fields.size()-1].c_str();

    if(isErrorMsg(type))
      printf("[ERROR %s %d] %s\n", file_name, line, msg.c_str());
    else if(isWarningMsg(type))
      printf("[WARNING %s %d] %s\n", file_name, line, msg.c_str());
    else if(isInfoMsg(type))
      printf("[INFO %s %d] %s\n", file_name, line, msg.c_str());
    else
      printf("[DEBUG %s %d] %s\n", file_name, line, msg.c_str());

    fflush(stdout); // Prints to screen or whatever your standard out is
  }

  bool Logger::setTimeSample(double t)
  {
    if(t<=0.0)
      return false;
    m_timeSample = t;
    return true;
  }

  bool Logger::setStreamPrintPeriod(double s)
  {
    if(s<=0.0)
      return false;
    m_streamPrintPeriod = s;
    return true;
  }

} // namespace robust_equilibrium

