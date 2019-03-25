/*
 * Copyright 2015, LAAS-CNRS
 * Author: Andrea Del Prete
 */

#ifndef __centroidal_dynamics_lib_logger_H__
#define __centroidal_dynamics_lib_logger_H__

/* --------------------------------------------------------------------- */
/* --- INCLUDE --------------------------------------------------------- */
/* --------------------------------------------------------------------- */

#include <hpp/centroidal-dynamics/local_config.hh>
#include <sstream>
#include <Eigen/Dense>
#include <map>
#include "boost/assign.hpp"

namespace centroidal_dynamics
{

//#define LOGGER_VERBOSITY_ERROR
//#define LOGGER_VERBOSITY_WARNING_ERROR
//#define LOGGER_VERBOSITY_INFO_WARNING_ERROR
//#define LOGGER_VERBOSITY_ALL
#define LOGGER_VERBOSITY_ALL

#define SEND_MSG(msg,type)         getLogger().sendMsg(msg,type,__FILE__,__LINE__)

#ifdef LOGGER_VERBOSITY_ERROR
#define SEND_DEBUG_MSG(msg)
#define SEND_INFO_MSG(msg)
#define SEND_WARNING_MSG(msg)
#define SEND_ERROR_MSG(msg)           SEND_MSG(msg,MSG_TYPE_ERROR)
#define SEND_DEBUG_STREAM_MSG(msg)
#define SEND_INFO_STREAM_MSG(msg)
#define SEND_WARNING_STREAM_MSG(msg)
#define SEND_ERROR_STREAM_MSG(msg)    SEND_MSG(msg,MSG_TYPE_ERROR_STREAM)
#endif

#ifdef LOGGER_VERBOSITY_WARNING_ERROR
#define SEND_DEBUG_MSG(msg)
#define SEND_INFO_MSG(msg)
#define SEND_WARNING_MSG(msg)         SEND_MSG(msg,MSG_TYPE_WARNING)
#define SEND_ERROR_MSG(msg)           SEND_MSG(msg,MSG_TYPE_ERROR)
#define SEND_DEBUG_STREAM_MSG(msg)
#define SEND_INFO_STREAM_MSG(msg)\
#define SEND_WARNING_STREAM_MSG(msg)  SEND_MSG(msg,MSG_TYPE_WARNING_STREAM)
#define SEND_ERROR_STREAM_MSG(msg)    SEND_MSG(msg,MSG_TYPE_ERROR_STREAM)
#endif

#ifdef LOGGER_VERBOSITY_INFO_WARNING_ERROR
#define SEND_DEBUG_MSG(msg)
#define SEND_INFO_MSG(msg)            SEND_MSG(msg,MSG_TYPE_INFO)
#define SEND_WARNING_MSG(msg)         SEND_MSG(msg,MSG_TYPE_WARNING)
#define SEND_ERROR_MSG(msg)           SEND_MSG(msg,MSG_TYPE_ERROR)
#define SEND_DEBUG_STREAM_MSG(msg)
#define SEND_INFO_STREAM_MSG(msg)     SEND_MSG(msg,MSG_TYPE_INFO_STREAM)
#define SEND_WARNING_STREAM_MSG(msg)  SEND_MSG(msg,MSG_TYPE_WARNING_STREAM)
#define SEND_ERROR_STREAM_MSG(msg)    SEND_MSG(msg,MSG_TYPE_ERROR_STREAM)
#endif

#ifdef LOGGER_VERBOSITY_ALL
#define SEND_DEBUG_MSG(msg)           SEND_MSG(msg,MSG_TYPE_DEBUG)
#define SEND_INFO_MSG(msg)            SEND_MSG(msg,MSG_TYPE_INFO)
#define SEND_WARNING_MSG(msg)         SEND_MSG(msg,MSG_TYPE_WARNING)
#define SEND_ERROR_MSG(msg)           SEND_MSG(msg,MSG_TYPE_ERROR)
#define SEND_DEBUG_STREAM_MSG(msg)    SEND_MSG(msg,MSG_TYPE_DEBUG_STREAM)
#define SEND_INFO_STREAM_MSG(msg)     SEND_MSG(msg,MSG_TYPE_INFO_STREAM)
#define SEND_WARNING_STREAM_MSG(msg)  SEND_MSG(msg,MSG_TYPE_WARNING_STREAM)
#define SEND_ERROR_STREAM_MSG(msg)    SEND_MSG(msg,MSG_TYPE_ERROR_STREAM)
#endif

  /** Enum representing the different kind of messages.
       */
  enum CENTROIDAL_DYNAMICS_DLLAPI MsgType
  {
    MSG_TYPE_DEBUG          =0,
    MSG_TYPE_INFO           =1,
    MSG_TYPE_WARNING        =2,
    MSG_TYPE_ERROR          =3,
    MSG_TYPE_DEBUG_STREAM   =4,
    MSG_TYPE_INFO_STREAM    =5,
    MSG_TYPE_WARNING_STREAM =6,
    MSG_TYPE_ERROR_STREAM   =7
  };

  template<typename T>
  std::string toString(const T& v)
  {
    std::stringstream ss;
    ss<<v;
    return ss.str();
  }

  template<typename T>
  std::string toString(const std::vector<T>& v, const std::string separator=", ")
  {
    std::stringstream ss;
    for(int i=0; i<v.size()-1; i++)
      ss<<v[i]<<separator;
    ss<<v[v.size()-1];
    return ss.str();
  }

  template<typename T, int n>
  std::string toString(const Eigen::MatrixBase<T>& v, const std::string separator=", ")
  {
    if(v.rows()>v.cols())
      return toString(v.transpose(), separator);
    std::stringstream ss;
    ss<<v;
    return ss.str();
  }

  enum CENTROIDAL_DYNAMICS_DLLAPI LoggerVerbosity
  {
    VERBOSITY_ALL,
    VERBOSITY_INFO_WARNING_ERROR,
    VERBOSITY_WARNING_ERROR,
    VERBOSITY_ERROR,
    VERBOSITY_NONE
  };

  /** A simple class for logging messages
      */
  class CENTROIDAL_DYNAMICS_DLLAPI Logger
  {
  public:

    /** Constructor */
    Logger(double timeSample=0.001, double streamPrintPeriod=1.0);

    /** Destructor */
    ~Logger(){}

    /** Method to be called at every control iteration
           * to decrement the internal Logger's counter. */
    void countdown();

    /** Print the specified message on standard output if the verbosity level
         * allows it. The file name and the line number are used to identify
         * the point where sendMsg is called so that streaming messages are
         * printed only every streamPrintPeriod iterations.
         */
    void sendMsg(std::string msg, MsgType type, const char* file="", int line=0);

    /** Set the sampling time at which the method countdown()
           * is going to be called. */
    bool setTimeSample(double t);

    /** Set the time period for printing of streaming messages. */
    bool setStreamPrintPeriod(double s);

    /** Set the verbosity level of the logger. */
    void setVerbosity(LoggerVerbosity lv);

  protected:
    LoggerVerbosity m_lv;                /// verbosity of the logger
    double          m_timeSample;        /// specify the period of call of the countdown method
    double          m_streamPrintPeriod; /// specify the time period of the stream prints
    double          m_printCountdown;    /// every time this is < 0 (i.e. every _streamPrintPeriod sec) print stuff

    /** Pointer to the dynamic structure which holds the collection of streaming messages */
    std::map<std::string, double> m_stream_msg_counters;

    bool isStreamMsg(MsgType m)
    { return m==MSG_TYPE_ERROR_STREAM || m==MSG_TYPE_DEBUG_STREAM || m==MSG_TYPE_INFO_STREAM || m==MSG_TYPE_WARNING_STREAM; }

    bool isDebugMsg(MsgType m)
    { return m==MSG_TYPE_DEBUG_STREAM || m==MSG_TYPE_DEBUG; }

    bool isInfoMsg(MsgType m)
    { return m==MSG_TYPE_INFO_STREAM || m==MSG_TYPE_INFO; }

    bool isWarningMsg(MsgType m)
    { return m==MSG_TYPE_WARNING_STREAM || m==MSG_TYPE_WARNING; }

    bool isErrorMsg(MsgType m)
    { return m==MSG_TYPE_ERROR_STREAM || m==MSG_TYPE_ERROR; }
  };

  /** Method to get the logger (singleton). */
  Logger& getLogger();

} // namespace centroidal_dynamics

#endif // #ifndef __sot_torque_control_trajectory_generators_H__
