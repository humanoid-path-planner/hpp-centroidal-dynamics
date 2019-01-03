/*
 * Copyright 2015, LAAS-CNRS
 * Author: Andrea Del Prete
 */

#ifndef _CENTROIDAL_DYNAMICS_LIB_CONFIG_HH
#define _CENTROIDAL_DYNAMICS_LIB_CONFIG_HH

// Package version (header).
# define CENTROIDAL_DYNAMICS_VERSION "UNKNOWN"

// Handle portable symbol export.
// Defining manually which symbol should be exported is required
// under Windows whether MinGW or MSVC is used.
//
// The headers then have to be able to work in two different modes:
// - dllexport when one is building the library,
// - dllimport for clients using the library.
//
// On Linux, set the visibility accordingly. If C++ symbol visibility
// is handled by the compiler, see: http://gcc.gnu.org/wiki/Visibility
# if defined _WIN32 || defined __CYGWIN__
// On Microsoft Windows, use dllimport and dllexport to tag symbols.
#  define CENTROIDAL_DYNAMICS_DLLIMPORT __declspec(dllimport)
#  define CENTROIDAL_DYNAMICS_DLLEXPORT __declspec(dllexport)
#  define CENTROIDAL_DYNAMICS_DLLLOCAL
# else
// On Linux, for GCC >= 4, tag symbols using GCC extension.
#  if __GNUC__ >= 4
#   define CENTROIDAL_DYNAMICS_DLLIMPORT __attribute__ ((visibility("default")))
#   define CENTROIDAL_DYNAMICS_DLLEXPORT __attribute__ ((visibility("default")))
#   define CENTROIDAL_DYNAMICS_DLLLOCAL  __attribute__ ((visibility("hidden")))
#  else
// Otherwise (GCC < 4 or another compiler is used), export everything.
#   define CENTROIDAL_DYNAMICS_DLLIMPORT
#   define CENTROIDAL_DYNAMICS_DLLEXPORT
#   define CENTROIDAL_DYNAMICS_DLLLOCAL
#  endif // __GNUC__ >= 4
# endif // defined _WIN32 || defined __CYGWIN__

# ifdef CENTROIDAL_DYNAMICS_STATIC
// If one is using the library statically, get rid of
// extra information.
#  define CENTROIDAL_DYNAMICS_DLLAPI
#  define CENTROIDAL_DYNAMICS_LOCAL
# else
// Depending on whether one is building or using the
// library define DLLAPI to import or export.
#  ifdef CENTROIDAL_DYNAMICS_EXPORTS
#   define CENTROIDAL_DYNAMICS_DLLAPI CENTROIDAL_DYNAMICS_DLLEXPORT
#  else
#   define CENTROIDAL_DYNAMICS_DLLAPI CENTROIDAL_DYNAMICS_DLLIMPORT
#  endif // CENTROIDAL_DYNAMICS_EXPORTS
#  define CENTROIDAL_DYNAMICS_LOCAL CENTROIDAL_DYNAMICS_DLLLOCAL
# endif // CENTROIDAL_DYNAMICS_STATIC

#endif //_CENTROIDAL_DYNAMICS_LIB_CONFIG_HH
