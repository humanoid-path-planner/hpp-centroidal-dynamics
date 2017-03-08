/*
 * Copyright 2015, LAAS-CNRS
 * Author: Andrea Del Prete
 */

#ifndef _ROBUST_EQUILIBRIUM_LIB_CONFIG_HH
#define _ROBUST_EQUILIBRIUM_LIB_CONFIG_HH

// Package version (header).
# define ROBUST_EQUILIBRIUM_VERSION "UNKNOWN"

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
#  define ROBUST_EQUILIBRIUM_DLLIMPORT __declspec(dllimport)
#  define ROBUST_EQUILIBRIUM_DLLEXPORT __declspec(dllexport)
#  define ROBUST_EQUILIBRIUM_DLLLOCAL
# else
// On Linux, for GCC >= 4, tag symbols using GCC extension.
#  if __GNUC__ >= 4
#   define ROBUST_EQUILIBRIUM_DLLIMPORT __attribute__ ((visibility("default")))
#   define ROBUST_EQUILIBRIUM_DLLEXPORT __attribute__ ((visibility("default")))
#   define ROBUST_EQUILIBRIUM_DLLLOCAL  __attribute__ ((visibility("hidden")))
#  else
// Otherwise (GCC < 4 or another compiler is used), export everything.
#   define ROBUST_EQUILIBRIUM_DLLIMPORT
#   define ROBUST_EQUILIBRIUM_DLLEXPORT
#   define ROBUST_EQUILIBRIUM_DLLLOCAL
#  endif // __GNUC__ >= 4
# endif // defined _WIN32 || defined __CYGWIN__

# ifdef ROBUST_EQUILIBRIUM_STATIC
// If one is using the library statically, get rid of
// extra information.
#  define ROBUST_EQUILIBRIUM_DLLAPI
#  define ROBUST_EQUILIBRIUM_LOCAL
# else
// Depending on whether one is building or using the
// library define DLLAPI to import or export.
#  ifdef ROBUST_EQUILIBRIUM_EXPORTS
#   define ROBUST_EQUILIBRIUM_DLLAPI ROBUST_EQUILIBRIUM_DLLEXPORT
#  else
#   define ROBUST_EQUILIBRIUM_DLLAPI ROBUST_EQUILIBRIUM_DLLIMPORT
#  endif // ROBUST_EQUILIBRIUM_EXPORTS
#  define ROBUST_EQUILIBRIUM_LOCAL ROBUST_EQUILIBRIUM_DLLLOCAL
# endif // ROBUST_EQUILIBRIUM_STATIC

#endif //_ROBUST_EQUILIBRIUM_LIB_CONFIG_HH
