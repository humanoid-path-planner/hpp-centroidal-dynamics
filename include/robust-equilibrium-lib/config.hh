#ifndef _ROBUST_EQUILIBRIUM_LIB_CONFIG_HH
#define _ROBUST_EQUILIBRIUM_LIB_CONFIG_HH

#include <Eigen/Dense>
#include <Eigen/src/Core/util/Macros.h>

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

namespace robust_equilibrium
{

//#define USE_FLOAT 1;
#ifdef USE_FLOAT
typedef float value_type;
#else
typedef double value_type;
#endif

typedef Eigen::Matrix <value_type, 2, 1>                                            Vector2;
typedef Eigen::Matrix <value_type, 1, 2>                                            RVector2;
typedef Eigen::Matrix <value_type, 3, 1>                                            Vector3;
typedef Eigen::Matrix <value_type, 1, 3>                                            RVector3;
typedef Eigen::Matrix <value_type, Eigen::Dynamic, 1>                               VectorX;
typedef Eigen::Matrix <value_type, 1, Eigen::Dynamic>                               RVectorX;
typedef Eigen::Matrix <value_type, 3, 3, Eigen::RowMajor>                           Rotation;
typedef Eigen::Matrix <value_type, Eigen::Dynamic, 3, Eigen::RowMajor>              MatrixX3;
typedef Eigen::Matrix <value_type, 4, 3, Eigen::RowMajor>                           Matrix43;
typedef Eigen::Matrix <value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXX;

//define Eigen ref if available
#if EIGEN_VERSION_AT_LEAST(3,2,0)
typedef Eigen::Ref<Vector2>     Ref_vector2;
typedef Eigen::Ref<Vector3>     Ref_vector3;
typedef Eigen::Ref<VectorX>     Ref_vectorX;
typedef Eigen::Ref<Rotation>    Ref_rotation;
typedef Eigen::Ref<MatrixX3>    Ref_matrixX3;
typedef Eigen::Ref<Matrix43>    Ref_matrix43;
typedef Eigen::Ref<MatrixXX>    Ref_matrixXX;

typedef const Eigen::Ref<const Vector2>     & Cref_vector2;
typedef const Eigen::Ref<const Vector3>     & Cref_vector3;
typedef const Eigen::Ref<const VectorX>     & Cref_vectorX;
typedef const Eigen::Ref<const Rotation>    & Cref_rotation;
typedef const Eigen::Ref<const MatrixX3>    & Cref_matrixX3;
typedef const Eigen::Ref<const Matrix43>    & Cref_matrix43;
typedef const Eigen::Ref<const MatrixXX>    & Cref_matrixXX;
#else
typedef vector2_t   & Ref_vector2;
typedef vector3_t   & Ref_vector3;
typedef vector_t    & Ref_vector;
typedef rotation_t  & Ref_rotation;
typedef T_rotation_t& Ref_matrixX3;
typedef T_rotation_t& Ref_matrix43;
typedef matrix_t    & Ref_matrixXX;

typedef const vector2_t   & Cref_vector2;
typedef const vector3_t   & Cref_vector3;
typedef const vector_t    & Cref_vector;
typedef const rotation_t  & Cref_rotation;
typedef const T_rotation_t& Cref_matrixX3;
typedef const T_rotation_t& Cref_matrix43;
typedef const matrix_t    & Cref_matrixXX;
#endif

} //namespace robust_equilibrium

#endif //_ROBUST_EQUILIBRIUM_LIB_CONFIG_HH
