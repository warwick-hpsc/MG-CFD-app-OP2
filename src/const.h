#ifndef CONST_H
#define CONST_H

#ifndef OP_FUN_PREFIX
    #ifdef __CUDACC__
        #define OP_FUN_PREFIX __host__ __device__
    #else
        #define OP_FUN_PREFIX
    #endif
#endif

#define VALUE_TO_STRING(x) #x
#define VALUE(x) VALUE_TO_STRING(x)
#define VAR_NAME_VALUE(var) #var "=" VALUE( (var) )
#define XMACRO_TO_STR(s) MACRO_TO_STR( s )
#define MACRO_TO_STR(s) #s
#define DO_PRAGMA(x) _Pragma ( #x )

#ifdef __clang__
#define restrict __restrict__
#endif

/*
 * Options
 *
 */
#define GAMMA 1.4f

#define NDIM 3

#define RK 3	// 3rd order RK
#define ff_mach 1.2f
#define deg_angle_of_attack 0.0f

/*
 * not options
 */
#define VAR_DENSITY 0
#define VAR_MOMENTUM  1
#define VAR_DENSITY_ENERGY (VAR_MOMENTUM+NDIM)
#define NVAR (VAR_DENSITY_ENERGY+1)

#define PI 3.1415926535897931f

#define MG_UP 0
#define MG_DOWN 1

#define MESH_FVCORR 0
#define MESH_LA_CASCADE 1
#define MESH_ROTOR37 2
#define MESH_M6WING 3

#endif
