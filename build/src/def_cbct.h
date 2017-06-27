#pragma once


#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <math.h>
#include "helper_math.h"
#include "cuMemMgr.h"

#include <string>

#include "Matrices.h"

#ifndef SQR
#define SQR(x) (x)*(x)
#endif

#ifndef RAD2DEG 
#define RAD2DEG 180.0f/(float)M_PI //57.295779513082   
#endif

#ifndef DEG2RAD
#define DEG2RAD (float)M_PI/180.0f //0.01745329252
#endif

//#define USEREG
#define MEDTRONIC

#define ORDER_TS
#define ORDER_ZXY

// image geometry
using cbct_ig = struct cbct_ig
{
	int3   nxyz;
	float3 dxyz;
	float3 invdxyz;
	float3 offset;
	float3 wxyz;
};

// cone-beam geometry
using cbct_cg = struct cbct_cg
{
	float dso;
	float dsd;
	int2  nxy;
	float2  dxy;
	float2 offset;
	float2 wxy;
	float2 invdxy;
};

using lut_area = struct lut_area
{
	float* lut;
	int2 dim;
	float2 delta;
};

using lut_height = struct lut_height
{
	float* lut;
	int3 dim;
	float3 delta;
};

typedef struct
{
	float4 m[3];
} float3x4;

struct Ray
{
	float3 o; // origin
	float3 d; // direction
};


// CUDA utils (def_cbct.cpp)
void CUDA_CHECK_ERROR(const char* name);
std::string CUDA_SET_BEST_DEVICE();

void show_cbct_ig(const cbct_ig& ig);
void show_cbct_cg(const cbct_cg& cg);

int iDivUp(int a, int b);
unsigned int nextPowerOf2(unsigned int n);
float GUPS(int3 voldim, int na, float sec);

float2 rotateCCW_z(float2 v, const float cs, const float sn);

namespace cbct
{
	// def_cbct.cu
	void initCg(const cbct_cg& cg);
	void initIg(const cbct_ig& ig);
	void showGeometry();

	// pre-defined geometry
	void Medtronic(cbct_cg& cg);
	void TimeGeometry(cbct_cg& cg, cbct_ig& ig, float scale);
	void TimeGeometry2(cbct_cg& cg, cbct_ig& ig);


	// cuda, global
	void bindProjToTexLayered(cudaArray* arr);
	void bindImgToTex3D(cudaArray* arr);

} // cbct

namespace pm
{
	Matrix4 extrinsic(const float sad, float deg, float dThetaX, float dThetaY, float dThetaZ);
	Matrix4 invExtrinsic(const float sad, float deg, float dThetaX, float dThetaY, float dThetaZ);
	Matrix4 intrinsic(const float sdd, float ox, float oy);
	Matrix4 PM(const float sdd, const float sad,
		const float ox, const float oy, const float deg,
		const float dThetaX, const float dThetaY, const float dThetaZ,
		const float dox, const float doy,
		const float dsdd, const float dsad);
}



