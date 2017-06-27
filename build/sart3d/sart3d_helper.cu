#include "def_cbct.h"

extern __constant__ cbct_cg g_cg;
extern __constant__ cbct_ig g_ig;

extern texture<float, cudaTextureType2D, cudaReadModeElementType> geodivTex;
extern texture<float, cudaTextureType2D, cudaReadModeElementType> amp1Tex;

__global__ static void
correction_kernel(
	float* d_proj,
	const float* d_measurement,
	const float* d_norm,
	const int pitch, 
	const float beta, // only for SF
	const int method
)
{
	const int is = blockIdx.x*blockDim.x + threadIdx.x;
	const int it = blockIdx.y*blockDim.y + threadIdx.y;

	if (is >= g_cg.nxy.x || it >= g_cg.nxy.y)
		return;

	const int index = is*pitch + it;

	const float p = d_proj[index];
	const float m = d_measurement[index];
	const float n = d_norm[index];

	float sc;
	if (method == 0 || method == 1 || method == 2) // LUT
	{
		sc = tex2D(geodivTex, is + 0.5f, it + 0.5f);
	}
	else // SF, A1
	{
		//const float s = (is - g_cg.wxy.x)*g_cg.dxy.x;
		//const float t = (it - g_cg.wxy.y)*g_cg.dxy.y;

		//float l_theta = atanf(t / (sqrtf(SQR(s) + SQR(g_cg.dsd))));
		//l_theta = 1.0f / fabsf(cosf(l_theta));

		//float l_phi = atan2f(s, g_cg.dsd) - beta;
		//l_phi = g_ig.dxyz.x / fmaxf(fabsf(cosf(l_phi)), fabsf(sinf(l_phi)));

		//sc = l_theta*l_phi;
		sc = tex2D(amp1Tex, is + 0.5f, it + 0.5f);
	}

	const float out = (n > 0.0f) ? sc*((m - sc*p) / (sc*n)) : 0.0f;
	d_proj[index] = out;
}

__global__ static void
update_kernel(
	float* d_recon,
	const float* d_back,
	const float* d_norm,
	const int pitch, const float lambda
)
{
	const int ix = blockIdx.y;
	const int iy = blockIdx.z;
	const int iz = blockIdx.x*blockDim.x + threadIdx.x;

	if (ix >= g_ig.nxyz.x || iy >= g_ig.nxyz.y || iz >= g_ig.nxyz.z)
		return;

	//(ix + iy*g_ig.nxyz.x)*pitch + iz
	const int index = (ix + iy*g_ig.nxyz.x)*pitch + iz;
	const float reval = d_recon[index];
	const float bpval = d_back[index];
	const float nval = d_norm[index];

	const float out = (nval > 0.0f) ? lambda*(bpval / nval) : 0.0f;
	d_recon[index] = fmaxf(reval + out, 0.0f);
}



extern "C"
void sart3d_correct(
	float* d_proj,
	float* d_measurement,
	float* d_norm,
	const int pitch, const int2 nxy, 
	const float beta,
	const int method
)
{
	dim3 blk(16, 16);
	dim3 grd(iDivUp(nxy.x, blk.x), iDivUp(nxy.y, blk.y));

	correction_kernel<<<grd,blk>>>(d_proj, d_measurement, d_norm, pitch, beta, method);
}

extern "C"
void sart3d_update(
	float* d_recon,
	const float* d_back,
	const float* d_norm,
	const int pitch, const int3 nxyz, const float lambda
)
{
	dim3 blk(nxyz.z, 1, 1);
	dim3 grd(iDivUp(nxyz.z, blk.x), nxyz.x, nxyz.y);

	update_kernel<<<grd,blk>>>(d_recon, d_back, d_norm, pitch, lambda);
}