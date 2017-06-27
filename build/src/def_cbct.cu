#include "def_cbct.h"

#include <stdio.h>

__constant__ cbct_cg g_cg;
__constant__ cbct_ig g_ig;

__constant__ float c_dt2;
__constant__ float c_ds2;
__constant__ float c_dz2;

texture<float, cudaTextureType2DLayered, cudaReadModeElementType> texRefLayered;
texture<float, 3, cudaReadModeElementType> tex3DRef;

__global__ void show_kernel()
{
	printf("\nImage geometry:\n");
	printf("[nx, ny, nz] : [%d, %d, %d]\n", g_ig.nxyz.x, g_ig.nxyz.y, g_ig.nxyz.z);
	printf("[dx, dy, dz] : [%.3f, %.3f, %.3f] mm\n", g_ig.dxyz.x, g_ig.dxyz.y, g_ig.dxyz.z);
	printf("offset_[x, y, z]: [%.3f, %.3f, %.3f] mm\n", g_ig.offset.x, g_ig.offset.y, g_ig.offset.z);

	printf("\nCone-beam geometry:\n");
	printf("[dso, dsd] : [%.3f, %.3f] mm\n", g_cg.dso, g_cg.dsd);
	printf("%s-detector\n", "Flat");
	printf("[ns, nt] : [%d, %d]\n", g_cg.nxy.x, g_cg.nxy.y);
	printf("[ds, dt] : [%.3f, %.3f]\n", g_cg.dxy.x, g_cg.dxy.y);
	printf("offset_[s, t]: [%.3f, %.3f] mm\n\n", g_cg.offset.x, g_cg.offset.y);
}

void cbct::initCg(const cbct_cg& cg)
{
	cudaMemcpyToSymbol(g_cg, &cg, sizeof(cbct_cg));

	float  dt2 = fabsf(cg.dxy.y*0.5f);
	float  ds2 = cg.dxy.x*0.5f;
	cudaMemcpyToSymbol(c_dt2, &dt2, sizeof(float));
	cudaMemcpyToSymbol(c_ds2, &ds2, sizeof(float));

	CUDA_CHECK_ERROR("cbct::initCg");
}

void cbct::initIg(const cbct_ig& ig)
{
	//cbct_ig cp = ig;
	cudaMemcpyToSymbol(g_ig, &ig, sizeof(cbct_ig));

	float  dz2 = ig.dxyz.z*0.5f;
	cudaMemcpyToSymbol(c_dz2, &dz2, sizeof(float));

	CUDA_CHECK_ERROR("cbct::initIg");
}

void cbct::showGeometry()
{
	show_kernel<<<1,1>>>();
	CUDA_CHECK_ERROR("cbct::showGeometry");
}

void cbct::bindProjToTexLayered(cudaArray* arr)
{
	texRefLayered.addressMode[0] = cudaAddressModeClamp;
	texRefLayered.addressMode[1] = cudaAddressModeClamp;
	texRefLayered.filterMode = cudaFilterModeLinear;
	texRefLayered.normalized = false;

	cudaChannelFormatDesc desc = cudaCreateChannelDesc<float>();
	cudaBindTextureToArray(texRefLayered, arr, desc);

	CUDA_CHECK_ERROR("cbct::bindProjToTexLayered");
}

void cbct::bindImgToTex3D(cudaArray* arr)
{
	tex3DRef.addressMode[0] = cudaAddressModeClamp;
	tex3DRef.addressMode[1] = cudaAddressModeClamp;
	tex3DRef.addressMode[2] = cudaAddressModeClamp;
	tex3DRef.filterMode = cudaFilterModeLinear;
	tex3DRef.normalized = false;

	cudaChannelFormatDesc desc = cudaCreateChannelDesc<float>();
	cudaBindTextureToArray(tex3DRef, arr, desc);

	CUDA_CHECK_ERROR("cbct::bindImgToTex3D");
}



