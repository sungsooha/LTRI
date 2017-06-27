#include "cbctProjector.h"
#include "mexMat.h"

extern __constant__ cbct_cg g_cg;
extern __constant__ cbct_ig g_ig;

extern __constant__ float c_dz2;
extern __constant__ float c_dt2;

extern texture<float, cudaTextureType3D, cudaReadModeElementType> tex3DRef;
extern texture<float, cudaTextureType2DLayered, cudaReadModeElementType> texRefLayered;

texture<float, cudaTextureType2D, cudaReadModeElementType> amp1Tex;

texture<float, cudaTextureType2D, cudaReadModeElementType> geodivTex;
texture<float, cudaTextureType2D, cudaReadModeElementType> areaTex;
texture<float, cudaTextureType3D, cudaReadModeElementType> heightTex;


texture<float4, cudaTextureType1D, cudaReadModeElementType> lineTex;
texture<float4, cudaTextureType1D, cudaReadModeElementType> planeTex;
texture<float , cudaTextureType1D, cudaReadModeElementType> polarTex;

//__constant__ float c_polar[1024];
//__constant__ float4 c_plane[1024];
__constant__ float2 c_sn_cs_view;

__constant__ float2 c_alutOffset;
__constant__ float  c_voxBase;
__constant__ int    c_nlines;

__constant__ float3 c_hlutOffset;
__constant__ int    c_nplanes;

__constant__ float2 c_src;
__constant__ float2 c_det;
__constant__ float2 c_uv_s;
__constant__ float2 c_uv_t;
__constant__ float  c_viewang;




#include "lut_helper_kernels.cuh"
#include "lut_proj_kernels.cuh"
#include "lut_back_kernels.cuh"
#include "lut_geodiv_kernels.cuh"

__global__ static void
update_lines_pos_kernel(float4* lines)
{
	const int is = blockIdx.x*blockDim.x + threadIdx.x;
	if (is >= c_nlines)
		return;

	const float s0 = (0 - g_cg.wxy.x - .5f)*g_cg.dxy.x;
	const float s = s0 + is*g_cg.dxy.x;

	const float2 Ps = c_det + s*c_uv_s;
	const float2 rayvec = Ps - c_src;

	float ang = atan2f(rayvec.y, rayvec.x) * (360.0f / (2.0f*(float)M_PI));
	if (ang < 0.0f)
		ang += 360.0f;

	const float A = Ps.y - c_src.y;
	const float B = c_src.x - Ps.x;
	const float C = c_src.x*Ps.y - Ps.x*c_src.y;
	lines[is] = make_float4(ang, A, B, C);

	//const float cs = cosf(ang*float(DEG2RAD));
	//const float sn = sinf(ang*float(DEG2RAD));
	//const float t = sn*c_src.x + cs*c_src.y;
	//lines[is] = make_float4(ang, sn, cs, t);
}


#include "sf_helper_kernels.cuh"
#include "sf_proj_kernels.cuh"
#include "sf_back_kernels.cuh"
#include "sf_amp1_kernels.cuh"





void cbctProjector::_bindGeoDiv()
{
	geodivTex.normalized = false;
	geodivTex.filterMode = cudaFilterModeLinear;
	geodivTex.addressMode[0] = cudaAddressModeClamp;
	geodivTex.addressMode[1] = cudaAddressModeClamp;

	cudaBindTextureToArray(geodivTex, _mm.getArray("geodiv"));
}

void cbctProjector::_bindAmp1()
{
	amp1Tex.normalized = false;
	amp1Tex.filterMode = cudaFilterModeLinear;
	amp1Tex.addressMode[0] = cudaAddressModeClamp;
	amp1Tex.addressMode[1] = cudaAddressModeClamp;

	cudaBindTextureToArray(amp1Tex, _mm.getArray("amparr"));
}

void cbctProjector::_bindAreaLut(const lut_area* alut)
{
	_nlines = p_cg->nxy.x + 1 + 2;

	_mm.allocArray<float>("alut", alut->dim);
	_mm.copyToArray("alut", alut->lut);

	_mm.allocLinear<float4>("lines", _nlines);

	areaTex.normalized = false;
	areaTex.filterMode = cudaFilterModeLinear;
	areaTex.addressMode[0] = cudaAddressModeClamp;
	areaTex.addressMode[1] = cudaAddressModeClamp;

	cudaBindTextureToArray(areaTex, _mm.getArray("alut"));
	cudaBindTexture(0, lineTex, _mm.getLinear<float4>("lines"), lineTex.channelDesc);

	float2 alutOffset = make_float2(1.0f / alut->delta.x, 1.0f / alut->delta.y);
	float  voxBase = fabsf(p_ig->dxyz.x * p_ig->dxyz.y);

	cudaMemcpyToSymbol(c_alutOffset, &alutOffset, sizeof(float2));
	cudaMemcpyToSymbol(c_voxBase, &voxBase, sizeof(float));
	cudaMemcpyToSymbol(c_nlines, &_nlines, sizeof(int));
}

void cbctProjector::_bindAreaLut_pos(const lut_area* alut)
{
	_nlines = p_cg->nxy.x + 1;

	_mm.allocArray<float>("alut", alut->dim);
	_mm.copyToArray("alut", alut->lut);

	_mm.allocLinear<float4>("lines", _nlines);

	areaTex.normalized = false;
	areaTex.filterMode = cudaFilterModeLinear;
	areaTex.addressMode[0] = cudaAddressModeClamp;
	areaTex.addressMode[1] = cudaAddressModeClamp;

	cudaBindTextureToArray(areaTex, _mm.getArray("alut"));
	cudaBindTexture(0, lineTex, _mm.getLinear<float4>("lines"), lineTex.channelDesc);

	float2 alutOffset = make_float2(1.0f / alut->delta.x, 1.0f / alut->delta.y);
	float  voxBase = fabsf(p_ig->dxyz.x * p_ig->dxyz.y);

	cudaMemcpyToSymbol(c_alutOffset, &alutOffset, sizeof(float2));
	cudaMemcpyToSymbol(c_voxBase, &voxBase, sizeof(float));
	cudaMemcpyToSymbol(c_nlines, &_nlines, sizeof(int));
}

void cbctProjector::_bindHeightLut(const lut_height* hlut)
{
	_nplanes = p_cg->nxy.y + 1 + 2;

	//printf("\n\
	//	Height lookup table\n\
	//	- Dimension: %d x %d x %d\n\
	//	- Delta: %f x %f x %f\n", 
	//	hlut->dim.x, hlut->dim.y, hlut->dim.z,
	//	hlut->delta.x, hlut->delta.y, hlut->delta.z);

	if (hlut)
	{
		_mm.alloc3DArray<float>("hlut", hlut->dim);
		_mm.copyTo3DArray("hlut", hlut->lut);

		heightTex.normalized = false;
		heightTex.filterMode = cudaFilterModeLinear;
		heightTex.addressMode[0] = cudaAddressModeClamp;
		heightTex.addressMode[1] = cudaAddressModeClamp;
		heightTex.addressMode[2] = cudaAddressModeClamp;

		cudaBindTextureToArray(heightTex, _mm.getArray("hlut"));

		float3 hlutOffset = make_float3(1.0f / hlut->delta.x, 1.0f / hlut->delta.y, 1.0f / hlut->delta.z);
		cudaMemcpyToSymbol(c_hlutOffset, &hlutOffset, sizeof(float3));


		_mm.allocLinear<float4>("planes", _nplanes);
		_mm.allocLinear<float>("polar", _nplanes);

		cudaBindTexture(0, planeTex, _mm.getLinear<float4>("planes"), planeTex.channelDesc);
		cudaBindTexture(0, polarTex, _mm.getLinear<float>("polar"), polarTex.channelDesc);

		// Pre-compute polar angles and copy to device memory & bind it to texture
		std::vector<float> polar(_nplanes);
		{
			float t0 = (0 - p_cg->wxy.y - 1.5f)*p_cg->dxy.y;
			for (int i = 0; i < _nplanes; i++, t0 += p_cg->dxy.y)
				polar[i] = atanf(t0 / p_cg->dsd) * (float)RAD2DEG;
		}
		_mm.copyToLinear("polar", &polar[0]);
		//cudaMemcpyToSymbol(c_polar, &polar[0], _nplanes * sizeof(float));

		// constants
		cudaMemcpyToSymbol(c_nplanes, &_nplanes, sizeof(int));
	}
}

__global__ void check_hlut(float* inout, int N)
{
	const int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i >= N)
		return;

	const float pos = inout[i];

	const float h = tex3D(heightTex, fabsf(pos)*c_hlutOffset.x + 0.5f, 0.5f, 0.5f);

	inout[i] = h;
}

void cbctProjector::_bindHeightLut_pos(const lut_height* hlut)
{
	_nplanes = p_cg->nxy.y + 1;

	if (hlut)
	{
		_mm.alloc3DArray<float>("hlut", hlut->dim);
		_mm.copyTo3DArray("hlut", hlut->lut);

		heightTex.normalized = false;
		heightTex.filterMode = cudaFilterModeLinear;
		heightTex.addressMode[0] = cudaAddressModeClamp;
		heightTex.addressMode[1] = cudaAddressModeClamp;
		heightTex.addressMode[2] = cudaAddressModeClamp;

		cudaBindTextureToArray(heightTex, _mm.getArray("hlut"));


		float3 hlutOffset = make_float3(1.0f / hlut->delta.x, 1.0f / hlut->delta.y, 1.0f / hlut->delta.z);
		cudaMemcpyToSymbol(c_hlutOffset, &hlutOffset, sizeof(float3));

		//if (1)
		//{
		//	printf("hlut dim: %d x %d x %d\n", hlut->dim.x, hlut->dim.y, hlut->dim.z);
		//	printf("hlut delta: %f x %f x %f\n", hlut->delta.x, hlut->delta.y, hlut->delta.z);
		//	std::vector<float> posSet(hlut->dim.x+10, 0);
		//	for (int i = 0; i < hlut->dim.x+10; ++i)
		//	{
		//		posSet[i] = i*hlut->delta.x;
		//	}

		//	_mm.allocLinear<float>("posSet", hlut->dim.x+10);
		//	_mm.copyToLinear("posSet", &posSet[0]);
		//	check_hlut<<<iDivUp(hlut->dim.x+10, 256), 256>>>(_mm.getLinear<float>("posSet"), hlut->dim.x+10);

		//	_mm.copyLinearToHost("posSet", &posSet[0]);
		//	mexMat::write("dbg_hlut.mat", "h", &posSet[0], hlut->dim.x+10, 1, 1);
		//}



		_mm.allocLinear<float4>("planes", _nplanes);
		_mm.allocLinear<float>("polar", _nplanes);

		cudaBindTexture(0, planeTex, _mm.getLinear<float4>("planes"), planeTex.channelDesc);
		cudaBindTexture(0, polarTex, _mm.getLinear<float>("polar"), polarTex.channelDesc);

		// Pre-compute polar angles and copy to device memory & bind it to texture
		std::vector<float> polar(_nplanes);
		{
			float t0 = (0 - p_cg->wxy.y - .5f)*p_cg->dxy.y;
			for (int i = 0; i < _nplanes; i++, t0 += p_cg->dxy.y)
				polar[i] = atanf(t0 / p_cg->dsd) * (float)RAD2DEG;
		}
		_mm.copyToLinear("polar", &polar[0]);
		//cudaMemcpyToSymbol(c_polar, &polar[0], _nplanes * sizeof(float));

		// constants
		cudaMemcpyToSymbol(c_nplanes, &_nplanes, sizeof(int));
	}
}

void cbctProjector::_updateRays(
	const float2 uv_s, const float2 uv_t, 
	const float2 src, const float2 det, 
	float ang)
{
	cudaMemcpyToSymbol(c_det, &det, sizeof(float2));
	cudaMemcpyToSymbol(c_src, &src, sizeof(float2));
	cudaMemcpyToSymbol(c_uv_s, &uv_s, sizeof(float2));
	cudaMemcpyToSymbol(c_uv_t, &uv_t, sizeof(float2));

	ang = fmodf(ang, 360.0f);
	if (ang < 0.0f)
		ang += 360.0f;

	cudaMemcpyToSymbol(c_viewang, &ang, sizeof(float));

	if (_method == LUT_LUT || _method == LUT_REG || _method == LUT_DD)
	{
		//_updateLines();
		_updateLines_pos();

		if (_method == LUT_REG || _method == LUT_LUT)
		{
			// why?
			const float2 sn_cs_view = make_float2(
				sinf((ang + 180.0f)*float(DEG2RAD)),
				cosf((ang + 180.0f)*float(DEG2RAD))
			);
			cudaMemcpyToSymbol(c_sn_cs_view, &sn_cs_view, sizeof(float2));
		}

		if (_method == LUT_LUT)
		{
			_updatePlanes_pos();

			// do I need to add 180??
			ang = fmodf(ang+180.0f, 360.0f);
			if (ang < 0.0f)
				ang += 360.0f;

			if (ang <= 45.0f) { ; }
			else if (ang <= 90.0f) { ang = 90.0f - ang; }
			else if (ang <= 135.0f) { ang = ang - 90.0f; }
			else if (ang <= 180.0f) { ang = 180.0f - ang; }
			else if (ang <= 225.0f) { ang = ang - 180.0f; }
			else if (ang <= 270.0f) { ang = 270.0f - ang; }
			else if (ang <= 315.0f) { ang = ang - 270.0f; }
			else { ang = 360.0f - ang; }
			cudaMemcpyToSymbol(c_viewang, &ang, sizeof(float));
		}


	}
}

void cbctProjector::_updateAmp()
{
	dim3 blk(16, 16);
	dim3 grd(iDivUp(p_cg->nxy.x, blk.x), iDivUp(p_cg->nxy.y, blk.y));

	int4 dim = _mm.getDim("amp");

	update_amp1_kernel<<<grd, blk>>>(_mm.getLinear<float>("amp"), dim.w);

	// todo: use cuMemMgr
	cudaMemcpy2DToArray(
		_mm.getArray("amparr"), 0, 0,
		_mm.getLinear<float>("amp"), dim.w*sizeof(float), dim.x*sizeof(float), dim.y,
		cudaMemcpyDeviceToDevice
	);
}

void cbctProjector::_updateLines()
{
	dim3 blk(256, 1, 1);
	dim3 grd(iDivUp(_nlines, blk.x), 1, 1);

	update_lines_kernel<<<grd, blk>>>(_mm.getLinear<float4>("lines"));
}

void cbctProjector::_updateLines_pos()
{
	dim3 blk(256, 1, 1);
	dim3 grd(iDivUp(_nlines, blk.x), 1, 1);

	update_lines_pos_kernel<<<grd, blk>>>(_mm.getLinear<float4>("lines"));

	//std::vector<float4> lines(_nlines);
	//_mm.copyLinearToHost("lines", &lines[0]);
	//mexMat::write("lines.mat", "lines", (float*)&lines[0], 4, _nlines, 1);
}

void cbctProjector::_updatePlanes()
{
	dim3 blk(256, 1, 1);
	dim3 grd(iDivUp(_nplanes, blk.x), 1, 1);

	update_plane_kernel<<<grd,blk>>>(_mm.getLinear<float4>("planes"));

	//cudaMemcpyToSymbol(c_plane, _mm.getLinear<float4>("planes"), _nplanes * sizeof(float4));
}

void cbctProjector::_updatePlanes_pos()
{
	dim3 blk(256, 1, 1);
	dim3 grd(iDivUp(_nplanes, blk.x), 1, 1);

	update_plane_pos_kernel<<<grd, blk>>>(_mm.getLinear<float4>("planes"));

	//std::vector<float4> planes(_nplanes);
	//_mm.copyLinearToHost("planes", &planes[0]);
	//mexMat::write("planes.mat", "planes", (float*)&planes[0], 4, _nplanes, 1);
}

void cbctProjector::_init_footprint()
{
	dim3 blk(16, 16, 1);
	dim3 grd(iDivUp(p_ig->nxyz.x, blk.x), iDivUp(p_ig->nxyz.y, blk.y));

	//todo: switch for each case according to nfoot
	if (_method == LUT_LUT)
	{
		switch (_nfoot)
		{
		case 5: lut_init_foot_kernel<5, 0><<<grd, blk>>>(_mm.getLinear<float>("footinfo")); break;
		case 7: lut_init_foot_kernel<7, 0><<<grd, blk>>>(_mm.getLinear<float>("footinfo")); break;
		}
		
	}
	else if (_method==LUT_REG)
	{
		switch (_nfoot)
		{
		case 5: lut_init_foot_kernel<5, 1><<<grd, blk>>>(_mm.getLinear<float>("footinfo")); break;
		case 7: lut_init_foot_kernel<7, 1><<<grd, blk>>>(_mm.getLinear<float>("footinfo")); break;
		}
		
	}
	else if (_method == LUT_DD)
	{
		switch (_nfoot)
		{
		case 5: lut_init_foot_kernel<5, 2><<<grd, blk>>>(_mm.getLinear<float>("footinfo")); break;
		case 7: lut_init_foot_kernel<7, 2><<<grd, blk>>>(_mm.getLinear<float>("footinfo")); break;
		}
		
	}
	else if (_method==SF_TR)
	{
		switch (_nfoot)
		{
		case 5: sf_init_foot_kernel<5, true><<<grd, blk>>>(_mm.getLinear<float>("footinfo")); break;
		case 7: sf_init_foot_kernel<7, true><<<grd, blk>>>(_mm.getLinear<float>("footinfo")); break;
		}
	}
	else
	{
		switch (_nfoot)
		{
		case 5: sf_init_foot_kernel<5, false><<<grd, blk>>>(_mm.getLinear<float>("footinfo")); break;
		case 7: sf_init_foot_kernel<7, false><<<grd, blk>>>(_mm.getLinear<float>("footinfo")); break;
		}
		
	}

	
}

void cbctProjector::_init_footprint_pos()
{
	dim3 blk(16, 16, 1);
	dim3 grd(iDivUp(p_ig->nxyz.x, blk.x), iDivUp(p_ig->nxyz.y, blk.y));

	//todo: switch for each case according to nfoot
	if (_method == LUT_LUT)
	{
		switch (_nfoot)
		{
		case 5: lut_init_foot_pos_kernel<5, 0><<<grd, blk>>>(_mm.getLinear<float>("footinfo")); break;
		case 7: lut_init_foot_pos_kernel<7, 0><<<grd, blk>>>(_mm.getLinear<float>("footinfo")); break;
		}

	}
	else if (_method == LUT_REG)
	{
		switch (_nfoot)
		{
		case 5: lut_init_foot_pos_kernel<5, 1><<<grd, blk>>>(_mm.getLinear<float>("footinfo")); break;
		case 7: lut_init_foot_pos_kernel<7, 1><<<grd, blk>>>(_mm.getLinear<float>("footinfo")); break;
		}

	}
	else if (_method == LUT_DD)
	{
		switch (_nfoot)
		{
		case 5: lut_init_foot_pos_kernel<5, 2><<<grd, blk>>>(_mm.getLinear<float>("footinfo")); break;
		case 7: lut_init_foot_pos_kernel<7, 2><<<grd, blk>>>(_mm.getLinear<float>("footinfo")); break;
		}

	}
	else if (_method == SF_TR)
	{
		switch (_nfoot)
		{
		case 5: sf_init_foot_kernel<5, true><<<grd, blk>>>(_mm.getLinear<float>("footinfo")); break;
		case 7: sf_init_foot_kernel<7, true><<<grd, blk>>>(_mm.getLinear<float>("footinfo")); break;
		}
	}
	else
	{
		switch (_nfoot)
		{
		case 5: sf_init_foot_kernel<5, false><<<grd, blk>>>(_mm.getLinear<float>("footinfo")); break;
		case 7: sf_init_foot_kernel<7, false><<<grd, blk>>>(_mm.getLinear<float>("footinfo")); break;
		}

	}
}

void cbctProjector::_do_projection(float* d_proj, const int pitch, const float val)
{
	dim3 blk(p_ig->nxyz.z <= 32 ? 32: p_ig->nxyz.z , 1, 1);
	dim3 grd(iDivUp(p_ig->nxyz.z,blk.x), p_ig->nxyz.x, p_ig->nxyz.y);

	//todo: add more cases based on nfoot
	if (_method == LUT_LUT)
	{
		switch (_nfoot)
		{
		case 5: lut_fp_ll_kernel<5><<<grd, blk>>>(d_proj, pitch, _mm.getLinear<float>("footinfo"), val); break;
		case 7: lut_fp_ll_kernel<7><<<grd, blk>>>(d_proj, pitch, _mm.getLinear<float>("footinfo"), val); break;
		default:
			;
		}
	}
	else if (_method == LUT_REG)
	{
		switch (_nfoot)
		{
		case 5:	lut_fp_lr_kernel<5><<<grd, blk>>>(d_proj, pitch, _mm.getLinear<float>("footinfo"), val); break;
		case 7:	lut_fp_lr_kernel<7><<<grd, blk>>>(d_proj, pitch, _mm.getLinear<float>("footinfo"), val); break;
		default:
			;
		}
	}
	else if (_method == LUT_DD)
	{
		switch (_nfoot)
		{
		case 5: lut_fp_ld_kernel<5><<<grd, blk>>>(d_proj, pitch, _mm.getLinear<float>("footinfo"), val); break;
		case 7: lut_fp_ld_kernel<7><<<grd, blk>>>(d_proj, pitch, _mm.getLinear<float>("footinfo"), val); break;
		default:
			;
		}
	}
	else if (_method == SF_TR)
	{
		switch (_nfoot)
		{
		case 5: sf_fp_tr_kernel<5><<<grd, blk>>>(d_proj, pitch, _mm.getLinear<float>("footinfo"), val); break;
		case 7: sf_fp_tr_kernel<7><<<grd, blk>>>(d_proj, pitch, _mm.getLinear<float>("footinfo"), val); break;
		default:
			;
		}
	}
	else
	{
		switch (_nfoot)
		{
		case 5: sf_fp_tt_kernel<5><<<grd, blk>>>(d_proj, pitch, _mm.getLinear<float>("footinfo"), val); break;
		case 7: sf_fp_tt_kernel<7><<<grd, blk>>>(d_proj, pitch, _mm.getLinear<float>("footinfo"), val); break;
		default:
			;
		}
	}
}

void cbctProjector::_do_backprojection(float* d_image, const int pitch, const int pid, const float val)
{
	dim3 blk;
	dim3 grd;

	blk.x = p_ig->nxyz.z; 
	blk.y = 1; blk.z = 1;
	grd.x = iDivUp(p_ig->nxyz.z, blk.x); grd.y = p_ig->nxyz.x; grd.z = p_ig->nxyz.y;

	//todo: add more cases based on nfoot
	if (val <= 0.0f)
	{
		if (_method == LUT_LUT)
		{
			lut_bp_ll_kernel<5, false><<<grd, blk>>>(d_image, pitch, _mm.getLinear<float>("footinfo"), pid);
		}
		else if (_method == LUT_REG)
		{
			lut_bp_lr_kernel<5, false><<<grd, blk>>>(d_image, pitch, _mm.getLinear<float>("footinfo"), pid);
		}
		else if (_method == LUT_DD)
		{
			lut_bp_ld_kernel<5, false><<<grd, blk>>>(d_image, pitch, _mm.getLinear<float>("footinfo"), pid);
		}
		else if (_method == SF_TR)
		{
			sf_bp_tr_kernel<5, false><<<grd, blk>>>(d_image, pitch, _mm.getLinear<float>("footinfo"), pid);
		}
		else
		{
			sf_bp_tt_kernel<5, false><<<grd, blk>>>(d_image, pitch, _mm.getLinear<float>("footinfo"), pid);
		}
	}
	else
	{
		if (_method == LUT_LUT)
		{
			lut_bp_ll_kernel<5, true><<<grd, blk>>>(d_image, pitch, _mm.getLinear<float>("footinfo"), pid);
		}
		else if (_method == LUT_REG)
		{
			lut_bp_lr_kernel<5, true><<<grd, blk>>>(d_image, pitch, _mm.getLinear<float>("footinfo"), pid);
		}
		else if (_method == LUT_DD)
		{
			lut_bp_ld_kernel<5, true><<<grd, blk>>>(d_image, pitch, _mm.getLinear<float>("footinfo"), pid);
		}
		else if (_method == SF_TR)
		{
			sf_bp_tr_kernel<5, true><<<grd, blk>>>(d_image, pitch, _mm.getLinear<float>("footinfo"), pid);
		}
		else
		{
			sf_bp_tt_kernel<5, true><<<grd, blk>>>(d_image, pitch, _mm.getLinear<float>("footinfo"), pid);
		}

	}
}

void cbctProjector::_scale_projection(float* d_proj, const int pitch) const
{
	dim3 blk(16, 16);
	dim3 grd(iDivUp(p_cg->nxy.x, blk.x), iDivUp(p_cg->nxy.y, blk.y));

	if (_method == LUT_LUT || _method == LUT_REG || _method == LUT_DD)
	{
		apply_geodiv_kernel<<<grd, blk>>>(d_proj, pitch);
	}
	else
	{
		apply_amp1_kernel<<<grd,blk>>>(d_proj, pitch);
	}
}

void cbctProjector::_scale_projection(float* d_proj, const int pitch, const float beta) const
{
	dim3 blk(16, 16);
	dim3 grd(iDivUp(p_cg->nxy.x, blk.x), iDivUp(p_cg->nxy.y, blk.y));

	if (_method == LUT_LUT || _method == LUT_REG || _method == LUT_DD)
	{
		//if (_order == TS) 
			apply_geodiv_kernel<<<grd, blk>>>(d_proj, pitch);
		//else
		//	apply_geodiv_kernel_st<<<grd, blk>>>(d_proj, pitch);
	}
	else
	{
		apply_amp1_beta_kernel<<<grd,blk>>>(d_proj, pitch, beta);
	}
}



