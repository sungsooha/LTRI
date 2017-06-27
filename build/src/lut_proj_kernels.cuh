template<int nfoot, int method>
__global__ static void
lut_init_foot_kernel(float* footinfo)
{
	const int ix = blockIdx.x*blockDim.x + threadIdx.x;
	const int iy = blockIdx.y*blockDim.y + threadIdx.y;

	if (ix >= g_ig.nxyz.x || iy >= g_ig.nxyz.y)
		return;

	const int stride = g_ig.nxyz.x*g_ig.nxyz.y;
	footinfo += (iy*g_ig.nxyz.x + ix);

	const int nfoot2 = nfoot / 2;
	const float2 xy = make_float2(
		(ix - g_ig.wxyz.x)*g_ig.dxyz.x,
		(iy - g_ig.wxyz.y)*g_ig.dxyz.y
	);

	const float div = SQR(xy.x - c_src.x) + SQR(xy.y - c_src.y);
	*footinfo = div; footinfo += stride;

	const float mag = g_cg.dsd / (dot(c_uv_t, xy) + g_cg.dso);
	*footinfo = mag;
	footinfo +=stride;

	const float s = mag*dot(c_uv_s, xy);
	const int s_bin = (int)floorf(s*g_cg.invdxy.x + g_cg.wxy.x) - nfoot2;

	*footinfo = float(s_bin);
	footinfo += stride;

	int is = s_bin + 1;
	float area0 = fetchAreaLut(is, xy); ++is;
	for (int ifoot = 0; ifoot < nfoot; ++ifoot, ++is, footinfo+=stride)
	{
		const float area1 = fetchAreaLut(is, xy);
		*footinfo = fabsf(area0 - area1);
		//const float area = fabsf(area1 - area0);
		//*footinfo = area < 0.0001f ? 0.0f: area;
		area0 = area1;
	}

	if (method==0 || method==1) // LL, LR
	{
		*footinfo = xy.x; footinfo += stride;
		*footinfo = xy.y; //footinfo += stride;
	}
	else 
	{
		// LD
		*footinfo = 1.0f / mag;
	}
}

template<int nfoot, int method>
__global__ static void
lut_init_foot_pos_kernel(float* footinfo)
{
	const int ix = blockIdx.x*blockDim.x + threadIdx.x;
	const int iy = blockIdx.y*blockDim.y + threadIdx.y;

	if (ix >= g_ig.nxyz.x || iy >= g_ig.nxyz.y)
		return;

	const int stride = g_ig.nxyz.x*g_ig.nxyz.y;
	footinfo += (iy*g_ig.nxyz.x + ix);

	const int nfoot2 = nfoot / 2;
	const float2 xy = make_float2(
		(ix - g_ig.wxyz.x)*g_ig.dxyz.x,
		(iy - g_ig.wxyz.y)*g_ig.dxyz.y
	);

	const float div = SQR(xy.x - c_src.x) + SQR(xy.y - c_src.y);
	*footinfo = div; footinfo += stride;

	const float mag = g_cg.dsd / (dot(c_uv_t, xy) + g_cg.dso);
	*footinfo = mag;
	footinfo += stride;

	const float s = mag*dot(c_uv_s, xy);
	const int s_bin = (int)floorf(s*g_cg.invdxy.x + g_cg.wxy.x) - nfoot2;

	*footinfo = float(s_bin);
	footinfo += stride;

	int is = s_bin;
	float area0 = fetchAreaLut_pos(is, xy); ++is;
	for (int ifoot = 0; ifoot < nfoot; ++ifoot, ++is, footinfo += stride)
	{
		const float area1 = fetchAreaLut_pos(is, xy);
		*footinfo = (area0 - area1);
		//const float area = fabsf(area1 - area0);
		//*footinfo = area < 0.0001f ? 0.0f: area;
		area0 = area1;
	}

	//if (ix==128 && iy==128)
	//{
	//	is = s_bin;
	//	area0 = fetchAreaLut_pos(is, xy); 
	//	printf("%d: %f\n", is, area0);
	//	++is;
	//	for (int ifoot=0; ifoot<nfoot; ++ifoot, ++is)
	//	{
	//		printf("%d: %f\n", is, fetchAreaLut_pos(is, xy));
	//	}
	//}

	if (method == 0 || method == 1) // LL, LR
	{
		*footinfo = xy.x; footinfo += stride;
		*footinfo = xy.y; //footinfo += stride;
	}
	else
	{
		// LD
		*footinfo = 1.0f / mag;
	}
}

template<int nfoot>
__global__ static void
lut_fp_ll_kernel(float* proj, const int pitch, const float* footinfo, const float val)
{
	//extern __shared__ float shdata[];
	__shared__ float shdata[3 + nfoot + 2];

	const int ix = blockIdx.y;
	const int iy = blockIdx.z;
	const int iz = blockIdx.x*blockDim.x + threadIdx.x;

	const float* p_div = shdata + 0;
	const float* p_mag = shdata + 1;
	const float* p_sbin = shdata + 2;
	const float* p_foot = shdata + 3;
	const float* p_x = shdata + 3 + nfoot;
	const float* p_y = shdata + 3 + nfoot + 1;

	if (threadIdx.x < nfoot + 5)
	{
		const int stride = g_ig.nxyz.x*g_ig.nxyz.y;
		footinfo += iy*g_ig.nxyz.x + ix;
		shdata[threadIdx.x] = footinfo[threadIdx.x*stride];
	}
	__syncthreads();

	if (ix >= g_ig.nxyz.x || iy >= g_ig.nxyz.y || iz >= g_ig.nxyz.z)
		return;

	const float att = (val>0.0f) ? val: tex3D(tex3DRef, iz + 0.5f, ix + 0.5f, iy + 0.5f); // fetching

	if (att == 0.0f)
		return;

	const float x = *p_x;
	const float y = *p_y;
	const float z = (iz - g_ig.wxyz.z)*g_ig.dxyz.z;
	const float div = att / (*p_div + SQR(z));
	const float tp = (z + c_dz2)*(*p_mag);
	const float tm = (z - c_dz2)*(*p_mag);

	const int tmax = min((int)ceilf(tp*g_cg.invdxy.y + g_cg.wxy.y + .5f), g_cg.nxy.y - 1);
	const int tmin = max((int)floorf(tm*g_cg.invdxy.y + g_cg.wxy.y - .5f), 0);

	//int tidx = tmin;//tmin +1;
	//float h0 = fetchHeightLut_pos(tidx, x, y, z); 
	//++tidx;
	float t = (tmin - g_cg.wxy.y)*g_cg.dxy.y;
	float h0 = fetchHeightLut_pos_t(x, y, z, t - c_dt2);
	for (int it = tmin; it <= tmax; ++it, /*++tidx*/ t+=g_cg.dxy.y)
	{
		//const float h1 = fetchHeightLut_pos(tidx, x, y, z);
		const float h1 = fetchHeightLut_pos_t(x, y, z, t + c_dt2);
		const float effh = fabsf(h0 - h1);//fabsf(h1 - h0);
		h0 = h1;

		if (effh <= 0.0f)
			continue;

		int is = int(*p_sbin);
		for (int ifoot = 0; ifoot<nfoot; ++is, ++ifoot)
		{
			const float area = p_foot[ifoot];

			if (is < 0 || is >= g_cg.nxy.x || area <= 0.0f)
				continue;

			atomicAdd(proj + is*pitch + it, div*area*effh);
		}
	}
}


template<int nfoot>
__global__ static void
lut_fp_lr_kernel(float* proj, const int pitch, const float* footinfo, const float val)
{
	//extern __shared__ float shdata[];
	__shared__ float shdata[3 +nfoot + 2];

	const int ix = blockIdx.y;
	const int iy = blockIdx.z;
	const int iz = blockIdx.x*blockDim.x + threadIdx.x;

	const float* p_div = shdata + 0;
	const float* p_mag = shdata + 1;
	const float* p_sbin = shdata + 2;
	const float* p_foot = shdata + 3;
	const float* p_x = shdata + 3 + nfoot;
	const float* p_y = shdata + 3 + nfoot + 1;

	if (threadIdx.x < nfoot + 5)
	{
		const int stride = g_ig.nxyz.x*g_ig.nxyz.y;
		footinfo += iy*g_ig.nxyz.x + ix;
		shdata[threadIdx.x] = footinfo[threadIdx.x*stride];
	}
	__syncthreads();

	if (ix >= g_ig.nxyz.x || iy >= g_ig.nxyz.y || iz >= g_ig.nxyz.z)
		return;

	const float att = (val>0.0f) ? val : tex3D(tex3DRef, iz + 0.5f, ix + 0.5f, iy + 0.5f); // fetching

	if (att == 0.0f)
		return;

	const float x = *p_x;
	const float y = *p_y;
	const float z = (iz - g_ig.wxyz.z)*g_ig.dxyz.z;
	const float div = att / (*p_div + SQR(z));
	const float tp = (z + c_dz2)*(*p_mag);
	const float tm = (z - c_dz2)*(*p_mag);

	const int tmax = min((int)ceilf(tp*g_cg.invdxy.y + g_cg.wxy.y + .5f), g_cg.nxy.y - 1);
	const int tmin = max((int)floorf(tm*g_cg.invdxy.y + g_cg.wxy.y - .5f), 0);

	float t = (tmin - g_cg.wxy.y)*g_cg.dxy.y;
	float h0 = cmpHeight(x, y, z, t-c_dt2); 
	for (int it = tmin; it <= tmax; ++it, t+=g_cg.dxy.y)
	{
		const float h1 = cmpHeight(x, y, z, t+c_dt2);
		const float effh = fabsf(h0 - h1);
		h0 = h1;

		if (effh <= 0.0f)
			continue;

		int is = int(*p_sbin);
		for (int ifoot = 0; ifoot<nfoot; ++is, ++ifoot)
		{
			const float area =  p_foot[ifoot];

			if (is < 0 || is >= g_cg.nxy.x || area <= 0.0f || att == 0.0f)
				continue;

			atomicAdd(proj + is*pitch + it, div*area*effh);
		}
	}
}


template<int nfoot>
__global__ static void
lut_fp_ld_kernel(float* proj, const int pitch, const float* footinfo, const float val)
{
	//extern __shared__ float shdata[];
	__shared__ float shdata[3 + nfoot + 1];

	const int ix = blockIdx.y;
	const int iy = blockIdx.z;
	const int iz = blockIdx.x*blockDim.x + threadIdx.x;

	const float* p_div = shdata + 0;
	const float* p_mag = shdata + 1;
	const float* p_sbin = shdata + 2;
	const float* p_foot = shdata + 3;
	const float* p_invmag = shdata + 3 + nfoot;

	if (threadIdx.x < nfoot+4)
	{
		const int stride = g_ig.nxyz.x*g_ig.nxyz.y;
		footinfo += iy*g_ig.nxyz.x + ix;

		shdata[threadIdx.x] = footinfo[threadIdx.x*stride];
	}
	__syncthreads();

	if (ix >= g_ig.nxyz.x || iy >= g_ig.nxyz.y || iz >= g_ig.nxyz.z)
		return;

	const float att = (val>0.0f) ? val : tex3D(tex3DRef, iz + 0.5f, ix + 0.5f, iy + 0.5f); // fetching

	if (att == 0.0f)
		return;

	const float z = (iz - g_ig.wxyz.z)*g_ig.dxyz.z;
	const float div = att / (*p_div + SQR(z));
	const float tp = (z + c_dz2)*(*p_mag);
	const float tm = (z - c_dz2)*(*p_mag);
	const float invmag = *p_invmag;

	//const int tmax = min((int)ceilf(tp*g_cg.invdxy.y + g_cg.wxy.y + .5f), g_cg.nxy.y - 1);
	const int tmin = max((int)floorf(tm*g_cg.invdxy.y + g_cg.wxy.y - .5f), 0);

	float t = (tmin - g_cg.wxy.y)*g_cg.dxy.y;
	//for (int it = tmin; it <= tmax; ++it, t += g_cg.dxy.y)
	for (int ifoot_y = 0; ifoot_y < nfoot; ++ifoot_y, t += g_cg.dxy.y)
	{
		const float ttp = t + c_dt2;
		const float ttm = t - c_dt2;
		const float effh = (ttm > tp || ttp < tm) ? 0.0f : invmag*(fminf(tp, ttp) - fmaxf(tm, ttm));

		if (effh <= 0.0f)
			continue;

		int is = int(*p_sbin);
		for (int ifoot = 0; ifoot<nfoot; ++is, ++ifoot)
		{
			const float area = p_foot[ifoot];

			if (is < 0 || is >= g_cg.nxy.x || area <= 0.0f || div == 0.0f)
				continue;

			atomicAdd(proj + is*pitch + ifoot_y + tmin, div*area*effh);
		}
	}
}
