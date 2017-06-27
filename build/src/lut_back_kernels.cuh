
template<int nfoot, bool allOnes>
__global__ static void
lut_bp_ll_kernel(float* image, const int pitch, const float* footinfo, const int pid)
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


	const float x = *p_x;
	const float y = *p_y;
	const float z = (iz - g_ig.wxyz.z)*g_ig.dxyz.z;
	const float div = 1.0f / (*p_div + SQR(z));
	const float tp = (z + c_dz2)*(*p_mag);
	const float tm = (z - c_dz2)*(*p_mag);

	const int tmax = min((int)ceilf(tp*g_cg.invdxy.y + g_cg.wxy.y + .5f), g_cg.nxy.y - 1);
	const int tmin = max((int)floorf(tm*g_cg.invdxy.y + g_cg.wxy.y - .5f), 0);

	float sum = 0.0f;

	//int tidx = tmin + 1;
	//float h0 = fetchHeightLut(tidx, x, y, z); 
	//++tidx;
	float t = (tmin - g_cg.wxy.y)*g_cg.dxy.y;
	float h0 = fetchHeightLut_pos_t(x, y, z, t - c_dt2);

	for (int it = tmin; it <= tmax; ++it,  /*++tidx*/ t += g_cg.dxy.y)
	{
		//const float h1 = fetchHeightLut(tidx,x, y, z);
		const float h1 = fetchHeightLut_pos_t(x, y, z, t + c_dt2);
		const float effh = fabsf(h0 - h1);
		h0 = h1;

		if (effh <= 0.0f)
			continue;

		int is = int(*p_sbin);
		for (int ifoot = 0; ifoot<nfoot; ++is, ++ifoot)
		{
			const float area = p_foot[ifoot];

			if (allOnes)
			{
				const float density = tex2D(geodivTex, is + 0.5f, it + 0.5f);
				sum += density*div*area*effh;
			}
			else
			{
				const float density = tex2DLayered(texRefLayered, it + 0.5f, is + 0.5f, pid);
				sum += density*div*area*effh;
			}
		}
	}

	atomicAdd(image + (ix + iy*g_ig.nxyz.x)*pitch + iz, sum);
}

template<int nfoot, bool allOnes>
__global__ static void
lut_bp_lr_kernel(float* image, const int pitch, const float* footinfo, const int pid)
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

	const float x = *p_x;
	const float y = *p_y;
	const float z = (iz - g_ig.wxyz.z)*g_ig.dxyz.z;
	const float div = 1.0f / (*p_div + SQR(z));
	//const float tp = (z + c_dz2)*(*p_mag);
	const float tm = (z - c_dz2)*(*p_mag);

	//const int tmax = min((int)ceilf(tp*g_cg.invdxy.y + g_cg.wxy.y + .5f), g_cg.nxy.y - 1);
	const int tmin = max((int)floorf(tm*g_cg.invdxy.y + g_cg.wxy.y - .5f), 0);

	float sum = 0.0f;
	float t = (tmin - g_cg.wxy.y)*g_cg.dxy.y;
	float h0 = cmpHeight(x, y, z, t - c_dt2);
	//for (int it = tmin; it <= tmax; ++it, t+=g_cg.dxy.y)
	for (int ifoot_y = 0; ifoot_y < nfoot; ++ifoot_y, t += g_cg.dxy.y)
	{
		const float h1 = cmpHeight(x, y, z, t + c_dt2);
		const float effh = fabsf(h0 - h1);
		h0 = h1;

		if (effh <= 0.0f)
			continue;

		int is = int(*p_sbin);
		for (int ifoot = 0; ifoot<nfoot; ++is, ++ifoot)
		{
			const float area = p_foot[ifoot];

			if (allOnes)
			{
				const float density = tex2D(geodivTex, is + 0.5f, ifoot_y + tmin + 0.5f);
				sum += density*div*area*effh;
			}
			else
			{
				const float density = tex2DLayered(texRefLayered, ifoot_y + tmin + 0.5f, is + 0.5f, pid);
				sum += density*div*area*effh;
			}
		}
	}

	atomicAdd(image + (ix + iy*g_ig.nxyz.x)*pitch + iz, sum);
}

template<int nfoot, bool allOnes>
__global__ static void
lut_bp_ld_kernel(float* image, const int pitch, const float* footinfo, const int pid)
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

	if (threadIdx.x < nfoot + 4)
	{
		const int stride = g_ig.nxyz.x*g_ig.nxyz.y;
		footinfo += iy*g_ig.nxyz.x + ix;

		shdata[threadIdx.x] = footinfo[threadIdx.x*stride];
	}
	__syncthreads();

	if (ix >= g_ig.nxyz.x || iy >= g_ig.nxyz.y || iz >= g_ig.nxyz.z)
		return;

	const float z = (iz - g_ig.wxyz.z)*g_ig.dxyz.z;
	const float div = 1.0f / (*p_div + SQR(z));
	const float tp = (z + c_dz2)*(*p_mag);
	const float tm = (z - c_dz2)*(*p_mag);
	const float invmag = *p_invmag;

	//const int tmax = min((int)ceilf(tp*g_cg.invdxy.y + g_cg.wxy.y + .5f), g_cg.nxy.y - 1);
	const int tmin = max((int)floorf(tm*g_cg.invdxy.y + g_cg.wxy.y - .5f), 0);

	float sum = 0.0f;
	float t = (tmin - g_cg.wxy.y)*g_cg.dxy.y;
	//for (int it = tmin; it <= tmax; ++it, t += g_cg.dxy.y)
	for (int ifoot_y = 0; ifoot_y < nfoot; ++ifoot_y, t += g_cg.dxy.y)
	{
		const float ttp = t + c_dt2;
		const float ttm = t - c_dt2;
		const float effh = (ttm >= tp || ttp <= tm) ? 0.0f : invmag*(fminf(tp, ttp) - fmaxf(tm, ttm));

		if (effh <= 0.0f)
			continue;

		int is = int(*p_sbin);
		for (int ifoot = 0; ifoot<nfoot; ++is, ++ifoot)
		{
			const float area = p_foot[ifoot];
			if (allOnes)
			{
				const float density = tex2D(geodivTex, is + 0.5f, ifoot_y + tmin + 0.5f);
				sum += density*div*area*effh;
			}
			else
			{
				const float density = tex2DLayered(texRefLayered, ifoot_y + tmin + 0.5f, is + 0.5f, pid);
				sum += density*div*area*effh;
			}
		}
	}

	atomicAdd(image + (ix + iy*g_ig.nxyz.x)*pitch + iz, sum);
}



