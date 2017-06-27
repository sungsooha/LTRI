__global__ static void
update_lines_kernel(float4* lines)
{
	const int is = blockIdx.x*blockDim.x + threadIdx.x;
	if (is >= c_nlines)
		return;

	const float s0 = (0 - g_cg.wxy.x - 1.5f)*g_cg.dxy.x;
	const float s = s0 + is*g_cg.dxy.x;

	const float2 Ps = c_det + s*c_uv_s;
	const float2 rayvec = Ps - c_src;

	float ang = atan2f(rayvec.y, rayvec.x) * (360.0f / (2.0f*(float)M_PI));
	if (ang < 0.0f)
		ang += 360.0f;

	const float A = Ps.y - c_src.y;
	const float B = c_src.x - Ps.x;
	const float C = Ps.x*c_src.y - c_src.x*Ps.y;
	const float Z = sqrtf(A*A + B*B);

	if (ang <= 45.0f) { ; }
	else if (ang <= 90.0f) { ang = 90.0f - ang; }
	else if (ang <= 135.0f) { ang = ang - 90.0f; }
	else if (ang <= 180.0f) { ang = 180.0f - ang; }
	else if (ang <= 225.0f) { ang = ang - 180.0f; }
	else if (ang <= 270.0f) { ang = 270.0f - ang; }
	else if (ang <= 315.0f) { ang = ang - 270.0f; }
	else { ang = 360.0f - ang; }

	lines[is] = make_float4(ang, A / Z, B / Z, C / Z);
}

__global__ static void
update_plane_kernel(float4* planes)
{
	const int it = blockIdx.x*blockDim.x + threadIdx.x;
	if (it >= c_nplanes)
		return;

	const float t0 = (0 - g_cg.wxy.y - 1.5f)*g_cg.dxy.y;
	const float t = t0 + it*g_cg.dxy.y;
	const float ang = atanf(t / g_cg.dsd);                      // always const.

	const float cs = cosf(ang);
	const float sn = sinf(ang);

	const float csb = cosf((c_viewang+180.0f)*float(DEG2RAD));
	const float snb = sinf((c_viewang+180.0f)*float(DEG2RAD));

	const float3 normal = make_float3(snb*sn, -csb*sn, cs);
	const float dist = -g_cg.dso*sn;

	planes[it] = make_float4(normal, dist);
}

__global__ static void
update_plane_pos_kernel(float4* planes)
{
	const int it = blockIdx.x*blockDim.x + threadIdx.x;
	if (it >= c_nplanes)
		return;

	const float t0 = (0 - g_cg.wxy.y - .5f)*g_cg.dxy.y;
	const float t = t0 + it*g_cg.dxy.y;
	const float ang = atanf(t / g_cg.dsd);                      // always const.

	const float cs = cosf(ang);
	const float sn = sinf(ang);

	const float csb = cosf((c_viewang + 180.0f)*float(DEG2RAD));
	const float snb = sinf((c_viewang + 180.0f)*float(DEG2RAD));

	const float3 normal = make_float3(snb*sn, -csb*sn, cs);
	const float dist = dot(normal, make_float3(c_src));

	planes[it] = make_float4(normal, dist);
}

__device__ __inline__ static float
fetchAreaLut(int sidx, const float2& xy)
{
	sidx = clamp(sidx, 0, c_nlines - 1);
	const float4 L = tex1Dfetch(lineTex, sidx);
	const float
		ang = L.x,
		pos = L.y*xy.x + L.z*xy.y + L.w;

#ifdef USEREG
	float value;
	if (ang > 0.001f)
	{
#else
	const float
#endif
		value = tex2D(
			areaTex,
			fabsf(pos)*c_alutOffset.x + 0.5f,
			ang*c_alutOffset.y + 0.5f
		);
#ifdef USEREG
	}
	else
	{
		const float dx2 = 0.5f*g_ig.dxyz.x;
		const float dx = 0.5f*g_ig.dxyz.x;
		value = fabsf(pos) < dx2 ? dx*(dx2 - fabsf(pos)) : 0.0f;
	}
#endif
	return pos < 0.0f ? c_voxBase - value : value;
}

__device__ __inline__ static float
fetchAreaLut_pos(int sidx, const float2& xy)
{
	sidx = clamp(sidx, 0, c_nlines - 1);
	const float4 L = tex1Dfetch(lineTex, sidx);

	//const float s = (sidx - g_cg.wxy.x)*g_cg.dxy.x;
	//const float2 Ps = c_det + s*c_uv_s;
	//const float2 rayvec = Ps - c_src;

	//float ang = atan2f(rayvec.y, rayvec.x) * (360.0f / (2.0f*(float)M_PI));
	//if (ang < 0.0f)
	//	ang += 360.0f;

	//const float A = Ps.y - c_src.y;
	//const float B = c_src.x - Ps.x;
	//const float C = c_src.x*Ps.y - Ps.x*c_src.y;
	//const float4 L = make_float4(ang, A, B, C);
	
	float ang = L.x;
	float pos = L.w - L.y*xy.x - L.z*xy.y;
	//float pos = L.w + L.y*xy.x + L.z*xy.y;
	float value;

	if (ang <= 45.0f)
	{
		pos = pos / L.z;
		value = tex2D(areaTex, fabsf(pos)*c_alutOffset.x + 0.5f, ang*c_alutOffset.y + 0.5f);
	}
	else if (ang <= 90.0f)
	{
		pos = pos / L.y; 
		ang = 90.0f - ang;
		value = c_voxBase - tex2D(areaTex, fabsf(pos)*c_alutOffset.x + 0.5f, ang*c_alutOffset.y + 0.5f);
	}
	else if (ang <= 135.0f)
	{
		pos = pos / L.y; 
		ang = ang - 90.0f;
		value = c_voxBase - tex2D(areaTex, fabsf(pos)*c_alutOffset.x + 0.5f, ang*c_alutOffset.y + 0.5f);
	}
	else if (ang <= 180.0f)
	{
		pos = pos / L.z; 
		ang = 180.0f - ang;
		value = c_voxBase - tex2D(areaTex, fabsf(pos)*c_alutOffset.x + 0.5f, ang*c_alutOffset.y + 0.5f);
	}
	else if (ang <= 225.0f)
	{
		pos = pos / L.z; 
		ang = ang - 180.0f;
		value = c_voxBase - tex2D(areaTex, fabsf(pos)*c_alutOffset.x + 0.5f, ang*c_alutOffset.y + 0.5f);
	}
	else if (ang <= 270.0f)
	{
		pos = pos / L.y; 
		ang = 270.0f - ang;
		value = tex2D(areaTex, fabsf(pos)*c_alutOffset.x + 0.5f, ang*c_alutOffset.y + 0.5f);
	}
	else if (ang <= 315.0f)
	{
		pos = pos / L.y; 
		ang = ang - 270.0f;
		value = tex2D(areaTex, fabsf(pos)*c_alutOffset.x + 0.5f, ang*c_alutOffset.y + 0.5f);
	}
	else
	{
		pos = pos / L.z; 
		ang = 360.0f - ang;
		value = tex2D(areaTex, fabsf(pos)*c_alutOffset.x + 0.5f, ang*c_alutOffset.y + 0.5f);
	}

	return pos < 0.0f ? c_voxBase - value : value;
}

__device__ __inline__ static float
fetchHeightLut(int tidx, const float x, const float y, const float z)
{
	const float4 Pl = tex1Dfetch(planeTex, tidx);
	const float polar = tex1Dfetch(polarTex, tidx);
	const float dist = -(Pl.x*x + Pl.y*y + Pl.z*z + Pl.w);

#ifdef USEREG
	float h;
	if (fabsf(polar) > 0.001f)
	{
#else
	const float
#endif
		h = tex3D(heightTex,
			fabsf(dist)*c_hlutOffset.x + 0.5f,
			fabsf(polar)*c_hlutOffset.y + 0.5f,
			c_viewang*c_hlutOffset.z + 0.5f
		);
#ifdef USEREG
	}
	else
	{
		const float dz2 = 0.5f*g_ig.dxyz.z;
		h = (fabsf(dist) < dz2) ? dz2 - fabsf(dist) : 0.0f;
	}
#endif

	
	return (dist > 0.0f) ? g_ig.dxyz.z - h : h;
}

__device__ __inline__ static float
fetchHeightLut_pos(int tidx, const float x, const float y, const float z)
{
	const float4 Pl = tex1Dfetch(planeTex, tidx);
	const float polar = tex1Dfetch(polarTex, tidx);
	const float pos = ( Pl.w - (Pl.x*x + Pl.y*y + Pl.z*z) )/Pl.z;

#ifdef USEREG
	float h;
	if (fabsf(polar) > 10000.1f)
	{
#else
	const float
#endif
		h = tex3D(heightTex,
			fabsf(pos)*c_hlutOffset.x + 0.5f,
			fabsf(polar)*c_hlutOffset.y + 0.5f,
			c_viewang*c_hlutOffset.z + 0.5f
		);

#ifdef USEREG
	}
	else
	{
		h = (fabsf(pos) < c_dz2) ? c_dz2 - fabsf(pos) : 0.0f;
	}
#endif

	return (pos < 0.0f) ? g_ig.dxyz.z - h : h;
}

__device__ __inline__ static float
fetchHeightLut_pos_t(const float x, const float y, const float z, const float t)
{
	const float L = sqrtf(t*t + SQR(g_cg.dsd));
	const float cs = g_cg.dsd / L;
	const float sn = t / L;

	const float3 normal = make_float3(c_sn_cs_view.x*sn, -c_sn_cs_view.y*sn, cs);
	const float w = dot(normal, make_float3(c_src));
	const float pos = (w - (normal.x*x + normal.y*y + normal.z*z)) / normal.z;

	const float polar = atanf(t/g_cg.dsd) * float(RAD2DEG);

#ifdef USEREG
	float h;
	if (fabsf(polar) > 0.001f)
	{
#else
	const float
#endif


		h = tex3D(heightTex,
			fabsf(pos)*c_hlutOffset.x + 0.5f, 
			fabsf(polar)*c_hlutOffset.y + 0.5f,
			c_viewang*c_hlutOffset.z + 0.5f
		);

#ifdef USEREG
	}
	else
	{
		h = (fabsf(pos) < c_dz2) ? c_dz2 - fabsf(pos) : 0.0f;
	}
#endif

	return (pos < 0.0f) ? g_ig.dxyz.z - h : h;
}

//__device__ __inline__ static float
//fetchHeightLut(const int tidx, const float x, const float y, const float z, const bool useReg)
//{
//	const float4 Pl = tex1Dfetch(planeTex, tidx);
//	//const float4 Pl = c_plane[tidx];// tex1Dfetch(planeTex, tidx);
//	const float dist = -(Pl.x*x + Pl.y*y + Pl.z*z + Pl.w);
//
//	float h;
//	if (useReg)
//	{
//		//const float dz2 = 0.5f*g_ig.dxyz.z;
//		h = (fabsf(dist) < c_dz2) ? c_dz2 - fabsf(dist) : 0.0f;
//	}
//	else
//	{
//		const float polar = tex1Dfetch(polarTex, tidx);
//		//const float polar = c_polar[tidx];// tex1Dfetch(polarTex, tidx);
//
//		h = tex3D(heightTex,
//			fabsf(dist)*c_hlutOffset.x + 0.5f,
//			fabsf(polar)*c_hlutOffset.y + 0.5f,
//			c_viewang*c_hlutOffset.z + 0.5f
//		);
//
//	}
//
//	return (dist > 0.0f) ? g_ig.dxyz.z - h : h;
//}

__device__ __inline__ static float
cmpHeight(const float x, const float y, const float z, const float t)
{
	const float L = sqrtf(t*t + SQR(g_cg.dsd));
	const float cs = g_cg.dsd / L;
	const float sn = t / L;

	const float3 normal = make_float3(c_sn_cs_view.x*sn, -c_sn_cs_view.y*sn, cs);
	const float w = -g_cg.dso*sn;
	const float dist = -(normal.x*x + normal.y*y + normal.z*z + w);

	const float h = (fabsf(dist) < c_dz2) ? c_dz2 - fabsf(dist) : 0.0f;

	return (dist > 0.0f) ? g_ig.dxyz.z - h : h;
}





