__global__ static void
apply_geodiv_kernel(float* proj, const int pitch)
{
	const int is = blockIdx.x*blockDim.x + threadIdx.x;
	const int it = blockIdx.y*blockDim.y + threadIdx.y;

	if (is >= g_cg.nxy.x || it >= g_cg.nxy.y)
		return;

	proj[is*pitch + it] *= tex2D(geodivTex, is + 0.5f, it + 0.5f);
}

__global__ static void
apply_geodiv_kernel_st(float* proj, const int pitch)
{
	const int is = blockIdx.x*blockDim.x + threadIdx.x;
	const int it = blockIdx.y*blockDim.y + threadIdx.y;

	if (is >= g_cg.nxy.x || it >= g_cg.nxy.y)
		return;

	proj[it*pitch + is] *= tex2D(geodivTex, is + 0.5f, it + 0.5f);
}