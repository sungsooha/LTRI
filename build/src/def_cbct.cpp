#include "def_cbct.h"

#include <iostream>
#include <string>

void CUDA_CHECK_ERROR(const char* str)
{
	cudaDeviceSynchronize();
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err)
	{
#ifdef USEMATLAB
		char buffer[2048];
		sprintf(buffer, "[%s ERROR]: %s\n", str, cudaGetErrorString(err));
		mexErrMsgTxt(buffer);
#else
		std::cerr << "[" << str << " ERROR]: " << cudaGetErrorString(err) << std::endl;
		exit(EXIT_FAILURE);
#endif // USEMATLAB

	}
}

std::string CUDA_SET_BEST_DEVICE()
{
	int num_devices;
	cudaGetDeviceCount(&num_devices);
	std::string deviceName;
	int max_multiprocessors = 0, max_device = 0;
	for (int device = 0; device<num_devices; device++)
	{
		cudaDeviceProp properties;
		cudaGetDeviceProperties(&properties, device);
		if (max_multiprocessors < properties.multiProcessorCount)
		{
			max_multiprocessors = properties.multiProcessorCount;
			max_device = device;
			deviceName = std::string(properties.name);
		}
	}
	cudaSetDevice(max_device);
//#ifdef USEMATLAB
//	mexPrintf("%s\n", deviceName.c_str());
//#else
//	std::cout << deviceName << std::endl;
//#endif

	return deviceName;
}

void show_cbct_ig(const cbct_ig& ig)
{
	printf("Image geometry:\n");
	printf("[nx, ny, nz] : [%d, %d, %d]\n", ig.nxyz.x, ig.nxyz.y, ig.nxyz.z);
	printf("[dx, dy, dz] : [%.3f, %.3f, %.3f] mm\n", ig.dxyz.x, ig.dxyz.y, ig.dxyz.z);
	printf("offset_[x, y, z]: [%.3f, %.3f, %.3f] mm\n", ig.offset.x, ig.offset.y, ig.offset.z);
}

void show_cbct_cg(const cbct_cg& cg)
{
	printf("Cone-beam geometry:\n");
	printf("[dso, dsd] : [%.3f, %.3f] mm\n", cg.dso, cg.dsd);
	printf("%s-detector\n", "Flat");
	printf("[ns, nt] : [%d, %d]\n", cg.nxy.x, cg.nxy.y);
	printf("[ds, dt] : [%.3f, %.3f]\n", cg.dxy.x, cg.dxy.y);
	printf("offset_[s, t]: [%.3f, %.3f] mm\n", cg.offset.x, cg.offset.y);
}

void cbct::Medtronic(cbct_cg& cg)
{
	cg.dsd = 1147.7f;
	cg.dso = 647.7f;
	cg.nxy = make_int2(1024, 384);
	cg.dxy = make_float2(393.432f/1024.0f, -290.224f/384.0f);
	cg.offset = make_float2(0.0f, 0.0f);
	cg.wxy = make_float2(
		0.5f*(cg.nxy.x - 1),
		0.5f*(cg.nxy.y - 1)
	);
	cg.invdxy = 1.0f / cg.dxy;
}

void cbct::TimeGeometry(cbct_cg& cg, cbct_ig& ig, float scale)
{
	cg.dsd = 949.075f;
	cg.dso = 647.7f;
	cg.nxy = make_int2(512, 512);
	cg.dxy = make_float2(1.0279f, -1.0964f);
	cg.offset = make_float2(0.0f, 0.0f);
	cg.wxy = make_float2(
		0.5f*(cg.nxy.x - 1) + cg.offset.x,
		0.5f*(cg.nxy.y - 1) + cg.offset.y
	);
	cg.invdxy = 1.0f / cg.dxy;

	ig.nxyz = make_int3(512, 512, 512);
	ig.dxyz = make_float3(0.4883f, -0.4883f, 0.6250f);
	ig.offset = make_float3(0.0f, 0.0f, 0.0f);

	ig.nxyz.x = (int)floorf(ig.nxyz.x * scale);
	ig.nxyz.y = (int)floorf(ig.nxyz.y * scale);
	ig.nxyz.z = (int)floorf(ig.nxyz.z * scale);

	ig.dxyz.x = ig.dxyz.x / scale;
	ig.dxyz.y = ig.dxyz.y / scale;
	ig.dxyz.z = ig.dxyz.z / scale;

	ig.wxyz = make_float3(
		0.5f*(ig.nxyz.x - 1) + ig.offset.x,
		0.5f*(ig.nxyz.y - 1) + ig.offset.y,
		0.5f*(ig.nxyz.z - 1) + ig.offset.z
	);
}


int iDivUp(int a, int b) {
	return (a % b != 0) ? (a / b + 1) : (a / b);
}

unsigned int nextPowerOf2(unsigned int n)
{
	unsigned int p = 1;
	if (n && !(n & (n - 1)))
		return n;

	while (p < n)
		p <<= 1;

	return p;
}

// Giga Updates Per Second 
float GUPS(int3 voldim, int na, float sec)
{
	return
		float(voldim.x*voldim.y*voldim.z)*float(na) / float(1024 * 1024 * 1024) / sec;
}

float2 rotateCCW_z(float2 v, const float cs, const float sn)
{
	return make_float2(v.x*cs - v.y*sn, v.x*sn + v.y*cs);
}
