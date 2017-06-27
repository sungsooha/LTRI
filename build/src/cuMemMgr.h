#pragma once

#include "def_cbct.h"
#include <string>
#include <unordered_map>

class cuMemMgr
{
public:
	cuMemMgr();
	~cuMemMgr();

	void freeAll();
	void freeArray(std::string name);
	void freeLinear(std::string name);

	void show();

	template<class T> void alloc3DArray(std::string name, int3 dim)
	{
		freeArray(name);

		cudaArray* d_array;
		cudaChannelFormatDesc desc = cudaCreateChannelDesc<T>();
		cudaMalloc3DArray(&d_array, &desc, make_cudaExtent(dim.x, dim.y, dim.z));
		CUDA_CHECK_ERROR("ctMem::alloc3DArray");

		//printf("%s: %d x %d x %d\n", name.c_str(), dim.x, dim.y, dim.z);

		cuArrayMap[name] = d_array;
		dimMap[name] = make_int4(dim, dim.x);
	}
	template<class T> void alloc3DArraySurface(std::string name, int3 dim)
	{
		freeArray(name);

		cudaArray* d_array;
		cudaChannelFormatDesc desc = cudaCreateChannelDesc<T>();
		cudaMalloc3DArray(&d_array, &desc, make_cudaExtent(dim.x, dim.y, dim.z), cudaArraySurfaceLoadStore);
		CUDA_CHECK_ERROR("ctMem::alloc3DArraySurface");

		cuArrayMap[name] = d_array;
		dimMap[name] = make_int4(dim, dim.x);
	}
	template<class T> void alloc3DArrayLayer(std::string name, int3 dim)
	{
		freeArray(name);

		cudaArray* d_array;
		cudaChannelFormatDesc desc = cudaCreateChannelDesc<T>();
		cudaMalloc3DArray(&d_array, &desc, make_cudaExtent(dim.x, dim.y, dim.z), cudaArrayLayered);
		CUDA_CHECK_ERROR("ctMem::alloc3DArrayLayer");

		cuArrayMap[name] = d_array;
		dimMap[name] = make_int4(dim, dim.x);
	}
	template<class T> void allocArray(std::string name, int2 dim)
	{
		freeArray(name);

		cudaArray* d_array;
		cudaChannelFormatDesc desc = cudaCreateChannelDesc<T>();
		cudaMallocArray(&d_array, &desc, dim.x, dim.y);
		CUDA_CHECK_ERROR("ctMem::allocArray");

		cuArrayMap[name] = d_array;
		dimMap[name] = make_int4(dim.x, dim.y, 1, dim.x);
	}

	template<class T> void copyTo3DArray(std::string arrname, T* h_data)
	{
		auto it = cuArrayMap.find(arrname);
		if (it == cuArrayMap.end())
		{
			printf("Cannot find allocated memory with name of %s.\n", arrname.c_str());
			exit(-1);
		}

		cudaArray* d_array = cuArrayMap[arrname];
		int4 dim = dimMap[arrname];

		//printf("%s: %d x %d x %d\n", arrname.c_str(), dim.x, dim.y, dim.z);

		cudaMemcpy3DParms cp = { 0 };
		cp.srcPtr = make_cudaPitchedPtr((void*)h_data, dim.x * sizeof(T), dim.x, dim.y);
		cp.dstArray = d_array;
		cp.extent = make_cudaExtent(dim.x, dim.y, dim.z);
		cp.kind = cudaMemcpyHostToDevice;
		cudaMemcpy3D(&cp);

		CUDA_CHECK_ERROR("ctMem::copyTo3DArray");
	}
	template<class T> void copyTo3DArray(std::string src, std::string dstarr)
	{
		//printf("copyTo3DArray(src, dstarr)\n");

		int4 dstdim = getDim(dstarr);
		int4 srcdim = getDim(src);

		if (dstdim.x != srcdim.x || dstdim.y*dstdim.z != srcdim.y*srcdim.z )
		{
			printf("Mismatched dimension between %s (src) and %s (dst)\n", src.c_str(), dstarr.c_str());
			exit(-1);
		}

		cudaMemcpy3DParms cp = { 0 };
		cp.srcPtr = make_cudaPitchedPtr(cuLinearMap[src], srcdim.w * sizeof(T), dstdim.x, dstdim.y);
		cp.dstArray = cuArrayMap[dstarr];
		cp.extent = make_cudaExtent(dstdim.x, dstdim.y, dstdim.z);
		cp.kind = cudaMemcpyDeviceToDevice;
		cudaMemcpy3D(&cp);
	}

	template<class T> void copyToArray(std::string name, T* h_data)
	{
		auto it = cuArrayMap.find(name);
		if (it == cuArrayMap.end())
		{
			printf("Cannot find allocated memory with name of %s.\n", name.c_str());
			exit(-1);
		}

		cudaArray* d_array = cuArrayMap[name];
		int4 dim = dimMap[name];
		cudaMemcpyToArray(
			d_array, 0, 0,
			h_data, dim.x*dim.y * sizeof(T),
			cudaMemcpyHostToDevice);
		CUDA_CHECK_ERROR("ctMem::copyToArray");
	}
	template<class T> void copyToArray(std::string src, std::string dstarr, int2 dim)
	{
		cudaMemcpyToArray(
			cuArrayMap[dstarr], 0, 0,
			cuLinearMap[src], dim.x*dim.y * sizeof(T),
			cudaMemcpyDeviceToDevice);
	}
	template<class T> void copyFrom3DArray(std::string src, T* h_dst)
	{
		int4 dim = dimMap[src];

		cudaMemcpy3DParms cp = { 0 };
		cp.srcArray = cuArrayMap[src];
		cp.dstPtr = make_cudaPitchedPtr((void*)h_dst, dim.x * sizeof(T), dim.x, dim.y);
		cp.extent = make_cudaExtent(dim.x, dim.y, dim.z);
		cp.kind = cudaMemcpyDeviceToHost;
		cudaMemcpy3D(&cp);
	}
	cudaArray* getArray(std::string name, int4* dim = nullptr)
	{
		auto it = cuArrayMap.find(name);
		if (it == cuArrayMap.end())
		{
			printf("Cannot find allocated memory with name of %s.\n", name.c_str());
			exit(-1);
		}
		if (dim) *dim = dimMap[name];
		return it->second;
	}

	template<class T> int allocLinearPitch(std::string name, int2 dim)
	{
		freeLinear(name);

		//printf("%s: %d x %d\n", name.c_str(), dim.x, dim.y);

		T * d_mem;
		size_t pitch;
		cudaMallocPitch((void**)&d_mem, &pitch, dim.x * sizeof(T), (size_t)dim.y);
		cudaMemset2D(d_mem, pitch, 0, dim.x * sizeof(T), dim.y);

		CUDA_CHECK_ERROR("ctMem::allocLinearPitch");

		cuLinearMap[name] = (void*)d_mem;
		dimMap[name] = make_int4(dim.x, dim.y, 1, (int)pitch / sizeof(T));
		return (int)pitch / sizeof(T);
	}
	template<class T> void allocLinear(std::string name, int N)
	{
		freeLinear(name);

		T * d_mem;
		cudaMalloc((void**)&d_mem, N * sizeof(T));
		cudaMemset(d_mem, 0, N * sizeof(T));

		CUDA_CHECK_ERROR("ctMem::allocLinear");

		cuLinearMap[name] = (void*)d_mem;
		dimMap[name] = make_int4(N, 1, 1, N);
	}

	template<class T> void memsetLinearPitch(std::string name, int value = 0)
	{
		int4 dim = dimMap[name];
		cudaMemset2D(
			cuLinearMap[name], dim.w * sizeof(T), value,
			dim.x * sizeof(T), dim.y);
	}
	template<class T> void memsetLinear(std::string name)
	{
		int4 dim = dimMap[name];
		cudaMemset(cuLinearMap[name], 0, dim.x * sizeof(T));
	}
	template<class T> T* getLinear(std::string name, int4* dim = nullptr)
	{
		auto it = cuLinearMap.find(name);
		if (it == cuLinearMap.end())
		{
			printf("Cannot find allocated memory with name of %s.\n", name.c_str());
			exit(-1);
		}
		if (dim) *dim = dimMap[name];
		return (T*)it->second;
	}

	template<class T> void copyToLinearPitch(std::string name, T* src)
	{
		int4 dim = dimMap[name];
		cudaMemcpy2D(
			cuLinearMap[name], dim.w * sizeof(T),
			src, dim.x * sizeof(T),
			dim.x * sizeof(T), dim.y,
			cudaMemcpyHostToDevice);
		CUDA_CHECK_ERROR("ctMem::copyToLinearPitch");
	}
	template<class T> void copyToLinearPitch(std::string name, T* src, int2 srcdim)
	{
		int4 dim = dimMap[name];
		cudaMemcpy2D(
			cuLinearMap[name], dim.w * sizeof(T),
			src, srcdim.x * sizeof(T),
			srcdim.x * sizeof(T), srcdim.y,
			cudaMemcpyHostToDevice);
		CUDA_CHECK_ERROR("ctMem::copyToLinearPitch");
	}
	template<class T> void copyToLinear(std::string name, T* src)
	{
		int4 dim = dimMap[name];
		cudaMemcpy(
			cuLinearMap[name],
			src,
			dim.x * sizeof(T),
			cudaMemcpyHostToDevice);
		CUDA_CHECK_ERROR("ctMem::copyToLinear");
	}
	template<class T> void copyToLinear(std::string name, int offset, T* src, int N)
	{
		cudaMemcpy(
			(void*)(getLinear<T>(name) + offset),
			(void*)src,
			N * sizeof(T),
			cudaMemcpyHostToDevice);
		CUDA_CHECK_ERROR("ctMem::copyToLinear");
	}
	template<class T> void copyLinearPitchToHost(std::string name, T* dst)
	{
		int4 dim = dimMap[name];
		cudaMemcpy2D(
			dst, dim.x * sizeof(T),
			getLinear<T>(name), dim.w * sizeof(T),
			dim.x * sizeof(T), dim.y,
			cudaMemcpyDeviceToHost);
		//printf("%s\n", name.c_str());
		CUDA_CHECK_ERROR("ctMem::copyLinearPitchToHost");
	}
	template<class T> void copyLinearPitchToHost(std::string name, T* dst, int2 dstDim)
	{
		int4 dim = dimMap[name];
		cudaMemcpy2D(
			dst, dstDim.x * sizeof(T),
			getLinear<T>(name), dim.w * sizeof(T),
			dstDim.x * sizeof(T), dstDim.y,
			cudaMemcpyDeviceToHost);
		CUDA_CHECK_ERROR("ctMem::copyToLinearPitch");
	}
	template<class T> void copyLinearToHost(std::string name, T* dst)
	{
		int4 dim = dimMap[name];
		cudaMemcpy(
			dst,
			getLinear<T>(name),
			dim.x * sizeof(T),
			cudaMemcpyDeviceToHost);
		CUDA_CHECK_ERROR("ctMem::copyLinearToHost");
	}

	int4 getDim(std::string name)
	{
		return dimMap[name];
	}

private:

	std::unordered_map<std::string, cudaArray*> cuArrayMap;
	std::unordered_map<std::string, void*>      cuLinearMap;
	std::unordered_map<std::string, int4>       dimMap;
};

