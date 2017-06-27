#include "cuMemMgr.h"



cuMemMgr::cuMemMgr()
{
}


cuMemMgr::~cuMemMgr()
{
	freeAll();
}

void cuMemMgr::freeAll()
{
	if (!cuArrayMap.empty())
	{
		for (auto it = cuArrayMap.begin(); it != cuArrayMap.end(); ++it)
			cudaFreeArray(it->second);
		cuArrayMap.clear();
	}

	if (!cuLinearMap.empty())
	{
		for (auto it = cuLinearMap.begin(); it != cuLinearMap.end(); ++it)
			cudaFree(it->second);
		cuLinearMap.clear();
	}

	if (!dimMap.empty())
		dimMap.clear();
}

void cuMemMgr::freeArray(std::string name)
{
	auto it = cuArrayMap.find(name);
	if (it != cuArrayMap.end())
	{
		cudaFreeArray(it->second);
		cuArrayMap.erase(name);
		dimMap.erase(name);
	}
}

void cuMemMgr::freeLinear(std::string name)
{
	auto it = cuLinearMap.find(name);
	if (it != cuLinearMap.end())
	{
		cudaFree(it->second);
		cuLinearMap.erase(name);
		dimMap.erase(name);
	}
}

void cuMemMgr::show()
{
	std::unordered_map<std::string, int4>::iterator iit;

	printf("\nGPU memory usages:\n");
	if (!cuArrayMap.empty())
	{
		printf("- CUDA array:\n");
		for (auto it = cuArrayMap.begin(); it != cuArrayMap.end(); ++it)
		{
			std::string name = it->first;
			int4 dim = make_int4(-1, -1, -1, -1);
			if ((iit = dimMap.find(name)) != dimMap.end())
				dim = iit->second;
			printf("\t%s : %d x %d x %d (%d)\n", name.c_str(), dim.x, dim.y, dim.z, dim.w);
		}
		printf("\n");
	}

	if (!cuLinearMap.empty())
	{
		printf("- CUDA linear:\n");
		for (auto it = cuLinearMap.begin(); it != cuLinearMap.end(); ++it)
		{
			std::string name = it->first;
			int4 dim = make_int4(-1, -1, -1, -1);
			if ((iit = dimMap.find(name)) != dimMap.end())
				dim = iit->second;
			printf("\t%s : %d x %d x %d (%d)\n", name.c_str(), dim.x, dim.y, dim.z, dim.w);
		}
		printf("\n");
	}
}