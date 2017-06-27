#pragma once

#include "mex.h"
#ifdef printf
#undef printf
#endif

#include "helper_math.h"
#include <string>


class mexArg
{
public:
	static std::string getString(const mxArray* mx, const char* arg = nullptr)
	{
		char *buf;
		std::string str;

		int n = mxGetM(mx)*mxGetN(mx) + 1;

		if (!mxIsChar(mx) && arg)
		{
			mexPrintf("%s must be char array\n", arg);
			return str;
		}

		buf = (char*)mxCalloc(n, sizeof(char));
		mxGetString(mx, buf, n);

		str = std::string(buf);
		return str;
	}

	template<typename T>
	static T getNumber(const mxArray* mx, const int pos = 0)
	{
		T* ptr = (T*)mxGetPr(mx);
		return ptr[pos];
	}

	template<typename T>
	static float3 getFloat3(const mxArray* mx)
	{
		T* ptr = (T*)mxGetPr(mx);
		return make_float3(float(ptr[0]), float(ptr[1]), float(ptr[2]));
	}

	template<typename T>
	static float2 getFloat2(const mxArray* mx)
	{
		T* ptr = (T*)mxGetPr(mx);
		return make_float2(float(ptr[0]), float(ptr[1]));
	}

	template<typename T>
	static float getFloat(const mxArray* mx)
	{
		T* ptr = (T*)mxGetPr(mx);
		return (float)ptr[0];
	}

	template<typename T>
	static int getInt(const mxArray* mx)
	{
		T* ptr = (T*)mxGetPr(mx);
		return int(ptr[0]);
	}

	template<typename T>
	static int2 getInt2(const mxArray* mx)
	{
		T* ptr = (T*)mxGetPr(mx);
		return make_int2(int(ptr[0]), int(ptr[1]));
	}

	template<typename T>
	static int3 getInt3(const mxArray* mx)
	{
		T* ptr = (T*)mxGetPr(mx);
		return make_int3(int(ptr[0]), int(ptr[1]), int(ptr[2]));
	}

};

