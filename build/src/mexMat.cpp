#include "mexMat.h"

#include <vector>
#include <string>

#include "mexArg.h"
#include <iostream>


mexMat::mexMat()
{
}

mexMat::~mexMat()
{
}

int mexMat::fetch_cg(std::string fn, cbct_cg& cg)
{
	const std::vector<std::string> varnames{"dso", "dsd", "nst", "dst", "ost"};

	MATFile* pmat;
	mxArray* pa;

	// Open file to get directory
	pmat = matOpen(fn.c_str(), "r");
	if (pmat == nullptr)
	{
		printf("Error opening file %s\n", fn.c_str());
		return -1;
	}

	for (auto name : varnames)
	{
		pa = matGetVariable(pmat, name.c_str());
		if (pa == nullptr) {
			matClose(pmat);
			return -2;
		}

		if      (name == "dso")    cg.dso = mexArg::getFloat<double>(pa);
		else if (name == "dsd")    cg.dsd = mexArg::getFloat<double>(pa);
		else if (name == "nst")    cg.nxy = mexArg::getInt2<double>(pa);
		else if (name == "dst")    cg.dxy = mexArg::getFloat2<double>(pa);
		else if (name == "ost")    cg.offset = mexArg::getFloat2<double>(pa);

		mxDestroyArray(pa);
	}

	matClose(pmat);

	cg.wxy = make_float2
	(
		0.5f*(cg.nxy.x - 1) + cg.offset.x,
		0.5f*(cg.nxy.y - 1) + cg.offset.y
	);

	cg.invdxy = make_float2
	(
		1.0f / cg.dxy.x,
		1.0f / cg.dxy.y
	);

	return 0;
}

int mexMat::fetch_ig(std::string fn, cbct_ig& ig)
{
	const std::vector<std::string> varnames{ "nxyz", "dxyz", "oxyz" };

	MATFile* pmat;
	mxArray* pa;

	// Open file to get directory
	pmat = matOpen(fn.c_str(), "r");
	if (pmat == nullptr)
	{
		printf("Error opening file %s\n", fn.c_str());
		return -1;
	}

	for (auto name : varnames)
	{
		pa = matGetVariable(pmat, name.c_str());
		if (pa == nullptr) {
			matClose(pmat);
			return -2;
		}
		if      (name == "nxyz")   ig.nxyz = mexArg::getInt3<double>(pa);
		else if (name == "dxyz")   ig.dxyz = mexArg::getFloat3<double>(pa);
		else if (name == "oxyz")   ig.offset = mexArg::getFloat3<double>(pa);

		mxDestroyArray(pa);
	}

	matClose(pmat);

	ig.wxyz = make_float3
	(
		0.5f*(ig.nxyz.x - 1) + ig.offset.x,
		0.5f*(ig.nxyz.y - 1) + ig.offset.y,
		0.5f*(ig.nxyz.z - 1) + ig.offset.z
	);

	return 0;
}

int mexMat::fetch_proj(const std::string fn, float*& proj, float*& betas, int& na)
{
	const std::vector<std::string> varnames{ "proj", "betas", "na" };

	MATFile* pmat;
	mxArray* pa;

	// Open file to get directory
	pmat = matOpen(fn.c_str(), "r");
	if (pmat == nullptr)
	{
		printf("Error opening file %s\n", fn.c_str());
		return -1;
	}

	for (auto name : varnames)
	{
		pa = matGetVariable(pmat, name.c_str());
		if (pa == nullptr) {
			matClose(pmat);
			return -2;
		}
		if (name == "proj") {
			proj = (float*)mxGetData(pa);
		}
		else if (name == "betas") {
			betas = (float*)mxGetData(pa);
		}
		else if (name == "na") {
			na = mexArg::getInt<double>(pa);
			mxDestroyArray(pa);
		}
	}

	matClose(pmat);
	return 0;
}

int mexMat::fetch_orig_proj(const std::string fn, float*& proj, float*& betas, int& na)
{
	const std::vector<std::string> varnames{ "orig_proj", "betas", "na" };

	MATFile* pmat;
	mxArray* pa;

	// Open file to get directory
	pmat = matOpen(fn.c_str(), "r");
	if (pmat == nullptr)
	{
		printf("Error opening file %s\n", fn.c_str());
		return -1;
	}

	for (auto name : varnames)
	{
		pa = matGetVariable(pmat, name.c_str());
		if (pa == nullptr) {
			matClose(pmat);
			return -2;
		}
		if (name == "orig_proj") {
			proj = (float*)mxGetData(pa);
		}
		else if (name == "betas") {
			betas = (float*)mxGetData(pa);
		}
		else if (name == "na") {
			na = mexArg::getInt<double>(pa);
			mxDestroyArray(pa);
		}
	}

	matClose(pmat);
	return 0;
}

int mexMat::fetch_proj_ad(const std::string fn, float*& proj, float*& ad, int& na)
{
	const std::vector<std::string> varnames{ "proj", "ad", "na", "orbit_start" };

	MATFile* pmat;
	mxArray* pa;
	float orbit_start = 0.0f;

	// Open file to get directory
	pmat = matOpen(fn.c_str(), "r");
	if (pmat == nullptr)
	{
		printf("Error opening file %s\n", fn.c_str());
		return -1;
	}

	for (auto name : varnames)
	{
		pa = matGetVariable(pmat, name.c_str());
		if (pa == nullptr) {
			matClose(pmat);
			return -2;
		}
		if (name == "proj") {
			proj = (float*)mxGetData(pa);
		}
		else if (name == "ad") {
			ad = (float*)mxGetData(pa);
		}
		else if (name == "na") {
			na = mexArg::getInt<double>(pa);
			mxDestroyArray(pa);
		}
		else if (name == "orbit_start")
		{
			orbit_start = mexArg::getFloat<double>(pa);
			mxDestroyArray(pa);
		}
	}

	for (int ia=0; ia<na; ++ia)
	{
		ad[ia] -= orbit_start;
	}

	matClose(pmat);
	return 0;
}

int mexMat::fetch_image(const std::string fn, float*& image)
{
	MATFile* pmat;
	mxArray* pa;

	// Open file to get directory
	pmat = matOpen(fn.c_str(), "r");
	if (pmat == nullptr)
	{
		printf("Error opening file %s\n", fn.c_str());
		return -1;
	}

	pa = matGetVariable(pmat, "image");
	if (pa == nullptr) {
		image = nullptr;
		matClose(pmat);
		return -2;
	}
	
	image = (float*)mxGetData(pa);

	matClose(pmat);
	return 0;
}

int mexMat::fetch_proj(const std::string fn, float*& proj)
{
	MATFile* pmat;
	mxArray* pa;

	// Open file to get directory
	pmat = matOpen(fn.c_str(), "r");
	if (pmat == nullptr)
	{
		printf("Error opening file %s\n", fn.c_str());
		return -1;
	}

	pa = matGetVariable(pmat, "proj");
	if (pa == nullptr) {
		matClose(pmat);
		return -2;
	}

	proj = (float*)mxGetData(pa);

	matClose(pmat);
	return 0;
}

int mexMat::fetch_betas(const std::string fn, float*& betas, int& na)
{
	const std::vector<std::string> varnames{ "betas", "na" };

	MATFile* pmat;
	mxArray* pa;

	// Open file to get directory
	pmat = matOpen(fn.c_str(), "r");
	if (pmat == nullptr)
	{
		printf("Error opening file %s\n", fn.c_str());
		return -1;
	}

	for (auto name : varnames)
	{
		pa = matGetVariable(pmat, name.c_str());
		if (pa == nullptr) {
			matClose(pmat);
			return -2;
		}
		if (name == "betas") {
			betas = (float*)mxGetData(pa);
		}
		else if (name == "na") {
			na = mexArg::getInt<double>(pa);
			mxDestroyArray(pa);
		}
	}

	matClose(pmat);
	return 0;
}

int mexMat::fetch_ad(const std::string fn, float*& betas, int& na)
{
	const std::vector<std::string> varnames{ "ad", "na",  "orbit_start"};

	MATFile* pmat;
	mxArray* pa;
	float orbit_start = 0.0f;

	// Open file to get directory
	pmat = matOpen(fn.c_str(), "r");
	if (pmat == nullptr)
	{
		printf("Error opening file %s\n", fn.c_str());
		return -1;
	}

	for (auto name : varnames)
	{
		pa = matGetVariable(pmat, name.c_str());
		if (pa == nullptr) {
			matClose(pmat);
			return -2;
		}
		if (name == "ad") {
			betas = (float*)mxGetData(pa);
		}
		else if (name == "na") {
			na = mexArg::getInt<double>(pa);
			mxDestroyArray(pa);
		}
		else if (name == "orbit_start")
		{
			orbit_start = mexArg::getFloat<double>(pa);
			mxDestroyArray(pa);
		}
	}

	orbit_start = 0.0f;
	for (int ia=0; ia<na; ++ia)
	{
		betas[ia] = float(DEG2RAD)* (betas[ia] - orbit_start);
	}

	matClose(pmat);
	return 0;
}

int mexMat::fetch_alut(const std::string fn, lut_area& lut)
{
	const std::vector<std::string> varnames{ "lut", "nd", "na", "ddist", "dang"};

	MATFile* pmat;
	mxArray* pa;

	// Open file to get directory
	pmat = matOpen(fn.c_str(), "r");
	if (pmat == nullptr)
	{
		printf("Error opening file %s\n", fn.c_str());
		return -1;
	}

	for (auto name : varnames)
	{
		pa = matGetVariable(pmat, name.c_str());
		if (pa == nullptr) {
			matClose(pmat);
			return -2;
		}
		if (name == "lut") {
			lut.lut = (float*)mxGetData(pa);
		}
		else if (name == "nd") {
			lut.dim.x = mexArg::getInt<double>(pa);
			mxDestroyArray(pa);
		}
		else if (name == "na") {
			lut.dim.y = mexArg::getInt<double>(pa);
			mxDestroyArray(pa);
		}
		else if (name == "ddist") {
			lut.delta.x = mexArg::getFloat<double>(pa);
			mxDestroyArray(pa);
		}
		else if (name == "dang") {
			lut.delta.y = mexArg::getFloat<double>(pa);
			mxDestroyArray(pa);
		}

	}

	matClose(pmat);
	return 0;
}

int mexMat::fetch_alut_pos(const std::string fn, lut_area& lut)
{
	const std::vector<std::string> varnames{ "alut", "np", "na", "step_pos", "step_ang" };

	MATFile* pmat;
	mxArray* pa;

	// Open file to get directory
	pmat = matOpen(fn.c_str(), "r");
	if (pmat == nullptr)
	{
		printf("Error opening file %s\n", fn.c_str());
		return -1;
	}

	for (auto name : varnames)
	{
		pa = matGetVariable(pmat, name.c_str());
		if (pa == nullptr) {
			matClose(pmat);
			return -2;
		}
		if (name == "alut") {
			lut.lut = (float*)mxGetData(pa);
		}
		else if (name == "np") {
			lut.dim.x = mexArg::getInt<double>(pa);
			mxDestroyArray(pa);
		}
		else if (name == "na") {
			lut.dim.y = mexArg::getInt<double>(pa);
			mxDestroyArray(pa);
		}
		else if (name == "step_pos") {
			lut.delta.x = mexArg::getFloat<double>(pa);
			mxDestroyArray(pa);
		}
		else if (name == "step_ang") {
			lut.delta.y = mexArg::getFloat<double>(pa);
			mxDestroyArray(pa);
		}

	}

	matClose(pmat);
	return 0;
}

int mexMat::fetch_hlut(const std::string fn, lut_height& lut)
{
	const std::vector<std::string> varnames{ "lut", "nd", "np", "na", "ddist", "dazi", "dpol" };
	
	MATFile* pmat;
	mxArray* pa;

	// Open file to get directory
	pmat = matOpen(fn.c_str(), "r");
	if (pmat == nullptr)
	{
		printf("Error opening file %s\n", fn.c_str());
		return -1;
	}

	for (auto name : varnames)
	{
		pa = matGetVariable(pmat, name.c_str());
		if (pa == nullptr) {
			matClose(pmat);
			return -2;
		}
		if (name == "lut") {
			lut.lut = (float*)mxGetData(pa);
		}
		else if (name == "nd") {
			lut.dim.x = mexArg::getInt<double>(pa);
			mxDestroyArray(pa);
			//printf("nd: %d\n", lut.dim.x);
		}
		else if (name == "np") {
			lut.dim.y = mexArg::getInt<double>(pa);
			mxDestroyArray(pa);
			//printf("np: %d\n", lut.dim.y);
		}
		else if (name == "na") {
			lut.dim.z = mexArg::getInt<double>(pa);
			mxDestroyArray(pa);
			//printf("na: %d\n", lut.dim.z);
		}
		else if (name == "ddist") {
			lut.delta.x = mexArg::getFloat<double>(pa);
			mxDestroyArray(pa);
			//printf("ddist: %f\n", lut.delta.x);
		}
		else if (name == "dazi") {
			lut.delta.z = mexArg::getFloat<double>(pa);
			mxDestroyArray(pa);
			//printf("dazi: %f\n", lut.delta.z);
		}
		else if (name == "dpol") {
			lut.delta.y = mexArg::getFloat<double>(pa);
			mxDestroyArray(pa);
			//printf("dpol: %f\n", lut.delta.y);
		}
	}

	matClose(pmat);
	return 0;
}

int mexMat::fetch_hlut_pos(const std::string fn, lut_height& lut)
{
	const std::vector<std::string> varnames{ "hlut", "npos", "npol", "nazi", "step_pos", "step_pol", "step_azi" };

	MATFile* pmat;
	mxArray* pa;

	// Open file to get directory
	pmat = matOpen(fn.c_str(), "r");
	if (pmat == nullptr)
	{
		printf("Error opening file %s\n", fn.c_str());
		return -1;
	}

	for (auto name : varnames)
	{
		pa = matGetVariable(pmat, name.c_str());
		if (pa == nullptr) {
			matClose(pmat);
			return -2;
		}
		if (name == "hlut") {
			lut.lut = (float*)mxGetData(pa);
		}
		else if (name == "npos") {
			lut.dim.x = mexArg::getInt<double>(pa);
			mxDestroyArray(pa);
			//printf("nd: %d\n", lut.dim.x);
		}
		else if (name == "npol") {
			lut.dim.y = mexArg::getInt<double>(pa);
			mxDestroyArray(pa);
			//printf("np: %d\n", lut.dim.y);
		}
		else if (name == "nazi") {
			lut.dim.z = mexArg::getInt<double>(pa);
			mxDestroyArray(pa);
			//printf("na: %d\n", lut.dim.z);
		}
		else if (name == "step_pos") {
			lut.delta.x = mexArg::getFloat<double>(pa);
			mxDestroyArray(pa);
			//printf("ddist: %f\n", lut.delta.x);
		}
		else if (name == "step_azi") {
			lut.delta.z = mexArg::getFloat<double>(pa);
			mxDestroyArray(pa);
			//printf("dazi: %f\n", lut.delta.z);
		}
		else if (name == "step_pol") {
			lut.delta.y = mexArg::getFloat<double>(pa);
			mxDestroyArray(pa);
			//printf("dpol: %f\n", lut.delta.y);
		}
	}

	matClose(pmat);
	return 0;
}

int mexMat::fetch_var(const std::string fn, const std::string name, mxArray*& pa)
{
	MATFile* pmat;

	// Open file to get directory
	pmat = matOpen(fn.c_str(), "r");
	if (pmat == nullptr)
	{
		printf("Error opening file %s\n", fn.c_str());
		return -1;
	}

	pa = matGetVariable(pmat, name.c_str());

	matClose(pmat);

	if (pa == nullptr)
		return -2;

	return 0;
}


void mexMat::write(std::string fn, std::string dataname, float* data, int nx, int ny, int nz)
{
	MATFile* pmat;
	mxArray* pa;
	matError status;

	pmat = matOpen(fn.c_str(), "w7.3");
	if (pmat == nullptr)
	{
		printf("Error creating file %s\n", fn.c_str());
		printf("(Do you have write permission in this directory?)\n");
		return;
	}
	//printf("matOpen\n");

	mwSize dims[3] = {nx, ny, nz};
	pa = mxCreateNumericArray(3, dims, mxSINGLE_CLASS, mxREAL);
	if (pa == nullptr)
	{
		printf("%s : Out of memory on line %d\n", __FILE__, __LINE__);
		printf("Unable to create mxArray.\n");
		return;
	}
	//printf("matAlloc\n");

	memcpy(mxGetData(pa), (void*)data, nx * ny * nz * sizeof(float));
	//printf("matCopy\n");

	status = matPutVariable(pmat, dataname.c_str(), pa);
	if (status != 0)
	{
		printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
		return;
	}
	//printf("matWrite\n");

	mxDestroyArray(pa);
	if (matClose(pmat) != 0)
	{
		printf("Error closing file %s\n", fn.c_str());
		return;
	}
	//printf("matClose\n");
}

void mexMat::write(std::string fn, std::string dataname, int* data, int nx, int ny, int nz)
{
	MATFile* pmat;
	mxArray* pa;
	matError status;

	pmat = matOpen(fn.c_str(), "w7.3");
	if (pmat == nullptr)
	{
		printf("Error creating file %s\n", fn.c_str());
		printf("(Do you have write permission in this directory?)\n");
		return;
	}
	//printf("matOpen\n");

	mwSize dims[3] = { nx, ny, nz };
	pa = mxCreateNumericArray(3, dims, mxINT32_CLASS, mxREAL);
	if (pa == nullptr)
	{
		printf("%s : Out of memory on line %d\n", __FILE__, __LINE__);
		printf("Unable to create mxArray.\n");
		return;
	}
	//printf("matAlloc\n");

	memcpy(mxGetData(pa), (void*)data, nx * ny * nz * sizeof(float));
	//printf("matCopy\n");

	status = matPutVariable(pmat, dataname.c_str(), pa);
	if (status != 0)
	{
		printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
		return;
	}
	//printf("matWrite\n");

	mxDestroyArray(pa);
	if (matClose(pmat) != 0)
	{
		printf("Error closing file %s\n", fn.c_str());
		return;
	}
	//printf("matClose\n");
}

void mexMat::write(std::string fn, std::vector<std::string> dataname, std::vector<float*> data, std::vector<int3> dim)
{
	MATFile* pmat;
	mxArray* pa;
	matError status;

	pmat = matOpen(fn.c_str(), "w7.3");
	if (pmat == nullptr)
	{
		printf("Error creating file %s\n", fn.c_str());
		printf("(Do you have write permission in this directory?)\n");
		return;
	}

	const int ndata = dataname.size();
	for (int i=0; i<ndata; ++i)
	{
		const int3 nxyz = dim[i];
		mwSize dims[3] = { nxyz.x, nxyz.y, nxyz.z };
		pa = mxCreateNumericArray(3, dims, mxSINGLE_CLASS, mxREAL);
		if (pa == nullptr)
		{
			printf("%s : Out of memory on line %d\n", __FILE__, __LINE__);
			printf("Unable to create mxArray.\n");
			//return;
		}

		memcpy(mxGetData(pa), (void*)data[i], nxyz.x * nxyz.y * nxyz.z * sizeof(float));

		//std::cout << dataname[i].c_str() << std::endl;

		status = matPutVariable(pmat, dataname[i].c_str(), pa);
		if (status != 0)
		{
			printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
			//return;
		}

		mxDestroyArray(pa);
	}

	if (matClose(pmat) != 0)
	{
		printf("Error closing file %s\n", fn.c_str());
		return;
	}
}
