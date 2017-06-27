#pragma once

#include "mat.h"

#include <unordered_map>
#include "def_cbct.h"



class mexMat
{
public:
	mexMat();
	~mexMat();

	static int fetch_cg(const std::string fn, cbct_cg& cg);
	static int fetch_ig(const std::string fn, cbct_ig& ig);
	static int fetch_proj(const std::string fn, float*& proj, float*& betas, int& na);
	static int fetch_orig_proj(const std::string fn, float*& proj, float*& betas, int& na);
	static int fetch_proj_ad(const std::string fn, float*& proj, float*& ad, int& na);
	static int fetch_image(const std::string fn, float*& image);
	static int fetch_proj(const std::string fn, float*& proj);
	static int fetch_betas(const std::string fn, float*& betas, int& na);
	static int fetch_ad(const std::string fn, float*& betas, int& na);

	static int fetch_alut(const std::string fn, lut_area& lut);
	static int fetch_alut_pos(const std::string fn, lut_area& lut);
	static int fetch_hlut(const std::string fn, lut_height& lut);
	static int fetch_hlut_pos(const std::string fn, lut_height& lut);

	static int fetch_var(const std::string fn, const std::string name, mxArray*& pa);


	static void write(std::string fn, std::string dataname, float* data, int nx, int ny, int nz);
	static void write(std::string fn, std::string dataname, int* data, int nx, int ny, int nz);
	static void write(std::string fn,
		std::vector<std::string> dataname,
		std::vector<float*> data,
		std::vector<int3> dim);
};

