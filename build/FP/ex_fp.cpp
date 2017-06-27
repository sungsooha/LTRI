#include "cbctProjector.h"
#include "mexMat.h"
#include "mexArg.h"
#include <iostream>
#include <fstream>

#define USE_CMDARG

static void Usage()
{
	printf("\n\
		Forward projection\n\
		-r || --root : root (or working) directory path\n\
		-i || --input : input mat filename (including extension, relative to root directory)\n\
		-m || --method : integer number (0: LUT_LUT, 1: LUT_REG, 2: LUT_DD, 3: SF_TT, 4: SF_TR)\n\
		-o || --output [optional]: output mat filename (including extension, relative to root directory, 'no' will do nothing)\n\
		-t || --time [optional]: measure time performance (0: no, 1: yes) [default: 0]\n\
		-s || --scale [optional]: 1: scaling on-the-fly, 0: scaling later [default: 1]\n\
		-a || --alut: area look-up table mat file name with full path relative to root\n\
		-v || --hlut: height look-up table mat file name with full path relative to root \n\
		-ost || --offset_st: use user provided detector offset for each views \n\
		-h || --help : show this message\n\
		\n\
		Note: for the LL, LR and LD, it expects look-up table(s) is in lut sub-directory of the root.\n");
	exit(-1);
}

void load_test_setting(
	cbct_cg& cg,
	cbct_ig& ig,
	float*&  image,
	int&     na,
	float*&  betas,
	int&     nfoot,
	lut_area& alut,
	lut_height& hlut
)
{
	cbct::Medtronic(cg);
	ig.nxyz = make_int3(512, 512, 512);
	ig.dxyz = make_float3(0.415f, -0.415f, 0.312f);
	ig.offset = make_float3(0.0f, 0.0f, 0.0f);
	ig.wxyz = make_float3(
		0.5f*(ig.nxyz.x - 1.0f),
		0.5f*(ig.nxyz.y - 1.0f),
		0.5f*(ig.nxyz.z - 1.0f)
	);

	image = new float[ig.nxyz.x*ig.nxyz.y*ig.nxyz.z];
	memset(image, 0, ig.nxyz.x*ig.nxyz.y*ig.nxyz.z * sizeof(float));
	for (int iz = ig.nxyz.z / 2 - 100; iz < ig.nxyz.z / 2 + 100; ++iz)
		for (int iy = ig.nxyz.y / 2 - 50; iy < ig.nxyz.y / 2 + 50; ++iy)
			for (int ix = ig.nxyz.x / 2 - 50; ix < ig.nxyz.x / 2 + 50; ++ix)
				image[(ix + iy*ig.nxyz.x)*ig.nxyz.z + iz] = 1.0f;

	na = 1;
	betas = new float[na];
	betas[0] = 230.0f * float(DEG2RAD);

	nfoot = 5;

	mexMat::fetch_alut("./test/test_alut.mat", alut);
	mexMat::fetch_hlut("./test/test_hlut.mat", hlut);
}

void load_cuboid_setting(
	std::string matfn,
	cbct_cg& cg,
	cbct_ig& ig,
	float*&  image,
	int&     na,
	float*&  betas,
	int&     nfoot,
	std::string alutfn, lut_area& alut,
	std::string hlutfn, lut_height& hlut
)
{
	mxArray* pa;

	mexMat::fetch_cg(matfn, cg);
	mexMat::fetch_ig(matfn, ig);
	mexMat::fetch_image(matfn, image);
	mexMat::fetch_ad(matfn, betas, na);


	mexMat::fetch_var(matfn, "nfoot", pa);
	nfoot = mexArg::getInt<double>(pa);
	mxDestroyArray(pa);

	//if (!alutfn.empty()) mexMat::fetch_alut(alutfn, alut);
	if (!alutfn.empty()) mexMat::fetch_alut_pos(alutfn, alut);
	//if (!hlutfn.empty()) mexMat::fetch_hlut(hlutfn, hlut);
	if (!hlutfn.empty()) mexMat::fetch_hlut_pos(hlutfn, hlut);

	if (nfoot % 2 == 0) nfoot++;
}

void load_offset_st(
	std::string matfn,
	float*& offset_s, 
	float*& offset_t
)
{
	mxArray* pa;

	if (mexMat::fetch_var(matfn, "all_offset_s", pa) < 0)
	{
		;
	}
	offset_s = (float*)mxGetData(pa);
	//mxDestroyArray(pa);

	mexMat::fetch_var(matfn, "all_offset_t", pa);
	offset_t = (float*)mxGetData(pa);
	//mxDestroyArray(pa);
}

int main(int argc, char** argv)
{
	cuMemMgr              mm;
	
	cbctProjector         sys;
	cbctProjector::method method = cbctProjector::SF_TR;
	std::string           methodName = "TR";
	cbctProjector::order  order = cbctProjector::TS;

	cbct_cg               cg;
	cbct_ig               ig;

	float*                proj = nullptr;
	int                   pitch;

	float*                image = nullptr;

	float*                betas = nullptr;
	int                   na;

	int                   nfoot;

	lut_area              alut;
	lut_height            hlut;

	std::string           root;
	std::string           matFn;
	std::string           outFn;
	std::string           gpuname;

	std::string           alutFn = "alut.mat";
	std::string           hlutFn = "hlut.mat";

	bool                  useMat = true;
	bool                  measureTime = false;
	bool                  makeOutput = true;

	cudaEvent_t           start, stop;
	float                 elapsedTime = 0.0f;

	bool                  doScaling = true;

	bool                  useOffset_st = false;
	float*                offset_s = nullptr;
	float*                offset_t = nullptr;

	// --------------------------------------------------------------------------
	// Test setting
	// --------------------------------------------------------------------------
	//load_test_setting(cg, ig, image, na, betas, nfoot, alut, hlut);

	//useMat = true;
	//load_cuboid_setting( "./test/cuboid/center/center.mat",
	//	cg, ig, image, na, betas, nfoot,
	//	"./test/cuboid/alut.mat", alut, 
	//	"./test/cuboid/hlut.mat", hlut
	//);
	//nfoot = 5;
	//outFn = "./test/cuboid/center/proj";

	// --------------------------------------------------------------------------
	// Command line
	// --------------------------------------------------------------------------
#ifdef USE_CMDARG
	{
		for (int i=1; i<argc; ++i)
		{
			std::string cmd = argv[i];
			if (cmd[0] != '-')
			{
				Usage();
			}

			if (cmd == "-r" || cmd == "--root")
			{
				root = argv[++i]; 
				if (root.back() != '/')
					root.append("/");
			}
			else if (cmd == "-i" || cmd == "--input")
			{
				matFn = argv[++i];
			}
			else if (cmd == "-m" || cmd == "--method")
			{
				switch(atoi(argv[++i]))
				{
				case 0: method = cbctProjector::LUT_LUT; methodName = "LL"; break;
				case 1: method = cbctProjector::LUT_REG; methodName = "LR"; break;
				case 2: method = cbctProjector::LUT_DD;  methodName = "LD"; break;
				case 3: method = cbctProjector::SF_TT;   methodName = "TT"; break;
				case 4: method = cbctProjector::SF_TR;   methodName = "TR"; break;
				default:
					Usage();
				}
			}
			else if (cmd == "-o" || cmd == "--output")
			{
				outFn = argv[++i];
			}
			else if (cmd == "-d" || cmd == "--order")
			{
				switch (atoi(argv[++i]))
				{
				case 0: order = cbctProjector::TS; break;
				case 1: order = cbctProjector::ST; break;
				default:
					Usage();
				}
			}
			else if (cmd == "-t" || cmd == "--time")
			{
				measureTime = atoi(argv[++i]) > 0;
			}
			else if (cmd == "-s" || cmd == "--scale")
			{
				doScaling = atoi(argv[++i]) == 1;
			}
			else if (cmd == "-a" || cmd == "-alut")
			{
				alutFn = argv[++i];
			}
			else if (cmd == "-v" || cmd == "-hlut")
			{
				hlutFn = argv[++i];
			}
			else if (cmd == "-ost" || cmd == "--offset_st")
			{
				useOffset_st = atoi(argv[++i]) > 0;
			}
			else if (cmd == "-h" || cmd == "--help")
			{
				Usage();
			}
			else
			{
				Usage();
			}
		}

		if (root.empty() || matFn.empty())
		{
			Usage();
		}

		if (outFn.empty())
		{
			outFn = root + "fp_" + methodName + ".mat";
		}
		else
		{
			if (outFn == "no")
				makeOutput = false;
			outFn = root + outFn;
		}

		alutFn = root + alutFn;
		hlutFn = root + hlutFn;

		if (method != cbctProjector::LUT_LUT)
		{
			hlutFn.clear();
		}
		if (method == cbctProjector::SF_TR || method == cbctProjector::SF_TT)
		{
			alutFn.clear();
			hlutFn.clear();
		}

	}
#else
	{
		root = "./test/cuboid_tmi/";
		matFn = "outcen/coutcen.mat";
		method = cbctProjector::SF_TR;
		outFn = "outcen/proj";
	}
#endif


	printf("\n\
		Input arguments:\n\
		- root: %s\n\
		- input: %s\n\
		- method: %s\n\
		- output: %s\n\
		- time: %s\n\
		- alut: %s\n\
		- hlut: %s\n\
		- offset: %s\n",
		root.c_str(), matFn.c_str(), methodName.c_str(), outFn.c_str(), measureTime ? "yes": "no",
		alutFn.c_str(), hlutFn.c_str(),
		useOffset_st ? "dynamic": "static");


	load_cuboid_setting(root + matFn,
		cg, ig, image, na, betas, nfoot,
		alutFn, alut,
		hlutFn, hlut
	);
	nfoot = 5;

	if (useOffset_st)
	{
		load_offset_st(root + matFn, offset_s, offset_t);
	}

	if (measureTime)
	{
		cudaEventCreate(&start);
		cudaEventCreate(&stop);
		elapsedTime = 0.0f;

		doScaling = false;
	}

	// --------------------------------------------------------------------------
	// Forward projection
	// --------------------------------------------------------------------------
	gpuname = CUDA_SET_BEST_DEVICE();
	//cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);

	cbct::initCg(cg);
	cbct::initIg(ig);
	//cbct::showGeometry();

	sys.setGeometry(&cg, &ig, method, order);
	sys.setFootsize(nfoot, &alut, &hlut);

	//CUDA_CHECK_ERROR("after setFootsize");

	//printf("%d, %d\n", cg.nxy.x, cg.nxy.y);

	pitch = 
	mm.allocLinearPitch<float>("proj", make_int2(cg.nxy.y, cg.nxy.x*na));        // [nt, ns]
	mm.alloc3DArray<float>("image", make_int3(ig.nxyz.z, ig.nxyz.x, ig.nxyz.y)); // [nz, nx, ny]
	mm.copyTo3DArray("image", image);

	cbct::bindImgToTex3D(mm.getArray("image"));

	for (int ia = 0; ia<na; ++ia)
	{
		//ia = 9;

		const float beta = betas[ia];
		const float   cs = cosf(beta);
		const float   sn = sinf(beta);

		const float2 uv_t = rotateCCW_z(make_float2(0, -1), cs, sn); // unit vector of ray
		const float2 uv_s = rotateCCW_z(make_float2(1, 0), cs, sn);
		const float2 src = rotateCCW_z(make_float2(0, cg.dso), cs, sn);
		const float2 det = rotateCCW_z(make_float2(0, -(cg.dsd - cg.dso)), cs, sn);

		if (useOffset_st)
		{
			cg.offset = make_float2(offset_s[ia], offset_t[ia]);
			cg.wxy = make_float2
			(
				0.5f*(cg.nxy.x - 1) + cg.offset.x,
				0.5f*(cg.nxy.y - 1) + cg.offset.y
			);
			cbct::initCg(cg);
			sys.setGeometry(&cg, &ig, method, order);
		}

		sys.setViews(uv_s, uv_t, src, det, beta*float(RAD2DEG));
		sys.initFootprint();

		//cudaThreadSynchronize();
		if (measureTime) cudaEventRecord(start);
		sys.Ax(mm.getLinear<float>("proj") + ia*pitch*cg.nxy.x, pitch, doScaling);
		//sys.Ax(mm.getLinear<float>("proj") + ia*pitch*cg.nxy.x, pitch, false);
		if (measureTime) cudaEventRecord(stop);

		if (measureTime)
		{
			float ms = 0.0f;
			cudaEventSynchronize(stop);
			cudaEventElapsedTime(&ms, start, stop);
			elapsedTime += ms;
		}
		//if (ia == 20) break;
		//break;
	}

	if (!doScaling && !measureTime)
	{
		for (int ia = 0; ia < na; ++ia)
		{
			const float beta = betas[ia];
			sys.scaleProjection(
				mm.getLinear<float>("proj") + ia*pitch*cg.nxy.x,
				pitch, beta
			);
		}
	}


	// --------------------------------------------------------------------------
	// Output
	// --------------------------------------------------------------------------
	if (makeOutput)
	{
		proj = new float[cg.nxy.x*cg.nxy.y*na];
		mm.copyLinearPitchToHost("proj", proj);

		mexMat::write(outFn, "proj", &proj[0], cg.nxy.y, cg.nxy.x, na);
	}

	if (measureTime)
	{
		std::string fn = root + "fp_time_" + methodName + ".txt";
		std::ofstream fout;

		fout.open(fn.c_str());
		fout << "FP with " << methodName << std::endl;;
		fout << "Using " << matFn << " , " << gpuname << std::endl;
		fout << "# projections: " << na << std::endl;
		fout << "Total time: " << elapsedTime << " msec" << std::endl;
		fout << "Average time: " << elapsedTime / float(na) << " msec" << std::endl;
		fout << "GUPS: " << GUPS(ig.nxyz, na, elapsedTime / 1000.0f) << std::endl;
		fout.close();

		std::cout << std::endl;
		std::cout << "FP with " << methodName << std::endl;;
		std::cout << "Using " << matFn << " , " << gpuname << std::endl;
		std::cout << "# projections: " << na << std::endl;
		std::cout << "Total time: " << elapsedTime << " msec" << std::endl;
		std::cout << "Average time: " << elapsedTime / float(na) << " msec" << std::endl;
		std::cout << "GUPS: " << GUPS(ig.nxyz, na, elapsedTime / 1000.0f) << std::endl;
		std::cout << std::endl;

		cudaEventDestroy(start);
		cudaEventDestroy(stop);
	}

	// --------------------------------------------------------------------------
	// Clear
	// --------------------------------------------------------------------------
	if (proj) delete[] proj;
	if (useMat)
	{
		if (image) mxFree(image);
		if (betas) mxFree(betas);
	}
	else
	{
		if (image) delete[] image;
		if (betas) delete[] betas;
	}

	if (method == cbctProjector::LUT_LUT)
	{
		mxFree(alut.lut);
		mxFree(hlut.lut);
	}
	if (method == cbctProjector::LUT_REG || method == cbctProjector::LUT_DD)
	{
		mxFree(alut.lut);
	}

	if (useOffset_st)
	{
		if (offset_s) mxFree(offset_s);
		if (offset_t) mxFree(offset_t);
	}

	mm.freeAll();
	sys.freeAll();
	cudaDeviceReset();

	return 0;
}



