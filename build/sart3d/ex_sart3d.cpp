#include "cbctProjector.h"
#include "mexMat.h"
#include "mexArg.h"
#include <algorithm>
#include <map>
#include <iostream>

#define USE_CMDARG

extern "C"
void sart3d_correct(
	float* d_proj,
	float* d_measurement,
	float* d_norm,
	const int pitch, const int2 nxy,
	const float beta,
	const int method
);

extern "C"
void sart3d_update(
	float* d_recon,
	const float* d_back,
	const float* d_norm,
	const int pitch, const int3 nxyz, const float lambda
);


static void Usage()
{
	printf("\n\
		Back-projection\n\
		-r || --root : root (or working) directory path\n\
		-i || --input : input mat filename (including extension, relative to root directory)\n\
		-m || --method : integer number (0: LUT_LUT, 1: LUT_REG, 2: SF_TT, 3: SF_TR)\n\
		-o || --output [optional]: output mat filename (including extension, relative to root directory, 'no' will do nothing)\n\
		-a || --alut: area look-up table mat file name with full path relative to root\n\
		-v || --hlut: height look-up table mat file name with full path relative to root \n\
		-T || --temp: interval to save intermediate recon result\n\
		-I || --iter: the max. number of iterations\n\
		-L || --lambda: relaxation parameter [default: 0.0001]\n\
		-S || --stop: stop value in HU (only valid when ref. is available)\n\
		-b || --verbose: 0: none, 1: level 1, 2: level 2\n\
		-h || --help : show this message\n\
		\n\
		Note: for the LUT_LUT and LUT_REG, it expects look-up table(s) is in lut sub-directory of the root.\n");
	exit(-1);
}

// shepp-logan pantom
void load_slp_setting(
	std::string matfn,
	cbct_cg& cg,
	cbct_ig& ig,
	float*&  proj,
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

	mexMat::fetch_proj(matfn, proj);
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

// assuming 'na' projections collected over 360 degrees uniformly
std::vector<int> get_recon_order(const float* betas, const int na)
{
	std::vector<int> order(na);
	//std::vector<bool> check(na, false);

	std::map<float, std::vector<int>> tb;

	for (int ia = 0; ia < na; ++ia)
	{
		float deg = fmodf(betas[ia] * float(RAD2DEG), 90.0f);
		if (tb.find(deg) == tb.end())
		{
			std::vector<int> orderset{ ia };
			tb[deg] = orderset;
		}
		else
		{
			tb[deg].push_back(ia);
		}
	}

	int wpos = 0;
	for (auto it=tb.begin(); it!=tb.end(); ++it)
	{
		std::vector<int> orderset = it->second;
		for (int i=0; i<orderset.size(); ++i)
		{
			//float deg = betas[orderset[i]] * float(RAD2DEG);
			//printf("%f ", deg);
			order[wpos] = orderset[i];
			wpos++;
		}
		//printf("\n");
	}

	return order;

}

int main(int argc, char** argv)
{
	cuMemMgr                       mm;

	cbctProjector                  sys;
	cbctProjector::method          method = cbctProjector::SF_TT;
	std::string                    methodName = "TR";
	cbctProjector::order           order = cbctProjector::TS;

	cbct_ig                        ig;
	cbct_cg                        cg;

	float *                        proj = nullptr;
	int                            proj_pitch;

	float *                        image = nullptr;
	int                            image_pitch;

	float *                        betas = nullptr;
	int                            na;

	int                            nfoot;

	lut_area                       alut;
	lut_height                     hlut;

	std::string                    root;
	std::string                    matFn;
	std::string                    outFn;
	std::string                    gpuname;

	std::string                    alutFn;
	std::string                    hlutFn;

	bool                           useMat = true;
	bool                           measureTime = false;
	bool                           makeOutput = true;


	int                            ImageInterval = 0;

	std::vector<float>             metric;
	float *                        ref_image = nullptr;

	int                            maxIter = 2000;
	float                          lambda = 0.0001f;
	float                          stopHU = 75.0f;

	int                            verbose = 0;

	// --------------------------------------------------------------------------
	// Test setting (command line argument)
	// --------------------------------------------------------------------------
#ifdef USE_CMDARG
	{
		for (int i = 1; i<argc; ++i)
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
				switch (atoi(argv[++i]))
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
			else if (cmd == "-T" || cmd == "--temp")
			{
				ImageInterval = atoi(argv[++i]);
			}
			else if (cmd == "-a" || cmd == "-alut")
			{
				alutFn = argv[++i];
			}
			else if (cmd == "-v" || cmd == "-hlut")
			{
				hlutFn = argv[++i];
			}
			else if (cmd == "-I" || cmd == "--iter")
			{
				maxIter = atoi(argv[++i]);
			}
			else if (cmd == "-L" || cmd == "--lambda")
			{
				lambda = (float)atof(argv[++i]);
			}
			else if (cmd == "-S" || cmd == "--stop")
			{
				stopHU = (float)atof(argv[++i]);
			}
			else if (cmd == "-b" || cmd == "--verbose")
			{
				verbose = atoi(argv[++i]);
			}
			else if (cmd == "-h" || cmd == "--help")
			{
				Usage();
			}
			else
			{
				printf("Unknown argument: %s\n", cmd.c_str());
				Usage();
			}
		}

		if (root.empty() || matFn.empty())
		{
			Usage();
		}

		if (outFn.empty())
		{
			outFn = root + "recon_" + methodName + ".mat";
		}
		else
		{
			if (outFn == "no")
				makeOutput = false;
			outFn = root + outFn;
		}

		mexMat::fetch_image(root + matFn, ref_image);
		if (ref_image == nullptr)
			ImageInterval = 0;

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
		- alut: %s\n\
		- hlut: %s\n\
		- max. # iter: %d\n\
		- lambda: %f\n\
		- stop HU: %f\n\
		- intermediate: %s (%d)\n",
		root.c_str(), matFn.c_str(), methodName.c_str(), outFn.c_str(), alutFn.c_str(), hlutFn.c_str(),
		maxIter, lambda, stopHU, ImageInterval>0 ? "yes" : "no", ImageInterval);

	if (verbose) printf("Load data ... ");
	load_slp_setting(root + matFn,
		cg, ig, proj, na, betas, nfoot,
		alutFn, alut,
		hlutFn, hlut
	);
	nfoot = 5;
	if (verbose) printf("done.\n");
	

	// --------------------------------------------------------------------------
	// SART3D
	// --------------------------------------------------------------------------
	gpuname = CUDA_SET_BEST_DEVICE();
	if (verbose) printf("Run on %s.\n", gpuname.c_str());

	cbct::initIg(ig);
	cbct::initCg(cg);
	//cbct::showGeometry();

	sys.setGeometry(&cg, &ig, method, order, true);
	sys.setFootsize(nfoot, &alut, &hlut);

	image = new float[ig.nxyz.x*ig.nxyz.y*ig.nxyz.z];

	image_pitch =
	mm.allocLinearPitch<float>("recon", make_int2(ig.nxyz.z, ig.nxyz.x*ig.nxyz.y));
	mm.allocLinearPitch<float>("bp_norm", make_int2(ig.nxyz.z, ig.nxyz.x*ig.nxyz.y));
	mm.allocLinearPitch<float>("cur_image", make_int2(ig.nxyz.z, ig.nxyz.x*ig.nxyz.y));
	mm.alloc3DArray<float>("arr_recon", make_int3(ig.nxyz.z, ig.nxyz.x, ig.nxyz.y));

	proj_pitch =
	mm.allocLinearPitch<float>("proj", make_int2(cg.nxy.y, cg.nxy.x*na));
	mm.allocLinearPitch<float>("fp_norm", make_int2(cg.nxy.y, cg.nxy.x));
	mm.allocLinearPitch<float>("cur_proj", make_int2(cg.nxy.y, cg.nxy.x));
	mm.alloc3DArrayLayer<float>("arr_proj", make_int3(cg.nxy.y, cg.nxy.x, 1));

	//mm.show();

	mm.copyToLinearPitch("proj", proj);
	
	cbct::bindImgToTex3D(mm.getArray("arr_recon"));
	cbct::bindProjToTexLayered(mm.getArray("arr_proj"));

	if (verbose) printf("[SART] Reconstructing\n\n");
	for (int iter=0; iter<maxIter; ++iter)
	{
		if (verbose) printf("\n[%d iteration]\n", iter+1);

		for (int ia=0; ia<na; ++ia)
		{
			const float beta = betas[ia];

			const float cs = cosf(beta);
			const float sn = sinf(beta);

			const float2 uv_t = rotateCCW_z(make_float2(0, -1), cs, sn); // unit vector of ray
			const float2 uv_s = rotateCCW_z(make_float2(1, 0), cs, sn);
			const float2 src = rotateCCW_z(make_float2(0, cg.dso), cs, sn);
			const float2 det = rotateCCW_z(make_float2(0, -(cg.dsd - cg.dso)), cs, sn);

			sys.setViews(uv_s, uv_t, src, det, beta*float(RAD2DEG));
			sys.initFootprint();

			mm.copyTo3DArray<float>("recon", "arr_recon"); 

			mm.memsetLinearPitch<float>("cur_proj");
			mm.memsetLinearPitch<float>("fp_norm");

			sys.Ax(mm.getLinear<float>("cur_proj"), proj_pitch, false);
			sys.Ax(mm.getLinear<float>("fp_norm"), proj_pitch, false, 1.0f);

			sart3d_correct(
				mm.getLinear<float>("cur_proj"),
				mm.getLinear<float>("proj") + ia*proj_pitch*cg.nxy.x,
				mm.getLinear<float>("fp_norm"),
				proj_pitch, cg.nxy, beta, int(method)
			);

			//{
			//	std::vector<float> pp(cg.nxy.x*cg.nxy.y);
			//	mm.copyLinearPitchToHost("cur_proj", &pp[0]);
			//	mexMat::write("./test/debug/amp.mat", "amp", &pp[0], cg.nxy.y, cg.nxy.x, 1);
			//}

			mm.copyTo3DArray<float>("cur_proj", "arr_proj"); // check

			mm.memsetLinearPitch<float>("cur_image");
			mm.memsetLinearPitch<float>("bp_norm");

			sys.At(mm.getLinear<float>("cur_image"), image_pitch, 0);
			sys.At(mm.getLinear<float>("bp_norm"), image_pitch, 0, 1.0f);

			//{
			//	std::vector<float> bp_norm(ig.nxyz.x*ig.nxyz.y*ig.nxyz.z);
			//	mm.copyLinearPitchToHost("bp_norm", &bp_norm[0]);
			//	mexMat::write("./test/debug/bp_norm.mat", "norm", &bp_norm[0], ig.nxyz.z, ig.nxyz.x, ig.nxyz.y);

			//	//mm.copyLinearPitchToHost("cur_image", &bp_norm[0]);
			//	//mexMat::write("./test/debug/cur_image.mat", "image", &bp_norm[0], ig.nxyz.z, ig.nxyz.x, ig.nxyz.y);
			//}

			sart3d_update(
				mm.getLinear<float>("recon"),
				mm.getLinear<float>("cur_image"),
				mm.getLinear<float>("bp_norm"),
				image_pitch, ig.nxyz, lambda
			);

			if (verbose == 2 && (ia == 0 || (ia + 1) % 100 == 0 || ia + 1 == na))
				printf("%d ", ia + 1);
			//break;

		} // ia

		if (verbose) printf("done. ");

		if (ref_image)
		{
			//printf("error in bmetric\n");
			float me = 0;
			mm.copyLinearPitchToHost("recon", image);

			int iz = ig.nxyz.z / 2;
			for (int iy=0; iy<ig.nxyz.y; ++iy)
				for (int ix=0; ix<ig.nxyz.x; ++ix)
				{
					int idx = iz + ig.nxyz.z*(ix + iy*ig.nxyz.x);
					me += SQR(1000.0f*(ref_image[idx] - image[idx]));
				}
			me = sqrtf(me) / float(ig.nxyz.x*ig.nxyz.y);

			//for (int i = 0; i < ig.nxyz.x*ig.nxyz.y*ig.nxyz.z; ++i)
			//	me += SQR(1000.0f*(ref_image[i] - image[i]));
			//me = sqrtf(me) / float(ig.nxyz.x*ig.nxyz.y*ig.nxyz.z);

			metric.push_back(me);
			if (verbose) printf("%f", 1000.0f*me);

			if (iter && metric[iter-1] <= me)
			{
				printf("\nTerminated as it is diverging.\n");
				break;
			}

			if ( 1000.0f*me <= stopHU)
			{
				printf("\nTerminated as it reaches stop HU (%f).\n", me * 1000);
				break;
			}
		}

		if (iter + 1 == maxIter)
		{
			printf("\nTerminated as it reaches max. # iterations.\n");
			break;
		}

		if (verbose) printf("\n");


		
		if (iter && ImageInterval &&  (iter + 1)<maxIter &&  (iter+1)%ImageInterval == 0)
		{
			if (verbose) printf("output so far ... ");
			mm.copyLinearPitchToHost("recon", image);
			std::string tmpFn = outFn.substr(0, outFn.size()-4) + "_" + std::to_string(iter+1) + ".mat";

			//std::cout << tmpFn.c_str() << std::endl;

			std::vector<std::string> dataname{ "recon" };
			std::vector<float*> data{ image };
			std::vector<int3> dim{ make_int3(ig.nxyz.z, ig.nxyz.x, ig.nxyz.y) };

			if (!metric.empty())
			{
				dataname.push_back("metric");
				data.push_back(&metric[0]);
				dim.push_back(make_int3(metric.size(), 1, 1));
			}

			mexMat::write(tmpFn, dataname, data, dim);
			if (verbose) printf("done.\n");
		}
		
		//printf("next iter\n");
	} // iter
	if (verbose) printf("Reconstruction is done.\n");

	// --------------------------------------------------------------------------
	// Output
	// --------------------------------------------------------------------------

	if (verbose) printf("Output ... ");
	//printf("check\n");
	mm.copyLinearPitchToHost("recon", image);

	//printf("check\n");
	//mexMat::write(outFn, "recon", image, ig.nxyz.z, ig.nxyz.x, ig.nxyz.y);
	{
		std::vector<std::string> dataname{"recon"};
		std::vector<float*> data{image};
		std::vector<int3> dim{make_int3(ig.nxyz.z, ig.nxyz.x, ig.nxyz.y)};

		if (!metric.empty())
		{
			dataname.push_back("metric");
			data.push_back(&metric[0]);
			dim.push_back(make_int3(metric.size(), 1, 1));
		}

		mexMat::write(outFn, dataname, data, dim);
	}
	if (verbose) printf("done.\n");

	// --------------------------------------------------------------------------
	// Clear
	// --------------------------------------------------------------------------
	if (verbose) printf("clean ... ");
	if (image) delete[] image;
	if (useMat)
	{
		if (proj) mxFree(proj);
		if (betas) mxFree(betas);
		if (ref_image) mxFree(ref_image);
	}
	else
	{
		if (proj) delete[] proj;
		if (betas) delete[] betas;
		if (ref_image) delete[] ref_image;
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

	mm.freeAll();
	sys.freeAll();
	cudaDeviceReset();

	if (verbose) printf("done. Goodbye~~~!!!\n");
	return 0;
}








