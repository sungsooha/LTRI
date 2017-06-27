#include "cbctProjector.h"
#include "mexMat.h"


cbctProjector::cbctProjector() 
: p_cg(nullptr), p_ig(nullptr), _method(SF_TR), _order(TS), _nfoot(5),
  _nlines(0), _nplanes(0)
{
}


cbctProjector::~cbctProjector()
{
}

void cbctProjector::freeAll()
{
	_mm.freeAll();
}

void cbctProjector::setGeometry(const cbct_cg* cg, const cbct_ig* ig, method m, order o, bool useamptex)
{
	p_cg = cg;
	p_ig = ig;
	_method = m;
	_order = o;

	_useamptex = useamptex;

	if (_method == LUT_LUT || _method == LUT_REG || _method == LUT_DD)
	{
		_cmpGeoDiv();
	}

	if ( (_method == SF_TR || _method == SF_TT) && _useamptex)
	{
		_mm.allocLinearPitch<float>("amp", make_int2(p_cg->nxy.x, p_cg->nxy.y));
		_mm.allocArray<float>("amparr", make_int2(p_cg->nxy.x, p_cg->nxy.y));
		_bindAmp1();
	}
}

void cbctProjector::setFootsize(const int nfoot, const lut_area* alut, const lut_height* hlut)
{
	_nfoot = nfoot;

	switch(_method)
	{
	case LUT_LUT: case LUT_REG: _mm.allocLinear<float>("footinfo", p_ig->nxyz.x*p_ig->nxyz.y*(3 + _nfoot + 2)); break;
	case LUT_DD:  _mm.allocLinear<float>("footinfo", p_ig->nxyz.x*p_ig->nxyz.y*(3 + _nfoot + 1)); break;
	case SF_TR:   _mm.allocLinear<float>("footinfo", p_ig->nxyz.x*p_ig->nxyz.y*(2 + _nfoot + 0)); break;
	case SF_TT:   _mm.allocLinear<float>("footinfo", p_ig->nxyz.x*p_ig->nxyz.y*(2 + _nfoot + 4)); break;
	default:
		;
	}

	if ((int)_method < 3 && alut)
	{
		//_bindAreaLut(alut);
		_bindAreaLut_pos(alut);
	}

	if (_method == LUT_LUT && hlut)
	{
		//_bindHeightLut(hlut);
		_bindHeightLut_pos(hlut);
	}
}

void cbctProjector::setViews(const float2 uv_s, const float2 uv_t, const float2 src, const float2 det, const float deg)
{
	_updateRays(uv_s, uv_t, src, det, deg);

	if (_useamptex && (_method==SF_TR || _method==SF_TT))
	{
		_updateAmp();
	}
}

void cbctProjector::scaleProjection(float* d_proj, const int pitch, const float beta) const
{
	_scale_projection(d_proj, pitch, beta);
}

void cbctProjector::initFootprint()
{
	//_init_footprint();
	_init_footprint_pos();

	//std::vector<float> footinfo(p_ig->nxyz.x*p_ig->nxyz.y*(_nfoot + 5),0);
	//_mm.copyLinearToHost("footinfo", &footinfo[0]);
	//mexMat::write("footinfo_LL.mat", "footinfo", &footinfo[0], p_ig->nxyz.x, p_ig->nxyz.y, _nfoot + 5);
}

void cbctProjector::Ax(float* d_proj, const int pitch, const bool doScaling, const float val)
{
	_do_projection(d_proj, pitch, val);

	if (doScaling)
		_scale_projection(d_proj, pitch); 
}

void cbctProjector::At(float* d_image, const int pitch, const int pid, const float val)
{
	_do_backprojection(d_image, pitch, pid, val);
}

void cbctProjector::_cmpGeoDiv()
{
	const float    dsd2 = SQR(p_cg->dsd);
	const float    dsd  = p_cg->dsd;
	const int2     nst  = p_cg->nxy;
	const float2   dst  = p_cg->dxy;
	const float2   wst  = p_cg->wxy;
	const float2   dst2 = 0.5f*dst;

	const int n_sample = 16;
	const float delta = dst.y / float(n_sample);
	const float delta2 = 0.5f*delta;

	std::vector<float> aziang(nst.x), geodiv(nst.x*nst.y);// , upward(nst.y);
	bool baziang = false;

	//for (int j=0; j<nst.y; ++j)
	//{
	//	const float t = (j - wst.y - 0.5f)*dst.y + delta2;

	//	float sum = 0.0f;
	//	for (int sample_i=0; sample_i<n_sample; ++sample_i)
	//	{
	//		const float tt = t + sample_i*delta;
	//		sum += dsd / sqrtf(SQR(tt) + dsd2);
	//	}
	//	upward[j] = sum / float(n_sample);
	//}

	//float l_theta = atanf(t / (sqrtf(SQR(s) + SQR(g_cg.dsd))));
	//l_theta = 1.0f / fabsf(cosf(l_theta));

	//float l_phi = atan2f(s, g_cg.dsd) + c_viewang*float(DEG2RAD);
	//l_phi = g_ig.dxyz.x / fmaxf(fabsf(cosf(l_phi)), fabsf(sinf(l_phi)));

	for (int j=0; j<nst.y; ++j)
	{
		const float t = (j - wst.y)*dst.y;
		
		const float upward = 1.0f;// dsd / sqrtf(SQR(t) + dsd2);
		const float polang = fabsf(atanf((t - dst2.y) / dsd) - atanf((t + dst2.y) / dsd));
		
		for (int i=0; i<nst.x; ++i)
		{
			const float s = (i - wst.x)*dst.x;
			float l_theta = atanf(t / (sqrtf(SQR(s) + dsd2)));
			l_theta = 1.0f / fabsf(cosf(l_theta));
			float l_phi = atan2f(s, dsd);// +c_viewang*float(DEG2RAD);
			l_phi = 1.0f/ fmaxf(fabsf(cosf(l_phi)), fabsf(sinf(l_phi)));

			//const float L = sqrtf(SQR(s) + dsd2);
			//const float upward = L / sqrtf(SQR(L) + SQR(t));
			if (!baziang)
			{	
				aziang[i] = fabsf(atanf((s - dst2.x) / dsd) - atanf((s + dst2.x) / dsd));

				// trick.. why?
				//const float sc = dsd / sqrtf(SQR(s) + dsd2);
				//aziang[i] *= sc;
			}
			geodiv[j*nst.x + i] = (1.0f / (upward*polang*aziang[i]))*(l_theta*l_phi);
		}
		baziang = true;
	}

	//mexMat::write("./test/debug/geodiv.mat", "geodiv", &geodiv[0], nst.x, nst.y, 1);

	_mm.allocArray<float>("geodiv",nst);
	_mm.copyToArray("geodiv", &geodiv[0]);
	_bindGeoDiv();
}



