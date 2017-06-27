#pragma once
#include "def_cbct.h"

class cbctProjector
{
public:
	enum method{LUT_LUT = 0, LUT_REG = 1, LUT_DD = 2, SF_TT=3, SF_TR=4};
	enum order{TS=0, ST=1};

	/**
	 * \brief constructor
	 */
	cbctProjector();
	/**
	 * \brief destructor
	 */
	~cbctProjector();

	/**
	 * \brief free all device memories
	 * \note required to call cudaDeviceReset() at the end of program
	 */
	void freeAll();

	/**
	 * \brief set cone-beam geometry
	 * \param cg cone-beam geometry
	 * \param ig image geometry
	 * \param m projection method
	 * \param d projection direction
	 * \param o data ordering
	 */
	void setGeometry(const cbct_cg* cg, const cbct_ig* ig, method m, order o, bool useamptex = false);

	/**
	 * \brief set footprint size (and set lut for LUT)
	 * \param nfoot footprint size
	 * \param alut  area lut
	 * \param hlut  height lut
	 */
	void setFootsize(const int nfoot, const lut_area* alut=nullptr, const lut_height* hlut=nullptr);

	/**
	 * \brief set projection views
	 * \param uv_s unit vector in s-direction
	 * \param uv_t unit vector in t-direction
	 * \param src source position
	 * \param det the center of the detector position
	 * \param ang projection angle in degree
	 */
	void setViews(const float2 uv_s, const float2 uv_t, const float2 src, const float2 det, const float ang);

	/**
	 * \brief scaling projection data with amplitude for SF or with geometric spreading term for LUT
	 * \param d_proj projection array in device memory
	 * \param pitch pitch of the projection array
	 * \param beta projection view in radian
	 */
	void scaleProjection(float* d_proj, const int pitch, const float beta) const;

	void initFootprint();

	/**
	 * \brief Forward projection, Ax
	 * \param d_proj projection array in device memory
	 * \param pitch pitch of the projection array
	 * \param doScaling boolean to apply amplitude term for SF or geometric spreading term for LUT
	 * \param val constant image, must positive number
	 */
	void Ax(float* d_proj, const int pitch, const bool doScaling = false, const float val=-1.0f);

	/**
	 * \brief Back projection
	 * \param d_image image array in device memory
	 * \param pitch pitch of the image array
	 * \param pid projection ID to refer projection view
	 * \param val constant projection, must positive number
	 */
	void At(float* d_image, const int pitch, const int pid, const float val=-1.0f);

private:

	/**
	 * \brief compute geometric spreading term and store in device array for LUT
	 */
	void _cmpGeoDiv();

	/**
	 * \brief bind geometric spreading term to texture (call by cmpGeoDiv())
	 */
	void _bindGeoDiv();
	
	void _bindAmp1();

	/**
	 * \brief bind area lookup table to texture
	 * \param alut area lookup table
	 */
	void _bindAreaLut(const lut_area* alut);
	void _bindAreaLut_pos(const lut_area* alut);

	/**
	 * \brief bind height lookup table to texture
	 * \param hlut height lookup table
	 */
	void _bindHeightLut(const lut_height* hlut);
	void _bindHeightLut_pos(const lut_height* hlut);

	/**
	 * \brief update cone-beam volumetric ray information at a view
	 * \param uv_s unit vector in s-direction
	 * \param uv_t unit vector in t-direction
	 * \param src source position
	 * \param det the center of the detector position
	 * \param ang projection angle in degree
	 * \note needed for LUT
	 */
	void _updateRays(const float2 uv_s, const float2 uv_t, const float2 src, const float2 det, float ang);

	void _updateAmp();

	/**
	 * \brief update ray in xy-plane, required for LUT
	 * \note call by _updateRays(...)
	 */
	void _updateLines();
	void _updateLines_pos();

	/**
	 * \brief update ray in z-direction, required for LUT
	 * \note call by _updateRays(...)
	 */
	void _updatePlanes();
	void _updatePlanes_pos();

	/**
	 * \brief pre-compute footprint in xy-plane (also other sharable information in z-direction)
	 */
	void _init_footprint();
	void _init_footprint_pos();

	/**
	 * \brief do projection
	 * \param d_proj projection array in device memory
	 * \param pitch pitch of the array
	 * \param val constant image, must positive number to activate
	 */
	void _do_projection(float* d_proj, const int pitch, const float val);

	void _do_backprojection(float* d_image, const int pitch, const int pid, const float val=-1.0f);

	/**
	 * \brief scale projection as post-processing
	 * \param d_proj projection array in device memory
	 * \param pitch pitch of the array
	 * \note scale with amplitude (A1) for SF and with geometric divergence for LUT
	 */
	void _scale_projection(float* d_proj, const int pitch) const; // amp or spreading

	void _scale_projection(float* d_proj, const int pitch, const float beta) const;


	const cbct_cg   *    p_cg;
	const cbct_ig   *    p_ig;

	cuMemMgr              _mm;

	method                _method;
	order                 _order;
	int                   _nfoot;

	int                   _nlines;
	int                   _nplanes;

	bool                  _useamptex;
};

