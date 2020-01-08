#include "chroma.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"
#include "meas/sources/z2_src.h"
#include "meas/sources/zN_src.h"

#include "gamma_ops.h"  // definition of gamma_ops() ftn
#include <iomanip>      // std::setw
#include <algorithm>
#include <fstream>

//constexpr int num_links = 1; // calculate links forward and backward displacement
//constexpr int link_dir = 2; // only z direction
//constexpr int num_disp = 8*num_links + 1;

// 1/sqrt(2)
#define INV_SQRT2 0.70710678118654752440


// Number of gamma operators (Ns * Ns = 16)
#define NUM_G     16

//#define CALC_ERR_ERR

#define POW2(a) ((a) * (a))
#define POW3(a) ((a) * (a) * (a))
#define POW4(a) ((a) * (a) * (a) * (a))

using namespace QDP;
using namespace Chroma;

// Anonymous namespace
namespace TSF {
	// Standard Time Slicery
	class TimeSliceFunc: public SetFunc {
	public:
		TimeSliceFunc(int dir) :
		dir_decay(dir) {
		}
		int operator()(const multi1d<int>& coordinate) const {
			return coordinate[dir_decay];
		}
		int numSubsets() const {
			return Layout::lattSize()[dir_decay];
		}
		int dir_decay;
	private:
		TimeSliceFunc() {
		}  // hide default constructor
	};
}

//====================================================================
// Structures for input parameters
//====================================================================

struct Inverter_t {
	GroupXML_t invParamHP;
	GroupXML_t invParamLP;
	GroupXML_t invParamLP_1;
	GroupXML_t fermactLP;
	GroupXML_t fermactHP;
	GroupXML_t fermactHPE;
	Real mass;
	bool setupParamLP_1;
};

struct Displacement_t {
	multi1d<int> link_dirs;
	int link_max;
};

struct NoiseSource_t {
	int version;
	
	int Max_Nr_LP;  // Minimum Nr_LP to terminate
	int Min_Nr_LP;  // Minimum Nr_LP to terminate
	
	int Nr_LPHP_ratio;
	
	bool dilute;
	multi1d<int> timeslices;
	
	std::string NoiseSrcType;
	
	// for version 2; restart if condition satisfied
	int restart_NrLP;  // Nr_LP to check if restart is needed
	double restart_factor;  // Restart criterion: err_Cr / err_LP > restart_factor
	GroupXML_t invParamLP_r;  // Restart LP loops with new inverter parameters
};

struct Checkpoint_t {
	int version;
	multi1d<double> checkpoints;
	std::string OutFileName;
	std::string LaMETOutFileName;
	
	// below variables are for the program process, not for read
	std::vector<double> cp;
	int chkout_order;
};

struct Params_t {
	multi1d<int> nrow;
	Inverter_t inverter;
	NoiseSource_t noise_src;
	Displacement_t displacement;
	multi1d<Checkpoint_t> checkpoint;
	bool use_HPE;
};

struct Inline_input_t {
	Params_t param;
	GroupXML_t cfg;
	QDP::Seed rng_seed;
};

//====================================================================
// Read input parameters
//====================================================================
void read(XMLReader& xml, const std::string& path, Inverter_t& p) {
	
	XMLReader paramtop(xml, path);
	p.invParamHP = readXMLGroup(paramtop, "InvertParamHP", "invType");
	p.invParamLP = readXMLGroup(paramtop, "InvertParamLP", "invType");
	read(paramtop, "Mass", p.mass);
	
	// Fermion action for Low Precision Inverter
	if (paramtop.count("FermionActionLP") > 0)
		p.fermactLP = readXMLGroup(paramtop, "FermionActionLP", "FermAct");
	else
		p.fermactLP = readXMLGroup(paramtop, "FermionAction", "FermAct");
	
	// Fermion action for Low Precision Inverter
	if (paramtop.count("FermionActionHP") > 0)
		p.fermactHP = readXMLGroup(paramtop, "FermionActionHP", "FermAct");
	else
		p.fermactHP = readXMLGroup(paramtop, "FermionAction", "FermAct");
	
	// Fermion action for Hopping Parameter Expansion; must be unpreconditioned
	p.fermactHPE = readXMLGroup(paramtop, "FermionActionHPE", "FermAct");
	
	// Inverter parameters for the first run of the inversions
	// For Multigrid inverter, the first inversion is setting up the inverter
	if (paramtop.count("InvertParamLP_1") > 0) {
		p.invParamLP_1 = readXMLGroup(paramtop, "InvertParamLP_1", "invType");
		p.setupParamLP_1 = true;
	} else {
		p.setupParamLP_1 = false;
	}
}

void read(XMLReader& xml, const std::string& path, Displacement_t& p) {
	
	XMLReader paramtop(xml, path);
	
	// Fermion action for Low Precision Inverter
	if (paramtop.count("link_dirs") > 0){
		read(paramtop, "link_dirs", p.link_dirs);
		for (int idir=0 ; idir<p.link_dirs.size() ; idir++){
			if(p.link_dirs[idir]<0 || p.link_dirs[idir]>7){
				QDPIO::cerr << "Error! Link direction should be within range [0,7]."<< std::endl;
				QDP_abort(1);
			}
		}
	}
	else{
		p.link_dirs.resize(8);
		for(int i=0;i<8;i++)
			p.link_dirs[i]=i;
	}
	
	if (paramtop.count("link_max") > 0){
		read(paramtop, "link_max", p.link_max);
		if(p.link_max < 0){
			QDPIO::cerr << "Error! Link length can't be negative."<< std::endl;
			QDP_abort(1);
		}
	}
	else
		p.link_max = 2;
}

void read(XMLReader& xml, const std::string& path, NoiseSource_t& p) {
	XMLReader paramtop(xml, path);
	
	// Read input style version
	if (paramtop.count("version") > 0) {
		read(paramtop, "version", p.version);
	} else {
		p.version = 1;
	}
	
	switch (p.version) {
		case 1:
			break;
			
		case 2:
			read(paramtop, "MaxNrLP", p.Max_Nr_LP);
			read(paramtop, "MinNrLP", p.Min_Nr_LP);
			read(paramtop, "RestartNrLP", p.restart_NrLP);
			read(paramtop, "RestartFactor", p.restart_factor);
			p.invParamLP_r = readXMLGroup(paramtop, "InvertParamLP_r", "invType");
			break;
			
		default:
			QDPIO::cerr << "Error! NoiseSource input parameter version of " << p.version
			<< " is not supported." << std::endl;
			QDP_abort(1);
	}
	
	read(paramtop, "NrLPHP_RATIO", p.Nr_LPHP_ratio);
	
	// Read Dilution information
	if (paramtop.count("Dilution") > 0) {
		p.dilute = true;
		read(paramtop, "Dilution/Timeslice", p.timeslices);
	} else {
		p.dilute = false;
	}
	
	// Read noise source type
	if (paramtop.count("NoiseSrcType") > 0) {
		read(paramtop, "NoiseSrcType", p.NoiseSrcType);
	} else {
		p.NoiseSrcType = "GAUSSIAN";
	}
}

void read(XMLReader& xml, const std::string& path, Checkpoint_t& p) {
	XMLReader paramtop(xml, path);
	
	// Read input style version
	if (paramtop.count("version") > 0) {
		read(paramtop, "version", p.version);
	} else {
		p.version = 1;
	}
	
	read(paramtop, "Checkpoints", p.checkpoints);
	read(paramtop, "OutFile", p.OutFileName);
	if(paramtop.count("LaMETOutFile") > 0){
		read(paramtop, "LaMETOutFile", p.LaMETOutFileName);
	} else {
		p.LaMETOutFileName=p.OutFileName+".LaMET";
	}
	
}

void read(XMLReader& xml, const std::string& path, Params_t& p) {
	XMLReader paramtop(xml, path);
	read(paramtop, "nrow", p.nrow);
	
	if (paramtop.count("UseHPE") > 0)
		read(paramtop, "UseHPE", p.use_HPE);
	else
		p.use_HPE = true;
	
	read(paramtop, "Inverter", p.inverter);
	read(paramtop, "NoiseSource", p.noise_src);
	read(paramtop, "Displacement", p.displacement);
	read(paramtop, "Checkpoint", p.checkpoint);
}

void read(XMLReader& xml, const std::string& path, Inline_input_t& p) {
	try {
		XMLReader paramtop(xml, path);
		
		read(paramtop, "Param", p.param);
		p.cfg = readXMLGroup(paramtop, "Cfg", "cfg_type");
		
		if (paramtop.count("RNG") > 0)
			read(paramtop, "RNG", p.rng_seed);
		else
			p.rng_seed = 11;     // default seed
	} catch (const std::string& e) {
		QDPIO::cerr << "Error reading XML : " << e << std::endl;
		QDP_abort(1);
	}
}

bool linkageHack(void) {
	bool success = true;
	
	success &= GaugeInitEnv::registerAll();
	success &= WilsonTypeFermActsEnv::registerAll();
	success &= LinOpSysSolverEnv::registerAll();
	
	return success;
}

//====================================================================
// Structures for error analysis
//====================================================================
struct ErrAnlyVars {
	multi3d<DComplex> TrM_inv_sum_LP;
	multi3d<DComplex> TrM_inv_sum_C;
	multi3d<double> TrM_inv_LP_sqsum_r;
	multi3d<double> TrM_inv_LP_sqsum_i;
	multi3d<double> TrM_inv_C_sqsum_r;
	multi3d<double> TrM_inv_C_sqsum_i;
	
#ifdef CALC_ERR_ERR
	multi3d<std::vector<DComplex>> TrM_inv_est;
#endif
	
	ErrAnlyVars() {
	}
	
	ErrAnlyVars(int num_disp,int NumTs) {
		TrM_inv_sum_LP.resize(num_disp, NUM_G, NumTs);
		TrM_inv_sum_C.resize(num_disp, NUM_G, NumTs);
		TrM_inv_LP_sqsum_r.resize(num_disp, NUM_G, NumTs);
		TrM_inv_LP_sqsum_i.resize(num_disp, NUM_G, NumTs);
		TrM_inv_C_sqsum_r.resize(num_disp, NUM_G, NumTs);
		TrM_inv_C_sqsum_i.resize(num_disp, NUM_G, NumTs);
		
#ifdef CALC_ERR_ERR
		TrM_inv_est.resize(num_disp, NUM_G, NumTs);
#endif
		
		for (int d = 0; d < num_disp; ++d)
			for (int i = 0; i < NUM_G; ++i)
				for (int t = 0; t < NumTs; ++t) {
					TrM_inv_sum_LP[d][i][t] = zero;
					TrM_inv_sum_C[d][i][t] = zero;
					
					TrM_inv_LP_sqsum_r[d][i][t] = 0.0;
					TrM_inv_LP_sqsum_i[d][i][t] = 0.0;
					TrM_inv_C_sqsum_r[d][i][t] = 0.0;
					TrM_inv_C_sqsum_i[d][i][t] = 0.0;
				}
	}
};

//====================================================================
// Calculate statistical error of scalar channel
// and return maximum/average error among timeslices
//====================================================================
void check_acc(int Nr_LP, int Nr_HP, ErrAnlyVars &errAnly, std::vector<int> &link_dirs,
			   int link_max, std::vector<int> &timeslices, double &ratio_C_LP_err,
			   double &m_err_max, double &m_err_av) {
	int NumTs = timeslices.size();
	int d = 0; // no displacement
	int g = 0; // scalar
	
	// TSM estimate of Tr [ M^{-1} \gamma ]
	multi1d<DComplex> TrM_inv_av(NumTs);
	
	// Low precision estimate and correction to Tr [ G^{-1} \Gamma ]
	multi1d<DComplex> TrM_inv_av_LP(NumTs);
	multi1d<DComplex> TrM_inv_av_C(NumTs);
	
	// Calculate average
	if (Nr_LP != 0) for (int t = 0; t < NumTs; ++t)
		TrM_inv_av_LP[t] = errAnly.TrM_inv_sum_LP[d][g][t] / (double)Nr_LP;
	
	if (Nr_HP != 0) for (int t = 0; t < NumTs; ++t)
		TrM_inv_av_C[t] = errAnly.TrM_inv_sum_C[d][g][t] / (double)Nr_HP;
	
	//-----------------------------
	// Calculate statistical error
	//-----------------------------
	std::vector<double> TrM_inv_LP_err_r(NumTs);
	std::vector<double> TrM_inv_C_err_r(NumTs);
	std::vector<double> TrM_inv_err_r(NumTs);
	
	double TrM_inv_err_r_av = 0.0;
	
	for (int t = 0; t < NumTs; ++t) {
		if (Nr_LP > 1) {
			double TrM_inv_av_LP_sq_r = pow(
											(double)TrM_inv_av_LP[t].elem().elem().elem().real(), 2);
			TrM_inv_LP_err_r[t] = sqrt(
									   (errAnly.TrM_inv_LP_sqsum_r[d][g][t] / Nr_LP - TrM_inv_av_LP_sq_r) / (Nr_LP
																											 - 1.0));
		} else {
			TrM_inv_LP_err_r[t] = 0.0;
		}
		
		if (Nr_HP > 1) {
			double TrM_inv_av_C_sq_r = pow(
										   TrM_inv_av_C[t].elem().elem().elem().real(), 2);
			TrM_inv_C_err_r[t] = sqrt(
									  (errAnly.TrM_inv_C_sqsum_r[d][g][t] / Nr_HP - TrM_inv_av_C_sq_r) / (Nr_HP
																										  - 1.0));
		} else {
			TrM_inv_C_err_r[t] = 0.0;
		}
		
		TrM_inv_err_r[t] = sqrt(
								pow(TrM_inv_LP_err_r[t], 2) + pow(TrM_inv_C_err_r[t], 2));
		TrM_inv_err_r_av += TrM_inv_err_r[t] / (double)NumTs;
	}  // for (int t=0; t<NumTs; ++t)
	
	// Find maximum error among all timeslices; C++11
	double max_LP_err = *std::max_element(TrM_inv_LP_err_r.begin(),
										  TrM_inv_LP_err_r.end());
	double max_C_err = *std::max_element(TrM_inv_C_err_r.begin(),
										 TrM_inv_C_err_r.end());
	double max_err = *std::max_element(TrM_inv_err_r.begin(),
									   TrM_inv_err_r.end());
	
	// Calculate ratio between error of correction term and error of LP term;
	// this will be used for the restart criterion
	if (max_LP_err != 0) ratio_C_LP_err = max_C_err / max_LP_err;
	
	m_err_max = max_err;
	m_err_av = TrM_inv_err_r_av;
	
	return;
}

void link_pattern(std::vector<int> &link_patterns, std::vector<int> &link_dirs, int link_max){
	link_patterns.push_back(0);
	if(link_max==1)
		for(int j=0;j<link_dirs.size();j++)
			link_patterns.push_back(10+link_dirs[j]);
	if(link_max==2)
		for(int i=0;i<link_dirs.size();i++){
			int dir1=link_dirs[i];
			link_patterns.push_back(10+dir1);
			for(int j=0;j<link_dirs.size();j++){
				int dir2=link_dirs[j];
				if(abs(dir1-dir2)!=4)
					link_patterns.push_back(200+10*dir1+dir2);
			}
		}
	if(link_max>=3){
		for(int i=0;i<8;i++){
			link_patterns.push_back(10+i);
			for(int j=0;j<8;j++){
				if(abs(i-j)!=4)
					link_patterns.push_back(200+10*i+j);
			}
		}
		for(int j=0;j<link_dirs.size();j++)
			for(int i=1;i<=link_max;i++)
				link_patterns.push_back(link_dirs[j]*100+i);
	}
}
//====================================================================
// Calculate statistical error, and save results
//====================================================================
void checkout(int Nr_LP, int Nr_HP, ErrAnlyVars &errAnly, std::string out_fname, std::string lamet_out_fname, std::vector<int> &link_dirs, int link_max,
			  std::vector<int> &timeslices, int chkout_order, bool Restarted) {
	int NumTs = timeslices.size();
	std::vector<int> link_patterns;
	link_pattern(link_patterns,link_dirs,link_max);
	int num_disp = link_patterns.size();
	int num_disp_mom= pow(8,2)+1;
	// TSM estimate of Tr [ M^{-1} \gamma ]
	multi3d<DComplex> TrM_inv_av(num_disp, NUM_G, NumTs);
	
	// Low precision estimate and correction to Tr [ G^{-1} \Gamma ]
	multi3d<DComplex> TrM_inv_av_LP(num_disp, NUM_G, NumTs);
	multi3d<DComplex> TrM_inv_av_C(num_disp, NUM_G, NumTs);
	
	// Calculate average
	if (Nr_LP != 0) for (int d = 0; d < num_disp; ++d)
		for (int g = 0; g < NUM_G; ++g)
			for (int t = 0; t < NumTs; ++t)
				TrM_inv_av_LP[d][g][t] = errAnly.TrM_inv_sum_LP[d][g][t] / (double)Nr_LP;
	
	if (Nr_HP != 0) for (int d = 0; d < num_disp; ++d)
		for (int g = 0; g < NUM_G; ++g)
			for (int t = 0; t < NumTs; ++t)
				TrM_inv_av_C[d][g][t] = errAnly.TrM_inv_sum_C[d][g][t] / (double)Nr_HP;
	
	//-----------------------------
	// Calculate statistical error
	//-----------------------------
	multi3d<double> TrM_inv_LP_err_r(num_disp, NUM_G, NumTs);
	multi3d<double> TrM_inv_LP_err_i(num_disp, NUM_G, NumTs);
	multi3d<double> TrM_inv_C_err_r(num_disp, NUM_G, NumTs);
	multi3d<double> TrM_inv_C_err_i(num_disp, NUM_G, NumTs);
	multi3d<double> TrM_inv_err_r(num_disp, NUM_G, NumTs);
	multi3d<double> TrM_inv_err_i(num_disp, NUM_G, NumTs);
	
	for (int d = 0; d < num_disp; ++d)
		for (int g = 0; g < NUM_G; ++g)
			for (int t = 0; t < NumTs; ++t) {
				if (Nr_LP > 1) {
					double TrM_inv_av_LP_sq_r = pow(
													(double)TrM_inv_av_LP[d][g][t].elem().elem().elem().real(), 2);
					double TrM_inv_av_LP_sq_i = pow(
													(double)TrM_inv_av_LP[d][g][t].elem().elem().elem().imag(), 2);
					TrM_inv_LP_err_r[d][g][t] =
					sqrt(
						 (errAnly.TrM_inv_LP_sqsum_r[d][g][t] / Nr_LP - TrM_inv_av_LP_sq_r) / (Nr_LP
																							   - 1.0));
					TrM_inv_LP_err_i[d][g][t] =
					sqrt(
						 (errAnly.TrM_inv_LP_sqsum_i[d][g][t] / Nr_LP - TrM_inv_av_LP_sq_i) / (Nr_LP
																							   - 1.0));
				} else {
					TrM_inv_LP_err_r[d][g][t] = 0.0;
					TrM_inv_LP_err_i[d][g][t] = 0.0;
				}
				
				if (Nr_HP > 1) {
					double TrM_inv_av_C_sq_r = pow(
												   TrM_inv_av_C[d][g][t].elem().elem().elem().real(), 2);
					double TrM_inv_av_C_sq_i = pow(
												   TrM_inv_av_C[d][g][t].elem().elem().elem().imag(), 2);
					TrM_inv_C_err_r[d][g][t] =
					sqrt(
						 (errAnly.TrM_inv_C_sqsum_r[d][g][t] / Nr_HP - TrM_inv_av_C_sq_r) / (Nr_HP
																							 - 1.0));
					TrM_inv_C_err_i[d][g][t] =
					sqrt(
						 (errAnly.TrM_inv_C_sqsum_i[d][g][t] / Nr_HP - TrM_inv_av_C_sq_i) / (Nr_HP
																							 - 1.0));
				} else {
					TrM_inv_C_err_r[d][g][t] = 0.0;
					TrM_inv_C_err_i[d][g][t] = 0.0;
				}
				
				TrM_inv_err_r[d][g][t] = sqrt(
											  pow(TrM_inv_LP_err_r[d][g][t], 2) + pow(TrM_inv_C_err_r[d][g][t], 2));
				TrM_inv_err_i[d][g][t] = sqrt(
											  pow(TrM_inv_LP_err_i[d][g][t], 2) + pow(TrM_inv_C_err_i[d][g][t], 2));
			}  // for (int d = 0; d < num_disp; ++d), for (int g=0; g<NUM_G; ++g), for (int t=0; t<NumTs; ++t)
	
	//------------------------------------------
	// Print LP estimate and correction term
	//------------------------------------------
	if (chkout_order == -1)  // -1 means final checkout
		QDPIO::cout << std::endl << "Checkout - fn; " << out_fname << std::endl;
	else
		QDPIO::cout << std::endl << "Checkout - " << chkout_order << "; "
		<< out_fname << std::endl;
	
	QDPIO::cout << std::endl;
	
//	for (int d = 0; d < num_disp; ++d) {
//		QDPIO::cout << "Displacement = " << link_patterns[d] << std::endl;
//		for (int t = 0; t < NumTs; ++t) {
//			QDPIO::cout << "Timeslice = " << timeslices[t] << std::endl;
//
//			for (int g = 0; g < NUM_G; ++g)
//				QDPIO::cout << "Tr [ M^{-1} g" << g << " ]_LP = "
//				<< TrM_inv_av_LP[d][g][t] << " err( " << TrM_inv_LP_err_r[d][g][t]
//				<< ", " << TrM_inv_LP_err_i[d][g][t] << " )" << std::endl;
//
//			QDPIO::cout << std::endl;
//
//			for (int g = 0; g < NUM_G; ++g)
//				QDPIO::cout << "Tr [ M^{-1} g" << g << " ]_Cr = "
//				<< TrM_inv_av_C[d][g][t] << " err( " << TrM_inv_C_err_r[d][g][t] << ", "
//				<< TrM_inv_C_err_i[d][g][t] << " )" << std::endl;
//
//			QDPIO::cout << std::endl;
//			QDPIO::cout << std::endl;
//		}  //for (int t=0; t<NumTs; ++t)
//	}  // for disp
	
	//--------------------------------------------------------------------
	// Calculate TSM estimate of TrM_inv
	//--------------------------------------------------------------------
	for (int d = 0; d < num_disp; ++d)
		for (int g = 0; g < NUM_G; ++g)
			for (int t = 0; t < NumTs; ++t)
				TrM_inv_av[d][g][t] = TrM_inv_av_LP[d][g][t] + TrM_inv_av_C[d][g][t];
	
	//--------------------------------------------------------------------
	// Calculate Another estimate of total error and error of error
	//--------------------------------------------------------------------
#ifdef CALC_ERR_ERR
	multi3d<double> TrM_inv_est_av_r(num_disp, NUM_G, NumTs);
	multi3d<double> TrM_inv_est_av_i(num_disp, NUM_G, NumTs);
	multi3d<double> TrM_inv_est_err_r(num_disp, NUM_G, NumTs);
	multi3d<double> TrM_inv_est_err_i(num_disp, NUM_G, NumTs);
	multi3d<double> TrM_inv_est_err_err_r(num_disp, NUM_G, NumTs);
	multi3d<double> TrM_inv_est_err_err_i(num_disp, NUM_G, NumTs);
	
	for (int d = 0; d < num_disp; ++d)
		for (int g=0; g<NUM_G; ++g)
			for (int t=0; t<NumTs; ++t) {
				TrM_inv_est_av_r[d][g][t] = 0.0;
				TrM_inv_est_av_i[d][g][t] = 0.0;
				
				int N_e = errAnly.TrM_inv_est[d][g][t].size();
				
				// Calculate average
				double tmp_av_r = 0.0;
				double tmp_av_i = 0.0;
				for (int i=0; i<N_e; ++i) {
					tmp_av_r += errAnly.TrM_inv_est[d][g][t][i].elem().elem().elem().real();
					tmp_av_i += errAnly.TrM_inv_est[d][g][t][i].elem().elem().elem().imag();
				}
				tmp_av_r /= N_e;
				tmp_av_i /= N_e;
				
				TrM_inv_est_av_r[d][g][t] = tmp_av_r;
				TrM_inv_est_av_i[d][g][t] = tmp_av_i;
				
				// Calculate error and error of error
				double tmp_p2sum_r = 0.0;
				double tmp_p2sum_i = 0.0;
				double tmp_p4sum_r = 0.0;
				double tmp_p4sum_i = 0.0;
				
				for (int i=0; i<N_e; ++i) {
					double m_r = errAnly.TrM_inv_est[d][g][t][i].elem().elem().elem().real();
					double m_i = errAnly.TrM_inv_est[d][g][t][i].elem().elem().elem().imag();
					tmp_p2sum_r += POW2(m_r - tmp_av_r);
					tmp_p2sum_i += POW2(m_i - tmp_av_i);
					tmp_p4sum_r += POW4(m_r - tmp_av_r);
					tmp_p4sum_i += POW4(m_i - tmp_av_i);
				}
				
				// error
				double tmp_err_r = sqrt(tmp_p2sum_r / (N_e*(N_e-1.0)));
				double tmp_err_i = sqrt(tmp_p2sum_i / (N_e*(N_e-1.0)));
				
				TrM_inv_est_err_r[d][g][t] = tmp_err_r;
				TrM_inv_est_err_i[d][g][t] = tmp_err_i;
				
				// error of error
				double tmp_err_var_r = sqrt( tmp_p4sum_r/(N_e*POW3(N_e-1.0)) - POW4(tmp_err_r)/(N_e-1.0) );
				double tmp_err_var_i = sqrt( tmp_p4sum_i/(N_e*POW3(N_e-1.0)) - POW4(tmp_err_i)/(N_e-1.0) );
				TrM_inv_est_err_err_r[d][g][t] = tmp_err_var_r / (2*tmp_err_r);
				TrM_inv_est_err_err_i[d][g][t] = tmp_err_var_i / (2*tmp_err_i);
			}  // Loop over d, g and t
#endif  // End of CALC_ERR_ERR
	
	//-----------------------------
	// Print results
	//-----------------------------
//	for (int d = 0; d < num_disp; ++d) {
//		QDPIO::cout << "Displacement = " << link_patterns[d] << std::endl;
//		for (int t = 0; t < NumTs; ++t) {
//			QDPIO::cout << "Timeslice = " << timeslices[t] << std::endl;
//
//			for (int g = 0; g < NUM_G; ++g)
//				QDPIO::cout << "Tr [ M^{-1} g" << g << " ]    = " << TrM_inv_av[d][g][t]
//				<< " std( " << TrM_inv_err_r[d][g][t] << ", "
//				<< TrM_inv_err_i[d][g][t] << " )" << std::endl;
//
//#ifdef CALC_ERR_ERR
//			for (int g=0; g<NUM_G; ++g)
//				QDPIO::cout << "Tr [ M^{-1} g" << g << " ]_est= ( "
//				<< TrM_inv_est_av_r[d][g][t] << ", "
//				<< TrM_inv_est_av_i[d][g][t] << " ) "
//				<< " ( " << TrM_inv_est_err_r[d][g][t] << ", "
//				<< TrM_inv_est_err_i[d][g][t] << " )"
//				<< " ( " << TrM_inv_est_err_err_r[d][g][t] << ", "
//				<< TrM_inv_est_err_err_i[d][g][t] << " )" << std::endl;
//#endif
//			QDPIO::cout << std::endl;
//		}  //for (int t=0; t<NumTs; ++t)
//	} // for disp
	
	//-----------------------------
	// Save results
	//-----------------------------
	char buffer[250];
	char buffer_LaMET[250];
	char buffer_cr[250];
	if (chkout_order == -1){  // -1 means that this checkout is the final
		sprintf(buffer, "%s_fn", out_fname.c_str());
		sprintf(buffer_cr, "%s_cr_fn", out_fname.c_str());
		sprintf(buffer_LaMET, "%s_fn", lamet_out_fname.c_str());
	} else{
		sprintf(buffer, "%s_%02d", out_fname.c_str(), chkout_order);
		sprintf(buffer_cr, "%s_cr_%02d", out_fname.c_str(), chkout_order);
		sprintf(buffer_LaMET, "%s_%02d", lamet_out_fname.c_str(), chkout_order);
	}
	
	std::string out_fname_c(buffer);
	std::string out_fname_cr_c(buffer_cr);
	std::string lamet_out_fname_c(buffer_LaMET);
	
	//write moment into file
	TextFileWriter fout(out_fname_c);
	
	for (int d = 0; d < num_disp_mom; ++d) {
		
		for (int t = 0; t < NumTs; ++t) {
			
#ifdef CALC_ERR_ERR
			fout << "# d t   g  Tr[M^-1 g_i]_re  Tr[M^-1 g_i]_im   StatErr_re      StatErr_im       StatErrErr_re      StatErrErr_im" << "\n";
			
			for (int g=0; g<NUM_G; ++g) {
				char buffer[250];
				sprintf(buffer, "%d %3d %2d %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n",
						link_patterns[d],
						timeslices[t],
						g,
						TrM_inv_av[d][g][t].elem().elem().elem().real(),
						TrM_inv_av[d][g][t].elem().elem().elem().imag(),
						TrM_inv_err_r[d][g][t],
						TrM_inv_err_i[d][g][t],
						TrM_inv_est_err_err_r[d][g][t],
						TrM_inv_est_err_err_i[d][g][t]);
				std::string oline(buffer);
				fout << oline;
			}  // for (int g=0; g<NUM_G; ++g)
#else
			fout << "#d  t   g  Tr[M^-1 g_i]_re  Tr[M^-1 g_i]_im   StatErr_re      StatErr_im"
			<< "\n";
			
			for (int g = 0; g < NUM_G; ++g) {
				char buffer[250];
				sprintf(buffer, "%d %3d %2d %16.8e %16.8e %16.8e %16.8e\n",
						link_patterns[d] , timeslices[t], g,
						TrM_inv_av[d][g][t].elem().elem().elem().real(),
						TrM_inv_av[d][g][t].elem().elem().elem().imag(),
						TrM_inv_err_r[d][g][t], TrM_inv_err_i[d][g][t]);
				std::string oline(buffer);
				fout << oline;
			}  // for (int g=0; g<NUM_G; ++g)
#endif
			
		}  // for (int t=0; t<NumTs; ++t)
	}  // for disp
	
	fout.close();
	
	// write the LP and Correction into separate files.
		TextFileWriter fout_cr(out_fname_cr_c);
		
		for (int d = 0; d < num_disp_mom; ++d) {
			
			for (int t = 0; t < NumTs; ++t) {
				
	
				fout_cr << "# d t   g  Tr[M^-1 g_i]_LP_re  Tr[M^-1 g_i]_LP_im   StatErr_re      StatErr_im       Correction_re      Correction_im	Correction_Err_re	Correction_Err_im" << "\n";
				
				for (int g=0; g<NUM_G; ++g) {
					char buffer[250];
					sprintf(buffer, "%d %3d %2d %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n",
							link_patterns[d],
							timeslices[t],
							g,
							TrM_inv_av_LP[d][g][t].elem().elem().elem().real(),
							TrM_inv_av_LP[d][g][t].elem().elem().elem().imag(),
							TrM_inv_LP_err_r[d][g][t],
							TrM_inv_LP_err_i[d][g][t],
							TrM_inv_av_C[d][g][t].elem().elem().elem().real(),
							TrM_inv_av_C[d][g][t].elem().elem().elem().imag(),
							TrM_inv_C_err_r[d][g][t],
							TrM_inv_C_err_i[d][g][t]);
					std::string oline(buffer);
					fout_cr << oline;
				}  // for (int g=0; g<NUM_G; ++g)
	
				
			}  // for (int t=0; t<NumTs; ++t)
		}  // for disp
		
		fout_cr.close();
		
		
	// write LaMET type operators into a separate file
	TextFileWriter fout_lamet(lamet_out_fname_c);
		
		for (int d = 0; d < num_disp; ++d) {
			if(d>0 && d<num_disp_mom) continue;
			for (int t = 0; t < NumTs; ++t) {
				
	#ifdef CALC_ERR_ERR
				fout << "# d t   g  Tr[M^-1 g_i]_re  Tr[M^-1 g_i]_im   StatErr_re      StatErr_im       StatErrErr_re      StatErrErr_im" << "\n";
				
				for (int g=0; g<NUM_G; ++g) {
					char buffer[250];
					sprintf(buffer, "%d %3d %2d %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n",
							link_patterns[d],
							timeslices[t],
							g,
							TrM_inv_av[d][g][t].elem().elem().elem().real(),
							TrM_inv_av[d][g][t].elem().elem().elem().imag(),
							TrM_inv_err_r[d][g][t],
							TrM_inv_err_i[d][g][t],
							TrM_inv_est_err_err_r[d][g][t],
							TrM_inv_est_err_err_i[d][g][t]);
					std::string oline(buffer);
					fout_lamet << oline;
				}  // for (int g=0; g<NUM_G; ++g)
	#else
				fout_lamet << "#d  t   g  Tr[M^-1 g_i]_re  Tr[M^-1 g_i]_im   StatErr_re      StatErr_im"
				<< "\n";
				
				for (int g = 0; g < NUM_G; ++g) {
					char buffer[250];
					sprintf(buffer, "%d %3d %2d %16.8e %16.8e %16.8e %16.8e\n",
							link_patterns[d] , timeslices[t], g,
							TrM_inv_av[d][g][t].elem().elem().elem().real(),
							TrM_inv_av[d][g][t].elem().elem().elem().imag(),
							TrM_inv_err_r[d][g][t], TrM_inv_err_i[d][g][t]);
					std::string oline(buffer);
					fout_lamet << oline;
				}  // for (int g=0; g<NUM_G; ++g)
	#endif
				
			}  // for (int t=0; t<NumTs; ++t)
		}  // for disp
		
		fout_lamet.close();
	
	//-----------------------------
	// Save number of iterations
	//-----------------------------
	if (Layout::primaryNode()) {
		// only the master node do the job
		char buffer[250];
		sprintf(buffer, "%s_iters", out_fname.c_str());
		std::string out_fname_c(buffer);
		
		std::ofstream fout;
		fout.open(out_fname_c, std::ios::out | std::ios::app);
		
		if (chkout_order == -1)  // -1 means that this checkout is the final
			fout << " fn";
		else
			fout << " " << std::setfill('0') << std::setw(2) << chkout_order;
		
		fout << "  " << Nr_LP;
		
		if (Restarted)
			fout << "  R" << std::endl;
		else
			fout << std::endl;
	}
}  // double checkout()

//===============================================
// Hopping parameter expansion (HPE)
//===============================================
// Let us write the clover operator M into
//   M = (1/2kappa) (1-kappa D)
// where D includes both the hopping term and the clover term.
// For clover action, x=y case, the trace of the M^-1 \Gamma is
//   Tr(M^-1 \Gamma) = Tr[ (2kappa*I + kappa^2 D^2 M^-1) \Gamma ]
// Except for the scalar case, the first term in Tr[] canceled
// so what we need to calculate is kappa^2 D^2 M^-1
// where D =(1/kappa -2M); kappa D = (1 - 2*kappa*M)
//
void do_HPE(Real kappa, LatticeFermion &psi,
			Handle<LinearOperator<LatticeFermion>> &M) {
	Real mtwo_kappa = -2.0 * kappa;
	
	LatticeFermion psi0, psi1;
	
	psi0 = zero;
	psi1 = psi;
	
	(*M)(psi0, psi, PLUS);  // psi = M psi
	psi0 *= mtwo_kappa;     // psi = -2*kappa*M psi
	psi = psi1 + psi0;      // psi = (1 - 2*kappa*M)psi
	
	psi0 = zero;
	psi1 = psi;
	
	(*M)(psi0, psi, PLUS);  // psi = M psi
	psi0 *= mtwo_kappa;     // psi = -2*kappa*M psi
	psi = psi1 + psi0;      // psi = (1 - 2*kappa*M)psi
}
//====================================================================
// Shift Links: For first moment, do normal shift. For second moment,
// do shift based on the first moment operators. For LaMET moment, do
// separate shifts on each direction.
//====================================================================

void shift_link(int d, LatticeFermion &chi,LatticeFermion &psi, LatticeFermion &shift_psi, int num_disp_mom, int num_disp, int disp, multi1d<LatticeColorMatrix> U){
	if(d<num_disp_mom){
		if(disp==0){
			chi = psi;
			shift_psi=chi;
		} else if(disp>9 && disp<18){
			chi = psi;
			int mu=disp%10;
			if (mu < 4) {
				shift_psi = U[mu] * shift(chi, FORWARD, mu);
			} else {
				shift_psi = shift(adj(U[mu-4])*chi, BACKWARD, mu-4);
			}
			chi=shift_psi;
		} else if(disp>199 && disp <278){
			int mu2=disp%10;
			if (mu2 < 4) {
				shift_psi = U[mu2] * shift(chi, FORWARD, mu2);
			} else {
				shift_psi = shift(adj(U[mu2-4])*chi, BACKWARD, mu2-4);
			}
		}
	} else{
		int mu=disp/100;
		int link_length=disp%100;
		if(link_length ==1) chi=psi;
		else chi=shift_psi;
		if (mu < 4) {
			shift_psi = U[mu] * shift(chi, FORWARD, mu);
		} else {
			shift_psi = shift(adj(U[mu-4])*chi, BACKWARD, mu-4);
		}
	}
}

//====================================================================
// Main program
//====================================================================
int main(int argc, char **argv) {
	// Put the machine into a known state
	Chroma::initialize(&argc, &argv);
	
	// Put this in to enable profiling etc.
	START_CODE();
	
	QDPIO::cout << "Linkage = " << linkageHack() << std::endl;
	
	StopWatch snoop;
	snoop.reset();
	snoop.start();
	
	// Instantiate xml reader for DATA
	// if you used -i input file, Chroma::getXMLInputFileName() returns
	//  the filename you specified
	// otherwise it will return "DATA"
	XMLReader xml_in;
	Inline_input_t input;
	
	try {
		xml_in.open(Chroma::getXMLInputFileName());
		read(xml_in, "/disco", input);
	} catch (const std::string& e) {
		QDPIO::cerr << "DISCO: Caught Exception reading XML: " << e << std::endl;
		QDP_abort(1);
	} catch (std::exception& e) {
		QDPIO::cerr << "DISCO: Caught standard library exception: " << e.what()
		<< std::endl;
		QDP_abort(1);
	} catch (...) {
		QDPIO::cerr << "DISCO: caught generic exception reading XML" << std::endl;
		QDP_abort(1);
	}
	
	XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
	push(xml_out, "disco");
	
	// Write out the input
	write(xml_out, "Input", xml_in);
	
	Layout::setLattSize(input.param.nrow);
	Layout::create();
	
	proginfo(xml_out);    // Print out basic program info
	
	// Initialize the RNG
	QDP::RNG::setrn(input.rng_seed);
	write(xml_out, "RNG", input.rng_seed);
	
	// Initialize stop watch
	StopWatch swatch;
	
	// Start up the config
	swatch.reset();
	multi1d<LatticeColorMatrix> U(Nd);
	XMLReader gauge_file_xml, gauge_xml;
	
	// Start up the gauge field
	QDPIO::cout << "Attempt to read gauge field" << std::endl;
	swatch.start();
	try {
		std::istringstream xml_c(input.cfg.xml);
		XMLReader cfgtop(xml_c);
		QDPIO::cout << "Gauge initialization: cfg_type = " << input.cfg.id
		<< std::endl;
		
		Handle<GaugeInit> gaugeInit(
									TheGaugeInitFactory::Instance().createObject(input.cfg.id, cfgtop,
																				 input.cfg.path));
		(*gaugeInit)(gauge_file_xml, gauge_xml, U);
	} catch (std::bad_cast) {
		QDPIO::cerr << "DISCO: caught cast error" << std::endl;
		QDP_abort(1);
	} catch (std::bad_alloc) {
		// This might happen on any node, so report it
		QDPIO::cerr << "DISCO: caught bad memory allocation" << std::endl;
		QDP_abort(1);
	} catch (const std::string& e) {
		QDPIO::cerr << "DISCO: Caught Exception: " << e << std::endl;
		QDP_abort(1);
	} catch (std::exception& e) {
		QDPIO::cerr << "DISCO: Caught standard library exception: " << e.what()
		<< std::endl;
		QDP_abort(1);
	} catch (...) {
		// This might happen on any node, so report it
		QDPIO::cerr << "DISCO: caught generic exception during gaugeInit"
		<< std::endl;
		QDP_abort(1);
	}
	swatch.stop();
	
	QDPIO::cout << "Gauge field successfully read: time= "
	<< swatch.getTimeInSeconds() << " secs" << std::endl;
	
	XMLBufferWriter config_xml;
	config_xml << gauge_xml;
	
	// Write out the config header
	write(xml_out, "Config_info", gauge_xml);
	
	//====================================================================
	// Prepare variables
	//====================================================================
	Real kappa = massToKappa(input.param.inverter.mass);
	
	QDPIO::cout << "kappa = " << kappa << std::endl;
	
	bool use_HPE = input.param.use_HPE;
	if (use_HPE) QDPIO::cout << "Use Hopping Parameter Expansion" << std::endl;
	
	bool setupParamLP_1 = input.param.inverter.setupParamLP_1;
	
	std::string NoiseSrcType = input.param.noise_src.NoiseSrcType;
	QDPIO::cout << "NoiseSrcType = " << NoiseSrcType << std::endl;
	
	// if noise_input_version == 1, do not use restart
	// if noise_input_version == 2, restart the loop if std of correction is large
	int noise_input_version = input.param.noise_src.version;
	int Nr_LPHP_ratio = input.param.noise_src.Nr_LPHP_ratio;
	int Max_Nr_LP = input.param.noise_src.Max_Nr_LP;
	int Min_Nr_LP = input.param.noise_src.Min_Nr_LP;
	int restart_NrLP = input.param.noise_src.restart_NrLP;
	double restart_factor = input.param.noise_src.restart_factor;
	
	bool dilute = input.param.noise_src.dilute;
	int NumTs = input.param.noise_src.timeslices.size();
	std::vector<int> timeslices;
	for (int i = 0; i < NumTs; ++i)
		timeslices.push_back(input.param.noise_src.timeslices[i]);
	
	
	if (dilute) {
		QDPIO::cout << "Calculate disconnected contribution on timeslice = ";
		for (int i = 0; i < timeslices.size(); ++i)
			QDPIO::cout << timeslices[i] << " ";
		QDPIO::cout << std::endl;
	}
	
	int num_checkp = input.param.checkpoint.size();
	std::vector<Checkpoint_t> checkp;
	for (int i = 0; i < num_checkp; ++i) {
		checkp.push_back(input.param.checkpoint[i]);
		
		// Index of the order of the checkout; initialization
		checkp.back().chkout_order = 1;
		
		//sort checkpoints
		for (int j = 0; j < checkp[i].checkpoints.size(); ++j)
			checkp[i].cp.push_back(checkp[i].checkpoints[j]);
		
		if (checkp[i].version == 1)
			// sort checkpoints in descending order
			std::sort(checkp[i].cp.begin(), checkp[i].cp.end(),
					  std::greater<double>());
		else
			// sort checkpoints in ascending order
			std::sort(checkp[i].cp.begin(), checkp[i].cp.end());
	}
	
	// Link length and directions
	std::vector<int> link_dirs;
	for (int i = 0; i < input.param.displacement.link_dirs.size(); ++i)
		link_dirs.push_back(input.param.displacement.link_dirs[i]);
	int link_max=input.param.displacement.link_max;
	std::vector<int> link_patterns;
	link_pattern(link_patterns,link_dirs,link_max);
	int num_disp=link_patterns.size();
	int num_disp_mom=pow(8,2)+1;
	
	// Noise source (eta) and Solution (psi)
	LatticeFermion eta, psi;
	
	// Variables for error analysis
	ErrAnlyVars errAnly(num_disp,NumTs);
	// Temp variables
	Complex TrM_inv;
	LatticeFermion chi, shift_psi;
#ifdef CALC_ERR_ERR
	multi3d<DComplex> TrM_inv_est_LP_sum(num_disp, NUM_G, NumTs);
	for (int d=0; d < num_disp; ++d)
		for (int g=0; g < NUM_G; ++g)
			for (int t=0; t < NumTs; ++t)
				TrM_inv_est_LP_sum[d][g][t] = zero;
#endif
	
	LatticeComplex corr_fn;
	multi1d<DComplex> corr_fn_t(Layout::lattSize()[Nd - 1]);
	
	LatticeBoolean mask = false;
	
	// Machinery to do timeslice sums with
	Set TS;
	TS.make(TSF::TimeSliceFunc(Nd - 1));
	
	//====================================================================
	// Do Inversion
	//====================================================================
	try {
		typedef LatticeFermion T;
		typedef multi1d<LatticeColorMatrix> P;
		
		//
		// Initialize fermion action
		//
		std::istringstream xml_sLP(input.param.inverter.fermactLP.xml);
		XMLReader fermacttopLP(xml_sLP);
		QDPIO::cout << "FermAct_LP = " << input.param.inverter.fermactLP.id
		<< std::endl;
		
		std::istringstream xml_sHP(input.param.inverter.fermactHP.xml);
		XMLReader fermacttopHP(xml_sHP);
		QDPIO::cout << "FermAct_HP = " << input.param.inverter.fermactHP.id
		<< std::endl;
		
		std::istringstream xml_sHPE(input.param.inverter.fermactHPE.xml);
		XMLReader fermacttopHPE(xml_sHPE);
		QDPIO::cout << "FermAct_HPE = " << input.param.inverter.fermactHPE.id
		<< std::endl;
		
		// Generic Wilson-type fermion action handles
		Handle<WilsonTypeFermAct<T, P, P>> S_f_LP(
												  TheWilsonTypeFermActFactory::Instance().createObject(
																									   input.param.inverter.fermactLP.id, fermacttopLP,
																									   input.param.inverter.fermactLP.path));
		
		Handle<WilsonTypeFermAct<T, P, P>> S_f_HP(
												  TheWilsonTypeFermActFactory::Instance().createObject(
																									   input.param.inverter.fermactHP.id, fermacttopHP,
																									   input.param.inverter.fermactHP.path));
		
		Handle<WilsonTypeFermAct<T, P, P>> S_f_HPE(
												   TheWilsonTypeFermActFactory::Instance().createObject(
																										input.param.inverter.fermactHPE.id, fermacttopHPE,
																										input.param.inverter.fermactHPE.path));
		
		Handle<FermState<T, P, P>> stateLP(S_f_LP->createState(U));
		Handle<FermState<T, P, P>> stateHP(S_f_HP->createState(U));
		Handle<FermState<T, P, P>> stateHPE(S_f_HPE->createState(U));
		
		// Solvers
		Handle<SystemSolver<LatticeFermion>> PP_LP;
		Handle<SystemSolver<LatticeFermion>> PP_LP_1;
		Handle<SystemSolver<LatticeFermion>> PP_HP;
		
		PP_LP = S_f_LP->qprop(stateLP, input.param.inverter.invParamLP);
		PP_LP_1 = S_f_LP->qprop(stateLP, input.param.inverter.invParamLP_1);
		// PP_HP will be initialized in the loop
		
		// For the multigrid inverter, the only way to re-initialize the inverter
		// parameters is to call S->qprop(). Hence, in order to switch the
		// parameters  LP <--> HP, we need to call qprop() again.
		// The flag "HP_inv_called" tells us the HP inverter is called
		// and we need to re-initialize the LP parameters.
		bool HP_inv_called = false;
		
		// Boolean that indicates whether the loop is restarted
		bool Restarted = false;
		
		// Linear Operator for Hopping Parameter Expansion
		Handle<LinearOperator<T>> M(S_f_HPE->linOp(stateHPE));
		
		//====================================================================
		// Truncated Solver
		//====================================================================
		// Here we use "Truncated Solver Method (TSM)" to estimate
		// inverse all-to-all propagator. It averages over NrLP
		// times of low-precision inversions, and corrects it by
		// using the NrHP times of high-precision inversions.
		
		//--------------------------------------------------------------------
		// Calculate estimation of propagator
		// M^-1 = (1 / Nr_LP) sum[ |psi_i>_LP <eta_i| ]
		// by using low precision calculation
		//--------------------------------------------------------------------
		
		int count_hp = 0;
		for (int count_lp = 1; count_lp <= Max_Nr_LP; ++count_lp) {
			swatch.reset();
			swatch.start();
			
			QDPIO::cout << "TSM Low Precision Estimation loop; iter = " << count_lp
			<< std::endl;
			
			// Make noise source
			if (NoiseSrcType == "Z2") {
				z2_src(eta);
				
				// gaussian() and z2_src require normalization factor (1/sqrt(2))
				eta *= INV_SQRT2;
			} else if(NoiseSrcType ==  "POINT") {
				LatticeFermion eta00;
				zN_src(eta00,1);
				// Make a source only at the origin (0,0,0,0)
				eta=where((Layout::latticeCoordinate(0)+ Layout::latticeCoordinate(1) + Layout::latticeCoordinate(2)+Layout::latticeCoordinate(3)) == 0,eta00,LatticeFermion(zero));
			} else if (NoiseSrcType == "GAUSSIAN") {
				gaussian(eta);
				
				// gaussian() and z2_src require normalization factor (1/sqrt(2))
				eta *= INV_SQRT2;
			} else if (NoiseSrcType == "Z2REAL") {
				zN_src(eta, 2);
			} else if (NoiseSrcType == "Z4") {
				zN_src(eta, 4);
			} else {
				QDPIO::cerr << "Error! Unknown noise source type " << NoiseSrcType
				<< std::endl;
				QDPIO::cerr << "Allowed types are Z2, Z2REAL, Z4 and GAUSSIAN."
				<< std::endl;
				QDP_abort(1);
			}
			
			// Calculate only on the timeslices
			if (dilute) {
				for (int t = 0; t < NumTs; ++t)
					mask |= Layout::latticeCoordinate(3) == timeslices[t];
				eta = where(mask, eta, LatticeFermion(zero));
			}
			
			psi = zero;  // Initialize psi
			
			// Since qprop trashes the source, it should be copied safely
			chi = eta;
			
			// Low precision solver
			// This is called every iteration because this is the only way to
			// (re)initialize the inverter parameters of multigrid solver
			if (HP_inv_called) {
				// if allocating new pointer to a Handle, it deletes old pointer
				if (Restarted)
					PP_LP = S_f_LP->qprop(stateLP, input.param.noise_src.invParamLP_r);
				else
					PP_LP = S_f_LP->qprop(stateLP, input.param.inverter.invParamLP);
				HP_inv_called = false;
			}
			
			// Calculate psi by using Dirac Inverter with Low Precision
			SystemSolverResults_t res;
			if (setupParamLP_1 && count_lp == 1 && !Restarted) {
				// When the inverter runs at the first time, it may use different input
				// parameters, especially in the multigrid solver
				res = (*PP_LP_1)(psi, chi);
			} else {
				res = (*PP_LP)(psi, chi);
			}
			
			// Hopping parameter expansion (HPE)
			if (use_HPE) do_HPE(kappa, psi, M);
			
			// Calculate Tr(M^-1) = (1/N) sum_i <eta_i| \gamma |psi_i>
			// Displace the source (covariantly) to form a new operator
			
			for (int d=0; d<num_disp; ++d){
				int disp=link_patterns[d];
				//if(Layout::primaryNode()) std::cout << "calculating link "<< d << " in direction "<< disp <<std::endl;
				
				shift_link(d,chi,psi,shift_psi,num_disp_mom,num_disp,disp,U);
				

//				if(d<num_disp_mom){
//					chi = psi;
//					if(disp==0){
//						shift_psi=chi;
//					}
//					//Calculate one link of pattern 1_mu
//					else if (disp>9 && disp<18){
//						int mu=disp%10;
//						if (mu < 4) {
//							shift_psi = U[mu] * shift(chi, FORWARD, mu);
//						} else {
//							shift_psi = shift(adj(U[mu-4])*chi, BACKWARD, mu-4);
//						}
//					}
//					//Calculate two link of pattern 2_mu_nu
//					else if(disp>199 && disp <278){
//						int mu1=(disp/10)%10;
//						int mu2=disp%10;
//						if (mu1 < 4) {
//							shift_psi = U[mu1] * shift(chi, FORWARD, mu1);
//						} else {
//							shift_psi = shift(adj(U[mu1-4])*chi, BACKWARD, mu1-4);
//						}
//						chi=shift_psi;
//						if (mu2 < 4) {
//							shift_psi = U[mu2] * shift(chi, FORWARD, mu2);
//						} else {
//							shift_psi = shift(adj(U[mu2-4])*chi, BACKWARD, mu2-4);
//						}
//					}
//				}else if(link_max>2){
//					int mu=disp/100;
//					int link_length=disp%100;
//
//					if(link_length ==1) chi=psi;
//					else chi=shift_psi;
//					if (mu < 4) {
//						shift_psi = U[mu] * shift(chi, FORWARD, mu);
//					} else {
//						shift_psi = shift(adj(U[mu-4])*chi, BACKWARD, mu-4);
//					}
//
//				}
				for (int g = 0; g < NUM_G; ++g) {
					corr_fn = localInnerProduct(eta, gamma_ops(g) * shift_psi);
					corr_fn_t = sumMulti(corr_fn, TS);
					
					for (int t = 0; t < NumTs; ++t) {
						TrM_inv = corr_fn_t[timeslices[t]];
						
						// For scalar case, HPE should correct Tr(2 kappa I) = 24*kappa*L^3
						if (g == 0 && use_HPE)
							TrM_inv += 24.0 * kappa * pow(Layout::lattSize()[0], 3);
						
						errAnly.TrM_inv_sum_LP[d][g][t] += TrM_inv;
						
						// Statistical error estimation
						errAnly.TrM_inv_LP_sqsum_r[d][g][t] += pow(
																									   (double)TrM_inv.elem().elem().elem().real(), 2);
						errAnly.TrM_inv_LP_sqsum_i[d][g][t] += pow(
																									   (double)TrM_inv.elem().elem().elem().imag(), 2);
						
#ifdef CALC_ERR_ERR
						TrM_inv_est_LP_sum[d][g][t] += TrM_inv;
#endif
					}//multi3d<DComplex> TrM_inv_C_HP(num_disp, NUM_G, NumTs);
				}//for (int g = 0; g < NUM_G; ++g)
			}//for (int d=0; d<num_disp; ++d)
			
			
			
			// Print time used
			swatch.stop();
			
			QDPIO::cout << "TSM LP loop: time= " << swatch.getTimeInSeconds()
			<< " secs" << std::endl;
			QDPIO::cout << std::endl;
			
			//--------------------------------------------------------------------
			// Calculate correction to the above low-precision estimation
			// Correction = (1 / Nr_HP) sum [ (|psi_j>_HP - |psi_j>_LP>) <eta_j| ]
			// by using LP and HP results
			//
			// This part runs every Nr_LPHP_ratio times of LP iterations
			//--------------------------------------------------------------------
			if (Nr_LPHP_ratio > 0 && count_lp % Nr_LPHP_ratio == 0) {
				swatch.reset();
				swatch.start();
				
				QDPIO::cout << "TSM Correction Estimation Loop; iter = " << count_hp + 1
				<< std::endl;
				
				// Temp variables
				multi3d<DComplex> TrM_inv_C_HP(num_disp, NUM_G, NumTs);
				multi3d<DComplex> TrM_inv_C_LP(num_disp, NUM_G, NumTs);
				for (int d = 0; d < num_disp; ++d)
					for (int g = 0; g < NUM_G; ++g)
						for (int t = 0; t < NumTs; ++t) {
							TrM_inv_C_HP[d][g][t] = zero;
							TrM_inv_C_LP[d][g][t] = zero;
						}
				
				// Make noise source
				if (NoiseSrcType == "Z2") {
					z2_src(eta);
					
					// gaussian() and z2_src require normalization factor (1/sqrt(2))
					eta *= INV_SQRT2;
				} else if(NoiseSrcType ==  "POINT") {
					LatticeFermion eta00;
					zN_src(eta00,1);
					// Make a source only at the origin (0,0,0,0)
					eta=where((Layout::latticeCoordinate(0)+ Layout::latticeCoordinate(1) + Layout::latticeCoordinate(2)+Layout::latticeCoordinate(3)) == 0,eta00,LatticeFermion(zero));
				} else if (NoiseSrcType == "GAUSSIAN") {
					gaussian(eta);
					
					// gaussian() and z2_src require normalization factor (1/sqrt(2))
					eta *= INV_SQRT2;
				} else if (NoiseSrcType == "Z2REAL") {
					zN_src(eta, 2);
				} else if (NoiseSrcType == "Z4") {
					zN_src(eta, 4);
				} else {
					QDPIO::cerr << "Error! Unknown noise source type " << NoiseSrcType
					<< std::endl;
					QDPIO::cerr << "Allowed types are Z2, Z2REAL, Z4 and GAUSSIAN."
					<< std::endl;
					QDP_abort(1);
				}
				
				// Calculate only on the timeslices
				if (dilute) {
					for (int t = 0; t < NumTs; ++t)
						mask |= Layout::latticeCoordinate(3) == timeslices[t];
					eta = where(mask, eta, LatticeFermion(zero));
				}
				
				//-----------------------------------------------------
				// Low precision calculation with the same source
				//-----------------------------------------------------
				psi = zero;  // Initialize psi
				
				// Since qprop trashes the source, it should be copied safely
				chi = eta;
				
				// Calculate psi by using Dirac Inverter with Low Precision
				SystemSolverResults_t res2 = (*PP_LP)(psi, chi);
				
				// Hopping parameter expansion (HPE)
				if (use_HPE) do_HPE(kappa, psi, M);
				
				// Calculate Tr(M^-1) = (1/N) sum_i <eta_i| \gamma |psi_i>
				// Displace the source (covariantly) to form a new operator
				
				for (int d=0; d<num_disp; ++d){
					int disp=link_patterns[d];
					//if(Layout::primaryNode()) std::cout << "calculating link "<< d << " in direction "<< disp <<std::endl;
					shift_link(d,chi,psi,shift_psi,num_disp_mom,num_disp,disp,U);
					
//					if(d<num_disp_mom){
//						chi = psi;
//						if(disp==0) shift_psi=chi;
//						// link one with patter 1_mu
//						else if(disp>9 && disp < 18){
//							int mu=disp%10;
//							if (mu<4) {
//								shift_psi = U[mu] * shift(chi, FORWARD, mu);
//							} else {
//								shift_psi = shift(adj(U[mu-4])*chi, BACKWARD, mu-4);
//							}
//						}
//						// link two with patter 2_mu_nu
//						else if(disp>199 && disp < 278){
//							int mu1=(disp/10)%10;
//							int mu2=disp%10;
//							if (mu1<4) {
//								shift_psi = U[mu1] * shift(chi, FORWARD, mu1);
//							} else {
//								shift_psi = shift(adj(U[mu1-4])*chi, BACKWARD, mu1-4);
//							}
//							chi=shift_psi;
//							if (mu2<4) {
//								shift_psi = U[mu2] * shift(chi, FORWARD, mu2);
//							} else {
//								shift_psi = shift(adj(U[mu2-4])*chi, BACKWARD, mu2-4);
//							}
//						}
//					}else if(link_max>2){
//						int mu=disp/100;
//						int link_length=disp%100;
//						if(link_length ==1) chi=psi;
//						else chi=shift_psi;
//						if (mu < 4) {
//							shift_psi = U[mu] * shift(chi, FORWARD, mu);
//						} else {
//							shift_psi = shift(adj(U[mu-4])*chi, BACKWARD, mu-4);
//						}
//					}
					
					for (int g = 0; g < NUM_G; ++g) {
						corr_fn = localInnerProduct(eta, gamma_ops(g) * shift_psi);
						corr_fn_t = sumMulti(corr_fn, TS);
						
						// For the correction term, we don't need to add constant part Tr[2 kappa I]
						for (int t = 0; t < NumTs; ++t)
							TrM_inv_C_LP[d][g][t] = corr_fn_t[timeslices[t]];
					}
				
				}//loop for (int d=0; d<num_disp; ++d)
				
				//-----------------------------------------------------
				// High precision calculation
				//-----------------------------------------------------
				// Result (psi) of Low precision inversion is used as a
				// initial guess for this High precision calculation
				
				// Since qprop trashes the source, it should be copied safely
				chi = eta;
				
				// High precision solver
				// This is called every iteration because this is the only way to
				// (re)initialize the inverter parameters of multigrid solver
				// if allocating new pointer to a Handle, it deletes old pointer
				PP_HP = S_f_HP->qprop(stateHP, input.param.inverter.invParamHP);
				HP_inv_called = true;
				
				// Calculate psi by using Dirac Inverter with High Precision
				res = (*PP_HP)(psi, chi);
				
				// Hopping parameter expansion (HPE)
				if (use_HPE) do_HPE(kappa, psi, M);
				
				// Calculate Tr(M^-1) = (1/N) sum_i <eta_i| \gamma |psi_i>
				// Displace the source (covariantly) to form a new operator
				
				for (int d=0; d<num_disp; ++d){
					int disp=link_patterns[d];
					
					//if(Layout::primaryNode()) std::cout << "calculating link "<< d << " in direction "<< link_patterns[d] <<std::endl;
					
					shift_link(d,chi,psi,shift_psi,num_disp_mom,num_disp,disp,U);
					
//					if(d<num_disp_mom){
//						chi = psi;
//						if (disp==0) shift_psi=chi;
//						// link one with patter 1_mu
//						else if(disp>9 && disp < 18){
//							int mu=disp%10;
//							if (mu<4) {
//								shift_psi = U[mu] * shift(chi, FORWARD, mu);
//							} else {
//								shift_psi = shift(adj(U[mu-4])*chi, BACKWARD, mu-4);
//							}
//						}
//						// link two with patter 2_mu_nu
//						else if(disp>199 && disp < 278){
//							int mu1=(disp/10)%10;
//							int mu2=disp%10;
//							if (mu1<4) {
//								shift_psi = U[mu1] * shift(chi, FORWARD, mu1);
//							} else {
//								shift_psi = shift(adj(U[mu1-4])*chi, BACKWARD, mu1-4);
//							}
//							chi=shift_psi;
//							if (mu2<4) {
//								shift_psi = U[mu2] * shift(chi, FORWARD, mu2);
//							} else {
//								shift_psi = shift(adj(U[mu2-4])*chi, BACKWARD, mu2-4);
//							}
//						}
//					}else if(link_max>2){
//						int mu=disp/100;
//						int link_length=disp%100;
//						if(link_length ==1) chi=psi;
//						else chi=shift_psi;
//						if (mu < 4) {
//							shift_psi = U[mu] * shift(chi, FORWARD, mu);
//						} else {
//							shift_psi = shift(adj(U[mu-4])*chi, BACKWARD, mu-4);
//						}
//					}
					
					for (int g = 0; g < NUM_G; ++g) {
						corr_fn = localInnerProduct(eta, gamma_ops(g) * shift_psi);
						corr_fn_t = sumMulti(corr_fn, TS);
						
						// For the correction term, we don't need to add constant part Tr[2 kappa I]
						for (int t = 0; t < NumTs; ++t)
							TrM_inv_C_HP[d][g][t] = corr_fn_t[timeslices[t]];
					}
				
				}
				
				//-----------------------------------------------------
				// Calculate correction term
				//-----------------------------------------------------
				for (int d = 0; d < num_disp; ++d)
					for (int g = 0; g < NUM_G; ++g)
						for (int t = 0; t < NumTs; ++t) {
							TrM_inv = TrM_inv_C_HP[d][g][t] - TrM_inv_C_LP[d][g][t];
							errAnly.TrM_inv_sum_C[d][g][t] += TrM_inv;
							
							// Statistical error estimation
							errAnly.TrM_inv_C_sqsum_r[d][g][t] += pow(
																	  (double)TrM_inv.elem().elem().elem().real(), 2);
							errAnly.TrM_inv_C_sqsum_i[d][g][t] += pow(
																	  (double)TrM_inv.elem().elem().elem().imag(), 2);
							
#ifdef CALC_ERR_ERR
							DComplex TrM_inv_est_tmp = TrM_inv_est_LP_sum[d][g][t]/Nr_LPHP_ratio + TrM_inv;
							errAnly.TrM_inv_est[d][g][t].push_back(TrM_inv_est_tmp);
							
							// Reinitialize
							TrM_inv_est_LP_sum[d][g][t] = zero;
#endif
						}
				
				// Print time used
				swatch.stop();
				
				QDPIO::cout << "TSM Correction loop (LP + HP): time= "
				<< swatch.getTimeInSeconds() << " secs" << std::endl;
				QDPIO::cout << std::endl;
				
				// Increase number of hp iteration counter
				++count_hp;
				
			}  // if ( Nr_LPHP_ratio > 0  &&  i%Nr_LPHP_ratio == 0)
			
			//------------------------------------------
			// Restart
			//------------------------------------------
			if (noise_input_version == 2) {
				// Check whether restart is needed
				if (count_lp == restart_NrLP && !Restarted) {
					double ratio_C_LP_err, s_err_max, s_err_av;
					check_acc(count_lp, count_hp, errAnly, link_dirs, link_max, timeslices, ratio_C_LP_err,
							  s_err_max, s_err_av);
					
					// Restart the loop with new LP inverter parameters if the error of
					// correction term is much larger than the error of LP term
					if (ratio_C_LP_err > restart_factor) {
						QDPIO::cout << "Restart with new LP inverter parameters!"
						<< std::endl;
						QDPIO::cout << "count_lp = " << count_lp << ", Cr_err / LP_err = "
						<< ratio_C_LP_err << ", Err_scalar = " << s_err_max << std::endl;
						
						count_lp = 0;  // for-loop will increase this by 1 in the end of this block
						count_hp = 0;
						HP_inv_called = false;
						Restarted = true;
						errAnly = ErrAnlyVars(num_disp,NumTs);
						PP_LP = S_f_LP->qprop(stateLP, input.param.noise_src.invParamLP_r);
						// Note that checkout order is not reinitialized; it will be proceed from where it is
					}  // if (ratio_C_LP_err > restart_factor)
				}  // if (count_lp == restart_NrLP && !Restarted)
			}  // if (noise_input_version == 2)
			
			//------------------------------------------
			// Checkout or Finish
			//------------------------------------------
			if (count_lp >= Min_Nr_LP) {
				// Calculate error of scalar channel (maximum value among timeslices)
				double ratio_C_LP_err, s_err_max, s_err_av;
				check_acc(count_lp, count_hp, errAnly,link_dirs,link_max, timeslices, ratio_C_LP_err,
						  s_err_max, s_err_av);
				
				QDPIO::cout << "count_lp = " << count_lp << ", Cr_err / LP_err = "
				<< ratio_C_LP_err << ", Err_scalar = " << s_err_max << std::endl;
				
				bool all_empty = true;
				
				// Loop over all types of checkpoints
				for (int idx_cp = 0; idx_cp < num_checkp; ++idx_cp) {
					while (1) {
						// break while(1) if checkpoints is empty
						if (checkp[idx_cp].cp.empty()) {
							all_empty &= true;
							break;
						}
						
						double next_checkpoint = checkp[idx_cp].cp.back();
						bool chkout = false;
						
						// if this checkout is the last checkpoint, set it final (-1)
						if (checkp[idx_cp].cp.size() == 1) checkp[idx_cp].chkout_order = -1;
						
						// Determine whether checkout
						switch (checkp[idx_cp].version) {
								// Version 1: save results if the number of iterations reached to a checkpoint
							case 1:
								if (count_lp >= next_checkpoint) chkout = true;
								break;
								
							case 2:
								if (next_checkpoint > s_err_max) chkout = true;
								break;
								
							case 3:
								if (next_checkpoint > s_err_av) chkout = true;
								break;
						}  // switch(checkp[idx_cp].version)
						
						// If it reaches to the Maximum iterations, checkout
						if (count_lp == Max_Nr_LP) {
							chkout = true;
							checkp[idx_cp].chkout_order = -1;
						}
						
						if (chkout) {
							// checkout
							checkout(count_lp, count_hp, errAnly, checkp[idx_cp].OutFileName,checkp[idx_cp].LaMETOutFileName,link_dirs,link_max,
									 timeslices, checkp[idx_cp].chkout_order, Restarted);
							checkp[idx_cp].chkout_order++;
							
							// remove the passed checkpoint (the largest checkout accuracy)
							checkp[idx_cp].cp.pop_back();
						}
						// break while(1) if there is no checkpoints
						else {
							all_empty &= false;  // if one of &= operates false, it will false
							break;
						}
						
						// break while(1) if checkpoints is empty
						if (checkp[idx_cp].cp.empty()) {
							all_empty &= true;
							break;
						}
					}  // end of while(1)
					
				}  // for (int idx_cp = 0; idx_cp < num_checkp; ++idx_cp)
				
				if (all_empty) break;
			}  // if (count_lp >= Min_Nr_LP)
			
		}  // for (int count_lp=1; count_lp <= Max_Nr_LP; ++count_lp)
		
		xml_out.flush();
		pop(xml_out);
		
	}  // try
	catch (const std::string &e) {
		QDPIO::cerr << "Error in inversion: " << e << std::endl;
		QDP_abort(1);
	}
	
	//-----------------------------
	// End program
	//-----------------------------
	snoop.stop();
	QDPIO::cout << "DISCO: total time = " << snoop.getTimeInSeconds() << " secs"
	<< std::endl;
	
	QDPIO::cout << "DISCO: ran successfully" << std::endl;
	
	END_CODE();
	
	Chroma::finalize();
	exit(0);
}
