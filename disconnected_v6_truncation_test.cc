#include "chroma.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"
#include "meas/sources/z2_src.h"
#include "meas/sources/zN_src.h"

#include "gamma_ops.h"  // definition of gamma_ops() ftn
#include <iomanip>      // std::setw
#include <algorithm>
#include <string>
#include <fstream>

// This version of code changes the noise source generation to site by site, generate a random number but skip the noise source when the site in outside of the sublattice. Use QLA random number generator

// This version of code works for the pdfchroma where the SftMom function can take mom2_list as input.

//constexpr int num_links = 1; // calculate links forward and backward displacement
//constexpr int link_dir = 2; // only z direction
//constexpr int NumDisp = 8*num_links + 1;

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
	
	int t_src;
	
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
	
	// below variables are for the program process, not for read
	std::vector<double> cp;
	int chkout_order;
};

struct Params_t {
	multi1d<int> nrow;
	Inverter_t inverter;
	NoiseSource_t noise_src;
	Displacement_t displacement;
	multi1d<int> mom2_list;
	multi1d<Checkpoint_t> checkpoint;
	bool use_HPE;
};

struct Inline_input_t {
	Params_t param;
	GroupXML_t cfg;
	int rng_seed; // random seed for random seed (update the seed manually)
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
	
	// Read time source location
	read(paramtop, "t_src", p.t_src);
	
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
	
}

void read(XMLReader& xml, const std::string& path, Params_t& p) {
	XMLReader paramtop(xml, path);
	read(paramtop, "nrow", p.nrow);
	
	if (paramtop.count("mom2_list") > 0)
		read(paramtop, "mom2_list", p.mom2_list);
	else {
		p.mom2_list.resize(1);
		p.mom2_list[0]=0;
	}
	
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
void checkout(int Nr_LP, multi4d<Complex> &TrM_inv, std::string out_fname, std::vector<int> &link_dirs, int link_max,int chkout_order, bool Restarted, SftMom* &phases) {
	std::vector<int> link_patterns;
	link_pattern(link_patterns,link_dirs,link_max);
	int NumDisp = link_patterns.size();
	int NumDisp_mom;
	if(link_max>2) NumDisp_mom= pow(8,2)+1;
	else NumDisp_mom= pow(link_dirs.size(),link_max)+1;
	int NumMom;
	NumMom=phases->numMom();
	multi2d<int> mom_list(NumMom,3);
	for(int p=0; p<NumMom; ++p) {
		multi1d<int> pp = phases->numToMom(p);
		mom_list[p][0]=pp[0];
		mom_list[p][1]=pp[1];
		mom_list[p][2]=pp[2];
	}
	
	
	


	//------------------------------------------
	// Print all LP measurements
	//------------------------------------------
	
	QDPIO::cout << std::endl << "Checkout - fn; " << out_fname << std::endl;
	


	//-----------------------------
	// Save results
	//-----------------------------

	for (int p = 0; p < NumMom; ++p){
		char buffer[250];
		if (chkout_order == -1){  // -1 means that this checkout is the final
			sprintf(buffer, "%s_qx%d_qy%d_qz%d_fn", out_fname.c_str(),mom_list[p][0],mom_list[p][1],mom_list[p][2]);
		} else{
			sprintf(buffer, "%s_qx%d_qy%d_qz%d_%02d", out_fname.c_str(),mom_list[p][0],mom_list[p][1],mom_list[p][2], chkout_order);
		}
		
		std::string out_fname_c(buffer);
		
		//write moment into file
		TextFileWriter fout(out_fname_c);
		fout << "#LP_idx d g  Tr[M^-1 g_i]_re  Tr[M^-1 g_i]_im "
		<< "\n";
		for (int lp_idx = 0; lp_idx < Nr_LP; ++lp_idx) {
			for (int d = 0; d < NumDisp_mom; ++d) {
			for (int g = 0; g < NUM_G; ++g) {
				char buffer[250];
				sprintf(buffer, "%d %3d %2d %16.8e %16.8e\n",
						lp_idx,link_patterns[d] , g,
						TrM_inv[lp_idx][d][p][g].elem().elem().elem().real(),
						TrM_inv[lp_idx][d][p][g].elem().elem().elem().imag()
						);
				std::string oline(buffer);
				fout << oline;
			}  // for (int g=0; g<NUM_G; ++g)

			
		}  // for (int d = 0; d < 1; ++d)
		}// for LP
		fout.close();
	}  //for (int p = 0; p < NumMom; ++p)

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

void shift_link(int d, LatticeFermion &chi,LatticeFermion &psi, LatticeFermion &shift_psi, int NumDisp_mom, int NumDisp, int disp, multi1d<LatticeColorMatrix> U){
	if(d<NumDisp_mom){
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
	
	// Initialize the RNG seed
	srand(input.rng_seed);
	QDP::Seed qdp_rng_seed(11);
	for(int idx_seed=0;idx_seed<4;idx_seed++) qdp_rng_seed.elem().elem().elem(idx_seed)=rand();
	QDP::RNG::setrn(qdp_rng_seed); // initialize a qdp seed
	write(xml_out, "RNG", input.rng_seed);
	
	int t_src = input.param.noise_src.t_src;
	// Initialize stop watch
	StopWatch swatch;
	
	// Start up the config
	swatch.reset();
	multi1d<LatticeColorMatrix> U(Nd), U_tmp(Nd);
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
		U_tmp=U;
		QDPIO::cout << "Checking lattice shifts:" << std::endl;
		QDPIO::cout << "Before shift (0,0,0,"<< t_src <<"):" << std::endl;
		multi1d<int> coord(4);
		coord[0]=coord[1]=coord[2]=0; coord[3]=t_src;
		QDPIO::cout << peekSite(U,coord) << std::endl;
		if(t_src<=Layout::lattSize()[Nd-1]/2){
			for(int t_shift=0;t_shift<t_src;t_shift++){
				for(int mu=0;mu<Nd-1;mu++){
				U[mu]=shift(U_tmp[mu], BACKWARD, Nd-1);// shift the whole lattice FWD so that the source is located at t=0.
				U_tmp[mu]=U[mu];
				}
			}
		}
		else{
			for(int t_shift=0;t_shift<=Layout::lattSize()[Nd-1]-t_src;t_shift++){
				for(int mu=0;mu<Nd-1;mu++){
				U[mu]=shift(U_tmp[mu], FORWARD, Nd-1);// shift the whole lattice BWD so that the source is located at t=0.
				U_tmp[mu]=U[mu];
				}
			}
		}
		QDPIO::cout << "Checking lattice shifts:" << std::endl;
		QDPIO::cout << "After shift (0,0,0,0):" << std::endl;
		multi1d<int> coord(4);
		coord[0]=coord[1]=coord[2]=coord[3]=0;
		QDPIO::cout << peekSite(U,coord) << std::endl;
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
	
	QDPIO::cout << "Calculate disconnected contribution on timeslice = 0";
	
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
	int NumDisp=link_patterns.size();
	
	int NumDisp_mom;
	if(link_max>2) NumDisp_mom= pow(8,2)+1;
	else NumDisp_mom= pow(link_dirs.size(),link_max)+1;
	SftMom* phases;
	
	
	
	// some dummy variables used for mom2_list functions
	multi1d<SftMomSrcPos_t> origin_offs(1);
	origin_offs[0].src_pos.resize(Nd);
	origin_offs[0].src_pos=0;
	origin_offs[0].t_min = 0;
	int j_decay=Nd-1;
    origin_offs[0].t_max = Layout::lattSize()[j_decay] - 1;
	multi1d<int> mom_offset;
	mom_offset.resize(Nd-1);  /*!< Origin for the momentum */
    mom_offset = 0;
	//end
	
	
	
	phases=new SftMom(input.param.mom2_list, origin_offs, mom_offset, false, j_decay);
	int NumMom=phases->numMom();
	
	// Noise source (eta) and Solution (psi)
	LatticeFermion eta, psi;
	
	// Temp variables
	multi4d<Complex> TrM_inv(Max_Nr_LP,NumDisp,NumMom,NUM_G);
	LatticeFermion chi, shift_psi;
	
	LatticeComplex corr_fn;
	multi2d<DComplex> corr_fn_t;
	
	
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
		
		
		for (int count_lp = 0; count_lp < Max_Nr_LP; ++count_lp) {
			swatch.reset();
			swatch.start();
			
			QDPIO::cout << "TSM Low Precision Estimation loop; iter = " << count_lp
			<< std::endl;
			
			// Make noise source
			if (NoiseSrcType == "Z4") {
				LatticeReal rnd1, theta;
				Real twopiN = Chroma::twopi / 4;   // twopi defined in chroma/lib/chromabase.h
				LatticeComplex c;
				LatticeColorVector colorvec = zero;
				for(int spin_index= 0; spin_index < Ns; ++spin_index)
				  for(int color_index= 0; color_index < Nc; ++color_index)
				  {
				random(rnd1);
/* Tests for random numbers, print the Z4 random number at different sites for comparison.
				QDPIO::cout << "Test random seed"
									  << std::endl;
				for(int n=0;n<1;n++){
					  QDPIO::cout << "Random number "<< n <<": " << std::endl;
				multi1d<int> coord(4);
				for(int i=0;i<4;i++){
					  coord[0]=coord[1]=coord[2]=coord[3]=i;
					  QDPIO::cout <<"(";
					  for(int mu=0;mu<4;mu++)
							  QDPIO::cout << coord[mu] << ", ";
					  QDPIO::cout <<"): "<< peekSite(rnd1,coord) <<std::endl;
					  coord[3]=0;
					  QDPIO::cout <<"(";
					  for(int mu=0;mu<4;mu++)
							  QDPIO::cout << coord[mu] << ", ";
					  QDPIO::cout <<"): "<< peekSite(rnd1,coord) <<std::endl;
				}
				}
 */
				for (int idx_seed=0;idx_seed<4;idx_seed++) {qdp_rng_seed.elem().elem().elem(idx_seed)=rand();}
				QDP::RNG::setrn(qdp_rng_seed);
				theta = twopiN * floor(4*rnd1);
				c = cmplx(cos(theta),sin(theta));

				colorvec = peekSpin(eta,spin_index);

				pokeSpin(eta,pokeColor(colorvec,c,color_index),spin_index);
				  }
			} else {
				QDPIO::cerr << "Error! Unknown noise source type " << NoiseSrcType
				<< std::endl;
				QDPIO::cerr << "Allowed types is Z4."
				<< std::endl;
				QDP_abort(1);
			}
			
			// Calculate only on the timeslices t=0
			eta = where(Layout::latticeCoordinate(3) == 0, eta, LatticeFermion(zero));
			
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
			
			for (int d=0; d<NumDisp; ++d){
				int disp=link_patterns[d];
				//if(Layout::primaryNode()) std::cout << "calculating link "<< d << " in direction "<< disp <<std::endl;
				
				shift_link(d,chi,psi,shift_psi,NumDisp_mom,NumDisp,disp,U);
				
				for (int g = 0; g < NUM_G; ++g) {
					corr_fn = localInnerProduct(eta, gamma_ops(g) * shift_psi);
					corr_fn_t = phases->sft(corr_fn);
					for (int p=0; p<NumMom; ++p){
						TrM_inv[count_lp][d][p][g] = corr_fn_t[p][0];
						
						// For scalar case, HPE should correct Tr(2 kappa I) = 24*kappa*L^3
						if (g == 0 && use_HPE)
							TrM_inv[count_lp][d][p][g] += 24.0 * kappa * pow(Layout::lattSize()[0], 3);
						
					}//for (int p = 0; p < NumMom; ++p)
				}//for (int g = 0; g < NUM_G; ++g)
			}//for (int d=0; d<NumDisp; ++d)
			
			
			
			// Print time used
			swatch.stop();
			
			QDPIO::cout << "TSM LP loop: time= " << swatch.getTimeInSeconds()
			<< " secs" << std::endl;
			QDPIO::cout << std::endl;
			
			

			//------------------------------------------
			// Checkout or Finish
			//------------------------------------------
			if (count_lp >= Min_Nr_LP) {
				// Calculate error of scalar channel (maximum value among timeslices)
				
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
						}  // switch(checkp[idx_cp].version)
						
						// If it reaches to the Maximum iterations, checkout
						if (count_lp == Max_Nr_LP) {
							chkout = true;
							checkp[idx_cp].chkout_order = -1;
						}
						
						if (chkout) {
							// checkout
							checkout(count_lp,TrM_inv, checkp[idx_cp].OutFileName,link_dirs,link_max,
									  checkp[idx_cp].chkout_order, Restarted,phases);
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
	
	delete phases;
	
	//-----------------------------
	// End program
	//-----------------------------
	snoop.stop();
	QDPIO::cout << "DISCO: total time = " << snoop.getTimeInSeconds() << " secs"
	<< std::endl;
	
	QDPIO::cout << "DISCO: ran successfully" << std::endl;
	QDPIO::cerr << "Successfully Run " << std::endl;
	QDP_abort(1);
	END_CODE();
	
	Chroma::finalize();
	exit(0);
}
