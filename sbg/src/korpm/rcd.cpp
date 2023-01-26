/***********************************************************************
 *
 * COPYRIGHT CSIC:  2011-2022
 *
 ***********************************************************************
 *
 * FILE: 				RCD.CPP
 *
 * DATE OF CREATION: 	August 31th, 2012
 * LAST MODIFIED:	 	Septmber 19th, 2022
 * AUTHOR (S):  		- Pieter Jozef Chys (2012-2015)
 * 						- Jose R. Lopez-Blanco (2015-2019)
 *                      - Pablo Chacon (2012-2022)
 *
 ***********************************************************************
 *
 * DESCRIPTION:	Random Coordinate Descent (RCD) loop closure
 * 				source code
 *
 **********************************************************************/

#include "parameters.h"

#ifdef MPI_RELEASE
#include <mpi.h>
#endif
// Notes to properly compile MPI release in sherpa machine:
// RCD program
// $mpicom = "mpirun -np "; // this one should work if everything was properly configured...
// $mpicom = "/usr/bin/orterun -np "; // for GNU's gcc
// IMPORTANT NOTE FOR SHERPA:
// Add -I/home/mon/mpi/include to Intel-MPI release and Compile it with /usr/local/openmpi/bin/mpiCC

#include "libga/include/rcd_stuff.h" // Mon: Now all includes are in rcd_stuff.h header...
#include "libga/include/timer.h"	 // I don't know why this must be here... otherwise multi-declaration errors appear.

#define LENGTH_LINE 500
#define LENGTH_SEQ 100
#define LENGTH_FILE 100

using namespace TCLAP;

int counter;
char ramafile[LENGTH_LINE]; // Input MAP file
char h3file[LENGTH_LINE]; // Input MAP file

char file_energy[LENGTH_LINE]; // Input Energy file (ICOSA or KORP)
char file_pdb[LENGTH_LINE]; // Input PDB file
char file_name[LENGTH_LINE]; // dummy text
bool load_ramap = false; // = true, Loads Ramachandran maps (Dunbrack), sequences, and related stuff...
float rama_thr =  0.9;     // rama cutoff not anchors or interkink
float rama_thrA = 0.9;     // rama cutoff anchors
float rama_thrK = 0.9;     // rama cutoff kink residues

float kink_mergeA = 0.99;
float kink_merge = 0.7;
float nokink_merge = 0.9;


float postcheck = -1.0;
int size=72; // 2D PDF side length
float step=5; // increment --> Depends on Dunbrack PDFs data set
float *****pdfs; // Pointer to pointer array
float *acc=NULL;
int cen,sid,lef,rig;
float constant = ((180.0/M_PI)/step);
float ratio = 0.0; // Ratio to verify the Dunbrack inside checking. Reinitialize in each PDB.
float native_energy = 0.0; // Native loop energy (Dunbrack-derived energy)
float native_loco = 0.0; // Native loop energy (LOCO-ICOS)
bool opt_condition = false;
int ramachandran_inside_aux = 0;
int ramachandran_check_aux = 0;
bool changed_thr = false;
bool output_debug = false;
bool output_dihedrals = true; // Output dihedral angles for PyRosetta...
bool bench_rmsd = false; // true= Considers non-native workflow but just using native-coordinates for rmsd benchmarking
bool mutate = false; // true= Uses input file sequence
bool self_clash = true;
bool flat_switch = false; // true= Uses a flat distribution [-1:1] for random number generator ("normRand")
bool exact_rama = true; // Enables/Disables exact Ramachandran maps instead of a binarized version
bool save_closed = true; // true= save all Closed loops
bool save_initial = false; // true= save all Initial loops
bool inside_rand = true; // Random Phi or Psi selection (according to given PDF distribution) after non-allowed Rama in inside RCD
int sidechain = 1; // 2= Include Side-chain (CB, H, HA, CG, SG) in the clash map of environment.
float pass_factor = 10000; // Factor to adjust which structures will be accepted after each RCD trial.
bool server = false; // Enables some server-related output
float omega = 179.10; // Omega dihedral angle (average value from Fig.1 in Dunbrack's paper: Berkholz PNAS 2012 http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3258596/)
bool delHydrogens_switch = true; // Delete hydrogens
bool delHeteros_switch = false; // Delete heteroatoms
bool delWaters_switch = true; // Delete waters
bool logbootstrap_switch = false; // Experimental log-min_rmsd representation
// Linear closure stuff
bool debug_linear = false;
bool linear_closure = false; // Enables/Disables the Analytic Linear Closure method
int linear_iters = 10; // Maximum number of linear closure iterations
double linear_cutoff = 0.001; // Below this threshold the loop is considered closed in the Analytic Linear Closure method
double linear_delta = 1.0; // Increment (ratio) to apply the linear displacement "X" to the dihedrals.
double linear_maxang = 0.017453293; // 0.087266465; // Maximum angular displacement [rad]. Amplitude will be reduced to match this value if any angle wants to surpass this maximum.
// EM-maps stuff
bool em_switch = false; // Enable/Disable EM-map loops scoring (by cross-correlation)
float em_thr = 0.0; // EM-map threshold (only voxel values ABOVE this threshold will be considered)

CRandomMersenne *rg; // Mersenne Twister global object
float unitsx;
float unitsy;
float unitsz;
int limitx;
int limity;
int limitz;
int stepy;
int stepz;

// Compute the derivatives matrix ("dr/d_theta") of the positions (r) of the 3 anchor atoms (N,CA,C) w.r.t. all dihedrals (theta)
//  der --> The derivatives matrix (#theta x 9 --> 9= 3 anchor atoms x 3 x,y,z coordinates). Set der=NULL to allocate memory.
//  co --> The N,CA,C coordinates (or C,CA,N in backward)
//  num_res --> number of residues including anchors
//  reverse --> (OPTIONAL) =0 forward, =1 backward
int der_anchor(double **p_der, double **co, int num_res, int reverse = 0);

// Compute the displacement vector B between the mobile Ct-anchor and the fixed one (Nt-anchor in reverse).
// It is the same for "forward" and for "backward" (reverse independent)
// B = [ b(i) ] = r_fix(i) - r_mob(i)      <--- Where "i" represents the N,CA,C atom indices of the Ct-anchor
//   b --> Displacement vector B (size= 9 x 1). Memory must be already allocated.
//   co --> The N,CA,C coordinates (or C,CA,N in backward) of the mobile loop
//   co_i --> The N,CA,C coordinates (or C,CA,N in backward) of the N/Ct-anchor (just 3 atoms coordinates)
//   num_res --> number of residues including anchors
double disp_b(double *b, double **co, double **co_i, int num_res);

// Get the maximum value
float getmax(float *a, int n);

// Get the minimum value
float getmin(float *a, int n);

// Load "n" random samples from "in" into "out" array without repetition
void randomSample(double *in, int n_in, float *out, int n_out);

// Compute the average of a float array "data" of "ndata" elements
float average(float *data, int ndata);

// Compute the sigma of a float array "data" of "ndata" elements. Optionally, the average can be provided if already known.
float sigma(float *data, int ndata, float avg=0.0);

// Compute the variance of a float array "data" of "size" elements. Optionally, the "mean" can be provided if already known.
float variance(float *data, int size, float mean=0.0);

// Linear regression: y = b + m·x  (correlation = r)
void linear_regression(float *x, float *y, int n, float *b, float *m, float *r);


extern "C" {

//SUBROUTINE DGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,
//$                   WORK, LWORK, IWORK, INFO )
//*  -- LAPACK driver routine (version 3.2.2) --
//INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
//DOUBLE PRECISION   RCOND
//INTEGER            IWORK( * )
//DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), S( * ), WORK( * )
void dgelsd_(int *nrows, int *Nmod, int *nrhs, double *Aop, int *lda,
		double *wop, int *ldb, double *svals, double *rcond, int *rank,
		double *work, int *lwork, int *iwork, int *info);
}


int main(int argc, char** argv)
{
	int32 seed = time( 0 );
	//int32 seed = 10;
	init_aminoacids(ICM,pdb);

	// ------------------------------------------------
	// PROGRAM VARIABLES
	// ------------------------------------------------
	int map;
	int reverse;	    // to set algorithm in forward or backward mode
	int reverse_input;  // to set algorithm in forward or backward mode

	bool rmsd_O = true; // final rmsd with native loop include O
	bool verb = false;
	bool save_movie = false;
	bool sim_text = false; // =false, output similarity files in plain-text, otherwise binary format to save disk space.
	bool simIC_switch = false; // =true, enable IC-RMSD similarity computation
	bool simCC_switch = false; // =true, enable CC-RMSD similarity computation
	bool no_ntpsi = true; // =false, optimize first Psi (Nt anchor's Psi angle)
	bool no_ctphi = false; // Enable (false) or Disable (true) Ct-anchor Phi optimization (watch out in Rosetta it is disabled)

	inbetween	  = 1; /* all atoms coordinates between the anchors are normally present; if not this       will become zero*/
	profile_counter = 0;
	randomize	  = 1; /*set to zero: inactive; set to one: active and randomization at each loop closure*/
	/*by default inactive according to reference experiment*/
	switching_side  = 1; /*by default the orientation should be normally forward in the algorithm and not combined forward-reverse */
	switching_side2 = 1; /*by default the orientation should be normally forward in the algorithm and not combined forward-reverse (2nd protocol)*/
	randbond        = 1; /* randbond = 1: bond lengths are being randomized; = 0 no randomization of bond lengths */
	randval	  = 1; /* randval = 1:  valence angles are being randomized; = 0 no randomization of valence angles */
	float omega_sigma = 0; // Standard deviation for Omega angles
	float Ctomega[4], Ntomega[4];

	cutoff 	   = 3;
	cutoffp	   = 3;		/*standard cutoff as like cutoff*/
	kindex 	   = 3;
	kindex2          = 3;
	counter_clash    = 0;
	atomshift	     = 3; // 3

	double factor_co = 2.89; // VdW radius^2
	double factor_ocb = 2.25; // VdW radius^2
	double factor_clashes = 1.0;
	double factor_cacb = 1.0;
	float kick_factor = 0.3; // Factor to trigger "kick" (the number of unsuccessful trials should be greater than factor times the number of loop atoms)
	float kick_thr = 0.001; // RMSD improvement factor to trigger "kick"
	float kick_cut = 0.35; // RMSD cutoff that enables increasing the kick counter
	//	float lc_msd = 0.25; // Analytic Loop-Closure option
	bool enable_kick = true; // Enable "kick" stage
	bool enable_lc = true; // Enable analytic (exact) loop-closure
	bool enable_vib = true;
	bool rasp_switch = false; // Enable RASP repacking to obtain heavy-atoms loop models
	float rasp_time=0.0; // total timer for RASP
	float bump_cutoff; // Bump filtering cutoff sigma
	float bump_rate = 1.0; // Bump filtering cutoff rate
	bool bump_filter = false; // Enable bumps filtering using PD2 bump energy
	bool presampling_switch = false; // Enable pre-sampling to estimate energy distribution
	bool presampling_continue = true; // To continue standard iterations...
	int presampling = 200; // Number of iterations of pre-sampling to estimate energy distribution
	// PD2 bump filter (MacDonald et al. 2013)
	double energy = 0.0; // energy
	double energy_intra = 0.0;
	double energy_extra = 0.0;
	//	bool loco_switch = false; // =true to enable LOCO-ICOS energy calculation
	int emodel = 1; // Potential energy model (=1 ICOSA, =5 KORP)
	float bonding_factor = -1.0; // Bonding energy factor (typically 0.3)
	bool notbest_switch = false;  // =true to enable random selection of not-best LOCO-ICOS energy loops
	float *loop_coord=NULL; // linear array of coordinates of the mobile N, CA, C loop atoms.
	int loco_best = 0; // 0= --> disabled, !=0 --> number of best solutions to retain
	int loco_notbest = 0; // 0= --> disabled, !=0 --> number of not-best solutions to retain
	bool trama_switch = false; // Enable (true) or Disable (false) the Terminal Ramachandran threshold relax
	int nrama = 0; // Number of terminal residues to be considered in Rama. relaxing...
	bool rama_switch = false; // Set "true" to activate Rama potential (added as bonding contribution to the selected "emodel" potential)
	bool h3_switch = false; // Set "true" to activate Rama h3 potential


	float rama_factor = 10.0; // Ramachandran energy factor (typically 10)
	int rama_model = 3; // Ramachandran energy model	reverse             = 0;// Direction initialization: default to forward mode:on = 1; off
	ramachandran_check  = 2; /*	Ramachandran initialization: on = 1; off = 0 */
	ramachandran_inside = 2; /*	Ramachandran inside protocol:     		    on = 1; off = 0 */
	include_hydrogens   = 0; /*	Include hydrogens in selected loop environment:     on = 1; off = 0 */
	post_check	        = 0; /*	Mainchain check backbone with grid		  	            */
	no_loop	            = 0; /*	Input with only anchors:                           on =  1; off = 0 */
	float rmsd_thr = 9999; // RMSD threshold for "begin" conformations (all rmsds should be less than rmsd_thr*nres)
	int present_seq = 0 ;
	int cg_mode = 3; // mode of operation: 0=(N,CA,C), 1=(N,CA,C,O), 2=(N,CA,C,CB), 3=(N,CA,C,O,CB)
	int protocol; // mode of coordinates: 0= Native, 1=Native+Random, 2=No loop coordinates
	switch_stop    = 0;   // default (0) = no stopping during program execution; (1) stopping during the experiments
	switch_file    = 1;   // default (1) = saving of output and result files on
	switch_solutions = 1; // default (1) = save the loop closure solutions
	switch_situs	 = 0; // default (1) = save situs map of the protein grid
	grid_backbone    = 0; // default (0) = all atoms put on the grid; (1) only backbone atoms used
	seq_initgrid     = 0; // default (0) = putting the loop closure on grid during initialization
	switch_dof       = 1; // default (0) = extra DoF off
	switch_bond      = 1; // default (1) = random bond choice selection
	switch_perturbation = 0; // def (1) = switching on of a perturbation cycle at start of algorithm
	native 	= 0;   	//	default (0) = not using native loop; (1) = using native loop conformation
	double shift_rev[3], rdump;

	randcount     =	4;
	cycles 	=    1000;  //	number of cycles to be computed for one loop structure to test it
	//	Make it large enough: cycles = 1 gives floating point exception!!!!!!
	//	Default cycles: 1000
	number_path   =	cycles - 5;
	max_path     = 	5000; // -m
	path_slope     = 	1000; // -p
	max_perturbation = 	100;
	rmsd_crit 	    = 	0.1; // -d
	stepsize 	    =   0.25;
	softness	    =	0.9; // --softness
	bondptb 	    = 	1.0; // multiplicative factor for normal distribution. N(bondptb,sigma) = sigma * N(bondptb,1)
	valptb  	    = 	1.0; // 0.9625 (old value)   N(valptb,sigma) = sigma * N(valptb,1)
	// from Berkholz paper (see chain_manip.cpp); other alternative is value 1.0 to account for 95% confidence interval
	int naa_fix = 0;   // number of aa to be fixed

	// CREATE ARRAYS WITH VALID PHI AND PSI ANGLES
	int size_rama;			/*number of Ramachandran angles in valid range*/
	int index_rama;			/*index to array with Ramachandran valid angles*/
	double range_phi;		/*valid numerical range PHI angles*/
	double range_psi;		/*valid numerical range PSI angles*/
	double step_phi;		/*computed step size in PHI range*/
	double step_psi;		/*computed step size in PSI range*/
	double current;			/*auxiliary current dihedral angle used in Ramachandran inside-filter approach*/

	size_rama 	= 	300;	/*resolution of 1 degree seems over e.g. more or less 180 degrees is interesting*/
	current		=	0;		/*extra initializiation prior to actual initialization*/
	range_phi	= 	UPPERPHI 	- 	LOWERPHI;
	range_psi	= 	UPPERPSI	-	LOWERPSI;
	step_phi	=	range_phi/double(size_rama - 1);
	step_psi	=	range_psi/double(size_rama - 1);
	double *angle_phi 	= 	new double[size_rama];	/*create array Ramachandran valid PHI angles*/
	double *angle_psi	=	new double[size_rama];	/*create array Ramachandran valid PSI angles*/

	for(int i = 0; i < size_rama; i++)				/*create actual Ramachandran angles*/
	{
		angle_phi[i]	=	LOWERPHI + i*step_phi;
		angle_psi[i]	=	LOWERPSI + i*step_psi;
	}

	// PRINT SETTINGS
	char    *testfile;			/*char array that will contain name input file with pdb names and loop residue range*/
	char folder_later[LENGTH_LINE];
	folder_later[0]='\0';
	bool folder_name = false;
	string line(80,'-');
	testfile = 	new char[LENGTH_LINE];


	// H3 stuff
	int h3CtShift=0; // Ct-term shift Ct-anchor is Cys in pos 92
	int h3NtShift=0; // h3ct_shift", "Ct-term shift Nt-anchor

	std::string temp;
	CmdLine cmd(argv[0],"RCD program is a fast loop-closure modeling tool based on an improved version of our RCD method "
			"[Chys and Chacón (2013)]. Accurate loop ensembles can be easily generated in a few minutes. You just need the "
			"atomic coordinates (PDB file v3.x, chain ID required), chain ID, and indices of the start/end residues. "
			"The presence in the PDB file of the N- and C-terminal anchor residues (indices Start-1 "
			"and End+1) coordinates is mandatory. Notice that all residues from \"Start\" to \"End\" indices (inclusive) "
			"will be modeled from scratch.\n"
			"Note that most of the above options are developer's options for fine tuning and debugging of the algorithm. "
			"Only a few options are usually required for optimal performance. "
			"Feel free to contact us to get some additional help!\n"
			"*If you would like to refer our work, please use:\n"
			"Server and improved method:\n\tLópez-Blanco JR, Canosa-Valls AJ, Li Y, and Chacón P (2016). "
			"RCD+: Fast loop modeling server. NAR (DOI: 10.1093/nar/gkw395).\n"
			"Original method:\n\tChys P and Chacón P (2013). "
			"Random coordinate descent with spinor-matrices and geometric filters for efficient loop closure. "
			"J. Chem. Theory Comput. 9:1821-1829.\n"
			"Thanks for using RCD!", "1.60" );

	try {
		//
		// Define required arguments no labeled
		//
		UnlabeledValueArg<string> pdb_in("txt","Plain text input file. Use the following format, one loop per line:\n<PDB_file> [Start_index] [End_index] [Chain_ID] [Loop_sequence]\nFor example: \"135l.pdb    84  91 X LSSDITAS\".","default","input_list");
		cmd.add( pdb_in );

		ValueArg<float> EMmap("","em", "EM-map threshold (values above will be considered, default= 0.0). It also enables using an EM-map to score loops. "
				"The base name of the MRC map must correspond to <PDB_file> basename.", false, 0.0,"float");
		cmd.add( EMmap );

		SwitchArg DelHeteros("","delete_heteros", "Delete Hetero-atoms, including waters (default=disabled).", true);
		cmd.add( DelHeteros );
		SwitchArg KeepWaters("","keep_waters", "Disables Water molecules deletion (default=disabled).", true);
		cmd.add( KeepWaters );
		SwitchArg KeepHydrogens("","keep_hydrogens", "Disables Hydrogen atoms deletion (default=disabled).", true);
		cmd.add( KeepHydrogens );

		ValueArg<float> frmsd("d","rmsd", "Minimum anchor RMSD distance in Angstroms (default=0.1)",false,0.1,"float");
		cmd.add( frmsd);

		ValueArg<int> fn("n","nloops", "Number of loops sampled.", false,1000,"int");
		cmd.add( fn );

		ValueArg<int> fm("m","max_path", "Basal maximum number of RCD iterations. (default = 5000)", false, 5000,"int");
		cmd.add( fm );

		ValueArg<float> slope("p","path_slope", "Slope of the maximum number of RCD iterations (Max number = max_path + (number_of_residues-4) * slope). (default = 1000)", false,1000,"float");
		cmd.add( slope );

		SwitchArg LogBootstrap("","logbootstrap", "Experimental log-min_rmsd representation. (default= disabled)", false);
		cmd.add( LogBootstrap );

		ValueArg<float> LinearMaxang("","linear_maxang", "Maximum angular displacement [deg]. Amplitude will be reduced to match this value if any angle wants to surpass this maximum. (default = 5)", false, 5.0,"float");
		cmd.add( LinearMaxang );

		ValueArg<float> LinearDelta("","linear_delta", "Increment (ratio) to apply the linear displacement X to the dihedrals. (default = 1.0)", false, 1.0,"float");
		cmd.add( LinearDelta );

		ValueArg<float> LinearCutoff("","linear_cutoff", "Below this threshold the loop is considered closed in the analytic linear closure method [A] (default = 0.001)", false, 0.001,"float");
		cmd.add( LinearCutoff );

		ValueArg<int> LinearIters("","linear_iters", "Maximum number of linear closure iterations (default = 20)", false, 20,"int");
		cmd.add( LinearIters );

		ValueArg<int> Movie_iter("","movie_iter", " Iteration to save a trajectory", false, 0,"int");
		cmd.add( Movie_iter );

		SwitchArg linear("","linear", "Enable the analytic linear closure method. (default= disabled)", false);
		cmd.add( linear );

		SwitchArg nat("f","native", "Loop already present in the pdb (maintain native original angles and lengths).", false);
		cmd.add( nat );

		SwitchArg Randomize("r","randomize", "Randomize bond lengths and valence angles. (small perturbation)", false);
		cmd.add( Randomize );

		SwitchArg Bench("","bench", "Considers non-native workflow but just using native-coordinates for RMSD benchmarking. (disabled by default)", false);
		cmd.add( Bench );

		SwitchArg Output_Debug("","output_debug", "Output debugging files.", false);
		cmd.add( Output_Debug );

		SwitchArg FlatRand("","flat_rand", "Use a flat random distribution in random number generator. (disabled by default)", false);
		cmd.add( FlatRand );

		ValueArg<float> OmeSigma("","omega_sigma", "Sigma (s) for omega dihedral angles (default=0 deg, but 6.3 deg in Berkholz et al. PNAS 2012 109:449-53). Note: 1s=68.3%, 2s=95.5%, 3s=99.7%", false, 0,"float");
		cmd.add( OmeSigma );

		// Omega dihedral angle (average value from Fig.1 in Dunbrack's paper: Berkholz et al. PNAS 2012 http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3258596/)
		ValueArg<float> Omega("","omega", "Omega dihedral angle (default=179.1 degrees, taken from: Berkholz et al. PNAS 2012 109:449-53)", false, 179.1,"float");
		cmd.add( Omega );

		ValueArg<float> Post("l","post", "Post-checking threshold. (Dunbrack)", false, -1,"float");
		cmd.add( Post );

		ValueArg<float> RMSD_thr("","rmsd_thr", "RMSD threshold for \"Begin\" conformations (all RMSDs should be less than rmsd_thr*nres). (default=disabled)", false, 9999,"float");
		cmd.add( RMSD_thr);

		ValueArg<float>Softness("","softness", "VdW radius factor for clashes map generation. (default=0.9)", false, 0.9,"float");
		cmd.add( Softness);

		ValueArg<float> Stepsize("","stepsize", "Voxel size in clashes map. (default=0.25)", false, 0.25,"float");
		cmd.add( Stepsize);

		ValueArg<int> switching("s","switching", "Switching protocol \n\t 0: forward \n\t 1: reverse \n\t 2: random  \n\t 3: failure based (Default)", false, 3, "int");
		cmd.add( switching );

		ValueArg<int> CG_mode("c","", "Coarse-graining level for atomic model\n\t 0: N,CA,C\n\t 1: N,CA,C,O\n\t 2: N,CA,C,CB\n\t 3: N,CA,C,O,CB (default)", false, 3, "int");
		cmd.add( CG_mode );

		ValueArg<float> Rama_thr("t","threshold", "Threshold for neighbor-dependent Ramachandran filter. For example, using 0.90 RCD accounts for the 90% of probability (default=0.98)", false,0.98,"float");
		cmd.add( Rama_thr );

		ValueArg<int> Ins("i","inside", "Inside filter: \n\t 0: Free \n\t 1: Original \n\t 2: Dunbrack (default) \n\t", false,2,"int");
		cmd.add( Ins );

		ValueArg<int> Begin("b","begin", "Pre filter: \n\t 0: Free \n\t 1: Original \n\t 2: Dunbrack (default) \n\t", false,2,"int");
		cmd.add( Begin );

		ValueArg<float> Factor_CO("","factor_co", "Backbone atoms intra-clash factor (square of inter-atomic distance). (default=2.89)", false, 2.89,"float");
		cmd.add( Factor_CO );

		ValueArg<float> Factor_OCB("","factor_ocb", "O and CB atoms intra-clash factor (square of inter-atomic distance). (default=2.25)", false, 2.25,"float");
		cmd.add( Factor_OCB );

		ValueArg<float> Factor_Clashes("","factor_clashes", "Factor that multiplies CO and CB factors for intra-clash detection. (default=1.0)", false, 1.0,"float");
		cmd.add( Factor_Clashes );

		ValueArg<float> Factor_CACB("","factor_cacb", "CA-CB distance multiplicative factor to define atomic model (set to <=0 for constant standard distance). (default=1.0)", false, 1.0,"float");
		cmd.add( Factor_CACB );

		ValueArg<int> Sidechain("","sidechain", "Side-chain atoms considered in the clash detection map: \n\t 0: N, CA, C and O \n\t 1: N, CA, C, O and CB  (default)\n\t 2: CB, H, HA, CG and SG\n\t 3: All heavy atoms present in PDB", false, 1,"int");
		cmd.add( Sidechain );

		SwitchArg SimilarityIC("","simic", "Outputs IC similarities between loops.", false);
		cmd.add( SimilarityIC );

		SwitchArg SimilarityCC("","simcc", "Outputs CC similarities between loops.", false);
		cmd.add( SimilarityCC );

		SwitchArg SimText("","sim_text", "Enables plain text format for similarities between loops instead of binary. (default=disabled)", false);
		cmd.add( SimText );

		SwitchArg RaspPost("","rasp", "Enable RASP repacking to obtain heavy-atoms loop models. WARNING: links to bbdep11.bin and RASP.ini files are mandatory! (default=disabled)", false);
		cmd.add( RaspPost );

		SwitchArg Rmsd_O("","rmsdO", "Disable O in the computation of rmsd with native. (default=enables)", true);
		cmd.add( Rmsd_O );

		ValueArg<float> Pass_Factor("","pass_factor", "Factor to adjust which structures will be accepted after each RCD trial (the higher the easier). (default=10000)", false, 10000,"float");
		cmd.add( Pass_Factor );

		ValueArg<float> Kick_Factor("","kick_factor", "Factor to trigger \"kick\" (the number of unsuccessful trials should be greater than factor times the number of loop atoms). (default=0.3)", false, 0.3,"float");
		cmd.add( Kick_Factor );

		ValueArg<float> Kick_Thr("","kick_thr", "RMSD improvement factor to trigger \"kick\". (default=0.001)", false, 0.001,"float");
		cmd.add( Kick_Thr );

		ValueArg<float> Kick_Cut("","kick_cut", "RMSD cutoff that enables increasing the kick counter. (default=0.35)", false, 0.35,"float");
		cmd.add( Kick_Cut );

		//		ValueArg<float> LC_MSD("","lc_msd", "Mean Square Deviation to accept an analytic Loop-Closure conformation. (default=0.25)", false, 0.25,"float");
		//		cmd.add( LC_MSD );

		SwitchArg Disable_Kick("","no_kick", "Disable \"kick\" stage. (default=disabled)", false);
		cmd.add( Disable_Kick );

		ValueArg<int> Naa_fix("","nfix", " Number of aa randomly fix. (default=0)", false, 0,"int");
		cmd.add( Naa_fix );

		SwitchArg Server("","server", "Enable RCD-Server output. (default=disabled)", false);
		cmd.add( Server );

		ValueArg<float> ValFactor("","valence_factor", "Multiplicative factor for the sigma of valence and omega angles. (default=1)", false, 1.0,"float");
		cmd.add( ValFactor );

		ValueArg<float> BondFactor("","bond_factor", "Multiplicative factor for the sigma of bond lengths. (default=1)", false, 1.0,"float");
		cmd.add( BondFactor );

		ValueArg<float> BumpRate("","bump_rate", "Bump filtering cutoff rate. (default=1.0)", false, 1.0,"float");
		cmd.add( BumpRate );

		SwitchArg BumpFilter("","bump_filter", "Enable bumps filtering using PD2 bump energy. (default=disabled)", true);
		cmd.add( BumpFilter );

		ValueArg<int> PreSampling("","presampling", "Number of iterations of pre-sampling to estimate energy distribution (default=200)", false, 200,"int");
		cmd.add( PreSampling );

		ValueArg<float> BondingFactor("","bf", "Bonding energy factor, typically 0.3 goes fine. (Default=0.3)",false,0.3,"float");
		cmd.add( BondingFactor );

		ValueArg<int> LocoBest("","loco_best", "Number of lowest energy loops to be further stored and processed (<nloops). (disabled by default)", false, 0,"int");
		cmd.add( LocoBest );

		ValueArg<int> LocoNotBest("","loco_notbest", "Number of Not-lowest ICOSA-energy solutions that will be randomly selected for storage and processing (<nloops). (disabled by default)", false, 0,"int");
		cmd.add( LocoNotBest );

		SwitchArg NoSaveIC("","nosave_ic", "Disable Internal Coordinates output (*_dh.txt, *_val.txt, and *_len.txt files). (default=disabled)", true);
		cmd.add( NoSaveIC );

		SwitchArg NoSaveClosed("","nosave_closed", "Disable all closed loops output as a Multi-PDB file. (default=disabled)", true);
		cmd.add( NoSaveClosed );

		SwitchArg SaveInitial("","save_initial", "Save all initial loops in Multi-PDB format. (default=disabled)", true);
		cmd.add( SaveInitial );

		SwitchArg InsideBorder("","inside_border", "Phi or Psi are forced to the Rama border (according to PDF) after non-allowed Rama in inside loop. (default=disabled)", true);
		cmd.add( InsideBorder );

		SwitchArg NtPsi("","nt_psi", "Include the Psi angle of the Nt anchor in the conformational search if Nt's Oxigen is present. (default=disabled)", true);
		cmd.add( NtPsi );

		SwitchArg NoCtPhi("","no_ctphi", "Disable C-terminal Phi dihedral angle optimization in conformational search (to do it like Rosetta). (default=disabled)", true);
		cmd.add( NoCtPhi );


		ValueArg<int> AtomShift("","atomshift", "Number of contiguous mainchain (co) atoms to ignore intra-chain contact. Typically, 3= <=5 bonds, 4= <5 bonds. (default=4)", false, 3,"int");
		cmd.add( AtomShift );

		//	SwitchArg Disable_LC("","no_lc", "Disable \"Analytic Loop-Closure\" (default=disabled)", false);
		//	cmd.add( Disable_LC );


		ValueArg<float> RamaFactor("","rf", "Ramachandran energy factor. (Default=6.0)",false,6.0,"float");
		cmd.add( RamaFactor );

		ValueArg<int> RamaModel("","rama_model", "Ramachandran energy model: "
				"\n\t 0 = Pc/1, "
				"\n\t 1 = Pc/Pt (weighted Pt), "
				"\n\t 5 = Pc/Pt1, "
				"\n\t 2 = Plcr/Pt (weighted Pt), "
				"\n\t 6 = Plcr/Pt1, "
				"\n\t 3 = Plcr/Pc (default)"
				"\n\t 4 = Plcr/1",false,3,"int");
		cmd.add( RamaModel );

		ValueArg<std::string> RamaFile("x","rama_file", "Dunbrack's Probability Density Function.",false,"pdf","string");
		cmd.add( RamaFile );

		ValueArg<std::string> h3File("","h3", "H3 Rama Probability Density Function.",false,"pdf","string");
		cmd.add( h3File );

		ValueArg<int> H3NtShift("","h3nt_shift", "Nt-term kink shift considering Ct-anchor is Cys in pos 92", false, 0,"int");
		cmd.add( H3NtShift );

		ValueArg<int> H3CtShift("","h3ct_shift", "Ct-term kink shift considering Nt-anchor is W in pos 104", false, 0,"int");
		cmd.add( H3CtShift );


		ValueArg<float> Kink_merge("","kink_merge", "Factor to merge rama maps between kink residues no anchors (default=0.90)", false, 0.90,"int");
		cmd.add( Kink_merge);

		ValueArg<float> Kink_mergeA("","kink_mergeA", "Factor to merge rama maps kink anchors (default=0.90)", false, 0.90,"int");
		cmd.add( Kink_mergeA);

		ValueArg<float> Nokink_merge("","nokink_merge", "Factor to merge rama maps betweens kink residues (default=0.1)", false, 0.1,"int");
		cmd.add( Nokink_merge);

		ValueArg<float> Rama_thrK("","tkink", "Kink threshold for neighbor-dependent Ramachandran distribution (default=0.70)", false,0.70,"float");
		cmd.add( Rama_thrK );

		ValueArg<float> Rama_thrA("","tkinkA", "anchor Kink threshold for neighbor-dependent Ramachandran distribution (default=0.70)", false,0.70,"float");
		cmd.add( Rama_thrA );

		ValueArg<int> Nterm_rama("","nterm_rama", "Number of N/C-terminal residues to relax Ramachandran filter. 0= Anchors-only, 1= Anchors and first/last loop residues, and so on. Threshold used in these residues should be set with --term_ramaA option. (default=1)", false, 1,"int");
		cmd.add( Nterm_rama );

		ValueArg<float> Term_ramaA("","term_ramaA", "Ramachandran threshold for the N/C-terminal anchors residues. (default=0.99)", false, 0.99,"float");
		cmd.add( Term_ramaA );

		//		ValueArg<float> Term_rama("","term_rama", "Ramachandran threshold except for the N/C-terminal residues. Only used if --nterm_rama is set. (default=0.99)", false, 0.99,"float");
		//		cmd.add( Term_rama );


		//		ValueArg<std::string> Loco("","loco", "Introducing a ICOSA energy data file (loco.score) enables ICOSA energy calculations.",false,"loco","string");
		ValueArg<std::string> EnergyFile("","energy_file", "Energy file data: ICOSA (loco.score) or KORP energy data file.",false,"energy","string");
		cmd.add( EnergyFile );

		ValueArg<int> Emodel("e","energy_model", "Energy model: \n\t 0 = None (just loop closure)"
				"\n\t  1 = ICOSA, "
				"\n\t  5 = KORP (default)"
				"\n\t  6 = Rama"
				"\n\t 56 = KORP + Rama",false,5,"int");
		cmd.add( Emodel );

		SwitchArg BinaryRama("","binary_rama", "Enables using binary Ramachandran maps instead of an exact version (default=disabled)", false);
		cmd.add( BinaryRama );

		ValueArg<std::string> Folder_name("o","name", "Output directory basename. By default, current date will be used as basename (e.g. run__11_May_2016__13h_23_06).",false,"./","string");
		cmd.add( Folder_name );
		//
		//*******************************************************

		// Parse the command line.
		cmd.parse(argc,argv);
		strcpy(testfile,((temp=pdb_in.getValue()).c_str()));


		//		fprintf(stderr,"ramafile= %s\n",ramafile);
		strcpy(file_energy,((temp=EnergyFile.getValue()).c_str()));
		strcpy(folder_later,((temp=Folder_name.getValue()).c_str())); // Gets output file name
		//
		//*******************************************************
		rmsd_crit=frmsd.getValue();
		cycles=fn.getValue();
		path_slope=slope.getValue();
		max_path = fm.getValue();
		ramachandran_inside = Ins.getValue();
		ramachandran_check = Begin.getValue();
		switching_mode = switching.getValue();
		rmsd_thr = RMSD_thr.getValue();
		postcheck = Post.getValue();
		softness = Softness.getValue();
		stepsize = Stepsize.getValue();
		cg_mode = CG_mode.getValue();
		output_debug = Output_Debug.getValue();
		factor_co = Factor_CO.getValue();
		factor_ocb = Factor_OCB.getValue();
		bench_rmsd = Bench.getValue();
		factor_clashes = Factor_Clashes.getValue();
		factor_co *= factor_clashes;
		factor_ocb	*= factor_clashes;
		factor_cacb= Factor_CACB.getValue();
		naa_fix=Naa_fix.getValue();
		bump_rate = BumpRate.getValue();
		presampling_switch = PreSampling.isSet();
		presampling = PreSampling.getValue();

		if(BumpFilter.isSet())
		{
			bump_filter = true;
			presampling_switch = true;
		}
		//		loco_switch = EnergyFile.isSet();
		emodel = Emodel.getValue();
		if( emodel == 0 && EnergyFile.isSet() ) // If old text loco.icos energy file is provided
			emodel = 1;

		if(Movie_iter.isSet())
		{
			number_path=Movie_iter.getValue();
			save_movie = true;
		}



		// Rama potential  stuff
		rama_model = RamaModel.getValue();
		rama_factor = RamaFactor.getValue();
		if(RamaFile.isSet())
		{
			strcpy(ramafile,((temp=RamaFile.getValue()).c_str()));
			rama_switch = true;
		}

		// Rama cut-off
		rama_thr = Rama_thr.getValue();  // -t option  or inter kink if H3
		rama_thrA = Rama_thrA.getValue(); // rama for the anchor with rama_thrA

		if( Nterm_rama.isSet() ) {
			nrama = Nterm_rama.getValue(); // Number of terminal rama residues to be relaxed
			trama_switch = true;
		}
		else {
			nrama=0;
			trama_switch =false;
		}

		// H3 stuff
		if(h3File.isSet())
		{
			strcpy(h3file,((temp=h3File.getValue()).c_str()));
			h3_switch = true;
			kink_merge=Kink_merge.getValue();
			kink_mergeA=Kink_mergeA.getValue();
			nokink_merge= Nokink_merge.getValue();
			rama_thrK = Rama_thrK.getValue();

		}

		if(H3CtShift.isSet())
			h3CtShift=H3CtShift.getValue();

		if(H3NtShift.isSet())
			h3NtShift=H3NtShift.getValue();




		if(emodel == 6 || emodel == 56)
			rama_switch = true;
		if( rama_switch && !RamaFile.isSet() ) // Rama requested but no Rama file found
		{
			fprintf(stdout,"Error! Please, introduce a valid Ramachandran PDFs file (dunbrack.bin) in --rama_file option.\n");
			exit(1);
		}

		if(BondingFactor.isSet())
			bonding_factor = BondingFactor.getValue(); // Bonding energy factor (typically 0.3)
		loco_best = LocoBest.getValue();
		if(loco_best > cycles)
		{
			fprintf(stderr,"Sorry, --loco_best (%d) can not be greater than -n (%d).\nForcing exit!\n",loco_best,cycles);
			exit(1);
		}
		if(loco_best == cycles)
		{
			fprintf(stderr,"Sorry, --loco_best (%d) is equal to -n (%d), thus --loco_best is automatically disabled!\n",loco_best,cycles);
			loco_best = 0;
		}
		notbest_switch = LocoNotBest.isSet();
		loco_notbest = LocoNotBest.getValue();
		if((loco_best+loco_notbest) > cycles)
		{
			fprintf(stderr,"Sorry, the sum (%d) of --loco_best (%d) and --loco_notbest (%d) can not exceed -n (%d).\nForcing exit!\n",loco_best+loco_notbest,loco_best,loco_notbest,cycles);
			exit(1);
		}

		if(Folder_name.isSet())
			folder_name = true;

		if (nat.isSet())
		{
			native=1;
			present_seq = 1;
		}
		else
		{
			native=0;
			present_seq = 0;
		}

		if(Randomize.isSet())
		{
			randbond = 1;
			randval = 1;
			randomize = 1;
		}
		else
		{
			randbond = 0;
			randval = 0;
			randomize = 0;
		}

		omega = Omega.getValue();
		omega_sigma = OmeSigma.getValue();

		simIC_switch = SimilarityIC.isSet();
		simCC_switch = SimilarityCC.isSet();
		sim_text = SimText.isSet();

		if(RamaFile.isSet())
			load_ramap = true; // Load rama-stuff...

		if( (ramachandran_check == 2 || ramachandran_inside == 2) )
		{
			if( !RamaFile.isSet() )
			{
				fprintf(stderr,"Please, introduce a Ramachandran map file!\n");
				exit(1);
			}
			if( !Rama_thr.isSet() && rama_thr<=0 && rama_thr>=1)
			{
				fprintf(stderr,"Please, introduce a valid threshold for Ramachandran map!\n");
				exit(1);
			}
		}

		kick_factor = Kick_Factor.getValue(); // Factor to trigger "kick" (the number of unsuccessful trials should be greater than factor times the number of loop atoms)
		kick_thr = Kick_Thr.getValue(); // RMSD improvement factor to trigger "kick"
		kick_cut = Kick_Cut.getValue(); // RMSD cutoff that enables increasing the kick counter
		enable_kick = !Disable_Kick.isSet(); // Enable "kick" stage
		// enable_lc = !Disable_LC.isSet(); // Enable analytic (exact) loop-closure
		//		lc_msd = LC_MSD.getValue(); // MSD threshold to accept a Loop-Closure conformation.
		rasp_switch = RaspPost.isSet(); // Enable RASP repacking to obtain heavy-atoms loop models
		server = Server.isSet();
		flat_switch = FlatRand.isSet();
		exact_rama = !BinaryRama.isSet();
		pass_factor = Pass_Factor.getValue();
		sidechain = Sidechain.getValue();

		if(NoSaveClosed.isSet())
			save_closed = false;
		else
			save_closed = true;

		if(NoSaveIC.isSet())
			output_dihedrals = false; // Disable output ICs for PyRosetta...
		else
			output_dihedrals = true; // Enable output ICs for PyRosetta...

		if(SaveInitial.isSet())
			save_initial = true;
		else
			save_initial = false;

		if(InsideBorder.isSet())
			inside_rand = false; // Random Phi or Psi selection (according to given PDF distribution) after non-allowed Rama in inside RCD

		if(Rmsd_O.isSet())
			rmsd_O = false; // Disable O-atom in RMSD computation

		if(cg_mode == 0 || cg_mode == 2)
		{
			fprintf(stdout,"rcd> Using N, CA and C atoms for RMSD computations. \nrcd>\n");
			rmsd_O = false; // Disable O-atom in RMSD computation
		}

		if( NtPsi.isSet() )
			no_ntpsi = false; // Enable first Psi (Nt anchor's Psi angle) optimization (default for all CGs)
		else
			no_ntpsi = true; // Disable first Psi (Nt anchor's Psi angle) optimization


		no_ctphi = NoCtPhi.isSet(); // Enable (false) or Disable (true) Ct-anchor Phi optimization (in Rosetta it is disabled)


		valptb = ValFactor.getValue(); // 0.9625 (old value)   N(valptb,sigma) = sigma * N(valptb,1)
		bondptb = BondFactor.getValue(); // multiplicative factor for normal distribution. N(bondptb,sigma) = sigma * N(bondptb,1)

		atomshift = AtomShift.getValue();

		linear_closure = linear.isSet(); // Enables/Disables the Analytic Linear Closure method
		linear_iters = LinearIters.getValue(); // Maximum number of linear closure iterations
		linear_cutoff = LinearCutoff.getValue(); // Below this threshold the loop is considered closed in the Analytic Linear Closure method
		linear_delta = LinearDelta.getValue(); // Increment (ratio) to apply the linear displacement "X" to the dihedrals.
		linear_maxang = LinearMaxang.getValue() * M_PI / 180; // 0.017453293 rad / deg: 5 deg = 0.087266465 rad   10 deg = 0.17453293

		logbootstrap_switch = LogBootstrap.isSet(); // Experimental log-min_rmsd representation

		if( EMmap.isSet() )
		{
			em_thr = EMmap.getValue(); // EM-map threshold (only voxel values ABOVE this threshold will be considered)
			em_switch = true; // Enable EM-map loops scoring (by cross-correlation)
		}

		if(DelHeteros.isSet())
			delHeteros_switch = true; // Delete heteroatoms
		if(KeepHydrogens.isSet())
			delHydrogens_switch = false; // Keep hydrogens
		if(KeepWaters.isSet())
			delWaters_switch = false; // Keep waters

	} catch ( ArgException& e )
	{
		cout << "  Error->" << e.error() << " " << e.argId() << endl;
	}


	//***********************************************************
	//
	int pid=0; // To be compatible with MPI

#ifdef MPI_RELEASE
	int np;
	MPI_Status status;
	MPI_Init (&argc, &argv);      /* starts MPI */
	MPI_Comm_rank (MPI_COMM_WORLD, &pid);        /* get current process id */
	MPI_Comm_size (MPI_COMM_WORLD, &np);        /* get number of processes */
	seed += pid;
	printf( "Process %d of %d awake! seed= %d\n", pid, np, seed );
#endif

	if(pid==0)
	{
		fprintf(stderr,"---------------------------------\n");
		fprintf(stderr,"USED PARAMETERS:\n");
		fprintf(stderr,"---------------------------------\n");
		fprintf(stderr,"\t COMMAND-LINE: ");
		for(int i=0;i<argc;i++)
			fprintf(stderr,"%s ", argv[i]);
		fprintf(stderr,"\n\t n: %d\n",cycles);
		fprintf(stderr,"\t max_path: %d\n",max_path);
		fprintf(stderr,"\t path_slope: %f\n",path_slope);
		fprintf(stderr,"\t t: %f\n",rama_thr);
		fprintf(stderr,"\t s: %d\n",switching_mode);
		fprintf(stderr,"\t c: %d\n",cg_mode);
		fprintf(stderr,"\t d: %f\n",rmsd_crit);
		fprintf(stderr,"\t rmsd_thr: %f\n",rmsd_thr);
		fprintf(stderr,"\t pass_factor: %f\n",pass_factor);
		fprintf(stderr,"\t kick_factor: %f\n",kick_factor);
		fprintf(stderr,"\t kick_thr: %f\n",kick_thr);
		fprintf(stderr,"\t kick_cut: %f\n",kick_cut);
		fprintf(stderr,"\t softness: %f\n",softness);
		fprintf(stderr,"\t bump_filter: %d\n",bump_filter);
		fprintf(stderr,"\t bump_rate: %f\n",bump_rate);
		fprintf(stderr,"\t presampling: %d\n",presampling);
		fprintf(stderr,"\t factor_clashes: %f\n",factor_clashes);
		fprintf(stderr,"\t factor_co: %f\n",factor_co);
		fprintf(stderr,"\t factor_ocb: %f\n",factor_ocb);
		fprintf(stderr,"\t enable_kick: %d\n",enable_kick);
		fprintf(stderr,"\t inside_rand: %d\n",inside_rand);
		fprintf(stderr,"\t exact_rama: %d\n",exact_rama);
		//		fprintf(stderr,"\t loco_switch: %d\n",loco_switch);
		fprintf(stderr,"\t emodel: %d\n",emodel);
		fprintf(stderr,"\t loco_best: %d\n",loco_best);
		fprintf(stderr,"\t loco_notbest: %d\n",loco_notbest);
		fprintf(stderr,"\t atomshift: %d\n",atomshift);
		fprintf(stderr,"\t rama_th: %f\n",rama_thr);
		fprintf(stderr,"\t rama_thrA: %f\n",rama_thrA);
		fprintf(stderr,"\t rama_thrK: %f\n",rama_thrK);

		fprintf(stderr,"\t kink_merge: %f\n",kink_merge);
		fprintf(stderr,"\t kink_mergeA: %f\n",kink_mergeA);
		fprintf(stderr,"\t nokink_merge: %f\n",nokink_merge);


		fprintf(stderr,"---------------------------------\n");
	}


	//***********************************************************
	// LOAD PROBABILITY DENSITY MAP: Dunbrack(TCB)
	// Load Dunbrack's neighbor-dependent PDFs for Ramachandran-based potentials MAP (Dunbrack_TCB)
	if(load_ramap) {
		pdfs = read_dunbrack(ramafile,&size,&step); // Reads binary
		if(pid==0)
			printf("rcd> reading %s bins %d step %f\n", ramafile, size, step );
	}
	float *****allmaps; // Complete Ramachandran potential (20x20x20 = 8000 72x72 maps)
	if(emodel == 6 || emodel == 56) // Rama potentials
	{
		// Rama stuff for protein modeling (many single runs)
		if(pid==0) {
			fprintf(stdout,"rcd> Generating Ramachandran potential model: %d (20x20x20 = 8000 72x72 maps)... ", rama_model);
			fflush(stdout);
		}
		allmaps = gen_allrama(pdfs,rama_model);
		if(pid==0) {
			fprintf(stdout,"Done!\n");
			fflush(stdout);
		}
	}

	// LOAD H3 PDF maps

	float ***h3maps;
	int h3size;
	float h3step;
	if (h3_switch)  // Set "true" to activate H3 Rama potential
	{
		// open file
		h3maps =  read_h3(h3file,&h3size,&h3step);
		if(pid==0) 
			printf("rcd> reading %s bins %d step %f\n", h3file, size, step );

	}




	int cycles0 = cycles; // required by "best energy" stuff
	int max_path0 = max_path; // required

	//CREATE FOLDER AND SYSTEM COMMANDS
	char    tmp_path[LENGTH_LINE];
	tmp_path[0]='\0';
	char    res_path[LENGTH_LINE];
	res_path[0]='\0';
	char    current_pdb_name[LENGTH_LINE];
	current_pdb_name[0]='\0';
	time_t rawtime_folder;				/*Retrieve current time*/
	struct tm * timeinfo_f_folder;

	if(!folder_name)
	{
		time ( &rawtime_folder );
		timeinfo_f_folder = localtime ( &rawtime_folder );
		strftime(folder_later,80,"run__%d_%b_%Y__%Hh_%M_%S",timeinfo_f_folder);	//:INFO: 	%S (seconds) field to prevent bugs of overwriting folders
	}

	if (switch_file == 1)
		mkdir(folder_later, 0777);

	//INPUT DATA RETRIEVAL
	marker_input	= 	1;


	// CALCULATE SAMPLE NUMBER FOR TRAJECTORY
	if (cycles < number_path)
		number_path = 0;								/*Improbable case but can occur-> set at first solution*/

	// Reverse
	reverse_input = reverse;

	rg = new CRandomMersenne( seed ); // this must be after MPI Init...

	// OPEN OVERALL RESULTS TXT-FILE for writing and creates header
	FILE *summary_file;										/*Create report file if wanted								*/
	if (pid==0 && switch_file==1)											/*File containing results over all different loops			*/
	{
		strcpy(res_path,folder_later);
		strcat(res_path,"/results.txt");

		if( !(summary_file = fopen(res_path,"w")) )
		{
			fprintf(stderr,"Sorry, unable to write file: %s --> Forcing exit!!!\n",res_path);
			exit(1);
		}

		fprintf(summary_file,"%-12s %-6s %-5s %-6s %6s %6s %6s %6s %6s %6s %5s %7s %3s %7s %7s %5s %3s %3s %-8s\n","PDB","RMSD","Sigma","RMSD","<0.5A","<1.0A","<1.5A","<2.0A","<2.5A","<3.0A","Iter","Time","Con","Trial","Success","RMSD","Int","Ext","Sequence(+Anchors)");
		fprintf(summary_file,"%-12s %6s %5s %6s %6s %6s %6s %6s %6s %6s %5s %7s %3s %7s %7s %5s %3s %3s %-8s\n\n","ID","Avg[A]","+-[A]","Min[A]","[%]","[%]","[%]","[%]","[%]","[%]","[#]","[s]","[%]","F/B [%]","F/B [%]","Anch","[%]","[%]"," ");
		fclose(summary_file);
	}

	//SCAN INPUT FILE
	typedef struct
	{									/*Struct to save data for loop closure each PDB entry				*/
		char string[LENGTH_FILE];				/*string for reading in the filename								*/
		int start_number;				/*the start residue number of the loop								*/
		int end_number;				/*the end residue number of the loop								*/
		char chainnr;				/*char array for containing the chain number (A,B,C,D...  )			*/				//:TODO: use strings instead?
		char sequence[LENGTH_SEQ];				/*char array for containing the sequence in single letter code		*/
	} filedata;

	//-------------- READING INPUT TEXT FILE FOR MULTI-LOOP ---------------//
	// (Reads the text file with the input PDB files and loop sequence definition)
	//
	FILE* p_file1;
	int ind 	= 0;					/*count number of pdb files	*/
	char myline[1024];

	p_file1	= fopen(testfile,"r");		/*Read input file*/
	if(p_file1 == NULL)					/*If problem with opening file...*/
	{
		fprintf(stdout,"rcd> Input file %s could not be opened...\nrcd>\n",testfile);
		exit(1);
	}									/*File opening check passed...*/
	while( fgets(myline,1024,p_file1) )
	{
		if(myline[0] != '#')
			ind++;
	}

	fclose(p_file1);

	//DEBUGGING PROTOCOL
	if(pid==0)
		printf("rcd> Reading input file ---> %d pdb cases\n", ind);


	filedata *pdbfile;
	if( !(pdbfile = (filedata*) malloc( sizeof(filedata)*ind )))	/*initial size*/
	{
		printf("Sorry, unable to allocate memory!!!\n");
		exit(1);
	}
	for(int i=0; i<ind; i++)
		pdbfile[i].sequence[0] = '\0'; // Initialize sequence

	p_file1	= fopen(testfile,"r");		/*Read input file*/
	if(p_file1 == NULL)					/*If problem with opening file...*/
	{
		fprintf(stdout,"rcd> Input file %s could not be opened...\nrcd>\n",testfile);
		exit(1);
	}									/*File opening check passed...*/



	// Mon: Simplify this!!!
	ind 	= 0;					/*count number of pdb files																*/
	char mystring[2];

	int nscan = 0;
	while( fgets(myline,1024,p_file1) )
	{
		//printf("%s\n",myline);
		if(myline[0] != '#')
		{
			nscan = sscanf(myline,"%s %d %d %s %s\n",pdbfile[ind].string,&pdbfile[ind].start_number,&pdbfile[ind].end_number,
					mystring, pdbfile[ind].sequence);
			// fprintf(stderr,"rcd> nscan=%d \n",nscan);
			if(nscan == 4 )
			{
				nscan = sscanf(myline,"%s %d %d %s\n",pdbfile[ind].string,&pdbfile[ind].start_number,&pdbfile[ind].end_number,mystring);
				pdbfile[ind].sequence[0] = '\0';
			}
			if( !(nscan == 4 || nscan == 5) )
			{
				fprintf(stderr,"rcd> Please, check input text file... nscan=%d \nForcing exit!\n",nscan);
				exit(1);
			}
			pdbfile[ind].chainnr = mystring[0]; // with %c this does not work...
			if(pid==0)
				printf("rcd> %2d %8s %8d %8d %c %s (%ld)\n",ind+1, pdbfile[ind].string,pdbfile[ind].start_number,
						pdbfile[ind].end_number, pdbfile[ind].chainnr,pdbfile[ind].sequence, (long int) strlen(pdbfile[ind].sequence));
			ind++;
		}


	}

	fclose(p_file1);
	// Finished reading input text file...

	// ALLOCATION OF GLOBAL OBJECTS

	// Define "cycles" for each worker...
#ifdef MPI_RELEASE
	if(pid != 0)
	{
		cycles = cycles / np;
		// loco_best = loco_best / np;
		if(!server)
			printf( "Worker %d number of loops: %d (loco_best= %d)\n", pid, cycles, loco_best);
	}
	else
	{
		printf( "Master (%d) number of loops: %d (loco_best= %d)\n", pid, cycles, loco_best);
		//		if(presampling_switch && presampling >= cycles/np)
		if(presampling_switch && presampling > cycles)
		{
			//			fprintf(stderr,"ERROR: The number of loops per process (%d) is lower or equal than presampling steps (%d) ... Forcing exit!\n",cycles/np,presampling);
			fprintf(stderr,"ERROR: The total number of loops (%d) is lower than presampling steps (%d) ... Forcing exit!\n",cycles,presampling);
			exit(1);
		}
	}
#endif


	//	TIMER OBJECT
	Timer *timer1 = new Timer();	//NEW ---- Maybe BUG PROBLEM

	//	STAT OBJECTS
	Stat *ftimer 		= new Stat(cycles);    //  Stat object for fast timings
	Stat *fcounter 		= new Stat(cycles);    //  Stat object for number of iterations/loop closure
	Stat *frmsd;    //  Stat object for RMSD all obtained loops versus native loop
	if(loco_best > 0)
		frmsd = new Stat(loco_best+loco_notbest);
	else
		frmsd = new Stat(cycles);
	Stat *energies 		= new Stat(cycles);    //  Stat object for Energies
	Stat *energies_intra= new Stat(cycles);    //  Stat object for Intra-Energies
	Stat *energies_extra= new Stat(cycles);    //  Stat object for Extra-Energies
	Stat *energies_presampling; //  Stat object for the PreSampling Energies
	if(presampling_switch)
		energies_presampling = new Stat(presampling);
	Stat *scores = new Stat(cycles);    //  Stat object for Scores (ods LOCO-ICOS energies), it may be Energy or Cross-correlation
	Stat *ccorrs;    //  Stat object for Cross-Correlations
	if(em_switch && (emodel == 1 || emodel == 5 || emodel == 56) ) // Only save ccorrs if Energies are also requested
		ccorrs = new Stat(cycles);    //  Stat object for Cross-Correlations
	Stat *dih_found 	= new Stat(cycles);    //  Stat object for dihedral found of loop
	Stat *rmsdprofile; //  Stat object to save a profile of rmsd evolution on one run
	Stat *counterprofile;//  Stat object to save a profile of counter perturbation related to rmsd evolution
	Stat *tot_rmsd_avg	= new Stat(ind); // Array to store averages of all PDB's
	Stat *tot_rmsd_sig	= new Stat(ind); // Array to store sigmas of all PDB's
	Stat *tot_rmsd_min	= new Stat(ind); // Array to store minima of all PDB's
	Stat *tot_rmsd_05	= new Stat(ind); // Array to store <0.5A of all PDB's
	Stat *tot_rmsd_10	= new Stat(ind); // Array to store <1.0A of all PDB's
	Stat *tot_rmsd_15	= new Stat(ind); // Array to store <1.5A of all PDB's
	Stat *tot_rmsd_20	= new Stat(ind); // Array to store <2.0A of all PDB's
	Stat *tot_rmsd_25	= new Stat(ind); // Array to store <2.5A of all PDB's
	Stat *tot_rmsd_30	= new Stat(ind); // Array to store <3.0A of all PDB's
	Stat *tot_counter	= new Stat(ind); // Array to store total counters of all PDB's
	Stat *tot_time	= new Stat(ind); // Array to store total times of all PDB's
	Stat *tot_con	= new Stat(ind); // Array to store Convergence rates of all PDB's
	Stat *tot_trialF	= new Stat(ind); // Array to store Trial Forward rates of all PDB's
	Stat *tot_trialB	= new Stat(ind); // Array to store Trial Backward rates of all PDB's
	Stat *tot_sucF	= new Stat(ind); // Array to store Success Forward rates of all PDB's
	Stat *tot_sucB	= new Stat(ind); // Array to store Success Backward rates of all PDB's
	Stat *anchor_rmsd = new Stat(ind); // Array to store the final anchor rmsd
	Stat *intra_clashes = new Stat(ind); // Array to store the number of intra-clashes
	Stat *extra_clashes = new Stat(ind); // Array with the number of times there is an extra-clash (clashes with map)

	//	MATRIX OBJECTS
	Matrix *matrix2 	= new Matrix(3,3);								// 	rotation matrix used
	Matrix *matrix  	= new Matrix(3,3);								// 	coordinates of the loop end atoms
	Matrix *polpoint	= new Matrix(3,3);								// 	coordinates of the fixed anchor to go
	Matrix *anchorpoint	= new Matrix(3,3);								//  rotated loop end coordinates
	Matrix *result1 	= new Matrix(3,3);
	Matrix *result2 	= new Matrix(3,3);								//	rotated anchor coordinates
	Matrix *result3 	= new Matrix(3,3);								//	rotation matrix in rotate
	Matrix *anchor_end 	= new Matrix(3,3);								//	Active end anchor loop coordinates
	Matrix *anchor_end_i 	= new Matrix(3,3);								//	Active end anchor loop coordinates forward direction
	Matrix *anchor_end_r 	= new Matrix(3,3);								//	Active end anchor loop coordinates backward direction

	// ALEX: RMSD GENERATION (BEGIN) VARIABLES AND COUNTERS
	float rmsd_pdbs_gen = 0.0, rmsd_pdbs_sol = 0.0;
	float rmsd_pdbs_avg_gen = 0.0, rmsd_pdbs_avg_sol = 0.0;
	int total_gen = 0, total_sols = 0, counter_no_filter_tot = 0;
	int self_clash_count = 0;
	int self_clash_count_failed = 0;
	int n_almost_closed = 0;
	int n_LoopClosure = 0;
	int n_bad_closures = 0;
	int n_good_closures = 0;

	// RASP's variables
	//	 rmodel rasp_in,rasp_out; // RASP's atomic model of protein
	string rasp_mut; // repacking mask (the sequence in 1-char per aminoacid format)

	// ICOS initialization stuff...
	double icosVertices[NVERT][3]; // Icosaedron vertices (size = 1.0)
	double ****loco; // Loco energy data: loco[NAAS][NAAS][NDIS][NTRI]
	if(emodel == 1) // ICOSA
	{
		// Initialize icosahedron data
		init_icos(icosVertices); // Load an icosaedron of size = 1.0

		//	Allocate dynamic array of dimensions: loco[NAAS][NAAS][NDIS][NTRI]
		loco = (double ****) malloc( sizeof(double ***) * NAASloco );
		for(int i=0; i<NAASloco; i++)
		{
			loco[i] = (double ***) malloc( sizeof(double **) * NAASloco );
			for(int j=0; j<NAASloco; j++)
			{
				loco[i][j] = (double **) malloc( sizeof(double *) * NDIS );
				for(int k=0; k<NDIS; k++)
					loco[i][j][k] = (double *) malloc( sizeof(double) * NTRI );
			}
		}
		read_loco(loco, file_energy); // Load loco energy
	}

	// KORP initialization stuff...
	korp *korpmap; // KORP's map structure with all stuff required...
	frame *frames; // Environment frames
	frame *framesloop; // Loop frames
	contact *contacts = NULL; // Contacts array
	int ncont; // Number of contacts
	int *resnums = NULL; // Environment residue numbers (loop-less)
	int *resnumsloop = NULL; // Loop residue numbers (only mobile)
	char *reschainids = NULL; // Environment chain-ids (loop-less)
	char *reschainidsloop = NULL; // Loop chain-ids (only mobile)
	if( emodel == 5 || emodel == 56) // KORP or Rama+KORP
	{
		korpmap = readKORP(file_energy); // Wrapper to read a whole KORP energy map (3/4/6D)

		// Overwrite map's bonding_factor with user-defined bonding_factor if parser's bonding_factor has been provided
		if(bonding_factor >= 0.0)
		{
			fprintf(stderr,"rcd> Overwrite map's bonding_factor (%f) with user's: %f\n",korpmap->bonding_factor,bonding_factor);
			korpmap->bonding_factor = bonding_factor; // Use parser's bonding_factor
		}
	}

	//------------------------------------------------------------------------------
	//
	//                       LOOP FOR ALL THE PDBS
	//
	//------------------------------------------------------------------------------
	//
	// START PDB ITERATIONS
	for (int i_pdb = 0; i_pdb < ind; i_pdb++)
	{
		// ALEX : This 4 variables allow a better control of threshold and filters applied to each pdb file.
		ramachandran_inside_aux = ramachandran_inside;
		ramachandran_check_aux = ramachandran_check;

		// Counters for final statistics...
		int n_rama_fail = 0; // Number of conformations with any bad dihedral angle
		int avg_rama_fail = 0; // Average number of aminoacids with bad dihedral angles
		int n_infinite =0;
		int n_max_path = 0;
		int n_bondselect_fail = 0;
		int n_max_rmsd_count = 0;
		int n_max_grid_clash = 0;
		int n_rand_repairs = 0;
		int for_phi_counter = 0;
		int for_psi_counter = 0;
		int back_phi_counter = 0;
		int back_psi_counter = 0;
		// INITIALIZING some required Pieter's variables
		//
		int start   =	pdbfile[i_pdb].start_number-1;
		int end     =	pdbfile[i_pdb].end_number+1;
		//	start and end residue position of loop-anchors in protein main chain
		//		- start residue anchor starts at nitrogen; anchor = N, CA & CO
		//		- end residue anchor ends at CO; anchor = N, CA & CO
		//	mind implementation single point

		float* co1;	//  pointer co1 to work in function readinPDB & nr_atoms = number main chain
		int nr_atoms;	//  atoms in protein chain
		int startres;	//  start position of first amino acid residue in PDB-file: at first initialized at 1
		//  but determined from reading in PDB-file
		double 	rmsd;	//  Computed RMSD value each iteration/loop closure*/
		long int counter_cycles_total;	//counter number executed cycles in total
		long int counter_cycles_single;	//counter number executed cycles per loop closure -locally
		long int total_counter;		//counter number executed cycles per loop closure in total
		long int counter_path;		// counter loop closure trajectory
		long int counter_failed;		// counter failed loop closures

		int ibond;       		// rotation bond in active computation following order actual bonds
		int ibond_r;		// rotation bond in randomization following order actual bonds
		double 	angle_r;	// dihedral angle value in randomization procedure
		double omega_r;
		int bond_previous;	// Track of the bond that was previously rotated: previously rotated bond in algorithm itself
		int size_path;		// Size of the array to save the loop closuretrajectory
		int range_mod;		// number of rotatable bond modified for variable choice of DoF

		//	Loop variables
		int nr_loop;		// numbers of residues of protein loop without end anchors
		int last;			// index last atom of the protein loop itself (first atom of Ct anchor)
		int	range;		// number of rotatable bonds from which one can choose
		int	nr_atoms_loop;	// total number of atoms in the loop
		int	start_co_prot;	// start index loop in backbone complete protein

		ibond 		        	= 2;
		range_mod 				= 5;
		counter_cycles_total 	= 0;
		counter_cycles_single 	= 0;
		total_counter			= 0;
		counter_path			= 0;
		counter_failed			= 0;

		nr_loop		= end - start - 1; // number of intermediate residues (without anchors)

		// fprintf(stderr,"max_path= %d\n",max_path);
		if(save_movie)
		{
			rmsdprofile 	= new Stat(max_path+1);//  Stat object to save a profile of rmsd evolution on one run
			counterprofile 	= new Stat(max_path+1);//  Stat object to save a profile of counter perturbation related to rmsd evolution
		}

		//************************************************
		int *iaa=NULL; // AA index for each loop aminoacid
		float **accs=NULL; // Accumulated profile for each loop aminoacid
		float ***maps=NULL; // PDF-maps for each loop aminoacid
		float ***ramaps=NULL; // Ramachandran PDFs for each loop aminoacid (Rama energy model)
		float ***binmaps=NULL;
		float ***accs_phi=NULL; // 1D probabilities for fixed Phi
		float ***accs_psi=NULL; // 1D probabilities for fixed Psi
		//		float ***binmaps_post=NULL;
		float *rama_cutoff=NULL; // Cutoff values for selected Ramachandran threshold

		char **secu;
		if( !( secu = (char **) malloc(sizeof(char *) * (nr_loop+2)) ) ) // Mon: free memory later...
		{
			printf("Sorry, unable to allocate secu memory!!!\n");
			exit(1);
		}

		int *seq_num;
		seq_num=(int *) malloc(sizeof(int) * (nr_loop+2));

		max_path = max_path0 + (int) path_slope * (nr_loop - 4);
		if(save_movie)
			if (max_path > 8000) max_path=8000;
		if(pid==0)
			printf("rcd> Current max_path= %d\n", max_path);

		nr_atoms_loop 	= 3 * (nr_loop + 2);
		size_path		= max_path*nr_atoms_loop; /*Size set to allow for maximum number of iterations*/
		last		= 3 * nr_loop + 3;
		range		= 2 * nr_loop + 1;
		start_co_prot 	= 3 * (start - 1);
		startres 		= 1;

		//CHANGE INTEGER TO CORRECT SIZE FOR ALGORITHM
		if (switch_dof == 1)
			range_mod	= range + 1;
		else if (switch_dof == 0)
			range_mod  	= range;

		//INTERNAL COORDINATE ARRAYS
		double *bond_length 	; // = new double[nr_atoms_loop]; //	Array bond lengths (Angstrom) for internal coordinate set (ic set)
		double *valence_angle   ; //= new double[nr_atoms_loop]; //	Array valence angles (rad) for internal coordinate set (ic set)
		double *dihedral_angle  ; //= new double[nr_atoms_loop]; //	Array dihedral angles for internal coordinate set (ic set)

		bond_length = (double *) malloc( sizeof(double) * nr_atoms_loop);
		valence_angle = (double *) malloc( sizeof(double) * nr_atoms_loop);
		dihedral_angle = (double *) malloc( sizeof(double) * nr_atoms_loop);


		double *bond_length_back ; //	= new double[nr_atoms_loop]; //	Array bond lengths (Angstrom) for internal coordinate set (ic set)
		double *valence_angle_back; //  = new double[nr_atoms_loop]; //	Array valence angles (rad) for internal coordinate set (ic set)

		bond_length_back = (double *) malloc( sizeof(double) * nr_atoms_loop);
		valence_angle_back = (double *) malloc( sizeof(double) * nr_atoms_loop);

		double *bond_length2 	 ;// = new double[nr_atoms_loop]; //	Array bond lengths (Angstrom) for ic set no_loop case
		double *valence_angle2   ; //= new double[nr_atoms_loop]; //	Array valence angles (rad) for ic set no_loop case
		double *dihedral_angle2  ;//= new double[nr_atoms_loop]; //	Array dihedral angles for ic set no_loop case

		bond_length2 = (double *) malloc( sizeof(double) * nr_atoms_loop);
		valence_angle2 = (double *) malloc( sizeof(double) * nr_atoms_loop);
		dihedral_angle2 = (double *) malloc( sizeof(double) * nr_atoms_loop);


		double *dihedral_angle_i   ; // = new double[nr_atoms_loop]; //Initialization array dihedral angles (ic set) in forward direction
		double *dihedral_angle_aux ; // = new double[nr_atoms_loop]; //Auxiliary array for dihedral angles
		double *dihedral_angle_r   ; // = new double[nr_atoms_loop]; //Initialization array dihedral angles (ic set) in backward direction
		double *dihedral_angle_f   ; //  = new double[nr_atoms_loop]; //Initialization array dihedral angles (ic set) in forward direction

		dihedral_angle_i = (double *) malloc( sizeof(double) * nr_atoms_loop);
		dihedral_angle_aux = (double *) malloc( sizeof(double) * nr_atoms_loop);
		dihedral_angle_r = (double *) malloc( sizeof(double) * nr_atoms_loop);
		dihedral_angle_f = (double *) malloc( sizeof(double) * nr_atoms_loop);




		//		 int    *hybridisation		= new int[nr_atoms_loop];    //      switch for tree algorithm
		//
		//		for (int i = 0; i < nr_atoms_loop/3; i++)   /* hybridisation marker index crucial */
		//		{
		//			hybridisation[3*i]   	= 1;			/*sp2 at N atom in loop itself*/
		//			hybridisation[3*i+1] 	= 2;			/*sp3 at CA atom*/
		//			hybridisation[3*i+2] 	= 1;			/*sp2 at CO atom*/
		//		}
		//

		for (int i = 0; i <  nr_atoms_loop; i++)
		{
			bond_length[i] 	= 0.0;/*initialize ic arrays*/
			bond_length2[i]	= 0.0;
			bond_length_back[i] = 0.0;
			valence_angle[i]   	= 0.0;
			valence_angle2[i]   	= 0.0;
			valence_angle_back[i] = 0.0;
			dihedral_angle[i]	= 0.0;/*initialize dihedral changes*/
			dihedral_angle2[i]  	= 0.0;
			dihedral_angle_i[i]  	= 0.0;
			dihedral_angle_r[i]  	= 0.0;
			dihedral_angle_f[i]  	= 0.0;
		}

		if(emodel == 1 || emodel == 5 || emodel == 56) // ICOSA or KORP
			// Allocate linear array of coordinates of the mobile N, CA, C loop atoms.
			loop_coord = (float *) malloc(sizeof(float) * (nr_atoms_loop-6)*3);

		// Mon: Change "rand()" by "rg->Random()" in GA libraries!!!
		srand( time(NULL) );	/*seed for random number generation*/

		if(pid==0)
		{
			printf("rcd>\nrcd>________________________________________________________________________________\n");
			printf("rcd>\nrcd> Processing %s                                                      %d/%d \nrcd>\n", pdbfile[i_pdb].string, i_pdb+1, ind );
		}
		int lon=0;
		current_pdb_name[0] = ' ';
		lon=strlen(pdbfile[i_pdb].string)-4;
		strncpy(current_pdb_name,pdbfile[i_pdb].string, lon);
		current_pdb_name[lon] = '\0';
		// END INITIALIZATION of Pieter's variables


		//
		// READING PDB COORDINATES
		//

		// READING PDB COORDINATES
		Macromolecule *protein = new Macromolecule("pdb");
		Macromolecule *mol = new Macromolecule("pdb1");	//  macromolecule object to obtain protein coordinates

		char *seq;
		protein->readPDB(pdbfile[i_pdb].string);
		//protein->info(stdout);
		if(delHydrogens_switch)
		{
			if(pid==0) printf( "rcd> Deleting Hydrogen atoms (if any)...\n" );
			protein->deleteHYDS();
		}
		if(delHeteros_switch) // CA model can't deal with HETATM's... (TO DO)
		{
			if(pid==0) printf( "rcd> Deleting Hetero-atoms (if any)...\n" );
			protein->delete_heteros();
		}
		if(delWaters_switch)
		{
			if(pid==0) printf( "rcd> Deleting Water molecules (if any)...\n" );
			protein->delete_waters();
		}
		protein->format_residues(false,-1); // Mon 4/11/2016: Not-sorted atoms bug fixed!
		// protein->writePDB("formated.pdb");

		// Select chain
		Chain *ch;
		pdbIter *iter_ch;
		int sel_Ch=-1;
		char in_chain;
		in_chain=pdbfile[i_pdb].chainnr;
		iter_ch = new pdbIter( protein ); // iters current protein

		// Getting chain number
		//		fprintf(stderr,"Requested chain: %c\n",in_chain);
		if(!(in_chain == '-' || in_chain == '*')) // By default, an '*' or '-' must be provided ALWAYS in the chain field!
			for( iter_ch->pos_chain = 0; !iter_ch->gend_chain(); iter_ch->next_chain() )
			{
				ch = iter_ch->get_chain();
				if (in_chain==ch->getName()[0])
					sel_Ch=iter_ch->pos_chain ;
			}
		if (sel_Ch==-1)
		{
			iter_ch->pos_chain = 0;
			ch = iter_ch->get_chain();
			in_chain = ch->getName()[0];
			if(pid==0)
				printf("rcd> Select first chain: %c\n",in_chain);
		}
		else
			if(pid==0)
				printf("rcd> Select chain %c\n", in_chain);

		// KORP, initializing arrays resnumsloop and reschainidsloop (residue numbers and chain ids for current loop)
		if(emodel == 5 || emodel == 56)
		{
			resnumsloop = (int *) malloc( sizeof(int) * nr_loop );
			reschainidsloop = (char *) malloc( sizeof(char) * nr_loop );
			for(int i=0; i<nr_loop; i++) // "start" is now the Nt-anchor index (thus +1)
			{
				resnumsloop[i] = start + i + 1;
				reschainidsloop[i] = in_chain; // Loop chain-ids (only mobile)
			}
		}
		// Get complete loop (N,CA,C,O)
		Conditions *conds = new Conditions();
		Condition *cond = new Condition(-1,-1,-1,-1,sel_Ch,sel_Ch,start,end,-1,-1);
		cond->add(" N  ");
		cond->add(" CA ");
		cond->add(" C  ");
		if(rmsd_O) // O-atom is considered
		{
			cond->add(" O  ");
			//			cond->add(" CB "); // Mon: Ideally, CB must be considered just for self-consistency in checking...
			// WARNING: if CB is added, the number of atoms present must be checked and reconsidered...
		}
		conds->add(cond);
		Macromolecule *loop = protein->select(conds);

		// Mon:
		if( (!rmsd_O && loop->num_atoms() == (nr_loop+2)*3) || (rmsd_O && loop->num_atoms() == (nr_loop+2)*4) )
		{
			present_seq = 1;
			if(pid==0)
				fprintf(stderr,"rcd> Loop present in PDB.\n");
		}
		else
		{
			present_seq = 0;
			if(pid==0)
			{
				if(!rmsd_O)
					fprintf(stderr,"rcd> Loop missing in PDB. (%d atoms expected but %d found)\n",(nr_loop+2)*3,loop->num_atoms());
				else
					fprintf(stderr,"rcd> Loop missing in PDB. (%d atoms expected but %d found)\n",(nr_loop+2)*4,loop->num_atoms());
			}

			if(bench_rmsd)
			{
				bench_rmsd = false;
				fprintf(stderr,"rcd> Disabling --bench option.\n");
			}
		}

		// BUG fixed! Mon: What the chapa!
		if( !((native == 1 || present_seq == 1) || bench_rmsd) ) // If loop not present in PDB
		{
			char *myseq = loop->get_sequence(); // 1-letter code for anchors
			seq = (char *) malloc(sizeof(char)*nr_atoms_loop/3+1); // +1 for termination char '\0'
			seq[0] = myseq[0]; // Anchor N-t
			seq[nr_atoms_loop/3-1] = myseq[1]; // Anchor C-t
			for(int x=1; x<(nr_atoms_loop/3-1); x++)
				seq[x] = pdbfile[i_pdb].sequence[x-1]; // Mobile residues
			seq[nr_atoms_loop/3] = '\0'; // termination char '\0'
			if(pid==0)
				fprintf(stderr,"rcd> Loop not present in PDB - Get sequence from input file.\n");
			mutate = true;
		}
		else // If loop present in PDB
		{
			seq=loop->get_sequence(); // 1-letter code
			if(pid==0)
				fprintf(stderr,"rcd> Sequence from PDB: seq= %s\n",seq);
			if(strlen(pdbfile[i_pdb].sequence) == 0) // If sequence was not inserted in Input file, then copy it from PDB
			{
				mutate = false; // sets requested sequence in <name>_loop1.pdb
				for(int x=1; x<(nr_atoms_loop/3-1); x++)
					pdbfile[i_pdb].sequence[x-1] = seq[x]; // Mobile residues
				if(pid==0)
				{
					fprintf(stderr,"rcd> Checking... pdbfile[i_pdb].sequence = %s\n",pdbfile[i_pdb].sequence);
					fprintf(stderr,"rcd> Using sequence from PDB: seq= %s\n",seq);
				}
			}
			else
			{
				mutate = true; // sets requested sequence in <name>_loop1.pdb
				if(pid==0)
					fprintf(stderr,"rcd> Loop present in PDB but using custom sequence from Input: %s\n",pdbfile[i_pdb].sequence);
				for(int x=1; x<(nr_atoms_loop/3-1); x++)
					seq[x] = pdbfile[i_pdb].sequence[x-1]; // Mobile residues
			}
		}

		if(pid==0)
			fprintf(stderr,"rcd> Loop %d-%d including anchors (sequence: %s, %d loops requested)\n",start,end,seq, cycles);
		// fprintf(stderr,"rcd> %d %d\n",loop->get_num_atoms(),3*(end-start+1));

		int *iseql=NULL; // mobile loop sequence in integer indices
		int *iseqe=NULL; // non-mobile (environment) sequence in integer indices
		float *coorde; // array with the coordinates of the non-mobile region (environment)
		int nseqe; // Number of environment residues
		int anchorNt,anchorCt;
		if(emodel == 1 || emodel == 5 || emodel == 56) // ICOSA or KORP per-PDB stuff...
		{
			// Get non-mobile region (environment)
			// Main PDB selection of N,CA,C atoms (and, optionally, non-loop residues)
			Condition *myncac;
			Conditions *myncac2= new Conditions();
			myncac = new Condition(-1,-1,-1,-1,-1,-1,start+1,end-1,-1,-1); // Excluding mobile region
			myncac2->add(myncac);
			if(pid==0)
				printf( "rcd> Selecting non-mobile residues (ichain= %d), i.e. excluding residues from %d to %d\n", sel_Ch, start+1, end-1 );
			Macromolecule *mole = protein->select_cpy(myncac2,false); // inverse selection
			Conditions *conds= new Conditions();
			Condition *cond= new Condition(-1,-1,-1,-1,-1,-1,-1,-1,-1,-1);
			cond->add(" N  ");
			cond->add(" CA ");
			cond->add(" C  ");
			conds->add(cond);
			Macromolecule *mole2 = mole->select_cpy(conds);
			delete mole;

			if(emodel == 5 || emodel == 56) // KORP
			{
				resnums = getResNums(mole2);
				reschainids = getResChainIds(mole2);
			}

			char *seqe = mole2->get_sequence(); // 1-letter code sequence of Environment (without mobile loop)
			nseqe = mole2->get_num_fragments(); // number of residues
			seq2iseq(seqe,&iseqe,nseqe);
			free(seqe);
			mole2->coordMatrix(&coorde); // Get environment coordinates

			if(pid==0)
				fprintf(stderr, "rcd> loop %d aminoacids sequence: %s\n", (nr_atoms_loop/3-2), seq );
			seq2iseq(seq+1,&iseql,nr_atoms_loop/3-2); // Get loop sequence indices array

			// Get residue indices from main pdb
			Chain *ch;
			pdbIter *iter_chain = new pdbIter(mole2);
			pdbIter *iter_frag;
			char molchain;
			Fragment *frag;
			int icont=0; // index counter
			bool exitnow = true;
			for( iter_chain->pos_chain = 0; !iter_chain->gend_chain() && exitnow; iter_chain->next_chain() ) // screen chains
			{
				ch=iter_chain->get_chain();
				molchain = ch->getName()[0];
				iter_frag = new pdbIter(ch);
				if(molchain == in_chain)
				{
					for( iter_frag->pos_fragment = 0; !iter_frag->gend_fragment() && exitnow; iter_frag->next_fragment() ) // screen residues
					{
						frag = iter_frag->get_fragment();
						if(start == frag->getIdNumber() ) // ri-1 --> Nt anchor (not mobile)
						{
							anchorNt = icont;
							if(pid==0)
								fprintf(stdout,"rcd> Nt-1 anchor residue %d (%c chain) found in main PDB: %d (internal index)\n",start,in_chain,anchorNt);
						}
						if(end == frag->getIdNumber() ) // rf+1 --> Ct anchor (not mobile)
						{
							anchorCt = icont;
							if(pid==0)
								fprintf(stdout,"rcd> Ct+1 anchor residue %d (%c chain) found in main PDB: %d (internal index)\n",end,in_chain,anchorCt);
							exitnow = false;
						}
						icont++;
					}
				}
				else
					icont += iter_frag->num_fragment();
				delete iter_frag;
			}
			delete iter_chain;
			delete mole2; // this copy is not needed anymore
		}

		// Get complete loop2 (N,CA,C,O,CB) // Mon: Chapa for Native Bump energy computation, do not remove...
		Conditions *conds2= new Conditions();
		cond->add(" CB ");
		conds2->add(cond);
		Macromolecule *loop2 = protein->select(conds2);

		if(rmsd_O) // O-atom is considered
		{
			if(	loop->get_num_atoms() != 4*(end-start+1)) // if all native-loop backbone coordinates are NOT present in PDB...
			{
				if(native == 1 || present_seq == 1)
				{
					if(pid==0)
						printf("rcd> Error: Missing loop backbone (N,CA,C,O) atoms, only %d found but %d expected!\n", loop->get_num_atoms(), 4*(end-start+1) );
					//exit(1);
				}
				if(bench_rmsd)
				{
					if(pid==0)
						printf("rcd> Warning: Missing loop backbone (N,CA,C,O) atoms, only %d found but %d expected! Disabling --bench options.\n", loop->get_num_atoms(), 4*(end-start+1) );
					bench_rmsd = false;
				}
			}
		}
		else
		{
			if(	loop->get_num_atoms() != 3*(end-start+1)) // if all native-loop backbone coordinates are NOT present in PDB...
			{
				if(native == 1 || present_seq == 1)
				{
					if(pid==0)
						printf("rcd> Error: Missing loop backbone (N,CA,C) atoms, only %d found but %d expected!\n", loop->get_num_atoms(), 3*(end-start+1) );
				}
				if(bench_rmsd)
				{
					if(pid==0)
						printf("rcd> Warning: Missing loop backbone (N,CA,C) atoms, only %d found but %d expected! Disabling --bench options.\n", loop->get_num_atoms(), 3*(end-start+1) );
					bench_rmsd = false;
				}
			}
		}

		// If input-file sequence and pdb-found sequence do not match, then "bench" stuff is not fair...
		if( bench_rmsd && server && strncmp(seq+1,pdbfile[i_pdb].sequence,nr_loop) != 0)
		{
			if(pid==0)
				printf("rcd> Remark: PDB sequence (%s) and input sequence (%s) mismatch! Disabling --bench options.\n", seq, pdbfile[i_pdb].sequence);
			bench_rmsd = false;
		}

		// Check missing atoms....

		//  N-terminal
		Condition *ncond= new Condition(-1,-1,-1,-1,sel_Ch,sel_Ch,start,start,-1,-1);
		Conditions *nconds= new Conditions();
		char *seq_nt;
		ncond->add(" N  ");
		ncond->add(" CA ");
		ncond->add(" C  ");
		nconds->add(ncond);
		Macromolecule * ntloop=protein->select(nconds);
		if(output_debug)
			ntloop->writePDB("loop_nt.pdb");
		if (ntloop->get_num_atoms()!=3)
		{
			if(pid==0)
				printf("rcd> Error: Missing coordinates at Nt anchor\n");
			exit(1);
		}
		seq_nt=ntloop->get_sequence();

		// C-terminal
		Condition *tcond= new Condition(-1,-1,-1,-1,sel_Ch,sel_Ch,end,end,-1,-1);
		Conditions *tconds= new Conditions();
		char *seq_ct;
		tcond->add(" N  ");
		tcond->add(" CA ");
		tcond->add(" C  ");
		tconds->add(tcond);
		Macromolecule * ctloop=protein->select(tconds);
		if(output_debug)
			ctloop->writePDB("loop_ct.pdb");
		if (ctloop->get_num_atoms()!=3)
		{
			if(pid==0)
				printf("rcd> Error: Missing coordinates at Ct anchor\n");
			exit(1);
		}
		seq_ct=ctloop->get_sequence();
		no_loop=0;

		// Check missing atoms....
		//  O N-terminal - - 1
		pdbIter *iter2;
		Condition *ncondOX= new Condition(-1,-1,-1,-1,sel_Ch,sel_Ch,start,start,-1,-1);
		Conditions *ncondsOX= new Conditions();
		ncondOX->add(" O  ");
		ncondsOX->add(ncondOX);
		Macromolecule * ntOX=protein->select(ncondsOX);
		if(output_debug)
			ntOX->writePDB("loop_O.pdb");
		if (pid==0 && ntOX->get_num_atoms()!=1)
		{
			fprintf(stderr,"rcd> Warning: Missing O coordinates at Nt anchor\n");
			if(no_ntpsi) // If Nt-Oxigen optimization is disabled
			{
				fprintf(stderr,"rcd> Error: O atom required for current options. Forcing exit!\n");
				exit(1);
			}
		}

		// Mon (13/7/2015) Get Nt's Oxygen position
		Tcoor NtOf;
		if(no_ntpsi) // If Nt-Oxygen optimization is disabled
		{
			pdbIter *iter4 = new pdbIter( ntOX );
			iter4->pos_atom = 0;
			(iter4->get_atom())->getPosition(NtOf);
			delete iter4;
		}

		//************************************************************
		// ALEX: Get AA-1 AA+1 next to C-term Anchor and N-term Anchor
		char *seq_nt1;
		char *seq_ct1;
		// REVERSE: N-term: Anchor - 1
		Condition *ncond1= new Condition(-1,-1,-1,-1,sel_Ch,sel_Ch,start-1,start-1,-1,-1);
		Conditions *nconds1= new Conditions();
		ncond1->add(" CA "); // two atoms are mandatory, otherwise it crashes...
		ncond1->add(" C  ");
		nconds1->add(ncond1);
		Macromolecule *ntloop1 = protein->select(nconds1);
		if (ntloop1->get_num_atoms() != 2)
		{
			if(pid==0)
				printf("rcd> Error: Missing C or CA atom coordinates at Nt anchor -1 \n");
			exit(1);
		}
		seq_nt1=ntloop1->get_sequence();

		// Nt-1  atom C
		Tcoor NtCf,NtCAf,NtNf;
		char *at_name;
		pdbIter *iter3 = new pdbIter( ntloop1 );
		Atom *atom;
		for ( iter3->pos_atom = 0; !iter3->gend_atom(); iter3->next_atom() )  // screens all-atoms
		{
			atom = ( Atom * ) iter3->get_atom();
			at_name = atom->getName();
			if ( strcmp(at_name," C  ") == 0  )
				atom->getPosition(NtCf);
		}
		delete iter3;

		double NtC[3];
		NtC[0] = (double) NtCf[0];
		NtC[1] = (double) NtCf[1];
		NtC[2] = (double) NtCf[2];
		if(pid==0)
			fprintf(stderr,"rcd> NtC= %f %f %f\n",NtC[0],NtC[1],NtC[2]);

		// FORWARD: C-terminal -> Anchor+1
		Condition *tcond1= new Condition(-1,-1,-1,-1,sel_Ch,sel_Ch,end+1,end+1,-1,-1);
		Conditions *tconds1= new Conditions();
		tcond1->add(" N  ");
		tcond1->add(" CA "); // two atoms are mandatory, otherwise it crashes...
		tconds1->add(tcond1);
		Macromolecule * ctloop1 = protein->select(tconds1);
		if (ctloop1->get_num_atoms() != 2)
		{
			if(pid==0)
				fprintf(stderr,"rcd> Error: Missing N or CA atom coordinates at Ct anchor +1\n");
			exit(1);
		}
		seq_ct1=ctloop1->get_sequence();

		// Ct+1  atom N
		Tcoor CtNf,CtCAf,CtCf;
		iter3 = new pdbIter( ctloop1 );
		for ( iter3->pos_atom = 0; !iter3->gend_atom(); iter3->next_atom() )  // screens all-atoms
		{
			atom = ( Atom * ) iter3->get_atom();
			at_name = atom->getName();
			if ( strcmp(at_name," N  ") == 0  )
				atom->getPosition(CtNf);
			// Load here Ct coordinates into coaux
		}
		delete iter3;

		double CtN[3];
		CtN[0] = (double) CtNf[0];
		CtN[1] = (double) CtNf[1];
		CtN[2] = (double) CtNf[2];
		if(pid==0)
			fprintf(stderr,"rcd> CtN= %f %f %f\n",CtN[0],CtN[1],CtN[2]);

		//************************************************************

		// Load "secu" array (contains sequence)
		if ( present_seq==1 ) // All loop atoms present (get sequence input PDB)
		{
			for (int i = 0; i < nr_loop+2;i++)
			{
				pdbfile[i_pdb].sequence[i]=seq[i];
				for (int j = 0; j < 20;j++)
					if(pdbfile[i_pdb].sequence[i]==AA[j].aa_name1 || pdbfile[i_pdb].sequence[i]==AA[j].aa_name1+32) // check both Uppercase and Lowercase
					{
						secu[i]=AA[j].aa_name3; // 3-letter code
						break;
					}
			}
			pdbfile[i_pdb].sequence[nr_loop+2]='\0';
			if(pid==0)
				printf("rcd> Loop coordinates present %s\n", pdbfile[i_pdb].sequence);

		}
		else // Only anchors present (sequence already provided in input file)
		{
			if(pid==0)
				printf("rcd> Loop coordinates NOT present\n");

			// Mon: Zurrarse esto...
			no_loop = 0;
			Atom *at=NULL;
			Chain *cad=NULL;
			Segment *seg=NULL;
			Fragment *frag=NULL;
			Fragment *fragm=NULL;
			cad = ( Chain * ) ntloop->getCurrent();
			seg = ( Segment * ) cad->getCurrent();
			fragm= ( Fragment * ) seg->getCurrent();

			char letter=' ';
			int cont_res=0,aa;
			int atom_no=3;
			float occ, tfact;
			//fprintf(stderr,"mysequence= %s\n",pdbfile[i_pdb].sequence);
			cont_res=1;

			//			for (int i = 0; i < nr_loop+2;i++) // <--- Watch out this!
			for (int i = 0; i < nr_loop;i++) // <--- Watch out this!
				for (int j = 0; j < 20;j++)
					if(pdbfile[i_pdb].sequence[i]==AA[j].aa_name1 || (j==12 && pdbfile[i_pdb].sequence[i]=='p')) // Proline is j=12. The 'p' required for cis-Proline
					{
						secu[i+1]=AA[j].aa_name3; // Loading "secu" array (sequence)
						// iaa[i+1]=resnum_from_resname
						//printf("--> %d %c %s\n", cont_res + start, AA[j].aa_name1, secu[i+1]);

						// Building backbone
						aa=resnum_from_resname( secu[i+1] );
						frag=new Residue(secu[i+1],start+cont_res, start+cont_res, letter);
						// N
						at = new Atom( atom_types[AA[aa].atom[0].fullatom_type-1], AA[aa].atom[0].atom_name, AA[aa].atom[0].icoor, i*3+1 );
						frag->add(at);
						// CA
						at = new Atom( atom_types[AA[aa].atom[1].fullatom_type-1], AA[aa].atom[1].atom_name, AA[aa].atom[1].icoor, i*3+2 );
						frag->add(at);
						// C
						at = new Atom( atom_types[AA[aa].atom[2].fullatom_type-1], AA[aa].atom[2].atom_name, AA[aa].atom[2].icoor, i*3+3 );
						frag->add(at);
						fragm->add(frag);
						cont_res++;
						break;
					}

			//			ntloop->info(stdout);
			//			ntloop->writePDB("loop_i.pdb");

			// Nt aminoacid id
			for (int j = 0; j < 20;j++)
				if(seq_nt[0]==AA[j].aa_name1)
				{
					secu[0]=AA[j].aa_name3;
					break;
				}

			// Ct aminoacid id
			for (int j = 0; j < 20;j++)
				if(seq_ct[0]==AA[j].aa_name1)
				{
					secu[nr_loop+1]=AA[j].aa_name3;
					break;
				}

			// Add C-term
			pdbIter *iter;
			iter = new pdbIter( ctloop );
			frag=new Residue(secu[nr_loop+1],start+cont_res, start+cont_res, letter);
			for(iter->pos_atom=0;!iter->gend_atom();iter->next_atom())
			{
				at=iter->get_atom();
				frag->add(at);
			}
			fragm->add(frag);
			cont_res++;
			delete iter;

			//ntloop->info(stdout);
			if(output_debug)
				ntloop->writePDB("nloop.pdb");

			char *seqb; // <-- Whatch out this (seqb != seq)
			seq=ntloop->get_sequence(); // "seq" contains loop+anchors sequence, and "pdbfile[i_pdb].sequence" the loop only sequence.

			for (int i = nr_loop; i >= 0; i--)
				pdbfile[i_pdb].sequence[i+1] = pdbfile[i_pdb].sequence[i];
			pdbfile[i_pdb].sequence[0] = seq[0];
			pdbfile[i_pdb].sequence[nr_loop+1] = seq[nr_loop+1];
			pdbfile[i_pdb].sequence[nr_loop+2] = '\0';

			if(pid==0)
			{
				printf("rcd> Rebuilding %d res anchors seq: %s\n", nr_loop, seq);
				printf("rcd> input file sequence: %s\n", pdbfile[i_pdb].sequence);
			}

			if (cont_res-2!=nr_loop)
			{
				if(pid==0)
					printf("rcd> Error: Mismatch loop size (%d) with seq residues (%d)\nrcd>\n",cont_res-2, nr_loop );
				exit(1);
			}
		}
		// END NOT-present coordinates

		//		ntloop->info(stdout);
		//		ntloop->writePDB("loop_i.pdb");
		//		fprintf(stderr,"pdbfile.sequence= %s\n",pdbfile[i_pdb].sequence);
		//		for (int i = 0; i < nr_loop+2;i++) // <--- Watch out this!
		//			fprintf(stderr,"secu[%d]= %s\n",i,secu[i]);


		// *******************************
		// RASP stuff initialization
		// *******************************
		//		int rasp_start,rasp_end;
		//		if(rasp_switch)
		//		{
		//			// macromolecule2rasp(protein,rasp_in,rasp_mut); // Interface between our "Macromolecule" and RASP's "rmodel" (Zhichao Miao)
		//			// read_pdb(pdbfile[i_pdb].string,rasp_in_array[i_pdb],rasp_mut[i_pdb]); // Zhichao's routine to read pdb into rmodel.
		//			// fprintf(stderr,"RASP's sequence: %s\n",rasp_mut.c_str());
		//
		//			// fprintf(stdout,"Loop seq: %s\n",seq);
		//
		//			// This overloaded function inserts the sequence of the loop ("loop_seq") between "start" and "end" residues (PDB indexed)
		//			macromolecule2rasp(protein,rasp_in,rasp_mut,seq,start,end); // Interface between our "Macromolecule" and RASP's "rmodel" (Zhichao Miao)
		////			fprintf(stdout,"rasp_mut (%d): %s\n",rasp_mut.size(),rasp_mut.c_str());
		////			write_mdl("myraspa.pdb",rasp_in);
		////			exit(0);
		//
		//			// Screen residues to get RASP's start/end indices
		//			int rasp_nres = protein->get_num_fragments();
		//			pdbIter *myiter = new pdbIter(protein);
		//			Fragment *res;
		//			for ( myiter->pos_fragment = 0; !myiter->gend_fragment(); myiter->next_fragment() )
		//			{
		//				res = ( Fragment * ) myiter->get_fragment();
		//				if(start==res->getIdNumber())
		//					rasp_start = myiter->pos_fragment;
		//				else
		//					if(end==res->getIdNumber())
		//						rasp_end = myiter->pos_fragment;
		//			}
		//			delete myiter;
		//			fprintf(stderr,"start= %d  rasp_start= %d  end= %d  rasp_end= %d  rasp_nres= %d\n",start,rasp_start,end,rasp_end,rasp_nres);
		//
		//			// Setting re-packing mask
		//			for(int i=0; i<=rasp_start; i++)
		//				rasp_mut[i] = rasp_mut[i] + 32; // upper-case means selected
		//			for(int i=rasp_end; i<rasp_nres; i++)
		//				rasp_mut[i] = rasp_mut[i] + 32;
		//			fprintf(stdout,"Selected loop: %s\n",rasp_mut.c_str());
		//		}

		//		// Zhichao's repacking routine
		//		rasp(rasp_in,rasp_out,rasp_mut);
		//
		//		// Save the complete structure (protein + loop)
		//		char myname[100];
		//		sprintf(myname,"raspa_%02d.pdb",i_pdb+1);
		//		write_mdl(myname,rasp_out);
		//		fprintf(stderr,"RASP re-built loop sidechains in: %s\n",myname);
		//
		//		// Destroy RASP's rmodels
		//		destroy_rmodel(rasp_in);
		//		destroy_rmodel(rasp_out);
		//		rasp_mut.clear();


		// *******************************
		// New RASP stuff initialization
		// *******************************
		CLD *rasp;
		int rasp_start,rasp_end;

		if(rasp_switch && pid==0) // I'm the Master
		{
			// Interface between SBG's "Macromolecule" (ours) and RASP's "CLD object" (Zhichao Miao)
			// Converts the "protein" Macromolecule into a new CLD object, inserting dummy residues in the loop region:
			// ("loop_seq") between "start" and "end" residues (PDB indexed) of selected chain (in_chain)
			rasp = macromol2rasp(protein, rasp_mut, seq+1, start, end, in_chain);
			// dump_rasp( rasp );

			pdbIter *myiter;
			pdbIter *myiter2 = new pdbIter(protein);
			Fragment *res;
			Chain *ch;
			int myresi = 0; // RASP's residue index
			int rasp_nres = protein->get_num_fragments();

			// Screen Chains
			for( myiter2->pos_chain = 0; !myiter2->gend_chain(); myiter2->next_chain() )
			{
				ch = (Chain *) myiter2->get_chain();
				myiter = new pdbIter( ch );

				if( in_chain == ch->getName()[0] ) // Select the input requested chain only
					// Screen residues to get RASP's start/end indices
					for ( myiter->pos_fragment = 0; !myiter->gend_fragment(); myiter->next_fragment() )
					{
						res = ( Fragment * ) myiter->get_fragment();
						if(start == res->getIdNumber())
							rasp_start = myresi;
						else
							if(end == res->getIdNumber())
								rasp_end = myresi;
						myresi++; // count residues
					}
				else
					myresi += myiter->num_fragment(); // count residues of the whole non-matching chain

				delete myiter;
			}
			delete myiter2;

			// RASP start/end indices must be of the first flanking moving residues, not anchors (since anchors are fully known by definition)
			//			rasp_start++;
			//			rasp_end--;

			fprintf(stderr,"start= %d  rasp_start= %d  end= %d  rasp_end= %d  rasp_nres= %d\n",start,rasp_start,end,rasp_end,rasp_nres);

			// Setting re-packing mask
			for(int i=0; i<=rasp_start; i++)
				rasp_mut[i] = rasp_mut[i] + 32; // this sets to lower-case (upper-case means selected for repacking)
			for(int i=rasp_end; i<rasp_nres; i++)
				rasp_mut[i] = rasp_mut[i] + 32;
			fprintf(stdout,"Selected loop for RASP: %s\n",rasp_mut.c_str());

			// Set repackeable residues from string
			//			rasp->GetSub(rasp_mut);
			//			rasp->GetPar();
			//			rasp->GetLib();
		}


		/*ALLOCATE & INITIALIZE LOCAL COORDINATE ARRAYS*/
		Matrix *co      	= new Matrix(nr_atoms_loop,4);/*Loop coordinates read in from loopdata.h*/
		Matrix *coaux     = new Matrix(nr_atoms_loop,4);/*Active loop coordinates read in with 2nd procedure*/
		Matrix *co_dummy     = new Matrix(nr_atoms_loop,4);/*Active loop coordinates read in with 2nd procedure*/
		Matrix *co_i 	= new Matrix(nr_atoms_loop,4);/*Native loop coordinates-3th coordinate would be Van der Waals radius*/

		Matrix *part	= new Matrix(nr_atoms_loop,4);/*array for containing segment coordinates*/
		Matrix *part1	= new Matrix(nr_atoms_loop,4);/*array for containing segment coordinates*/
		Matrix *part2	= new Matrix(nr_atoms_loop,4);/*array for containing segment coordinates*/
		Matrix *conoloop  = new Matrix(nr_atoms_loop,4);/*Active loop coordinates*/
		Matrix *p1coaux   = new Matrix(nr_atoms_loop,4);/*Auxiliary loop coordinates parallel chain 1*/
		Matrix *p2coaux   = new Matrix(nr_atoms_loop,4);/*Auxiliary loop coordinates parallel chain 2*/
		Matrix *co2      	= new Matrix(nr_atoms_loop,4);/*Copy active loop coordinates*/
		Matrix *co_r      = new Matrix(nr_atoms_loop,4);/*Reverse loop coordinates*/
		Matrix *p1co_dummy    = new Matrix(nr_atoms_loop,4);/*Reverse loop coordinates parallel chain 1*/
		Matrix *p1co_r    = new Matrix(nr_atoms_loop,4);/*Reverse loop coordinates parallel chain 1*/
		Matrix *p2co_r    = new Matrix(nr_atoms_loop,4);/*Reverse loop coordinates parallel chain 2*/
		Matrix *co_sol  	= new Matrix(nr_atoms_loop*cycles,3);/*Solution coordinates: size equal to all loops ...*/
		Matrix *co_soli  	= new Matrix(nr_atoms_loop*cycles,3);/*Solution coordinates: size equal to all loops ...*/
		Matrix *p1co_i    = new Matrix(nr_atoms_loop,4);/*Parallel side chain coordinates 1 for loop: initial ones*/
		Matrix *p2co_i    = new Matrix(nr_atoms_loop,4);/*Parallel side chain coordinates 2 for loop: initial ones*/
		Matrix *co_f 	= new Matrix(nr_atoms_loop,4);/*Forward loop coordinates-coordinate would be Van der Waals radius*/
		Matrix *p1co_f    = new Matrix(nr_atoms_loop,4);/*Forward parallel side chain coordinates 1 for loop: initial ones*/
		Matrix *p2co_f    = new Matrix(nr_atoms_loop,4);/*Forward parallel side chain coordinates 2 for loop: initial ones*/
		Matrix *p1co      = new Matrix(nr_atoms_loop,4);/*Parallel side chain coordinates 1 for loop*/
		Matrix *p1co2     = new Matrix(nr_atoms_loop,4);/*Parallel side chain coordinates 1 for loop*/
		Matrix *p1co_sol  = new Matrix(nr_atoms_loop*cycles,3);/*All p1co coordinates of solutions*/
		Matrix *p1co_soli = new Matrix(nr_atoms_loop*cycles,3);/*All p1co coordinates of solutions*/
		Matrix *p2co      = new Matrix(nr_atoms_loop,4);/*Parallel side chain coordinates 2 for loop*/
		Matrix *p2co2     = new Matrix(nr_atoms_loop,4);/*Parallel side chain coordinates 2 for loop*/
		Matrix *p2co_sol  = new Matrix(nr_atoms_loop*cycles,3);/*All p1co coordinates of solutions*/
		Matrix *p2co_soli = new Matrix(nr_atoms_loop*cycles,3);/*All p1co coordinates of solutions*/
		Matrix *path 	= new Matrix(size_path,3);	     /*Array to save an actual loop closure*/
		Matrix *s1co      = new Matrix(nr_atoms_loop+3,4);     /*Parallel hydrogens-on-side chain-carbon coordinates*/

		for(int i = 0; i < nr_atoms_loop; i++)	/*Initialize arrays to zero*/
			for(int j = 0; j < 3; j++)
			{
				part->el[i][j] 	= 0.0;
				part1->el[i][j]  	= 0.0;
				part2->el[i][j]  	= 0.0;

				co_i->el[i][j] 	= 0.0;
				p1co_i->el[i][j] 	= 0.0;

				co_f->el[i][j] 	= 0.0;
				p1co_f->el[i][j] 	= 0.0;

				coaux->el[i][j]	= 0.0;
				p1coaux->el[i][j] = 0.0;

				co_dummy->el[i][j]	= 0.0;
				p1co_dummy->el[i][j]	= 0.0;

				co_r->el[i][j]   	= 0.0;
				p1co_r->el[i][j]  = 0.0;

				co->el[i][j]   	= 0.0;
				p1co->el[i][j] 	= 0.0;

				conoloop->el[i][j]= 0.0;
				co2->el[i][j]    	= 0.0;

				s1co->el[i][j]	= 0.0;	/*hydrogens on side chain-carbon*/
			}

		for(int i = nr_atoms_loop; i < nr_atoms_loop+3; i++)	/*Initialize arrays to zero*/
			for(int j = 0; j < 3; j++)
				s1co->el[i][j]	= 0.0;			        /*hydrogens on side chain-carbon*/

		int *marker	= new int[nr_atoms_loop];     // 1 = bond can rotate, 0 = bond cannot rotate 2= gly*/
		int *marker_r	= new int[nr_atoms_loop];
		int *marker2	= new int[nr_atoms_loop];  // 1 phi y psi 0 // Mon: Watch out this memory allocation!
		int *marker2_r	= new int[nr_atoms_loop];
		int *marker2_i	= new int[nr_atoms_loop];
		int *psimarker	= new int[nr_atoms_loop];	/*marker indicating which bonds are psi bonds*/
		int *psimarker_i	= new int[nr_atoms_loop];	/*marker indicating which bonds are psi bonds*/
		int *psimarker_r	= new int[nr_atoms_loop];	/*marker indicating which bonds are psi bonds*/
		int *marker_i    = new int[nr_atoms_loop];

		int *residuemarker= new int[nr_atoms_loop/3];	//  Marker in from PDB-file for residues:
		//  0 = regular residue, 1 = glycine residue, 2 = proline residue
		for(int i = 0; i < nr_atoms_loop/3; i++)
			residuemarker[i] 	= 0; //  Marker initialization

		// Loading loop coordinates from pdb
		if (present_seq == 1) // If "-f" option enabled
		{
			loop_co(loop, coaux->el,p1co_i->el, residuemarker, seq_num, rmsd_O);
			if(pid==0 && !rmsd_O) // O-coordinates disabled in RMSD
				printf("rcd> Remark: Oxygen coordinates of backbone atoms are ignored! \n");
		}
		else
		{
			// loop_co(ntloop, coaux->el,p1co_i->el, residuemarker, seq_num);
			residuemarker_seq( ntloop, residuemarker, seq_num);

			// Mon: load Nt and Ct directly into "coaux" here...
			// Load Nt coordinates into "coaux"
			Tcoor pos;
			pdbIter *iter;
			iter = new pdbIter( ntloop );
			iter->pos_atom=0; // N
			(iter->get_atom())->getPosition(pos);
			coaux->el[0][0] = pos[0];
			coaux->el[0][1] = pos[1];
			coaux->el[0][2] = pos[2];
			iter->pos_atom=1; // CA
			(iter->get_atom())->getPosition(pos);
			coaux->el[1][0] = pos[0];
			coaux->el[1][1] = pos[1];
			coaux->el[1][2] = pos[2];
			iter->pos_atom=2; // C
			(iter->get_atom())->getPosition(pos);
			coaux->el[2][0] = pos[0];
			coaux->el[2][1] = pos[1];
			coaux->el[2][2] = pos[2];
			delete iter;
			// Other "coaux" values are zero initalized...

			// Load Ct coordinates into "coaux"
			iter = new pdbIter( ctloop );
			iter->pos_atom=0; // N
			(iter->get_atom())->getPosition(pos);
			coaux->el[nr_atoms_loop-3][0] = pos[0];
			coaux->el[nr_atoms_loop-3][1] = pos[1];
			coaux->el[nr_atoms_loop-3][2] = pos[2];
			iter->pos_atom=1; // CA
			(iter->get_atom())->getPosition(pos);
			coaux->el[nr_atoms_loop-2][0] = pos[0];
			coaux->el[nr_atoms_loop-2][1] = pos[1];
			coaux->el[nr_atoms_loop-2][2] = pos[2];
			iter->pos_atom=2; // C
			(iter->get_atom())->getPosition(pos);
			coaux->el[nr_atoms_loop-1][0] = pos[0];
			coaux->el[nr_atoms_loop-1][1] = pos[1];
			coaux->el[nr_atoms_loop-1][2] = pos[2];
			delete iter;
			// Other "coaux" values are zero initalized...

		}
		// What happens when loop-coordiantes are not present ???????????? Check this....

		//		// Getting First-Phi and Last-Psi dihedral angles (forward) --> Phi/Psi from anchors
		double FirstPhi = get_dihedral(NtC,coaux->el[0],coaux->el[1],coaux->el[2]);
		correctdihedral(&FirstPhi);
		double FirstPsi;
		double LastPsi = get_dihedral(coaux->el[nr_atoms_loop-3],coaux->el[nr_atoms_loop-2],coaux->el[nr_atoms_loop-1],CtN);
		correctdihedral(&LastPsi);

		if(pid==0)
			fprintf(stderr,"rcd> FirstPhi= %f (%f)  LastPsi= %f (%f)\n",FirstPhi,FirstPhi*180/M_PI,LastPsi,LastPsi*180/M_PI);

		//
		// MARKERS
		//
		for (int i = 0; i <  nr_atoms_loop; i++)
		{
			marker[i]   	=	0;
			marker_r[i]	= 0;
			marker_i[i]	= 0;
			psimarker[i]	=   0;
			psimarker_r[i]	=   0;

			marker2[i]	= 0;
			marker2_r[i]	= 0;
			marker2_i[i]	= 0;
		}

		// MAKE FORWARD ARRAYS
		for (int i = 0; i < nr_atoms_loop/3; i++)
		{
			marker[3*i] 	=	0;						/*Assignment of OMEGA angle-inactive*/
			marker[3*i + 1]	=	1;						/*Assignment of PHI angle-active*/
			marker[3*i + 2]	=   1;						/*Assignment of PSI angle-active*/
			marker2[3*i + 1]	=	1;						/*Assignment of PHI angle-active*/
			marker2[3*i + 2]	=   1;						/*Assignment of PSI angle-active*/
			psimarker[3*i + 2]= 1;

			if (residuemarker[i]==1) // gly
			{
				marker2[3*i + 1]	=	2;
				marker2[3*i + 2]	=   2;
			}
			else if (residuemarker[i]==2) // pro
			{
				if(ramachandran_inside == 2 && ramachandran_check == 2)
					marker[3*i + 1]	=	1; // We move PHI angle of PRO (Dunbrack maps have info about PRO's PHI)
				else
					marker[3*i + 1]	=	0; // PRO has no PHI... (without Dunbrack)
			}
		}

		// Mon (13/7/2015) Bug...
		if(no_ntpsi) // If Nt-Oxigen optimization is disabled
		{
			marker[2] = 0; // First Psi is not mobile !!! (Psi of Nt anchor is (approximately) known given Nt-Oxigen)
			if(pid==0)
				fprintf(stderr,"rcd> Nt-Oxigen optimization is disabled\n");
		} else if(pid==0)
			fprintf(stderr,"rcd> Nt-Oxigen optimization is enable\n");



		if(no_ctphi) // Disable Ct-anchor Phi angle to freeze CA and C of last loop residue and do it like Rosetta.
		{
			marker[nr_atoms_loop-2] = 0; // Ct-anchor Phi
			if(pid==0)
				fprintf(stderr,"rcd> Ct-Oxigen optimization is disabled\n");
		} else if(pid==0)
			fprintf(stderr,"rcd> Ct-Oxigen optimization is enable\n");


		// MAKE REVERSE ARRAYS BASED ON FORWARD ARRAYS
		for (int z = 1; z < nr_atoms_loop; z++)		/*Markers initialization in reverse direction*/
		{
			marker_r[z] 	=	marker[nr_atoms_loop-z];		//:TRICKY:	MUST BE SHIFTED ONE HIGHER THE REVERSED ELEMENTS
			marker2_r[z]	= 	marker2[nr_atoms_loop-z];
			psimarker_r[z]	= 	psimarker[nr_atoms_loop-z];	//:HOW ABOUT THE NATIVE PROTOCOL BELOW?
		}

		//marker_r[nr_atoms_loop-2] = 1;
		//marker_r[2] = 0;


		for(int i = 0; i < nr_atoms_loop; i++)					/*Initialize & backup values*/
		{
			marker_i[i]=marker[i];  // backup
			psimarker_i[i]=psimarker[i];
			marker2_i[i]=marker2[i];
		}

		//NOW USING NEW READING PROTOCOL
		for(int i = 0; i < nr_atoms_loop; i++)
			for(int j = 0; j < 3; j++)
				co_i->el[i][j] =	coaux->el[i][j];				/*Initialize co_i array from new loading function   */

		//CREATING ANCHORS
		anchor_end->copy(0,0,&co_i->el[last],3,3);			//:INFO: Initialization active fixed end anchor (set to forward anchor initially) <-> must be in agreement with reverse =  0 setting initially?

		anchor_end_i->copy(0,0,&co_i->el[last],3,3);		// Initialization fixed forward end anchor

		for (int i = 0; i < 3; i++)							// Initialization fixed backward end anchor
			for(int j = 0; j < 3; j++)
				anchor_end_r->el[i][j] = co_i->el[2-i][j];  // Reversed Nt anchor (to calculate rmsd in backward mode)

		// ***************************
		// CLASHES MAP STUFF... (GRID)
		// ***************************
		Tcoor centroid;
		double xam;	      /*geometrical center of backbone*/
		double max1[3], max2[3];   /*distance with low end loop*/
		double Anchor_dist=0;
		for(int i = 0; i < 3; i++)
		{
			centroid[i] =  coaux->el[0][i] + coaux->el[nr_atoms_loop-1][i];
			Anchor_dist +=   (coaux->el[0][i] - coaux->el[nr_atoms_loop-1][i])*(coaux->el[0][i] - coaux->el[nr_atoms_loop-1][i]);
		}
		for(int i = 0; i < 3;i++)
		{
			centroid[i] /= 2.0;
			max1[i] =	fabs(coaux->el[0][i] - centroid[i]);
			max2[i] =	fabs(coaux->el[nr_atoms_loop-1][i] - centroid[i]);
		}

		xam = max1[0];	/*select greatest maximum*/
		for(int i = 0; i < 3; i++)
		{
			if(max1[i] > xam)
				xam = max1[i];

			if (max2[i] > xam)
				xam = max2[i];
		}
		xam *= 2.2;
		if(pid==0)
		{
			printf("rcd> Creating clash map with stepsize %f\n",stepsize);// stepsize
			printf("rcd> Anchor Distance %8.3f %8.3f\n",sqrt(Anchor_dist), xam);
		}

		// distance between aa 3...
		if (nr_atoms_loop/3.0*1.8 > xam)
			xam += (nr_atoms_loop/3.0*1.8 - xam)/2.0 ;

		if(pid==0)
			printf("rcd> Box added distance %8.3f\n",xam);

		// CLASHES MAP (PENDING REMOVE around NT & CT ????
		Conditions *condC= new Conditions();
		Condition *condc= new Condition(-1,-1,-1,-1,sel_Ch,sel_Ch,start,end,-1,-1); // Non-loop selection (it does not include anchors)
		condC->add(condc);
		mol=protein->select_cpy(condC, false);  // full coordinates without anchors (except mobile loop)
		Macromolecule *box = mol->select_Box(centroid,xam);

		Macromolecule *box2;
		vlVolume *vol;
		if(sidechain == 3) // Whatever sidechain present in input PDB
		{
			vol = box->project_radiusVDW(stepsize,1,false,softness); // softness --> VdW radius factor
		}
		else
		{
			Condition *ncondT= new Condition(-1,-1,-1,-1,-1,-1,-1,-1,-1,-1);
			Conditions *ncondsT= new Conditions();
			ncondT->add(" N  ");
			ncondT->add(" CA ");
			ncondT->add(" C  ");
			ncondT->add(" O  ");
			if(sidechain > 0) // including CB
				ncondT->add(" CB ");
			if(sidechain > 1) // including some other atoms
			{
				ncondT->add(" H "); // With Hydrogen atoms it seems to improve...
				ncondT->add(" HA "); // With Hydrogen atoms it seems to improve...
				ncondT->add(" CG "); // With Gamma atoms it seems to improve...
				ncondT->add(" SG "); // With Gamma atoms it seems to improve...
			}
			ncondsT->add(ncondT);
			box2 = box->select(ncondsT,true,true); // positive selection of protein/nucleic but including ALL Heteros
			vol = box2->project_radiusVDW(stepsize,1,false,softness); // softness --> VdW radius factor
		}

		if(pid==0 && output_debug)
		{
			FOPS::write(vol,"box.sit");
			if(sidechain == 3) // Whatever sidechain present in input PDB
				box->writePDB("box.pdb");
			else
				box2->writePDB("box.pdb");
		}

		sizegrid = pow(vol->dim().x()*vol->dim().y()*vol->dim().z(),(1/3.0));
		if(pid==0)
			printf("rcd> Clashes map dimensions: %d %d %d vol %d\n", vol->dim().x(),vol->dim().y(),vol->dim().z(),sizegrid);

		unitsx = vol->units().x();
		unitsy = vol->units().y();
		unitsz = vol->units().z();
		vlDim m_dimLimit = vol->dim();
		limitx = m_dimLimit.x();
		limity = m_dimLimit.y();
		limitz = m_dimLimit.z();
		vlStep m_step = vol->stepping();
		stepy = m_step.y();
		stepz = m_step.z();
		float posV[3];
		int stepV[3];
		float unitV[3];
		int limitV[3];
		unitV[0] = unitsx;
		unitV[1] = unitsy;
		unitV[2] = unitsz;
		limitV[0] = limitx;
		limitV[1] = limity;
		limitV[2] = limitz;
		stepV[0] = m_step.x();
		stepV[1] = stepy;
		stepV[2] = stepz;
		vlPoint3f posVol;
		vol->getPosition(&posVol);
		posV[0] = posVol.x();
		posV[1] = posVol.y();
		posV[2] = posVol.z();




		//
		// Cross-correlation with EM-maps as scoring function
		//
		int *workmask; // Integer mask for fast cross-correlation computation
		vlVolume *workmap; // Cropped and masked input EM-map
		vlVolume *loopmap; // Simulated loop map
		if(em_switch) // EM-map loops scoring enabled
		{
			// Current EM-map file name
			strncpy(tmp_path, pdbfile[i_pdb].string, lon);
			tmp_path[lon]='\0';
			strcat(tmp_path,".mrc");
			tmp_path[lon+4]='\0';

			fprintf(stdout,"rcd> Reading EM-map for Cross-correlation computations: %s\n", tmp_path);
			vlVolume *emmap = FOPS::readFile(tmp_path);
			// FOPS::writeFile(emmap,"emmap.mrc");

			// Resize EM-map to match "vol" map (clashes map) dimensions.
			vlDim dim;
			dim = vol->dim(); // Map sizes
			vlUnit uni;
			uni = vol->units(); // Samplings (voxel size)
			vlPoint3f ori; // Origin (minimum corner)
			vol->getPosition(&ori); // Origins
			vlPoint3f corn; // Corner (maximum) position (real space)
			corn.x( (float)dim.x() * (float)uni.x() + ori.x() );
			corn.y( (float)dim.y() * (float)uni.y() + ori.y() );
			corn.z( (float)dim.z() * (float)uni.z() + ori.z() );

			// EM-map resize (now, the motion should be contained...)
			emmap = FOPS::resize(emmap,ori,corn,true);
			// FOPS::writeFile(emmap,"emmap_box.mrc");
			// FOPS::writeFile(vol,"vol_box.mrc");

			// Copy "vol" map (clashes map)
			vlVolume *maskmap = new vlVolume( vol );

			// Create a mask map of the same size as EM-map
			maskmap = FOPS::interpolate_map(maskmap, emmap->units(), true);
			// FOPS::writeFile(maskmap,"maskmap.mrc");

			// Compute mask of indices (array of integers) with two maps. Values above or below the respective thresholds can be selected.
			workmask = FOPS::indices_mask(maskmap, 0.5, false, emmap, em_thr, true);

			// Returns the dot product between two maps (faster than cross-correlation). Mon made (24/8/2018)
			// float dot = FOPS::dotprod_mask(emmap, maskmap, imask);
			// fprintf(stderr, "dot= %f\n", dot);

			// Create a masked output map from the input map and an integers-array mask
			workmap = FOPS::apply_mask(emmap, workmask);
			delete(emmap);
			// FOPS::writeFile(workmap,"workmap.mrc");

			// Create an empty loop map compatible (same size) with "workmap"
			loopmap = new vlVolume( workmap ); // loop simulated map
			// loopmap->clear(0); // Set all voxels to zero?
			// FOPS::writeFile(loopmap,"loopmap.mrc");

			// exit(0);
		}
		delete box; // "box" not needed anymore

		// PD2 Bump filtering stuff
		float *boxcoord;
		char *boxtype = NULL;
		int *boxres = NULL;
		int loopindex = -1; // index of the first moving residue of the loop (in iterator numeration scheme)
		int npocket;
		if(bump_filter)
		{
			Conditions *myconds = new Conditions();
			Condition *condWithAnchors= new Condition(-1,-1,-1,-1,sel_Ch,sel_Ch,start+1,end-1,-1,-1);
			myconds->add(condWithAnchors);
			Macromolecule *noloop=protein->select_cpy(myconds, false);  // full coordinates with anchors (except mobile loop)
			Macromolecule *pocket = noloop; // Pocket and anchors (whole protein)

			Conditions *myconds2= new Conditions();
			Condition *mycond2= new Condition(-1,-1,-1,-1,-1,-1,-1,-1,-1,-1);
			mycond2->add(" N  ");
			mycond2->add(" CA ");
			mycond2->add(" C  ");
			mycond2->add(" O  ");
			mycond2->add(" CB ");
			myconds2->add(mycond2);

			Macromolecule *pocket2 = pocket->select_cpy(myconds2); // N, CA, C, O, CB selection
			if(pid==0 && output_debug)
				pocket2->writePDB("pocket2.pdb");
			npocket = pocket2->get_num_atoms();

			pd2_type(pocket2, &boxtype); // Creates array of atomic types (N,CA,etc...) for the loop pocket
			if(pid==0 && output_debug)
			{
				for(int i=0; i<pocket2->get_num_atoms(); i++)
					fprintf(stderr,"%d ",boxtype[i]);
				fprintf(stderr,"\n");
			}
			pd2_res(pocket2, &boxres); // Creates array with the residue index of the loop pocket
			if(pid==0 && output_debug)
			{
				for(int i=0; i<pocket2->get_num_atoms(); i++)
					fprintf(stderr,"%3d ",boxres[i]);
				fprintf(stderr,"\n");
			}

			pocket2->coordMatrix(&boxcoord); // Get coordinates
			// pd2_bump(float *coord, char *type, int *res, int natoms); // The "PD2" bump energy calculation for loops and pockets

			loopindex = 1 + resindex(noloop, start); // index of the first moving residue of the loop (in iterator numeration scheme)
			if(pid==0 && output_debug)
				fprintf(stderr,"Residue number %d corresponds to %d \n",start+1,loopindex);

			delete pocket2;
		}

		//ALLOCATE SPINORS
		double *spinor1;											/*1st auxiliary spinor*/
		double *spinor2;											/*2nd auxiliary spinor*/
		double *spinor3;											/*3rd auxiliary spinor*/
		double *spinor4;											/*4th auxiliary spinor*/
		double *spinor5;
		double *spinor6;
		double *spinor7;
		double *full;
		double *fullr;
		double *posespinor; 										/*spinor for correct start orientation						*/
		double *posespinor2;										/*auxiliary spinor for correct start orientation reversely 	*/
		double *between;
		double *bond;												/*Optimization routine allocation*/			/*Allocation integrated approach matrix bestanglesm function*/
		double *rotor;
		double *result;
		double *s3tot;

		spinor1 	= 	(double*) malloc(sizeof(double)*4); 		/*in principle only size 2 integrated approach*/
		spinor2	    = 	(double*) malloc(sizeof(double)*4);  		/*in principle only size 2 integrated approach*/
		spinor3 	= 	(double*) malloc(sizeof(double)*4);
		spinor4 	= 	(double*) malloc(sizeof(double)*4);
		spinor5 	= 	(double*) malloc(sizeof(double)*4);
		spinor6  	= 	(double*) malloc(sizeof(double)*4); 		/*in principle only size 2 integrated approach*/
		spinor7     = 	(double*) malloc(sizeof(double)*4);  		/*in principle only size 2 integrated approach*/
		full    	= 	(double*) malloc(sizeof(double)*4);  		/*in principle only size 2 integrated approach*/
		fullr    	= 	(double*) malloc(sizeof(double)*4);  		/*in principle only size 2 integrated approach*/
		posespinor  = 	(double*) malloc(sizeof(double)*4);
		posespinor2 = (double*) malloc(sizeof(double)*4);
		between 	= 	(double*) malloc(sizeof(double)*4);
		bond  	= 	(double*) malloc(sizeof(double)*4);			/*Optimization routine allocation*/			/*one more than formally needed*/
		rotor 	= 	(double*) malloc(sizeof(double)*4);
		result	= 	(double*) malloc(sizeof(double)*4);
		s3tot		= 	(double*) malloc(sizeof(double)*4);

		for(int i = 0; i < 4; i++)
		{
			spinor1[i] 	  = 	0.0;
			spinor2[i]    = 	0.0;
			spinor3[i]    = 	0.0;
			spinor4[i]    = 	0.0;
			spinor5[i]    = 	0.0;
			spinor6[i]    = 	0.0;
			spinor7[i]    = 	0.0;
			full[i]       = 	0.0;
			between[i]    = 	0.0;
			posespinor[i] = 	0.0;
			posespinor2[i]= 	0.0;
			bond[i]       = 	0.0;
			rotor[i]      = 	0.0;
			result[i]     = 	0.0;
			s3tot[i]	  =     0.0;
		}

		//:TODO:	rename arrays with appropriate name? e.g. precos, presin?
		double *prearray1   		= new double[nr_atoms_loop];	/*Array 1 for precalculated cosines valence angles*/
		double *prearray2   		= new double[nr_atoms_loop];	/*Array 2 for precalculated sines valence angles*/
		double *prearray3   		= new double[nr_atoms_loop];	/*Array 3 for precalculated cosines dihedral angles*/
		double *prearray4   		= new double[nr_atoms_loop];	/*Array 4 for precalculated sines dihedral angles*/
		double *prearray5   		= new double[nr_atoms_loop];	/*Array 4 for precalculated sines dihedral angles*/
		double *prearray6   		= new double[nr_atoms_loop];	/*Array 4 for precalculated sines dihedral angles*/

		//CREATE THE POSE SPINOR
		double **coord;				/*local coordinates first 3 atoms			*/
		double **coord2;			/*auxiliary local coordinates first 3 atoms	*/
		double term1;				/*auxiliary scalar*/
		pointeralloc(&coord,3,3);
		pointeralloc(&coord2,3,3);

		if(pid==0 && output_debug)
			multiPDB_cg(co_i->el,p1co->el,1,nr_loop+2,"co_i.pdb",secu,1,0,in_chain);	    /*full backbone: PDB format?*/ /*All solutions*/

		sprintf(tmp_path,"%s/%s_native.txt",folder_later,current_pdb_name); // pdbfile[i_pdb].string);

		if(native == 1) // Native
		{
			// Taking ICs from PDB
			findlength(co_i->el,bond_length,nr_atoms_loop);
			valencefind(co_i->el,valence_angle,nr_atoms_loop);
			finddihedral(co_i->el,dihedral_angle_i,nr_atoms_loop);
			//			if(pid==0 && output_debug)
			//				dump_IC(bond_length, valence_angle, dihedral_angle_i, nr_atoms_loop, secu, FirstPhi, LastPsi);

			dihedral_angle_i[2] = FirstPhi;
			dihedral_angle_i[0] = LastPsi;
			save_dihedrals(dihedral_angle_i, nr_atoms_loop, tmp_path, pdbfile[i_pdb].sequence); // Save all dihedral angles into text file (Phi, Psi and Omega, in degrees)
		}
		else // No-Native
		{
			if (bench_rmsd) // The middle-loop coordinates from "co_i" are not used if native==0, only for final rmsd computations.
			{
				loop_co(loop2, co_i->el,p1co_i->el, residuemarker, seq_num, rmsd_O);
				//				multiPDB_cg(co_i->el,p1co->el,1,nr_loop+2,"co_i2.pdb",secu,1,0,in_chain);	    /*full backbone: PDB format?*/ /*All solutions*/

				findlength(co_i->el,bond_length2,nr_atoms_loop);
				valencefind(co_i->el,valence_angle2,nr_atoms_loop);
				finddihedral(co_i->el,dihedral_angle_i,nr_atoms_loop);

				dihedral_angle_i[2] = FirstPhi;
				dihedral_angle_i[0] = LastPsi;
				if(pid==0 && !server)
				{
					dump_dihedrals("Native bond_length",bond_length2,nr_atoms_loop);
					dump_dihedrals("Native valence_angle",valence_angle2,nr_atoms_loop);
					dump_dihedrals("Native dihedral_angle",dihedral_angle_i,nr_atoms_loop);
				}
				save_dihedrals(dihedral_angle_i, nr_atoms_loop, tmp_path, pdbfile[i_pdb].sequence); // Save all dihedral angles into text file (Phi, Psi and Omega, in degrees)

			}

			findlength(coaux->el,bond_length2,nr_atoms_loop);
			valencefind(coaux->el,valence_angle2,nr_atoms_loop);
			//	dump_dihedrals("backup (bond_length):",bond_length2,nr_atoms_loop);
			//	dump_dihedrals("backup (valence_angle):",valence_angle2,nr_atoms_loop);
			// Use "loopdata" directly...

			// MON: 8/4/2016 I finally decided to remove this shit:
			//
			//			// MON: fixing this shit...(13/05/2015)
			//			// This should be removed some day... see standard Omegas stuff below and watch out indices before removal!
			//			for (int i = 0; i < nr_atoms_loop; i++)
			//				for (int j = 0; j < 3; j++)
			//					conoloop->el[i][j]	= loopdata[i][j];
			//			finddihedral(conoloop->el,dihedral_angle_i,nr_atoms_loop); // To get "omegas" form 1cruA (do something wiser...)

			dihedral_angle_i[2] = FirstPhi;
			dihedral_angle_i[0] = LastPsi;

			// Mon (13/7/2015) Bug
			double NtO[3];
			if(no_ntpsi) // If Nt-Oxygen optimization is disabled (Psi-Nt)
			{
				NtO[0] = (double) NtOf[0];
				NtO[1] = (double) NtOf[1];
				NtO[2] = (double) NtOf[2];
				FirstPsi = get_dihedral(coaux->el[0],coaux->el[1],coaux->el[2],NtO); // Estimate First-Psi angle from Nt anchor "O" atom
				FirstPsi += M_PI;
				correctdihedral(&FirstPsi);
				if(pid==0)
					fprintf(stderr,"rcd> Remark: Nt-Oxigen optimization is disabled, considering FirstPsi= %f (%f)\n",FirstPsi,FirstPsi*180/M_PI);
				dihedral_angle_i[3] = FirstPsi;
			}

			// Set standard values for Non-Omega ICs

			standic(bond_length,valence_angle,nr_atoms_loop,seq,seq_ct1[0]);			//:INFO: use the default peptide geometry
			//			standic(bond_length,valence_angle,nr_atoms_loop,residuemarker);			//:INFO: use the default peptide geometry

			//			findlength(co_i->el,bond_length,nr_atoms_loop);
			//			valencefind(co_i->el,valence_angle,nr_atoms_loop);

			// RESTORE ANCHORs BOND-LENGTHS AND BOND-ANGLES:
			//
			//  AMINO-TERMINAL ONE
			bond_length[1]		= 	bond_length2[1];	//:TRICKY:   THESE FIRST COORDINATES SHOULD ALWAYS HAVE REMAINED THE SAME
			bond_length[2]		= 	bond_length2[2];
			valence_angle[2]	=   valence_angle2[2];
			//  CARBOXY TERMINAL ONE
			bond_length[nr_atoms_loop -1]		= 	bond_length2[nr_atoms_loop -1];
			bond_length[nr_atoms_loop -2]		= 	bond_length2[nr_atoms_loop -2];
			valence_angle[nr_atoms_loop -1]	=   valence_angle2[nr_atoms_loop -1];

			// Here all standard bond lengths and angles should be ok...
			if (bench_rmsd && pid==0 && !server) // The middle-loop coordinates from "co_i" are not used if native==0, only for final rmsd computations.
			{
				dump_dihedrals("Standard bond_length",bond_length,nr_atoms_loop);
				dump_dihedrals("Standard valence_angle",valence_angle,nr_atoms_loop);
			} // MON: I'm "previous"


			Ctomega[0]=-178; // N-92
			Ctomega[1]= 178; // N-93
			Ctomega[2]= 178; // N-94
			Ctomega[3]= 178; // N-95
			Ntomega[0]=-176; // N-92
			Ntomega[1]= 179; // N-93
			Ntomega[2]=-178; // N-94
			Ntomega[3]= 171; // N-95

			// Set standard values for Omega dihedrals
			for(int i = 3; i < nr_atoms_loop-1; i+=3 ) // i+1 --> OMEGA
			{
				if(pdbfile[i_pdb].sequence[i/3]=='p') // If cis-Proline
					dihedral_angle_i[i+1]	=	0.0; // M_PI/2; // Omega for cis-Pro
				else
				{
					dihedral_angle_i[i+1] = -omega * M_PI / 180; // Omega for trans-Aminoacids
			//		fprintf(stdout,"omega %f %f \n", omega, dihedral_angle_i[i+1] );
			//		if (h3_switch) {
			//	        if ( (i/3 + h3NtShift) < 4) dihedral_angle_i[i+1] =  Ctomega[(i/3 + h3NtShift)-1]*M_PI/180.0;
			//		    	else if ((i/3 + h3CtShift +1) >=nr_loop-2)   dihedral_angle_i[i+1]  = Ntomega[(i/3 -nr_loop +1 + h3NtShift)]*M_PI/180.0;
			//		}

				}
			}
			if(pid==0 && !server) // MON: merge this with previous...
			{
				dump_dihedrals("Standard dihedral_angle",dihedral_angle_i,nr_atoms_loop);
				sprintf(tmp_path,"%s/%s_std.txt",folder_later,current_pdb_name); // pdbfile[i_pdb].string);
				save_dihedrals(dihedral_angle_i, nr_atoms_loop, tmp_path); // Save all dihedral angles into text file (Phi, Psi and Omega, in degrees)
			}

		} // COMPUTED LOOP COORDINATES IF STANDARD COORDINATES

		if(pid==0 && output_debug)
		{
			fprintf(stdout,"\n");
			dump_IC(bond_length, valence_angle, dihedral_angle_i, nr_atoms_loop, secu, FirstPhi, LastPsi);
			fflush(stdout);
		}

		// AUXILIARY ARRAYS FOR POSE SPINOR COMPUTATION
		for (int i = 0; i < nr_atoms_loop-1; i++)		/*native loop initialization*/
		{
			prearray1[i] 	    = cos( (M_PI - valence_angle[i+1])/2.0 );
			prearray2[i] 	    = sin( (M_PI - valence_angle[i+1])/2.0 );
			// fprintf(stderr,"%2d  O atom valence angle  %f\n",(180/M_PI)*(M_PI - valence_angle[i+1])/2.0);
		}
		prearray1[nr_atoms_loop-1] 	    = cos( (M_PI - valence_angle[nr_atoms_loop-3])/2.0 );
		prearray2[nr_atoms_loop-1] 	    = sin( (M_PI - valence_angle[nr_atoms_loop-3])/2.0 );

		// Backup ICs (bond lengths and valence angles) for further randomization (randomizeIC)
		for(int i=0;i<nr_atoms_loop;i++)
		{
			valence_angle_back[i] = valence_angle[i];
			bond_length_back[i] = bond_length[i];
		}

		// CALCULATED transformation (full) for forward
		//
		coord[0][0] = 0.0;			/*First atom (N)*/
		coord[0][1] = 0.0;
		coord[0][2] = 0.0;
		spinor4[0] = 1.0;			/*Initialize start spinor*/
		spinor4[1] = 0.0;
		spinor4[2] = 0.0;
		spinor4[3] = 0.0;
		// Generating CA for Nt anchor
		spinor1[0]    = 1;	// First omega is zero ( cos(0/2) = 1)
		spinor1[1]    = 0;	// First omega is zero ( sin(0/2) = 0)
		spinor2[0]    = 0; // First C-N-CA angle is zero ( cos(M_PI/2) = 0)
		spinor2[1]    = 1; // First C-N-CA angle is zero ( sin(M_PI/2) = 1)
		spinor3prod1(spinor4,spinor1,between);
		spinor3prod2(between,spinor2,spinor4);
		term1 = 2*bond_length[1];
		coord[1][0] = coord[0][0] + term1*(spinor4[1]*spinor4[1] + spinor4[0]*spinor4[0] - 0.5);
		coord[1][1] = coord[0][1] + term1*(spinor4[1]*spinor4[2] + spinor4[0]*spinor4[3]);
		coord[1][2] = coord[0][2] + term1*(spinor4[1]*spinor4[3] - spinor4[0]*spinor4[2]);

		// Generating C for Nt anchor
		spinor1[0]    = 1; // cos(FirstPhi/2.0);	// 1;	// First Phi is considered zero ( cos(0/2) = 1)
		spinor1[1]    = 0; // cos(FirstPhi/2.0);	// 0;	// First Phi is considered zero ( sin(0/2) = 0)
		spinor2[0]    = cos((M_PI-valence_angle[2])/2);
		spinor2[1]    = sin((M_PI-valence_angle[2])/2);
		spinor3prod1(spinor4,spinor1,between);
		spinor3prod2(between,spinor2,spinor4);
		term1 = 2*bond_length[2];
		coord[2][0] = coord[1][0] + term1*(spinor4[1]*spinor4[1] + spinor4[0]*spinor4[0] - 0.5);
		coord[2][1] = coord[1][1] + term1*(spinor4[1]*spinor4[2] + spinor4[0]*spinor4[3]);
		coord[2][2] = coord[1][2] + term1*(spinor4[1]*spinor4[3] - spinor4[0]*spinor4[2]);

		// Computing initial pose ("full") to begin chain construction from Nt by isarotatetree...
		pose(coord,co_i->el,posespinor);
		full[0] = 1.0;				/*Add hydrogen at first nitrogen*/
		full[1] = 0.0;
		full[2] = 0.0;
		full[3] = 0.0;
		quatproduct(posespinor,full,full);
		// REBUILDING THE NATIVE LOOP IF ALL DATA AVAILABLE --- SIMPLE DO APPROACH

		// RECOMPUTATION START LOOP WITH FULLY UP-TO-DATE INTERNAL COORDINATES
		dihedral_angle_i[2] = FirstPhi;
		dihedral_angle_i[0] = LastPsi;
		// dihedral_angle_i[3] = FirstPsi;

		//		if(pid==0)
		//		{
		//			dump_dihedrals("bond_length",bond_length,nr_atoms_loop);
		//			dump_dihedrals("valence_angle",valence_angle,nr_atoms_loop);
		//			dump_dihedrals("dihedral_angle",dihedral_angle_i,nr_atoms_loop);
		//		}

		co->copy(0,0,&co_i->el[0],nr_atoms_loop,3);				/*Need for first coordinate to have in co matrix*/
		isafrotatetreef_cg(co->el,p1co->el,nr_atoms_loop,bond_length,dihedral_angle_i,prearray3,prearray4,prearray5,prearray6,prearray1,prearray2,seq,full,cg_mode,factor_cacb);
		//		dump_dihedrals("rcd> Initial dihedral angles",dihedral_angle_i,nr_atoms_loop);

		if(output_debug)
			multiPDB_cg(co->el,p1co->el,1,nr_loop+2,"regenerated_loop_cg.pdb",secu,1,cg_mode,in_chain);	    /*full backbone: PDB format?*/ /*All solutions*/


		// CALCULATED transformation (full) for backward



		//RESUBSTITUTE THE INITIAL COORDINATES
		//	:TRICKY: OXYGEN ATOMS MUST BE CREATED FOR SURE
		//
		if(pid==0 && (native == 1 || present_seq == 1))
		{
			double sumRMSDO, sumRMSDBB;
			sumRMSDO=0.0;
			for (int j = 0; j < nr_loop+2; j++)	/*iterate # residues*/ /*local j used and override global j*/
				for(int i = 0; i < 3; i++)
					sumRMSDO += pow( p1co->el[j*3+2][i] - p1co_i->el[j*3+2][i], 2 );

			sumRMSDBB=0.0;
			for(int i = 0; i < nr_atoms_loop; i++)
				for(int j = 0; j < 3; j++)
				{
					//			co_i->el[i][j] 	=	  co->el[i][j];			//NOT NEEDED FORMALLY
					sumRMSDBB += pow( co->el[i][j] -co_i->el[i][j], 2 );
					if(cg_mode>0)  {
						// p1co_i->el[i][j] 	= 	p1co->el[i][j];
					}
				}
			fprintf(stderr,"rcd> RMSD native-generated BB %8.4f +O %8.4f\n",sqrt(sumRMSDBB/nr_atoms_loop),sqrt(sumRMSDO/(nr_loop+2)));

			bool debug_clash = false;

			if (clash_cg(vol,co_i->el,p1co_i->el,nr_atoms_loop,marker2,0,cg_mode) == 1) {
				fprintf(stderr,"rcd> Warning: Clash in PDB %s loop\n",pdbfile[i_pdb].string);
				debug_clash = true;
			}
			if (intra_clash_cg(co_i->el,p1co_i->el,factor_co,factor_ocb,nr_atoms_loop,atomshift,marker2,0,cg_mode) == 1) {
				fprintf(stderr,"rcd> Warning: Intra-clash in PDB %s loop\n",pdbfile[i_pdb].string);
				debug_clash = true;
			}
			if (clash_cg(vol,co->el,p1co->el,nr_atoms_loop,marker2,0,cg_mode) == 1) {
				fprintf(stderr,"rcd> Warning: Clash in regenerated %s loop\n",pdbfile[i_pdb].string);
				debug_clash = true;
			}
			if (intra_clash_cg(co->el,p1co->el,factor_co,factor_ocb,nr_atoms_loop,atomshift,marker2,0,cg_mode) == 1) {
				fprintf(stderr,"rcd> Warning: Intra-clash in regenerated %s loop\n",pdbfile[i_pdb].string);
				debug_clash = true;
			}
			if(output_debug || debug_clash)
			{
				sprintf(file_name,"%s/%s_loop_cg_re.pdb",folder_later,current_pdb_name); // pdbfile[i_pdb].string);
				multiPDB_cg(co->el,p1co->el,1,nr_loop+2,file_name,secu,1,cg_mode);
				sprintf(file_name,"%s/%s_loop_cg_i.pdb",folder_later,current_pdb_name); // pdbfile[i_pdb].string);
				multiPDB_cg(co_i->el,p1co_i->el,1,nr_loop+2,file_name,secu,1,cg_mode);
			}
		}

		// Do something elegant with "omegas"
		for (int z = 0; z < nr_atoms_loop; z++)						/*Initialize forward dihedral array*/
			dihedral_angle_f[z] = dihedral_angle_i[z];

		for (int z = 1; z < nr_atoms_loop-2; z++)					/*Initialize reverse dihedral array*/
			dihedral_angle_r[z+2] = dihedral_angle_i[nr_atoms_loop-z];

		// Loading missing dihedral angles: FirstPhi and LastPsi (in anchors)
		dihedral_angle_f[2] = FirstPhi; // forward
		dihedral_angle_f[0] = LastPsi; // forward
		dihedral_angle_r[0] = FirstPhi; // backward
		dihedral_angle_r[2] = LastPsi; // backward

		//CREATE CORRECT INITIALIZED ARRAYS
		for (int z = 0; z < nr_atoms_loop; z++)
			for (int x = 0; x < 3; x++)
			{
				co_f->el[z][x] 	= 	  co->el[z][x];
				co_r->el[z][x] 	= 	  co->el[nr_atoms_loop-1-z][x];

				if(cg_mode>0)
				{
					p1co_f->el[z][x]	= 	p1co->el[z][x];
					p1co_r->el[z][x]	= 	p1co->el[nr_atoms_loop-1-z][x];
				}
			}

		if(load_ramap)
		{
			//    Dunbrack based PDFs
			if( !(iaa = (int *) malloc(sizeof(int)*(nr_loop+2))) )
			{
				if(pid==0)
					printf("rcd> Error: Sorry, unable to allocate iaa memory!!!\n");
				exit(1);
			}

			int seq_ct1ID, seq_nt1ID;

			for (int j = 0; j < 20;j++)
				if(*(seq_ct1)==AA[j].aa_name1)
				{
					seq_ct1ID= j;
					break;
				}
			for (int j = 0; j < 20;j++)
				if(*(seq_nt1)==AA[j].aa_name1)
				{
					seq_nt1ID= j;
					break;
				}

			// Loop AA ID (including anchors)
			for (int i = 0; i < nr_loop+2;i++)
			{
				if(pdbfile[i_pdb].sequence[i]=='p') // If cis-Proline
				{
					if(pid==0)
						printf("rcd> Remark: cis-Proline requested by user in residue %d\n",i+1);
					iaa[i] = 20; // cis-Proline index
				}
				else // other aminoacids
					for (int j = 0; j < 20;j++)
						if(pdbfile[i_pdb].sequence[i]==AA[j].aa_name1)
						{
							iaa[i]= j; // Array with ID number of AA except C-anchor +1 or N-anchor -1 (added previously in line 952)
							break;
						}
			}

			if( !(maps = (float ***) malloc(sizeof(float **)*(nr_loop+2)) ))
			{
				printf("Sorry, unable to allocate maps memory!!!\n");
				exit(1);
			}
			if( !(ramaps = (float ***) malloc(sizeof(float **)*(nr_loop+2)) ))
			{
				printf("Sorry, unable to allocate maps memory!!!\n");
				exit(1);
			}
			if( !(binmaps = (float ***) malloc(sizeof(float **)*(nr_loop+2)) ))
			{
				printf("Sorry, unable to allocate maps memory!!!\n");
				exit(1);
			}
			//			if( !(binmaps_post = (float ***) malloc(sizeof(float **)*(nr_loop+2)) ))
			//			{
			//				printf("Sorry, unable to allocate binmaps_post memory!!!\n");
			//				exit(1);
			//			}

			// NT anchor AA
			int ileft,iright; // left and right indices (required for cis-Pro)
			if(seq_nt1ID==20) // If cis-Pro to the left... (all Prolines are merged when Left or Right side, see Dunbrack's paper...)
				ileft = 12; // left-side cis-Pro is set back to trans-Pro (code 12)
			else
				ileft = seq_nt1ID;
			if(iaa[1]==20) // If cis-Pro to the right... (all Prolines are merged when Left or Right side, see Dunbrack's paper...)
				iright = 12; // right-side cis-Pro is set back to trans-Pro (code 12)
			else
				iright = iaa[1];
			// fprintf(stderr,"iaa[0]= %d  seq_nt1ID= %d  iaa[1]= %d\n",iaa[0],seq_nt1ID,iaa[1]);
			//			maps[0] = map_gen(iaa[0],ileft,iright, pdfs, size); // PDF map of anchor Nt
			maps[0] = map_gen(iaa[0],ileft,iright, pdfs, size, pdfs[iaa[0]][1][21]); // PDF map of anchor Nt
			norm_map(maps[0],size);

			if( emodel == 6 || emodel == 56 ) // Rama energy models
				ramaps[0] = allmaps[iaa[0]][ileft][iright]; // Rama energy map for Nt-anchor

			// CT anchor AA
			if(iaa[nr_loop]==20) // If cis-Pro to the left... (all Prolines are merged when Left or Right side, see Dunbrack's paper...)
				iright = 12; // left-side cis-Pro is set back to trans-Pro (code 12)
			else
				iright = iaa[nr_loop];
			if(seq_ct1ID==20) // If cis-Pro to the right... (all Prolines are merged when Left or Right side, see Dunbrack's paper...)
				iright = 12; // right-side cis-Pro is set back to trans-Pro (code 12)
			else
				iright = seq_ct1ID;
			// fprintf(stderr,"Ct iaa[nr_loop+1]= %d  iaa[nr_loop]=%d  seq_ct1ID=%d   size=%d\n",iaa[nr_loop+1],iaa[nr_loop],seq_ct1ID, size);
			//			maps[nr_loop+1] = map_gen(iaa[nr_loop+1],ileft,iright, pdfs, size); // PDF map of anchor Ct
			maps[nr_loop+1] = map_gen(iaa[nr_loop+1],ileft,iright, pdfs, size, pdfs[iaa[nr_loop+1]][1][21]); // PDF map of anchor Ct
			norm_map(maps[nr_loop+1],size);

			if( emodel == 6 || emodel == 56 ) // Rama energy models
				ramaps[nr_loop+1] = allmaps[iaa[nr_loop+1]][ileft][iright]; // Rama energy map for Ct-anchor

			// LOOP AAs
			for(int i=1;i<nr_loop+1;i++)
			{
				// fprintf(stderr,"AA=%d  iaa[i]=%d iaa[i-1]=%d  iaa[i+1]=%d  size=%d\n",i,iaa[i],iaa[i-1],iaa[i+1],size);
				if(iaa[i-1]==20) // If cis-Pro to the left... (all Prolines are merged when Left or Right side, see Dunbracks paper...)
					ileft = 12; // left-side cis-Pro is set back to trans-Pro (code 12)
				else
					ileft = iaa[i-1];
				if	(iaa[i+1]==20) // If cis-Pro to the right... (all Prolines are merged when Left or Right side, see Dunbracks paper...)
					iright = 12; // right-side cis-Pro is set back to trans-Pro (code 12)
				else
					iright = iaa[i+1];
				//				maps[i] = map_gen(iaa[i],ileft,iright, pdfs, size);
				maps[i] = map_gen(iaa[i],ileft,iright, pdfs, size, pdfs[iaa[i]][1][21]);
				norm_map(maps[i],size);

				if( emodel == 6 || emodel == 56 ) // Rama energy models
					ramaps[i] = allmaps[iaa[i]][ileft][iright]; // Rama energy map for mobile aminoacids

				//				char kk[100];
				//				sprintf(kk,"mierda_%02d.txt",i);
				//				show_map(maps[i],size,kk);
			}

			// COMPUTING NATIVE ENERGY
			// MON: Consider move this to outside Rama region...
			if(native==1 || bench_rmsd)
			{
				if(bump_filter)
				{
					// Warning: this only works with "-f", check why?
					float extra = pd2_extrabump(boxcoord, boxtype, boxres, npocket, co_i->el, p1co_i->el, nr_atoms_loop, loopindex, cg_mode, residuemarker);
					float intra = pd2_intrabump(co_i->el, p1co_i->el, nr_atoms_loop, cg_mode, residuemarker);
					native_energy = intra + extra;
					if(pid==0)
						fprintf(stderr,"rcd> Native loop Bump energy: %f (extra: %f  intra: %f)\n",native_energy,extra,intra);
				}
				//				else // OLD RAMA ENERGY
				//				{
				//					native_energy = rama_energy(dihedral_angle_i, nr_atoms_loop, maps, size);
				//					if(pid==0)
				//						fprintf(stderr,"rcd> Native loop Rama energy: %f\n",native_energy);
				//				}

				if(emodel == 1) // ICOSA
				{
					// Convert "co->el" 2D-matrix into linear array (coordMatrix-like)
					co2coord(co_i->el + 3, loop_coord,  nr_atoms_loop-6);

					// Mon: use "nr_loop"
					native_loco = loco_energy(loop_coord, iseql, nr_atoms_loop/3-2, coorde, iseqe, nseqe , anchorCt, loco, icosVertices);
				}

				if( emodel == 5 || emodel == 56 ) // KORP energy models
				{
					// Convert "co->el" 2D-matrix into linear array (coordMatrix-like)
					co2coord(co_i->el + 3, loop_coord,  nr_atoms_loop-6);

					// Loop vs. Loop
					framesloop = frameCoord(loop_coord, nr_loop, resnumsloop, reschainidsloop);
					ncont = contactCoord(&contacts, korpmap->cutoff, framesloop, nr_loop, iseql);
					native_loco = korp6D(contacts,ncont,korpmap); // Non-bonding = (nb2 < |i-j|)
					free(contacts);
					contacts = NULL; // Enables automatic memory allocation on next iteration

					// Loop vs. Environment
					frames = frameCoord(coorde, nseqe, resnums, reschainids);
					ncont = contactCoord(&contacts, korpmap->cutoff, frames, nseqe, iseqe, anchorNt, framesloop, nr_loop, iseql);
					free(framesloop); // Free loop frames
					native_loco += korp6D(contacts,ncont,korpmap); // Non-bonding = (nb2 < |i-j|)
					free(contacts);
					contacts = NULL; // Enables automatic memory allocation on next iteration
					if(pid==0)
						fprintf(stderr,"rcd> Native loop KORP energy: %f\n",native_loco);
				}

				if( emodel == 6 || emodel == 56 ) // Rama energy models
				{

					native_energy = rama_factor * rama_energy2(dihedral_angle_i, nr_atoms_loop/3, ramaps, size);
					if(pid==0)
						if(emodel == 6)
							fprintf(stderr,"rcd> Native loop Rama energy: %f\n",native_energy);
						else
							fprintf(stderr,"rcd> Native loop KORP+RAMA energy: %f\n",native_energy + native_loco);

				}
			}
			else
			{
				// Just compute Environment Frames for KORP stuff...
				if( emodel == 5 || emodel == 56 ) // KORP energy models
					frames = frameCoord(coorde, nseqe, resnums, reschainids);
			}

			if( !(accs = (float **) malloc(sizeof(float *)*(nr_loop+2))) ) // We reserve memory (including anchors)
			{
				printf("Sorry, unable to allocate accs memory!!!\n");
				exit(1);
			}

			if( !(rama_cutoff = (float *) malloc(sizeof(float)*(nr_loop+2))))
			{
				printf("Sorry, unable to allocate rama_cutoff memory!!!\n");
				exit(1);
			}

			// Selecting the appropriate Ramachandran threshold for Terminal residues
			float rama_thr_new;


			if(trama_switch)
				rama_thr_new = rama_thrA;
			else {
				 rama_thr_new = rama_thr;
				 if (h3_switch) rama_thr_new = rama_thrK;
			}



			////////////////////////////////////////////////
			// Building "binmaps" for First and Last AA
			///////////////////////////////////////////////
			if(exact_rama)
			{
				// NT anchor AA /////////////////////////////////////
				if (h3_switch && h3NtShift<4)
				{
					binmaps[0] = map_mergeh3(h3maps, maps[0], (21+h3NtShift),  size, kink_mergeA); // pos 92
				}
				else {
					binmaps[0] = map_copy(maps[0], size); // Copy Nt map
				}
				rama_cutoff[0] = map_cutoff(binmaps[0], rama_thr_new, size);  // NT anchor AA
				map_setzero(binmaps[0], size, rama_cutoff[0], true); // Set to zero those elements with values lower than some threshold

				// CT anchor AA
				if (h3_switch && h3CtShift<4)
				{
					binmaps[nr_loop+1] = map_mergeh3(h3maps, maps[nr_loop+1], (28-h3CtShift),  size, kink_mergeA); // pos 92
					// binmaps[nr_loop+1] = map_copyh3(h3maps, (28-h3CtShift),  size);
				}
				else {
					binmaps[nr_loop+1] = map_copy(maps[nr_loop+1], size); // Copy Ct map
				}
				rama_cutoff[nr_loop+1] = map_cutoff(binmaps[nr_loop+1], rama_thr_new, size);  // CT anchor AA
				map_setzero(binmaps[nr_loop+1], size, rama_cutoff[nr_loop+1], true); // Set to zero those elements with values lower than some threshold
			}
			else
			{
				binmaps[0] = map_bin2(maps[0], rama_thr_new, size); // NT anchor AA
				rama_cutoff[0] = 0.5;
				binmaps[nr_loop+1] = map_bin2(maps[nr_loop+1], rama_thr_new, size); // CT anchor AA
				rama_cutoff[nr_loop+1] = 0.5;
			}

			int nbad,nbad2,iphi,ipsi,ibest,ibest2;





			// Checking whether first Psi and last Phi have any valid value...
			//
			// Nt anchor
			//
			float mysum;
			nbad = check1angle(binmaps[0], dihedral_angle_i[2], true, size, 0, true, rama_cutoff[0]); // Nt anchor
			if(pid==0)
				fprintf(stderr,"rcd> Remark: Phi angle in Nt-anchor (%f) has %d valid Psi angles.\n",dihedral_angle_i[2]*180/M_PI, size-nbad);
			if(nbad >= size) // Fix the Psi anchor angle
			{
				// Get the most likely Phi, if bestphi=true, (or Psi if bestphi=false) for a given Psi (or Phi) "index" from a PDF distribution "map" of side "size".
				iphi = angle2index(dihedral_angle_i[2],size);
				if (h3_switch && h3NtShift<4) {
					binmaps[0] = map_mergeh3(h3maps, maps[0], (21+h3NtShift),  size, kink_mergeA); // pos 92 C
				}
				else {
					binmaps[0] = map_copy(maps[0], size); // Copy Nt map
				}
				ibest = best_index(binmaps[0],false,iphi,size); // res 92 H3

				dihedral_angle_i[3] = index2angle( ibest, size ); // Set best Psi (fixed)
				dihedral_angle_f[3] = dihedral_angle_i[3]; // this required given loop generation does not change angles if marker==1
				dihedral_angle_r[nr_atoms_loop-1] = dihedral_angle_i[3]; // this required given loop generation does not change angles if marker==1
				if(pid==0)
					fprintf(stderr,"rcd> Warning: First Psi angle auto-fixed to the most likely: %.1lf +- 2.5 deg! (%lf rad) (ibest=%d)\n", dihedral_angle_i[3]*180/M_PI, dihedral_angle_i[3], ibest);

				// NT anchor AA

				rama_cutoff[0] = binmaps[0][iphi][ibest];

				// this Rama cutoff must be re-defined
				map_setzero(binmaps[0], size, rama_cutoff[0], true); // Set to zero those elements with values lower than some threshold
			}
			accs[0] = NULL;
			accu_1D(binmaps[0],size,FirstPhi,true,&accs[0]); // Nt moves PSI (PHI is fixed, then isphi=true)
			//
			// Ct anchor
			//
			nbad2 = check1angle(binmaps[nr_loop+1], dihedral_angle_i[0], false, size, 0, true, rama_cutoff[nr_loop+1]); // Ct anchor
			if(pid==0)
				fprintf(stderr,"rcd> Remark: Psi angle in Ct-anchor (%f) has %d valid Phi angles.\n",dihedral_angle_i[0]*180/M_PI, size-nbad2);
			if(nbad2 >= size) // Fix the phi anchor angle
			{

				// Get the most likely Phi, if bestphi=true, (or Psi if bestphi=false) for a given Psi (or Phi) "index" from a PDF distribution "map" of side "size".
				ipsi = angle2index(dihedral_angle_i[0],size);
				if (h3_switch && h3CtShift<4)
					binmaps[nr_loop+1] = map_mergeh3(h3maps, maps[nr_loop+1], (28-h3CtShift),  size, kink_mergeA); // pos 92
				else
					binmaps[nr_loop+1] = map_copy(maps[nr_loop+1], size); // Copy Ct map

				ibest2 = best_index(binmaps[nr_loop+1],true,ipsi,size); // Get best Phi (to be fixed)

				dihedral_angle_i[nr_atoms_loop-1] = index2angle( ibest2, size );
				if(pid==0)
					fprintf(stderr,"rcd> Warning: Last Phi angle auto-fixed to the most likely: %.1lf +- 2.5 deg! (%lf rad) (ibest2=%d)\n", dihedral_angle_i[nr_atoms_loop-1]*180/M_PI, dihedral_angle_i[nr_atoms_loop-1], ibest2);

				// CT anchor AA

				rama_cutoff[nr_loop+1] = binmaps[nr_loop+1][ibest2][ipsi];
				map_setzero(binmaps[nr_loop+1], size, rama_cutoff[nr_loop+1], true); // Set to zero those elements with values lower than some threshold
			}
			accs[nr_loop+1] = NULL;
			accu_1D(binmaps[nr_loop+1],size,LastPsi,false,&accs[nr_loop+1]); // Ct moves PHI (PSI is mobile, then isphi=false)


			if(no_ntpsi)
			{
				int imypsi,imyphi;
				imyphi = angle2index(dihedral_angle_i[2],size); // Nt's Phi angle
				imypsi = angle2index(dihedral_angle_i[3],size); // Nt's Psi angle
				if(binmaps[0][imyphi][imypsi] < rama_cutoff[0])
				{
					binmaps[0][imyphi][imypsi]=rama_cutoff[0];
					if (pid==0)
						fprintf(stderr,"rcd> Warning: N-terminal Phi/Psi angle pair (%f/%f) is outside Ramachandran ...\n",dihedral_angle_i[2],dihedral_angle_i[3]);

				} else
				{
					if (pid==0)
					  fprintf(stderr,"rcd> Nt Phi/Psi angle pair (%f/%f) is ok\n",dihedral_angle_i[2],dihedral_angle_i[3]);

				}

			}

			if(no_ctphi)
			{
				int imypsi,imyphi;
//				imyphi = angle2index(dihedral_angle_i[nr_atoms_loop-1],size); // Ct's Phi angle
//				imypsi = angle2index(LastPsi,size); // Ct's Psi angle
//
//				if(binmaps[0][imyphi][imypsi] < rama_cutoff[0])
//				{
//					binmaps[0][imyphi][imypsi]=rama_cutoff[0];
//					if (pid==0)
//						fprintf(stderr,"rcd> Warning: C-terminal Anchor Phi/Psi angle pair (%f/%f) is outside Ramachandran ...\n",dihedral_angle_i[nr_atoms_loop-1],LastPsi);
//
//				}


				imyphi = angle2index(dihedral_angle_i[nr_atoms_loop-4],size); // Ct's Phi angle
				imypsi = angle2index(dihedral_angle_i[nr_atoms_loop-3],size); // Ct's Psi angle


				if(binmaps[0][imyphi][imypsi] < rama_cutoff[0])
				{
					binmaps[0][imyphi][imypsi]=rama_cutoff[0];
					if (pid==0)
					fprintf(stderr,"rcd> Warning: C-terminal Anchor Phi/Psi angle pair (%f/%f) is outside Ramachandran ...\n",dihedral_angle_i[nr_atoms_loop-4],dihedral_angle_i[nr_atoms_loop-3]);

				} else
					{
					fprintf(stderr,"rcd> Phi/Psi Ct angle pair (%f/%f) ok\n",dihedral_angle_i[nr_atoms_loop-4],dihedral_angle_i[nr_atoms_loop-3]);

					}




			}


			// LOOP rest of AAs
			for(int i=1;i<nr_loop+1;i++)
			{
				// Mon (23/5/2016): Bug solved! RCD versions 1.10 and 1.11 used a Rama. thr. 99% in the first and last loop residues.

				if(trama_switch)
				{
					if(i <= nrama)
						rama_thr_new = rama_thrA;     // put rama_thrA for nrama
					else if(i >= nr_loop+1 - nrama)
						rama_thr_new = rama_thrA;
					else {
						rama_thr_new = rama_thr; // put rama_thr for all or inter kink
						if (h3_switch) {
							if ( (i + h3NtShift) < 4) rama_thr_new = rama_thrK; // inter kink
							else if ((i + h3CtShift) >=nr_loop-2) rama_thr_new = rama_thrK; // inter kink
						}
					}
				}
				else {
					rama_thr_new = rama_thr;
					if (h3_switch) {
						if ( (i + h3NtShift) < 4) rama_thr_new = rama_thrK;
						else if ((i + h3CtShift) >=nr_loop-2) rama_thr_new = rama_thrK;
					}
				}
				// fprintf(stderr,"Residue %2d using %8f as Ramachandran threshold.\n",i,rama_thr_new);

				if (h3_switch)
				{
					// get corresponding H3 maps
					/*  21 pos 92 C
				        22 pos 93 A*
				        23 pos 94 R
				        24 pos 95 X
				        28 pos 100
				        27 pos 101 D
				        26 pos 102 Y
				        25 pos 103 W?
					 */
					if ( (i + h3NtShift) < 4)
					{
						binmaps[i] = map_mergeh3(h3maps, maps[i], (21+i+h3NtShift),  size, kink_merge);
						//rama_cutoff[i] = map_cutoff(binmaps[i], rama_thrK, size);  // cutoff of i-th AA

					}
					else if ((i + h3CtShift) >=nr_loop-2)
					{
						binmaps[i] = map_mergeh3(h3maps, maps[i], (27+i-nr_loop-h3CtShift),  size, kink_merge);
						//rama_cutoff[i] = map_cutoff(binmaps[i], rama_thrK, size);  // cutoff of i-th AA

					}
					else
					{
						// binmaps[i] = map_copyh3(h3maps, iaa[i],  size);
						binmaps[i] = map_mergeh3(h3maps, maps[i], iaa[i],  size, nokink_merge);
					}

					//fprintf(stderr,"Residue %2d aa %d %s %f \n",i,iaa[i], secu[i], rama_cutoff[i]);
					rama_cutoff[i] = map_cutoff(binmaps[i], rama_thr_new, size);  // cutoff of i-th AA
					map_setzero(binmaps[i], size, rama_cutoff[i], true); // Set to zero those elements with values lower than some threshold
				}

				else
				{

					if(exact_rama)
					{
						binmaps[i] = map_copy(maps[i], size); // Copy i-th map
						// fprintf(stderr,"binmap %d\n",i);
						rama_cutoff[i] = map_cutoff(binmaps[i], rama_thr_new, size);  // cutoff of i-th AA
						map_setzero(binmaps[i], size, rama_cutoff[i], true); // Set to zero those elements with values lower than some threshold
					}
					else // if binary maps
					{
						binmaps[i] = map_bin2(maps[i], rama_thr_new, size);
						rama_cutoff[i] = 0.5;
						// NOt sure is wright 2022
					}
				}

				accs[i] = accunorm (binmaps[i],size); //  Acc probability density in loop AA (between anchors).
			}

			//			fprintf(stderr,"Rama_cutoff:\n");
			//			for(int i=0;i<nr_loop+2;i++)
			//				fprintf(stderr, "%f ", rama_cutoff[i]);
			//			fprintf(stderr,"\n");
			//			fprintf(stderr,"Check_angles. Number of bad angles: %d\n",check_angles(binmaps, dihedral_angle_i, size, nr_atoms_loop, 0, true, marker, rama_cutoff));

			if(inside_rand)
			{
				if( !(acc = (float *) malloc(sizeof(float)*size)) ) // We reserve memory for one profile
				{
					printf("Sorry, unable to allocate acc memory!!!\n");
					exit(1);
				}

				if( !(accs_phi= (float ***) malloc(sizeof(float **)*(nr_loop+2))))
				{
					printf("Sorry, unable to allocate accs_phi memory!!!\n");
					exit(1);
				}
				if( !(accs_psi= (float ***) malloc(sizeof(float **)*(nr_loop+2))))
				{
					printf("Sorry, unable to allocate accs_psi memory!!!\n");
					exit(1);
				}

				// Some pre-calculations... (otherwise up to 8x slower!)
				for(int j = 0; j < nr_loop+2; j++) // screen loop aminoacids
				{
					if( !(accs_phi[j]= (float **) malloc(sizeof(float *)*size)))
					{
						printf("Sorry, unable to allocate accs_phi memory!!!\n");
						exit(1);
					}
					if( !(accs_psi[j]= (float **) malloc(sizeof(float *)*size)))
					{
						printf("Sorry, unable to allocate accs_phi memory!!!\n");
						exit(1);
					}

					for( int i=0; i<size; i++) // rows or cols
					{
						if( !(accs_phi[j][i]= (float *) malloc(sizeof(float)*size)))
						{
							printf("Sorry, unable to allocate accs_phi memory!!!\n");
							exit(1);
						}
						if( !(accs_psi[j][i]= (float *) malloc(sizeof(float)*size)))
						{
							printf("Sorry, unable to allocate accs_phi memory!!!\n");
							exit(1);
						}
						accu_1D(binmaps[j],size,i,true,accs_phi[j][i]); // Fixes PHI if isphi=true (PSI is moved)
						accu_1D(binmaps[j],size,i,false,accs_psi[j][i]); // Fixes PSI if isphi=false (PHI is moved)
					}
				}

				//				if(nbad >= size) // Fix the First Psi anchor angle
				//					for( int i=0; i<size; i++)
				//						accs_psi[0][i] = accs[0];
				//
				//				if(nbad2 >= size) // Fix the Last Phi anchor angle
				//					for( int i=0; i<size; i++)
				//						accs_phi[nr_loop+1][i] = accs[nr_loop+1];
			}

			if(pid==0 && output_debug)
			{
				// SHOW DUNBRACK MAPS: Normal & bin
				char name[100];
				char num[10];
				char basename[]="pdf_";
				char ext[]=".map";
				for(int i=0; i<nr_loop+2; i++)
				{
					strcpy(name,basename);
					sprintf(num,"%02d",i);
					strcat(name,num);
					strcat(name,ext);
					show_map(maps[i],size,name);
				}
				if(rama_thr >= 0)
				{
					char basename2[]="pdf_bin_";
					for(int i=0; i<nr_loop+2; i++)
					{
						strcpy(name,basename2);
						sprintf(num,"%02d",i);
						strcat(name,num);
						strcat(name,ext);
						show_map(binmaps[i],size,name);
					}
				}
				// SHOW DUNBRACK MAPS:  bin_post
				//				if(postcheck >= 0)
				//				{
				//					char basename2[]="pdf_binpost_";
				//					for(int i=0; i<nr_loop+2; i++)
				//					{
				//						strcpy(name,basename2);
				//						sprintf(num,"%02d",i);
				//						strcat(name,num);
				//						strcat(name,ext);
				//						show_map(binmaps_post[i],size,name);
				//					}
				//				}
			}

			//			// Mon: Consider removal...
			//			if(postcheck >= 0)
			//			{
			//				binmaps_post[0] = map_bin2(maps[0], postcheck, size); // MON
			//				binmaps_post[nr_loop+1] = map_bin2(maps[nr_loop+1], postcheck, size); //MON
			//				for(int i=1;i<nr_loop+1;i++)
			//					binmaps_post[i] = map_bin2(maps[i], postcheck, size); //MON
			//			}

			// Checking FORWARD dihedral_angles. Ramachandran?
			if(output_debug)
			{
				nbad = check_angles(binmaps, dihedral_angle_i, size, nr_atoms_loop, 0, true, NULL, rama_cutoff);
				if(pid==0)
					fprintf(stderr,"rcd> Remark: %d dihedral_angle_i (forward) are not in Ramachandran at start (check starti files)\n",nbad);
				paint_maps_points(binmaps, size, nr_loop+2, "starti", dihedral_angle_i, 0);
			}
		}

		//		if(postcheck >= 0)
		//		{
		//			binmaps_post[0] = map_bin2(maps[0], postcheck, size); // MON
		//			binmaps_post[nr_loop+1] = map_bin2(maps[nr_loop+1], postcheck, size); //MON
		//			for(int i=1;i<nr_loop+1;i++)
		//			{
		//				binmaps_post[i] = map_bin2(maps[i], postcheck, size); //MON
		//			}
		//		}

		// ----------------------------------------------------------------------------------------------------
		//
		//				START LOOP CLOSURE ITERATIONS
		//
		//------------------------------------------------------------------------------------------------------
		int failureloop;
		int trials[2];
		int fails[2];
		double stat_n[2];
		stat_n[0]=0;
		stat_n[1]=0;
		fails[0]=0;
		fails[1]=0;
		trials[0]=0;
		trials[1]=0;
		failureloop = 0;
		save_counter = 0;			/*reset for new PDB loop*/
		k 	   = 0;			/*reset for new PDB loop*/
		totalk 	   = 0;			/*reset for new PDB loop*/
		ibond_r 	   = 0;			/*just some value to make it work*/
		timer1->start();							/*Start slow timer*/

		float rmsd_total_sum = 0.0, rmsd_avg_solution = 0.0, rmsd_solution_sum = 0.0, rmsd_avg_gen = 0.0;
		float ratio = 0.0;
		int counter_per_gen = 0,counter_solution_rmsd = 0;
		float anchor_rmsd_avg; // Current average of anchor RMSD (for linear closure)
		int intra_clashes_avg; // Current average of "intra-clashes"
		int extra_clashes_avg; // Current average of "extra-clashes" (clashes with map)

		double t_time = 0.0;
		double timer_0 = 0.0; // pre-begin
		double timer_1 = 0.0; // begin
		double timer_2 = 0.0; // inside
		double timer_3 = 0.0; // infinite
		double timer_4 = 0.0; // RCD

		double *dihedrals_sol = (double *) malloc(sizeof(double) * cycles * nr_atoms_loop);
		double *valence_sol = (double *) malloc(sizeof(double) * cycles * nr_atoms_loop);
		double *length_sol = (double *) malloc(sizeof(double) * cycles * nr_atoms_loop);

		Matrix *co_best; // Best "co" for energy filtration
		Matrix *p1co_best; // Best "p1co" for energy filtration
		double *dihedrals_best; // Best dihedrals for energy filtration
		double *valence_best; // Best valence angles for energy filtration
		double *length_best; // Best bond lengths for energy filtration

		self_clash_count=0;
		self_clash_count_failed=0;
		n_infinite=0;
		n_max_grid_clash = 0;
		if(bump_filter)
			presampling_continue = false; // true, to continue standard iterations...
		bump_cutoff = 999999; // Bump filtering cutoff sigma

#ifdef MPI_RELEASE
		// Mon: Bug 6/4/2016 (--bump_filter option could not be removed...)
		// After allocating "full cycles" memory for results, the master process only should do the corresponding "cycles/np" work...
		// ... unless it has to perform "presampling"
		if(pid==0 && !presampling_switch) // If Master
			cycles = cycles0 / np + cycles0 % np;

		if(pid!=0 && presampling_switch) // If Worker
		{
			MPI_Bcast(&bump_cutoff, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
			if(!server)
				fprintf(stderr,"Worker %d received \"bump_cutoff\" (%f) from Master...\n",pid,bump_cutoff);
			presampling_continue = true;
		}
#endif

		anchor_rmsd_avg = 0.0;
		intra_clashes_avg = 0;
		extra_clashes_avg = 0; // Current average of "extra-clashes" (clashes with map)

		while(k < cycles)
		{

			timer1->fstart(); // Mon: Start Pieter's fast timer

			counter_cycles_single 	 = 		0;
			total_counter		 	 = 		0;
			bond_previous		 	 =		0;
			//:TODO:	check the relevance of this !

			if (k == number_path)
				counter_path = 0; //	reset counter_path to zero for new protein loop to be saved

			switch (switching_mode)
			{
			case 0: // FORWARD
				reverse = 0;
				break;
			case 1: // REVERSE
				reverse = 1;
				break;
			case 2: // RANDOM
				reverse = rg->IRandom(0,1);
				break; //	CHANGING PROTOCOL REVERSE-FORWARD 50-50 case
			case 3: // CONMUTE
				if (failureloop == 1)
				{
					if (reverse == 0)
						reverse = 1;
					else
						reverse = 0;
					failureloop = 0;
				}
				else
				{
					if (k<50)
						reverse = rg->IRandom(0,1); //	CHANGING PROTOCOL REVERSE-FORWARD 50-50 case to get stats
					else
						if (rg->Random() < 1-(fails[0]/(trials[0]+1.0)+0.01)/(fails[1]/(trials[1]+1.0)+fails[0]/(trials[0]+1.0)+0.02) )
							reverse=0;
						else
							reverse=1;
				}
				break;
			case 4: // Aggressive forward/backward choice (fastest)
				if(k < (cycles/50) )
					reverse = k%2; // 50%/50% always!
				else
					if (rg->Random() < (trials[0]-fails[0])/(k+1.0) )
						reverse=0;
					else
						reverse=1;
				break;
			case 5:
				// Switching between forward and backward (50%/50% always)
				if (reverse == 0)
					reverse = 1;
				else
					reverse = 0;
				break;
			default:
				if(pid==0)
					fprintf(stderr,"Unknown option...forcing exit!\n");
				exit(1);
				break;
			}

			if (reverse == 0)
				trials[0]++;
			else
				trials[1]++;

			//REINITIALIZE LOOP
			if (reverse == 1)  // backward
			{
				for (int z = 0; z < nr_atoms_loop; z++)
				{
					// MON: consider removal "dihedral_angle", warning but anchor dihedrals must be ok here!
					dihedral_angle[z]   =   dihedral_angle_r[z];  /*Initialize dihedrals reversely*/
					marker[z] 			=	marker_r[z];      /*Reinitialization markers reversely*/
					marker2[z]			= 	marker2_r[z];
					psimarker[z]		= 	psimarker_r[z];				//:TODO:	HOW ABOUT THE NATIVE PROTOCOL BELOW
				}

				for (int i = 0; i < 3; i++)								/*Reinitialize backward anchor*/
					for(int j = 0; j < 3; j++)
						anchor_end->el[i][j] = anchor_end_r->el[i][j];

			}
			else // forward
			{
				for (int z = 0; z < nr_atoms_loop; z++)
				{
					dihedral_angle[z] = dihedral_angle_f[z]; 	  /*Initialize dihedrals*/
					marker[z] 			=	marker_i[z];      /*Reinitialize markers*/
					marker2[z]			= 	marker2_i[z];
					psimarker[z]	    = 	psimarker_i[z];
				}

				for (int i = 0; i < 3; i++)						/*Reinitialize forward anchor*/
					for(int j = 0; j < 3; j++)
						anchor_end->el[i][j] = anchor_end_i->el[i][j];

			}

			//RANDOMIZE THE INTERNAL COORDINATES
			if (randomize == 1)
			{
				//DEFAULT PEPTIDE GEOMETRY FOR NO_LOOP CASE
				randomizeIC(bond_length,valence_angle,bond_length_back,valence_angle_back,nr_atoms_loop,bondptb,valptb,seq,seq_ct1[0],flat_switch);

				// INITIALIZE SOME ANCHOR INTERNAL COORDINATES
				for (int i = 0; i < nr_atoms_loop-1; i++)		// native loop initialization
				{
					prearray1[i] = cos( (M_PI - valence_angle[i+1])/2.0 );
					prearray2[i] = sin( (M_PI - valence_angle[i+1])/2.0 );
				}
				prearray1[nr_atoms_loop-1] = cos( (M_PI - valence_angle[nr_atoms_loop-3])/2.0 );
				prearray2[nr_atoms_loop-1] = sin( (M_PI - valence_angle[nr_atoms_loop-3])/2.0 );


			}

			// RANDOMIZE OMEGA ANGLE
			if(omega_sigma > 0)
				for(int i = 3; i < nr_atoms_loop-1; i+=3 ) // i+1 --> OMEGA
				{
					if(pdbfile[i_pdb].sequence[i/3]=='p') // If cis-Proline
						dihedral_angle[i+1]	=  (-0.0 + omega_sigma*normRand(flat_switch)) *M_PI/180.0;
					else
					{
						dihedral_angle[i+1]	=  (-omega + omega_sigma*normRand(flat_switch)) *M_PI/180.0;
//						if (h3_switch) {
//									        if ( (i/3 + h3NtShift) < 4) dihedral_angle_i[i+1] =  (Ctomega[(i/3 + h3NtShift)-1]+ omega_sigma*normRand(flat_switch)) *M_PI/180.0;
//										    	else if ((i/3 + h3CtShift +1) >=nr_loop-2)   dihedral_angle_i[i+1]  = (Ntomega[(i/3 -nr_loop +1 + h3NtShift)]+ omega_sigma*normRand(flat_switch))*M_PI/180.0;
//										}

					}
				correctdihedral(&dihedral_angle[i + 1]);


				}

			// end randomize protocol



			rmsd = rmsd_crit - 0.01;	// Reinitialize RMSD < RMSD critical
			timer1->fstart();			// Start fast timer

			// ***************************************************
			// BEGIN LOOP (INITIALIZE DIHEDRALS)
			// ***************************************************
			do {
				check2=0;
				counter_per_gen++; // Counter of trials (current PDB)
				total_gen++; // Counter of trials (all PDBs)

				if (ramachandran_check_aux == 0) // NO rama check
				{
					for(int i = 2; i < nr_atoms_loop-1; i++ )	/*randomize 4x all torsions*/
						if(marker[i] == 1)
						{
							angle_r = 2.0*M_PI*(rg->Random() -0.5); // from -PI to +PI radians
							dihedral_angle[i + 1] +=  angle_r; // Mon: Check this! shit...
							correctdihedral(&dihedral_angle[i + 1]);
						}
				}
				else if (ramachandran_check_aux == 1) /*check for Ramachandran validity*/
				{
					// RANDOMIZE DIHEDRAL
					for(int i = 2; i < nr_atoms_loop-1; i++ )
					{
						if(i%3 != 0) // if i+1 is not OMEGA, i.e. it is Phi or Psi
						{
							// GET PHI AND PSI WITHIN RAMACHANDRAN LIMITS
							index_rama 	=rg->IRandom(0,size_rama);
							if(marker2[i] == 1)		//:INFO: non-glycine residue
							{
								if(psimarker[i] == 0) // PHI
								{
									angle_r 			=	angle_phi[index_rama] - dihedral_angle[i+1];
									dihedral_angle[i+1]	=	angle_phi[index_rama];
								}
								else if(psimarker[i] == 1) // PSI
								{
									angle_r				=	angle_psi[index_rama] - dihedral_angle[i+1];
									dihedral_angle[i+1]	=	angle_psi[index_rama];
								}
							}
							else if (marker2[i] == 2)	//glycine residue pablo
							{
								angle_r = 2.0*M_PI*(rg->Random() -0.5);
								dihedral_angle[i + 1] +=  angle_r;	/*the correct and logical code*/ /*update angle*/
							}
							correctdihedral(&dihedral_angle[i+1]);
						}
					}
				}
				else if (ramachandran_check_aux == 2) // Dunbrack ramachandran maps
				{
					int j;
					float phi,psi;

					// RANDOMIZE DIHEDRAL (atom-fashion...)
					//
					// Randomize phi or psi of First anchor
					if(reverse==0) // forward
					{
						if(marker[2] == 1) // mobile First Psi (forward)
						{
							dihedral_angle[3] = rang_index1D_dun(accs[0], size) + 0.087*(rg->Random() - 0.5); // Psi, radians (+- 2.5 degrees)
							correctdihedral(&dihedral_angle[3]);
						}
					}
					else // backward
					{
						if(marker[nr_atoms_loop-2] == 1) // mobile First Psi (backward)
						{
							dihedral_angle[nr_atoms_loop-1] = rang_index1D_dun(accs[0], size) + 0.087*(rg->Random() - 0.5); // Psi, radians (+- 2.5 degrees)
							correctdihedral(&dihedral_angle[nr_atoms_loop-1]);
						}
					}

					// Randomize phi or psi of Last anchor
					if(reverse==0) // forward
					{
						if(marker[nr_atoms_loop-2] == 1) // mobile Last Phi (forward)
						{
							dihedral_angle[nr_atoms_loop-1] = rang_index1D_dun(accs[nr_loop+1], size) + 0.087*(rg->Random() - 0.5); // Psi, radians (+- 2.5 degrees)
							correctdihedral(&dihedral_angle[nr_atoms_loop-1]);
						}
					}
					else // backward
					{
						if(marker[2] == 1) // mobile Last Phi (backward)
						{
							dihedral_angle[3] = rang_index1D_dun(accs[nr_loop+1], size) + 0.087*(rg->Random() - 0.5); // Psi, radians (+- 2.5 degrees)
							correctdihedral(&dihedral_angle[3]);
						}
					}

					for(int i = 3; i < nr_atoms_loop-3; i++ )
					{
						if(i%3==0) // if OMEGA
						{
							if (reverse == 0)
								rang_index2D_dun(accs[i/3],size,step,&phi,&psi); // Randomize phi/psi in all AA with Dunbrack values
							else
							{
								j = (nr_loop+2) - ((i/3)+1); // check this...
								rang_index2D_dun(accs[j],size,step,&phi,&psi);
							}
						}
						else // PHI or PSI (mobile)
						{
							if(psimarker[i] == 0) // PHI
							{
								if(reverse == 0) // forward
								{
									if(marker[i+1] == 1) // if the corresponding PSI is mobile
										angle_r = (double)phi; // moving PHI
								}
								else // backward
									if(marker[i-1] == 1) // if the corresponding PSI is mobile
										angle_r = (double)phi; // moving PHI
							}
							else // PSI
							{
								if(reverse == 0) // forward
								{
									if(marker[i-1] == 1) // if the corresponding PSI is mobile
										angle_r = (double)psi; // moving PSI
								}
								else // backward
									if(marker[i+1] == 1) // if the corresponding PSI is mobile
										angle_r = (double)psi; // moving PSI
							}

							// GET PHI AND PSI WITHIN RAMACHANDRAN LIMITS
							dihedral_angle[i+1]	= angle_r  + 0.087*(rg->Random() - 0.5);
							correctdihedral(&dihedral_angle[i + 1]);
						}
					}
				}

				// fprintf(stderr,"BEGIN Number of bad angles: %d\n",check_angles(binmaps, dihedral_angle, size, nr_atoms_loop, reverse, true, marker, rama_cutoff));
				// if(reverse==0)
				//  fprintf(stderr,"reverse=%d  First Psi: %lf    Last Phi: %lf\n", reverse, dihedral_angle[3], dihedral_angle[nr_atoms_loop-1]);
				// else
				//  fprintf(stderr,"reverse=%d  First Psi: %lf    Last Phi: %lf\n", reverse, dihedral_angle[nr_atoms_loop-1], dihedral_angle[3]);

				//				dihedral_angle_i[2] = FirstPhi;
				//				dihedral_angle_i[0] = LastPsi;

				// MON: REMOVE THIS! TESTING SOMETHING
				//				if(reverse == 0) // forward
				//					dihedral_angle[2] = dihedral_angle_i[2];
				//				else // backward
				//					dihedral_angle[0] = dihedral_angle_i[0];
				//				dihedral_angle[3] = dihedral_angle_i[3];
				//				dihedral_angle[4] = dihedral_angle_i[4];
				////				dihedral_angle[5] = dihedral_angle_i[5];
				////				dihedral_angle[6] = dihedral_angle_i[6];
				// MON: END TESTING SOMETHING
				// MON: Put here the "randomize" if isarotatetree always enabled
				co->copy(0,0,&co_i->el[0],nr_atoms_loop,3);								// The 3 first atoms are required in "co" matrix

				if(reverse==0) // Forward
				{
					//FORWARD COORDINATES ARE NOW CREATED
					isafrotatetreef_cg(co->el,p1co->el,nr_atoms_loop,bond_length,dihedral_angle,prearray3,prearray4,prearray5,prearray6,prearray1,prearray2,seq,full,cg_mode,factor_cacb);
				}
				else
				{
					// backward --> to forward
					for (int z = 1; z < nr_atoms_loop-2; z++)	// Initialize reverse dihedrals
						dihedral_angle2[z+2]	=	dihedral_angle[nr_atoms_loop-z];
					dihedral_angle2[2] = FirstPhi; // forward
					dihedral_angle2[0] = LastPsi; // forward

					isafrotatetreef_cg(co->el,p1co->el,nr_atoms_loop,bond_length,dihedral_angle2,prearray3,prearray4,prearray5,prearray6,prearray1,prearray2,seq,full,cg_mode,factor_cacb);

					for (int x = 0; x < 3; x++)
						shift_rev[x] 	= 	co->el[nr_atoms_loop-1][x] - co_i->el[nr_atoms_loop-1][x];


					int start=0,end=nr_atoms_loop-1;
					while(start<=end)
					{
						for (int x = 0; x < 3; x++)  {
							rdump=co->el[start][x]-shift_rev[x];
							co->el[start][x] = co->el[end][x]-shift_rev[x];
							co->el[end][x]=rdump;
							if(cg_mode>0)
							{
								rdump=p1co->el[start][x]-shift_rev[x];
								p1co->el[start][x] = p1co->el[end][x]-shift_rev[x];
								p1co->el[end][x]=rdump;
							}
						}
						start++;
						end--;
					}

					//ROTATE REVERSE ANCHOR IN CORRECT PLACE
					for (int z = 0; z < 3; z++)
						for (int x = 0; x < 3; x++) {
							coord2[z][x] = co_i->el[nr_atoms_loop - 1 - z][x]; // Get Ct-anchor for alignment...
						}

					pose(co->el,coord2,posespinor2);
					quatrotate_cg(co->el,p1co->el,p2co->el,nr_atoms_loop,posespinor2,cg_mode);

					// old stuf.... PABLO 2021
					//					for (int z = 0; z < nr_atoms_loop; z++)
					//					{
					//						for (int x = 0; x < 3; x++)
					//						{
					//							co_r->el[z][x] 	= 	  co->el[nr_atoms_loop-1-z][x] - shift_rev[x]; // co->el[nr_atoms_loop-1][x] + co_i->el[nr_atoms_loop-1][x];
					//							if(cg_mode>0)
					//								p1co_r->el[z][x] = 	p1co->el[nr_atoms_loop-1-z][x] - shift_rev[x]; // co->el[nr_atoms_loop-1][x] + co_i->el[nr_atoms_loop-1][x] ;
					//						}
					//
					//					}
					//
					//
					//					//ROTATE REVERSE ANCHOR IN CORRECT PLACE
					//					for (int z = 0; z < 3; z++)
					//						for (int x = 0; x < 3; x++) {
					//							coord2[z][x] = co_i->el[nr_atoms_loop - 1 - z][x]; // Get Ct-anchor for alignment...
					//						}
					//
					//					pose(co_r->el,coord2,posespinor2);
					//					quatrotate_cg(co_r->el,p1co_r->el,p2co_r->el,nr_atoms_loop,posespinor2,cg_mode);
					//
					//					copyco(co_r->el,co->el,nr_atoms_loop); // Copy reverse coordinates "co_r" into "co"
					//					if(cg_mode>0)
					//						copyco(p1co_r->el,p1co->el,nr_atoms_loop);

					// multiPDB_cg(co->el,p1co_r->el,1,nr_loop+2,"reverseTF.pdb",secu,1,cg_mode,in_chain);
					//exit(1);

				}
				check2 = clash_cg(vol,co->el,p1co->el,nr_atoms_loop,marker2,reverse,cg_mode);
				if(check2 == 0 && self_clash)
					check2 = intra_clash_cg(co->el,p1co->el,factor_co,factor_ocb,nr_atoms_loop,atomshift,marker2,reverse, cg_mode);

				rmsd = rmsd_anchor3(co->el+last,anchor_end->el);
				rmsd_pdbs_gen += rmsd;
				rmsd_total_sum += rmsd; // Summatory of trials
			} while( check2 == 1 || rmsd < rmsd_crit || rmsd > rmsd_thr*(nr_loop+2) ); // This while is checking proposes only...
			// ********************************************************
			// END RANDOMIZE NATIVE LOOP (generated structures are OK)
			// ********************************************************

			// ALEX: Counters and accumulators for generation statistics
			total_sols++;
			counter_solution_rmsd++; // Counter of good RMSD
			rmsd_solution_sum += rmsd; // Sumatory of good RMSD
			rmsd_pdbs_sol += rmsd;
			rmsd_avg_solution = rmsd_solution_sum/counter_solution_rmsd; // Avg of solutions
			rmsd_avg_gen = rmsd_total_sum/counter_per_gen; // Avg of trials
			ratio = ( counter_solution_rmsd*1.0 / counter_per_gen )*100.0; // Ratio solutions/trials
			// fprintf(stderr,"rmsd_pdbs_gen = %f rmsd_total_sum = %f rmsd_pdbs_sol = %f rmsd_avg_solution = %f rmsd_avg_gen = %f ratio = %f  \n",rmsd_pdbs_gen, rmsd_pdbs_sol,rmsd_solution_sum,rmsd_avg_solution,rmsd_avg_solution,rmsd_avg_gen,ratio);

			// Storing the generated loops
			for(int i = 0; i < nr_atoms_loop ; i++)
				for(int j = 0; j < 3; j++)
					if (reverse == 1)
					{
						co_soli->el[i   + k*nr_atoms_loop][j]   =    co->el[nr_atoms_loop-1-i][j];
						if(cg_mode>0)
							p1co_soli->el[i + k*nr_atoms_loop][j]   =  p1co->el[nr_atoms_loop-1-i][j];
					}
					else
					{
						co_soli->el[i   + k*nr_atoms_loop][j]   =    co->el[i][j];
						if(cg_mode>0)
							p1co_soli->el[i + k*nr_atoms_loop][j]   =  p1co->el[i][j];
					}

			if(pid==0 && output_debug)
				multiPDB_cg(co->el,p1co->el,1,nr_loop+2,"nofunciona_2.pdb",secu,1,cg_mode,in_chain);

			save_counter++;

			if (save_movie && k == number_path)
			{
				profile_counter 						= 0;
				rmsdprofile->array[profile_counter] 	= rmsd;
				counterprofile->array[profile_counter] 	= counter_cycles_single;
				profile_counter++;
			}

			// START ACTUAL ALGORITHM
			check2			  = 0;	//:TRICKY: must be zero to allow for without-grid approach
			counter_perturbed = 0;
			counter_optimized = 0;
			counter_passed    = 0;
			self_marker		  = 0;
			atom_marker       = 3;	// new code --- must be checked first---initialization: problems?*/ /*anyway must be set at start

			// -----------------------------------------
			//        KERNEL LOOP CLOSURE
			//------------------------------------------
			float rmsd_old=9999;
			int rmsd_count=0;
			int check_count=0;
			int check_bond=0;
			int check_old;

			// *************
			// FIX DIHEDRALS
			// *************
			int aa_fix;
			for (int z =0 ; z < naa_fix ; z++)
			{
				aa_fix=	rg->IRandom(3,nr_atoms_loop/3-3);
				if (reverse == 1)
				{
					marker[aa_fix*3+1]=0; // psi
					marker[aa_fix*3+2]=0; // phi
				}
				else
				{
					marker[nr_atoms_loop-aa_fix*3+1]=0; // psi
					marker[nr_atoms_loop-aa_fix*3+2]=0; // phi
				}
			}

			// ****************************************************
			// INSIDE LOOP
			// ****************************************************
			int exited_rmsd = 0;
			int exited_infinite = 0;
			int exited_check = 0;

			int check2alex = 0;
			int indexphi,indexpsi; // Index in PDF
			int aux_index,opt_index, ind;
			double ang_phi,ang_psi;
			double opt_angle;

			rmsd=300000000;
			ibond = 0;

			while( rmsd > rmsd_crit && total_counter < max_path )
			{
				//  SELECT A VALID BOND
				ibond = bondselect_marker(nr_atoms_loop, marker,&bond_previous,reverse); // SELECT RANDOM BOND
				check_bond=0;

				// OPTIMIZE
				angle = bestangle_anchor_new(co->el,ibond,nr_atoms_loop,anchor_end->el,matrix->el,polpoint->el,anchorpoint->el,result1->el,result2->el,result3->el);

				// RAMACHANDRAN INSIDE-FILTER
				if ( ramachandran_inside == 1 && marker2[ibond] == 1 )
				{
					current = dihedral_angle[ibond+1] + angle;
					correctdihedral(&current);

					if (psimarker[ibond] == 0) // Is PHI
					{
						if( current < LOWERPHI) // -175-57
							angle = (LOWERPHI - dihedral_angle[ibond+1])/1.00; // <-- Mon: this is not consistent...
						else if (current > UPPERPHI + 0.00000 ) // -40+57=+17 *************************
							angle = (UPPERPHI - dihedral_angle[ibond+1])/1.00;
					}
					else if (psimarker[ibond] == 1) // Is PSI
					{
						if ( current < LOWERPSI -  0.00000 )
							angle = (LOWERPSI - dihedral_angle[ibond+1])/1.00;
						else if  (current > UPPERPSI )
							angle = (UPPERPSI - dihedral_angle[ibond+1])/1.00;
					}
				}

				// RAMACHANDRAN DUNBRACK INSIDE-FILTER
				else if (ramachandran_inside_aux == 2)
				{
					if(!inside_rand)
						aux_index = angle2index(dihedral_angle[ibond+1], size); // get the old angle index

					current = dihedral_angle[ibond+1] + angle;
					correctdihedral(&current);

					if(reverse==0) // FORWARD
					{
						if (psimarker[ibond] == 0) // Phi
						{
							indexphi = angle2index(current,size);
							if(ibond == nr_atoms_loop-2)
								indexpsi = angle2index(LastPsi,size);
							else
								indexpsi = angle2index(dihedral_angle[ibond+2],size);

							// Phi --> // WE ROTATE PHI AND FIX PSI
							if (binmaps[ibond/3][indexphi][indexpsi] < rama_cutoff[ibond/3]) // Outside allowed Rama...
							{
								if(inside_rand)
								{
									// accu_1D(binmaps[ibond/3],size,indexpsi,false,acc); // Fixes PHI if isphi=true (PSI is moved) or fixes PSI if isphi=false (PHI is moved)
									// angle = rang_index1D_dun(acc, size) + 0.087*(rg->Random() -0.5); // radians (+- 2.5 degrees)
									// 0.087 * 180 / PI = 4.984737028 (less than 5.0, then OK)
									angle = rang_index1D_dun(accs_psi[ibond/3][indexpsi], size) + 0.087*(rg->Random()-0.5) - dihedral_angle[ibond+1]; // radians (+- 2.5 degrees)
								}
								else
									angle = getoptphi(indexphi, indexpsi, aux_index, binmaps, dihedral_angle, ibond/3, ibond, size, &opt_index, rama_cutoff[ibond/3]);
							}
						}
						else // Psi
						{
							indexpsi = angle2index(current,size);

							if(ibond == 2)
								indexphi = angle2index(FirstPhi,size);
							else
								indexphi = angle2index(dihedral_angle[ibond],size);

							// Psi --> WE ROTATE PSI AND FIX PHI
							if (binmaps[ibond/3][indexphi][indexpsi] < rama_cutoff[ibond/3])
								if(inside_rand)
									angle = rang_index1D_dun(accs_phi[ibond/3][indexphi], size) + 0.087*(rg->Random() -0.5) - dihedral_angle[ibond+1]; // radians (+- 2.5 degrees)
								else
									angle = getoptpsi(indexphi, indexpsi, aux_index, binmaps, dihedral_angle, ibond/3, ibond, size, &opt_index, rama_cutoff[ibond/3]);
						}
					}
					else // REVERSE == 1
					{
						int j = ((nr_atoms_loop/3) - (ibond/3) - 1);

						if (psimarker[ibond] == 0) // Phi
						{
							indexphi = angle2index(current,size); // Bug 2 !!!

							if(ibond == 2)
								indexpsi = angle2index(LastPsi,size);
							else
								indexpsi = angle2index(dihedral_angle[ibond],size);

							// Phi --> // WE ROTATE PHI AND FIX PSI
							if (binmaps[j][indexphi][indexpsi] < rama_cutoff[j])
								if(inside_rand)
									angle = rang_index1D_dun(accs_psi[j][indexpsi], size)+0.087*(rg->Random() -0.5) - dihedral_angle[ibond+1]; // radians (+- 2.5 degrees)
								else
									angle = getoptphi(indexphi, indexpsi, aux_index, binmaps, dihedral_angle, j, ibond, size, &opt_index, rama_cutoff[j]);
						}
						else // Psi
						{
							indexpsi = angle2index(current,size); // BUG 1 !!!

							if(ibond == nr_atoms_loop-2)
								indexphi = angle2index(FirstPhi,size);
							else
								indexphi = angle2index(dihedral_angle[ibond+2],size);

							// Psi --> WE ROTATE PSI AND FIX PHI
							if (binmaps[j][indexphi][indexpsi] < rama_cutoff[j])
								if(inside_rand)
									angle = rang_index1D_dun(accs_phi[j][indexphi], size)+0.087*(rg->Random() -0.5) - dihedral_angle[ibond+1]; // radians (+- 2.5 degrees)
								else
									angle = getoptpsi(indexphi, indexpsi, aux_index, binmaps, dihedral_angle, j, ibond, size, &opt_index, rama_cutoff[j]);
						}
					}
					correctdihedral(&angle);
				}  // end ramachandrans...

				//
				//  ROTATION and COLLISION TESTS
				//

				// Current structure backup
				copyco(co->el,co2->el,ibond,nr_atoms_loop);     /*make copy backbone before rotation*/
				if(cg_mode>0)
					copyco(p1co->el,p1co2->el,ibond,nr_atoms_loop);	/*make copy backbone before rotation*/

				// slow
				// check2 =rotatetree_cg_check(co->el,p1co->el,ibond,nr_atoms_loop,angle,reverse,cg_mode,vol,posV, stepV, unitV, limitV, marker2);


				// Single bond ROTATION
				rotatetree_cg(co->el,p1co->el,ibond,nr_atoms_loop,angle,reverse,cg_mode);

				// CLASH TESTs:
				// LOOP vs. MAP TEST
				// check2 = clash_cg(vol,co->el,p1co->el,nr_atoms_loop,marker2,reverse,cg_mode);
				check2 = clash_cg_ibond(vol,co->el,p1co->el,ibond,nr_atoms_loop,marker2,reverse,cg_mode);

				// LOOP vs. LOOP (un-comment this to slow down 2x)
				//if(self_clash && check2 == 0)
				//check2 = intra_clash_cg(co->el,p1co->el,factor_co,factor_ocb,nr_atoms_loop,atomshift,marker2,reverse,0);
				// check2 = intra_clash_cg_ibond(co->el,p1co->el,factor_co,factor_ocb,ibond,nr_atoms_loop,atomshift,marker2,reverse,0);

				if( check2 == 0) // IT IS OK (NO CLASH)
				{
					rmsd_old = rmsd; // MON: this would be removed...
					rmsd = rmsd_anchor3(co->el+last,anchor_end->el);		//:TRICKY:	also variant using three atoms -> rmsd_anchor3(double **, double **)
					if(rmsd  < rmsd_old * pass_factor)
					{
						// Accept conformation
						dihedral_angle[ibond + 1] += angle;
						correctdihedral(&dihedral_angle[ibond + 1]);
						counter_passed++;
						check_count=0;
					}
					else
					{
						// Reject conformation
						copyco(co2->el,co->el,ibond,nr_atoms_loop); // paste the original conformation
						if(cg_mode>0)
							copyco(p1co2->el,p1co->el,ibond,nr_atoms_loop);
						rmsd=rmsd_old;
					}
				}
				else // Reject conformation
				{
					n_max_grid_clash++;
					copyco(co2->el,co->el,ibond,nr_atoms_loop); // paste the original conformation
					if(cg_mode>0)
						copyco(p1co2->el,p1co->el,ibond,nr_atoms_loop);
					rmsd=rmsd_old; // MON: this would be removed...?
				}

				//  Exit.....
				if( (rmsd_old-rmsd < kick_thr) && (rmsd > kick_cut) ) // If stuck at some minima... (or clash...)
					rmsd_count++;  // Increases if both clash or worse rmsd
				else
					rmsd_count=0;

				// The "Kick" ...
				if( enable_kick && (rmsd_count > nr_atoms_loop*kick_factor) && (rmsd > rmsd_crit) )
				{
					//fprintf(stderr,"rmsd_count= %d  kick= %f\n",rmsd_count, nr_atoms_loop*kick_factor);
					n_rand_repairs++; // Number of "repair trials"

					check2=1;
					int check_rep_count=0;
					float phi=0,psi=0;
					while(check2 == 1) // Reparation loop (it only exits when either 1) the structure is OK, or 2) if too many clash-checks
					{
						copyco(co->el,co2->el,nr_atoms_loop); // make copy backbone before rotation
						if(cg_mode>0)
							copyco(p1co->el,p1co2->el,nr_atoms_loop); // make copy backbone before rotation

						// select bond
						ibond = bondselect_marker2(nr_atoms_loop,marker,&bond_previous, reverse);

						int aa_ibond; // Index of the aminoacid given the "ibond"
						if(ramachandran_inside_aux == 2) // Get angles "in Ramachandran"
						{

							// aa_ibond is forward
							if(reverse == 0) // forward
								aa_ibond = ibond/3; // aa index
							else // reverse
								aa_ibond = ((nr_atoms_loop/3) -1) - ibond/3;

							// Randomize phi or psi of first anchor
							if(ibond/3 == 0) // First anchor
							{
								if(reverse==0)
									psi = rang_index1D_dun(accs[0], size) - dihedral_angle[ibond+1]; // radians
								else
									phi = rang_index1D_dun(accs[nr_loop+1], size) - dihedral_angle[ibond+1]; // radians
							}
							else if(ibond/3 == nr_loop+1) // Last anchor
							{
								if(reverse==0)
									phi = rang_index1D_dun(accs[nr_loop+1], size) - dihedral_angle[ibond+1]; // radians
								else
									psi = rang_index1D_dun(accs[0], size) - dihedral_angle[ibond+1]; // radians
							}
							else // mobile aminoacids
							{
								rang_index2D_dun(accs[aa_ibond],size,step,&phi,&psi); // Randomize phi/psi in all AA with Dunbrack values
								if(reverse == 0) // forward
								{
									phi -= dihedral_angle[aa_ibond*3+2];
									psi -= dihedral_angle[aa_ibond*3+3];
								}
								else
								{
									psi -= dihedral_angle[(ibond/3)*3+2];
									phi -= dihedral_angle[(ibond/3)*3+3];
								}
							}
						}
						else if (ramachandran_check_aux == 1) /*check for Ramachandran validity*/
						{
							// OJO REVERSE ??? BUG PABLO 2020
							// GET PHI AND PSI WITHIN RAMACHANDRAN LIMITS
							index_rama 	=rg->IRandom(0,size_rama);
							if(marker2[ibond/3] == 1)		//:INFO: non-glycine residue
							{
								if(psimarker[ibond/3] == 0) // PHI
									angle_r =	angle_phi[index_rama] - dihedral_angle[ibond/3+1];
								else if(psimarker[ibond/3] == 1) // PSI
									angle	=	angle_psi[index_rama] - dihedral_angle[ibond/3+1];
							}
							else if (marker2[ibond/3] == 2)	//glycine residue pablo
								angle = 2.0*M_PI*(rg->Random() -0.5);
						}
						else // NO rama check
							angle = 2.0*M_PI*(rg->Random() -0.5);

						if(ramachandran_inside_aux == 2) // Get angles "in Ramachandran"
						{
							// fprintf(stderr,"ibond= %d  aa_ibond= %d  phi= %f  psi= %f\n",ibond,aa_ibond,phi,psi);
							if(aa_ibond == 0) // Nt aminoacid
							{
								rotatetree_cg(co->el,p1co->el,ibond,nr_atoms_loop,psi,reverse,cg_mode);
							}
							else if(aa_ibond == nr_loop+1) // Ct aminoacid
							{
								rotatetree_cg(co->el,p1co->el,ibond,nr_atoms_loop,phi,reverse,cg_mode);
							}
							else
							{
								if(reverse == 0)
								{
									rotatetree_cg(co->el,p1co->el,aa_ibond*3+1,nr_atoms_loop,phi,reverse,cg_mode);
									rotatetree_cg(co->el,p1co->el,aa_ibond*3+2,nr_atoms_loop,psi,reverse,cg_mode);
								}
								else
								{
									rotatetree_cg(co->el,p1co->el,(ibond/3)*3+1,nr_atoms_loop,psi,reverse,cg_mode);
									rotatetree_cg(co->el,p1co->el,(ibond/3)*3+2,nr_atoms_loop,phi,reverse,cg_mode);
								}
							}
						}
						else
							rotatetree_cg(co->el,p1co->el,ibond,nr_atoms_loop,angle,reverse,cg_mode);

						rmsd = rmsd_anchor3(co->el+last,anchor_end->el); // Mon: put below...

						// Checking whether perturbed structure is OK
						check2 = clash_cg(vol,co->el,p1co->el,nr_atoms_loop,marker2,reverse,cg_mode);

						if( (self_clash == 1) && (check2 == 0) )
							check2 = intra_clash_cg(co->el,p1co->el,factor_co,factor_ocb,nr_atoms_loop,atomshift,marker2,reverse,cg_mode);

						if(check2 > 0) // If clashes detected, restoring un-perturbed structure
						{
							check_rep_count++;
							copyco(co2->el,co->el,nr_atoms_loop);		// paste the original conformation
							if(cg_mode>0)
								copyco(p1co2->el,p1co->el,nr_atoms_loop);

							if (check_rep_count>nr_atoms_loop*0.666) // if too many clash-checks...
							{
								total_counter = max_path; // forces "main-loop" exit and full regeneration (restart, rebegin, ...)
								check2=0;
							}
						}
						else // If structure is OK, then update dihedral_angles and continue...
						{
							rmsd_old=rmsd;
							rmsd_count=0;
							if(ramachandran_inside_aux == 2) // Get angles "in Ramachandran"
							{

								if(aa_ibond == 0) // Nt aminoacid
								{
									dihedral_angle[ibond + 1] +=  psi;	// restore dihedral_angle....
									correctdihedral(&dihedral_angle[ibond + 1]);
								}
								else if(aa_ibond == nr_loop+1) // Ct aminoacid
								{
									dihedral_angle[ibond + 1] +=  phi;	// restore dihedral_angle....
									correctdihedral(&dihedral_angle[ibond + 1]);
								}
								else
								{
									if(reverse == 0) // forward
									{
										dihedral_angle[aa_ibond*3+2] +=  phi;	// restore dihedral_angle....
										correctdihedral(&dihedral_angle[aa_ibond*3+2]);
										dihedral_angle[aa_ibond*3+3] +=  psi;	// restore dihedral_angle....
										correctdihedral(&dihedral_angle[aa_ibond*3+3]);
									}
									else
									{
										dihedral_angle[(ibond/3)*3+2] +=  psi;	// restore dihedral_angle....
										correctdihedral(&dihedral_angle[(ibond/3)*3+2]);
										dihedral_angle[(ibond/3)*3+3] +=  phi;	// restore dihedral_angle....
										correctdihedral(&dihedral_angle[(ibond/3)*3+3]);
									}
								}
							}
							else
							{
								dihedral_angle[ibond + 1] +=  angle;	// restore dihedral_angle.... Bug del copon! Pero no vale pa na...
								correctdihedral(&dihedral_angle[ibond + 1]);
							}
						}
					}
				}

				// *******************************************************
				// Here it was the ANALYTIC LOOP-CLOSURE (Coutsias et al.)
				// *******************************************************

				// MOVIE
				if(save_movie && k == number_path)	/*Save RMSD profiling*/
				{
					rmsdprofile->array[profile_counter] 	= (double) rmsd;
					counterprofile->array[profile_counter] 	= counter_cycles_single;
					profile_counter++;
					if (counter_path < max_path)	/*Save closure trajectory one loop*/
					{
						if (reverse == 1)	//NEW: modified for reverse moving...	//:TODO: does it work properly?
						{
							for(int i = 0; i < nr_atoms_loop; i++)
								for(int j = 0; j < 3; j++)
									path->el[i + counter_path*nr_atoms_loop][j] = co->el[nr_atoms_loop-1-i][j];
						}
						else // forward
						{
							for(int i = 0; i < nr_atoms_loop; i++)
								for(int j = 0; j < 3; j++)
									path->el[i + counter_path*nr_atoms_loop][j] = co->el[i][j];
						}
						counter_path += 1;
					}
				}

				total_counter++;
			}	// END OF RMSD LOOP

			counter_cycles_single=total_counter;
			counter_cycles_total+=counter_cycles_single;

			// MON: "dgelsd_" routine must be declared for Intel releases in Sherpa server (add -mkl library)
			//      Alternatively, "linear_closure" stuff can be commented.
			//
			// *****************************************************
			// Linear Loop Closure
			// *****************************************************
			//
			// MON: Just un-comment this to work...(otherwise it will require additional declarations for Intel in Sherpa)
			//
			if( linear_closure && rmsd <= rmsd_crit )
			{
				double *der = NULL; // enables automatic memory allocation
				double *b; // 9x1 array
				double rmsd_zero;

				// Using SVD (more sophisticated: computes the rank of Aop)
				int info=-1;
				int M = 9; // Number of constraints ( 3(x,y,z) x 3(atoms) )
				int N = (nr_loop+2)*2; // Number of variables
				int NRHS = 1; // Number of Right-Hand-Side columns. One solution only...
				int LDA = M;
				int LDB = N;
				int RANK=0;
				int SMLSIZ = 40; // Try different values
				int NLVL = 1 + floor( log(M/(SMLSIZ+1))/log(2.0) );
				// int LWORK = 12*M + 2*M*SMLSIZ + 8*M*NLVL + M*NRHS + pow(SMLSIZ+1,2);
				int LWORK = 20000;
				// int LWORK = -1; // query
				// int LIWORK = 3 * M * NLVL + 11 * M;
				int LIWORK = 20000;
				double *work;
				int *iwork;
				double RCOND = 0.0001;
				double *svals;
				svals = (double *) malloc( sizeof(double) * M);
				work = (double *) malloc( sizeof(double) * LWORK );
				iwork = (int *) malloc( sizeof(int) * LIWORK );
				b = (double *) malloc(sizeof(double) * N); // "b" must be the same size as "X" to allocate output (X)

				double *val;
				val = (double *) malloc( sizeof(double) * M);
				double *der2;
				der2 = (double *) malloc( sizeof(double) * N*M);
				double *b2;
				b2 = (double *) malloc( sizeof(double) * N);

				if(debug_linear)
					fprintf(stderr,"M= %d\nN= %d\nNRHS= %d\nLDA= %d\nLDB= %d\nRANK= %d\nSMLSIZ= %d\nNLVL= %d\nLWORK= %d\nLIWORK= %d\nRCOND= %f\n"
							,M,N,NRHS,LDA,LDB,RANK,SMLSIZ,NLVL,LWORK,LIWORK,RCOND);

				for(int nn=0; nn<linear_iters; nn++)
				{
					if(debug_linear)
						fprintf(stderr, "LOOP %3d  ITER %2d\n", k, nn);

					// Compute the derivatives matrix ("dr/d_theta") of the positions (r) of the 3 anchor atoms (N,CA,C), i.e. 9 constraints,
					// w.r.t. all dihedrals (theta)
					der = NULL; // DGELSD destroys "der" array... (it is freed...)
					der_anchor(&der, co->el, nr_loop+2, reverse);

					for(int x=0; x<N; x++)
						b[x] = 0.0; // Reset displacement vector B for DGELSD

					// Compute the displacement vector B between the mobile Ct-anchor and the fixed one.
					if(nn==0)
						rmsd_zero = disp_b(b, co->el, anchor_end->el, nr_loop+2);
					else
						disp_b(b, co->el, anchor_end->el, nr_loop+2);

					//exit(0);

					//  SUBROUTINE DGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,
					// $                   WORK, LWORK, IWORK, INFO )
					//*
					//*  -- LAPACK driver routine (version 3.1) --
					//*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
					//*     November 2006
					//*
					//*     .. Scalar Arguments ..
					//  INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
					//  DOUBLE PRECISION   RCOND
					//*     ..
					//*     .. Array Arguments ..
					//  INTEGER            IWORK( * )
					//  DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), S( * ), WORK( * )
					//*     ..
					//*
					//*  Purpose
					//*  =======
					//*
					//*  DGELSD computes the minimum-norm solution to a real linear least
					//*  squares problem:       minimize 2-norm(| b - A*x |)
					//*  using the singular value decomposition (SVD) of A. A is an M-by-N
					//*  matrix which may be rank-deficient.
					//*
					//*  Several right hand side vectors b and solution vectors x can be
					//*  handled in a single call; they are stored as the columns of the
					//*  M-by-NRHS right hand side matrix B and the N-by-NRHS solution
					//*  matrix X.
					//*
					//*  The problem is solved in three steps:
					//*  (1) Reduce the coefficient matrix A to bidiagonal form with
					//*      Householder transformations, reducing the original problem
					//*      into a "bidiagonal least squares problem" (BLS)
					//*  (2) Solve the BLS using a divide and conquer approach.
					//*  (3) Apply back all the Householder tranformations to solve
					//*      the original least squares problem.
					//*
					//*  The effective rank of A is determined by treating as zero those
					//*  singular values which are less than RCOND times the largest singular
					//*  value.
					//*
					//*  The divide and conquer algorithm makes very mild assumptions about
					//*  floating point arithmetic. It will work on machines with a guard
					//*  digit in add/subtract, or on those binary machines without guard
					//*  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
					//*  Cray-2. It could conceivably fail on hexadecimal or decimal machines
					//*  without guard digits, but we know of none.
					//*
					//*  Arguments
					//*  =========
					//*
					//*  M       (input) INTEGER
					//*          The number of rows of A. M >= 0. (Mon: Constraints)
					//*
					//*  N       (input) INTEGER
					//*          The number of columns of A. N >= 0. (Mon: Variables)
					//*
					//*  NRHS    (input) INTEGER
					//*          The number of right hand sides, i.e., the number of columns
					//*          of the matrices B and X. NRHS >= 0.
					//*
					//*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
					//*          On entry, the M-by-N matrix A.
					//*          On exit, A has been destroyed.
					//*
					//*  LDA     (input) INTEGER
					//*          The leading dimension of the array A.  LDA >= max(1,M).
					//*
					//*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
					//*          On entry, the M-by-NRHS right hand side matrix B.
					//*          On exit, B is overwritten by the N-by-NRHS solution
					//*          matrix X.  If m >= n and RANK = n, the residual
					//*          sum-of-squares for the solution in the i-th column is given
					//*          by the sum of squares of elements n+1:m in that column.
					//*
					//*  LDB     (input) INTEGER
					//*          The leading dimension of the array B. LDB >= max(1,max(M,N)).
					//*
					//*  S       (output) DOUBLE PRECISION array, dimension (min(M,N))
					//*          The singular values of A in decreasing order.
					//*          The condition number of A in the 2-norm = S(1)/S(min(m,n)).
					//*
					//*  RCOND   (input) DOUBLE PRECISION
					//*          RCOND is used to determine the effective rank of A.
					//*          Singular values S(i) <= RCOND*S(1) are treated as zero.
					//*          If RCOND < 0, machine precision is used instead.
					//*
					//*  RANK    (output) INTEGER
					//*          The effective rank of A, i.e., the number of singular values
					//*          which are greater than RCOND*S(1).
					//*
					//*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
					//*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
					//*
					//*  LWORK   (input) INTEGER
					//*          The dimension of the array WORK. LWORK must be at least 1.
					//*          The exact minimum amount of workspace needed depends on M,
					//*          N and NRHS. As long as LWORK is at least
					//*              12*N + 2*N*SMLSIZ + 8*N*NLVL + N*NRHS + (SMLSIZ+1)**2,
					//*          if M is greater than or equal to N or
					//*              12*M + 2*M*SMLSIZ + 8*M*NLVL + M*NRHS + (SMLSIZ+1)**2,
					//*          if M is less than N, the code will execute correctly.
					//*          SMLSIZ is returned by ILAENV and is equal to the maximum
					//*          size of the subproblems at the bottom of the computation
					//*          tree (usually about 25), and
					//*             NLVL = MAX( 0, INT( LOG_2( MIN( M,N )/(SMLSIZ+1) ) ) + 1 )
					//*          For good performance, LWORK should generally be larger.
					//*
					//*          If LWORK = -1, then a workspace query is assumed; the routine
					//*          only calculates the optimal size of the WORK array, returns
					//*          this value as the first entry of the WORK array, and no error
					//*          message related to LWORK is issued by XERBLA.
					//*
					//*  IWORK   (workspace) INTEGER array, dimension (MAX(1,LIWORK))
					//*          LIWORK >= 3 * MINMN * NLVL + 11 * MINMN,
					//*          where MINMN = MIN( M,N ).
					//*
					//*  INFO    (output) INTEGER
					//*          = 0:  successful exit
					//*          < 0:  if INFO = -i, the i-th argument had an illegal value.
					//*          > 0:  the algorithm for computing the SVD failed to converge;
					//*                if INFO = i, i off-diagonal elements of an intermediate
					//*                bidiagonal form did not converge to zero.
					//*
					//*  Further Details
					//*  ===============
					//*
					//*  Based on contributions by
					//*     Ming Gu and Ren-Cang Li, Computer Science Division, University of
					//*       California at Berkeley, USA
					//*     Osni Marques, LBNL/NERSC, USA
					//*
					//*  =====================================================================


					// DGELSD EXAMPLE from: http://www.nag.com/lapack-ex/node46.html
					//				double A[6*5] = { // Row-major
					//				 -0.09,  -1.56,  -1.48,  -1.09,   0.08,  -1.59,
					//				  0.14,   0.20,  -0.43,   0.84,   0.55,  -0.72,
					//				 -0.46,   0.29,   0.89,   0.77,  -1.13,   1.06,
					//				  0.68,   1.09,  -0.71,   2.11,   0.14,   1.24,
					//				  1.29,   0.51,  -0.96,  -1.27,   1.74,   0.34
					//				};
					//				double A[6*5] = { // Column-major (LAPACK)
					//						-0.09,  0.14, -0.46,  0.68,  1.29,
					//						-1.56,  0.20,  0.29,  1.09,  0.51,
					//						-1.48, -0.43,  0.89, -0.71, -0.96,
					//						-1.09,  0.84,  0.77,  2.11, -1.27,
					//						 0.08,  0.55, -1.13,  0.14,  1.74,
					//						-1.59, -0.72,  1.06,  1.24,  0.34
					//				};
					//				double B[5] = { 7.4, 4.3, -8.1, 1.8, 8.7 };
					//				DGELSD Example Program Results
					//				 Least squares solution:      1.5938    -0.1180    -3.1501     0.1554     2.5529    -1.6730
					//				 Tolerance used to estimate the rank of A:       1.00E-02
					//				 Estimated rank of A:     4
					//				 Singular values of A:    3.9997     2.9962     2.0001     0.9988     0.0025

					//				int info=0;
					//				int M = 5; // Number of constraints
					//				int N = 6; // Number of variables
					//				int NRHS = 1;
					//				int LDA = M;
					//				int LDB = N;
					//				int RANK=0;
					//				int SMLSIZ = 40; // Try different values
					//				int NLVL = 1 + floor( log(M/(SMLSIZ+1))/log(2.0) );
					////				int LWORK = 12*M + 2*M*SMLSIZ + 8*M*NLVL + M*NRHS + pow(SMLSIZ+1,2);
					//				int LWORK = 10000000;
					////				int LWORK = -1; // query
					////				int LIWORK = 3 * M * NLVL + 11 * M;
					//				int LIWORK = 1000000;
					//				double *work;
					//				int *iwork;
					//				double RCOND = 1.00E-02;
					//				double *svals;
					//				svals = (double *) malloc( sizeof(double) * M);
					//				work = (double *) malloc( sizeof(double) * LWORK );
					////				work = (double *) malloc( sizeof(double) * 100000 );
					//				iwork = (int *) malloc( sizeof(int) * LIWORK );


					for(int i=0; i<N*M; i++)
						der2[i] = der[i]; // Backup (A-matrix) --> "A" is destroyed by DGELSD
					for(int i=0; i<M; i++)
						b2[i] = b[i]; // Backup (B-matrix)

					//  SUBROUTINE DGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,
					// $                   WORK, LWORK, IWORK, INFO )
					dgelsd_(&M, &N, &NRHS, der, &LDA, b, &LDB, svals, &RCOND, &RANK, work, &LWORK, iwork, &info);
					free(der); // Despite "der" (A-matrix) is destroyed by DGELSD... it seem good to free it here!

					if(debug_linear)
						fprintf(stderr,"DGELSD OUTPUT: info=%d\n",info);

					if(LWORK<0)
					{
						fprintf(stderr,"WORK(0)= %f\n",work[0]);
						exit(0);
					}

					double mymax=0.0; // maximum angular displacement
					if(info==0)
					{
						if(debug_linear)
						{
							fprintf(stderr, "RANK= %d\n",RANK);
							fprintf(stderr, "SVALS= ");
							for(int i=0; i<M; i++)
								fprintf(stderr, "%8.6f ",svals[i]);
							fprintf(stderr, "\nX(rad)= ");
						}

						double myavg=0.0;
						int myimax=0;
						for(int i=0; i<N; i++)
						{
							if(debug_linear)
								fprintf(stderr, "%8.6f ",b[i]);

							myavg += abs(b[i]);
							//						fprintf(stderr, "%8.6f ",B[i]);
							if(abs(b[i])>mymax)
							{
								mymax = abs(b[i]); // store maximum dihedral increment
								myimax = i; // index of the maximum
							}
						}
						myavg /= N;

						if(debug_linear)
						{
							fprintf(stderr, "\nX(deg)= ");
							for(int i=0; i<N; i++)
								fprintf(stderr, "%8.6f ",b[i]*180/M_PI);
							fprintf(stderr,"\n");
							fprintf(stderr,"X_max(deg)= %8.6f (%d-th variable)\n",mymax*180/M_PI,myimax);
							fprintf(stderr,"X_avg(deg)= %8.6f\n",myavg*180/M_PI);
						}
					}
					else
					{
						fprintf(stderr,"dgelsd_'s error in linear closure, info= %d. Forcing exit!\n",info);
						exit(3);
					}

					// Thus if you have no extra elements in your matrix with M rows and N columns, then:
					//	- In row-major case: row i lays in memory right after row i-1 and thus LDA=N - number of elements in row.
					//	- In column-major case: column i lays in memory right after column i-1 and thus LDA=M - number of elements in column. (LAPACK's style)

					// Compute A*X
					for(int j=0; j<M; j++) // screen cols (constraints)
					{
						val[j]=0.0;
						for(int i=0; i<N; i++) // screen rows (variables)
							val[j] += der2[i*M+j] * b[i]; // A*X
					}

					if(debug_linear)
					{
						fprintf(stderr, "B=     ");
						for(int j=0; j<M; j++) // screen cols
							fprintf(stderr,"%9.6f ",b2[j]); // "b2" because on exit "b" contains the results "X"

						fprintf(stderr, "\nA*X=   ");
						for(int j=0; j<M; j++) // screen cols
							fprintf(stderr,"%9.6f ",val[j]);
					}

					for(int j=0; j<M; j++) // screen cols
						val[j] -= b2[j]; // -= B

					if(debug_linear)
					{
						fprintf(stderr, "\nA*X-B= ");
						for(int j=0; j<M; j++) // screen cols
							fprintf(stderr,"%9.6f ",val[j]);
						fprintf(stderr, "\n");
					}

					double dummy=0.0;
					dummy += sqrt(val[0]*val[0] + val[1]*val[1] + val[2]*val[2]);
					dummy += sqrt(val[3]*val[3] + val[4]*val[4] + val[5]*val[5]);
					dummy += sqrt(val[6]*val[6] + val[7]*val[7] + val[8]*val[8]);
					dummy /= 3;

					if(debug_linear)
						fprintf(stderr, "Final RMSD aprox. (%02d): %f --> final/initial: %f [%]\n", nn, dummy, 100*dummy/rmsd_zero);

					// Preventing large angular displacements
					if( mymax > linear_maxang )
						linear_delta = linear_maxang / mymax; // set new delta to prevent angular displacements greater than linear_maxang

					if(reverse == 1)
						linear_delta *= -1.0;

					for(int n=1; n<=1; n++)
					{
						// Apply angular increments to the "dihedral_angle" array (the same is valid for forward and backward)
						for(int a=0; a<(nr_loop+2); a++) // screen residues
						{
							if(a>0) // First Phi does not move
							{
								dihedral_angle[a*3 + 2] += linear_delta* b[2*a]; // Phi increment
								correctdihedral(&dihedral_angle[a*3 + 2]);
							}
							if(a<nr_loop+1) // Last Psi is at dihedral_angle[0] position
							{
								dihedral_angle[a*3 + 3] += linear_delta* b[2*a+1]; // Psi increment
								correctdihedral(&dihedral_angle[a*3 + 3]);
							}
						}

						// Rebuild the whole loop from scratch (many dihedrals are modified at the same time!)
						co->copy(0,0,&co_i->el[0],nr_atoms_loop,3);	// The 3 first atoms are required in "co" matrix
						if(reverse==0) // Forward
						{
							//FORWARD COORDINATES ARE NOW CREATED
							isafrotatetreef_cg(co->el,p1co->el,nr_atoms_loop,bond_length,dihedral_angle,prearray3,prearray4,prearray5,prearray6,prearray1,prearray2,seq,full,cg_mode,factor_cacb);
						}
						else // Backward
						{

							// backward --> to forward
							for (int z = 1; z < nr_atoms_loop-2; z++)	// Initialize reverse dihedrals
								dihedral_angle2[z+2]	=	dihedral_angle[nr_atoms_loop-z];
							dihedral_angle2[2] = FirstPhi; // forward
							dihedral_angle2[0] = LastPsi; // forward

							isafrotatetreef_cg(co->el,p1co->el,nr_atoms_loop,bond_length,dihedral_angle2,prearray3,prearray4,prearray5,prearray6,prearray1,prearray2,seq,full,cg_mode,factor_cacb);

							for (int x = 0; x < 3; x++)
								shift_rev[x] 	= 	co->el[nr_atoms_loop-1][x] - co_i->el[nr_atoms_loop-1][x];


							int start=0,end=nr_atoms_loop-1;
							while(start<=end)
							{
								for (int x = 0; x < 3; x++)  {
									rdump=co->el[start][x]-shift_rev[x];
									co->el[start][x] = co->el[end][x]-shift_rev[x];
									co->el[end][x]=rdump;
									if(cg_mode>0)
									{
										rdump=p1co->el[start][x]-shift_rev[x];
										p1co->el[start][x] = p1co->el[end][x]-shift_rev[x];
										p1co->el[end][x]=rdump;
									}
								}
								start++;
								end--;
							}

							//ROTATE REVERSE ANCHOR IN CORRECT PLACE
							for (int z = 0; z < 3; z++)
								for (int x = 0; x < 3; x++) {
									coord2[z][x] = co_i->el[nr_atoms_loop - 1 - z][x]; // Get Ct-anchor for alignment...
								}

							pose(co->el,coord2,posespinor2);
							quatrotate_cg(co->el,p1co->el,p2co->el,nr_atoms_loop,posespinor2,cg_mode);

							//							// OLD STUFF PABLO 2020
							//							for (int z = 1; z < nr_atoms_loop-2; z++)	// Initialize reverse dihedrals
							//								dihedral_angle2[z+2]	=	dihedral_angle[nr_atoms_loop-z];
							//							dihedral_angle2[2] = FirstPhi; // forward
							//							dihedral_angle2[0] = LastPsi; // forward
							//
							//							isafrotatetreef_cg_rev(co->el,p1co->el,nr_atoms_loop,bond_length,dihedral_angle2,prearray3,prearray4,prearray5,prearray6,prearray1,prearray2,seq,full,cg_mode,factor_cacb);
							//
							//							//TRANSLATE AND COPY COORDINATES TO REVERSE ARRAY
							//							for (int z = 0; z < nr_atoms_loop; z++)
							//								for (int x = 0; x < 3; x++)
							//								{
							//									co_r->el[z][x] 	= 	  co->el[nr_atoms_loop-1-z][x] 	- co->el[nr_atoms_loop-1][x] + co_i->el[nr_atoms_loop-1][x];
							//									if(cg_mode>0)
							//										p1co_r->el[z][x]	= 	p1co->el[nr_atoms_loop-1-z][x]  - co->el[nr_atoms_loop-1][x] + co_i->el[nr_atoms_loop-1][x] ;
							//								}
							//
							//							//ROTATE REVERSE ANCHOR IN CORRECT PLACE
							//							for (int z = 0; z < 3; z++)
							//								for (int x = 0; x < 3; x++)
							//									coord2[z][x] = co_i->el[nr_atoms_loop - 1 - z][x]; // Get Ct-anchor for alignment...
							//
							//							pose(co_r->el,coord2,posespinor2);
							//							quatrotate_cg(co_r->el,p1co_r->el,p2co_r->el,nr_atoms_loop,posespinor2,cg_mode);
							//
							//							copyco(co_r->el,co->el,nr_atoms_loop); // Copy reverse coordinates "co_r" into "co"
							//							if(cg_mode>0)
							//								copyco(p1co_r->el,p1co->el,nr_atoms_loop);


						}

						rmsd = rmsd_anchor3(co->el+last,anchor_end->el);

						if(debug_linear)
							fprintf(stdout, "True-rotation RMSD (%02d): %f --> final/initial: %f [%]\n",nn,rmsd, 100*rmsd/rmsd_zero);
					}

					if(rmsd < linear_cutoff)
					{
						if(debug_linear)
							fprintf(stderr, "Next loop?  rmsd= %f  dummy= %f\n",rmsd,dummy);
						break;
					}
				}

				free(svals);
				free(work);
				free(iwork);
				free(val);
				free(der2);
				free(b2);
				free(b);
			}
			// Linear closure end

			// Post-RCD (+linear) checking


			// INTRA_CLASH POST CHECK
			//
			check2 = 0;
			if(self_clash)
			{
				check2 = intra_clash_cg(co->el,p1co->el,factor_co,factor_ocb,nr_atoms_loop,atomshift,marker2,reverse,cg_mode);
				if ((check2 == 1) && (rmsd < rmsd_crit))
					self_clash_count++;
				if ((check2 == 1) && (rmsd > rmsd_crit))
					self_clash_count_failed++;
				if(check2 > 0)
				{
					// fprintf(stderr,"intra-clashes= %d\n",check2);
					intra_clashes_avg++; // count the number of times there is an intra-clash
				}
			}

			// Checking whether linear closure structure is OK
			//			if(check2 == 0 && linear_closure)
			if(check2 == 0)
			{
				// check2 = clash_cg(vol,co->el,p1co->el,nr_atoms_loop,marker2,reverse,cg_mode, true);
				check2 = clash_cg(vol,co->el,p1co->el,nr_atoms_loop,marker2,reverse,cg_mode);
				if(check2 > 0)
				{
					// fprintf(stderr,"linear closure map-clashes= %d\n",check2);
					extra_clashes_avg++; // count the number of "extra-clashes" (clashes with map)
				}
			}

			// Additional check for Analytical Loop-Closure...
			// if(check2==0 && enable_lc && rmsd < rmsd_crit )

			// Mon: When Psi-Nt is fixed, Ramachandran checks should be disabled.
			// linear closure add tons of rama outliers PABLO 2020!!
			//		if(!no_ntpsi && ramachandran_inside_aux == 2 && check2==0 && rmsd < rmsd_crit )
			if(!linear_closure && ramachandran_inside_aux == 2 && check2==0 && rmsd < rmsd_crit )
				//			if( ramachandran_inside_aux == 2 && check2==0 && rmsd < rmsd_crit )
			{
				check2 = check_angles(binmaps, dihedral_angle, size, nr_atoms_loop, reverse, false, marker, rama_cutoff);
				if( check2 == 1 )
					n_rama_fail++;

			}

			if( check2==0 && rmsd < rmsd_crit && ramachandran_inside_aux > 2 && load_ramap )
			{
				int curr_rama_fail = 0;
				if( (curr_rama_fail = check_angles(binmaps, dihedral_angle, size, nr_atoms_loop, reverse, true, marker, rama_cutoff)) != 0 )
				{
					// getchar();
					n_rama_fail++;
					avg_rama_fail += curr_rama_fail;
				}
				//if (curr_rama_fail>0) fprintf( stderr, "Warning %d ramafails\n", curr_rama_fail);

			}

			// PD2 bump filter (MacDonald et al. 2013)
			if(bump_filter && check2==0)
			{
				if (reverse == 1)
				{
					for(int i=0; i<nr_atoms_loop; i++)
					{
						co_dummy->el[i][0] = co->el[nr_atoms_loop-1-i][0];
						co_dummy->el[i][1] = co->el[nr_atoms_loop-1-i][1];
						co_dummy->el[i][2] = co->el[nr_atoms_loop-1-i][2];
					}
					if(cg_mode>0)
						for(int i=0; i<nr_atoms_loop; i++)
						{
							p1co_dummy->el[i][0] = p1co->el[nr_atoms_loop-1-i][0];
							p1co_dummy->el[i][1] = p1co->el[nr_atoms_loop-1-i][1];
							p1co_dummy->el[i][2] = p1co->el[nr_atoms_loop-1-i][2];
						}
					energy_extra = pd2_extrabump(boxcoord, boxtype, boxres, npocket, co_dummy->el, p1co_dummy->el, nr_atoms_loop, loopindex, cg_mode, residuemarker);
					energy_intra = pd2_intrabump(co_dummy->el, p1co_dummy->el, nr_atoms_loop, cg_mode, residuemarker);
				}
				else
				{
					energy_extra = pd2_extrabump(boxcoord, boxtype, boxres, npocket, co->el, p1co->el, nr_atoms_loop, loopindex, cg_mode, residuemarker);
					energy_intra = pd2_intrabump(co->el, p1co->el, nr_atoms_loop, cg_mode, residuemarker);
				}
				energy = energy_extra + energy_intra;
				if( energy > bump_cutoff)
					check2 = 1; // clashed...
			}

			timer1->fstop();			         //	STOP FAST TIMER
			timer1->felapsed();			        //	calculate elapsed time
			ftimer->array[k] = timer1->fcorrect();	//	correct for negative timings
			fcounter->array[k] = total_counter;													//	STORAGE full counter each cycle

			// Checking solution for convergence and incrementing counters
			totalk++;		// add total counter for converged loop
			if( total_counter >= max_path || check2 > 0 )
			{
				// fprintf(stderr,"Failed k %d rmsd %f total_counter %d Check2 %d \n",k, rmsd, total_counter, check2);
				failureloop = 1;
				fails[reverse]++;
				if(check2 == 0)
					n_max_path++; // Count the number of times that the maximum-iterations-per-loop are exceeded
			}
			else // Solution is OK !!!
			{
				anchor_rmsd_avg += rmsd; // Store the final RMSD for statistics

				// SAVE LOOP in SOLUTION array
				stat_n[reverse] += total_counter;
				for(int i = 0; i < nr_atoms_loop; i++)						/*Store coordinates current loop solution		*/
					for(int j = 0; j < 3; j++)
						if (reverse == 1)
						{
							co_sol->el[i   + k*nr_atoms_loop][j]   =    co->el[nr_atoms_loop-1-i][j];
							if(cg_mode>0)
								p1co_sol->el[i + k*nr_atoms_loop][j]   =  p1co->el[nr_atoms_loop-1-i][j];	/*store parallel side chain*/
						}
						else // if (reverse == 0)
						{
							co_sol->el[i   + k*nr_atoms_loop][j]   =    co->el[i][j];
							if(cg_mode>0)
								p1co_sol->el[i + k*nr_atoms_loop][j]   =  p1co->el[i][j];	/*store parallel side chain*/
						}

				if(factor_cacb > 0.0 && cg_mode > 0) // If detailed CA-CB distances were selected...
					// Set back CA-CB distances to standard length (bondCACB2 = 3.0604 A)
					set_cacb_dist2(co_sol->el + k*nr_atoms_loop, p1co_sol->el + k*nr_atoms_loop, nr_atoms_loop, bondCACB2/2, seq);

				if (reverse == 1)
				{
					for (int z = 1; z < nr_atoms_loop-2; z++)								// reverse dihedrals
						dihedrals_sol[k*nr_atoms_loop + z+2] = dihedral_angle[nr_atoms_loop-z];
					dihedrals_sol[k*nr_atoms_loop] = dihedral_angle[2]; // Last Psi
					dihedrals_sol[k*nr_atoms_loop+2] = dihedral_angle[0]; // First Phi
				}
				else
					for (int z = 0; z < nr_atoms_loop; z++)
						dihedrals_sol[k*nr_atoms_loop + z] = dihedral_angle[z];

				for (int z = 0; z < nr_atoms_loop; z++)
				{
					valence_sol[k*nr_atoms_loop + z] = valence_angle[z];
					length_sol[k*nr_atoms_loop + z] = bond_length[z];
				}

				// Store Energy statistics for OK solutions
				// MON: this is PD2 energy stuff...
				energies->array[k] = energy;
				energies_extra->array[k] = energy_extra;
				energies_intra->array[k] = energy_intra;

				// Compute LOCO-ICOS energy for OK solution
				if(emodel == 1 && presampling_continue) // ICOSA
				{
					// Convert "co->el" 2D-matrix into linear array (coordMatrix-like)
					co2coord(co_sol->el + k*nr_atoms_loop + 3, loop_coord,  nr_atoms_loop-6);

					// ICOSA energy calculation
					scores->array[k] = loco_energy(loop_coord, iseql, nr_atoms_loop/3-2, coorde, iseqe, nseqe , anchorCt, loco, icosVertices);
				}

				// Compute KORP energy for OK solution
				if( (emodel == 5 || emodel == 56 ) && presampling_continue) // KORP energy models
				{
					// Convert "co->el" 2D-matrix into linear array (coordMatrix-like)
					co2coord(co_sol->el + k*nr_atoms_loop + 3, loop_coord,  nr_atoms_loop-6);

					// Loop vs. Loop
					framesloop = frameCoord(loop_coord, nr_loop, resnumsloop, reschainidsloop);
					ncont = contactCoord(&contacts, korpmap->cutoff, framesloop, nr_loop, iseql);
					scores->array[k] = korp6D(contacts,ncont,korpmap); // Non-bonding = (nb2 < |i-j|)
					free(contacts);
					contacts = NULL; // Enables automatic memory allocation on next iteration

					// Loop vs. Environment
					ncont = contactCoord(&contacts, korpmap->cutoff, frames, nseqe, iseqe, anchorNt, framesloop, nr_loop, iseql);
					free(framesloop); // Free loop frames
					scores->array[k] += korp6D(contacts,ncont,korpmap); // Non-bonding = (nb2 < |i-j|)
					free(contacts);
					contacts = NULL; // Enables automatic memory allocation on next iteration
				}

				if( emodel == 6 || emodel == 56 ) // Rama energy models
				{
					if(emodel == 6)
						scores->array[k] = 0.0;

					scores->array[k] += rama_factor * rama_energy2(dihedrals_sol + k*nr_atoms_loop, nr_atoms_loop/3, ramaps, size);
				}

				// Compute Cross-Correlation between loop simulated map and input EM-map for OK solution
				if(em_switch && presampling_continue) // KORP
				{
					// Clear the loop simulated map
					loopmap->clear(0);

					// Project a "co" minimal backbone representation (just N, CA, and C atoms) into a pre-existent map.
					co2map_real( loopmap, loopmap->units().x() * 3, nr_atoms_loop, co_sol->el + k*nr_atoms_loop, p1co_sol->el + k*nr_atoms_loop, cg_mode, residuemarker);

					// FOPS::writeFile(loopmap,"loopmap.mrc");
					// exit(0);

					// scores->array[k] = 1 / FOPS::dotprod_mask(loopmap, workmap, workmask); // 1/dot --> The lower the better...
					// scores->array[k] = 1 - FOPS::dotprod_mask(loopmap, workmap, workmask); // 1-dot --> The lower the better...

					// Cross-correlation calculation
					if(emodel == 1 || emodel == 5 || emodel == 56 || emodel == 6) // Only save ccorrs if Energies are also requested
						ccorrs->array[k] = 1 - FOPS::correlation_mask(loopmap, workmap, workmask); // 1-dot --> The lower the better...
					else
						scores->array[k] = 1 - FOPS::correlation_mask(loopmap, workmap, workmask); // 1-dot --> The lower the better...
				}

				if(pid==0 && bump_filter && presampling_switch && !presampling_continue)
				{
					energies_presampling->array[k] = energy;
					if(k == presampling-1)
					{
						// After allocating "full cycles" memory for results, the master process only should do the corresponding "cycles/np" work...
						presampling_continue = true;
						total_counter = 0;
						totalk = 0;
						fails[0] = 0;
						fails[1] = 0;
						stat_n[0] = 0;
						stat_n[1] = 0;
						n_max_path = 0;
						save_counter = 0;
						k = -1; // -1 --> Because we do k++ later... (see a few lines below)
						energies_presampling->compute();
						if(pid==0)
							energies_presampling->print2("rcd> Pre-sampling energy ",6,stdout);
						bump_cutoff = bump_rate * energies_presampling->average;

#ifdef MPI_RELEASE
						cycles = cycles0 / np + cycles0 % np; // Tricky, but it should work...
						// everyone calls bcast, data is taken from root and ends up in everyone's buffer...
						MPI_Bcast(&bump_cutoff, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
						if(!server)
							fprintf(stderr,"Master (%d) broadcasted \"bump_cutoff\" (%f) to all workers...\n",pid,bump_cutoff);
#endif

						if(pid==0)
							fprintf(stdout,"rcd> Automatic bump cutoff (presampling= %d) = %f  (pid= %d)\n",presampling,bump_cutoff,pid);
					}
				}

				//  KEEP TRACK OF COUNTERS
				save_counter++;
				k++; // add counter actual solution
			}
			fflush(stdout);

			if(pid==0)
			{
				if(server || !presampling_continue)
				{
					if(k == cycles) // if server, only once...
						indicator("rcd> Loop closure ",k,cycles);
				}
				else
				{
					if(cycles > 50) // if more than 50 cycles...
					{
						if(k%(cycles/50)==0)
							indicator("rcd> Loop closure ",k,cycles);
					}
					else
						indicator("rcd> Loop closure ",k,cycles);
				}
			}
		}
		// ******************************
		// END ALL LOOP-CLOSURES
		// ******************************
		//#ifdef MPI_RELEASE
		//	MPI_Finalize();
		//#endif
		//exit(0);

#ifdef MPI_RELEASE
		// Send results from workers to master:  CC-results: co_sol and p1co_sol   IC-results: dihedrals_sol, valence_sol and length_sol
		if(pid==0) // I'm the Master
		{
			// Mon: Bug 6/4/2016 (--bump_filter option could not be removed...)
			cycles = cycles0 / np;
			int offset;
			offset = (cycles0 / np + cycles0 % np) * nr_atoms_loop; // accounting for the "extra" Master's closed loops

			// Translate the contiguous array into the complex "Matrix" object... (Fucking shit...)
			double *co_all = (double *) malloc(sizeof(double) * cycles * nr_atoms_loop * 3);
			double *p1co_all;
			if(cg_mode != 0) // if "p1co" exists
				p1co_all = (double *) malloc(sizeof(double) * cycles * nr_atoms_loop * 3);

			for(int id = 1; id < np; id++)
			{
				MPI_Recv( co_all, cycles*nr_atoms_loop*3, MPI_DOUBLE, id, 1, MPI_COMM_WORLD, &status);
				if(!server)
					fprintf(stderr,"Contiguous \"co\" array for pdb %d received from worker %d\n", i_pdb, status.MPI_SOURCE);
				for(int x = 0; x < cycles*nr_atoms_loop; x++) // screen all atoms for current coordinates chunk
				{
					co_sol->el[offset + x][0] = co_all[3*x]; // x-coord
					co_sol->el[offset + x][1] = co_all[3*x + 1]; // y-coord
					co_sol->el[offset + x][2] = co_all[3*x + 2]; // z-coord
				}

				if(cg_mode != 0) // if "p1co" exists
				{
					MPI_Recv( p1co_all, cycles*nr_atoms_loop*3, MPI_DOUBLE, id, 2, MPI_COMM_WORLD, &status);
					if(!server)
						fprintf(stderr,"Contiguous \"p1co\" array for pdb %d received from worker %d\n", i_pdb, status.MPI_SOURCE);
					for(int x = 0; x < cycles*nr_atoms_loop; x++) // screen all atoms for current coordinates chunk
					{
						p1co_sol->el[offset + x][0] = p1co_all[3*x]; // x-coord
						p1co_sol->el[offset + x][1] = p1co_all[3*x + 1]; // y-coord
						p1co_sol->el[offset + x][2] = p1co_all[3*x + 2]; // z-coord
					}
				}
				MPI_Recv( dihedrals_sol + offset, cycles*nr_atoms_loop, MPI_DOUBLE, id, 3, MPI_COMM_WORLD, &status);
				if(!server)
					fprintf(stderr,"\"dihedrals\" array for pdb %d received from worker %d\n", i_pdb, status.MPI_SOURCE);
				MPI_Recv( valence_sol + offset, cycles*nr_atoms_loop, MPI_DOUBLE, id, 4, MPI_COMM_WORLD, &status);

				if(!server)
					fprintf(stderr,"\"valence\" array for pdb %d received from worker %d\n", i_pdb, status.MPI_SOURCE);
				MPI_Recv( length_sol + offset, cycles*nr_atoms_loop, MPI_DOUBLE, id, 5, MPI_COMM_WORLD, &status);

				if(!server)
					fprintf(stderr,"\"length\" array for pdb %d received from worker %d\n", i_pdb, status.MPI_SOURCE);

				MPI_Recv( energies->array + offset/nr_atoms_loop, cycles, MPI_DOUBLE, id, 6, MPI_COMM_WORLD, &status);

				if(!server)
					fprintf(stderr,"\"energies\" array for pdb %d received from worker %d\n", i_pdb, status.MPI_SOURCE);
				MPI_Recv( energies_extra->array + offset/nr_atoms_loop, cycles, MPI_DOUBLE, id, 7, MPI_COMM_WORLD, &status);

				if(!server)
					fprintf(stderr,"\"energies_extra\" array for pdb %d received from worker %d\n", i_pdb, status.MPI_SOURCE);
				MPI_Recv( energies_intra->array + offset/nr_atoms_loop, cycles, MPI_DOUBLE, id, 8, MPI_COMM_WORLD, &status);

				if(!server)
					fprintf(stderr,"\"energies_intra\" array for pdb %d received from worker %d\n", i_pdb, status.MPI_SOURCE);
				MPI_Recv( scores->array + offset/nr_atoms_loop, cycles, MPI_DOUBLE, id, 9, MPI_COMM_WORLD, &status);

				if(!server)
					fprintf(stderr,"\"scores\" array for pdb %d received from worker %d\n", i_pdb, status.MPI_SOURCE);

				if(em_switch && (emodel == 1 || emodel == 5 || emodel == 56 || emodel == 6) ) // Only save ccorrs if Energies are also requested
				{
					MPI_Recv( ccorrs->array + offset/nr_atoms_loop, cycles, MPI_DOUBLE, id,10, MPI_COMM_WORLD, &status);
					if(!server)
						fprintf(stderr,"\"ccorrs\" array for pdb %d received from worker %d\n", i_pdb, status.MPI_SOURCE);
				}

				offset += cycles*nr_atoms_loop;
			}
			free(co_all);
			if(cg_mode != 0) // if "p1co" exists
				free(p1co_all);

			// Once results were sent to master, the master process can do things as usual...
			cycles = cycles0; // Setting back "cycles" to the initial value
		}
		else // I'm a Worker
		{
			// Translate the complex "Matrix" object into a contiguous array... (Fucking shit...)
			double *co_all = (double *) malloc(sizeof(double) * cycles * nr_atoms_loop * 3);
			double *p1co_all;
			for(int x=0; x < cycles * nr_atoms_loop; x++) // screen all atoms
			{
				co_all[3*x]     = co_sol->el[x][0]; // x-coord
				co_all[3*x + 1] = co_sol->el[x][1]; // y-coord
				co_all[3*x + 2] = co_sol->el[x][2]; // z-coord
			}

			// Send data chunk to master process (pid=0) from current worker (pid)
			if(!server)
				fprintf(stderr,"Worker %d sending \"co\" for pdb %d\n",pid,i_pdb);
			MPI_Send( co_all, cycles * nr_atoms_loop * 3, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD); // Sending data chunk...
			free(co_all); // release memory

			if(cg_mode != 0) // if "p1co" exists
			{
				p1co_all = (double *) malloc(sizeof(double) * cycles * nr_atoms_loop * 3);
				for(int x=0; x < cycles * nr_atoms_loop; x++) // screen all atoms
				{
					p1co_all[3*x]     = p1co_sol->el[x][0];
					p1co_all[3*x + 1] = p1co_sol->el[x][1];
					p1co_all[3*x + 2] = p1co_sol->el[x][2];
				}
				if(!server)
					fprintf(stderr,"Worker %d sending \"p1co\" for pdb %d\n",pid,i_pdb);
				MPI_Send( p1co_all, cycles * nr_atoms_loop * 3, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD); // Sending data chunk...
				free(p1co_all); // release memory
			}

			if(!server)
				fprintf(stderr,"Worker %d sending \"dihedrals\" for pdb %d\n",pid,i_pdb);
			MPI_Send( dihedrals_sol, cycles * nr_atoms_loop, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD); // Sending data chunk...
			free(dihedrals_sol);

			if(!server)
				fprintf(stderr,"Worker %d sending \"valence\" for pdb %d\n",pid,i_pdb);
			MPI_Send( valence_sol, cycles * nr_atoms_loop, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD); // Sending data chunk...
			free(valence_sol);

			if(!server)
				fprintf(stderr,"Worker %d sending \"length\" for pdb %d\n",pid,i_pdb);
			MPI_Send( length_sol, cycles * nr_atoms_loop, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD); // Sending data chunk...
			free(length_sol);

			if(!server)
				fprintf(stderr,"Worker %d sending \"energies\" for pdb %d\n",pid,i_pdb);
			MPI_Send( energies->array, cycles, MPI_DOUBLE, 0, 6, MPI_COMM_WORLD); // Sending data chunk...
			//delete energies;

			if(!server)
				fprintf(stderr,"Worker %d sending \"energies_extra\" for pdb %d\n",pid,i_pdb);
			MPI_Send( energies_extra->array, cycles, MPI_DOUBLE, 0, 7, MPI_COMM_WORLD); // Sending data chunk...
			//delete energies_extra;

			if(!server)
				fprintf(stderr,"Worker %d sending \"energies_intra\" for pdb %d\n",pid,i_pdb);
			MPI_Send( energies_intra->array, cycles, MPI_DOUBLE, 0, 8, MPI_COMM_WORLD); // Sending data chunk...
			//delete energies_intra;


			if(!server)
				fprintf(stderr,"Worker %d sending \"scores\" for pdb %d\n",pid,i_pdb);
			MPI_Send( scores->array, cycles, MPI_DOUBLE, 0, 9, MPI_COMM_WORLD); // Sending data chunk...
			//delete scores;

			if(em_switch && (emodel == 1 || emodel == 5 || emodel == 56 || emodel == 6) ) // Only save ccorrs if Energies are also requested
			{
				if(!server)
					fprintf(stderr,"Worker %d sending \"ccorrs\" for pdb %d\n",pid,i_pdb);
				MPI_Send( ccorrs->array, cycles, MPI_DOUBLE, 0,10, MPI_COMM_WORLD); // Sending data chunk...
				//delete ccorrs;

			}
		}
#endif

#ifdef MPI_RELEASE
		MPI_Barrier(MPI_COMM_WORLD);
#endif




		// Energy filtration
		if(pid==0 && (emodel == 1 || emodel == 5 || emodel == 56 || emodel == 6 || em_switch) && loco_best > 0)
		{
			// fprintf(stderr,"MON: cycles= %d  loco_best= %d\n", cycles, loco_best);
			if(pid==0)
				fprintf(stderr,"rcd> Energy filtration: cycles= %d  loco_best= %d  loco_notbest= %d\n",cycles,loco_best,loco_notbest);

			// Allocate memory for the sorted array of energies
			double **se = (double **) malloc(sizeof(double *) * cycles);
			// Initialize sorted array of energies
			float se_min=1000000;
			float se_max=-1000000;
			for(int u = 0; u < cycles; u++) {
				se[u] = &scores->array[u];
				if (*se[u]<se_min) se_min=*se[u];
				if (*se[u]>se_max) se_max=*se[u];
			}
			// Sort array of energies
			quicksort(se, 0, cycles-1); // Quick sort algorithm

			// Energy cutoff to get the "loco_best" energies
			double ecut = *se[loco_best];
			free(se); // not needed anymore...

			if(pid==0)
				fprintf(stderr,"rcd> Energy filtration: Max %f Min %f cutoff %f\n",se_min,se_max,ecut);


			// Generate an indices array with best energy loops
			int *ibest = (int *) malloc(sizeof(int) * (loco_best+loco_notbest));
			int ei = 0; // best energy index
			for(int u = 0; u < cycles; u++)
				if(scores->array[u] < ecut)
				{
					ibest[ei] = u; // store best indices
					ei++; // increase best energy counter
				}
			// Random and efficiently selecting loops of not-best LOCO-ICOSA energy
			int rnb = loco_notbest; // Number of remaining not-best loops
			for(int u = 0; u < cycles && rnb > 0; u++)
				if(scores->array[u] >= ecut)
				{
					ibest[ei] = u; // store best indices
					ei++; // increase best energy counter
					rnb--; // decrease the number of remaining not-best loops
				}

			// Load best energy loops into "best" arrays
			Stat *energies2 	= new Stat((loco_best+loco_notbest));    //  Stat object for Energies
			Stat *energies_intra2	= new Stat((loco_best+loco_notbest));    //  Stat object for Intra-Energies
			Stat *energies_extra2	= new Stat((loco_best+loco_notbest));    //  Stat object for Extra-Energies
			Stat *scores2           = new Stat((loco_best+loco_notbest));    //  Stat object for Scores (energies or cross-correlations)
			for(int u = 0; u < loco_best+loco_notbest; u++)
			{
				energies2->array[u] = energies->array[ibest[u]];
				energies_intra2->array[u] = energies_intra->array[ibest[u]];
				energies_extra2->array[u] = energies_extra->array[ibest[u]];
				scores2->array[u] = scores->array[ibest[u]];
				// fprintf(stderr,"u= %d  scores2->array[u]= %f   scores->array[ibest[u]]= %f \n",u,scores2->array[u],scores->array[ibest[u]]);
			}
			for(int u = 0; u < loco_best+loco_notbest; u++)
			{
				energies->array[u] = energies2->array[u];
				energies_intra->array[u] = energies_intra2->array[u];
				energies_extra->array[u] = energies_extra2->array[u];
				scores->array[u] = scores2->array[u];
				// fprintf(stderr,"u= %d  scores2->array[u]= %f   scores->array[ibest[u]]= %f \n",u,scores2->array[u],scores->array[ibest[u]]);
			}
			// delete swap data
			delete energies2;
			delete energies_intra2;
			delete energies_extra2;
			delete scores2;

			//
			// swap sol array & best arrays
			//
			// swap co_best
			co_best  = new Matrix(nr_atoms_loop * (loco_best+loco_notbest),3);
			for(int u = 0; u < loco_best+loco_notbest; u++)
			{
				for(int i = 0; i < nr_atoms_loop; i++)
				{
					co_best->el[u*nr_atoms_loop + i][0] = co_sol->el[ibest[u]*nr_atoms_loop + i][0];
					co_best->el[u*nr_atoms_loop + i][1] = co_sol->el[ibest[u]*nr_atoms_loop + i][1];
					co_best->el[u*nr_atoms_loop + i][2] = co_sol->el[ibest[u]*nr_atoms_loop + i][2];
				}
			}
			for(int u = 0; u < loco_best+loco_notbest; u++)
			{
				for(int i = 0; i < nr_atoms_loop; i++)
				{
					co_sol->el[u*nr_atoms_loop + i][0] = co_best->el[u*nr_atoms_loop + i][0];
					co_sol->el[u*nr_atoms_loop + i][1] = co_best->el[u*nr_atoms_loop + i][1];
					co_sol->el[u*nr_atoms_loop + i][2] = co_best->el[u*nr_atoms_loop + i][2];
				}
			}
			delete co_best;

			// swap p1co_best
			p1co_best  = new Matrix(nr_atoms_loop * (loco_best+loco_notbest),3);
			for(int u = 0; u < loco_best+loco_notbest; u++)
			{
				for(int i = 0; i < nr_atoms_loop; i++)
				{
					p1co_best->el[u*nr_atoms_loop + i][0] = p1co_sol->el[ibest[u]*nr_atoms_loop + i][0];
					p1co_best->el[u*nr_atoms_loop + i][1] = p1co_sol->el[ibest[u]*nr_atoms_loop + i][1];
					p1co_best->el[u*nr_atoms_loop + i][2] = p1co_sol->el[ibest[u]*nr_atoms_loop + i][2];

				}
			}
			for(int u = 0; u < loco_best+loco_notbest; u++)
			{
				for(int i = 0; i < nr_atoms_loop; i++)
				{
					p1co_sol->el[u*nr_atoms_loop + i][0] = p1co_best->el[u*nr_atoms_loop + i][0] ;
					p1co_sol->el[u*nr_atoms_loop + i][1] = p1co_best->el[u*nr_atoms_loop + i][1];
					p1co_sol->el[u*nr_atoms_loop + i][2] = p1co_best->el[u*nr_atoms_loop + i][2];

				}
			}
			delete p1co_best;

			// swap dihedrals
			dihedrals_best = (double *) malloc(sizeof(double) * (loco_best+loco_notbest) * nr_atoms_loop);
			for(int u = 0; u < loco_best+loco_notbest; u++)
				for(int i = 0; i < nr_atoms_loop; i++)
					dihedrals_best[u*nr_atoms_loop + i] = dihedrals_sol[ibest[u]*nr_atoms_loop + i];
			for(int u = 0; u < loco_best+loco_notbest; u++)
				for(int i = 0; i < nr_atoms_loop; i++)
					dihedrals_sol[u*nr_atoms_loop + i] = dihedrals_best[u*nr_atoms_loop + i];
			free(dihedrals_best);

			// swap valence_best
			valence_best = (double *) malloc(sizeof(double) * (loco_best+loco_notbest) * nr_atoms_loop);
			for(int u = 0; u < loco_best+loco_notbest; u++)
				for(int i = 0; i < nr_atoms_loop; i++)
					valence_best[u*nr_atoms_loop + i] = valence_sol[ibest[u]*nr_atoms_loop + i];
			for(int u = 0; u < loco_best+loco_notbest; u++)
				for(int i = 0; i < nr_atoms_loop; i++)
					valence_sol[u*nr_atoms_loop + i] = valence_best[u*nr_atoms_loop + i];
			free(valence_best);

			// swap length_best
			length_best = (double *) malloc(sizeof(double) * (loco_best+loco_notbest) * nr_atoms_loop);
			for(int u = 0; u < loco_best+loco_notbest; u++)
				for(int i = 0; i < nr_atoms_loop; i++)
					length_best[u*nr_atoms_loop + i] = length_sol[ibest[u]*nr_atoms_loop + i];
			for(int u = 0; u < loco_best+loco_notbest; u++)
				for(int i = 0; i < nr_atoms_loop; i++)
					length_sol[u*nr_atoms_loop + i] = length_best[u*nr_atoms_loop + i];
			free(length_best);

			if(em_switch && (emodel == 1 || emodel == 5 || emodel == 56 || emodel == 6) ) // Only save ccorrs if Energies are also requested
			{
				Stat *ccorrs2 = NULL; //  Stat object for Cross-correlations
				ccorrs2           = new Stat((loco_best+loco_notbest));    //  Stat object for Cross-correlations
				for(int u = 0; u < loco_best+loco_notbest; u++)
					ccorrs2->array[u] = ccorrs->array[ibest[u]]; // Save only the best cross-correlations
				delete ccorrs; // this must be here to fix some Bug
				ccorrs = ccorrs2;
			}

			free(ibest); // not needed anymore...
			cycles = loco_best+loco_notbest; // Warning, from now on "cycles" is different! (sorry, but this makes things much easier...)
		}



		// Mon: this should be an independent I/O function... pid==0????
		if(pid==0)
		{
			// Boundary Internal Coordinates must be recalculated to obtain the true ICs (for PyRosetta)
			//
			double bond1[3],bond2[3];
			for(int k=0; k<cycles; k++)
			{

		//		printf("cycles %d %d %d\n", cycles, loco_best, loco_notbest);
				for(int i=0; i<3; i++) { // screen residues
			    co_sol->el[k*nr_atoms_loop+0][i]=coaux->el[0][i];
			    co_sol->el[k*nr_atoms_loop+1][i]=coaux->el[1][i];
			    co_sol->el[k*nr_atoms_loop+2][i]=coaux->el[2][i];
			    co_sol->el[k*nr_atoms_loop+nr_atoms_loop-3][i]=coaux->el[nr_atoms_loop-3][i];
			    co_sol->el[k*nr_atoms_loop+nr_atoms_loop-2][i]=coaux->el[nr_atoms_loop-2][i];
			    co_sol->el[k*nr_atoms_loop+nr_atoms_loop-1][i]=coaux->el[nr_atoms_loop-1][i];
				}
			    //getchar();

			    // recompute everything
			    findlength(&co_sol->el[k*nr_atoms_loop],&length_sol[k*nr_atoms_loop],nr_atoms_loop);
			    valencefind(&co_sol->el[k*nr_atoms_loop],&valence_sol[k*nr_atoms_loop],nr_atoms_loop);
			    finddihedral(&co_sol->el[k*nr_atoms_loop],&dihedrals_sol[k*nr_atoms_loop], nr_atoms_loop);

			    dihedrals_sol[k*nr_atoms_loop+2] = FirstPhi;
			   	dihedrals_sol[k*nr_atoms_loop+0] = LastPsi;


				// Dihedral angles

//				dihedrals_sol[k*nr_atoms_loop + 3] = get_dihedral(coaux->el[0],coaux->el[1],coaux->el[2],co_sol->el[k*nr_atoms_loop+3]); // Psi(1)
//				dihedrals_sol[k*nr_atoms_loop + 4] = get_dihedral(coaux->el[1],coaux->el[2],co_sol->el[k*nr_atoms_loop+3],co_sol->el[k*nr_atoms_loop+4]); // Omega(2)
//				dihedrals_sol[k*nr_atoms_loop + 5] = get_dihedral(coaux->el[2],co_sol->el[k*nr_atoms_loop+3],co_sol->el[k*nr_atoms_loop+4],co_sol->el[k*nr_atoms_loop+5]); // Phi(2)
//
//				dihedrals_sol[(k+1)*nr_atoms_loop - 3] = get_dihedral(co_sol->el[(k+1)*nr_atoms_loop-6],co_sol->el[(k+1)*nr_atoms_loop-5],co_sol->el[(k+1)*nr_atoms_loop-4],coaux->el[nr_atoms_loop-3]); // Psi(N-1)
//				dihedrals_sol[(k+1)*nr_atoms_loop - 2] = get_dihedral(co_sol->el[(k+1)*nr_atoms_loop-5],co_sol->el[(k+1)*nr_atoms_loop-4],coaux->el[nr_atoms_loop-3],coaux->el[nr_atoms_loop-2]); // Omega(N)
//				dihedrals_sol[(k+1)*nr_atoms_loop - 1] = get_dihedral(co_sol->el[(k+1)*nr_atoms_loop-4],coaux->el[nr_atoms_loop-3],coaux->el[nr_atoms_loop-2],coaux->el[nr_atoms_loop-1]); // Phi(N)
//				// Bond angles (valence)
//				subtract(coaux->el[1],coaux->el[2],bond1); // Bond vector CA(1)-C(1)
//				subtract(coaux->el[2],co_sol->el[k*nr_atoms_loop+3],bond2); // Bond vector C(1)-N(2)
//				valence_sol[k*nr_atoms_loop+3] = M_PI- acos(dotprod(bond1,bond2)/(norm(bond1)*norm(bond2)));	// CA(1)^C(N)^N(2) bond angle (in radians)
//				subtract(co_sol->el[k*nr_atoms_loop+3],co_sol->el[k*nr_atoms_loop+4],bond1); // Bond vector N(2)-CA(2)
//				valence_sol[k*nr_atoms_loop+4] = M_PI- acos(dotprod(bond2,bond1)/(norm(bond2)*norm(bond1)));	// C(1)^N(2)^CA(2) bond angle (in radians)
//
//				subtract(co_sol->el[(k+1)*nr_atoms_loop-5],co_sol->el[(k+1)*nr_atoms_loop-4],bond1); // Bond vector CA(N-1)-C(N-1)
//				subtract(co_sol->el[(k+1)*nr_atoms_loop-4],coaux->el[nr_atoms_loop-3],bond2); // Bond vector C(N-1)-N(N)
//				valence_sol[(k+1)*nr_atoms_loop - 3] = M_PI- acos(dotprod(bond1,bond2)/(norm(bond1)*norm(bond2)));	// CA(N-1)^C(N-1)^N(N) bond angle (in radians)
//				subtract(coaux->el[nr_atoms_loop-3],coaux->el[nr_atoms_loop-2],bond1); // Bond vector N(N)-CA(N)
//				valence_sol[(k+1)*nr_atoms_loop - 2] = M_PI- acos(dotprod(bond2,bond1)/(norm(bond2)*norm(bond1)));	// C(0)^N(1)^CA(1) bond angle (in radians)
//				// Bond lengths
//				length_sol[k*nr_atoms_loop + 3] = distance(coaux->el[2],co_sol->el[k*nr_atoms_loop+3]); // Distance C(1)-N(2)
//				length_sol[k*nr_atoms_loop + 2] = distance(coaux->el[1],coaux->el[2]); //
//				length_sol[k*nr_atoms_loop + 1] = distance(coaux->el[0],coaux->el[1]); //
//
//
//				length_sol[(k+1)*nr_atoms_loop - 3] = distance(co_sol->el[(k+1)*nr_atoms_loop-4],coaux->el[nr_atoms_loop-3]); // Distance C(N-1)-N(N)
//				length_sol[(k+1)*nr_atoms_loop - 2] = distance(coaux->el[nr_atoms_loop-3],coaux->el[nr_atoms_loop-2]); // Distance C(N-1)-N(N)
//				length_sol[(k+1)*nr_atoms_loop - 1] = distance(coaux->el[nr_atoms_loop-2],coaux->el[nr_atoms_loop-1]); // Distance C(N-1)-N(N)

			}

			if(output_dihedrals && pid==0)
			{
				char title[20];
				FILE *file;

				// DIHEDRAL ANGLES
				//
				sprintf(file_name,"%s/%s_dh.txt",folder_later,current_pdb_name);
				fprintf(stderr,"rcd> Writing dihedral angles file (Phi,Psi): %s\n",file_name);
				file = fopen(file_name,"w");

				// Write header for dihedral angles
				for(int i=0; i<nr_loop+2; i++) // screen residues
				{
					sprintf(title,"%s%d","Omega_",i+1);
					if(i==0)
						fprintf(file,"#%16s",title);
					else
						fprintf(file," %16s",title);
					sprintf(title,"%s%d","Phi_",i+1);
					fprintf(file," %16s",title);
					sprintf(title,"%s%d","Psi_",i+1);
					fprintf(file," %16s",title);
				}
				fprintf(file,"\n");

				for(int k=0; k<cycles; k++)
				{
					for(int i=0; i<nr_loop+2; i++) // screen residues
					{
						if(i==nr_loop+1)
							fprintf(file," %16.12lf %16.12lf %16.12lf", dihedrals_sol[k*nr_atoms_loop + i*3 + 1], dihedrals_sol[k*nr_atoms_loop + i*3 + 2], dihedrals_sol[k*nr_atoms_loop] ); // Omega_n, Phi_n, Psi_n
						else
							fprintf(file," %16.12lf %16.12lf %16.12lf", dihedrals_sol[k*nr_atoms_loop + i*3 + 1], dihedrals_sol[k*nr_atoms_loop + i*3 +2], dihedrals_sol[k*nr_atoms_loop + i*3 + 3] ); //  // Omega_i, Phi_i, Psi_i
					}
					fprintf(file,"\n");
				}
				fclose(file);

				// VALENCE ANGLES
				//
				sprintf(file_name,"%s/%s_val.txt",folder_later,current_pdb_name); // pdbfile[i_pdb].string);
				fprintf(stderr,"rcd> Writing valence angles file: %s\n",file_name);
				file = fopen(file_name,"w");

				// Write header for valence angles
				for(int i=0; i<nr_loop+2; i++) // screen residues
				{
					sprintf(title,"%s%d","NCAC_",i+1);
					if(i==0)
						fprintf(file,"# %16s",title);
					else
						fprintf(file," %16s",title);
					if(i != nr_loop+1)
					{
						sprintf(title,"%s%d","CACN_",i+1);
						fprintf(file," %16s",title);
						sprintf(title,"%s%d","CNCA_",i+1);
						fprintf(file," %16s",title);
					}
				}
				fprintf(file,"\n");

				for(int k=0; k<cycles; k++)
				{
					// fprintf(file,"  %9f %9f", 0.0, 0.0);
					fprintf(file," ");
					for(int i=2; i<nr_atoms_loop; i++) // screen residues
						fprintf(file," %16.12lf", valence_sol[k*nr_atoms_loop + i]);
					fprintf(file,"\n");
				}
				fclose(file);

				// BOND LENGTHS
				//
				sprintf(file_name,"%s/%s_len.txt",folder_later,current_pdb_name); // pdbfile[i_pdb].string);
				fprintf(stderr,"rcd> Writing bond lengths file: %s\n",file_name);
				file = fopen(file_name,"w");

				// Write header for bond lengths
				for(int i=0; i<nr_loop+2; i++) // screen residues
				{
					sprintf(title,"%s%d","N-CA_",i+1);
					if(i==0)
						fprintf(file,"#%16s",title);
					else
						fprintf(file," %16s",title);
					sprintf(title,"%s%d","CA-C_",i+1);
					fprintf(file," %16s",title);
					if(i != nr_loop+1)
					{
						sprintf(title,"%s%d","C-N_",i+1);
						fprintf(file," %16s",title);
					}
				}
				fprintf(file,"\n");

				for(int k=0; k<cycles; k++)
				{
					// fprintf(file,"  %9f %9f", 0.0, 0.0);
					fprintf(file," ");
					for(int i=1; i<nr_atoms_loop; i++) // screen residues
						fprintf(file," %16.12lf", length_sol[k*nr_atoms_loop + i]);
					fprintf(file,"\n");
				}
				fclose(file);
			}
		}

		// Debugging dihedrals for server
		if(pid==0 && (output_debug || server))
		{
			// Log-file name
			strcpy(file_name,folder_later);
			strcat(file_name,"/");
			strcat(file_name,current_pdb_name);
			strcat(file_name,"_dhs_");
			paint_all_dihedrals(dihedrals_sol, size, nr_loop+2, 0, file_name, cycles );
		}

		// Debugging Ramachandran maps for server
		if(pid==0 && (server || output_debug) && load_ramap)
		{
			// Plot the Ramachandran maps for all aminoacids
			strcpy(file_name,folder_later);
			strcat(file_name,"/");
			strcat(file_name,current_pdb_name);
			strcat(file_name,"_map_");
			paint_maps_points(binmaps, size, nr_loop+2, file_name);
		}

		// For clustering...
		FILE *f_clust;
		float simil;
		if(pid==0 && simIC_switch)
		{
			if(!sim_text)
			{
				sprintf(tmp_path,"%s/%s_simic.bin",folder_later,current_pdb_name); // pdbfile[i_pdb].string);
				f_clust = fopen(tmp_path,"wb");
				fprintf(stderr,"Writing similarity IC-RMSD file: %s\n",tmp_path);
				fwrite(&cycles,sizeof(int),1,f_clust); // Number of loops
				// Computing similarity IC-RMSD
				for(int i=0; i<cycles; i++)
					for(int j=i+1; j<cycles; j++)
					{
						simil = 0.0;
						for(int k=0; k<nr_atoms_loop; k++)
							simil += pow(dihedrals_sol[nr_atoms_loop*i+k] - dihedrals_sol[nr_atoms_loop*j+k],2);
						simil = sqrtf(simil/nr_atoms_loop);
						fwrite(&simil,sizeof(float),1,f_clust); // Similarity
					}
			}
			else
			{
				sprintf(tmp_path,"%s/%s_simic.txt",folder_later,current_pdb_name); // pdbfile[i_pdb].string);
				f_clust = fopen(tmp_path,"w");
				fprintf(stderr,"Writing similarity IC-RMSD file: %s\n",tmp_path);
				fprintf(f_clust,"# <i> <j> <similarity>\n"); // Matlab request (i begins in value 1)
				// Computing similarity IC-RMSD
				for(int i=0; i<cycles; i++)
					for(int j=i+1; j<cycles; j++)
					{
						simil = 0.0;
						for(int k=0; k<nr_atoms_loop; k++)
							simil += pow(dihedrals_sol[nr_atoms_loop*i+k] - dihedrals_sol[nr_atoms_loop*j+k],2);
						simil = sqrt(simil/nr_atoms_loop);
						fprintf(f_clust,"%d %d %f\n",i+1,j+1,simil); // Matlab request (i begins in value 1)
					}
			}
			fclose(f_clust);
		}

		if(pid==0 && simCC_switch)
		{
			fprintf(stderr,"BEGIN similarity computation\n");

			if(!sim_text) // binary
			{
				sprintf(tmp_path,"%s/%s_simcc.bin",folder_later,current_pdb_name); // pdbfile[i_pdb].string);
				f_clust = fopen(tmp_path,"wb");
				fprintf(stderr,"Writing similarity CC-RMSD file: %s\n",tmp_path);
				fwrite(&cycles,sizeof(int),1,f_clust); // Number of loops
				// Computing similarity CC-RMSD
				for(int i=0; i<cycles; i++)
					for(int j=i+1; j<cycles; j++)
					{
						simil = (float) rmsd_co_p1co(co_sol->el+i*nr_atoms_loop, p1co_sol->el+i*nr_atoms_loop, co_sol->el+j*nr_atoms_loop, p1co_sol->el+j*nr_atoms_loop,nr_atoms_loop,cg_mode);
						fwrite(&simil,sizeof(float),1,f_clust); // Similarity
					}
			}
			else // plain-text
			{
				sprintf(tmp_path,"%s/%s_simcc.txt",folder_later,current_pdb_name); // pdbfile[i_pdb].string);
				f_clust = fopen(tmp_path,"w");
				fprintf(stderr,"Writing similarity CC-RMSD file: %s\n",tmp_path);
				fprintf(f_clust,"# <i> <j> <similarity>\n"); // Matlab request (i begins in value 1)
				// Computing similarity CC-RMSD
				for(int i=0; i<cycles; i++)
					for(int j=i+1; j<cycles; j++)
					{
						simil = rmsd_co_p1co(co_sol->el+i*nr_atoms_loop, p1co_sol->el+i*nr_atoms_loop, co_sol->el+j*nr_atoms_loop, p1co_sol->el+j*nr_atoms_loop,nr_atoms_loop,cg_mode);
						fprintf(f_clust,"%d %d %f\n",i+1,j+1,simil); // Matlab request (i begins in value 1)
					}
			}
			fclose(f_clust);
			fprintf(stderr,"END similarity computation\n");
		}
		//exit(0);
		if(pid==0 && (bench_rmsd || native==1))
		{
			fprintf(stderr,"rcd> Computing %d RMSDs between loops and native\n",cycles);
			//		frmsd->array[i] = total_rmsd_co(co_i->el,&co_sol->el[nr_atoms_loop*i],nr_atoms_loop);
			if(rmsd_O) // including O-atom in rmsd computations
			{
				for(int i = 0; i < cycles; i++)		/*Calculate and store overall RMSDs*/
					// RMSD is computed using "moving" residues only (not anchors)
					frmsd->array[i] = total_rmsd_NCACO(co_i->el,&co_sol->el[nr_atoms_loop*i],p1co_i->el,&p1co_sol->el[nr_atoms_loop*i],nr_loop);
			}
			else
			{
				for(int i = 0; i < cycles; i++)		/*Calculate and store overall RMSDs*/
					// RMSD is computed using "moving" residues only (not anchors)
					frmsd->array[i] = total_rmsd_NCAC(co_i->el,&co_sol->el[nr_atoms_loop*i],nr_loop);
			}
			frmsd->compute(); // Mon: consider removal...
		}
		delete(loop);

		// Set "cycles" back...
		//		if(loco_switch)
		if(emodel == 1 || emodel == 5 || emodel == 56 || emodel == 6)
			cycles = cycles0; // required for "best energy" stuff...

		if(pid==0)
		{
			fprintf(stderr,"rcd> Cycles %d  trials %d F/R %d/%d  map_clashes %d  rand_repairs %d  infinite= %d\n"
					,cycles, totalk,trials[0],trials[1],n_max_grid_clash, n_rand_repairs,n_infinite);

			fprintf(stderr,"rcd> Fails %d  F/R %d/%d  max_path %d  Total_self_clashes %d  Open/Closed_self_clashes %d/%d\n"
					,(fails[0]+fails[1]),fails[0],fails[1],n_max_path, self_clash_count_failed, self_clash_count, (self_clash_count_failed-self_clash_count));
			// if(postcheck >= 0 || enable_lc)
			fprintf(stderr,"rcd> n_rama_fail= %d   avg_rama_fail= %d\n", n_rama_fail, avg_rama_fail);
		}
		//	if(enable_lc)
		//		fprintf(stderr,"rcd> LoopClosure= %d  almost_closed= %d  bad_closures= %d  good_closures= %d\n", n_LoopClosure, n_almost_closed, n_bad_closures, n_good_closures);

		// ---------------------------
		// END LOOP CLOSURE ITERATIONS
		// ---------------------------

		// STOP SLOW TIMER
		timer1->stop();
		timer1->elapsed();

		free(spinor1);
		free(spinor2);
		free(spinor3);
		free(spinor4);
		free(spinor5);
		free(spinor6);
		free(spinor7);
		free(full);
		free(posespinor);
		free(posespinor2);
		free(between);
		free(bond);
		free(rotor);
		free(result);
		free(s3tot);

		freepointer(coord,3);
		freepointer(coord2,3);

		//	// ULTIMATE CROSS-CHECK BETWEEN INTERNAL COORDINATES AND CARTESIAN COORDINATES
		//	dump_dihedrals("Lengths IC",bond_length,nr_atoms_loop);
		//	findlength(co->el,bond_length,nr_atoms_loop);
		//	dump_dihedrals("Lengths CC",bond_length,nr_atoms_loop);
		//	dump_dihedrals("Valences IC",valence_angle,nr_atoms_loop);
		//	valencefind(co->el,valence_angle,nr_atoms_loop);
		//	dump_dihedrals("Valences CC",valence_angle,nr_atoms_loop);
		//	dump_dihedrals("Dihedrals IC",dihedral_angle,nr_atoms_loop);
		//	finddihedral(co->el,dihedral_angle,nr_atoms_loop);
		//	dump_dihedrals("Dihedrals CC",dihedral_angle,nr_atoms_loop);
		delete(vol);

		//printf("rcd>Total elapsed time %f ms", timer1->felapsed();				//	TIMERS
		//timer1->print("",stdout);

		if(pid==0)
		{
			//	STATISTICS
			ftimer->compute();
			ftimer->print2("rcd> Timings",5,stdout);
			fcounter->compute();
			fcounter->print2("rcd> Counts ",5,stdout);
			frmsd->compute();
			frmsd->print2("rcd> Rmsd   ",6,stdout);
			energies->compute();
			energies->print2("rcd> Energy ",6,stdout);
			energies_extra->compute();
			energies_intra->compute();
			scores->compute();
			scores->print2("rcd> Score",6,stdout);
			if(em_switch && (emodel == 1 || emodel == 5 || emodel == 56 || emodel == 6) ) // Only do ccorrs stuff when energies computation was also required
			{
				ccorrs->compute();
				ccorrs->print2("rcd> Cross-correlations with EM-maps",6,stdout);
			}
			printf("rcd> Post-rama-check -> Failed-post-rama-loops: %d\n", counter);
			printf("rcd> Final solutions -> Trials (totalk): %d \n", totalk);
			printf("rcd> Solutions/trials %d/%d ratio %.3f   F/B %4.2f/%4.2f %4.2f/%4.2f %4.2f/%4.2f %4.2f/%4.2f\n"
					,k,totalk,  double (k)/ double (totalk),
					trials[0]/ double (totalk),
					trials[1]/ double (totalk),
					(trials[0]-fails[0])/(double) k,
					(trials[1]-fails[1])/(double) k,
					(trials[0]-fails[0])/double (totalk),
					(trials[1]-fails[1])/double (totalk),
					1-(fails[0]/(trials[0]+1.0)+0.01)/(fails[1]/(trials[1]+1.0)+fails[0]/(trials[0]+1.0)+0.02),
					1-(fails[1]/(trials[1]+1.0)+0.01)/(fails[1]/(trials[1]+1.0)+fails[0]/(trials[0]+1.0)+0.02)
			);
			printf("rcd> Solutions with clashes remaining %d ratio %.2f\n",counter_clash, double (counter_clash) /double(cycles) );
		}

		if (pid==0 && switch_file == 1)
		{
			char file_name_output[LENGTH_LINE];  /*Char array[] to safe created output filename	*/

			time_t rawtime;		/*Time object to retrieve time	*/
			struct tm * timeinfo_f;
			time ( &rawtime );
			timeinfo_f = localtime ( &rawtime );

			// Log-file name
			strcpy(tmp_path,folder_later);
			strcat(tmp_path,"/");
			strcat(tmp_path,current_pdb_name);
			strcpy(file_name_output,tmp_path);	/*PDB format*/
			strcat(file_name_output,".log");

			FILE*	p_file;
			p_file = fopen(file_name_output,"w");

			fprintf(p_file,"\n");		/*File output name & time details	*/
			fprintf(p_file,"______________________________________________________________________________________________\n");
			fprintf(p_file,"\n");
			fprintf(p_file,"%-30s %40s",file_name_output, asctime (timeinfo_f) );
			fprintf(p_file,"______________________________________________________________________________________________\n");
			fprintf(p_file,"\n");
			fprintf(p_file,"%-30s %40s\n","Input file:",pdbfile[i_pdb].string);
			fprintf(p_file,"\n");
			//ALEX: MI STATISTICS
			fprintf(p_file,"\n");
			fprintf(p_file,"SOLUTIONS/TRIALS IN GENERATION LOOP \n");
			fprintf(p_file,"___________________________________________________________________________________________________\n");
			fprintf(p_file,"rmsd_avg_solutions	rmsd_avg_trials		solutions	 trials 	solutions/trials (%%)\n");
			fprintf(p_file,"%15.5f	%20.5f	%15d	%15d	%20.5f \n", rmsd_avg_solution, rmsd_avg_gen, counter_solution_rmsd, counter_per_gen, ratio);
			fprintf(p_file,"___________________________________________________________________________________________________\n");
			fprintf(p_file,"\n");
			fprintf(p_file,"\n");
			fprintf(p_file,"MON's TIMERS \n");
			fprintf(p_file,"___________________________________________________________________________________________________\n");
			fprintf(p_file,"%8s %8s %8s %6s %2s %6s %2s %6s %2s %6s %2s %6s %2s\n","total","sum","loop[ms]","init","%","begin","%","inside","%","inf","%","RCD","%");
			fprintf(p_file,"%8.4lf %8.4lf %8.4lf %6.4lf %2.0f %6.4lf %2.0f %6.4lf %2.0f %6.4lf %2.0f %6.4lf %2.0f\n",t_time,timer_0+timer_1+timer_2, 1000*t_time/cycles,
					timer_0,100*timer_0/t_time,timer_1,100*timer_1/t_time,timer_2,100*timer_2/t_time,timer_3,100*timer_3/t_time,timer_4,100*timer_4/t_time);
			fprintf(p_file,"___________________________________________________________________________________________________\n");
			fprintf(p_file,"\n");
			//***************************************************************************************************************************************
			fprintf(p_file,"%s\n","DATA PROTEIN LOOP");		/*Protein loop details	*/
			fprintf(p_file,"%-50s %15s\n","PDB abbreviation protein:",current_pdb_name);
			fprintf(p_file,"%-50s %15i\n","Start residue loop itself",start+1);
			fprintf(p_file,"%-50s %15i\n","End residue loop itself",end-1);
			fprintf(p_file,"%-50s %15i\n","Begin anchor residue",start);
			fprintf(p_file,"%-50s %15i\n","End anchor residue",end);
			fprintf(p_file,"%-64s %c\n","Chain ID",in_chain);
			fprintf(p_file,"\n");

			fprintf(p_file,"\n");					/*Parameter settings	*/
			fprintf(p_file,"%s\n","PARAMETER SETTINGS");
			fprintf(p_file,"%-50s %15i\n","Number closures/pdb",cycles);
			fprintf(p_file,"%-50s %15i\n","Maximum iteration count/closure",max_path);
			fprintf(p_file,"%-50s %15i\n","Maximum iteration count/perturbation cycle",max_perturbation);
			fprintf(p_file,"%-50s %15i\n","Index for saving a trajectory path",number_path);
			fprintf(p_file,"\n");
			fprintf(p_file,"\n%s\n","SLOW TIMER");	/*Timing and statistics*/
			ftimer->printsize(p_file);
			ftimer->print2("FAST TIMER STATISTICS",2,p_file);
			fcounter->print2("COUNTER STATISTICS",3,p_file);
			frmsd->print2("RMSD STATISTICS LOOP SOLUTIONS",4,p_file);
			energies->print2("BUMP (PD2) ENERGY STATISTICS LOOP SOLUTIONS",4,p_file);
			energies_extra->print2("BUMP (PD2) EXTRA ENERGY STATISTICS LOOP SOLUTIONS",4,p_file);
			energies_intra->print2("BUMP (PD2) INTRA ENERGY STATISTICS LOOP SOLUTIONS",4,p_file);
			scores->print2("Score statistics for LOOP SOLUTIONS",4,p_file);
			if(em_switch && (emodel == 1 || emodel == 5 || emodel == 56 || emodel == 6) ) // Only do ccorrs stuff when energies computation was also required
				ccorrs->print2("Cross-correlations with EM-maps",4,p_file);
			fprintf(p_file,"\n");
			fprintf(p_file,"%-50s %15d/%d\n","Solutions/trials",k,totalk);
			fprintf(p_file,"%-50s %15.3f\n","Convergence fraction",  double(k)/ double(totalk) );
			fprintf(p_file,"%s\n"," ");

			fprintf(p_file,"\n%s\n","COUNTER ANALYSIS");		/*Counter analysis	*/
			fprintf(p_file,"%-50s %15li\n","Total iteration count",counter_cycles_total);
			fprintf(p_file,"%-50s %15li\n","Average count/LC", counter_cycles_total/cycles);
			fprintf(p_file,"%s\n"," ");

			fprintf(p_file,"\n%s\n","POST-FILTER CLASH ANALYSIS");	/*Post-filter clash analysis*/
			fprintf(p_file,"%-50s %15i\n","Loop closure solutions with clashes remaining",counter_clash);
			fprintf(p_file,"%-50s %15.3f\n","Fraction post-filter clashes", double (counter_clash) /double(cycles) );
			fprintf(p_file,"%s\n"," ");
			fclose(p_file);
			// PABLO 2020 bug...
			anchor_rmsd->array[i_pdb] = anchor_rmsd_avg / (double)cycles; // Store the final RMSD for statistics

			if((emodel == 1 || emodel == 5 || emodel == 56 || emodel == 6 || em_switch ) && loco_best > 0)
				cycles = loco_best+loco_notbest; // required for "best energy" stuff...

			tot_rmsd_avg->array[i_pdb]	= frmsd->average; 	// store averages of all PDB's
			tot_rmsd_sig->array[i_pdb]	= frmsd->deviation; // store sigmas of all PDB's
			tot_rmsd_min->array[i_pdb]	= frmsd->minimum; 	// store minima of all PDB's
			tot_rmsd_05->array[i_pdb]	= 100 * frmsd->nless(0.5)/(double)cycles; // store <0.5A of all PDB's
			tot_rmsd_10->array[i_pdb]	= 100 * frmsd->nless(1.0)/(double)cycles; // store <1.0A of all PDB's
			tot_rmsd_15->array[i_pdb]	= 100 * frmsd->nless(1.5)/(double)cycles; // store <1.5A of all PDB's
			tot_rmsd_20->array[i_pdb]	= 100 * frmsd->nless(2.0)/(double)cycles; // store <2.0A of all PDB's
			tot_rmsd_25->array[i_pdb]	= 100 * frmsd->nless(2.5)/(double)cycles; // store <2.5A of all PDB's
			tot_rmsd_30->array[i_pdb]	= 100 * frmsd->nless(3.0)/(double)cycles; // store <3.0A of all PDB's
			tot_counter->array[i_pdb]	= (float) counter_cycles_total/ (float) totalk; // Array to store total counters of all PDB's
			tot_time->array[i_pdb]		= timer1->elapsed(); // Array to store total times of all PDB's
			tot_con->array[i_pdb]		= 100 * (double) k/ (double) totalk; // Array to store Convergence rates of all PDB's
			tot_trialF->array[i_pdb]	= 100 * trials[0]/ (double) totalk; // Array to store Trial Forward rates of all PDB's
			tot_trialB->array[i_pdb]	= 100 * trials[1]/ (double) totalk; // Array to store Trial Backward rates of all PDB's
			tot_sucF->array[i_pdb]		= 100 * (trials[0]-fails[0])/(double) k; // Array to store Success Forward rates of all PDB's
			tot_sucB->array[i_pdb]		= 100 * (trials[1]-fails[1])/(double) k; // Array to store Success Backward rates of all PDB's


			intra_clashes->array[i_pdb] = 100 * intra_clashes_avg / (double)cycles; // Store the number of times there is an intra-clash (intra-loop)
			extra_clashes->array[i_pdb] = 100 * extra_clashes_avg / (double)cycles; // Store the number of times there is an extra-clash (clashes with map)

			// Save results single PDB's in results.txt file
			if( !(summary_file = summary_file =	fopen(res_path,"a+")) )
			{
				fprintf(stderr,"Sorry, unable to write file in a+ mode: %s --> Forcing exit!!!\n",res_path);
				exit(1);
			}

			fprintf(summary_file,"%-12s %6.3lf %5.3lf %6.3lf %6.2lf %6.2lf %6.2lf %6.2lf %6.2lf %6.2lf %5.0f %7.2lf %3.0f %3.0f/%-3.0f %3.0f/%-3.0f %5.3f %3.0f %3.0f %s\n",
					current_pdb_name,
					tot_rmsd_avg->array[i_pdb], 	// store averages of all PDB's
					tot_rmsd_sig->array[i_pdb], // store sigmas of all PDB's
					tot_rmsd_min->array[i_pdb], 	// store minima of all PDB's
					tot_rmsd_05->array[i_pdb], // store <0.5A of all PDB's
					tot_rmsd_10->array[i_pdb], // store <1.0A of all PDB's
					tot_rmsd_15->array[i_pdb], // store <1.5A of all PDB's
					tot_rmsd_20->array[i_pdb], // store <2.0A of all PDB's
					tot_rmsd_25->array[i_pdb], // store <2.5A of all PDB's
					tot_rmsd_30->array[i_pdb], // store <3.0A of all PDB's
					tot_counter->array[i_pdb], // Array to store total counters of all PDB's
					tot_time->array[i_pdb], // Array to store total times of all PDB's
					tot_con->array[i_pdb], // Array to store Convergence rates of all PDB's
					tot_trialF->array[i_pdb], // Array to store Trial Forward rates of all PDB's
					tot_trialB->array[i_pdb], // Array to store Trial Backward rates of all PDB's
					tot_sucF->array[i_pdb], // Array to store Success Forward rates of all PDB's
					tot_sucB->array[i_pdb], // Array to store Success Backward rates of all PDB's
					anchor_rmsd->array[i_pdb], // Store the final RMSD for statistics
					intra_clashes->array[i_pdb], // Number of times there is an intra-clash
					extra_clashes->array[i_pdb], // Number of times there is an extra-clash (clashes with map)
					pdbfile[i_pdb].sequence // Sequence
			);
			fclose(summary_file);

			if (i_pdb == ind - 1)		/*Calculate averages over all PDB loops*/
			{
				tot_rmsd_avg->compute(); // store averages of all PDB's
				tot_rmsd_sig->compute(); // store sigmas of all PDB's
				tot_rmsd_min->compute(); // store minima of all PDB's
				tot_rmsd_05->compute(); // store <0.5A of all PDB's
				tot_rmsd_10->compute(); // store <1.0A of all PDB's
				tot_rmsd_15->compute(); // store <1.5A of all PDB's
				tot_rmsd_20->compute(); // store <2.0A of all PDB's
				tot_rmsd_25->compute(); // store <2.5A of all PDB's
				tot_rmsd_30->compute(); // store <3.0A of all PDB's
				tot_counter->compute(); // store total counters of all PDB's
				tot_time->compute(); // store total times of all PDB's
				tot_con->compute(); // store convergence rates of all PDB's
				tot_trialF->compute(); // store Trial Forward rates of all PDB's
				tot_trialB->compute(); // store Trial Backward rates of all PDB's
				tot_sucF->compute(); // store Success Forward rates of all PDB's
				tot_sucB->compute(); // store Success Backward rates of all PDB's
				anchor_rmsd->compute(); // Store the final RMSD for statistics
				intra_clashes->compute(); // Store the number of times there is an intra-clash
				extra_clashes->compute(); // Store the number of times there is an extra-clash (clashes with map)

				summary_file =	fopen(res_path,"a+");
				fprintf(summary_file,"\n%-12s %6.3lf %5.3lf %6.3lf %6.2lf %6.2lf %6.2lf %6.2lf %6.2lf %6.2lf %5.0f %7.2lf %3.0f %3.0f/%-3.0f %3.0f/%-3.0f %5.3f %3.0f %3.0f\n",
						"Average",
						tot_rmsd_avg->average, // averages of all PDB's
						tot_rmsd_sig->average, // sigmas of all PDB's
						tot_rmsd_min->average, // minima of all PDB's
						tot_rmsd_05->average, // <0.5A of all PDB's
						tot_rmsd_10->average, // <1.0A of all PDB's
						tot_rmsd_15->average, // <1.5A of all PDB's
						tot_rmsd_20->average, // <2.0A of all PDB's
						tot_rmsd_25->average, // <2.5A of all PDB's
						tot_rmsd_30->average, // <3.0A of all PDB's
						tot_counter->average, // total counters of all PDB's
						tot_time->average, // total times of all PDB's
						tot_con->average, // convergence rates of all PDB's
						tot_trialF->average, // trial Forward rates of all PDB's
						tot_trialB->average, // trial Backward rates of all PDB's
						tot_sucF->average, // success Forward rates of all PDB's
						tot_sucB->average, // success Backward rates of all PDB's
						anchor_rmsd->average, // final RMSD between target and closed anchor
						intra_clashes->average, // Number of times there is an intra-clash (intra-loop clashes)
						extra_clashes->average // Number of times there is an extra-clash (clashes with map)
				);

				fprintf(summary_file,"%-12s %6.3lf %5.3lf %6.3lf %6.2lf %6.2lf %6.2lf %6.2lf %6.2lf %6.2lf %5.0f %7.2lf %3.0f %3.0f/%-3.0f %3.0f/%-3.0f %5.3f %3.0f %3.0f\n",
						"Sigma",
						tot_rmsd_avg->deviation, // deviations of all PDB's
						tot_rmsd_sig->deviation, // sigmas of all PDB's
						tot_rmsd_min->deviation, // minima of all PDB's
						tot_rmsd_05->deviation, // <0.5A of all PDB's
						tot_rmsd_10->deviation, // <1.0A of all PDB's
						tot_rmsd_15->deviation, // <1.5A of all PDB's
						tot_rmsd_20->deviation, // <2.0A of all PDB's
						tot_rmsd_25->deviation, // <2.5A of all PDB's
						tot_rmsd_30->deviation, // <3.0A of all PDB's
						tot_counter->deviation, // total counters of all PDB's
						tot_time->deviation, // total times of all PDB's
						tot_con->deviation, // convergence rates of all PDB's
						tot_trialF->deviation, // trial Forward rates of all PDB's
						tot_trialB->deviation, // trial Backward rates of all PDB's
						tot_sucF->deviation, // success Forward rates of all PDB's
						tot_sucB->deviation, // success Backward rates of all PDB's
						anchor_rmsd->deviation, // final RMSD between target and closed anchor
						intra_clashes->deviation, // Number of times there is an intra-clash
						extra_clashes->deviation // Number of times there is an extra-clash (clashes with map)
				);

				fprintf(summary_file,"\n\n");
				fclose(summary_file);
			}
		}

		time_t rawtime_path;	/*getting time data this PDB run*/
		struct tm * timeinfo_f_path;
		time ( &rawtime_path );
		timeinfo_f_path = localtime ( &rawtime_path );

		char 	pdb_name_path[5];	/*to store directories and filenames*/
		char 	file_name_output_path[LENGTH_LINE];									//	to store name of the output solutions pdb-file
		char 	file_name_output_path3[LENGTH_LINE];								//	to store name of the output solutions txt-file
		char 	file_name_output_path4[LENGTH_LINE];								//	to store name of the output initial loops pdb-file

		strcpy(tmp_path,folder_later);
		strcat(tmp_path,"/");

		if (pid==0 && save_movie)	/*save one closure trajectory*/
		{

			strcpy(tmp_path,folder_later);
			strcat(tmp_path,"/");
			strcat(tmp_path,current_pdb_name);
			strcpy(file_name_output_path,tmp_path);
			strcat(file_name_output_path,"_traj.pdb");
			multiPDB_write(path->el,counter_path,nr_loop+2,file_name_output_path);
		}

		if (pid==0 && switch_solutions == 1)	/*save all found solutions*/
		{
			strcpy(tmp_path,folder_later);
			strcat(tmp_path,"/");
			strcat(tmp_path,current_pdb_name);
			printf("rcd> Saved files in %s\n",folder_later);

			if(save_closed)
			{
				strcpy(file_name_output_path,tmp_path);	/*PDB format*/
				strcat(file_name_output_path,"_closed.pdb");
				multiPDB_cg(co_sol->el,p1co_sol->el,cycles,nr_loop+2,file_name_output_path,secu,start,cg_mode,in_chain);	    /*full backbone: PDB format?*/ /*All solutions*/
			}
			if(save_initial)
			{
				strcpy(file_name_output_path4,tmp_path);	/*TXT format*/
				strcat(file_name_output_path4,"_initial.pdb");
				multiPDB_cg(co_soli->el,p1co_soli->el,cycles,nr_loop+2,file_name_output_path4,secu,start,cg_mode,in_chain);	    /*full backbone: PDB format?*/ /*All solutions*/
			}
		}

		// Save first (or native) loop. It can be inserted into input structure with "insertloop.pl" to feed Rosetta...
		if(pid==0 && mutate)
		{
			strcpy(file_name_output_path3,tmp_path);	/*TXT format*/
			strcat(file_name_output_path3,"_loop1.pdb");

			// Mon: Bug in RCD+ server (6/10/2016 detected)
			//      Warning! "loop1" will be used in RMSD computations by Rosetta!!! It should contain Native coordinates, otherwise RMSDs are crap!
			//      Solution: Always output the PDB coordinates (backbone+CB) when they're present. Then, Rosetta's RMSDs will be computed fine!
			if( native == 1 || present_seq == 1 || bench_rmsd ) // If loop is present in PDB
				multiPDB_cg(co_i->el,p1co_i->el,1,nr_loop+2,file_name_output_path3,secu,start,cg_mode,in_chain);	    /*full backbone: PDB format?*/ /*All solutions*/
			else
				multiPDB_cg(co_sol->el,p1co_sol->el,1,nr_loop+2,file_name_output_path3,secu,start,cg_mode,in_chain);	    /*full backbone: PDB format?*/ /*All solutions*/
		}

		// RASP repacking
		//		if(rasp_switch)
		//		{
		//			timer1->start();			         //	START TIMER
		//
		//			char myname[LENGTH_LINE];
		//
		//			strcpy(myname,folder_later);
		//			strcat(myname,"/");
		//			strcat(myname,current_pdb_name);
		//			strcat(myname,"_rasp.txt");
		//			fprintf(stderr,"RASP energy will be stored in: %s\n",myname);
		//			FILE *myf2 = fopen(myname,"w"); // Delete previous file
		//			fprintf(myf2,"#%4s %11s\n","","RASP Energy");
		//
		//			strcpy(myname,folder_later);
		//			strcat(myname,"/");
		//			strcat(myname,current_pdb_name);
		//			strcat(myname,"_rasp.pdb");
		//			FILE *myf = fopen(myname,"w"); // Delete previous file
		//			fclose(myf);
		//
		//			for(int i=0; i<cycles; i++)
		//			{
		//				// Update coordinates of atoms ("co") between "start" and "end" residues (RASP indexed).
		//				co2rasp(co_sol->el + i*nr_atoms_loop, rasp_in, rasp_start, rasp_end);
		//
		//				// Zhichao's repacking routine
		//				fprintf(myf2,"%5d %11f\n",i+1,rasp(rasp_in,rasp_out,rasp_mut));
		//
		//				// Mon's routine to append a segment [start,end] of some "rmodel" into a Multi-PDB file
		//				append_mdl(myname, rasp_out, i+1, start, end);
		//
		//				// Save the complete structure (protein + loop)
		//				//	char myname[100];
		//				//	sprintf(myname,"raspa_%02d_%04d.pdb",i_pdb+1,i+1);
		//				//	write_mdl(myname,rasp_out);
		//				//	fprintf(stderr,"RASP re-built loop sidechains in: %s\n",myname);
		//
		//				destroy_rmodel(rasp_out);
		//			}
		//			fclose(myf2);
		//			fprintf(stderr,"RASP re-built loop sidechains in: %s\n",myname);
		//			// Destroy RASP's rmodels
		//			destroy_rmodel(rasp_in);
		//			rasp_mut.clear();
		//
		//			timer1->stop();			         //	STOP FAST TIMER
		//			fprintf(stderr,"rcd> RASP repacking time= %f s\n",timer1->elapsed());
		//
		//			rasp_time += timer1->elapsed();
		//		}

		// New RASP repacking stuff...
		if(rasp_switch && pid==0) // I'm the Master
		{
			timer1->start();			         //	START TIMER

			char myname[LENGTH_LINE];

			//			strcpy(myname,folder_later);
			//			strcat(myname,"/");
			//			strcat(myname,current_pdb_name);
			//			strcat(myname,"_rasp.txt");
			//			fprintf(stderr,"RASP energy will be stored in: %s\n",myname);
			//			FILE *myf2 = fopen(myname,"w"); // Delete previous file
			//			fprintf(myf2,"#%4s %11s\n","","RASP Energy");

			strcpy(myname,folder_later);
			strcat(myname,"/");
			strcat(myname,current_pdb_name);
			strcat(myname,"_rasp.pdb");
			FILE *myf = fopen(myname,"w"); // Delete previous file
			fclose(myf);

			// Get parameters from RASP.ini file
			rasp->GetPar();

			// Set re-packable residues from string
			rasp->GetSub(rasp_mut);

			rasp->GetTop();

			//			// Set repackeable residues from string
			//			rasp->GetSub(rasp_mut);
			// fprintf(stderr,"RASP cycles= %d\n",cycles);

			for(int i=0; i<cycles; i++)
			{

				// Update coordinates of atoms ("co") between "start" and "end" residues (RASP indexed).
				// fprintf(stderr,"i=%d Before co2rasp...\n",i);
				co2rasp(co_sol->el + i*nr_atoms_loop, rasp, rasp_start, rasp_end);
				// fprintf(stderr,"i=%d After co2rasp...\n",i);

				// Rebuild Oxygen atom in residues selected for repacking via "Oss"
				rebuildO(rasp);

				// Populate 4-atoms (N,CA,C,O) backbone "bcbn" (without any checking)
				rasp->Pdb2bcbn();

				// Compute all Phi/Psi dihedral angles from "bcbn"
				rasp->PhiPsi();

				// Load Dunbrack's data directly from file (bbdep11.bin) depending of current Phi and Psi...
				rasp->GetLib(); // Must be here...
				rasp->Build(); // Some initializations? Check whether this can go here...
				rasp->SelfEnergy();
				rasp->PairEnergy();
				rasp->Search();

				append_mdl(myname, rasp, i+1, rasp_start, rasp_end);

				rasp->bcbn.clear(); // Clear all Backbone
				rasp->sdcn.clear(); // Clear all side-chains
				eraseSC(rasp); // Erase Side-Chains of residues selected for repacking via "Oss" in "Pdb"

				//				exit(0);

				// Save the complete structure (protein + loop)
				//	char myname[100];
				//	sprintf(myname,"raspa_%02d_%04d.pdb",i_pdb+1,i+1);
				//	write_mdl(myname,rasp_out);
				//	fprintf(stderr,"RASP re-built loop sidechains in: %s\n",myname);

				// destroy_rmodel(rasp_out);

				indicator("rcd> RASP repacking ",i,cycles);
			}
			indicator("rcd> RASP repacking ",cycles,cycles);

			// fclose(myf2);
			fprintf(stderr,"RASP re-built loop sidechains in: %s\n",myname);

			// destroy_rmodel(rasp_in);
			rasp_mut.clear();

			// Destroy RASP's rmodels
			delete rasp;

			timer1->stop();			         //	STOP FAST TIMER
			fprintf(stderr,"rcd> RASP repacking time= %f s\n",timer1->elapsed());

			rasp_time += timer1->elapsed();
		}


		/*Print rmsd data to txt-file*/
		char 	file_name_output_path5[LENGTH_LINE];								//	to store name of the output solutions txt-file
		char    move_path5[LENGTH_LINE];

		if(pid==0)
		{
			fprintf(stderr,"rcd> Saving RMSD values\n");
			//tmp_path[0]=' ';
			file_name_output_path5[0]=' ';
			strcpy(tmp_path,folder_later);
			strcat(tmp_path,"/");
			strcat(tmp_path,current_pdb_name);
			strcpy(file_name_output_path5,tmp_path);	/*PDB format*/
			strcat(file_name_output_path5,"_rmsd.txt");

			FILE *file_rmsd;
			file_rmsd = fopen(file_name_output_path5,"w"); //:REMARK: is this correct output filename


			bool local_energy_switch =  emodel == 1 || emodel == 5 || emodel == 56 || emodel == 6 || em_switch; // Energy and/or Cross-corr were requested
			bool local_corr_switch =  (emodel == 1 || emodel == 5 || emodel == 56 || emodel == 6) && em_switch; // Both requested at the same time ???

			// Print header of scores file
			fprintf(file_rmsd,"%-8s %8s %8s %8s %8s","# Loop","RMSD[A]","Bump","BumpEx","BumpIn");
			//			if(local_energy_switch)
			//				fprintf(file_rmsd," %10s","Score1");
			if(emodel == 1 || emodel == 5 || emodel == 56 || emodel == 6) // energies computation was required
				fprintf(file_rmsd," %10s","Energy");
			if(em_switch) // EM-maps correlation war required
				fprintf(file_rmsd," %10s","1-Corr.");
			fprintf(file_rmsd,"\n#%7s %8s %8.3f %8s %8s","native","0.0",native_energy,"-","-");
			//			if(local_energy_switch)
			//				fprintf(file_rmsd," %10.3f",native_loco);
			if(emodel == 1 || emodel == 5 || emodel == 56 || emodel == 6 || em_switch) // energies or cross-correlation computation were required
				fprintf(file_rmsd," %10.3f",native_loco);
			//			else if(em_switch)
			//				fprintf(file_rmsd," %10s","NoComputed");

			fprintf(file_rmsd,"\n");


			//			if(emodel == 1 || emodel == 5 || em_switch)
			//			{
			if(loco_best > 0)
				cycles = loco_best+loco_notbest; // Warning, from now on "cycles" is different! (sorry, but this makes things much easier...)

			//				for (int u = 0; u < cycles; u++)
			//					fprintf(file_rmsd,"%8i %8.3lf %8.3f %8.3lf %8.3f %10.5f\n",u,
			//							frmsd->array[u],energies->array[u],energies_extra->array[u],energies_intra->array[u],scores->array[u]);

			for (int u = 0; u < cycles; u++)
			{
				fprintf(file_rmsd,"%8i %8.3lf %8.3f %8.3lf %8.3f",u,frmsd->array[u],
						energies->array[u],energies_extra->array[u],energies_intra->array[u]);

				if(local_energy_switch) // Energy or Cross-corr were requested
					fprintf(file_rmsd," %10.5f",scores->array[u]);

				if(local_corr_switch) // Both requested at the same time
					fprintf(file_rmsd," %10.5f",ccorrs->array[u]);

				fprintf(file_rmsd,"\n");
			}


			cycles = cycles0;
			//			}
			//			else
			//				for (int u = 0; u < cycles; u++)
			//					fprintf(file_rmsd,"%8i %8.3lf %8.3f %8.3lf %8.3f\n",u,
			//							frmsd->array[u],energies->array[u],energies_extra->array[u],energies_intra->array[u]);

			fclose(file_rmsd);

			// Mon: Experimental plot
			if(logbootstrap_switch)
			{
				int nreplica = 100; // Number of replicas
				int npoints = 100; // Number of points
				int reduceby = 4;
				int ndelta = (cycles/reduceby)/npoints; // Number of samples increment
				float *work = (float *) malloc( sizeof(float) * cycles );
				float *reps = (float *) malloc( sizeof(float) * nreplica );
				float *log_minrmsds = (float *) malloc( sizeof(float) * (cycles / ndelta) );
				float *log_ns = (float *) malloc( sizeof(float) * (cycles / ndelta) );
				float avg,sig;

				file_name_output_path5[0]='\0';
				strcpy(tmp_path,folder_later);
				strcat(tmp_path,"/");
				strcat(tmp_path,current_pdb_name);
				strcpy(file_name_output_path5,tmp_path);	/*PDB format*/
				strcat(file_name_output_path5,"_logboot.txt");

				FILE *file_logboot;
				file_logboot = fopen(file_name_output_path5,"w"); //:REMARK: is this correct output filename
				fprintf(file_logboot,"%6s %7s +- %6s  %6s %6s\n","#loops", "MinRMSD", "sigma", "log(N)", "log(min)");

				int i = 0;
				for(int n=ndelta; n <= cycles/reduceby; n += ndelta)
				{

					for(int r=0; r<nreplica; r++)
					{
						randomSample(frmsd->array,cycles,work,n); // Load "n" random samples into "work" from input array
						//						for (int i = 0; i < n; i++)
						//							fprintf(stderr," %6.3f",work[i]);
						reps[r] = getmin(work,n);
						//						fprintf(stderr,"%4d %6.3f\n",n,getmin(work,n));
					}
					// Compute the average of a float array "data" of "ndata" elements
					avg = average(reps, nreplica);

					// Compute the sigma of a float array "data" of "ndata" elements. Optionally, the average can be provided if already known.
					sig = sigma(reps, nreplica, avg);

					log_ns[i] = logf( (float) n );
					log_minrmsds[i] = logf( avg );

					fprintf(file_logboot,"%6d  %6.3f +- %6.3f  %6.3f %6.3f\n",n, avg, sig, logf((float)n), logf(avg));
					i++;
				}
				// Linear regression: y = b + m·x  (correlation = r)
				float b,m,r;
				linear_regression(log_ns+npoints/2, log_minrmsds+npoints/2, npoints/2, &b, &m, &r);
				fprintf(file_logboot,"# Linear regression: y = %11.6f %11.6f · x   r= %11.6f\n", b, m ,r);

				// log(min_RMSD) = b + m·log(N)
				// log(N) = (log(min_RMSD) - b) / m
				// N = e^( (log(min_RMSD) - b) / m )

				fprintf(file_logboot, "#");

				for(float r=0.5; r<3.0; r += 0.1 )
					fprintf(file_logboot, " %d", (int) expf( (logf(r)-b) / m) );
				fprintf(file_logboot, "\n");

				fclose(file_logboot);
				free(work);
				free(reps);
				free(log_minrmsds);
				free(log_ns);
				// exit(0);
			}
		}



		if(load_ramap)
		{
			free(iaa);
			iaa=NULL;
			for(int i=0; i<nr_loop+2; i++)
				free(maps[i]);
			free(maps);
			maps=NULL;
		}

		if (ramachandran_check_aux == 2 || ramachandran_inside_aux == 2)
		{
			for(int i=0; i<nr_loop+2; i++)
				free(accs[i]);
			free(accs);
			accs=NULL;
		}

		if(save_movie)
		{
			delete rmsdprofile;
			delete counterprofile;
		}

#ifdef MPI_RELEASE
		if(pid != 0)
			cycles = cycles0 / np;
		else
			cycles = cycles0; // required for "best energy" stuff...
#else
		// Mon bug (14/7/2015) (for non-mpi implementation)
		//		if(loco_switch && loco_best > 0)
		if((emodel == 1 || emodel == 5 || emodel == 56 || emodel == 6 || em_switch ) && loco_best > 0)
			cycles = cycles0; // required for "best energy" stuff...
#endif

		// MON: Warning, "loop_coord should be freed!!!"
		//		if(loco_switch)
		if(emodel == 1)
		{
			free(loop_coord);
			delete energies;
			delete energies_intra;
			delete energies_extra;
			delete scores;
			energies 		= new Stat(cycles);    //  Stat object for Energies
			energies_intra	= new Stat(cycles);    //  Stat object for Intra-Energies
			energies_extra	= new Stat(cycles);    //  Stat object for Extra-Energies
			scores		= new Stat(cycles);    //  Stat object for LOCO-ICOS energies
		}

		if(emodel == 5 || emodel == 56)
		{
			free(loop_coord);
			delete energies;
			delete energies_intra;
			delete energies_extra;
			delete scores;
			energies 		= new Stat(cycles);    //  Stat object for Energies
			energies_intra	= new Stat(cycles);    //  Stat object for Intra-Energies
			energies_extra	= new Stat(cycles);    //  Stat object for Extra-Energies
			scores		= new Stat(cycles);    //  Stat object for LOCO-ICOS energies
			free(resnums);
			free(reschainids);
			free(resnumsloop);
			free(reschainidsloop);
		}

		if(emodel == 6)
		{
			delete energies;
			delete energies_intra;
			delete energies_extra;
			delete scores;
			energies 		= new Stat(cycles);    //  Stat object for Energies
			energies_intra	= new Stat(cycles);    //  Stat object for Intra-Energies
			energies_extra	= new Stat(cycles);    //  Stat object for Extra-Energies
			scores		= new Stat(cycles);    //  Stat object for LOCO-ICOS energies
			//			free(resnums);
			//			free(reschainids);
			//			free(resnumsloop);
			//			free(reschainidsloop);
		}

		if(em_switch)
		{
			delete energies;
			delete energies_intra;
			delete energies_extra;
			delete scores;
			energies 		= new Stat(cycles);    //  Stat object for Energies
			energies_intra	= new Stat(cycles);    //  Stat object for Intra-Energies
			energies_extra	= new Stat(cycles);    //  Stat object for Extra-Energies
			scores		= new Stat(cycles);    //  Stat object for LOCO-ICOS energies

			delete(workmap); // Delete the cropped and masked EM-map
			delete(loopmap); // Delete the loop simulated map
			free(workmask); // Free integers mask for fast cross-correlation computation
		}

		free(bond_length2);
		free(valence_angle2);
		free(dihedral_angle2);

		free(dihedral_angle_i);
		free(dihedral_angle_f);
		free(dihedral_angle_r);
		free(dihedral_angle_aux);

		free(bond_length_back);
		free(valence_angle_back);

		free(bond_length);
		free(valence_angle);
		free(dihedral_angle);





		//free(length_sol);
		//free(valence_sol);
		//free(dihedrals_sol);

		delete part;			/*deallocate matrix objects*/
		delete part1;
		delete part2;

		delete co_sol;
		delete p1co_sol;
		delete p2co_sol;

		delete co_soli;
		delete p1co_soli;
		delete p2co_soli;

		delete co;
		delete p1co;
		delete p2co;

		delete co_i;
		delete p1co_i;
		delete p2co_i;

		delete co_f;
		delete p1co_f;
		delete p2co_f;

		delete co_r;
		delete p1co_r;
		delete p2co_r;

		delete coaux;
		delete p1coaux;
		delete p2coaux;

		delete co2;
		delete p1co2;
		delete p2co2;

		delete conoloop;

		delete path;

		delete s1co;

		delete[] residuemarker;

		//EXPERIMENTAL CODE
		delete[] marker;
		delete[] marker2;
		delete[] psimarker;

		delete[] marker_i;
		delete[] marker2_i;
		delete[] psimarker_i;

		delete[] marker_r;
		delete[] marker2_r;
		delete[] psimarker_r;

		delete[] prearray1;
		delete[] prearray2;
		delete[] prearray3;
		delete[] prearray4;
		delete[] prearray5;
		delete[] prearray6;



	}	// INFO: Single PDB run completed
	//exit(0);

#ifdef MPI_RELEASE
	MPI_Finalize();
#endif


	// ALEX: Info of RMSD avg in generation
	rmsd_pdbs_avg_gen = rmsd_pdbs_gen/total_gen*1.0;
	rmsd_pdbs_avg_sol = rmsd_pdbs_sol/total_sols*1.0;
	//	fprintf(stderr,"rmsd_pdbs_avg_gen = %d 	rmsd_pdbs_avg_sol = %d \n",rmsd_pdbs_avg_gen,rmsd_pdbs_avg_sol);

	//***************************************
	/*Saving date of the run*/		/*REMARK: MUST BE SAME AS DIRECTORY!!!!*/
	char	date_run[LENGTH_LINE];
	time_t 	date;
	struct tm * timeinfo_date;
	time ( &date );
	timeinfo_date = localtime ( &date );
	strftime(date_run,80,"%d %b %Y %H:%M:%S",timeinfo_date);	/*seconds added for same file*/


if (pid==0 && switch_file==1)	//	Write parameter settings to summary report
	{
		summary_file =	fopen(res_path,"a+");
		fprintf(summary_file," \n");
		fprintf(summary_file,"%s\n","RUN SPECIFICATIONS");
		fprintf(summary_file,"\n\n");
		fprintf(summary_file,"%-15s %50s\n","Date of run",date_run);
		fprintf(summary_file,"%-35s %30s\n","Input file",testfile);
		fprintf(summary_file,"%-50s %15i\n","Number of PDBs tested",ind);
		fprintf(summary_file,"%-50s %15i\n","Number closures/pdb",cycles);
		fprintf(summary_file,"%-50s %15i\n","Energy model",emodel);
		fprintf(summary_file,"%-50s %15i\n","Energy model",cg_mode);
		fprintf(summary_file,"%-50s %15i\n","Maximum iteration count/closure",max_path);
		fprintf(summary_file,"%-50s %15i\n","Maximum iteration count/perturbation cycle",max_perturbation);
		fprintf(summary_file,"%-50s %15i\n","Index for saving a trajectory path",number_path);
		fprintf(summary_file,"%-50s %15i\n","Intra-anchor degree of freedom:",switch_dof);
		fprintf(summary_file,"%-50s %15i\n","Using native loop internal coordinates (ics)",native);
		fprintf(summary_file,"%-50s %15i\n","Randomizing ics each closure ",randomize);
		fprintf(summary_file,"%-50s %15i\n","Randomizing bond lengths ",randbond);
		fprintf(summary_file,"%-50s %15i\n","Randomizing valence angles ",randval);
		fprintf(summary_file,"%-50s %15.3f\n","Omega sigma",omega_sigma);
		fprintf(summary_file,"%-50s %15i\n","Reverse loop closure ",reverse_input);
		fprintf(summary_file,"%-50s %15i\n","Switching direction mode:",switching_mode); // not operative
		fprintf(summary_file,"%-50s %15.3f\n","RMSD criterium for convergence",rmsd_crit);
		fprintf(summary_file,"%-50s %15.3f\n","Step size grid",stepsize);
		fprintf(summary_file,"%-50s %15.3f\n","Softness grid spheres",softness);
		fprintf(summary_file,"%-50s %15i\n","Number cutoff atoms loop in grid checking",cutoff);
		fprintf(summary_file,"%-50s %15i\n","Backbone grid",grid_backbone);
		fprintf(summary_file,"%-50s %15i\n","Sequential initialization on the grid",seq_initgrid);
		fprintf(summary_file,"%-50s %2.1i\n","Ramachandran initialization (dunbrack = 2 classic = 1; off = 0)",ramachandran_check);
		fprintf(summary_file,"%-50s %10i\n","Ramachandran inside (dunbrack = 2 classic = 1; off = 0)",ramachandran_inside);
		fprintf(summary_file,"%-50s %15.5f\n","Inside check threshold", rama_thr);
		fprintf(summary_file,"%-50s %15i\n","Ignore intra-chain contact shift",atomshift);
		fprintf(summary_file,"%-50s %15.5f\n","Post check threshold", postcheck);
		fprintf(summary_file,"%-50s %15i\n","Self-loop detection (on = 1; off = 0)",self_clash);
		fprintf(summary_file,"%-50s %15i\n","Post-algorithm grid check  (on = 1; off = 0)",post_check);
		fprintf(summary_file,"%-50s %15i\n","Exact_rama (on = 1; off = 0)",exact_rama);
		fprintf(summary_file,"%-50s %15i\n","Selected % of loops (loco_best)",loco_best);
		fprintf(summary_file,"%-50s %15i\n","Extra selected (loco_notbest)",loco_notbest);
		fprintf(summary_file,"%-50s %15i\n","Anchor residues to be relaxed",nrama);

		if (h3_switch)  // Set "true" to activate H3 Rama potential
		{
		fprintf(summary_file,"%-35s %30s\n","Kink maps file",h3file);
		fprintf(summary_file,"%-50s %15.5f\n","Rama cutoff inter-kink",rama_thr);
		fprintf(summary_file,"%-50s %15.5f\n","Rama cutoff anchors",rama_thrA);
		fprintf(summary_file,"%-50s %15.5f\n","Rama cutoff kink",rama_thrK);
		fprintf(summary_file,"%-50s %15.5f\n","Merge ratio kink",kink_merge);
		fprintf(summary_file,"%-50s %15.5f\n","Merge ratio anchors",kink_mergeA);
		fprintf(summary_file,"%-50s %15.5f\n","Merge ratio inter-kink",nokink_merge);
		}
		// ALEX
		fprintf(summary_file,"%-50s %15.3f\n","Avg RMSD (all pdbs) begin (trials)",rmsd_pdbs_avg_gen);
		fprintf(summary_file,"%-50s %15.3f\n","Avg RMSD (all pdbs) begin solution",rmsd_pdbs_avg_sol);
		//*********************************************************************************************************
		fprintf(summary_file,"%-50s %15i\n","Include hydrogens in protein map (on = 1; off = 0)",include_hydrogens);
		fprintf(summary_file,"\n");
		fclose(summary_file);

	}

	//	DEALLOCATION OUTSIDE PDB's LOOP --- OBJECTS NOT DEFINED BY PDB LOOP
	delete[] testfile;

	// delete timer1;
	delete ftimer;	//	DELETING STAT OBJECTS						Destructors of the form matrix2->Matrix::~Matrix(); are called when executing delete
	delete fcounter;
	delete frmsd;
	delete dih_found;

	delete matrix2;	//	DELETING MATRIX OBJECTS						Destructors of the form matrix2->Matrix::~Matrix(); are called when executing delete
	delete matrix;
	delete polpoint;
	delete anchorpoint;
	delete result1;
	delete result2;
	delete result3;
	delete anchor_end;
	delete anchor_end_i;
	delete anchor_end_r;

	delete[] angle_phi;
	delete[] angle_psi;

	free(pdbfile);		//:TODO: ASK WHETHER CORRECTLY DEALLOCATED?
	if(pid==0)
		printf("rcd> RCD execution completed\nrcd>\n");

	return 0;
}

// Compute the derivatives matrix ("dr/d_theta") of the positions (r) of the Ct-anchor atoms (N,CA,C) w.r.t. all dihedrals (theta)
//  der --> The derivatives matrix (#theta x 9 --> 9= 3 anchor atoms x 3 x,y,z coordinates). Set der=NULL to allocate memory.
//  co --> The N,CA,C coordinates (or C,CA,N in backward)
//  num_res --> number of residues including anchors
//  reverse --> (OPTIONAL) =0 forward (default), =1 backward
int der_anchor(double **p_der, double **co, int num_res, int reverse)
{
	bool debug = false; // set "true" to dump debug info
	double *der;
	double e[3],y[3],r[3]; // Bond vector
	double norm;
	int size,a,j=0;
	size = 9; // Number of constraints (3-atoms x 3 coordinates)

	if(*p_der==NULL)
	{
		if( !(der = (double *) malloc( sizeof(double) * 2*num_res * 9)) )
		{
			fprintf(stderr,"Msg(der_anchor): Unable to allocate memory for derivatives...\nForcing exit!\n");
			exit(1);
		}
		*p_der = der; // Output derivatives
	}

	for(int i=0; i < 2*num_res*9; i++)
		der[i] = 0.0; // Reset "der"

	if(reverse == 0) // Forward
	{
		// Screen residues
		for(int i=0; i < num_res; i++)
		{
			// Get some coordinates laying on the bond vector (e.g. CA coordinates)
			y[0] = co[3*i+1][0];
			y[1] = co[3*i+1][1];
			y[2] = co[3*i+1][2];

			//
			// PHI dihedral angle ------------------------------------------------
			//
			if(i > 0) // Not Nt-anchor residue
			{
				// Compute and normalize the bond vector:  e = r_CA - r_N
				e[0] = y[0] - co[3*i][0];
				e[1] = y[1] - co[3*i][1];
				e[2] = y[2] - co[3*i][2];
				// |e| = norm
				norm = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
				// |e| = 1
				e[0] /= norm;
				e[1] /= norm;
				e[2] /= norm;
				// Compute the derivatives w.r.t. the N,CA,C atoms of Ct-anchor
				for(int k = 0; k < 3; k++) // Screen N,CA,C atoms of Ct-anchor
				{
					a = 3*(num_res-1) + k; // Index of N,CA,C atoms of Ct-anchor
					// Coordinates of the N,CA,C atoms from Ct-anchor
					r[0] = co[a][0];
					r[1] = co[a][1];
					r[2] = co[a][2];
					// Derivative of the rotation of current atom ("r") around the vector "e" and center of rotation "y":   der = e x (r-y)
					// One constraint per row, and variables in columns.
					der[j*size+(3*k)]   = e[1]*(r[2]-y[2])-e[2]*(r[1]-y[1]);
					der[j*size+(3*k+1)] = e[2]*(r[0]-y[0])-e[0]*(r[2]-y[2]);
					der[j*size+(3*k+2)] = e[0]*(r[1]-y[1])-e[1]*(r[0]-y[0]);
				}
			}
			j++; // next dihedral

			//
			// PSI dihedral angle ------------------------------------------------
			//
			if(i < num_res-1) // Not Ct-anchor residue
			{
				// Compute and normalize the bond vector:  e = r_C - r_CA
				e[0] = co[3*i+2][0] - y[0];
				e[1] = co[3*i+2][1] - y[1];
				e[2] = co[3*i+2][2] - y[2];
				// |e| = norm
				norm = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
				// |e| = 1
				e[0] /= norm;
				e[1] /= norm;
				e[2] /= norm;
				// Compute the derivatives w.r.t. the N,CA,C atoms of Ct-anchor
				for(int k = 0; k < 3; k++) // Screen N,CA,C atoms of Ct-anchor
				{
					a = 3*(num_res-1) + k; // Index of N,CA,C atoms of Ct-anchor
					// Coordinates of the N,CA,C atoms from Ct-anchor
					r[0] = co[a][0];
					r[1] = co[a][1];
					r[2] = co[a][2];
					// Derivative of the rotation of current atom ("r") around the vector "e" and center of rotation "y":  der = e x (r-y)
					// One constraint per row, and variables in columns.
					der[j*size+(3*k)]   = e[1]*(r[2]-y[2])-e[2]*(r[1]-y[1]);
					der[j*size+(3*k+1)] = e[2]*(r[0]-y[0])-e[0]*(r[2]-y[2]);
					der[j*size+(3*k+2)] = e[0]*(r[1]-y[1])-e[1]*(r[0]-y[0]);
				}
			}
			j++; // Next dihedral
		}
	}
	else // Backward (atoms order in "co" array is reversed wrt forward)
	{
		// Screen residues
		for(int i=0; i < num_res; i++)
		{
			// Get some coordinates laying on the bond vector (e.g. CA coordinates)
			y[0] = co[3*i+1][0];
			y[1] = co[3*i+1][1];
			y[2] = co[3*i+1][2];

			//
			// PSI dihedral angle ------------------------------------------------
			//
			if(i > 0) // Not Ct-anchor residue (Ct-anchor's PSI does not move)
			{
				// Compute and normalize the bond vector:  e = r_C - r_CA
				e[0] = co[3*i][0] - y[0];
				e[1] = co[3*i][1] - y[1];
				e[2] = co[3*i][2] - y[2];
				// |e| = norm
				norm = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
				// |e| = 1
				e[0] /= norm;
				e[1] /= norm;
				e[2] /= norm;
				// Compute the derivatives w.r.t. the C,CA,N atoms of Nt-anchor
				for(int k = 0; k < 3; k++) // Screen C,CA,N atoms of Nt-anchor
				{
					a = 3*(num_res-1) + k; // Index of C,CA,N atoms of Nt-anchor
					// Coordinates of the C,CA,N atoms from Nt-anchor
					r[0] = co[a][0];
					r[1] = co[a][1];
					r[2] = co[a][2];
					// Derivative of the rotation of current atom ("r") around the vector "e" and center of rotation "y":  der = e x (r-y)
					// One constraint per row, and variables in columns.
					der[j*size+(3*k)]   = e[1]*(r[2]-y[2])-e[2]*(r[1]-y[1]);
					der[j*size+(3*k+1)] = e[2]*(r[0]-y[0])-e[0]*(r[2]-y[2]);
					der[j*size+(3*k+2)] = e[0]*(r[1]-y[1])-e[1]*(r[0]-y[0]);
				}
			}
			j++; // Next dihedral

			//
			// PHI dihedral angle ------------------------------------------------
			//
			if(i < num_res-1) // Not Nt-anchor residue (Nt-anchor's PHI does not move)
			{
				// Compute and normalize the bond vector:  e = r_CA - r_N
				e[0] = y[0] - co[3*i+2][0];
				e[1] = y[1] - co[3*i+2][1];
				e[2] = y[2] - co[3*i+2][2];
				// |e| = norm
				norm = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
				// |e| = 1
				e[0] /= norm;
				e[1] /= norm;
				e[2] /= norm;
				// Compute the derivatives w.r.t. the C,CA,N atoms of Nt-anchor
				for(int k = 0; k < 3; k++) // Screen C,CA,N atoms of Nt-anchor
				{
					a = 3*(num_res-1) + k; // Index of C,CA,N atoms of Nt-anchor
					// Coordinates of the C,CA,N atoms from Nt-anchor
					r[0] = co[a][0];
					r[1] = co[a][1];
					r[2] = co[a][2];
					// Derivative of the rotation of current atom ("r") around the vector "e" and center of rotation "y":   der = e x (r-y)
					// One constraint per row, and variables in columns.
					der[j*size+(3*k)]   = e[1]*(r[2]-y[2])-e[2]*(r[1]-y[1]);
					der[j*size+(3*k+1)] = e[2]*(r[0]-y[0])-e[0]*(r[2]-y[2]);
					der[j*size+(3*k+2)] = e[0]*(r[1]-y[1])-e[1]*(r[0]-y[0]);
				}
			}
			j++; // next dihedral
		}
	}

	if(debug)
	{
		fprintf(stderr,"Derivatives (A matrix) Column-major:\n");
		for(int i=0; i < size; i++)
			fprintf(stderr,"%9s%02d","X",i+1);
		fprintf(stderr,"\n");

		for(int k = 0; k < num_res*2; k++) // screen variables (dihedrals)
		{
			for(int j=0; j < size; j++) // Screen x,y,z components of N,CA,C atoms positions of Ct-anchor
				fprintf(stderr,"%11.6f",der[k*size+j]);
			fprintf(stderr,"\n");
		}
		fprintf(stderr,"\n");
	}

	return 0; // Everything fine...
}


// Compute the displacement vector B between the mobile Ct-anchor and the fixed one.
// It is the same for "forward" and for "backward" (reverse independent)
// B = [ b(i) ] = [ r_fix(i) - r_mob(i) ]     <--- Where "i" represents the N,CA,C atom indices of the Ct-anchor
//   b --> Displacement vector B (size= 9 x 1). Memory must be already allocated.
//   co --> The N,CA,C coordinates (or C,CA,N in backward) of the mobile loop
//   co_i --> The N,CA,C coordinates (or C,CA,N in backward) of the N/Ct-anchor (just 3 atoms coordinates)
//   num_res --> number of residues including anchors
double disp_b(double *b, double **co, double **co_i, int num_res)
{
	bool debug = false; // set "true" to dump debug info
	int last = 3*(num_res-1); // Index of the first atom of N/Ct-anchor
	double rmsd=0.0,temp;

	for(int i=0; i<3; i++)
	{
		temp=0.0;
		for(int x=0; x<3; x++)
		{
			b[i*3 + x] =  co_i[i][x] - co[last+i][x];
			temp += b[i*3 + x]*b[i*3 + x];
		}
		rmsd += sqrt(temp);
	}

	rmsd /= 3;

	if(debug)
	{
		fprintf(stderr,"Vector b (B matrix):\n");
		for(int i=0; i < 9; i++)
			fprintf(stderr,"%9.6f ",b[i]);
		fprintf(stderr,"\nRMSD= %f\n",rmsd);;
	}

	return rmsd; // everything was fine
}

// Get the maximum value
float getmax(float *a, int n)
{
	float max = a[0];
	for(int i=1; i<n; i++)
		if(a[i] > max)
			max = a[i];
	return max;
}

// Get the minimum value
float getmin(float *a, int n)
{
	float min = a[0];
	for(int i=1; i<n; i++)
		if(a[i] < min)
			min = a[i];
	return min;
}

// Load "n_out" random samples from "in" ("n_in" sized) into "out" array ("n_out" sized) without repetition
void randomSample(double *in, int n_in, float *out, int n_out)
{
	int *indices;
	indices = (int *) malloc( sizeof(int) * n_in );

	for(int i=0; i<n_in; i++)
		indices[i] = i; // Initialize array of indices

	int index;
	int i=0;
	while(i < n_out)
	{
		index = rg->IRandomX(0, n_in-1-i); // Get some random index (in-range and unused)
		out[i] = in[ indices[ index ] ]; // Output value without repetition
		indices[ index ] = indices[ n_in-1-i ]; // Replace the used index by an unused one
		i++; // effectively shorten the indices array
	}

	free(indices);
}


// Compute the average of a float array "data" of "ndata" elements
float average(float *data, int ndata)
{
	double dummy = 0.0;
	for(int i=0; i<ndata; i++)
		dummy += data[i];
	dummy /= ndata;
	return (float)dummy;
}

// Compute the sigma of a float array "data" of "ndata" elements. Optionally, the average can be provided if already known.
float sigma(float *data, int ndata, float avg)
{
	double dummy = 0.0;
	if(avg == 0.0)
		avg = average(data,ndata);
	for(int i=0; i<ndata; i++)
		dummy += pow(data[i]-avg,2);
	dummy /= ndata;
	return (float)sqrt(dummy);
}

// Compute the variance of a float array "data" of "size" elements. Optionally, the "mean" can be provided if already known.
float variance(float *data, int size, float mean)
{
	float square = 0.0;

	// Compute sum of squares
	for(int i = 0; i < size; i++)
		square += powf(data[i], 2);
	square /= size;

	// Compute average (if necessary)
	if(mean == 0.0)
		mean = average(data, size);

	return square - mean * mean;
}

// Linear regression: y = b + m·x  (correlation = r)
// C-code  --> https://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
//             http://www.code-in-c.com/linear-regression-fitting-a-line-to-data/
// m and b --> http://www.statisticshowto.com/probability-and-statistics/regression-analysis/find-a-linear-regression-equation/
// r       --> https://mathbits.com/MathBits/TISection/Statistics2/correlation.htm
void linear_regression(float *x, float *y, int n, float *b, float *m, float *r)
{
	//    float x_mean = average(x, n);
	//    float y_mean = average(y, n);
	//    float x_variance = variance(x, n, x_mean);
	//    float y_variance = variance(y, n, y_mean);
	//
	//    float xy_mean = 0.0;
	//    for(int i = 0; i < n; i++)
	//    	xy_mean += x[i] * y[i];
	//    xy_mean /= n;
	//
	//    *m = (xy_mean - (x_mean * y_mean) ) / x_variance;
	//
	//    *b = y_mean - (*m * x_mean);

	float sumx=0.0;
	float sumx2=0.0;
	float sumy=0.0;
	float sumy2=0.0;
	float sumxy=0.0;
	for (int i=0;i<n;i++)
	{
		sumx  += x[i];
		sumx2 += powf(x[i],2);
		sumxy += x[i] * y[i];
		sumy  += y[i];
		sumy2 += powf(y[i],2);
	}

	*m = ((sumxy/n) - (sumx*sumy)/(n*n)) / (sumx2/n - sumx*sumx/(n*n));
	*b = sumy/n - (*m * sumx/n);

	// compute correlation coeff
	*r = (sumxy - sumx * sumy / n) / sqrtf((sumx2 - powf(sumx,2)/n) * (sumy2 - powf(sumy,2)/n));
	//    *r = n*(xy_mean - x_mean * y_mean) / sqrtf( x_variance * y_variance );
}
