#include <libenergy/include/korpe.h> // KORP energy
#include <libenergy/include/rama.h> // Rama energy

#include <libpdb/include/pdbIter.h>
#include <libpdb/include/Macromolecule.h>

#include "time.h" // Santy's timer (parallelization independent)
#include <tclap/CmdLine.h>
using namespace TCLAP;


#define FILENAME 300
#define VERSION "v1.18"
char prog[]="korpm";


// AA list    0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19
char aa[] = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'};

// Return aminoacid index given its 1-letter code character.
int aa2index(char aa, char *list)
{
	for(int i=0; i<20; i++)
		if(aa == list[i] || aa-32 == list[i]) // this handles lowercase characters...
			return i;
	fprintf(stderr,"\naa2index> Sorry, %c aminoacid not found!!! Forcing exit!\n\n",aa);
	exit(1);
	return -1;
}


int main(int argc, char *argv[])
{
	bool debug = false;
	int cgmodel = 0; // Atomic model (for formatting)

	std::string temp;
	char *maindir; // Main Directory Path
	char *mainpdb;
	char *pre_pdbid;
	char *pdbid;
	char *input;
	char *mapfile; // Potential energy map filename
	char *ramafile; // Dunbrack's neighbor-dependent PDFs filename for Ramachandran-based potentials
	char *output; // Output file name
	char *weightstr; // Weights for potential energy map (a quoted array of whitespace separated values)
	maindir = (char *) malloc( sizeof(char) * FILENAME );
	mainpdb = (char *) malloc( sizeof(char) * FILENAME );
	pre_pdbid = (char *) malloc( sizeof(char) * FILENAME );
	pdbid = (char *) malloc( sizeof(char) * FILENAME );
	input = (char *) malloc( sizeof(char) * FILENAME );
	mapfile = (char *) malloc( sizeof(char) * FILENAME );
	ramafile = (char *) malloc( sizeof(char) * FILENAME );
	output = (char *) malloc( sizeof(char) * FILENAME );
	weightstr = (char *) malloc( sizeof(char) * FILENAME );
	mainpdb = (char *) malloc( sizeof(char) * FILENAME );

	double *W6D; // Weights array for KORP maps
	bool custom_weights = false; // Custom weights

	korp *map; // KORP's map structure with all stuff required...
	bool rama_switch = false; // Set "true" to activate Rama potential (added as bonding contribution to the selected "emodel" potential)
	float rama_factor = 10.0; // Ramachandran energy factor (typically 10)
	int rama_model = 3; // Ramachandran energy model
	int run_mode = -1;
	bool dump_sums = false; // Dump partial sums into file
	bool weighted = true; // To weighted (true) or not to weighted (false)
	bool single_chain = false; // Use just the mutated chain indicated in the mutation code (e.g. IB163A will only use chain B)
	bool binding = false; // Use ddG = ddGab - ddGa (where ddGab = dGab_Mut - dGab_WT, just valid for single mutation)
	bool pclase = false; // Homology class index
	bool rsa_switch = false; //Read Relative Surface Area (RSA) from input file. RSA value is expected at 4th or 5th (if --class) columns
	bool dexp_switch = false;  // Read DDG exp from file
	bool frev_switch = false;  // Read DDG exp from file

	// COMMAND-LINE PARSER:
	using namespace TCLAP;
	CmdLine cmd(argv[0],' ', VERSION );
	try {

		UnlabeledValueArg<std::string> Input("list","List of PDB/MUT. Format \"160L A A  1 -0.200  AA120M\" ",true, "default","mutations");
		cmd.add( Input );

		ValueArg<std::string> Output("o","output", "Output file name (default=\"out.txt\")",false,"out.txt","string");
		cmd.add( Output );

		ValueArg<std::string> Weights("","weights", "Weights for potential energy map (a quoted array of whitespace separated values), e.g. \"1.545 1.174 ... 2.181\". Typically, 20 values. (default=disabled)",false,"weights","string");
		cmd.add( Weights );

		ValueArg<std::string> ScoreFile("","score_file", "Potential energy map:"
				"\n\t - ICOSA's plain-text format: <AAi><AAj> <face> <dist> <energy> (for energy model 1)"
				//				"\n\t - Dunbrack's neighbor dependent PDFs for our Ramachandran potential (our dunbrack.bin file) (for energy model 6)"
				"\n\t - KORP's energy format (3D, 4D, or 6D, for energy model 5) (default)",false,"score_file","string");
		cmd.add( ScoreFile );

		ValueArg<float> RamaFactor("","rf", "Ramachandran energy factor, typically 10 goes fine. (Default=10)",false,10.0,"float");
		cmd.add( RamaFactor );

		ValueArg<std::string> RamaFile("","rama_file", "Neighbor-dependent PDFs must be provided here for any Ramachandran-based potential (our dunbrack.bin file)",false,"rama_file","string");
		cmd.add( RamaFile );

		ValueArg<std::string> MainDir("d","dir", "Main Directory Path (PDBs directory path), e.g. pdbsall",false,"pdbs","dir_path");
		cmd.add( MainDir );

		ValueArg<int> Runmode ("r","runmode","0  unweighted\n 5 20 Waa \n 21 20Waa+C", false, 5,"int");
		cmd.add( Runmode);

		SwitchArg Binding("","binding", "Use ddG = ddGab - ddGa (where ddGab = dGab_Mut - dGab_WT, just valid for single mutation)", false);
		cmd.add( Binding );

		SwitchArg SingleChain("","single_chain", "Use just the mutated chain indicated in the mutation code (e.g. IB163A will only use chain B)", false);
		cmd.add( SingleChain );

		SwitchArg Dump_sums("","opt", "Enable saving partial sums for optimization (NOT COMPATIBLE WITH -r)", false);
		cmd.add( Dump_sums );

		SwitchArg Pclase("","class", "Read write homology class info", false);
		cmd.add( Pclase );


		SwitchArg Frev("","frev", "Fake reverse mutation", false);
		cmd.add( Frev );

		SwitchArg RSA("","rsa", "Read Relative Surface Area (RSA) from input file. RSA value is expected at 4th or 5th (if --class) columns.", false);
		cmd.add( RSA);

		SwitchArg Dexp("","dexp", "Read  DDG experimental", false);
		cmd.add( Dexp);


		// Parse the command line.
		cmd.parse(argc,argv);

		strcpy(input,((temp=Input.getValue()).c_str())); // Input file name

		strcpy(output,((temp=Output.getValue()).c_str())); // Output file name

		strcpy(maindir,((temp=MainDir.getValue()).c_str())); // Main Directory Path

		strcpy(mapfile,((temp=ScoreFile.getValue()).c_str()));

		rama_factor = RamaFactor.getValue();

		run_mode = Runmode.getValue();;

		if(RamaFile.isSet())
		{
			strcpy(ramafile,((temp=RamaFile.getValue()).c_str()));
			rama_switch = true;
		}


		if(Frev.isSet())
		{
			frev_switch = true;
		}

		if (run_mode != 0)
			weighted = true; // Enable weighted workflow



		if( Weights.isSet() )
		{
			strcpy(weightstr,((temp=Weights.getValue()).c_str())); // Get weights weights for potential energy map (a quoted array of whitespace separated values)
			custom_weights = true; // Custom weights
			weighted = true; // Enable weighted workflow
		}

		single_chain = SingleChain.isSet(); // Use just the mutated chain indicated in the mutation code (e.g. IB163A will only use chain B)

		binding = Binding.isSet();  // Use ddG = ddGab - ddGa (where ddGab = dGab_Mut - dGab_WT, just valid for single mutation)

		if(Dump_sums.isSet())
		{
			dump_sums = Dump_sums.getValue();
			weighted = false;
		}

		if(single_chain and binding)
		{
			fprintf(stdout, "Don't be ridiculous! Single_chain and binding options are incompatible\n");
			exit(1);
		}

		if( binding and (run_mode == 6 || run_mode == 56) )
		{
			fprintf(stdout, "Work In Progress... to be done\n");
			exit(1);
		}

		if(Pclase.isSet()) {
			pclase = true;
			dexp_switch = true;
		}

		if(RSA.isSet()) {
			rsa_switch = true;
			dexp_switch = true;
		}

		if(Dexp.isSet())
			dexp_switch = true;

	}
	catch ( ArgException& e )
	{
		std::cout << "  Error->" << e.error() << " " << e.argId() << std::endl;
	}
	using namespace std;

	// Read weights
	if(custom_weights)
	{
		FILE *MYOUT = stdout;

		fprintf(MYOUT,"Read weights string: %s\n", weightstr);

		// Count whitespaces
		int len = strlen(weightstr); // Get weights string length

		// Get total number of weights
		int nw = 1; // Number of weights
		bool first = true;
		for(int i = 0; i < len; i++)
			if(weightstr[i] == ' ' and first)
			{
				nw++;
				first = false;
			}
			else
				first = true;

		W6D = (double *) malloc( sizeof(double) * nw ); // Allocate the required memory for weights


		// Store the indices of the first character for each weight
		first = true;
		int *iw; // Array with with indices of the first character for each weight
		iw = (int *) malloc( sizeof(int) * (nw + 1) ); // Allocate memory
		iw[0] = 0; // Store position
		int p = 1; // Array position
		for(int i = 0; i < len; i++)
			if(weightstr[i] == ' ' and first)
			{
				iw[p] = i; // Store position
				p++;
				first = false;
			}
			else
				first = true;

		// Load the weights into array
		for(int i = 0; i < nw; i++) // Screen all weigths
		{
			float myw = 0.0;
			sscanf( weightstr + iw[i], "%f", &myw );
			W6D[i] = myw;
		}

		fprintf(MYOUT, "Loaded %d weights: ", nw);
		for(int i = 0; i < nw; i++)
			fprintf(MYOUT, " %.3f", W6D[i]);
		fprintf(MYOUT, "\n");
	}


	// Initialize aminoacids (ADP)
	init_aminoacids(Rosseta,pdb);


	// Wrapper to read a whole KORP energy map (3/4/6D)
	map = readKORP(mapfile);

	fprintf(stdout, "Dimensions %d Cutoff %f Radial %d Cells %d Chis %d Maps %d\n", map->dimensions, map->cutoff,  map->nr, map->ncells, map->nchis, map->nsmaps);



	cgmodel = 0; // 0= N,CA,C model

	char *aa1,*aa2,*chn;
	int *posInChain, *aa2I, *aa1I;
	int nmut = 0;

	// esto no debia estar aqui...
	FILE *f, *fout;
	float deltaE;
	char mutstr[500];

	if ( (fout=fopen(output, "w"))==NULL)
	{
		fprintf(stderr, "\n  Error->Cannot open file %s\n", output);
		exit(1);
	}

	if ( (f=fopen(input, "r"))==NULL)
	{
		fprintf(stderr, "\n  Error->Cannot open file %s\n", input);
		exit(1);
	}

	// Dunbrack stuff
	int size = 72; // Number of bins per dimension of the 2D PDF
	float step = 5; // Angular increment [deg] (may depend on Dunbrack's data set)
	float *****pdfs; // Dunbrack's PDFs
	float *****allmaps; // Complete Ramachandran potential (20x20x20 = 8000 72x72 maps)
	float dihedrals[15];
	float ***ramaps = NULL; // Rama maps for Wild-type protein
	float ***ramaps2 = NULL; // Rama maps for Mutated protein
	Condition *ncac;
	Conditions *ncac2;

	if(run_mode == 6 || run_mode == 56)
	{

		// Reading Dunbrack's PDFs (Rama potential)
		// Load Dunbrack's neighbor-dependent PDFs for Ramachandran-based potentials MAP (Dunbrack_TCB)
		pdfs = read_dunbrack(ramafile,&size,&step); // Reads binary
		fprintf(stdout,"rama2\n" );

		// Rama stuff for protein modeling (many single runs)
		fprintf(stdout,"%s> Generating the complete Ramachandran potential (20x20x20 = 8000 72x72 maps)... ",prog);
		fflush(stdout);
		allmaps = gen_allrama(pdfs,rama_model);
		fprintf(stdout,"Done!\n");
		fflush(stdout);

		if( !(ramaps = (float ***) malloc(sizeof(float **) * 5) ) )
		{
			printf("Sorry, unable to allocate maps memory!!!\n");
			exit(1);
		}
		if( !(ramaps2 = (float ***) malloc(sizeof(float **) * 5) ) )
		{
			printf("Sorry, unable to allocate maps memory!!!\n");
			exit(1);
		}
	}

	float *fwt;
	float *fmut;
	fwt = (float *) malloc( sizeof(float) * 400 );
	fmut = (float *) malloc( sizeof(float) * 400 );

	//	float *fwtA;
	//	float *fmutA;
	//	if(binding)
	//	{
	//		fwtA = (float *) malloc( sizeof(float) * 400 );
	//		fmutA = (float *) malloc( sizeof(float) * 400 );
	//	}

	char myline[1024];
	int nscan = 0; deltaE=0;
	int Hclass; // Homology class index
	float rsa; // Relative Surface Area
	bool diferent_pdb;
	Macromolecule *mol; // Main PDB
	Macromolecule *molr;

	mol = new Macromolecule();
	molr = new Macromolecule();


	while( fgets(myline,1024,f) )
	{


		// Main loop
		//	while(!feof(f))
		//	{
		// Parsing input file
		//    fscanf(f,"%4s %c %c %d %f %s %s %f %f %f",&dumpf, &dumpc, &dumpc, &dumpi, &deltaE, &mutstr, &dumpc1, &dumpf2, &dumpf2, &dumpf2);
		// fscanf(f,"%s %c %c %d %f %s",&dumpf, &dumpc, &dumpc, &dumpi, &deltaE, &mutstr);
		//		fscanf(f,"%s %s %f",&dumpf, &mutstr, &deltaE);
		// fscanf(f,"%s %s %f %*[^\n]\n",&dumpf, &mutstr, &deltaE);

		if ((myline[0] == '#')|| strlen(myline)< 3)
			continue;

		Hclass=-1;
		rsa=-1;

		if (!dexp_switch)  // simple input
			nscan = sscanf(myline,"%s %s ",pdbid, mutstr);
		else if (pclase)
		{
			if(rsa_switch)
				nscan = sscanf(myline,"%s %s %f %d %f",pdbid, mutstr, &deltaE, &Hclass, &rsa);
			else
				nscan = sscanf(myline,"%s %s %f %d",pdbid, mutstr, &deltaE, &Hclass);
		}
		else
		{
			if(rsa_switch)
				nscan = sscanf(myline,"%s %s %f %f",pdbid, mutstr, &deltaE, &rsa);
			else
				nscan = sscanf(myline,"%s %s %f",pdbid, mutstr, &deltaE);
		}

		// fprintf(stdout,"%s %s %f %d %f\n",pdbid, mutstr, deltaE, Hclass, rsa);
		//	getchar();

		sprintf(mainpdb,"%s/%s.pdb",maindir,pdbid);

		diferent_pdb=false;
		if (strcmp(pre_pdbid,pdbid)) {
			// fprintf(stdout, "distinta %8s  %8s\n", pdbid, pre_pdbid);
			diferent_pdb = true;
		}
		strcpy(pre_pdbid,pdbid);


		// multiple mutations
		int lenmut = strlen(mutstr); // Get weights string length
		char list_mutstr[50][20];

		// Get total number of mutations
		nmut = 0; // Number of mutations
		int pstart =0;
		for(int i = 0; i < lenmut; i++)
			if(mutstr[i] == ',') {
				mutstr[i] = ' ';
				sscanf( mutstr + pstart, "%s", (char *) &list_mutstr[nmut]);
				pstart =i;
				//fprintf(stdout, "%s\n", list_mutstr[nmut]);
				nmut++;
			}
		sscanf( mutstr + pstart, "%s", (char *) &list_mutstr[nmut]);
		//fprintf(stdout, "%s\n", list_mutstr[nmut]);

		//fprintf(stdout, "%d\n", nmut);

		aa1 = (char *) malloc( sizeof(char) * nmut );
		aa2 = (char *) malloc( sizeof(char) * nmut );
		chn = (char *) malloc( sizeof(char) * nmut );
		aa1I = (int *) malloc( sizeof(int) * nmut );
		aa2I = (int *) malloc( sizeof(int) * nmut );
		posInChain = (int *) malloc( sizeof(int) * nmut );


		for(int i = 0; i <= nmut; i++) {
			sscanf(list_mutstr[i], "%c%c%d%c", &aa1[i], &chn[i], &posInChain[i], &aa2[i]);
			aa1I[i] = aa2index(aa1[i],aa);
			aa2I[i] = aa2index(aa2[i],aa);
			if (debug) fprintf(stdout, "%8s %7.3f %6s %c %c %4d %c\n", mainpdb, deltaE, mutstr, aa1[i], chn[i], posInChain[i], aa2[i]);
		}




		//fprintf(stdout, "\n\nPDB %s:%c pos %d mut %c(%d) -> %c(%d)\n", mainpdb, chn[0], posInChain[0], aa1[0], aa1I[0], aa2[0], aa2I[0]);

		// Reading main PDB
		if (diferent_pdb)
		{
			delete mol;
			mol = new Macromolecule();

			//printf( "%s> Reading PDB: %s\n", prog, mainpdb);
			mol->readPDB(mainpdb); // reading main PDB


			//	printf( "%s> Deleting Hydrogen atoms (if any)...\n", prog );
			mol->deleteHYDS();
			//	printf( "%s> Deleting Hetero-atoms (if any)...\n", prog );
			mol->delete_heteros();
			//	printf( "%s> Deleting Water molecules (if any)...\n", prog );
			mol->delete_waters();






		}

		// Get chain index (internal) for Reference PDB
		int *ichain;
		ichain = (int *) malloc( sizeof(int) * nmut );
		for(int i = 0; i <= nmut; i++)
			ichain[i] = -1;



		pdbIter *iter_mol = new pdbIter(mol);

		Fragment *frag;
		int r_left=2, r_right=2;

		int firstP=0;
		int lastP=100000;
		for(int i = 0; i <= nmut; i++)
		{

			for( iter_mol->pos_chain = 0; !iter_mol->gend_chain(); iter_mol->next_chain() ) // screen chains
			{
				Chain *ch = iter_mol->get_chain();
				//printf("\n Processing %c %c\n", ch->getName()[0], chn[i]);


				if( ch->getName()[0] == chn[i] )
				{
					ichain[i] = iter_mol->pos_chain; // chain internal index

					pdbIter *iter_seq = new pdbIter( ch ); // iter to screen fragments (residues)


					iter_seq->pos_fragment = 0;
					firstP = ((Residue *) iter_seq->get_fragment())->getIdNumber(); // Get first residue index in PDB numeration

					do  {
						frag =iter_seq->get_fragment();
						iter_seq->next_fragment();
						// printf("\n range %d %d %d\n", firstP, frag->getIdNumber() , posInChain[i]);

					} while(!iter_seq->gend_fragment());

					lastP =frag->getIdNumber();  // Get last residue index in PDB numeration

					r_left=r_right=2;

					if (firstP + 2 > posInChain[i])
						r_left= posInChain[i]-firstP;

					if (lastP - 2 < posInChain[i])
						r_right= lastP - posInChain[i];

					//                        if ((r_left!=2) or (r_right!=2)) {
					//                    	printf("\n range2 %d %d %d %d %d\n", firstP, lastP, posInChain[i], r_left, r_right);
					//  					    getchar();
					//                        }

					//printf("\n range2 %d %d %d %d %d %d %c\n", firstP, lastP, posInChain[i], r_left, r_right,ichain[i], chn[i] );

					// check mutation is in the chain
					if ( (firstP>posInChain[i])||(lastP<posInChain[i]) ) {
						r_right= 2;
						r_left=2;
						continue;
					}

					Condition *mut;
					Conditions *mut2;


					mut= new Condition(-1,-1,-1,-1,ichain[i],ichain[i],posInChain[i],posInChain[i],-1,-1);
					//mut= new Condition(-1,-1,-1,-1,-1,-1,posInChain[i],posInChain[i],-1,-1);
					mut2= new Conditions();
					mut2->add(mut);

					Macromolecule *mutmol2;

					mutmol2 = new Macromolecule();

					mutmol2 = mol->select_cpy(mut2);

					delete  iter_seq;

					if  ( mutmol2->get_num_atoms() == 0 ) {
						delete mutmol2;
						//getchar();

					} else
					{

						delete mutmol2;
						break;
					}


				}

			}



			if(ichain[i] < 0) // Some checking...
			{
				fprintf(stderr,"> Sorry, chain \"%c\" not found in %s (reference PDB)! Using all chains!\n",chn[0],mainpdb);
			}








		}

		delete iter_mol;


		// N,CA,C  Conditions
		ncac = NULL;
		ncac2 = NULL;

		//		Condition *ncacA;
		if(single_chain)
			ncac = new Condition(-1,-1,-1,-1,ichain[0],ichain[0],-1,-1,-1,-1);
		else
		{
			ncac = new Condition(-1,-1,-1,-1,-1,-1,-1,-1,-1,-1);
			//			if(binding)
			//				ncacA = new Condition(-1,-1,-1,-1,ichain,ichain,-1,-1,-1,-1);
		}

		// Select N, CA, C atoms from input PDB

		char *Nstr  = (char *)" N  ";
		char *CAstr = (char *)" CA ";
	  char *Cstr  = (char *)" C  ";


		ncac->add(Nstr);
		ncac->add(CAstr);
		ncac->add(Cstr);
		ncac2 = new Conditions();
		ncac2->add(ncac);

		molr = mol->select_cpy(ncac2);

		// Formatting main PDB
		//	printf( "%s> Formatting main PDB (receptor) residues order and checking for missing atoms at CG-model %d\n", prog, cgmodel );
		//	else	fprintf(stderr, "%s> No missing atom(s) in main PDB (receptor)!\n", prog );

		if(molr->format_residues(false,cgmodel) > 0)
		{
			fprintf( stdout, "%s> Error, missing atom(s) found in main PDB (receptor) according to %d CG-model!\n", prog, cgmodel );
			exit(2);
		}


		if (frev_switch) {
			molr->mutseq(aa2[0],posInChain[0]);
			char dump; int dumpi;
			dump=aa2[0];
			dumpi=aa2I[0];
			aa2[0]=aa1[0];
			aa2I[0]=aa1I[0];
			aa1[0]=dump;
			aa1I[0]=dumpi;

			sprintf(mutstr, "%c%c%d%c", aa1[0], chn[0], posInChain[0], aa2[0]);
			deltaE*=-1.0;
		}




		// Selecting 5 residues for Rama stuff
		int l56size;
		int *iaa = NULL; // AA index for each loop aminoacid
		if (!(run_mode == 6 || run_mode == 56))
				{
			r_right=0;
			r_left=0;
				}


		// Atomic selection Conditions
		Condition *mutzone= new Condition(-1,-1,-1,-1,ichain[0],ichain[0],posInChain[0]-r_left,posInChain[0]+r_right,-1,-1);
		//Condition *mutzone= new Condition(-1,-1,-1,-1,-1,-1,posInChain[0]-r_left,posInChain[0]+r_right,-1,-1);
		Conditions *mutzone2= new Conditions();
		mutzone2->add(mutzone);

		//mol->writePDB("mierda.pdb");
		Macromolecule *mutmol = molr->select_cpy(mutzone2);



		l56size=(r_left+r_right+1);


		// Checking
		if( mutmol->get_num_atoms() != (l56size)*3 )
		{
			fprintf(stdout, "Error, missing N, CA or C atoms %8s %7.3f %6s %c %c %4d %c, expected %d present %d\n", mainpdb, deltaE, mutstr, aa1[0], chn[0], posInChain[0], aa2[0], l56size*3, mutmol->get_num_atoms());
			exit(1);

		}

		if (run_mode == 6 || run_mode == 56) {

		char *mutseq = mutmol->get_sequence(); // PDB sequence in 1-letter code

		// Get the dihedrals for Rama stuff
		float *coord_mut; // 5 residues atomic coordinates
		mutmol->coordMatrix( &coord_mut ); // get coordinates vector
		for(int i = 0; i < 15; i++)
			dihedrals[i]=999;
		finddihedral( coord_mut, l56size*3, dihedrals ); // Get the dihedrals array for supplied coordinates
		free( coord_mut );

		// Setting cis-Prolines in sequence only for protein modeling and Rama potential
		if(setCisPro(dihedrals,mutseq,l56size) > 0)
			fprintf(stdout,"\n%s> Warning, cis-PRO(s) detected: %s\n",prog,mutseq);

		iaa = seq2iaa(mutseq,5); // Translate 1-letter code sequence into our indexing scheme, cis-Prolines included
		if(debug)
		{
			printf( "%s> %d aminoacids sequence from %s is mutseq= %s and iaa=", prog, 5, mainpdb, mutseq );
			for(int i = 0; i < l56size; i++)
				printf( " %d",iaa[i]);
			printf("\n");
			for(int i = 0; i < l56size; i++)
				printf( " %7.2f",dihedrals[i]*180/M_PI);
			printf("\n");
		}
		free( mutseq );
		}
		//			exit(0);
        delete mutmol;


		int ncontM; // Number of contacts

		float *coord = NULL; // Main-PDB coordinates
		int *iseq = NULL;
		contact *contactsM = NULL;
		char *reschainids = NULL; // Main-PDB residue chain-ids
		int *resnums =NULL;      // Main-PDB residue numbers
		char *seq = NULL;         // sequence array
		int nres;

		//		float *coordA = NULL; // Mutated-Chain PDB coordinates
		//		int *iseqA = NULL; // Mutated-Chain sequence indices
		//		contact *contactsMA = NULL;
		//		char *reschainidsA = NULL; // Main-PDB residue chain-ids
		//		int *resnumsA =NULL;      // Main-PDB residue numbers
		//		char *seqA = NULL;         // sequence array
		//		int nresA;


		molr->coordMatrix( &coord );
		reschainids = getResChainIds(molr);
		resnums = getResNums(molr);
		nres = molr->get_num_fragments();
		if (debug)
			printf ("iseq %d nres\n", nres);


		seq = molr->get_sequence();
		iseq = (int *) malloc( sizeof(int) * nres );
		for(int i=0; i<nres; i++)
			iseq[i] = aa2index(seq[i],aa); // numeric index for a given AA

		//sprintf(mainpdb,"%s.pdb",dumpf);
		//molr->writePDB(mainpdb);

		// Compute Environment frames
		frame *frames; // Environment frames
		//		frame *framesA; // Environment frames
		//printf ("frames\n");

		// Getting frames for all...
		frames = frameCoord(coord, nres, resnums, reschainids);
		// MON: make a Mutation version ASAP... to generate a framesM (just the mutated contacts)

		// Compute contacts
		ncontM = contactCoordM(&contactsM, map->cutoff, frames, nres, iseq, posInChain[0], chn[0]);
		//ncontM = contactCoordML(&contactsM, map->cutoff, frames, nres, iseq, posInChain, chn, nmut);

		//printf ("conts %d %d \n", ncont, ncontM);

		//		if(binding)
		//		{
		//			molrA->coordMatrix( &coordA );
		//			reschainidsA = getResChainIds(molrA);
		//			resnumsA = getResNums(molrA);
		//			nresA = molrA->get_num_fragments();
		//			seqA = molrA->get_sequence();
		//			iseqA = (int *) malloc( sizeof(int) * nresA );
		//			for(int i=0; i<nresA; i++)
		//				iseqA[i] = aa2index(seqA[i],aa); // numeric index for a given AA
		//
		//			framesA = frameCoord(coordA, nresA, resnumsA, reschainidsA);
		//			ncontMA = contactCoordM(&contactsMA, map->cutoff, framesA, nresA, iseqA, posInChain, chn);
		//		}

		float delta, Pwt, Pmut, rmse;

		// left/right 1-letter code aminoacids
		int ileft = -1;
		int iright = -1;
		int icenter = -1; // cis- (20) or trans- (12) Proline???
		float EramaWT[3] = { 0.0, 0.0, 0.0 };
		float EramaMut[3] = { 0.0, 0.0, 0.0 };

		// Compute KORP energy
		if(map->dimensions==3) // 3D
		{
			//printf("hola 3D\n");

			switch(run_mode)
			{
			case 0: // REGULAR KORP WITHOUT WEIGHTING
				// Prest = korp3D(contacts,ncont,map);
				Pwt = -korp3DM(contactsM,ncontM,map);
				Pmut = -korp3DM(contactsM, ncontM, map, posInChain[0], chn[0], aa2I[0]);
				delta = Pmut - Pwt;
				rmse = powf( deltaE - delta, 2 );
				break;

			case 1: // DEVELOPMENT SAVING PARTIAL AA CONTIBUTIONS
				fwt = (float *) malloc( sizeof(float) * 400 );
				fmut = (float *) malloc( sizeof(float) * 400 );
				Pwt = -korp3DM(contactsM, ncontM, map, posInChain[0], chn[0], aa1I[0], &(fwt));
				Pmut = -korp3DM(contactsM, ncontM, map, posInChain[0], chn[0], aa2I[0], &(fmut));
				delta = Pmut - Pwt;
				rmse = powf( deltaE - delta, 2 );
				break;

			default: // KORP WITH WEIGHTING
				// double W3D[20]={ 0.755954, 0.787728, 0.761371, 0.778867, 1.049757, 0.827210, 0.610316, 1.055043, 0.773767, 1.175402, 0.822197, 0.834593, 0.740355, 0.619559, 0.986451, 0.832175, 0.847291, 0.990947, 0.874341, 1.069103};
				double W3D[20]={0.85722,0.84165,0.86677,0.91016,1.16726,0.93063,0.61962,1.14223,0.85427,1.23540,0.91959,0.95848,0.82508,0.69785,1.05358,0.93565,0.90401,1.08846,0.99933,1.14419,};
				Pwt = -korp3DMW(contactsM, ncontM, map, posInChain[0], chn[0], aa1I[0], W3D);
				Pmut = -korp3DMW(contactsM, ncontM, map, posInChain[0], chn[0], aa2I[0], W3D);
				delta = Pmut - Pwt;
				rmse = powf( deltaE - delta, 2 );
			}

		}
		else // 6D
		{
			// printf("hola\n");
			switch(run_mode)
			{
			case 56:
			case 5: // KORP WITH WEIGHTING
			{
				if(dump_sums)
				{
					// Initialization
					//					Pwt = -korp6DM(contactsM, ncontM, map, posInChain, chn, aa1I, &(fwt)); // WT
					//					Pmut = -korp6DM(contactsM, ncontM, map, posInChain, chn, aa2I, &(fmut)); // Mut

					if(binding)
					{
						//						Pwt -= -korp6DM(contactsMA, ncontMA, map, posInChain, chn, aa1I, &(fwtA)); // WT (mutated chain)
						//						Pmut -= -korp6DM(contactsMA, ncontMA, map, posInChain, chn, aa2I, &(fmutA)); // Mut (mutated chain)
						//						Pwt -= -korp6DM(contactsMA, ncontMA, map, posInChain, chn, aa1I, fwtA); // WT (mutated chain)
						//						Pmut -= -korp6DM(contactsMA, ncontMA, map, posInChain, chn, aa2I, fmutA); // Mut (mutated chain)
						Pwt = -korp6DM_BIND(contactsM, ncontM, map, posInChain[0], chn[0], aa1I[0], fwt); // WT (mutated chain)
						Pmut = -korp6DM_BIND(contactsM, ncontM, map, posInChain[0], chn[0], aa2I[0], fmut); // Mut (mutated chain)

						if (debug) fprintf(stderr, " ncontM= %d   chn= %s  posInChain= %d\n",  ncontM,  chn, posInChain[0]);
					}
					else
					{

						Pwt = -korp6DM(contactsM, ncontM, map, posInChain[0], chn[0], aa1I[0], fwt); // WT
						Pmut = -korp6DM(contactsM, ncontM, map, posInChain[0], chn[0], aa2I[0], fmut); // Mut

					}
				}
				else
				{
					if(weighted)
					{
						double W6DR[20] = {1.875,0.718,0.876,0.876,2.440,0.973,1.030,2.512,1.030,2.512,1.643,1.113,1.411,1.113,1.608,1.113,1.113,2.512,1.411,2.440};
						double W6DK[20] = {1.875,0.718,0.876,0.876,2.440,0.973,1.030,2.512,1.030,2.512,1.643,1.113,1.411,1.113,1.608,1.113,1.113,2.512,1.411,2.440};

						if(!custom_weights) // Set default weights
						{	if(run_mode == 56) //  KORP + RAMA
								W6D = W6DR;
							else   // KORP
								W6D = W6DK;
						}
						if(binding)
						{
							Pwt = -korp6DMW_BIND(contactsM, ncontM, map, posInChain[0], chn[0], aa1I[0], W6D);
							Pmut = -korp6DMW_BIND(contactsM, ncontM, map, posInChain[0], chn[0], aa2I[0], W6D);
						}
						else
						{

							// Pwt = -korp6DMW(contactsM, ncontM, map, posInChain[0], chn[0], aa1I[0], W6D);
							//Pmut = -korp6DMW(contactsM, ncontM, map, posInChain[0], chn[0], aa2I[0], W6D);

							Pmut=0; Pwt=0;
							for(int mut = 0; mut <= nmut; mut++)
							{
								Pwt += -korp6DMW(contactsM, ncontM, map, posInChain[mut], chn[mut], aa1I[mut], W6D);
								Pmut += -korp6DMW(contactsM, ncontM, map, posInChain[mut], chn[mut], aa2I[mut], W6D); // Mut
							}
							Pwt /= (nmut+1.0);
							Pmut /= (nmut+1.0);
						}
					}
					else
					{
						Pwt = -korp6DM(contactsM, ncontM, map, posInChain[0], chn[0], aa1I[0]);
						Pmut = -korp6DM(contactsM,ncontM, map, posInChain[0], chn[0], aa2I[0]);
					}
				}
				delta = Pmut - Pwt;
				rmse = powf( deltaE - delta, 2 );

				if(run_mode != 56)
					break;

				//			case 56:
				//				fwt = (float *) malloc( sizeof(float) * 400 );
				//				fmut = (float *) malloc( sizeof(float) * 400 );
				//				Pwt = -korp6DM(contactsM, ncontM, map, posInChain, chn, aa1I, &(fwt));
				//				Pmut = -korp6DM(contactsM, ncontM, map, posInChain, chn, aa2I, &(fmut));

			}

			case 6: // Ramachandran alone
			{
				double *Wrama;
				//double WramaK[21] = {1.888,0.699,0.848,0.848,2.447,0.897,1.025,2.567,1.025,2.567,1.698,1.115,1.400,1.115,1.620,1.115,1.115,2.567,1.400,2.447,-0.239};
				double WramaR[21] = {1.888,0.699,0.848,0.848,2.447,0.897,1.025,2.567,1.025,2.567,1.698,1.115,1.400,1.115,1.620,1.115,1.115,2.567,1.400,2.447,-0.239};
				// double WramaK[20] = {1.896,0.721,0.854,0.854,2.436,0.920,1.022,2.545,1.022,2.545,1.667,1.118,1.418,1.118,1.622,1.118,1.118,2.545,1.418,2.436};
				double WramaK[20] = {1.875,0.718,0.876,0.876,2.440,0.973,1.030,2.512,1.030,2.512,1.643,1.113,1.411,1.113,1.608,1.113,1.113,2.512,1.411,2.440};

				if(run_mode == 56) // RAMA weights with KORP
					Wrama = WramaK;
				else // RAMA only weights
					Wrama = WramaR;

				//  Generating neighbor dependent distributions for current loop

				//				if ((r_left!=2) or (r_right!=2)) {
				//					printf("\n range2 %d %d %d %d\n", firstP, lastP, r_left, r_right);
				//					getchar();
				//				}

				int endi=r_right+r_left;
				for(int i=1; i < endi; i++) // Screen relevant AAs (first and last residues from 5-residues region are discarded)
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

					// show_map(allmaps[iaa[i]][ileft][iright],size,"deleteme.txt"); // Gnuplot command:  plot "mierda_01.txt" using 1:2:3 with image
					ramaps[i] = allmaps[iaa[i]][ileft][iright]; // assign the corresponding Rama-energy map to every residue
				}
				//icenter = iaa[2]; // cis/trans Proline?

				iaa[r_left] = aa2I[0]; // Mutate "central" residue of the 5-residues segment

				for(int i=1; i < endi; i++) // Screen relevant AAs (first and last residues from 5-residues region are discarded)
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

					// show_map(allmaps[iaa[i]][ileft][iright],size,"deleteme.txt"); // Gnuplot command:  plot "mierda_01.txt" using 1:2:3 with image
					ramaps2[i] = allmaps[iaa[i]][ileft][iright]; // assign the corresponding Rama-energy map to every residue
				}


				if ((r_left==0)||(r_right==0)) {
					r_left=2;
					endi=0; // pending
				}



				ileft = iaa[r_left-1];
				iright = iaa[r_left+1];

				// N- or C-terminal aminoacids do not have valid Ramachandran (Nt residue lacks Phi and Ct lacks Psi)

				// Ramachandran energy derived from neighbor-dependent PDFs from Dunbrack's paper.
				//				energy = rama_energy2(dihedrals, nres, maps, size);
				if(run_mode == 6) // Only RAMA
				{
					Pwt = 0.0;
					Pmut = 0.0;
				}

				// energy = rama_factor * rama_energy2(dihedrals, nres, ramaps, size);
				int iL,iR;
				for(int i = 0; i < 3; i++) // screen residues
					EramaWT[i]=0;

				for(int i = 1; i < endi; i++) // screen residues
				{
					iL = angle2index(dihedrals[i*3+2], size);
					iR = angle2index(dihedrals[i*3+3], size);
					EramaWT[i-1] = log10(ramaps[ i ][ iL ][ iR ]);
					EramaMut[i-1] = log10(ramaps2[ i ][ iL ][ iR ]);
					//					Pwt += EramaWT[i-1];
					//					Pmut += EramaMut[i-1];
				}

				// Weight correction...
				/*	if(weighted)
				{
					EramaWT[0] *= Wrama[ ileft ];
					EramaWT[1] *= Wrama[ icenter ];
					EramaWT[2] *= Wrama[ iright ];
					EramaMut[0] *= Wrama[ ileft ];
					EramaMut[1] *= Wrama[ iaa[2] ];
					EramaMut[2] *= Wrama[ iright ];
				}
				 */
				for(int i = 0; i < 3; i++) // screen residues
				{
					Pwt += W6D[20]*EramaWT[i];
					Pmut += W6D[20]*EramaMut[i];
				}

				delta = Pmut - Pwt;
				// delta *= 1.0;
				rmse = powf( deltaE - delta, 2 );

				free(iaa); // Not needed anymore
				break;
			}
			case 21:
			{
				// double W6D21[21] = {1.888,0.699,0.848,0.848,2.447,0.897,1.025,2.567,1.025,2.567,1.698,1.115,1.400,1.115,1.620,1.115,1.115,2.567,1.400,2.447,-0.239};
				// double W6D21[21] = {1.896,0.721,0.854,0.854,2.436,0.920,1.022,2.545,1.022,2.545,1.667,1.118,1.418,1.118,1.622,1.118,1.118,2.545,1.418,2.436,-0.239};
				double W6D21[21] =  {1.875,0.718,0.876,0.876,2.440,0.973,1.030,2.512,1.030,2.512,1.643,1.113,1.411,1.113,1.608,1.113,1.113,2.512,1.411,2.440, 0.1};

				if(!custom_weights) {
					W6D = W6D21;
				}

				//Pwt = -korp6DMW21(contactsM, ncontM, map, posInChain[0], chn[0], aa1I[0], W6D);
				//Pmut = -korp6DMW21(contactsM, ncontM, map, posInChain[0], chn[0], aa2I[0], W6D);

				Pwt = -korp6DMW(contactsM, ncontM, map, posInChain[0], chn[0], aa1I[0], W6D);
				Pmut = -(korp6DMW(contactsM, ncontM, map, posInChain[0], chn[0], aa2I[0], W6D)+W6D[20]);


				delta = Pmut - Pwt;
				rmse = powf( deltaE - delta, 2 );

				break;
			}
			case 40:
			{
				double W6DRSA[40] = {1.552,0.515,1.185,1.493,2.007,1.342,0.945,2.528,1.211,3.034,1.354,1.217,1.272,1.027,1.691,1.301,1.377,2.059,1.412,2.262,1.525,1.241,0.706,1.238,2.216,1.023,1.266,2.315,1.997,2.831,1.716,2.474,1.613,0.536,2.198,1.069,1.705,2.264,1.308,2.691};

				if(!custom_weights)
					W6D = W6DRSA;

				Pwt = -korp6DMWRSA(contactsM, ncontM, map, posInChain[0], chn[0], aa1I[0], W6D, rsa);
				Pmut = -korp6DMWRSA(contactsM, ncontM, map, posInChain[0], chn[0], aa2I[0], W6D, rsa);

				delta = Pmut - Pwt;
				rmse = powf( deltaE - delta, 2 );

				break;
			}
			case 99:
			{



				Pwt = W6D[aa1I[0]];
				Pmut = W6D[aa2I[0]];

				delta = Pmut - Pwt;

				rmse = powf( deltaE - delta, 2 );
				//fprintf(stdout, "%f %f %f\n", deltaE, sumaa,rmse);

				break;
			}

			}

			//			for(int c=0; c<20; c++)  {
			//				delta = korp6DM(contactsM,ncontM,map, posInChain, chn, c);
			//				printf("Score %c %f %f\n", aa[c], enat, delta);
			//				}

		}
		// fprintf(fout, "%s %c %d %c %c %d %d %f %f %s\n", mainpdb, chn, posInChain, aa1, aa2, aa1I, aa2I, deltaE, delta, mutstr);

		if (dump_sums)
		{
			fprintf(fout, "%s %-7s %c %c %2d %2d %2d %2d %2d %8.3f %12.6f %10.3f %10.3f %12.6f %12.6f %3d ", pdbid, mutstr, aa1[0], aa2[0], aa1I[0], aa2I[0], ileft, icenter, iright, deltaE, delta, rmse, sqrt(rmse), Pwt, Pmut, ncontM);
			fprintf(stdout, "%s %-7s %c %c %2d %2d %2d %2d %2d %8.3f %12.6f %10.3f %10.3f %12.6f %12.6f %3d ", pdbid, mutstr, aa1[0], aa2[0], aa1I[0], aa2I[0], ileft, icenter, iright, deltaE, delta, rmse, sqrt(rmse), Pwt, Pmut, ncontM);

		}
		else
		{
			if (dexp_switch) {
				fprintf(fout, "%s %-7s %c %c %2d %2d %2d %2d %2d %8.3f %12.6f %10.3f %10.3f %12.6f %12.6f %3d %3d %8.3f", pdbid, mutstr, aa1[0], aa2[0], aa1I[0], aa2I[0], ileft, icenter, iright, deltaE, delta, rmse, sqrt(rmse), Pwt, Pmut, ncontM, Hclass, rsa);
				fprintf(stdout, "%s %-7s %c %c %2d %2d %2d %2d %2d %8.3f %12.6f %10.3f %10.3f %12.6f %12.6f %3d %3d %8.3f", pdbid, mutstr, aa1[0], aa2[0], aa1I[0], aa2I[0], ileft, icenter, iright, deltaE, delta, rmse, sqrt(rmse), Pwt, Pmut, ncontM, Hclass, rsa);
			}
			else {
				fprintf(fout, "%s %-7s %10.3f", pdbid, mutstr, delta);
				fprintf(stdout, "%s %-7s %10.3f", pdbid, mutstr, delta);

			}
		}
		//		if( run_mode == 1 || run_mode == 56)
		if( dump_sums )
		{
			if( run_mode == 6 ) // Mega chapa
			{
				// Write partial sums into output file
				for(int j=0;j<40;j++)
					fprintf(fout," 0 0");
			}
			else
			{
				// Write partial sums into output file
				//				if(binding)
				//				{
				//					for(int j=0;j<20;j++)
				//						if (j!=aa1I)
				//							fprintf(fout," %10f %10f", -fwt[aa1I+20*j] +fwtA[aa1I+20*j], -fwt[j+20*aa1I] +fwtA[j+20*aa1I]);
				//						else
				//							fprintf(fout," %10f %10f", -fwt[aa1I+20*j] +fwtA[aa1I+20*j], 0.0);
				//
				//					for(int j=0;j<20;j++)
				//						if (j!=aa2I)
				//							fprintf(fout," %10f %10f", -fmut[aa2I+20*j] +fmutA[aa2I+20*j], -fmut[j+20*aa2I] +fmutA[j+20*aa2I]);
				//						else
				//							fprintf(fout," %10f %10f", -fmut[aa2I+20*j] +fmutA[aa2I+20*j], 0.0);
				//				}
				//				else
				//				{
				if (pclase)
				{
					fprintf(fout, "%3d ", Hclass);
					fprintf(stdout, "%3d ", Hclass);
				}

				if(rsa_switch)
				{
					fprintf(fout, "%5.1f ", rsa);
					fprintf(stdout, "%5.1f ", rsa);
				}


				for(int j=0;j<20;j++)
					if (j!=aa1I[0])
						fprintf(fout," %10f %10f", -fwt[aa1I[0]+20*j], -fwt[j+20*aa1I[0]]);
					else
						fprintf(fout," %10f %10f", -fwt[aa1I[0]+20*j], 0.0);

				for(int j=0;j<20;j++)
					if (j!=aa2I[0])
						fprintf(fout," %10f %10f", -fmut[aa2I[0]+20*j], -fmut[j+20*aa2I[0]]);
					else
						fprintf(fout," %10f %10f", -fmut[aa2I[0]+20*j], 0.0);
				//				}
			}

			// Ramachandran stuff partial sums
			for(int i=0; i<3; i++)
			{
				fprintf(fout, " %10f %10f", EramaWT[i], EramaMut[i] );
				// fprintf(stderr, " %10f %10f", EramaWT[i], EramaMut[i] );
			}
		}
		fprintf(fout, "\n");
		fprintf(stdout, "\n");

		free(coord);
		free(frames);
		free(contactsM);
		free(iseq);
		free(resnums);
		free(reschainids);
		delete molr;



		//		if(binding)
		//		{
		//			free(coordA);
		//			free(framesA);
		//			free(contactsMA);
		//			free(iseqA);
		//			free(resnumsA);
		//			free(reschainidsA);
		//			delete molrA;
		//		}
	}

	fprintf(stdout, "\n");
	free(fwt);
	free(fmut);

	//	if(binding)
	//	{
	//		free(fwtA);
	//		free(fmutA);
	//	}

	free(ramaps);
	//	free(dihedrals); // free dihedrals

	fclose(f);
	fclose(fout);
}
// END
