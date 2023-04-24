
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <list>

#include <libpdb/include/Macromolecule.h>
#include <libpdb/include/pdbIter.h>
#include <libtools/include/timer.h>
#include <libtools/include/Io.h>
#include <libpdb/include/SS_DISICL.h>
#include <tclap/CmdLine.h>

#define VERSION "v1.1"


using namespace TCLAP;

// Compute the RMSD between two arrays
float atomic_rmsd(float *array, float *array2, int natoms);

// Computes the Gaussian weights profile from inter-atomic distances between "mol" and "mol2"
// allocating profile memory if (*p_profile==NULL).
// Weights are computed according to ec.(7) from: w= exp(-(d^2)/c ) --> c=5A (for big motions)
// Damm & Carlson. "Gaussian-Weighted RMSD Superposition of Proteins: A Structural
// Comparison for Flexible Proteins and Predicted Protein Structures".
// Biophysical Journal. Volume 90, Issue 12, 15 June 2006, Pages 4558-4573
//void gaussian_weight(Macromolecule *mol, Macromolecule *mol2, double **p_profile, double c=5.0);

int main(int argc, char *argv[])
{


	CmdLine cmd(argv[0],' ', VERSION );

	char file_ref[500]; // pdb input reference
	char file_mob[500]; // pdb input
	char file_mob2[500]; // pdb input
	char file_out[500]; // pdb output
	char file_aln[500]; // Sequence alignment file name
	char file_sim[500]; // Similarity file name
	char tmp_path[500]; // Temporary path
	char dist_profiles[] = "dist_profiles.txt";
	float matrix4[4][4];
	std::string temp;
	int aln_level = 2; // Clustal alignment level  (1, 2 or 3 to match *, *: or *:. ,respectively)
	bool *maskres,*maskres2;
	bool *maskmol=NULL, *maskmol2=NULL;
	int nmatch,nali1,nali2;
	double *prof=NULL; // Distance profile before alignment
	double *prof2=NULL; // Distance profile after alignment
	Chain *ch; // Mon(7/1/2016): See below...
	bool sim_binary = true; // Set true to dump a binary similarity matrix (non-redundant), otherwise plain text format used.

	// Reset transformation matrix
	for(int i=0; i < 4; i++)
		for(int j=0; j < 4; j++)
			matrix4[i][j]=0;

	try {
		//
		// Define required arguments no labeled
		//
		UnlabeledValueArg<std::string> pdb_ref("ref_pdb","Reference PDB file",true, "default","pdb");
		cmd.add( pdb_ref );

		UnlabeledValueArg<std::string> pdb_mob("mob_pdb","Mobile PDB file (if Multi-PDB, Model index and initial and aligned RMSDs will be dumped.)",true, "default","pdb");
		cmd.add( pdb_mob );

		// optional arguments

		SwitchArg Dist("","dist", "Compute the inter-residue distance profiles (averaged over residue atoms) before and after alignment: dist_profiles.txt", false);
		cmd.add( Dist );

		ValueArg<int> AlnLevel("","aln_level", "Clustal alignment level (1, 2 or 3 to match *, *: or *:. ,respectively) (default=2)",false,2,"int");
		cmd.add( AlnLevel );

		ValueArg<std::string> AlnFile("a","aln_file", "Clustal's \".aln\" file defining alignment between initial and target atomic structures.\n",false,"alnstring","string");
		cmd.add( AlnFile );

		ValueArg<double> gauss_c("C","gauss_c", "Gaussian width (default=20)",false,20.0,"double");
		cmd.add( gauss_c );

		SwitchArg gauss("g","gauss", "Gaussian-weighted RMSD superposition", false);
		cmd.add( gauss );

		SwitchArg RemoveH("", "noH","Removing Hydrogen atoms from output files and computations.", false);
		cmd.add( RemoveH );

		SwitchArg BBone("b","bb", "only BackBone atoms (N,CA,C,O) in proteins", false);
		cmd.add( BBone );

		SwitchArg CAlpha("c","calpha", "only C alpha in proteins or P in nucleic acids", false);
		cmd.add( CAlpha );

		SwitchArg CB_add("","CB", "add CB to RMSD", false);
		cmd.add( CB_add );

		SwitchArg O_add("","O", "add O to RMSD", false);
		cmd.add( O_add );

		SwitchArg sim_text("","sim_text", "Force plain text format for similarity matrix output.", false);
		cmd.add( sim_text );

		ValueArg<std::string> sim_matrix("s","similarity","Compute the similarity matrix (RMSD). Only for Multi-PDB input.",false,"","string");
		cmd.add(sim_matrix);

		ValueArg<std::string> pdb_moved("o","output","output pdb file (mobile onto reference)",false,"","string");
		cmd.add(pdb_moved);

		ValueArg<std::string> pdb_mob2("","input2","extra PDB to move with transformation",false,"","string");
		cmd.add(pdb_mob2);

		ValueArg<int> Start("","start","N-terminal loop residue index, first mobile residue (PDB residue index). (default=disabled)",false,-1,"int");
		cmd.add( Start );

		ValueArg<int> End("","end","C-terminal loop residue index, last mobile residue (PDB residue index). (default=disabled)",false,-1,"int");
		cmd.add( End );

		ValueArg<char> Chain("","chain", "Chain ID (default=first chain of reference PDB)",false,1,"int");
		cmd.add( Chain );

		ValueArg<int> NFlanks("","flanks", "Number of loop flanking residues. E.g. \"2\" means two flanking residues at Nt and other two at Ct, \"0\" means disabled (default=0).",false,0,"int");
		cmd.add( NFlanks );

		ValueArg<int> ResOffset("","res_offset", "Residue index offset from both ends. Typically used to discard anchors (fixed) in Multi-PDB loops file. (default=0)",false,0,"int");
		cmd.add( ResOffset );

		SwitchArg save_rmsd("r","save_rmsd", " Save rmsd & dihedrals",false);
		cmd.add( save_rmsd );

		// Parse the command line.
		cmd.parse(argc,argv);

		init_aminoacids(Rosseta,pdb);


		// reading pdb
		Macromolecule *mol=new Macromolecule("pdb");
		Macromolecule *mol2=new Macromolecule("pdb2");
		Macromolecule *mol3=new Macromolecule("pdb3");
		Macromolecule *mol_ca;
		Macromolecule *mol2_ca;

		// stupid transformation from string to char
		strcpy(file_ref,((temp=pdb_ref.getValue()).c_str()));
		strcpy(file_mob,((temp=pdb_mob.getValue()).c_str()));
		mol->readPDB(file_ref);
		mol2->readPDB(file_mob);

		if(pdb_mob2.isSet()) {
			strcpy(file_mob2,((temp=pdb_mob2.getValue()).c_str()));
			mol3->readPDB(file_mob2);
		}

		// Removing Hydrogens
		if(RemoveH.isSet())
		{
			mol->delete_hydrogens(); // Remove Hydrogens
			mol2->delete_hydrogens(); // Remove Hydrogens

			//			printf( ">> Hydrogens were removed.\n" );
		}

		pdbIter *iter_molec;
		iter_molec = new pdbIter(mol2);
		int nmol = iter_molec->num_molecule();
		delete(iter_molec);

		int ri = Start.getValue(); // All residues by default
		int rf = End.getValue(); // All residues by default
		int ichain = -1; // All chains by default
		char chain; // Chain ID

		Condition *firstmolec; // First molecule conditions
		firstmolec = new Condition(0,0,-1,-1,-1,-1,-1,-1,-1,-1); // First loop selection
		Conditions *firstmolec2= new Conditions();
		firstmolec2->add(firstmolec);
		Macromolecule *molec0;
		molec0 = mol2->select_cpy(firstmolec2);
		iter_molec = new pdbIter(molec0);

		// Get loop limits and chain id from first model in Multi-PDB (Mobile)
		if(nmol > 1 && (ri < 0 || rf < 0))
		{

			iter_molec->pos_chain=0;
			ch = iter_molec->get_chain();
			pdbIter *iter_chain = new pdbIter(ch);
			if(ri < 0)
			{
				iter_chain->pos_fragment=0;
				ri = (iter_chain->get_fragment())->getIdNumber(); // get Nt anchor residue number (PDB)
			}
			if(rf < 0)
			{
				iter_chain->pos_fragment=iter_chain->num_fragment()-1;
				rf = (iter_chain->get_fragment())->getIdNumber(); // get Ct anchor residue number (PDB), i.e. the last that moves...
			}
		}

		// Chain *ch; // Mon(7/1/2016): Compilers (intel and gcc) say this is not declared within scope... WTF???
		iter_molec->pos_chain = 0;
		ch = iter_molec->get_chain();
		chain = ch->getName()[0];
		delete iter_molec;

		//		// Apply the selected number of residues offset (0 by default)
		//		ri += ResOffset.getValue();
		//		rf -= ResOffset.getValue();
		//		fprintf(stderr, "> Considered chain \"%c\" goes from residue %d (Nt) to %d (Ct), including an index offset of %d.\n", chain, ri, rf, ResOffset.getValue() );

		// Number of loop flanking residues (e.g. "2" means two flanking residues at Nt, and other two at Ct)
		int nflanks = NFlanks.getValue();

		if(nflanks <= 0) // No-flanks stuff...
		{
			// Apply the selected number of residues offset (0 by default)
			ri += ResOffset.getValue();
			rf -= ResOffset.getValue();
			fprintf(stderr, "> Considered chain \"%c\" goes from residue %d (Nt) to %d (Ct), including an index offset of %d.\n", chain, ri, rf, ResOffset.getValue() );
		}

		// Get chain index (internal) for Reference PDB
		pdbIter *iter_mol = new pdbIter(mol);
		for( iter_mol->pos_chain = 0; !iter_mol->gend_chain(); iter_mol->next_chain() ) // screen chains
		{
			ch=iter_mol->get_chain();
			if(ch->getName()[0] == chain)
			{
				ichain = iter_mol->pos_chain; // chain internal index
				break;
			}
		}
		delete iter_mol;
		if(ichain < 0) // Some checking...
		{
			fprintf(stderr,"> Sorry, chain \"%c\" not found in %s (reference PDB)! Using all chains!\n",chain,file_ref);
			//exit(1);
		}
		// Note: chain index "0" is assumed for Mobile PDB

		// Atomic selection Conditions
		Condition *calpha= new Condition(-1,-1,-1,-1,ichain,ichain,ri,rf,-1,-1);
		// Condition *calpha= new Condition(-1,-1,ichain,ichain,-1,-1,ri,rf,-1,-1); // Mon (9/4/2018): Bug, we must select by chain, not segment.
		// MON: Check above lines !!! Rare....

		Condition *calphaX= new Condition(-1,-1,-1,-1,-1,-1,ri,rf,-1,-1); // Selection for mobile Multi-PDB
		//		Condition *calphaX= new Condition(-1,-1,ichain,ichain,-1,-1,-1,-1,-1,-1); // Selection for mobile Multi-PDB

		//		          mobile = new Condition(-1,-1,-1,-1,-1,-1,ri+1,rf-1,-1,-1); // get residues from Nt+1 to Ct-1, i.e. only those mobile...

		if( CAlpha.isSet() )
		{
			calpha->add(" CA ");
			calpha->add(" P  ");
		}

		if( CB_add.isSet() )
		{
			calpha->add(" CB ");
		}

		if( O_add.isSet() )
		{
			calpha->add(" O  ");
		}


		if( BBone.isSet() )
		{
			calpha->add(" N  ");
			calpha->add(" CA ");
			calpha->add(" C  ");
			//	calpha->add(" O  ");
		}
		Conditions *calpha2= new Conditions();
		calpha2->add(calpha);

		if( CAlpha.isSet() )
		{
			calphaX->add(" CA ");
			calphaX->add(" P  ");
		}
		if( BBone.isSet() )
		{
			calphaX->add(" N  ");
			calphaX->add(" CA ");
			calphaX->add(" C  ");
			calphaX->add(" O  ");
		}
		Conditions *calpha2X= new Conditions();
		calpha2X->add(calphaX);


		// Compute RMSDs between some loop in the reference PDB and all the loops in the mobile Multi-PDB
		//		if(nmol > 1 || ri >= 0 || rf >= 0) // If mobile Multi-PDB found
		if(nmol > 1) // If mobile Multi-PDB found
		{
			printf( "# %d models found in input mobile Multi-PDB: %s\n", nmol, file_mob);


			// Select relevant atoms
			mol_ca = mol->select(calpha2);
			mol2_ca = mol2->select(calpha2X);


			// Copy trimmed reference Macromolecule (just mobile loop region) to load decoys coordinates later
			Macromolecule *mol2_decoy = new Macromolecule(mol_ca);

			//			}
			//			else
			//			{
			//				fprintf(stderr,"%d missing atoms in mol\n",mol->format_residues(false,2));
			//				fprintf(stderr,"%d missing atoms in mol2\n",mol2->format_residues(false,2));
			//				mol_ca = mol;
			//				mol2_ca = mol2;
			//			}

			int natom = mol_ca->get_num_atoms(); // Get number of atoms in reference loop
			int natom2 = mol2_ca->get_num_atoms(); // Get number of atoms in mobile loops

			// fprintf(stderr,"natom= %d  natom2= %d\n",natom,natom2);

			if(natom2 % nmol != 0) // Some checking
			{
				fprintf(stderr,"Error, different number of atoms between reference PDB and mobile Multi-PDB (natom2= %d, nmol= %d). Forcing exit!\n", natom2, nmol);
				exit(2);
			}

			fprintf(stderr,"natom= %d  natom2= %d  ichain= %d\n",natom,natom2,ichain);

			float *coord;
			float *coord2;
			mol_ca->coordMatrix( &coord );
			mol2_ca->coordMatrix( &coord2 );

			for(int i = 0; i < nmol; i++) // Screen all molecules
			{
				// Set current loop coordinates into dummy loop Macromolecule
				mol2_decoy->coordMatrixSet( coord2 + i * natom * 3 );

				// Compute the RMSD between two arrays
				fprintf(stdout,"%6d %f %f\n",i+1, mol_ca->rmsd(mol2_decoy), mol_ca->minRmsd(mol2_decoy,matrix4) );

				// OLD suff...
				// fprintf(stdout,"%6d %f\n",i+1, atomic_rmsd(coord, coord2 + i * natom * 3, natom) );
			}

			// Computing similarity matrix based on RMSD
			if(sim_matrix.isSet())
			{
				// Get similarity file name from parser
				strcpy( file_sim,( (temp = sim_matrix.getValue()).c_str() ) );

				FILE *f_sim; // Output similarity file handler
				fprintf(stdout,"rmsd> Writing similarity CC-RMSD into file %s\n",file_sim);
				float simil; // Similarity
				int cont = 0;

				if( !sim_text.isSet() ) // binary
				{
					f_sim = fopen(file_sim,"wb");
					fwrite(&nmol,sizeof(int),1,f_sim); // Store the total number of loops

					// Computing similarity CC-RMSD
					for(int i=0; i<nmol; i++)
						for(int j=i+1; j<nmol; j++)
						{
							simil = atomic_rmsd(coord2 + i * natom * 3, coord2 + j * natom * 3, natom);
							fwrite(&simil,sizeof(float),1,f_sim); // Dump similarity into file
							cont++;
						}
				}
				else // plain-text
				{
					f_sim = fopen(file_sim,"w");
					fprintf(f_sim,"# <i> <j> <similarity> (indices begin from 1)\n"); // Matlab request (i begins in value 1)
					// Computing similarity CC-RMSD
					for(int i=0; i<nmol; i++)
						for(int j=i+1; j<nmol; j++)
						{
							simil = atomic_rmsd(coord2 + i * natom * 3, coord2 + j * natom * 3, natom);
							fprintf(f_sim,"%5d %5d %10f\n",i+1,j+1,simil); // Matlab request (i begins in value 1)
							cont++;
						}
				}
				fclose(f_sim);
				fprintf(stdout,"rmsd> %d RMSD values written. (cont= %d)\n", nmol*(nmol-1)/2, cont);
			}

			delete mol2_decoy;
			exit(0);
		}

		// Single loop alignment (in no flanks stuff requested)
		if( nflanks <= 0 && (ri >= 0 || rf >= 0) )
		{
			// Select relevant atoms
			mol_ca = mol->select(calpha2);
			mol2_ca = mol2->select(calpha2X);

			printf("rmsd> Initial_RMSD: %8f", mol_ca->rmsd(mol2_ca));
			printf("  Aligned_RMSD: %8f\n", mol_ca->minRmsd(mol2_ca,matrix4));

			if(pdb_moved.isSet())
			{
				M4Rot *matrix4_op=new M4Rot(matrix4);
				mol2_ca->applyAtoms(matrix4_op); // superpose full-atom
				delete matrix4_op;

				strcpy(file_out,((temp=pdb_moved.getValue()).c_str()));
				mol2_ca->writePDB(file_out);
				printf(">> Saved in %s\n", file_out);
			}

			exit(0);
		}

		// Flanks rigid body alignment using Kabsch & Sander stuff...
		if(nflanks > 0) // Flanks
		{
			// Select relevant atoms (the whole chain)
			mol_ca = mol->select(calpha2);
			mol2_ca = mol2->select(calpha2X);

			// Condition(s) to select flanks
			Conditions *cflank= new Conditions();
			Condition *cflankN = new Condition(-1,-1,-1,-1,ichain,ichain,ri-nflanks,ri-1,-1,-1);
			Condition *cflankC = new Condition(-1,-1,-1,-1,ichain,ichain,rf+1,rf+nflanks,-1,-1);
			cflankN->add(" N  ");
			cflankN->add(" CA ");
			cflankN->add(" C  ");
			cflankC->add(" N  ");
			cflankC->add(" CA ");
			cflankC->add(" C  ");
			cflank->add(cflankN);
			cflank->add(cflankC);

			// Select relevant atoms (the whole chain)
			Macromolecule *mol_flanks = mol->select(cflank);
			mol_flanks->writePDB("mol_flanks.pdb");

			Conditions *cflank2= new Conditions();
			Condition *cflankN2 = new Condition(-1,-1,-1,-1,ichain,ichain,ri-nflanks,ri-1,-1,-1);
			Condition *cflankC2 = new Condition(-1,-1,-1,-1,ichain,ichain,rf+1,rf+nflanks,-1,-1);
			cflankN2->add(" N  ");
			cflankN2->add(" CA ");
			cflankN2->add(" C  ");
			cflankC2->add(" N  ");
			cflankC2->add(" CA ");
			cflankC2->add(" C  ");
			cflank2->add(cflankN2);
			cflank2->add(cflankC2);
			Macromolecule *mol2_flanks = mol2->select(cflank2);
			mol2_flanks->writePDB("mol2_flanks.pdb");

			// Compute transformation of target PDB and RMSD evaluation
			float matrix4[4][4];
			float flank_rmsd = mol_flanks->rmsd(mol2_flanks); // Compute un-aligned RMSD
			float flank_rmsd_min = mol_flanks->minRmsd(mol2_flanks, matrix4); // Compute aligned RMSD

			printf("rmsd> Flanks RMSD: Initial= %8f  Aligned= %8f (%d residues per flank)\n", flank_rmsd, flank_rmsd_min, nflanks);

			if(pdb_moved.isSet())
			{
				// Alignment of target PDB (minRMSD) (computed with flanks selection, but applied to full model)
				M4Rot *matrix4_op = new M4Rot(matrix4);
				mol2->applyAtoms(matrix4_op); // superpose
				delete matrix4_op;

				strcpy(file_out,((temp=pdb_moved.getValue()).c_str()));
				mol2->writePDB(file_out);
				printf("rmsd> Complete mobile PDB was aligned to flanks into %s\n", file_out);
			}
			mol2->writePDB("aliflanks.pdb");

			delete mol_ca;
			delete mol2_ca;
			delete mol_flanks;
			delete mol2_flanks;
			exit(0);
		}


		if(AlnFile.isSet())
		{
			strcpy(file_aln,((temp=AlnFile.getValue()).c_str())); // Gets Alignment file name
			aln_level = AlnLevel.getValue();

			printf( "%s> Reading Clustal's sequence alignment from: %s\n","",file_aln);
			read_aln(file_aln, &maskres, &nali1, &maskres2, &nali2, &nmatch, aln_level);
			printf( "%s> %d residues matched!  nali1= %d  nali2= %d\n","",nmatch,nali1,nali2);
		}

		// only selection is considered
		float rmsd_i,rmsd_f;
		if( CAlpha.isSet() || BBone.isSet() )
		{
			// if CA / BB selection
			mol_ca=mol->select(calpha2);
			mol2_ca=mol2->select(calpha2);
			//			mol->writePDB("mymol.pdb");
			//			mol_ca->writePDB("mymolca.pdb");
			//			mol2->writePDB("mymol2.pdb");
			//			mol2_ca->writePDB("mymol2ca.pdb");

			if(AlnFile.isSet())
			{
				maskmol = mol_ca->maskres2maskatom(maskres);
				maskmol2 = mol2_ca->maskres2maskatom(maskres2);
				if(save_rmsd.isSet())
				{				printf("> Saving initial.rmsd\n");

				rmsd_i = mol_ca->rmsd_file(mol2_ca,maskmol,maskmol2, "initial.rmsd");
				}
				rmsd_i = mol_ca->rmsd(mol2_ca,maskmol,maskmol2);

				// here save dihedral....
			}
			else
				rmsd_i = mol_ca->rmsd(mol2_ca);

			//			fprintf(stderr,"ichain=%d  ri=%d  rf=%d  maskmol %d:\n", ichain, ri, rf, mol_ca->get_num_atoms());
			//			for(int x=0; x<mol_ca->get_num_atoms(); x++)
			//				fprintf(stderr," %1d",maskmol[x]);
			//			fprintf(stderr,"\n");

			if(pdb_moved.isSet() || gauss.isSet())
				if(AlnFile.isSet())
					rmsd_f = mol_ca->minRmsd(mol2_ca,matrix4,maskmol,maskmol2);
				else
					rmsd_f = mol_ca->minRmsd(mol2_ca,matrix4);
			else
				if(AlnFile.isSet())
				{
					// rmsd_f = mol_ca->minRmsd(mol2_ca,maskmol,maskmol2);
					// fprintf(stderr,"Not implemented yet!\n");
					// exit(0);
					rmsd_f = mol_ca->minRmsd(mol2_ca,matrix4,maskmol,maskmol2); // A different function (without matrix4 computation) should be made...
				}
				else
					rmsd_f = mol_ca->minRmsd(mol2_ca);
		}
		else // All atoms
		{
			if(AlnFile.isSet())
			{
				maskmol = mol->maskres2maskatom(maskres);
				maskmol2 = mol2->maskres2maskatom(maskres2);
				if(save_rmsd.isSet())
				{
					printf("> Saving initial.rmsd\n");
					rmsd_i = mol->rmsd_file(mol2,maskmol,maskmol2,"initial.rmsd");
				}
				rmsd_i = mol->rmsd(mol2,maskmol,maskmol2);

				if(save_rmsd.isSet())
				{
					// output dihedrals && SS
					//
					float *di, *di2;
					int *ss, *ss2, dangi, dangi2;
					Fragment * res, *res2;
					int resn, Nres,Nres2, simple, simple2;
					Segment * seg;
					int chino;
					float * chis;
					pdbIter *iter1, *iter2;
					iter1 = new pdbIter( mol,  false, false, true, false); // r maskres2
					iter2 = new pdbIter( mol2,  false, false, true, false ); // t maskres1
					// first SS assign

					mol->all_dihedrals( &di);
					mol2->all_dihedrals( &di2);

					dihedrals2DISICL( di, &ss, (iter1->num_fragment()) );
					dihedrals2DISICL( di2, &ss2, (iter2->num_fragment()) );

					FILE *f_out=NULL;
					fprintf(stdout,">Saving initial.di\n");
					if( !(f_out = fopen("initial.di","w")) )
					{
						fprintf(stderr,"Sorry, unable to write output file! Forcing exit!\n");
						exit(2);
					}
					fprintf(f_out, "# AA     resn    PHI      PSI     OMEGA     DISICL  simple    chi1     chi2     chi3      chi4\n");
					iter2->pos_fragment = 0;
					dangi=0;
					Nres=0;
					Nres2=0;
					for ( iter1->pos_fragment = 0; !iter1->gend_fragment(); iter1->next_fragment() )
					{
						res = ( Fragment * ) iter1->get_fragment();

						if(maskres[iter1->pos_fragment])
						{
							while(!iter2->gend_fragment() && !maskres2[iter2->pos_fragment])
							{// this places the index into the corresponding pas
								iter2->next_fragment();
							}
							Nres=iter1->pos_fragment;
							Nres2=iter2->pos_fragment;
							dangi=Nres*7;
							dangi2=Nres2*7;

							res2 = ( Fragment * ) iter1->get_fragment();
							simple=DISICL_d2simple[ss[Nres]];  // change DISICL simple SS classes
							simple2=DISICL_d2simple[ss2[Nres2]];  // change DISICL simple SS classes

							fprintf( f_out, "%3s %8d  %8.3f %8.3f %8.3f  %3d %3s %3d %3s  %8.3f %8.3f %8.3f %8.3f   %3s %8d  %8.3f %8.3f %8.3f  %3d %3s %3d %3s  %8.3f %8.3f %8.3f %8.3f\n",
									res->getName(), res->getIdNumber(),
									di[dangi], di[dangi+1],  di[dangi+2], // (phi,psi,omega)
									ss[Nres],DISICL_d[ss[Nres]],	simple, DISICL_s[simple],
									di[dangi+3], di[dangi+4], di[dangi+5], di[dangi+6],  // (4x Chi)

									res2->getName(), res2->getIdNumber(),
									di2[dangi], di2[dangi+1],  di2[dangi+2], // (phi,psi,omega)
									ss2[Nres],DISICL_d[ss2[Nres2]],	simple2, DISICL_s[simple2],
									di2[dangi2+3], di2[dangi2+4], di2[dangi2+5], di2[dangi2+6]  // (4x Chi)

							);
							iter2->next_fragment();
						}
					}
					iter1->~pdbIter();
					iter2->~pdbIter();
					fclose(f_out);
				}
			}
			else
				rmsd_i = mol->rmsd(mol2);

			if (pdb_moved.isSet() || gauss.isSet())
				if(AlnFile.isSet())
					rmsd_f = mol->minRmsd(mol2,matrix4,maskmol,maskmol2);
				else
					rmsd_f = mol->minRmsd(mol2,matrix4);
			else
				if(AlnFile.isSet())
				{
					// rmsd_f = mol_ca->minRmsd(mol2_ca,maskmol,maskmol2);
					// fprintf(stderr,"Not implemented yet!\n");
					// exit(0);
					rmsd_f = mol->minRmsd(mol2,matrix4,maskmol,maskmol2);
				}
				else
					rmsd_f = mol->minRmsd(mol2);
		}
		printf("RMSD: init. %.4f ",rmsd_i);

		printf("\n");
		for(int i=0; i < 4; i++)
		{
			for(int j=0; j < 4; j++)
				printf("%6.3f ", matrix4[i][j]);
			printf("\n");
		}

		if(Dist.isSet())
		{
			if( CAlpha.isSet() || BBone.isSet() )
				mol_ca->dist_profile(mol2_ca,&prof,maskmol,maskmol2);
			else
				mol->dist_profile(mol2,&prof,maskmol,maskmol2);
		}

		if(gauss.isSet()) // YES-Gauss
		{
			int imax = 1000; // maximum number of iterations
			int iter=0;
			double *prof_wrmsd=NULL;
			float conv = 1e-6; // convergence criteria
			float delta = 1e6;
			float score_old = 1e6;
			float score;
			double sum_w=0.0;
			int num_atoms,num_res;

			// Initial alignment (minRMSD) (computed with selection, but applied to full model)
			M4Rot *matrix4_op=new M4Rot(matrix4);
			mol2->applyAtoms(matrix4_op); // superpose full-atom
			delete matrix4_op;

			if( CAlpha.isSet() || BBone.isSet() )
			{
				num_atoms = mol_ca->get_num_atoms();
				num_res = mol_ca->get_num_fragments();
			}
			else
			{
				num_atoms = mol->get_num_atoms();
				num_res = mol->get_num_fragments();
			}

			// Iterate until convergence
			for(iter=0; iter<imax && delta>conv; iter++)
			{
				score_old = score;
				if( CAlpha.isSet() || BBone.isSet() )
				{
					//			norm_profile(prof_wrmsd, mol_ca->get_num_atoms(), 0.0, 1.0);
					if(AlnFile.isSet())
					{
						mol_ca->gaussian_weight(mol2_ca,&prof_wrmsd,maskmol,maskmol2,gauss_c.getValue());
						score = mol_ca->minWRmsd(mol2_ca,matrix4,prof_wrmsd,maskmol,maskmol2); // computing transformation
					}
					else
					{
						mol_ca->gaussian_weight(mol2_ca,&prof_wrmsd,gauss_c.getValue());
						score = mol_ca->minWRmsd(mol2_ca,matrix4,prof_wrmsd); // computing transformation
					}
					matrix4_op = new M4Rot(matrix4);
					mol2->applyAtoms(matrix4_op);
				}
				else
				{
					//			norm_profile(prof_wrmsd, mol->get_num_atoms(), 0.0, 1.0);
					if(AlnFile.isSet())
					{
						mol->gaussian_weight(mol2,&prof_wrmsd,maskmol,maskmol2,gauss_c.getValue());
						score = mol->minWRmsd(mol2,matrix4,prof_wrmsd,maskmol,maskmol2); // computing transformation
					}
					else
					{
						mol->gaussian_weight(mol2,&prof_wrmsd,gauss_c.getValue());
						score = mol->minWRmsd(mol2,matrix4,prof_wrmsd); // computing transformation
					}
					matrix4_op = new M4Rot(matrix4);
					mol2->applyAtoms(matrix4_op);
				}
				delta = fabs(score_old - score);
				//			printf("iter= %d  score= %f  old= %f\n",iter,score,score_old);
			}
			if(AlnFile.isSet())
				rmsd_f = mol->rmsd(mol2,maskmol,maskmol2);
			else
				rmsd_f = mol->rmsd(mol2);

			for(int i=0; i < num_atoms; i++)
				sum_w += prof_wrmsd[i];
			//		  printf(">> rmsd_i %.4f -> rmsd_f %.4f  wrmsd= %f %%match= %.2f iters= %d C= %.1f\n",rmsd_i,rmsd_f,score,100*(sum_w/num_atoms),iter,gauss_c.getValue());
			printf(" final %.4f  wrmsd= %f  %%Matching: atoms= %.2f  res= %.2f  iter= %d  C= %.1f\n",rmsd_f,score,100*(sum_w/num_atoms),100*(sum_w/num_res),iter,gauss_c.getValue());
		}
		else
			printf(" final %.4f\n",rmsd_f);

		if(pdb_moved.isSet())
		{
			if(!gauss.isSet()) // NO-Gauss (with Gauss, transformation was already applied)
			{
				M4Rot *matrix4_op=new M4Rot(matrix4);
				mol2->applyAtoms(matrix4_op); // superpose full-atom
				delete matrix4_op;
			}

			if(Dist.isSet())
				//				mol->dist_profile(mol2,&prof2,maskmol,maskmol2);
			{
				if( CAlpha.isSet() || BBone.isSet() )
					mol_ca->dist_profile(mol2_ca,&prof2,maskmol,maskmol2);
				else
					mol->dist_profile(mol2,&prof2,maskmol,maskmol2);
			}

			strcpy(file_out,((temp=pdb_moved.getValue()).c_str()));
			if(pdb_mob2.isSet()) {
				M4Rot *matrix4_op=new M4Rot(matrix4);
				mol3->applyAtoms(matrix4_op); // superpose full-atom
				mol3->writePDB("temp2.pdb");
			}

			mol2->writePDB(file_out);
			printf(">> Saved in %s\n", file_out);
		}

		if(Dist.isSet())
		{
			FILE *f_prof;
			if((f_prof=fopen(dist_profiles,"w"))==NULL)
			{
				fprintf(stderr, "\n Error opening file: %s\n\n",dist_profiles);
				exit(1);
			}

			pdbIter *it = new pdbIter(mol,false); // without SMol
			int num_res = it->num_fragment();
			delete it;
			// printf("Number of fragments: %d\n",num_res);
			if(pdb_moved.isSet())
			{
				fprintf(f_prof,"# Res  Before[A]   After[A]\n");
				for(int i=0; i < num_res; i++)
					fprintf(f_prof,"%5d %10f %10f\n",i+1,prof[i],prof2[i]);
				free(prof2);
			}
			else
			{
				fprintf(f_prof,"# Res  Before[A]\n");
				for(int i=0; i < num_res; i++)
					fprintf(f_prof,"%5d %10f\n",i+1,prof[i]);
			}
			fclose(f_prof);
			printf(">> Distance profile(s) saved in %s\n", dist_profiles);
			free(prof);
		}


		// END
	}
	catch ( ArgException& e )
	{
		std::cout << "  Error->" << e.error() << " " << e.argId() << std::endl;
	}
}

// Compute the RMSD between two arrays
float atomic_rmsd(float *array, float *array2, int natoms)
{
	double accum = 0.0;
	for(int i = 0; i < natoms*3; i++)
	{
		accum += powf(array[i] - array2[i], 2);
	}

	return( (float) sqrt(accum / natoms) );
}
