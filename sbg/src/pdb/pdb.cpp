/***************************************************************************
                          main.cpp  -  description
                             -------------------
    begin                : Wed Sep 22 18:08:48 CEST 2004
    copyright            : (C) 2004 by Jose Ignacio Garzon
    email                :
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdlib.h>
#include <libpdb/include/Macromolecule.h>
#include <libpdb/include/pdbIter.h>
#include <libtools/include/timer.h>
#include <cmath>
#include <libpdb/include/SS_DISICL.h>
#include <string.h>


#include <cmdl/CmdLine.h>
using namespace TCLAP;
using namespace std;

char input_pdb[256],output_pdb[256], base_name[256], out_file[256];
float euler[3],offset[3];
int euler_convention;
bool only_back,only_ca,in_origin,in_rnumb,in_hy, in_wat, in_name=false, in_dihe, in_nma=false;
char in_chain, in_res[30]={}, in_lsms;
bool in_ncac;

int main(int argc, char *argv[])
{
	//Read parameters
	void parseOptions(int argc, char** argv);


	timer t;
	Tcoor mov;
	Tcoor original_center, invert_original_center;
	char *ext_in,   *ext_out;
	init_aminoacids(Rosseta,pdb);
	int type_in=0; int type_out=0;

	parseOptions(argc,argv);

	//Read Macromolecule
	Macromolecule *mol=new Macromolecule("pdb");




	ext_in=strrchr(argv[1],'.');
	if(ext_in!=NULL)
	{
		if((pdb_strcasecmp(ext_in,".pdb")==0)||(pdb_strcasecmp(ext_in,".ent")==0))
		{
			type_in=0;
		}
		if(pdb_strcasecmp(ext_in,".mol2")==0)
		{
			type_in=1;
		}
		if(pdb_strcasecmp(ext_in,".sdf")==0)
		{
			type_in=2;
		}
	}

	ext_out=strrchr(argv[2],'.');
	if(ext_out!=NULL)
	{
		if((pdb_strcasecmp(ext_out,".pdb")==0)||(pdb_strcasecmp(ext_out,".ent")==0))
		{
			type_out=0;
		}
		if(pdb_strcasecmp(ext_out,".mol2")==0)
		{
			type_out=1;
		}
		if(pdb_strcasecmp(ext_out,".sdf")==0)
		{
			type_out=2;
		}

	}


	switch(type_in)
	{
	case 0:
		fprintf(stdout,"emt_pdb>\nemt_pdb> Reading molecule in PDB format\n");
		break;
	case 1:
		fprintf(stdout,"emt_pdb>\nemt_pdb> Reading molecule in mol2 format\n");
		break;
	case 2:
		fprintf(stdout,"emt_pdb>\nemt_pdb> Reading molecule in sdf format\n");
		break;
	default:
		fprintf(stdout,"emt_pdb>\nemt_pdb> Warning input format %s not supported\n",ext_in);
	}

	if(!mol->readPDB(input_pdb))
	{
		printf("emt_pdb> Error: Cannot open file %s\n",input_pdb);
		return 0;
	}

	//Delete hydrogens if selected
	if(in_hy)
	{
		fprintf(stdout,"emt_pdb> Removing hydrogens\n");
		mol->deleteHYDS();

	}

	//Delete waters if selected
	if(in_wat)
	{
		fprintf(stdout,"emt_pdb> Removing waters\n");
		mol->delete_waters();
	}

	mol->info(stdout);


	if (in_nma) {
		pdbIter *iter_seg;
		pdbIter *iter_ch;
		pdbIter *iter_prot = new pdbIter( mol, false, true, true, true );

		Segment *seg;
		Molecule *pro;
		Chain *ch;


		for( iter_prot->pos_molecule = 0; !iter_prot->gend_molecule(); iter_prot->next_molecule() )
		{
			pro = (Molecule *) iter_prot->get_molecule();

			iter_ch = new pdbIter( pro );

			for( iter_ch->pos_chain = 0; !iter_ch->gend_chain(); iter_ch->next_chain() ) {
				ch = iter_ch->get_chain();
				iter_seg = new pdbIter(ch);
				for( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
				{
					seg = ( Segment * ) iter_seg->get_segment();

					if ((seg->getLimit()==1) )  {

						//printf ("%d %d %d  %d \n",  iter_seg->pos_segment, iter_ch->pos_chain, iter_prot->pos_molecule, ch->getLimit());


						ch->erase(iter_seg->pos_segment);

						// if (iter_seg->pos_segment ==0) iter_seg->pos_segment++;


						// pro->erase(iter_prot->pos_molecule);

						// ch->eraseAll();
						//if (pro->getLimit()==1) {
						// printf ("borar %d %d %d \n", iter_seg->pos_segment, iter_ch->pos_chain, iter_prot->pos_molecule);
						// pro->erase(iter_ch->pos_chain);
						//}
					}
				}
				delete iter_seg;
			}   delete iter_ch;
		}
		mol-> init();
		mol->info(stdout);

	}

	// select chain
	if (in_chain!=NULL)
	{
		Chain *ch;
		pdbIter *iter_ch;
		int sel_Ch=-1;
		iter_ch = new pdbIter( mol ); // iters current protein
		// getting chain number
		for( iter_ch->pos_chain = 0; !iter_ch->gend_chain(); iter_ch->next_chain() ) {
			ch = iter_ch->get_chain();
			if (in_chain==ch->getName()[0])
				sel_Ch=iter_ch->pos_chain ;
		}
		delete iter_ch;
		if (sel_Ch==-1) {
			printf("emt_pdb> WARNING chain %c not found\n", in_chain);
		} else
			printf("emt_pdb> Selecting chain %c\n", in_chain);
		Conditions *conds= new Conditions();
		Condition *cond= new Condition(-1,-1,-1,-1,sel_Ch,sel_Ch,-1,-1,-1,-1);
		conds->add(cond);
		Macromolecule * mol2=mol->select(conds);
		mol=mol2;
		//mol->writePDB("cadena.pdb",!in_rnumb,in_name);
	}

	// selecting residues
	if (in_res[0]!=NULL) {
		int initres,endres;
		sscanf(in_res,"%d-%d",&initres,&endres);
		printf("emt_pdb> Selecting residues from %d to %d\n", initres, endres);
		Conditions *conds= new Conditions();
		Condition *cond= new Condition(-1,-1,-1,-1,-1,-1,initres,endres,-1,-1);
		conds->add(cond);
		Macromolecule * mol2=mol->select(conds);
		mol=mol2;
		mol->writePDB("cadena.pdb",!in_rnumb,in_name);
	}

	if (in_dihe) {
		// dihedral angles....
		float *di;
		int *ss, dangi;
		Residue * res;
		int resn, Nres;
		Segment * seg;
		int chino;
		float * chis;
		int simple;

		chis = new float[4];

		FILE *f_out=NULL;

		strncpy(base_name, argv[1], strlen(argv[1])-4);
		sprintf( out_file, "%s.ss", base_name );


		if( !(f_out = fopen(out_file,"w")) )
		{
			fprintf(stderr,"Sorry, unable to write output file! Forcing exit!\n");
			exit(2);
		}

		pdbIter * iter1;
		iter1 = new pdbIter( mol, false, false, true, false );
		int num_res = ( int )iter1->num_fragment();
		//printf("Nres %d %d\n", num_res, mol->get_num_fragments());

		// mol->format_residues(4);
		mol->all_dihedrals( &di);


		// first SS assign
		dihedrals2DISICL( di, &ss, num_res );


		pdbIter  *iter2;
		dangi=0; Nres=0;
		Chain *ch;
		fprintf(f_out, "# AA   resn  CH     PHI      PSI     OMEGA     DISICL  simple    chi1     chi2     chi3      chi4\n");
		for ( iter1->pos_chain = 0; !iter1->gend_chain(); iter1->next_chain() ) // Iter chains to compute sequential distances fine
		{
			ch = (Chain *) iter1->get_chain();
			iter2 = new pdbIter( ch );
			for ( iter2->pos_fragment = 0; !iter2->gend_fragment(); iter2->next_fragment() )
			{
				res = ( Residue * ) iter2->get_fragment();
				simple=DISICL_d2simple[ss[Nres]];  // change DISICL simple SS classes

				fprintf( f_out, "%3s %8d %2s  %8.3f %8.3f %8.3f  %2d %3s %2d %3s  %8.3f %8.3f %8.3f %8.3f\n",
						res->getName(), res->getIdNumber(), ch->getName(),
						di[dangi], di[dangi+1],  di[dangi+2], // (phi,psi,omega)
						ss[Nres], DISICL_d[ss[Nres]],	simple, DISICL_s[simple],
						di[dangi+3], di[dangi+4], di[dangi+5], di[dangi+6]  // (4x Chi)
				);

				dangi+=7;
				Nres++;
			}
			iter2->~pdbIter();
		}
		iter1->~pdbIter();

		fclose(f_out);
		type_out=3;

	}









	//select backbone
	if ((only_back)&&(type_in==0))
	{
		fprintf(stdout,"emt_pdb> Restricted to backbone\n");
		Conditions *conds= new Conditions();
		Condition *cond= new Condition(-1,-1,-1,-1,-1,-1,-1,-1,-1,-1);
		cond->add(" C  ");
		cond->add(" CA ");
		cond->add(" N  ");
		cond->add(" O  ");
		cond->add(" CB ");
		cond->add(" HA ");
		cond->add(" H  ");
		conds->add(cond);
		Macromolecule * mol2=mol->select(conds);
		mol=mol2;
	}

	//Select Alfa Carbons
	if ((only_ca)&&(type_in==0))
	{
		fprintf(stdout,"emt_pdb> Restricted to CA atoms\n");
		Conditions *conds= new Conditions();
		Condition *cond= new Condition(-1,-1,-1,-1,-1,-1,-1,-1,-1,-1);
		cond->add(" CA ");
		conds->add(cond);
		Macromolecule * mol2=mol->select(conds);
		mol=mol2;
	}

	// select NCAC
	if (in_ncac)
		{

		    if(mol->format_residues(false,0) > 0)
			{
				printf( "pdb> Warning, missing atom(s) receptor found! Forcing exit!\n" );
			}

			fprintf(stdout,"emt_pdb> Restricted to NCAC\n");
			Conditions *conds= new Conditions();
			Condition *cond= new Condition(-1,-1,-1,-1,-1,-1,-1,-1,-1,-1);
			cond->add(" C  ");
			cond->add(" CA ");
			cond->add(" N  ");
			conds->add(cond);
			Macromolecule * mol2=mol->select(conds);
			mol=mol2;
		}





	if ((euler[0]==0)&&(euler[1]==0)&&(euler[2]==0))
		fprintf(stdout,"emt_pdb> No rotation...\n");
	else {


		//Compute center
		mol->geoCenter( original_center );
		invert_original_center[0] = original_center[0] * ( -1.0 );
		invert_original_center[1] = original_center[1] * ( -1.0 );
		invert_original_center[2] = original_center[2] * ( -1.0 );
		fprintf(stdout,"emt_pdb> Center at %f %f %f\n",original_center[0],original_center[1],original_center[2]);

		fprintf(stdout,"emt_pdb> Rotating PDB wait a moment...\n");

		//Movement to the center of coordinates (necessary for rotation)
		mol->moveAll( invert_original_center );
		// printf("%s> original center: %f %f %f\n",argv[0],original_center[0],original_center[1],original_center[2]);

		//Rotation
		euler[0]*=3.1415297/180.0;
		euler[1]*=3.1415297/180.0;
		euler[2]*=3.1415297/180.0;
		printf("emt_pdb> Euler: %f %f %f\n",euler[0],euler[1],euler[2]);

		EulerRot *rot=new EulerRot(euler[0],euler[1],euler[2],euler_convention);
		mol->applyAtoms(rot);

		//If the movement is relative, return the Macromolecule to its initial position
		if(!in_origin)
			mol->moveAll(original_center);
	}

	if ((offset[0]==0)&&(offset[1]==0)&&(offset[2]==0))
		fprintf(stdout,"emt_pdb> No translation...\n");
	else {
		//Apply movement
		fprintf(stdout,"emt_pdb> Moving PDB.wait a moment...\n");

		mol->moveAll(offset);

	}

	// save a surface for VMD
	if(in_lsms>0)
		mol->lsms(in_lsms, false);

	//write macromolecule


	switch(type_out)
	{
	case 0:
		fprintf(stdout,"emt_pdb> Saving molecule in PDB format\nemt_pdb>\n");
		mol->writePDB(argv[2],!in_rnumb,in_name);
		break;
	case 1:
		fprintf(stdout,"emt_pdb> Saving molecule in mol2 format\nemt_pdb>\n");
		mol->writeMol2(argv[2]);
		break;
	case 2:
		fprintf(stdout,"emt_pdb> Saving molecule in sdf format\nemt_pdb>\n");
		mol->writeSDF(argv[2]);
		break;
	case 3:
		fprintf(stdout,"emt_pdb> Saving %s file \nemt_pdb>\n",out_file);
		break;
	default:
		fprintf(stdout,"emt_pdb> Saving molecule in PDB format\nemt_pdb>\n");
		mol->writePDB(argv[2],!in_rnumb,in_name);

	}


}


void parseOptions(int argc, char** argv)
{

	string temp;
	CmdLine cmd(argv[0],"   ", "1.01" );

	try {
		//
		// Define required arguments no labeled
		//
		UnlabeledValueArg<string> pdb_in("pdb_in","pdb input file","default","pdb_in");
		cmd.add( pdb_in );

		UnlabeledValueArg<string> pdb_out("pdb_out","pdb output file","default","pdb_out");
		cmd.add( pdb_out );

		UnlabeledValueArg<float> E1("Euler_1","First Euler Angle",0,"Euler_1");
		cmd.add( E1 );

		UnlabeledValueArg<float> E2("Euler_2","Second Euler Angle",0,"Euler_2");
		cmd.add( E2 );

		UnlabeledValueArg<float> E3("Euler_3","Third Euler Angle",0,"Euler_3");
		cmd.add( E3 );

		UnlabeledValueArg<float> Ox("X_Offset","X offset",0,"X_offset");
		cmd.add( Ox );

		UnlabeledValueArg<float> Oy("Y_Offset","Y offset",0,"Y_offset");
		cmd.add( Oy );

		UnlabeledValueArg<float> Oz("Z_Offset","Z offset",0,"Z_offset");
		cmd.add( Oz );

		SwitchArg recept("r","relative", "rotation from current center of molecule (does not center the molecule)", false);
		cmd.add( recept );
		//
		// Define labeled arguments
		//

		ValueArg<int> Ec("","Ec", "(by default: 0) Convention of the Euler angles (0=x-convention 1=xyz-convention 2=y-convention)",false,0,"int");
		cmd.add( Ec );

		SwitchArg bck("","bck", "delete side chains", true);
		cmd.add( bck );

		SwitchArg ncac("","ncac", "minimal backbone", true);
		cmd.add( ncac );

		SwitchArg ca("c","ca", "only CAs", true);
		cmd.add( ca );

		SwitchArg hy("","hyd", "Erase hydrogens", true);
		cmd.add( hy );

		SwitchArg wat("","wat", "Erase waters", true);
		cmd.add( wat );


		SwitchArg rnumb("","rnumb", "Renumber the residues", true);
		cmd.add( rnumb );

		SwitchArg dihe("","dihedral", "Save dihedral & secondary structure", true);
		cmd.add( dihe );

		SwitchArg name("n","name", "Simplify atoms names", false);
		cmd.add( name );

		SwitchArg nma("","nma", "Remove single aa chains ", false);
		cmd.add( nma );

		ValueArg<char> sel_chain("","selchain", "select chain",false,NULL,"char");
		cmd.add( sel_chain );

		ValueArg<string> range("s","selrange", " Residue range e.g. -r 6-19",false,"","string");
		cmd.add( range );

		ValueArg<float> lsms("","ms", "save a molecular surface",false,NULL," float");
		cmd.add( lsms );



		// Parse the command line.
		cmd.parse(argc,argv);

		if (range.isSet()) {
			strcpy(in_res,((temp=range.getValue()).c_str()));
		} else in_res[0]=NULL;

		// selecting chain
		if (sel_chain.isSet())
			in_chain=sel_chain.getValue();
		else in_chain=NULL;

		if (nma.isSet()) {
			in_nma=true;
		}
		if (lsms.isSet())
			in_lsms=lsms.getValue();
		else in_lsms=0;

		strcpy(input_pdb,((temp=pdb_in.getValue()).c_str()));
		strcpy(output_pdb,((temp=pdb_out.getValue()).c_str()));

		euler[0]=E1.getValue();
		euler[1]=E2.getValue();
		euler[2]=E3.getValue();

		offset[0]=Ox.getValue();
		offset[1]=Oy.getValue();
		offset[2]=Oz.getValue();

		in_origin=!recept.isSet();

		euler_convention=Ec.getValue();
		if (bck.isSet())
			only_back=true;
		else
			only_back=false;

		if (ncac.isSet())
			in_ncac=true;
		else
			in_ncac=false;

		if (dihe.isSet())
			in_dihe=true;
		else
			in_dihe=false;


		if (ca.isSet())
		{
			only_ca=true;
			only_back=false;
		}
		else
			only_ca=false;

		if (rnumb.isSet())
			in_rnumb=true;
		else
			in_rnumb=false;

		if (hy.isSet())
			in_hy=true;
		else
			in_hy=false;

		if (wat.isSet())
			in_wat=true;
		else
			in_wat=false;


		if(name.isSet())
		{
			in_name=true;
		}
		else
		{
			in_name=false;
		}

	} catch ( ArgException& e )
	{
		cout << "  Error->" << e.error() << " " << e.argId() << endl;
	}
	//cmd.~CmdLine();

}
