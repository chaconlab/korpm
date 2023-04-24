#include "Macromolecule.h"

#include <stdio.h>
// #include "ResIni.h"

#define RESIDUE_DIST 25 // Mon: RESIDUE_DIST= d^2  , where "d" is the inter-CA distance between consecutive aminoacids
#define NUCLE_DIST 99999 // Mon: NUCLE_DIST= d^2  , where "d" is the inter-P distance between consecutive nucleotides
enum INPUT_MOL {input_null,input_nucleacid,input_protein,input_smol};



//Auxiliar Function
float sqrt_dist(float *d1, float *d2)
{
	return ( (d1[0]-d2[0])*(d1[0]-d2[0]) + (d1[1]-d2[1])*(d1[1]-d2[1]) + (d1[2]-d2[2])*(d1[2]-d2[2]));
}


bool Macromolecule::readPDB( char *name)
{
	char extension[10];
	char *point;
	point = strrchr(name, '.');
	int type=0;
	int i=0;

	if(point!=NULL)
	{
		while(point[i]!='\0')
		{
			extension[i]=point[i];
			i++;
		}
		extension[i]='\0';

		type=3;
		if(pdb_strcasecmp(extension,".pdb")==0)
		{
			type=0;
		}
		if(pdb_strcasecmp(extension,".ent")==0)
		{
			type=0;
		}
		if(pdb_strcasecmp(extension,".pqr")==0)
		{
			type=0;
		}
		if(pdb_strcasecmp(extension,".mol2")==0)
		{
			type=1;
		}
		if(pdb_strcasecmp(extension,".sdf")==0)
		{
			type=2;
		}
	}

	switch(type)
	{
	case 0:
		return readPDB_pdb(name);
	case 1:
		return readMol2(name);
	case 2:
		return readSDF(name);
	default:
		fprintf(stderr,"Unknown atomic structure datatype (%s). Please use .pdb, .ent, .mol or .sdf extension\n",extension);
		exit(1);
		break;
	}

}

bool Macromolecule::readPDB_pdb( char *name)
{
	FILE *file;

	char line[83];
	char x_c[9],y_c[9],z_c[9],occ_c[7], tfac_c[7];
	char * fmt;
	int atom_no, res_no, i,k, cont_CA=0,cont_res=0;
	double occ, tfac;
	char atom_name[5], res_name[4],old_res_name[4], rec_type[7];
	//	name_len = sizeof( atom_name );
	char number[6];


	char *ext;
	ext=strrchr(name,'.');

	if(ext!=NULL)
	{
		if(strcmp(ext,".pqr")!=0)
			fmt = ( char * ) "%6c%5c%*c%4c%c%3c%*c%c%4d%c%*3c%8c%8c%8c%6c%6c%c%3c";
		else
		{
			fmt = ( char * ) "%6c%5c%*c%4c%c%3c%*c%c%4d%c%*3c%8c%8c%8c%8c%8c%c%3c";
		}
	}
	else
	{
		fmt = ( char * ) "%6c%5c%*c%4c%c%3c%*c%c%4d%c%*3c%8c%8c%8c%8c%8c%c%3c";
	}

	for(i=0;i<9;i++)
	{
		x_c[i]=y_c[i]=z_c[i]='\0';
	}
	for(i=0;i<7;i++)
	{
		occ_c[i]=tfac_c[i]='\0';
	}

	//	for ( i = 0; i < name_len; atom_name[i++] = ' ' );
	//	for ( i = 0; i < 5; atom_name[i++] = ' ' );
	//	atom_name[4] = '\0';

	char blank, alt_loc, entity_id[2], insert_code,old_insert_code, no[4];
	entity_id[1] = '\0';

	for(i=0;i<4;i++)
	{
		no[i] = '\0';
		res_name[i] = '\0';
	}
	for(i=0;i<7;i++)
	{
		rec_type[i] = '\0';
	}

	number[5] = '\0';

	Molecule * mol=NULL;
	Chain * ch=NULL;
	Segment * seg=NULL;
	Fragment * frag=NULL;
	Atom * at=NULL;

	double x, y, z;
	float pos[3];
	char sym[3];
	sym[2] = '\0';
	char old_entity_id[2];
	int old_res_no;
	float old_pos[3],old_ca_pos[3];
	INPUT_MOL mol_type,old_mol_type;
	float dist_ca=-1.0;
	bool CA_previous=false;
	bool res_found;
	bool jump=true;
	bool inserted_frag=false,newsegment;
	bool hexadecimal=false;
	bool first_model=true;
	bool first_warning=true;

	old_mol_type=input_null;
	strcpy(old_res_name,"");
	old_insert_code='Z';
	strcpy( old_entity_id, "" );
	old_res_no = -1;
	old_pos[0]=old_pos[1]=old_pos[2]=-1000000;
	old_ca_pos[0]=old_ca_pos[1]=old_ca_pos[2]=-1000000;

	//Open the file
	file=fopen(name,"rt");
	if(file==NULL) {
		fprintf(stderr,"    Can't open -%s- for reading: No such pdb file\n\n",name);
		exit(1);
	}


	//Reads file line
	cont_res=0;
	while ( fgets( line, 82, file ) != 0 )
	{
		occ = 1.0;
		jump=true;
		//Interprets the read line
		if ( sscanf( line, fmt, rec_type, number, atom_name, & alt_loc, res_name, & entity_id[0], & res_no, & insert_code, x_c,
				y_c, z_c, occ_c, tfac_c, & blank, & no[0] ) == EOF )
		{
			fclose(file);
			return true;
		}
		x=atof(x_c);
		y=atof(y_c);
		z=atof(z_c);
		occ=atof(occ_c);
		tfac=atof(tfac_c);
		atom_name[4] = '\0';

		atom_no=atoi(number);
		// fprintf(stderr,"atom_name=%s--- atom_no= %d (number=%s)\n",atom_name,atom_no,number);
		//getchar();
		if( strcmp( rec_type, "ATOM  " ) == 0 || strcmp( rec_type, "HETATM" ) == 0)
		{
			jump=false;
			if(hexadecimal)

			{
				atom_no=0;
			}
			for(i=0;i<5;i++)
			{
				if((number[i]>=65 && number[i]<=70) || (number[i]>=97 && number[i]<=102))
				{
					//fprintf(stderr,"Hola: %c %s\n",number[i],number);
					//sscanf(number,"%x",&atom_no);
					atom_no=0;
					hexadecimal=true;
					i=6;
				}
				if(number[i]=='\0')
					i=5;
			}
		}


		//printf("%s %s %d\n",atom_name,res_name, CA_previous);

		//if(alt_loc!=' ' && alt_loc!='A' )
		if(alt_loc>65 && alt_loc<=90)
		{
			if( strcmp( rec_type, "ATOM  " ) == 0)
			{
				if (first_warning)
				{
					first_warning=false;
					// Mon: Hasta los cojones de ver estos warnings... con uno va que chuta!
					fprintf(stderr,"WARNING Alternate locations found in %s. E.g: Alt-loc[%d %s] (only \"A\" location considered)",name,atom_no,atom_name);
					//					 fprintf(stderr,"WARNING reading %s : Alt-loc[%d %s]",name,atom_no,atom_name);
				}
				//				 else
				//					 fprintf(stderr," Alt-loc[%d %s]",atom_no,atom_name);

			}

			jump=true;
		}

		if ( strcmp( rec_type, "     " ) == 0
				||  strcmp( rec_type, "" ) == 0
				||  strcmp( rec_type, "END   " ) == 0 )
		{
			jump=true;
		}
		//if ( strcmp( rec_type, "TER   " ) == 0 )
		if(rec_type[0]=='T' && rec_type[1]=='E' && rec_type[2]=='R')
		{

			old_mol_type=input_null;
			strcpy(old_res_name,"");
			old_insert_code='Z';
			strcpy( old_entity_id, "" );
			old_res_no = -1;
			old_pos[0]=old_pos[1]=old_pos[2]=-1000000;
			old_ca_pos[0]=old_ca_pos[1]=old_ca_pos[2]=-1000000;
			dist_ca=0;
			CA_previous=false;
			cont_CA=0;
			jump=true;
		}

		if ( strcmp( rec_type, "MODEL " ) == 0 )
		{
			if(first_model)
			{
				fprintf(stderr,"WARNING: ReadPDB. File %s has multiple models. Models are read as different Molecule sets\n",name);
				first_model=false;
			}
			old_mol_type=input_null;
			strcpy(old_res_name,"");
			old_insert_code='Z';
			strcpy( old_entity_id, "" );
			old_res_no = -1;
			old_pos[0]=old_pos[1]=old_pos[2]=-1000000;
			old_ca_pos[0]=old_ca_pos[1]=old_ca_pos[2]=-1000000;
			dist_ca=0;
			CA_previous=false;
			cont_CA=0;
			jump=true;
		}
		// fprintf(stderr,"arriba [%s] %d [%s]\n",res_name, mol_type, rec_type);
		// getchar();
		//IF VALID ATOM LINE
		if(!jump)
		{
			//DETERMINE ATOM TYPE

			if ((strcmp(res_name,old_res_name)==0))   mol_type=old_mol_type;

			else if ( strcmp( rec_type, "ATOM  " ) == 0 )
			{
				if(     (strcmp(res_name,"  A")==0) ||
						(strcmp(res_name,"  C")==0) ||
						(strcmp(res_name,"  G")==0) ||
						(strcmp(res_name,"  T")==0) ||
						(strcmp(res_name,"  U")==0) ||
						(strcmp(res_name," DA")==0) ||
						(strcmp(res_name," DG")==0) ||
						(strcmp(res_name," DT")==0) ||
						(strcmp(res_name," DC")==0) ||
						(strcmp(res_name," DU")==0)
				)
				{
					mol_type=input_nucleacid;
				}
				else
				{
					int ii;
					res_found=false;
					for (ii = 0; ii < 30; ii++ ) {
						// fprintf(stderr," [%s] [%s] %d\n",res_name, AA[ii].aa_name3, ii);
						if ( strcmp( res_name, AA[ii].aa_name3 ) == 0 )  { res_found=true; break; }
					}

					if (res_found) mol_type=input_protein;
					else { mol_type=input_smol;
					strcpy(rec_type,"HETATM");
					}

					//fprintf(stderr," [%s] %d %s %d\n",res_name, mol_type, rec_type, ii);
					//getchar();
				}
			}
			else
				if ( strcmp( rec_type, "HETATM" ) == 0 )
				{
					mol_type=input_smol;
				}
				else
				{
					mol_type=input_null;
				}


			//MODIFICATIONS IN READ VALUES
			//**************rename hydrogens
			if(strcmp(atom_name," HN ")==0) // ????
			{
				strcpy(atom_name, " H  "); // ????
			}

			//fprintf(stderr,"fuera [%s] %d %s\n",res_name, mol_type, rec_type);

			k=0;
			while (!isalpha(atom_name[k]) &&  (k<4) )
				k++; // Index "k" pointing to first letter

			if (k==4) {
				fprintf(stderr," Non atom name in line: \n %s \n\n",line);
				exit(1);
			}


			sym[0] = atom_name[k]; // atom identity for 1 Char (C,N,...)
			k++;
			sym[1] =' ';

			if  (isalpha(atom_name[k]) && (mol_type==input_smol) ) {

				if (islower(atom_name[k]))
					atom_name[k]=toupper(atom_name[k]);

				if ((sym[0] =='C')&&(atom_name[k] =='A')&&(strcmp( res_name, "CA ") == 0 )) sym[1] ='A';
				if ((sym[0] =='M')&&(atom_name[k] =='G')) sym[1] ='G';
				if ((sym[0] =='M')&&(atom_name[k] =='N')) sym[1] ='N';
				if ((sym[0] =='N')&&(atom_name[k] =='A')&&(strcmp( res_name, "NA ") == 0 )) sym[1] ='A';
				if ((sym[0] =='Z')&&(atom_name[k] =='N')) sym[1] ='N';
				if ((sym[0] =='F')&&(atom_name[k] =='E')) sym[1] ='E';
				if ((sym[0] =='N')&&(atom_name[k] =='I')) sym[1] ='I';
				if ((sym[0] =='C')&&(atom_name[k] =='U')&&(strcmp( res_name, "CU ") == 0 )) sym[1] ='U';
				if ((sym[0] =='A')&&(atom_name[k] =='L')) sym[1] ='L';
				if ((sym[0] =='C')&&(atom_name[k] =='L')&&(strcmp( res_name, "CL ") == 0 )) sym[1] ='L';
				if ((sym[0] =='C')&&(atom_name[k] =='O')&&(strcmp( res_name, "CO ") == 0 )) sym[1] ='O';
				if ((sym[0] =='B')&&(atom_name[k] =='R')) sym[1] ='R';
				if ((sym[0] =='A')&&(atom_name[k] =='U')) sym[1] ='U';
				if ((sym[0] =='P')&&(atom_name[k] =='T')&&(strcmp( res_name, "PT ") == 0 ) ) sym[1] ='T';
				if ((sym[0] =='H')&&(atom_name[k] =='G')&&(strcmp( res_name, "HG ") == 0 ) ) sym[1] ='G';


				// fprintf(stderr," [%c%c]->[%s] resname [%s]\n",sym[0],sym[1], atom_name, res_name);
			}

			// fprintf(stderr," [%c%c] %f\n",sym[0],sym[1], k);


			//*****************Pos
			pos[0] = ( float )x;
			pos[1] = ( float )y;
			pos[2] = ( float )z;

			//ATOM CREATION
			if(strcmp( entity_id, old_entity_id ) != 0  ||  mol_type!=old_mol_type)
			{
				cont_res=0;
			}

			// at = new Atom( Table_Elements::getElement( sym ), pos, 0.0, atom_name, atom_no, occ, tfac );
			// fprintf(stderr,"About to getElement(sym). sym= %s\n",sym);
			Element *elem = Table_Elements::getElement(sym);
			// fprintf(stderr,"elem.sym= %s\n",elem->sym);
			if (elem == NULL)
			{
				fprintf(stderr,"Warning! Unknown atom: [%c%c]->[%s] resname [%s], changed to C (carbon).\n",sym[0],sym[1], atom_name, res_name);
				at = new Atom( Table_Elements::getElement((char *)"C "), pos, 0.0, (char *)" C  ", atom_no, occ, tfac );
			}
			else
				at = new Atom( elem, pos, 0.0, atom_name, atom_no, occ, tfac );

			//CASE PROTEIN
			if(mol_type==input_protein)
			{

				//*****************CA
				if(strcmp( atom_name, " CA " ) == 0) // ????
				{
					cont_CA++;
					if(!CA_previous)
						dist_ca=0.0;
					else
						dist_ca=sqrt_dist(pos,old_ca_pos);
					old_ca_pos[0]=pos[0];
					old_ca_pos[1]=pos[1];
					old_ca_pos[2]=pos[2];
					CA_previous=true;
				}
				else
					dist_ca=0;

				//No New residue?
				if ( old_res_no == res_no
						&& insert_code==old_insert_code
						&& strcmp(res_name,old_res_name)==0
						&&  (cont_CA<2)
						&& strcmp( entity_id, old_entity_id ) == 0
						&& mol_type==old_mol_type)
				{
					frag->add(at);
				}
				else //New Residue
				{

					frag=new Residue(res_name,res_no,cont_res, insert_code );
					cont_res++;
					if(cont_res==9999)
						cont_res=0;

					inserted_frag=false;
					frag->add(at);
					if(strcmp( at->getName(), " CA " ) == 0) // ????
						cont_CA=1;
					else
						cont_CA=0;

				}

				if( (dist_ca>=RESIDUE_DIST) || !inserted_frag)
				{
					//No New Segment?
					if(res_no==(old_res_no+1)
							&& dist_ca<RESIDUE_DIST
							&& strcmp( entity_id, old_entity_id ) == 0
							&& mol_type==old_mol_type)
					{
						seg->add(frag);
						inserted_frag=true;
					}
					else //New Segment
					{
						newsegment=false;
						if( inserted_frag )
						{
							if(seg->getLimit()>1)
							{
								seg->remove(seg->getLimit()-1);
								seg=new Segment((char *)"Segment" );
								seg->add(frag);
								inserted_frag=true;
								newsegment=true;
							}
						}
						else
						{
							seg=new Segment((char *)"Segment" );
							seg->add(frag);
							inserted_frag=true;
							newsegment=true;
						}

						if(newsegment)
						{
							//No New Chain?
							if(strcmp( entity_id, old_entity_id ) == 0
									&& mol_type==old_mol_type)
							{
								ch->add(seg);
							}
							else //New Chain
							{

								ch=new Chain( entity_id );

								ch->add(seg);
								CA_previous=false;

								//No new Molecule?
								if(mol_type==old_mol_type)
								{
									mol->add(ch);
								}
								else //New Protein
								{
									mol=new Protein( ( char * ) "Protein" );
									mol->add(ch);
									this->add(mol);
								}
							}
						}
					}
				}
			}
			//END CASE PROTEIN


			//CASE NUCLEACID
			if(mol_type==input_nucleacid)
			{
				//*****************CA
				if(strcmp( atom_name, " P  " ) == 0)
				{
					cont_CA++;
					if(!CA_previous) {
						dist_ca=0.0;
					}
					else {
						dist_ca=sqrt_dist(pos,old_ca_pos);
					}
					old_ca_pos[0]=pos[0];
					old_ca_pos[1]=pos[1];
					old_ca_pos[2]=pos[2];
					CA_previous=true;
				}
				else
					dist_ca=0;


				//No New residue?
				if ( old_res_no == res_no
						&& insert_code==old_insert_code
						&& strcmp(res_name,old_res_name)==0 &&  (cont_CA<2)
						&& strcmp( entity_id, old_entity_id ) == 0
						&& mol_type==old_mol_type)
				{
					frag->add(at);
				}
				else //New Nucleotide
				{
					frag=new Nucleotide(res_name,res_no,cont_res, insert_code );
					inserted_frag=false;
					frag->add(at);
					if(strcmp( at->getName(), " P  " ) == 0) // ????
						cont_CA=1;
					else
						cont_CA=0;

					cont_res++;
				}

				if( dist_ca>=NUCLE_DIST || !inserted_frag)
				{
					//No New Segment?
					if(res_no==(old_res_no+1)
							&& strcmp( entity_id, old_entity_id ) == 0
							&& dist_ca<NUCLE_DIST
							&& mol_type==old_mol_type)
					{
						seg->add(frag);
						inserted_frag=true;
					}
					else //New Fragment
					{
						newsegment=false;
						if( inserted_frag )
						{
							if(seg->getLimit()>1)
							{
								seg->remove(seg->getLimit()-1);
								seg=new Segment((char *)"Segment" );
								seg->add(frag);
								inserted_frag=true;
								newsegment=true;
							}
						}
						else
						{
							seg=new Segment((char *)"Segment" );
							seg->add(frag);
							inserted_frag=true;
							newsegment=true;
						}


						if(newsegment)
						{
							//No New Chain?
							if(strcmp( entity_id, old_entity_id ) == 0
									&& mol_type==old_mol_type)
							{
								ch->add(seg);
							}
							else //New Chain
							{
								ch=new Chain( entity_id );
								ch->add(seg);
								CA_previous=false;

								//No new Molecule?
								if(mol_type==old_mol_type)
								{
									mol->add(ch);
								}
								else //New NAcid
								{
									mol=new NAcid( ( char * ) "NAcid" );
									mol->add(ch);
									this->add(mol);
								}
							}
						}
					}
				}
			}
			//END CASE NUCLEACID

			//CASE SMOL
			if(mol_type==input_smol)
			{
				//No New residue?
				if ( old_res_no == res_no
						&& insert_code==old_insert_code
						&& strcmp(res_name,old_res_name)==0
						&& strcmp( entity_id, old_entity_id ) == 0
						&& mol_type==old_mol_type)
				{
					mol->add(at);
				}
				else //New SMOL
				{
					mol=new SMol(res_name,res_no,cont_res );
					mol->add(at);
					this->add(mol);
					CA_previous=false;
					cont_res++;
				}
			}
			//END CASE SMOL


			if ( strcmp( rec_type, "ATOM  " ) == 0
					|| strcmp( rec_type, "HETATM" ) == 0 )
			{
				old_mol_type=mol_type;
				strcpy(old_res_name,res_name);
				old_insert_code=insert_code;
				strcpy( old_entity_id, entity_id );
				old_res_no = res_no;
				old_pos[0]=pos[0];
				old_pos[1]=pos[1];
				old_pos[2]=pos[2];
			}

		}//JUMP
	}//WHILE

	if (!first_warning)  fprintf(stderr,"\n");
	fclose(file);
	return true;
}


bool Macromolecule::writePDB( char *name, bool number, bool change_name)
{
	FILE *file;
	bool contP, contCh, contF, contR, contA;
	char line[83];
	char * fmt;
	char * fmt2;
	int res_no = 0;
	int cont_at=1;

	Segment * seg;
	Chain * cad;
	Molecule * p;
	Fragment * res;
	Atom * a;

	Tcoor coor;

	char name_simple[5];

	if(change_name)
	{
		strcpy(name_simple,"    ");
	}

	char *ext;
	ext=strrchr(name,'.');


	if(strcmp(ext,".pqr")!=0)
	{
		fmt=( char * ) "ATOM  %5d %-5s%3s %c%4d%c   %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n";
		fmt2= ( char * ) "HETATM%5d %-5s%3s  %4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n";
	}
	else
	{
		fmt=( char * ) "ATOM  %5d %-5s%3s %c%4d%c   %8.3lf%8.3lf%8.3lf%8.4lf%8.4lf\n";
		fmt2= ( char * ) "HETATM%5d %-5s%3s  %4d    %8.3lf%8.3lf%8.3lf%8.4lf%8.4lf\n";
	}


	if(this->get_num_atoms()<=0)
	{
		return false;
	}


	file=fopen(name,"wt");
	if(file==NULL)
		return false;
	//rename_residues();


	initAll();
	contP = true;
	while ( contP != false )
	{
		p = ( Molecule * ) this->getCurrent();

		if(p->getClass()==pdb_smol)
		{
			contA = true;
			while ( contA != false )
			{
				a = ( Atom * ) p->getCurrent();

				if (change_name)
				{
					name_simple[1]=(a->getElement())->sym[0];
					name_simple[2]=(a->getElement())->sym[1];
					a->setPdbName(name_simple);
					//fprintf(stderr,"hola: %s\n",a->getPdbName());
				}
				/* Escritura de un atomo */
				a->getPosition(coor);

				if(number)
					sprintf( line, fmt2, a->getPdbSerial(), a->getPdbName(), p->getName(),
							p->getIdNumber(), coor[0], coor[1], coor[2], a->getPdbocc(), a->getPdbfact() );
				else
					sprintf( line, fmt2, cont_at++, a->getPdbName(), p->getName(),
							((SMol*)p)->get_pos()+1, coor[0], coor[1], coor[2], a->getPdbocc(), a->getPdbfact() );

				fputs( line, file );
				if(cont_at==100000)
					cont_at=1;
				contA = p->next();
			}
		}
		else
		{
			contCh = true;
			while ( contCh != false )
			{
				cad = ( Chain * ) p->getCurrent();
				cont_at=1;

				contF = true;
				while ( contF != false )
				{
					seg = ( Segment * ) cad->getCurrent();
					contR = true;
					while ( contR != false )
					{
						res = ( Fragment * ) seg->getCurrent();
						res_no++;

						contA = true;
						while ( contA != false )
						{
							a = ( Atom * ) res->getCurrent();

							if (change_name)
							{
								name_simple[1]=(a->getElement())->sym[0];
								name_simple[2]=(a->getElement())->sym[1];
								a->setPdbName(name_simple);
							}
							/* Escritura de un atomo */
							a->getPosition(coor);
							if(number) {
								sprintf( line, fmt, a->getPdbSerial(), a->getPdbName(), res->getName(),
										cad->getName() [0],res->getIdNumber()/*res->get_pos()*/,res->get_letter(), coor[0], coor[1], coor[2], a->getPdbocc(), a->getPdbfact() );
							}
							else {
								sprintf( line, fmt, cont_at++, a->getPdbName(), res->getName(),
										cad->getName() [0],res->get_pos()+1,res->get_letter(), coor[0], coor[1], coor[2], a->getPdbocc(), a->getPdbfact() );
							}
							if(cont_at==100000) {
								cont_at=1;
							}

							fputs( line, file );
							contA = res->next();
						}
						contR = seg->next();
					}
					contF = cad->next();
				}
				contCh = p->next();
			}

			if(number) {
				sprintf( line, "TER    %4d      %s %c%4d                                                      \n",
						a->getPdbSerial() + 1, res->getName(),
						cad->getName() [0], res->getIdNumber());
			}
			else {
				sprintf( line, "TER    %4d      %s %c%4d                                                      \n",
						a->getPdbSerial() + 1, res->getName(),
						cad->getName() [0], res->get_pos()+1);
			}
			fputs( line, file );
		}
		contP = this->next();
	}

	fclose(file);
	return true;
}

bool Macromolecule::writeMPDB( char *name, int n, bool number)
{
	FILE *file;
	bool contP, contCh, contF, contR, contA;
	char line[83];
	char * fmt;
	char * fmt2;
	int res_no = 0;
	int cont_at=1;
	Segment * seg;
	Chain * cad;
	Molecule * p;
	Fragment * res;
	Atom * a;
	Tcoor coor;

	char *ext;
	ext=strrchr(name,'.');


	if(strcmp(ext,".pqr")!=0)
	{
		fmt=( char * ) "ATOM  %5d %-5s%3s %c%4d%c   %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n";
		fmt2= ( char * ) "HETATM%5d %-5s%3s  %4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n";
	}
	else
	{
		fmt=( char * ) "ATOM  %5d %-5s%3s %c%4d%c   %8.3lf%8.3lf%8.3lf%8.4lf%8.4lf\n";
		fmt2= ( char * ) "HETATM%5d %-5s%3s  %4d    %8.3lf%8.3lf%8.3lf%8.4lf%8.4lf\n";
	}

	file=fopen(name,"at");
	if(file==NULL)
		return false;

	//	COLUMNS       DATA TYPE      FIELD         DEFINITION
	//	----------------------------------------------------------------------
	//	 1 -  6       Record name    "MODEL "
	//	11 - 14       Integer        serial        Model serial number.

	sprintf( line, "MODEL   %6d\n",n);
	fputs( line, file );

	initAll();
	contP = true;
	while ( contP != false )
	{
		p = ( Molecule * ) this->getCurrent();

		if(p->getClass()==pdb_smol)
		{
			contA = true;
			while ( contA != false )
			{
				a = ( Atom * ) p->getCurrent();

				/* Escritura de un atomo */
				a->getPosition(coor);
				if(number)
					sprintf( line, fmt2, a->getPdbSerial(), a->getPdbName(), p->getName(),
							p->getIdNumber(), coor[0], coor[1], coor[2], a->getPdbocc(), a->getPdbfact() );
				else
					sprintf( line, fmt2, cont_at++, a->getPdbName(), p->getName(),
							((SMol*)p)->get_pos()+1, coor[0], coor[1], coor[2], a->getPdbocc(), a->getPdbfact() );
				if(cont_at==100000)
					cont_at=1;
				fputs( line, file );
				contA = p->next();
			}
		}
		else
		{
			contCh = true;
			while ( contCh != false )
			{
				cad = ( Chain * ) p->getCurrent();
				cont_at=1;
				contF = true;
				while ( contF != false )
				{
					seg = ( Segment * ) cad->getCurrent();

					contR = true;
					while ( contR != false )
					{
						res = ( Fragment * ) seg->getCurrent();
						res_no++;

						contA = true;
						while ( contA != false )
						{
							a = ( Atom * ) res->getCurrent();

							/* Escritura de un atomo */
							a->getPosition(coor);

							if(number)
								sprintf( line, fmt, a->getPdbSerial(), a->getPdbName(), res->getName(),
										cad->getName() [0],res->getIdNumber()/*res->get_pos()*/,res->get_letter(), coor[0], coor[1], coor[2], a->getPdbocc(), a->getPdbfact() );
							else
								sprintf( line, fmt, cont_at++, a->getPdbName(), res->getName(),
										cad->getName() [0],res->get_pos()+1,res->get_letter(), coor[0], coor[1], coor[2], a->getPdbocc(), a->getPdbfact() );

							if(cont_at==100000)
								cont_at=1;
							fputs( line, file );
							contA = res->next();
						}
						contR = seg->next();
					}
					contF = cad->next();
				}
				contCh = p->next();
			}

			if(number)
				sprintf( line, "TER    %4d      %s %c%4d                                                      \n",
						a->getPdbSerial() + 1, res->getName(),
						cad->getName() [0], res->getIdNumber());
			else
				sprintf( line, "TER    %4d      %s %c%4d                                                      \n",
						a->getPdbSerial() + 1, res->getName(),
						cad->getName() [0], res->get_pos()+1);
			fputs( line, file );
		}
		contP = this->next();
	}


	sprintf( line, "ENDMDL\n");
	fputs( line, file );

	fclose(file);
	return true;
}


bool Macromolecule::writePDB( char *name, int *list, int total, bool number)
{
	FILE *file;
	bool contP, contCh, contF, contR, contA;
	char line[83];
	char * fmt;
	char * fmt2;
	int res_no = 0;
	int cont=0;
	int i;
	bool Found;

	Segment * seg;
	Chain * cad;
	Molecule * p;
	Fragment * res;
	Atom * a;
	int cont_at=1;

	Tcoor coor;
	char *ext;
	ext=strrchr(name,'.');


	if(strcmp(ext,".pqr")!=0)
	{
		fmt=( char * ) "ATOM  %5d %-5s%3s %c%4d%c   %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n";
		fmt2= ( char * ) "HETATM%5d %-5s%3s  %4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n";
	}
	else
	{
		fmt=( char * ) "ATOM  %5d %-5s%3s %c%4d%c   %8.3lf%8.3lf%8.3lf%8.4lf%8.4lf\n";
		fmt2= ( char * ) "HETATM%5d %-5s%3s  %4d    %8.3lf%8.3lf%8.3lf%8.4lf%8.4lf\n";
	}


	file=fopen(name,"wt");
	if(file==NULL)
		return false;

	initAll();
	contP = true;
	while ( contP != false )
	{
		p = ( Molecule * ) this->getCurrent();

		if(p->getClass()==pdb_smol)
		{
			contA = true;
			while ( contA != false )
			{
				a = ( Atom * ) p->getCurrent();

				/* Escritura de un atomo */
				Found=false;
				for(i=0;i<total && !Found;i++)
					if(list[i]==cont)
						Found=true;

				if(Found)
				{

					/* Escritura de un atomo */
					a->getPosition(coor);
					if(number)
						sprintf( line, fmt2, a->getPdbSerial(), a->getPdbName(), p->getName(),
								p->getIdNumber(), coor[0], coor[1], coor[2], a->getPdbocc(), a->getPdbfact() );
					else
						sprintf( line, fmt2, cont_at++, a->getPdbName(), p->getName(),
								((SMol*)p)->get_pos()+1, coor[0], coor[1], coor[2], a->getPdbocc(), a->getPdbfact() );
					if(cont_at==100000)
						cont_at=1;

					fputs( line, file );
				}
				contA = p->next();
			}
			cont++;
		}
		else
		{

			contCh = true;
			while ( contCh != false )
			{
				cad = ( Chain * ) p->getCurrent();
				cont_at=1;
				contF = true;
				while ( contF != false )
				{
					seg = ( Segment * ) cad->getCurrent();

					contR = true;
					while ( contR != false )
					{
						res = ( Fragment * ) seg->getCurrent();
						res_no++;

						contA = true;
						while ( contA != false )
						{
							a = ( Atom * ) res->getCurrent();

							/* Escritura de un atomo */
							Found=false;
							for(i=0;i<total && !Found;i++)
								if(list[i]==cont)
									Found=true;

							if(Found)
							{
								a->getPosition(coor);

								if(number)
									sprintf( line, fmt, a->getPdbSerial(), a->getPdbName(), res->getName(),
											cad->getName() [0],res->getIdNumber()/*res->get_pos()*/,res->get_letter(), coor[0], coor[1], coor[2], a->getPdbocc(), a->getPdbfact() );
								else
									sprintf( line, fmt, cont_at++, a->getPdbName(), res->getName(),
											cad->getName() [0],res->get_pos()+1,res->get_letter(), coor[0], coor[1], coor[2], a->getPdbocc(), a->getPdbfact() );

								if(cont_at==100000)
									cont_at=1;
								fputs( line, file );

							}

							contA = res->next();
						}
						contR = seg->next();
						cont++;
					}
					contF = cad->next();
				}
				contCh = p->next();
			}


			if(number)
				sprintf( line, "TER    %4d      %s %c%4d                                                      \n",
						a->getPdbSerial() + 1, res->getName(),
						cad->getName() [0], res->getIdNumber());
			else
				sprintf( line, "TER    %4d      %s %c%4d                                                      \n",
						a->getPdbSerial() + 1, res->getName(),
						cad->getName() [0], res->get_pos()+1);
			fputs( line, file );
		}
		contP = this->next();
	}


	fclose(file);
	return true;
}

bool Macromolecule::writePDB_atoms( char *name, int *list, int total, bool sign,bool number)
{
	FILE *file;
	bool contP, contCh, contF, contR, contA;
	char line[83];
	char * fmt;
	char * fmt2;
	int res_no = 0;
	int cont=0;
	int i;
	bool Found;

	Segment * seg;
	Chain * cad;
	Molecule * p;
	Fragment * res;
	Atom * a;
	Tcoor coor;
	int cont_at=1;

	char *ext;
	ext=strrchr(name,'.');


	if(strcmp(ext,".pqr")!=0)
	{
		fmt=( char * ) "ATOM  %5d %-5s%3s %c%4d%c   %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n";
		fmt2= ( char * ) "HETATM%5d %-5s%3s  %4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n";
	}
	else
	{
		fmt=( char * ) "ATOM  %5d %-5s%3s %c%4d%c   %8.3lf%8.3lf%8.3lf%8.4lf%8.4lf\n";
		fmt2= ( char * ) "HETATM%5d %-5s%3s  %4d    %8.3lf%8.3lf%8.3lf%8.4lf%8.4lf\n";
	}


	file=fopen(name,"wt");
	if(file==NULL)
		return false;

	initAll();
	contP = true;
	while ( contP != false )
	{
		p = ( Molecule * ) this->getCurrent();

		if(p->getClass()==pdb_smol)
		{
			contA = true;
			while ( contA != false )
			{
				a = ( Atom * ) p->getCurrent();

				if(sign==true)
				{
					Found=false;
					for(i=0;i<total && !Found;i++)
						if(list[i]==cont)
							Found=true;
				}
				else
				{
					Found=true;
					for(i=0;i<total && Found;i++)
						if(list[i]==cont)
							Found=false;
				}

				if(Found)
				{

					/* Escritura de un atomo */
					a->getPosition(coor);
					if(number)
						sprintf( line, fmt2, a->getPdbSerial(), a->getPdbName(), p->getName(),
								p->getIdNumber(), coor[0], coor[1], coor[2], a->getPdbocc(), a->getPdbfact() );
					else
						sprintf( line, fmt2, cont_at++, a->getPdbName(), p->getName(),
								((SMol*)p)->get_pos()+1, coor[0], coor[1], coor[2], a->getPdbocc(), a->getPdbfact() );
					if(cont_at==100000)
						cont_at=1;
					fputs( line, file );
				}
				contA = p->next();
				cont++;
			}

		}
		else
		{

			contCh = true;
			while ( contCh != false )
			{
				cad = ( Chain * ) p->getCurrent();
				cont_at=1;

				contF = true;
				while ( contF != false )
				{
					seg = ( Segment * ) cad->getCurrent();

					contR = true;
					while ( contR != false )
					{
						res = ( Fragment * ) seg->getCurrent();
						res_no++;


						contA = true;
						while ( contA != false )
						{
							a = ( Atom * ) res->getCurrent();

							/* Escritura de un atomo */
							if(sign==true)
							{
								Found=false;
								for(i=0;i<total && !Found;i++)
									if(list[i]==cont)
										Found=true;
							}
							else
							{
								Found=true;
								for(i=0;i<total && Found;i++)
									if(list[i]==cont)
										Found=false;
							}

							if(Found)
							{
								a->getPosition(coor);

								if(number)
									sprintf( line, fmt, a->getPdbSerial(), a->getPdbName(), res->getName(),
											cad->getName() [0],res->getIdNumber()/*res->get_pos()*/,res->get_letter(), coor[0], coor[1], coor[2], a->getPdbocc(), a->getPdbfact() );
								else
									sprintf( line, fmt, cont_at++, a->getPdbName(), res->getName(),
											cad->getName() [0],res->get_pos()+1,res->get_letter(), coor[0], coor[1], coor[2], a->getPdbocc(), a->getPdbfact() );
								if(cont_at==100000)
									cont_at=1;
								fputs( line, file );

							}

							contA = res->next();
							cont++;
						}
						contR = seg->next();

					}
					contF = cad->next();
				}
				contCh = p->next();
			}
			if(number)
				sprintf( line, "TER    %4d      %s %c%4d                                                      \n",
						a->getPdbSerial() + 1, res->getName(),
						cad->getName() [0], res->getIdNumber());
			else
				sprintf( line, "TER    %4d      %s %c%4d                                                      \n",
						a->getPdbSerial() + 1, res->getName(),
						cad->getName() [0], res->get_pos()+1);
			fputs( line, file );    fputs( line, file );
		}
		contP = this->next();
	}

	fclose(file);
	return true;
}




bool Macromolecule::writePDB(Macromolecule *mol2, char *name, int **list, int total, bool number)
{
	FILE *file;
	bool contP, contCh, contF, contR, contA;
	char line[83];
	char * fmt;
	char * fmt2;
	int res_no = 0;
	int cont=0;
	int i;
	bool Found;

	Segment * seg;
	Chain * cad;
	Molecule * p;
	Fragment * res;
	Atom * a;
	Tcoor coor;
	int cont_at=1;
	char *ext;
	ext=strrchr(name,'.');

	if(strcmp(ext,".pqr")!=0)
	{
		fmt=( char * ) "ATOM  %5d %-5s%3s %c%4d%c   %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n";
		fmt2= ( char * ) "HETATM%5d %-5s%3s  %4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n";
	}
	else
	{
		fmt=( char * ) "ATOM  %5d %-5s%3s %c%4d%c   %8.3lf%8.3lf%8.3lf%8.4lf%8.4lf\n";
		fmt2= ( char * ) "HETATM%5d %-5s%3s  %4d    %8.3lf%8.3lf%8.3lf%8.4lf%8.4lf\n";
	}

	file=fopen(name,"wt");
	if(file==NULL)
		return false;

	//PRIMERA MOLECULA
	initAll();
	contP = true;
	while ( contP != false )
	{
		p = ( Molecule * ) this->getCurrent();
		fprintf(stderr,"%s\n",p->getName());
		if(p->getClass()==pdb_smol)
		{
			contA = true;
			while ( contA != false )
			{
				a = ( Atom * ) p->getCurrent();

				if(list!=NULL)
				{
					Found=false;
					for(i=0;i<total && !Found;i++)
						if(list[0][i]==cont)
							Found=true;
				}
				else Found=true;

				if(Found)
				{
					/* Escritura de un atomo */
					a->getPosition(coor);
					if(number)
						sprintf( line, fmt2, a->getPdbSerial(), a->getPdbName(), p->getName(),
								p->getIdNumber(), coor[0], coor[1], coor[2], a->getPdbocc(), a->getPdbfact() );
					else
						sprintf( line, fmt2, cont_at++, a->getPdbName(), p->getName(),
								((SMol*)p)->get_pos()+1, coor[0], coor[1], coor[2], a->getPdbocc(), a->getPdbfact() );

					fputs( line, file );
					if(cont_at==100000)
						cont_at=1;
				}
				contA = p->next();
			}
			cont++;
		}
		else
		{


			contCh = true;
			while ( contCh != false )
			{
				cad = ( Chain * ) p->getCurrent();
				cont_at=1;
				contF = true;
				while ( contF != false )
				{
					seg = ( Segment * ) cad->getCurrent();

					contR = true;
					while ( contR != false )
					{
						res = ( Fragment * ) seg->getCurrent();
						res_no++;

						contA = true;
						while ( contA != false )
						{
							a = ( Atom * ) res->getCurrent();

							/* Escritura de un atomo */
							if(list!=NULL)
							{
								Found=false;
								for(i=0;i<total && !Found;i++)
									if(list[0][i]==cont)
										Found=true;
							}
							else Found=true;

							if(Found)
							{
								a->getPosition(coor);

								if(number)
									sprintf( line, fmt, a->getPdbSerial(), a->getPdbName(), res->getName(),
											cad->getName() [0],res->getIdNumber()/*res->get_pos()*/,res->get_letter(), coor[0], coor[1], coor[2], a->getPdbocc(), a->getPdbfact() );
								else
									sprintf( line, fmt, cont_at++, a->getPdbName(), res->getName(),
											cad->getName() [0],res->get_pos()+1,res->get_letter(), coor[0], coor[1], coor[2], a->getPdbocc(), a->getPdbfact() );

								if(cont_at==100000)
									cont_at=1;
								fputs( line, file );

							}

							contA = res->next();
						}
						contR = seg->next();
						cont++;
					}
					contF = cad->next();
				}
				contCh = p->next();
			}

			if(number)
				sprintf( line, "TER    %4d      %s %c%4d                                                      \n",
						a->getPdbSerial() + 1, res->getName(),
						cad->getName() [0], res->getIdNumber());
			else
				sprintf( line, "TER    %4d      %s %c%4d                                                      \n",
						a->getPdbSerial() + 1, res->getName(),
						cad->getName() [0], res->get_pos()+1);
			fputs( line, file );
		}
		contP = this->next();
	}



	//SEGUNDA MOLECULA
	cont=0;
	cont_at=1;
	mol2->initAll();
	contP = true;
	while ( contP != false )
	{
		p = ( Molecule * ) mol2->getCurrent();
		fprintf(stderr,"%s\n",p->getName());

		if(p->getClass()==pdb_smol)
		{
			contA = true;
			while ( contA != false )
			{
				a = ( Atom * ) p->getCurrent();

				if(list!=NULL)
				{
					Found=false;
					for(i=0;i<total && !Found;i++)
						if(list[0][i]==cont)
							Found=true;
				}
				else Found=true;

				if(Found)
				{

					/* Escritura de un atomo */
					a->getPosition(coor);
					if(number)
						sprintf( line, fmt2, a->getPdbSerial(), a->getPdbName(), p->getName(),
								p->getIdNumber(), coor[0], coor[1], coor[2], a->getPdbocc(), a->getPdbfact() );
					else
						sprintf( line, fmt2, cont_at++, a->getPdbName(), p->getName(),
								((SMol*)p)->get_pos()+1, coor[0], coor[1], coor[2], a->getPdbocc(), a->getPdbfact() );
					if(cont_at==100000)
						cont_at=1;
					fputs( line, file );
				}
				contA = p->next();
			}
			cont++;

		}
		else
		{


			contCh = true;
			while ( contCh != false )
			{
				cad = ( Chain * ) p->getCurrent();
				cont_at=1;
				contF = true;
				while ( contF != false )
				{
					seg = ( Segment * ) cad->getCurrent();

					contR = true;
					while ( contR != false )
					{
						res = ( Fragment * ) seg->getCurrent();
						res_no++;

						contA = true;
						while ( contA != false )
						{
							a = ( Atom * ) res->getCurrent();

							/* Escritura de un atomo */
							if(list!=NULL)
							{
								Found=false;
								for(i=0;i<total && !Found;i++)
									if(list[1][i]==cont)
										Found=true;
							}
							else Found=true;

							if(Found)
							{
								a->getPosition(coor);

								if(number)
									sprintf( line, fmt, a->getPdbSerial(), a->getPdbName(), res->getName(),
											cad->getName() [0],res->getIdNumber()/*res->get_pos()*/,res->get_letter(), coor[0], coor[1], coor[2], a->getPdbocc(), a->getPdbfact() );
								else
									sprintf( line, fmt, cont_at++, a->getPdbName(), res->getName(),
											cad->getName() [0],res->get_pos()+1,res->get_letter(), coor[0], coor[1], coor[2], a->getPdbocc(), a->getPdbfact() );
								if(cont_at==100000)
									cont_at=1;

								fputs( line, file );

							}

							contA = res->next();
						}
						contR = seg->next();
						cont++;
					}
					contF = cad->next();
				}
				contCh = p->next();
			}
			if(number) {
				sprintf( line, "TER    %4d      %s %c%4d                                                      \n",
						a->getPdbSerial() + 1, res->getName(),
						cad->getName() [0], res->getIdNumber());
			}
			else {
				sprintf( line, "TER    %4d      %s %c%4d                                                      \n",
						a->getPdbSerial() + 1, res->getName(),
						cad->getName() [0], res->get_pos()+1);
			}
			fputs( line, file );

		}

		contP = mol2->next();
	}


	sprintf( line, "END\n");
	fputs( line, file );

	fclose(file);
	return true;
}

bool Macromolecule::readSDF( char *name)
{
	FILE *file;
	file=fopen(name,"rt");
	if(file==NULL) {
		fprintf(stderr,"Can't open %s for reading: No such sdf file\n",name);
		exit(1);
	}
	int step=0;
	char aux1[100],aux2[100],line[500];
	float pos[3];
	int cont=0, cont_at=0,cont_bonds=0,cont_mol=1;
	int num_atoms, num_bonds;
	char type_aux[3],type_aux2[6];
	int init,final,link,i;

	SMol *smol;
	Atom *at,*at1,*at2;
	Bond *bond;

	step=0;

	for(i=0;i<5;i++)
	{
		type_aux2[i]=' ';
	}
	type_aux2[5]='\0';

	while ( fgets( line, 82, file ) != 0 )
	{
		//fprintf(stderr,"{%s}",line);
		sscanf(line,"%s %s\n",aux1,aux2);
		//fprintf(stderr,"{%s} {%s}\n",aux1,aux2);

		if(strcmp(aux1,"$$$$")==0 || (strcmp(aux1,"M")==0 && strcmp(aux2,"END")==0  ) )
		{
			//fprintf(stderr,"hola\n");
			this->add(smol);
			step=0;


		}
		else
		{
			switch(step)
			{
			case 0:
				sscanf(line,"%s %*s\n",aux1);

				smol=new SMol(aux1, cont_mol, cont_mol);
				cont_mol++;
				cont=0;
				step=1;
				break;
			case 1:
				cont++;
				if(cont==2)
					step=2;
				break;
			case 2:
				sscanf(line," %d %d %*s\n",&num_atoms, &num_bonds);
				//fprintf(stderr,"num_atoms=%d num_bonds=%d\n",num_atoms,num_bonds);
				cont=0;
				step=3;
				break;
			case 3:
				sscanf(line,"   %f %f %f %s %*s\n",&(pos[0]),&(pos[1]),&(pos[2]),type_aux);
				fprintf(stderr,"%f|%f|%f %s\n",pos[0],pos[1],pos[2],type_aux);

				for(i=0;i<3;i++)
				{

					if(type_aux[i]!='\0')
						type_aux2[i]=type_aux[i];
					else
					{
						type_aux2[i]=' ';
						i=3;
					}


				}
				fprintf(stderr, "{%s}\n",type_aux2);
				at = new Atom( Table_Elements::getElement( type_aux ), pos, 0.0, type_aux2, cont_at, 0.0, 0.0 );
				smol->add(at);
				cont_at++;
				cont++;
				if(cont==num_atoms)
				{
					step=4;
					cont_bonds=0;
				}
				break;
			case 4:
				sscanf(line," %d %d %d %*d %*d %*d %*d\n",&init,&final,&link);
				at1=(Atom*)smol->getE(init-1);
				at2=(Atom*)smol->getE(final-1);

				bond=new Bond(at1,at2,link);
				at1->insertBond(bond);
				at2->insertBond(bond);
				smol->insertBond(bond);
				cont_bonds++;
				if(cont_bonds==num_bonds)
				{
					step=5;
					cont_bonds=0;
				}
				break;

			default:
				break;
			}
		}
	}
	fclose(file);

	return true;
}

int return_pos_aux(SMol *smol,Atom *at);

void Macromolecule::writeSDF( char *name)
{
	FILE *file;
	file=fopen(name,"wt");

	SMol *smol;
	Atom *at;
	Bond *bond;
	float pos[3];
	int i;

	initAll();

	do
	{
		if(((Molecule*)(this->getCurrent()))->getMolType()==tmol_smol)
		{
			smol=(SMol*)(this->getCurrent());
			fprintf(file,"%s\nFrom libPDB\n\n",smol->getName());
			fprintf(file,"%3d%3d  0  0  0  0  0  0  0\n",smol->getLimit(),smol->get_numBonds());

			do{
				at=(Atom*)smol->getCurrent();
				at->getPosition(pos);
				fprintf(file,"%10.4f%10.4f%10.4f %-4s 0  0  0  0  0  1  0  0  0  0  0  0\n",pos[0],pos[1],pos[2],at->getName());

			}while(smol->next());

			for(i=0;i<smol->get_numBonds();i++)
			{
				bond=smol->getBond(i);
				fprintf(file,"%3d%3d%3d  0  0  0  0 \n",return_pos_aux(smol,bond->getInit()),return_pos_aux(smol,bond->getFinal()),bond->getLink());
			}

			fprintf(file,"\n$$$$\n");

		}
	}while(next());
	fclose(file);
}

bool Macromolecule::readMol2(char *name)
{
	FILE * f=fopen(name,"rt");
	if(f==NULL) {
		fprintf(stderr,"Can't open %s for reading: No such mol2 file\n",name);
		exit(1);
	}
	char line[500],aux1[100],aux2[100];
	int step=0;
	int num_atoms,num_bonds,num_struct;
	int contMol=1,cont_at=0,i,j;
	SMol *smol;
	bool issmol=false;
	float pos[3],charge;
	char type_aux[3],type_aux2[7]; //estaba diferente, asi type_aux2[6]
	type_aux[2]='\0';
	Atom *at,*at1,*at2;
	int init,final,link;
	Bond *bond;

	for(i=0;i<5;i++)
	{
		type_aux2[i]=' ';
	}
	type_aux2[5]='\0';

	while ( fgets( line, 82, f ) != 0 )
	{
		step=0;
		if(line[0]=='@')
		{
			sscanf(line,"@<TRIPOS>%s\n",aux1);

			//fprintf(stderr,"TRIPOS->%s\n",aux1);




			if(strcmp(aux1,"MOLECULE")==0)
			{
				step=1;
				issmol=false;

				if (contMol!=1) {
					this->add(smol);
				}
				smol=new SMol(aux1,contMol,contMol);
				contMol++;
			}
			else
				if(strcmp(aux1,"ATOM")==0)
				{
					step=2;

				}
				else
					if(strcmp(aux1,"BOND")==0)
					{
						step=3;
					}
		}




		switch(step)
		{

		case 1:

			fscanf(f,"%s\n",aux1);
			// fprintf(stderr,"Molecule:%s ",aux1);


			fgets( line, 82, f );
			//fprintf(stderr,"Molecule:END\n");

			sscanf(line,"%d %d %d",&num_atoms,&num_bonds,&num_struct);
			fscanf(f,"%s\n%*s",aux2);

			//fprintf(stderr," n %d:%s %d %d %d %s\n",contMol,aux1,num_atoms,num_bonds,num_struct,aux2);


			if(strcmp(aux2,"SMALL")==0)
			{
				issmol=true;
			}
			else
			{
				issmol=false;
			}
			break;

		case 2:
			if(issmol)
			{
				cont_at=0;
				for(i=0;i<num_atoms;i++)
				{
					int k=0; //nuevo Erney
					fscanf(f,"%*s %*s %f %f %f  %s %d %s %f\n",  &(pos[0]),&(pos[1]),&(pos[2]),aux1,&num_struct,aux2,&charge);
					//fprintf(stderr,"Atom:\n\t%s %f %f %f %d %s %f\n", aux1, pos[0],pos[1],pos[2],num_struct,aux2,charge);
					type_aux[0]=aux1[0];
					if(!isalpha(aux1[1])) type_aux[1]=' ';
					else if(islower(aux1[1]))
						type_aux[1]=toupper(aux1[1]);
					else type_aux[1]=aux1[1];
					//fprintf(stderr,"%s %s\n", type_aux,type_aux2 );
					type_aux[2]='\0';
					for(j=0;j<6;j++) //puse 6 pero iba un 5 en el original, donde se puso la k iba la j Erney
					{
						if (aux1[j]!='.')
						{
							if (aux1[j]!='\0')
								type_aux2[k]=aux1[j];
							else
							{
								type_aux2[k]='\0';
								j=7;
							}
							k++;//nuevo Erney
						}

					}
					//fprintf(stderr,"[%s] [%s]\n", type_aux,type_aux2 );

					at = new Atom( Table_Elements::getElement( type_aux ), pos, 0.0, type_aux2, cont_at, 0.0, charge );
					//fprintf(stderr,"Element: %s :[%s] [%s]\n", at->getElement(), type_aux,type_aux2 );
					cont_at++;
					smol->add(at);
					// old_num_struct=num_struct;
				}


				//fprintf(stderr,"New mol %d %d\n", contMol,cont_at );


			}
			break;
		case 3:
			if(issmol)
			{
				for(i=0;i<num_bonds;i++)
				{
					fscanf(f,"%*d %d %d %s\n", &init,&final,aux1);

					if(strcmp(aux1,"1")==0)
						link=1;
					else
						if(strcmp(aux1,"2")==0)
							link=2;
						else
							if(strcmp(aux1,"3")==0)
								link=3;
							else
								if(strcmp(aux1,"am")==0)
									link=4;
								else
									if(strcmp(aux1,"ar")==0)
										link=5;
									else
										if(strcmp(aux1,"du")==0)
											link=6;
										else
											if(strcmp(aux1,"un")==0)
												link=7;
											else
												if(strcmp(aux1,"nc")==0)
													link=8;


					at1=(Atom*)smol->getE(init-1);
					at2=(Atom*)smol->getE(final-1);

					bond=new Bond(at1,at2,link);
					at1->insertBond(bond);
					at2->insertBond(bond);
					smol->insertBond(bond);

				}




			}
			break;



		}

	}

	// add the last one
	this->add(smol);

	return true;
	fclose(f);
}

void Macromolecule::writeMol2(char *name)
{
	FILE * f=fopen(name,"wt");
	bool contP = true;
	Molecule *p;
	Atom *at;
	Bond *bond;
	Tcoor pos;
	char aux_cad[3];
	int i,cont_at=1,cont_mol=1,cont_bond=1;
	this->initAll();

	while ( contP != false )
	{

		p = ( Molecule * ) this->getCurrent();
		if(p->getClass()==pdb_smol)
		{
			cont_at=1;
			fprintf(f,"@<TRIPOS>MOLECULE\n%s\n%d %d 0 0 0\nSMALL\nNO_CHARGES\n\n\n",p->getName(),p->getLimit(),((SMol*)p)->get_numBonds());
			fprintf(f,"@<TRIPOS>ATOM\n");
			for(i=0;i<p->getLimit();i++)
			{
				at=(Atom*)p->getCurrent();
				at->getPosition(pos);
				fprintf(f,"%6d  %-5s  %10.3f%10.3f%10.3f   %-5s%6d %s%12.3f\n",cont_at,at->getName(),pos[0],pos[1],pos[2],at->getPdbName(),cont_mol,p->getName(),at->getPdbfact());
				cont_at++;
				p->next();
			}


			if(((SMol*)p)->get_numBonds()>0)
			{
				fprintf(f,"@<TRIPOS>BOND\n");
				cont_bond=1;
				for(i=0;i<((SMol*)p)->get_numBonds();i++)
				{
					bond=((SMol*)p)->getBond(i);

					switch(bond->getLink())
					{
					case 1:
						strcpy(aux_cad,"1 ");
						break;
					case 2:
						strcpy(aux_cad,"2 ");
						break;
					case 3:
						strcpy(aux_cad,"3 ");
						break;
					case 4:
						strcpy(aux_cad,"am");
						break;
					case 5:
						strcpy(aux_cad,"ar");
						break;
					case 6:
						strcpy(aux_cad,"du");
						break;
					case 7:
						strcpy(aux_cad,"un");
						break;
					case 8:
						strcpy(aux_cad,"nc");
						break;
					}

					fprintf(f,"%5d%5d%5d %s\n",cont_bond, return_pos_aux(((SMol*)p),bond->getInit()),return_pos_aux(((SMol*)p),bond->getFinal()),aux_cad);
					cont_bond++;
				}
			}

			cont_mol++;


		}

		contP=this->next();
	}


	fclose(f);
}


int return_pos_aux(SMol *smol,Atom *at)
{

	int i;
	Atom *at2;
	for(i=0;smol->getLimit();i++)
	{
		at2=(Atom*)smol->getE(i);
		if (at2==at)
			return i+1;
	}
	return -1;
}

int pdb_strcasecmp(const char *s1, const char *s2)
{
	while ((*s1 != '\0')
			&& (tolower(*(unsigned char *) s1) ==
					tolower(*(unsigned char *) s2))) {
		s1++;
		s2++;
	}

	return tolower(*(unsigned char *) s1) - tolower(*(unsigned char *) s2);
}
