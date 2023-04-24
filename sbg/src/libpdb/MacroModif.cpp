
#include "Macromolecule.h"
#include "pdbIter.h"
// #include "ResIni.h"

void Macromolecule::centerBox()
{
  geoBox();
  macromoleculeDimension.geometricCenter[0] *= -1.0;
  macromoleculeDimension.geometricCenter[1] *= -1.0;
  macromoleculeDimension.geometricCenter[2] *= -1.0;

  moveAll( macromoleculeDimension.geometricCenter );
  geoBox();
}

bool Macromolecule::applyAtoms( AtomOperator * op )
{
  bool end = true;

  bool out = true;
  Atom * at;


  initAll();
  while ( end != false )
  {
    /* Obtengo un atomo */
    at = getCurrentAtom();
    /* Aplicacion de la operacion*/
    if ( !op->apply( at ) )
     out = false;

    end = nextAtom();
  }


  return out;
}

bool Macromolecule::applyAtoms(AtomOperator *op, Macromolecule *ref)
{
  bool end = true;
  bool out = true;
  Atom * at,* at2;

  initAll();
  ref->initAll();
  while ( end != false )
  {
    /* Obtengo un atomo del sistema a modificar*/
    at = getCurrentAtom();
    //Obtengo atomo de referencia
    at2=ref->getCurrentAtom();
    /* Aplicacion de la operacion*/
    if ( !op->apply( at2,at ) )
     out = false;

    end = nextAtom();
    ref->nextAtom();
  }

  return out;
}

bool Macromolecule::copy_coordinates(Macromolecule *ref)
{
  Atom * at,* at2;
  bool end=true;
  Tcoor pos;

  if(ref->get_num_atoms()!=get_num_atoms())
	  return false;

  initAll();
  ref->initAll();
  while ( end != false )
  {
    /* Obtengo un atomo del sistema a modificar*/
    at =ref->getCurrentAtom();
    at->getPosition(pos);

    at2 =getCurrentAtom();
    at2->setPosition(pos);


    end = nextAtom();
    ref->nextAtom();
  }

  return true;
}

/**
 * Sets atomic positions from the provided array.
 * @param coords: Coordinates array
 */
bool Macromolecule::copy_coordinates(double *coords)
{
	Tcoor pos;
	pdbIter *iter = new pdbIter(this);

	for(iter->pos_atom=0; !iter->gend_atom(); iter->next_atom())
	{
		pos[0] = coords[iter->pos_atom*3];
		pos[1] = coords[iter->pos_atom*3 + 1];
		pos[2] = coords[iter->pos_atom*3 + 2];
		(iter->get_atom())->setPosition(pos); // now "aver" is the initial average
	}
	delete iter;

	return true;
}

void Macromolecule::format_residues(int DETAIL)
{
 Fragment * res;
 int j,resn;
 char * at_name;
 Segment * seg;
 Atom * atom;
 pdbIter *iter1, *iter2;
 int atoms_for_checking;

  //Iterador para recorrer segmentos
  iter1 = new pdbIter( this );
  Chain *ch;
  //Bucle para recorrer segmentos
  for ( iter1->pos_segment = 0; !iter1->gend_segment(); iter1->next_segment() )
  {
    seg = ( Segment * ) iter1->get_segment();
    iter2 = new pdbIter( seg );
    ch = (Chain *) seg->getFather(); // ->get_chain();
    //Bucle para recorrer residuos
    for ( iter2->pos_fragment = 0; !iter2->gend_fragment(); iter2->next_fragment() )
    {
      res = ( Fragment * ) iter2->get_fragment();
      resn = resnum_from_resname( res->getName() );

      //fprintf( stderr, "Processing aa %s\n", res->getName());

      if ( resn > N_AMINO )
      {
        fprintf( stderr, "Warning %s residue is not an aa\n", res->getName());
        //getchar();
      }

      //Iterador sobre los atomos del residuo
      pdbIter * iter3 = new pdbIter( res );
      bool  atdetected; int menos;

       menos=0;
       switch(DETAIL)
     	{
     	case 0: // Level of detail: full atoms
     		atoms_for_checking=AA[resn].natoms;
     	break;
     	case 1: // Level of detail: backbone
     	    atoms_for_checking=4;
     	break;
     	case 2: // Level of detail: backbone + cb
     		atoms_for_checking=5;
     		if (resn==5) atoms_for_checking=4;

     	break;
     	case 3: // "Level of detail: full atoms without H
     	    atoms_for_checking=AA[resn].nheavyatoms;
     	break;
     	default:
     		atoms_for_checking=AA[resn].nheavyatoms;
     	};


       for ( j = 0; j < atoms_for_checking; j++ ) {
         atdetected=false;
         for ( iter3->pos_atom = 0; !iter3->gend_atom(); iter3->next_atom() )
          {
               atom = ( Atom * ) iter3->get_atom();
               at_name = atom->getName();

              if ( !strncmp( at_name, AA[resn].atom[j].atom_name, 4 ) ) {
                if ( iter3->pos_atom != j )
                res->exchange( iter3->pos_atom, j-menos );
                //fprintf( stderr, "  Warning atom %s of residue %s exchange position (%2d -> %2d) %d not detected\n",
                //at_name, res->getName(), iter3->pos_atom, j-menos, menos);

                atdetected=true;
                break;
             }
          } // iter 3
          iter3->~pdbIter();
          iter3 = new pdbIter( res );

          if (!atdetected) {
              fprintf( stderr, " Warning atom <%s> is not detected in aa %s res %d chain %s\n",
				  AA[resn].atom[j].atom_name, res->getName(), res->getIdNumber(), ch->getName());
				menos++;
            //getchar();
          }

       } // j-rosetta

    }
    iter2->~pdbIter();

  }
iter1->~pdbIter();

}

// Mon modified for iModS (2/10/2013)
// It returns the number of missing atoms (=0, NO missing atoms found)
int Macromolecule::format_residues(bool withH, int model)
{
	int j,resn;
	char * at_name;
	int atoms_for_checking;
	int missing=0; // Number of missing atoms detected
	Segment * seg;
	Fragment * res;
	Atom * atom;
	pdbIter *iter1, *iter2;

	//Iterador para recorrer segmentos
	iter1 = new pdbIter( this );

	//Bucle para recorrer segmentos
	for ( iter1->pos_segment = 0; !iter1->gend_segment(); iter1->next_segment() )
	{
		seg = ( Segment * ) iter1->get_segment();
		iter2 = new pdbIter( seg );

		//Bucle para recorrer residuos
		for ( iter2->pos_fragment = 0; !iter2->gend_fragment(); iter2->next_fragment() )
		{
			res = ( Fragment * ) iter2->get_fragment();
			// resn = resnum_from_resname( res->getName() ); // Now all fragments contain their ID number (Mon 21/02/2018)
			resn = res->fragid; // Directly get the ID from Fragment (it was set by constructor) (Mon 21/02/2018)
			//      fprintf( stderr, "Processing aa %s resn= %d\n", res->getName(), resn);

			if ( resn > N_AMINO )
			{
				fprintf( stderr, "Warning %s residue is not an aa\n", res->getName());
				//getchar();
			}

			if(model == 0) // CA-model
				atoms_for_checking = 3; // The only needed atoms in CA-model are: N, CA, C (indices 0,1,2)
			else if(model == 1) // C5-model
			{
				if(resn==5) // if GLY
					atoms_for_checking = 4; // The only needed atoms in GLY according to C5-model are: N, CA, C, O (indices 0,1,2,3)
				else
					atoms_for_checking = 5; // The only needed atoms in other aminoacids with C5-model are: N, CA, C, O, CB (indices 0,1,2,3,4)
			}
			else // HA-model
				if (withH==true)
					atoms_for_checking = AA[resn].natoms;
				else
					atoms_for_checking = AA[resn].nheavyatoms;

			//       fprintf( stderr, "resn= %d  nheavyatoms= %d  at_check= %d\n",resn,AA[resn].nheavyatoms,atoms_for_checking);
			//       fprintf( stderr, "resn= %d  natoms= %d  at_check= %d resatom= %d\n",resn,AA[resn].natoms,atoms_for_checking,iter3->num_atom());

			//Iterador sobre los atomos del residuo
			pdbIter * iter3 = new pdbIter( res );
			bool atdetected;
			int menos = 0;

			for ( j = 0; j < atoms_for_checking; j++ ) // Check all fragment atoms (from database) vs. current residue atoms
			{
				atdetected=false;
				for ( iter3->pos_atom = 0; !iter3->gend_atom(); iter3->next_atom() )
				{
					atom = ( Atom * ) iter3->get_atom();
					at_name = atom->getName();
					//               fprintf( stderr, " Atom \"%s\" of residue %s %d at_check= %d %d it3= %d\n",
					//                               at_name, res->getName(), resn, AA[resn].nheavyatoms, atoms_for_checking, iter3->pos_atom );

					if ( !strncmp( at_name, AA[resn].atom[j].atom_name, 4 ) )
					{
						if ( iter3->pos_atom != j )
							res->exchange( iter3->pos_atom, j-menos );
						//                fprintf( stderr, "  Warning atom %s of residue %s exchange position (%2d -> %2d), %d not detected -%s-\n",
						//                		at_name, res->getName(), iter3->pos_atom, j-menos, menos,AA[resn].atom[j].atom_name );

						atdetected=true;
						break;
					}
				} // iter 3
				//          iter3->~pdbIter();
				delete iter3; // CHECK THIS LATER, we should reset the iterator not re-creating it! NO, because atoms may have been exchanged...
				iter3 = new pdbIter( res );

				if (!atdetected)
				{
					//            fprintf( stderr, " Warning atom <%s> is not detected in aa %s, frag= %d\n",
					//            		AA[resn].atom[j].atom_name, res->getName(), iter2->pos_fragment );
					fprintf( stdout, "Warning, missing %s atom from residue %s %5d in chain %s\n", AA[resn].atom[j].atom_name, res->getName(), res->getIdNumber(), ((Chain *)seg->getFather())->getName() );
					menos++;
					missing++;
					//getchar();
				}

			} // j-rosetta
			delete iter3;
		}
		delete iter2;
		//    iter2->~pdbIter();
	}
	delete iter1;
	//  iter1->~pdbIter();
	return missing;
}

Dimensions * Macromolecule::getDimension()
{
  return & macromoleculeDimension;
}






bool Macromolecule::moveAll(Tcoor offset)
{
  int i;
  for(i=0;i<limit;i++)
      if(elements[i]->moveAll(offset)==false)
        return false;

 return true;
}

//Auxiliar function
void check_resnames(Fragment *r)
{
	int option=-1;
	r->initAll();
	Atom *a;
	bool atom_name1,atom_name2;
	if(strcmp(r->getName(),"ASP")==0)
		option=0;
	else
		if(strcmp(r->getName(),"CYS")==0)
			option=1;
		else
			if(strcmp(r->getName(),"GLU")==0)
				option=2;
			else
				if(strcmp(r->getName(),"HIS")==0)
					option=3;
				else
					if(strcmp(r->getName(),"LYS")==0)
						option=4;
					else
						if(strcmp(r->getName(),"TYR")==0)
							option=5;


	switch(option)
	{
		//ASH
		case 0:
		atom_name1=false;
		do{
			a=(Atom*)r->getCurrent();
			if(strcmp(a->getPdbName(),"1HD  ")==0
			|| strcmp(a->getPdbName(),"2HD  ")==0
		  || strcmp(a->getPdbName()," HD2 ")==0)
			{
				atom_name1=true;
				r->modifyName((char *)"ASH");
			}

		}while(r->next() && !atom_name1);
		break;

		//CYX
		case 1:
		atom_name1=false;
		do{
			a=(Atom*)r->getCurrent();
			if(strcmp(a->getPdbName(),(char *)" HG  ")==0)
			{
				atom_name1=true;
			}

		}while(r->next() && !atom_name1);
		if(!atom_name1)
			r->modifyName((char *)"CYX");
		break;

		//GLH
		case 2:
		atom_name1=false;
		do{
			a=(Atom*)r->getCurrent();
			if(strcmp(a->getPdbName(),"1HE  ")==0
			|| strcmp(a->getPdbName(),"2HE  ")==0
		  || strcmp(a->getPdbName()," HE2 ")==0)
			{
				atom_name1=true;
				r->modifyName((char *)"GLH");
			}

		}while(r->next() && !atom_name1);
		break;

	//HIP & HID
		case 3:

		atom_name1=false;
		atom_name2=false;
		do{
			a=(Atom*)r->getCurrent();
			if(strcmp(a->getPdbName()," HD1 ")==0)
			{
				atom_name1=true;
			}

			if(strcmp(a->getPdbName()," HE2 ")==0)
			{
				atom_name2=true;
			}

		}while(r->next());
		if(atom_name2 && atom_name1)
				r->modifyName((char *)"HIP");
		else
				if(!atom_name2 && atom_name1)
				{
						r->modifyName((char *)"HID");
				}

		break;

		//LYN
		case 4:
		atom_name1=false;
		do{
			a=(Atom*)r->getCurrent();
			if(strcmp(a->getPdbName()," HZ3 ")==0)
			{
				atom_name1=true;
			}

		}while(r->next() && !atom_name1);
		if(!atom_name1)
			r->modifyName((char *)"LYN");
		break;

		//LYN
		case 5:
		atom_name1=false;
		do{
			a=(Atom*)r->getCurrent();
			if(strcmp(a->getPdbName()," HH  ")==0)
			{
				atom_name1=true;
			}

		}while(r->next() && !atom_name1);
		if(!atom_name1)
			r->modifyName((char *)"TYM");
		break;
	}

}

void Macromolecule::rename_residues()
{
	pdbIter *iter = new pdbIter( this );

	for(iter->pos_fragment=0; !iter->gend_fragment(); iter->next_fragment())
	if((iter->get_fragment())->getMolType()!=tmol_smol)
	{
			check_resnames(iter->get_fragment());
	}
	delete iter;
}


// mut some sequence "seq" (in 1-letter format) from residue "index" (internal index) up to end of "seq"
void Macromolecule::mutseq(char seq, int index)
{
	pdbIter *iter = new pdbIter( this );
	Residue *res;


	int i; // Aminoacids table index



	for(iter->pos_fragment = 0; !iter->gend_fragment(); iter->next_fragment())
	{
		res = (Residue *) iter->get_fragment();

		if ((res->getIdNumber() == index) && (res->getMolType() != tmol_smol) )
		{
		    printf(" Warning res %s --> %c pos %d %d\n",res->getName(), seq, index, res->getIdNumber()  );

			for(i = 0; i < N_AMINO; i++)
			{
				if(AA[i].aa_name1 == seq)
				{
					res->modifyName( AA[i].aa_name3 );
					break;
				}
			}
		}
	}


	delete iter;
}



// Mon modified (5/4/2010)
// Convert to 3BB + 2R reduced model (Erases some atoms and changes values in others )
// If nomass = true, all masses set to 1.0 (default=false)
// If nocharge = true, all charges set to 1.0 (default=false)
int  Macromolecule::reducedmodel_3BBR2( bool nomass, bool nocharge)
{
  Residue * res;
  int j, resn, atoms, atoms_R, atom_res;
  char * at_name;
  Segment * seg;
  Atom * atom;
  Tcoor pos;
  pdbIter * iter1, * iter2;
  float cx, cy, cz; // center of masses
  int num_res=0;
  float masstotal, mass, massCH, massCH2,massCH3, massNH, massCO, imass;
  float chargetotal, charge, chargeCH, chargeCH2, chargeCH3, chargeNH, chargeCO, icharge;

  imass = 1.0; // when nomass enabled
  icharge = 1.0; // when nocharge enabled

  // Electronic Charges (more EM-like!)
  if(nocharge)
  {
	  chargeCH = 1.0;
	  chargeCH2 = 1.0;
	  chargeCH3 = 1.0;
	  chargeNH = 1.0;
	  chargeCO = 1.0;
  }
  else
  {
	  chargeCH = 6 + 1;
	  chargeCH2 = chargeCH + 1;
	  chargeCH3 = chargeCH2 + 1;
	  chargeNH = 7 + 1;
	  chargeCO = 6 + 8;
  }

  // Masses
  if(nomass)
  {
	  massCH = 1.0;
	  massCH2 = 1.0;
	  massCH3 = 1.0;
	  massNH = 1.0;
	  massCO = 1.0;
  }
  else
  {
	  massCH = 12.011 + 1.00797;
	  massCH2 = massCH + 1.00797;
	  massCH3 = massCH2 + 1.00797;
	  massNH = 14.00674 + 1.00797;
	  massCO = 12.011 + 15.9994;
  }

  //Iterador para recorrer segmentos
  iter1 = new pdbIter( this,false,false );

  //Bucle para recorrer segmentos
  for ( iter1->pos_segment = 0; !iter1->gend_segment(); iter1->next_segment() )
  {
    seg = ( Segment * ) iter1->get_segment();
    if(seg->getMolType()==tmol_protein)
    {

		iter2 = new pdbIter( seg );

		//Bucle para recorrer residuos
		for ( iter2->pos_fragment = 0; !iter2->gend_fragment(); iter2->next_fragment() )
		{
		  res = ( Residue * ) iter2->get_fragment();
		  resn = resnum_from_resname( res->getName() );
		  num_res++;
		  //Iterador sobre los atomos del residuo
		  pdbIter * iter3 = new pdbIter( res );

		  switch (resn)
		  {
			case GLY:
			if (iter3->num_atom() < 3)
			{
			  printf( " Residue  %s%d only %d atoms\n", res->getName(),num_res, iter3->num_atom());
			  exit(1);
			}
			iter3->pos_atom = 0; //NH
			(( Atom * ) iter3->get_atom())->setPdbocc(massNH);
			(( Atom * ) iter3->get_atom())->setPdbfact(chargeNH);

			iter3->next_atom(); // CA+HA
			(( Atom * ) iter3->get_atom())->setPdbocc(massCH);
			(( Atom * ) iter3->get_atom())->setPdbfact(chargeCH);

			iter3->next_atom(); // CO
			(( Atom * ) iter3->get_atom())->setPdbocc(massCO);
			(( Atom * ) iter3->get_atom())->setPdbfact(chargeCO);

			// Deleting remaining un-necessary atoms
			for ( iter3->pos_atom = 3; !iter3->gend_atom(); iter3->next_atom() )
				res->remove(3);
			break;

			case ALA:
			if (iter3->num_atom() < 4)
			{
			  printf( " Residue  %s%d only %d atoms\n", res->getName(),num_res, iter3->num_atom());
			  exit(1);
			}
			iter3->pos_atom = 0; //NH
			(( Atom * ) iter3->get_atom())->setPdbocc(massNH);
			(( Atom * ) iter3->get_atom())->setPdbfact(chargeNH);

			iter3->next_atom(); // CA+HA
			(( Atom * ) iter3->get_atom())->setPdbocc(massCH);
			(( Atom * ) iter3->get_atom())->setPdbfact(chargeCH);

			iter3->next_atom(); // CO
			(( Atom * ) iter3->get_atom())->setPdbocc(massCO);
			(( Atom * ) iter3->get_atom())->setPdbfact(chargeCO);

			iter3->next_atom(); // CB
			res->exchange( 3, 4 );   // exchange O for CB

			iter3->pos_atom = 4; // O in the Macromolecule struct
			iter3->get_atom()->setPdbocc(massCH3);
			iter3->get_atom()->setPdbfact(chargeCH3);

			// Deleting remaining un-necessary atoms
			for ( iter3->pos_atom = 4; !iter3->gend_atom(); iter3->next_atom() )
				res->remove(4);
			break;

			default:
			if (iter3->num_atom() < 5)
			{
			  printf( " Residue  %s%d only %d atoms\n", res->getName(),num_res, iter3->num_atom());
			  exit(1);
			}

			iter3->pos_atom = 0; //NH
			(( Atom * ) iter3->get_atom())->setPdbocc(massNH);
			(( Atom * ) iter3->get_atom())->setPdbfact(chargeNH);

			iter3->next_atom(); // CA+HA
			(( Atom * ) iter3->get_atom())->setPdbocc(massCH);
			(( Atom * ) iter3->get_atom())->setPdbfact(chargeCH);

			iter3->next_atom(); // CO
			(( Atom * ) iter3->get_atom())->setPdbocc(massCO);
			(( Atom * ) iter3->get_atom())->setPdbfact(chargeCO);

			iter3->next_atom(); // CB
			res->exchange( 3, 4 );   // exchange O for CB

			iter3->pos_atom = 4; // O in the Macromolecule struct
			if (( resn == ILE ) ||( resn == THR ) ||( resn == VAL ))
			{
			 (( Atom * ) iter3->get_atom())->setPdbocc(massCH); // Setting CH mass on O-atom ??? (pos 4)
			 (( Atom * ) iter3->get_atom())->setPdbfact(chargeCH); // Setting CH mass on O-atom ??? (pos 4)
			}
			else
			{
				(( Atom * ) iter3->get_atom())->setPdbocc(massCH2);
				(( Atom * ) iter3->get_atom())->setPdbfact(chargeCH2);
			}

			// 5th R2 --> center of mass of the rest of the chain..
			cx = 0; cy = 0; cz = 0; atoms = 0;

			masstotal=0;
			chargetotal=0;
			atom_res=5;
			atoms_R=0;
			for ( iter3->pos_atom = 5; !iter3->gend_atom(); iter3->next_atom() )
			{
			  atom = ( Atom * ) iter3->get_atom();
			  at_name = atom->getName();
			  atom->getPosition(pos);
			  atoms_R++;
			  //printf( " %d ++%s++ %f %f %f\n", iter3->pos_atom, at_name, pos[0],pos[1],pos[2]);
			  if ( strncmp( at_name, " H  ",4 ) )
			  if ( strncmp( at_name, " HN ",4 ) )
			  if ( strncmp( at_name, " HA ",4 ) )
			  if ( strncmp( at_name, " HB ",4 ) )
			  if ( strncmp( at_name, "1HB ",4 ) )
			  if ( strncmp( at_name, "2HB ",4 ) )
			  { // whether it is not those Hydrogens...
				  if(nomass)
					  mass= imass;
				  else
					  mass= ( atom->getElement() )->weight; // atomic weight
				  if(nocharge)
					  charge= icharge;
				  else
					  charge = ( atom->getElement() )->number; // atomic number (charge)
				masstotal+=mass;
				chargetotal+=charge;
				cx += pos[0]*mass;
				cy += pos[1]*mass;
				cz += pos[2]*mass;
				atoms++;
				//printf (" %s Mass %f pos  %f %f %f\n",at_name,  mass, pos[0],pos[1],pos[2]);
				//getchar();
			  }
			}
		  pos[0] = cx / masstotal;
		  pos[1] = cy / masstotal;
		  pos[2] = cz / masstotal;
		  //printf (" res %s masstotal %f atomos R %d\n",res->getName(), masstotal, atoms_R);

		  // 5th rest of the chain..
		  iter3->pos_atom = 3; // "O" ???
		  atom = ( Atom * ) iter3->get_atom();
		  atom->setPosition(pos);
		  atom->setPdbocc(masstotal);
		  atom->setPdbfact(chargetotal);
		  iter3->~pdbIter();
		  //printf( " mass %f cxyz %f %f %f atoms %d %d\n",
		  //masstotal, pos[0],pos[1],pos[2],atoms, atoms_R);

		  // MON: this is not needed!
//		  if (AA[resn].natoms!=atom_res+atoms_R)
//		  {
//			  printf(" Warning res %s have %d atoms instead of %d\n",res->getName(), atom_res+atoms_R, AA[resn].natoms);
//			  //getchar();
//		  }
		  for ( j = 0 ; j < atoms_R; j++)
			  res->remove(5);
		  break;
		 }
		}
		delete iter2;
	//    iter2->~pdbIter();
    }
  }
  delete iter1;
//  iter1->~pdbIter();

return num_res;

}

// Delete hydrogens from the macromolecule. IMPORTANT: Use this instead of delete_hydrogens()
void Macromolecule::deleteHYDS( )
{
	SMol *smol;
	int i;
	for(i=0;i<limit;i++)
	{
		if(((PDB_Container*)elements[i])->getMolType()==tmol_smol)
		{
			smol=(SMol*)elements[i];
			smol->deleteHydrogenBonds();
		}
	}

	delete_hydrogens();
}

// Delete waters from the macromolecule. IMPORTANT: Use this instead of delete_hydrogens()
void Macromolecule::delete_waters( )
{
	delete_water();
}

// Delete Hetero-atoms
void Macromolecule::delete_heteros( )
{
	delete_hetero();
}

// Delete duplicate atoms within a fragment
void Macromolecule::delete_duplicates( )
{
	Fragment *frag;
	pdbIter *iter_frags;

	iter_frags = new pdbIter(this);

	// Screen fragment-wise
	for ( iter_frags->pos_fragment = 0; !iter_frags->gend_fragment(); iter_frags->next_fragment() )
	{
		frag = ( Fragment * ) iter_frags->get_fragment();
		frag->delete_duplicates(); // Delete contained duplicate atoms preserving the first occurrence (Mon 21/02/2018).
	}

	delete iter_frags;
}

Macromolecule *Macromolecule::select_Box(Tcoor C, float xdim)
{

  bool contM, contCh, contS, contR, contA;
  bool insert;
  int cM,cCh,cS, cR, cA ;
	cM=cCh=cS= cR= cA=-1;

  Segment * seg;
  Chain * cad;
  Molecule * mol;
  Fragment * res;
  Atom * a;
  Tcoor pos;

  Segment * seg2;
  Chain * cad2;
  Molecule * mol2;
  Fragment * res2;

  float maxX, maxY, maxZ;
  float minX, minY, minZ;

  minX=C[0]-xdim;
  maxX=C[0]+xdim;
  minY=C[1]-xdim;
  maxY=C[1]+xdim;
  minZ=C[2]-xdim;
  maxZ=C[2]+xdim;


  Macromolecule * m = new Macromolecule( (char *)"selection" );

  initAll();

  contM = true;
  while ( contM != false )
  {
    mol = ( Molecule * ) this->getCurrent();

    cM++;
    switch(mol->getClass())
    {
    case pdb_protein:
    	mol2 = new Protein( (char *)"Selection" );
    	break;
    case pdb_nacid:
    	mol2 = new NAcid( (char *)"Selection" );
    	break;
    case pdb_smol:
    	mol2 = new SMol( mol->getName(),mol->getIdNumber() );
    	break;
    default :
    	break;
    }

    if(mol->getClass()==pdb_smol)
    {
    	contA = true;
    	while ( contA != false )
    	{
    		cA++;
    		a = ( Atom * ) mol->getCurrent();
    		insert = false;

    		a ->getPosition(pos);

    		if ( ((pos[0]>minX)&&(pos[0]<maxX)) &&
    			     ((pos[1]>minY)&&(pos[1]<maxY)) &&
    			     ((pos[2]>minZ)&&(pos[2]<maxZ))   ) {
    				insert=true;
    		}
    		else {
    				insert=false;
    		}
    			if ( insert == true )
    		{
    			mol2->add( a );
    		}

    		contA = mol->next();
    	}
    	if(mol2->getLimit()>0)
    		m->add(mol2);
    	else
        	delete(mol2);

    }
    else
    {
    	contCh = true;
    	while ( contCh != false )
    	{
    		cad = ( Chain * ) mol->getCurrent();

    		cCh++;
    		cad2 = new Chain( cad->getName() );

    		contS = true;
    		while ( contS != false )
    		{
    			seg = ( Segment * ) cad->getCurrent();

    			cS++;
    			seg2 = new Segment((char *)"Selection" );

    			contR = true;
    			while ( contR != false )
    			{
    				res = ( Fragment * ) seg->getCurrent();

    				//cR++;
    				cR=res->getIdNumber();
    				if(res->getClass()==pdb_residue)
    					res2 = new Residue( res->getName(), res->getIdNumber(),res->get_pos(),res->get_letter());
    				if(res->getClass()==pdb_nucleotide)
    					res2 = new Nucleotide( res->getName(), res->getIdNumber(),res->get_pos(),res->get_letter());

    				contA = true;
    				while ( contA != false )
    				{
    					cA++;
    					a = ( Atom * ) res->getCurrent();
    					insert = false;

    					a ->getPosition(pos);

    					    		if ( ((pos[0]>minX)&&(pos[0]<maxX)) &&
    					    			     ((pos[1]>minY)&&(pos[1]<maxY)) &&
    					    			     ((pos[2]>minZ)&&(pos[2]<maxZ))   ) {
    					    				insert=true;
    					    		}
    					    		else {
    					    				insert=false;
    					    		}

    					if ( insert == true )
    					{
    						// cout<<"Insertando atomo: "<<a->getPdbName()<<endl;
    						res2->add( a );
    					}

    					contA = res->next();
    				}

    				if ( res2->getLimit() > 0 )
    					seg2->add( res2 );
					else
						delete(res2);

    				contR = seg->next();
    			}

    			if ( seg2->getLimit() > 0 )
    				cad2->add( seg2 );
				else
					delete(seg2);

    			contS = cad->next();
    		}

    		if ( cad2->getLimit() > 0 )
    			mol2->add( cad2 );
    		else
    			delete(cad2);

    		contCh = mol->next();
    	}

    	if ( mol2->getLimit() > 0 )
    		m->add( mol2 );
    	else
    		delete(mol2);
    }

    contM = this->next();

  }
    return m;

}

// Select residues within the box (rectangular) defined by minimum (min) and maximum (max) corner points.
Macromolecule *Macromolecule::select_Box(float *min, float *max)
{

  bool contM, contCh, contS, contR, contA;
  bool insert;
  int cM, cCh, cS, cR, cA ;
  cM= cCh= cS= cR= cA=-1;

  Segment * seg;
  Chain * cad;
  Molecule * mol;
  Fragment * res;
  Atom * a;
  Tcoor pos;

  Segment * seg2;
  Chain * cad2;
  Molecule * mol2;
  Fragment * res2;

  float maxX, maxY, maxZ;
  float minX, minY, minZ;

  minX=min[0];
  maxX=max[0];
  minY=min[1];
  maxY=max[1];
  minZ=min[2];
  maxZ=max[2];


  Macromolecule * m = new Macromolecule((char *)"selection" );

  initAll();

  contM = true;
  while ( contM != false )
  {
    mol = ( Molecule * ) this->getCurrent();

    cM++;
    switch(mol->getClass())
    {
    case pdb_protein:
    	mol2 = new Protein((char *)"Selection" );
    	break;
    case pdb_nacid:
    	mol2 = new NAcid((char *)"Selection" );
    	break;
    case pdb_smol:
    	mol2 = new SMol( mol->getName(),mol->getIdNumber() );
    	break;
    default:
    	break;
    }

    if(mol->getClass()==pdb_smol)
    {
    	contA = true;
    	while ( contA != false )
    	{
    		cA++;
    		a = ( Atom * ) mol->getCurrent();
    		insert = false;

    		a ->getPosition(pos);

    		if ( ((pos[0]>minX)&&(pos[0]<maxX)) &&
    			     ((pos[1]>minY)&&(pos[1]<maxY)) &&
    			     ((pos[2]>minZ)&&(pos[2]<maxZ))   ) {
    				insert=true;
    		}
    		else {
    				insert=false;
    		}
    			if ( insert == true )
    		{
    			mol2->add( a );
    		}

    		contA = mol->next();
    	}
    	if(mol2->getLimit()>0)
    		m->add(mol2);
    	else
        	delete(mol2);

    }
    else
    {
    	contCh = true;
    	while ( contCh != false )
    	{
    		cad = ( Chain * ) mol->getCurrent();

    		cCh++;
    		cad2 = new Chain( cad->getName() );

    		contS = true;
    		while ( contS != false )
    		{
    			seg = ( Segment * ) cad->getCurrent();

    			cS++;
    			seg2 = new Segment((char *)"Selection" );

    			contR = true;
    			while ( contR != false )
    			{
    				res = ( Fragment * ) seg->getCurrent();

    				//cR++;
    				cR=res->getIdNumber();
    				if(res->getClass()==pdb_residue)
    					res2 = new Residue( res->getName(), res->getIdNumber(),res->get_pos(),res->get_letter());
    				if(res->getClass()==pdb_nucleotide)
    					res2 = new Nucleotide( res->getName(), res->getIdNumber(),res->get_pos(),res->get_letter());

    				contA = true;
    				while ( contA != false )
    				{
    					cA++;
    					a = ( Atom * ) res->getCurrent();
    					insert = false;

    					a ->getPosition(pos);

    					    		if ( ((pos[0]>minX)&&(pos[0]<maxX)) &&
    					    			     ((pos[1]>minY)&&(pos[1]<maxY)) &&
    					    			     ((pos[2]>minZ)&&(pos[2]<maxZ))   ) {
    					    				insert=true;
    					    		}
    					    		else {
    					    				insert=false;
    					    		}

    					if ( insert == true )
    					{
    						// cout<<"Insertando atomo: "<<a->getPdbName()<<endl;
    						res2->add( a );
    					}

    					contA = res->next();
    				}

    				if ( res2->getLimit() > 0 )
    					seg2->add( res2 );
					else
						delete(res2);

    				contR = seg->next();
    			}

    			if ( seg2->getLimit() > 0 )
    				cad2->add( seg2 );
				else
					delete(seg2);

    			contS = cad->next();
    		}

    		if ( cad2->getLimit() > 0 )
    			mol2->add( cad2 );
    		else
    			delete(cad2);

    		contCh = mol->next();
    	}

    	if ( mol2->getLimit() > 0 )
    		m->add( mol2 );
    	else
    		delete(mol2);
    }

    contM = this->next();

  }
    return m;

}


// MON made (23/5/2019)
// Write just the indicated loop (from the index of first residue "ifr" to index of last residue "lfr") into a Multi-PDB
bool Macromolecule::writeMloop( char *name, int n, int ifr, int ilr, char chain)
{
	FILE *file;
	char line[83];
	char *fmt;
	Fragment *res;
	Atom *a;

	fmt = ( char * ) "ATOM  %5d %-5s%3s %c%4d%c   %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n";

	file = fopen(name,"at");
	if(file == NULL)
		return false;

//	COLUMNS       DATA TYPE      FIELD         DEFINITION
//	----------------------------------------------------------------------
//	 1 -  6       Record name    "MODEL "
//	11 - 14       Integer        serial        Model serial number.

	sprintf( line, "MODEL   %6d\n",n);
	fputs( line, file );

	initAll();


	Tcoor pos;

//	fprintf(stderr,"About to create iterator\n");
	pdbIter *iter = new pdbIter( this, true, true, true, true ); // Iterator for residues
//	fprintf(stderr,"Iterator created\n");
//	exit(0);

	pdbIter *iter2;
	for( iter->pos_fragment = ifr; iter->pos_fragment <= ilr; iter->next_fragment() ) // screen current loop residue atoms
	{
		res = (Fragment *) iter->get_fragment(); // Get residue
		iter2 = new pdbIter( res ); // Iterator for atoms
		for( iter2->pos_atom = 0; !iter2->gend_atom(); iter2->next_atom() ) // screen current loop residue atoms
		{
			// Get atom position
			a = (Atom *) iter2->get_atom();
			a->getPosition(pos);

//			fprintf(stderr,"About to write atom\n");

//			sprintf( line, fmt, a->getPdbSerial(), a->getPdbName(), res->getName(),
//					cad->getName() [0],res->getIdNumber()/*res->get_pos()*/,res->get_letter(), pos[0], pos[1], pos[2], a->getPdbocc(), a->getPdbfact() );
			sprintf( line, fmt, a->getPdbSerial(), a->getPdbName(), res->getName(),
//					chain, res->getIdNumber()/*res->get_pos()*/,res->get_letter(), pos[0], pos[1], pos[2], a->getPdbocc(), a->getPdbfact() );
	        // remove  Occupancy annd Bfactor
					chain, res->getIdNumber()/*res->get_pos()*/,res->get_letter(), pos[0], pos[1], pos[2], 0.0, 0.0 );

			fputs( line, file );

//			fprintf(stderr,"Atom wrote\n");
		}
		delete iter2;
	}
	delete iter;

	sprintf( line, "TER    %4d      %s %c%4d                                                      \n",
			a->getPdbSerial() + 1, res->getName(),
			chain, res->getIdNumber()); // MON: substitute "0" by "cad->getName()[0]" when Chain ids are implemented
	fputs( line, file );

	sprintf( line, "ENDMDL\n");
	fputs( line, file );

	fclose(file);
	return true;
}


