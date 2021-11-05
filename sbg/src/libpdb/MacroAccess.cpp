#include "Macromolecule.h"

Atom * Macromolecule::getCurrentAtom()
{
  Atom * a;
  Fragment * frag;
  Segment * seg;
  Chain * cad;
  Molecule * p;

  p = ( Molecule * ) this->getCurrent();
  if(p->getClass()==pdb_smol)
  {
	  a=(Atom*)p->getCurrent();

  }
  else
  {
	  cad = ( Chain * ) p->getCurrent();
	  seg = ( Segment * ) cad->getCurrent();
	  frag = ( Fragment * ) seg->getCurrent();
	  a = ( Atom * ) frag->getCurrent();
  }
  return a;
}

bool Macromolecule::nextAtom()
{


  Segment *seg;
  Fragment * frag;
  Chain * cad;
  Molecule * p;
  bool more;

  p = ( Molecule * ) this->getCurrent();
  if(p->getClass()==pdb_smol)
  {
	 more= p->next();

	 if(more)
	 {
		 return true;
	 }
	 else
	 {
		 /* Se avanza a la siguiente molecula */
		 more = this->next();
		 if ( more )
		{
			 return ( true );
		}
		/* Se ha llegado al final del sistema, no quedan mas atomos por explorar */
		 else
		{
			 return ( false );
		}
	 }
  }
  else
  {
	  cad = ( Chain * ) p->getCurrent();
	  seg = ( Segment * ) cad->getCurrent();
	  frag = ( Fragment * ) seg->getCurrent();

	  /* Se avanza en el residuo */
	  more = frag->next();

	  if ( more  )
		  return ( true );

	  /* Se ha llegado al fin del residuo */
	  else
	  {

		  /* Se avanza al siguiente residuo de el segmento */
		  more = seg->next();
		  if ( more )
			  return ( true );

		  /* Se ha llegado al final del segmento */
		  else
		  {
			  /* Se avanza al siguiente segmento de la cadena */
			  more = cad->next();

			  if ( more )
				  return ( true );

			  /* Se ha llegado al final de la cadena */
			  else
			  {

				  /* Se avanza a la siguiente cadena de la molecula */
				  more = p->next();

				  if ( more  )
					  return ( true );

				  /* Se ha llegado al final de la molecula */
				  else
				  {

					  /* Se avanza a la siguiente molecula */
					  more = next();

					  if ( more )
						  return ( true );
					  /* Se hallegado al final del sistema, no quedan mas atomos por explorar */
					  else
						  return ( false );
				  }

			  }
		  }
	  }
  }
}
