/*
 * Nucleotide.h
 *
 *  Created on: Feb 5, 2010
 *      Author: nacho
 */

#ifndef NUCLEOTIDE_H_
#define NUCLEOTIDE_H_

 #include "Fragment.h"





/**
 *An aminoacid inside a protein. Stores a set of atoms
 *For more information about functionability see classes Contained and Container
 */
class Nucleotide : public Fragment
{
public:

	///Constructor
	///@param name: name of the object created
	Nucleotide( Tname name, int in_nid=0, int in_npos=0, char i_letter=' ' ) ;

	///Copy Constructor
	///@param old: Reference object to make the copy
	///@param with_elements: If true , it copies the array of elements in the buffer of the new object, otherwise the buffer is empty
	Nucleotide( Nucleotide * old ,bool with_elements=true);



  //Possibly neccessary
  //void get_achis( int aan, float * chis );
  //void rotamerize2( float * newchi );
  //void Align( Nucleotide * resi, COORDS positions );
  //void Align_backbone( float NH[3], float CA[3], float C[3], char atom[2] );
  //void Align_backbone( float NH[3], float CA[3], float C[3], char atom[2], float ***matrix_out, float **vector_out,float ***matrix_out2, float **vector_out2 );
  //void Rotate( float **matrix );
  //void Rotrans( float **matrix, float *vector );
  //void OriginCA();
  //void OriginCA( COORDS positions );
  //void OriginCA2( COORDS positions );
  //void OriginC();
  //void place_atom( int atom );

  ///Returns the class of the object.
  TElement getClass();

  /// Returns the Molecule type that contains the object
  TMOL getMolType();

};

#endif /* NUCLEOTIDE_H_ */
