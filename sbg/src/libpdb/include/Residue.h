#ifndef RESIDUE_H
  #define RESIDUE_H

  #include "Fragment.h"

typedef float( * COORDS ) [3];



/**
 *An aminoacid inside a protein. Stores a set of atoms
 *For more information about functionability see classes Contained and Container
 */
class Residue : public Fragment
{
public:

	///Constructor
	  ///@param name: name of the object created
	  Residue( Tname name, int in_nid=0, int in_npos=0, char i_letter=' ') ;

	  ///Copy Constructor
	  ///@param old: Reference object to make the copy
	  ///@param with_elements: If true , it copies the array of elements in the buffer of the new object, otherwise the buffer is empty
	  Residue( Residue * old ,bool with_elements=true);


  /** Constructor by template. Creates a residue using the list AA of standard
    * aminoacid.
    *THE AA LIST MUST BE INITIALIZED
    *@param aa: Position of the reference aminoacid in the AA list
    *@param in_nid: Residue Identifier
    *@param in_npos: Position of residue in Macromolecule
    *@param in_letter: Associated letter of residue in PDB file
    */
  Residue( int aa , int in_nid=0, int in_npos=0, char in_letter=' ',TMOL t=tmol_null);

  /// Gets all angle chis in an aminoacid
  ///THE AA LIST MUST BE INITIALIZED
  ///@param aan: Reference aminoacid in the AA list
  ///@param chis: List of angles. It must be allocated (4 floats)
  void get_achis( int aan, double * chis );



  /**
   * Calculates the coordinates of the atoms that required
   * chi angles to build it, given new angles. Modifies the
   * positions of the atoms in the residue
   *
   *@param: newchi: List of new chi angles
   */
  void rotamerize2( float * newchi );

  /**Aligns the NH, CA and C atoms of two residues. Returns the new
   * coordinates of the atoms
   *
   * @param resi: Second residue
   * @param positions: New coordinates of the residue. It must be allocated first
   */
  void Align( Residue * resi, COORDS positions );

  /**Aligns the residue backbone over the given coordinates. The atoms coordinates of the residue are modified
   *
   * @param NH: Coordinates of the NH atom
   * @param CA: Coordinates of the CA atom
   * @param C: Coordinates of the C atom
   * @param atom: Determines the vector used to align in the first step. Each vector produces different alignments
   */
  void Align_backbone( float NH[3], float CA[3], float C[3], char atom[2] );

  /**Aligns the residue backbone over the given coordinates. The atoms coordinates of the residue are modified
   * It also returns the rotational matrices and the translational vectors used in the two alignment steps
   * @param NH: Coordinates of the NH atom
   * @param CA: Coordinates of the CA atom
   * @param C: Coordinates of the C atom
   * @param atom: Determines the vector used to align in the first step. Each vector produces different alignments
   * @param matrix_out: Rotational matrix in first alignment step
   * @param vector_out: Translational Vector in the alignment first step
   * @param matrix_out2: Rotational matrix in second alignment step
   * @param vector_out2: Translational Vector in the second alignment step
   *
   */
  void Align_backbone( float NH[3], float CA[3], float C[3], char atom[2], float ***matrix_out, float **vector_out,float ***matrix_out2, float **vector_out2 );

 /**
  * Rotation of the residue
  *
  * @param matrix: Rotational matrix
  */
  void Rotate( float **matrix );

 // Rotation and Translation of the residue
 // @param matrix: Rotational matrix
 // @param vector: Traslational vector
  void Rotrans( float **matrix, float *vector );

  // Moves the residue so CA is in Origin point(0.0, 0.0, 0.0)
  void OriginCA();

  // Computes movement of the residue so CA is in Origin point(0.0, 0.0, 0.0)
  //@param positions: Traslation required
  void OriginCA( COORDS positions );

   // Computes movement in a position list so CA is in Origin point(0.0, 0.0, 0.0)
   // It does not affect to the residue
   //@param positions: Traslation required
     void OriginCA2( COORDS positions );

  // Moves the residue so CA is in Origin point(0.0, 0.0, 0.0)
  void OriginC();

   /**
   * Replaces an atom in the space in function of other atoms of the residue
   *   There are variables that tell how to build atoms that are missing ta.
   * A complete coordinate set of each amino acid is stored (icoor).
   * For each atom that must be placed, three 'template' atoms are defined.
   * The complete icoor set is
   * translated and rotated such that the template atoms are best fit.
   * the coordinates for the missing atom can be taken from the
   * processed icoor set.
   *
   *
   * @param atom: Atom to move
   */
  void place_atom( int atom );

  ///Returns the class of the object.
  TElement getClass();

  /// Returns the Molecule type that contains the object
  TMOL getMolType();

};



void lineup_bk( float * a1, float * a2, float * b1, float * b2, float * * mat, float * vec );

void align_bk( float * a1, float * a2, float * a3, float * b1, float * * mat, float * vec );

#endif


