
#ifndef MACROMOLECULE_H
#define MACROMOLECULE_H




#include <math.h>
#include <ctype.h>
#include "Protein.h"
#include "SMol.h"
#include "NAcid.h"
#include "ResIni.h"
#include "atomoperator.h"
#include "container.h"
#include "Condition.h"



#include <vector>
using namespace std;



#ifdef VOLUME_INCLUDED

#include "libvolume/include/vlvolume.h"
#include "libtools/include/Memstuff.h"



#ifdef USE_PTHREAD // Enables PThread parallel routines
typedef struct // required by project_RealNE_3BB2R_par() (parallel)
{
	//	Macromolecule *mol;
	//	pdbIter *iter;
	//	void *mol; // it will be casted into Macromolecule when needed (mandatory to prevent some "dichotomy")
	void *iter; // it will be casted into pdbIter when needed (mandatory to prevent some "dichotomy")
	vlVolume *vlm;
	float voxel;
	int stepy;
	int stepz;
	int limitx;
	int limity;
	int limitz;
	float unitsx;
	float unitsy;
	float unitsz;
	int first;
	int last;
	pthread_mutex_t *p_mutex_map; // mutex pointer (to protect result map access)
	pthread_mutex_t *p_mutex_begin; // mutex pointer
	pthread_cond_t *p_cond_begin; // condition pointer
	pthread_mutex_t *p_mutex_end; // mutex pointer
	pthread_cond_t *p_cond_end; // condition pointer
	int *p_nended; // pointer to the number of "ended" jobs
	bool begin; // bool telling whether the job started or not
} project_RealNE_data;

//#ifdef USE_PTHREAD // Enables PThread parallel routines

// Mon made (29/05/2013)
// Projects over an pre-existing vlVolume a Macromolecule with NE and without moving
// the reference system in parallel. (Thread routine)
// If fast=true, a fast method will be used instead of trilinear-interpolation (default=false)
void *project_RealNE_3BB2R_thread( void *threadarg );
#endif

#ifndef _floatops_hpp
namespace FOPS
{
extern bool writeFile(vlVolume *vol,char *filename,bool ccp4Endian=false);
extern vlVolume * padVolume(vlVolume *old, vlDim margin, bool self=false);
extern vlVolume * padVolume(vlVolume *old, vlDim margin,int dummy,bool self = false);

}
#endif
#endif

///Auxiliar function
int pdb_strcasecmp(const char *s1, const char *s2);

/**
 * Class that stores all the information about a molecular system,
 * usually read from a PDB file. In a system there could be one or more molecules
 * (proteins, etc.)
 * This class is a specialization of Class Container
 *
 * Basic structure:
 * It contains a array of Molecule objects (Proteins, etc.)
 * It also contains a array of Heteromolecule (Fragment specialization)
 * It provides direct access functions to the atoms
 * For more information about functionability see classes Contained and PDB_Container
 */
class Macromolecule: public PDB_Container {

private:
	/**Object that contains the spatial dimensions of the system.
	 * The center stored is computed as the middle point between the farthest atoms of the system (included hydrogens)
	 *
	 */
	Dimensions macromoleculeDimension;




public:

	//CONSTRUCTION AND DESTRUCTION METHODS
	/// Simple Constructor
	Macromolecule();
	/// Simple Constructor
	///@param name: Name of the macromolecule
	///@param i_nid: ID number of the macromolecule
	Macromolecule(Tname name,int i_nid=0);
	/// Constructor. Make a copy of the referenced object
	/// @param old: reference object to be copied
	Macromolecule(Macromolecule *old);
	/// Constructor. Combines diferent system in only one object.
	/// The lists of Molecules and heteromolecules are combined
	/// @param list: List of molecules of the macromolecule
	/// @param num_m: number of molecules in the list
	Macromolecule( Macromolecule **list, int num_m );
	/// Destructor
	virtual ~Macromolecule();
	///Returns the class of the object
	TElement getClass();
	/// Returns the Molecule type that contains the object
	TMOL getMolType();


	//ATOM ACCESS
	/// Return the current pointed atom
	Atom* getCurrentAtom();
	///Point to the atom next to the current one
	bool nextAtom();


	//MODIFICATION OF THE SYSTEM
	/** Move all to center (using the macromoleculeDimension structure information) the system in the origin of coordinates.
	 *  It employs the geoBox function, therefore The center stored is computed as the middle point between the farthest atoms of the system (included hydrogens)
	 */
	void centerBox();
	/// Apply an operation to all the atoms of the system
	/// @param op: Operator to be applied
	bool applyAtoms(AtomOperator *op);
	/**
	 *Apply an operation to all the atoms of the system
	 * The information of the atoms is taen from a reference system that is not
	 * modified. The reference system must be equal to the system to modify
	 * @param op: Operator to be applied
	 * @param ref: Reference Macromolecule
	 */
	bool applyAtoms(AtomOperator *op, Macromolecule *ref);
	/**
	 * Move the macromolecule to the position of a similar macromolecule.
	 * The reference system must be equal to the system to modify
	 * @param ref: Reference Macromolecule
	 */
	bool copy_coordinates(Macromolecule *ref);
	/**
	 * Sets atomic positions from the provided array.
	 * @param coords: Coordinates array
	 */
	bool copy_coordinates(double *coords);
	/**
	 * Sort the atoms in a standard way inside the residues
	 * FUNCTION PROBABLY OBSOLETE
	 *
	 * @param DETAIL: Atoms to sort. 0: All atoms 1: Backbone 2: Backbone+CB 3: All but hydrogens (Default)
	 *
	 */

	void format_residues(int DETAIL);
	/**
	 * Sort the atoms in a standard way inside the residues
	 * And checks whether they are present or not according to a given Coarse-Graining model
	 *
	 * @param withH: Sort the hydrogens (true). Do not sort the hydrogens (False. Default)
	 * @param model: If indicated (>=0), then atom presence checking will be enabled  (Default= -1)
	 *
	 */
	int format_residues(bool withH=false, int model=-1);

	/**
	 * Moves the macromolecule structure
	 *
	 * @param offset: Movement to apply
	 */
	bool moveAll(Tcoor offset);
	/**
	 * Searches and renames modified residues
	 */
	void rename_residues();

	//  some sequence "seq" (in 1-letter format) from residue "index" (internal index) up to end of "seq"
	void mutseq(char seq, int index);

	/// Returns an object that stores the spatial dimensions of the system
	/// The center stored is computed as the middle point between the farthest atoms of the system (included hydrogens)
	Dimensions* getDimension();

	//I/O OPERATIONS
	/// Reads a Molecule File to fill the system. Depending on the extension in the filename it calls one of the read funcion
	///@param name: Name of the file to be read
	bool readPDB(char*name);

	/// Reads a PDB File to fill the system
	///@param name: Name of the file to be read
	bool readPDB_pdb(char*name);
	bool readPDB_old(char*name);
	/// Reads a SDF File to fill the system
	///@param name: Name of the file to be read
	bool readSDF( char *name);
	/// Reads a MOL2 File to fill the system
	///@param name: Name of the file to be read
	bool readMol2(char *name);
	/// Writes a PDB file with the information of the system
	/// @param name: name of the file to be written
	/// @param number: if true it writes the read residue number, if it is false it writes the ordinal residue number
	bool writePDB( char *name, bool number=true, bool change_name=false);
	/// Writes a PDB file with only the marked aminoacids/fragments of the system
	/// @param name: name of the file to be written
	/// @param list: list of marked aminoacids (by position)
	/// @param total: Number of marked aminoacids
	/// @param number: if true it writes the read residue number, if it is false it writes the ordinal residue number
	bool writePDB( char *name, int *list, int total,bool number=true);
	/// Writes a PDB file with the marked/nonmarked atoms of the system
	/// @param name: name of the file to be written
	/// @param list: list of marked atoms (by position)
	/// @param total: Number of marked atoms
	/// @param sign: if true, marked atoms are written, elsewhere non marked atoms are written
	/// @param number: if true it writes the read residue number, if it is false it writes the ordinal residue number
	bool writePDB_atoms( char *name, int *list, int total, bool sign=true, bool number=true);
	/// Writes a PDB file with the markeds atoms of two systems
	/// @param mol2: the other macromolecule
	/// @param name: name of the file to be written
	/// @param list: lists of marked atoms (If null, all the residues are written)
	/// @param total: Number of valid atoms in lists
	/// @param number: if true it writes the read residue number, if it is false it writes the ordinal residue number
	bool writePDB(Macromolecule *mol2, char *name, int **list=NULL, int total=0, bool number=true);
	/// Writes at the end of a PDB file the information of the system as a mode
	/// @param name: name of the file to be written
	/// @param n: Number of the mode to write
	/// @param number: if true it writes the read residue number, if it is false it writes the ordinal residue number
	bool writeMPDB( char *name, int n, bool number=true);

	void writeSDF( char *name);
	void writeMol2(char *name);
	// Mon modified (5/4/2010)
	// Convert to 3BB + 2R reduced model (Erases some atoms and changes values in others )
	// If nomass = true, all masses set to 1.0 (default=false)
	// If nocharge = true, all charges set to 1.0 (default=false)
	int  reducedmodel_3BBR2( bool nomass=false, bool nocharge=false );
	// Delete hydrogens from the macromolecule. IMPORTANT: Use this instead of delete_hydrogens()
	void  deleteHYDS( );
	// Delete waters from the macromolecule
	void  delete_waters( );
	// Delete heteros from the macromolecule
	void  delete_heteros( );
	// Delete duplicate atoms within a fragment
	void delete_duplicates( );

	// Select inside the cube box
	Macromolecule *select_Box(Tcoor C, float xdim);
	// Select residues within the box (rectangular) defined by minimum (min) and maximum (max) corner points.
	Macromolecule *select_Box(float *min, float *max);

	// Write just the indicated loop (from the index of first residue "ifr" to index of last residue "lfr") into a Multi-PDB
	bool writeMloop( char *name, int n, int ifr, int ilr, char chain = 'A');

	void lsms( float PR, bool inner=false );

	//VOLUME OPERATIONS
	// Compilation only if the Volume library is accessible
#ifdef VOLUME_INCLUDED

	/**
	 * Creations of a volume that could contain the system using voxels of the size determined by the argument
	 * This volume is empty and is not moved over the macromolecule. It is only guaranteed that it can contain the full structure if the average center is instanced over the center of the volume.
	 * The average center is the average position of all atoms (included hydrogens)
	 *	@param unit: Voxel size of the volume to be created
	 */
	vlVolume *createVolume(float unit);
	/**
	 * Creations of a volume that could contain the system using voxels of the size determined by the argument
	 * This volume is empty and is moved over the macromolecule. The center of the volume is set over the average center of the Macromolecule
	 * The average center is the average position of all atoms (included hydrogens)
	 *	@param unit: Voxel size of the volume to be created
	 *	@param Hyd: Use (true, default) or not use (false) the hydrogens to compute the average center
	 */
	vlVolume *createVolumeCentered(float unit, bool Hyd=true);
	/**
	 * Creates and fills a volume that contains a projection of the system (Number of electrons)
	 * The center of the volume is set over the average center of the macromolecule
	 * The average center is the average position of all atoms (included hydrogens)
	 * @param unit: Voxel size of the volume to be created
	 * FITTED FORM-FACTORS
	 */
	vlVolume *fillVolumeBFS(float unit, float *bfs);
	/**
	 * Creates and fills a volume that contains a projection of the system (Number of electrons)
	 * The center of the volume is set over the average center of the macromolecule
	 * The average center is the average position of all atoms (included hydrogens)
	 * @param unit: Voxel size of the volume to be created
	 */
	vlVolume *fillVolumeNE(float unit);
	/**
	 * Creates and fills a volume that contains a projection of the system (Number of electrons.3BB2R reduced PDB model. In bfact column)
	 * The center of the volume is set over the average center of the macromolecule
	 * The average center is the average position of all atoms (included hydrogens)
	 * @param unit: Voxel size of the volume to be created
	 */
	vlVolume *fillVolumeNE_3BB2R(float unit);
	/**
	 * Creates and fills a volume that contains a projection of the system (Mass weight)
	 * The center of the volume is set over the average center of the macromolecule
	 * The average center is the average position of all atoms (included hydrogens)
	 * @param unit: Voxel size of the volume to be created
	 */
	vlVolume *fillVolumeMW(float unit);
	/**
	 * Creates and fills a volume that contains a projection of the system (size of the volume (Number of electrons)
	 * The center of the volume is set over the average center of the macromolecule
	 * The average center is the average position of all atoms (included hydrogens)
	 * @param maxLenght: size of the volume in voxels (maxLenghtXmaxLenghtXmaxLenght)
	 * @param unit: Voxel size of the volume to be created
	 */
	vlVolume *fillVolume(float maxLenght,float unit);
	/**
	 * Creates and fills a volume that contains a projection of the system (Beta. fact)
	 * The center of the volume is set over the average center of the macromolecule
	 * The average center is the average position of all atoms (included hydrogens)
	 * @param unit: Voxel size of the volume to be created
	 */
	vlVolume *fillVolumeBeta( float unit );

	/**
	 * Projects in a Volume the system (number of electrons). The Volume is moved to be centered over the system
	 * The center of the volume is set over the average center of the macromolecule
	 * The average center is the average position of all atoms (included hydrogens)
	 *  @param vol: Volume to be filled
	 */
	void project(vlVolume *vol);
	/**
	 *  Projects directly in a Volume the system (number of electrons). The Volume is not moved
	 *  @param vol: Volume to be filled
	 */
	void project_Real( vlVolume * vlm );
	/**
	 *  Projects directly in a Volume the system (1.0 per atom). The Volume is not moved
	 *  @param vol: Volume to be filled
	 */
	void project_unit( vlVolume * vlm );

	/**
	 *  Projects directly in a Volume the system (Number of electrons.3BB2R reduced PDB model. In bfact column). The Volume is not moved.
	 *  @param vol: Volume to be filled
	 *  @param fast: If fast=true, a fast method will be used instead of trilinear-interpolation (default=false)
	 */
	void project_RealNE_3BB2R( vlVolume * vlm, bool fast=false );

#ifdef USE_PTHREAD // Enables PThread parallel routines
	// Mon made (29/05/2013)
	// Projects over an pre-existing vlVolume a Macromolecule with NE and without moving
	// the reference system in parallel. (Initialization routine)
	// If fast=true, a fast method will be used instead of trilinear-interpolation (default=false)
	void project_RealNE_3BB2R_par_init(int nthreads, project_RealNE_data **p_threads_data, pthread_t **p_threads);

	// Mon made (29/05/2013)
	// Projects over an pre-existing vlVolume a Macromolecule with NE and without moving
	// the reference system in parallel. (Controller routine)
	// If fast=true, a fast method will be used instead of trilinear-interpolation (default=false)
	void project_RealNE_3BB2R_par( vlVolume *vlm, int nthreads, project_RealNE_data *threads_data );
#endif

	/**
	 * Creations of a volume that could contain the system using voxels of the size determined by the argument.
	 * This volume is empty and is not moved over the macromolecule. It is only guaranteed that it can contain the full structure if the average center is instanced over the center of the volume.
	 * Size of volume is increased
	 * The average center is the average position of all atoms (included hydrogens)
	 *	@param unit: Voxel size of the volume to be created
	 *	@param pad: Volume padding (not duplicate)
	 */
	vlVolume *createVolume_pad(float unit, int pad);
	/**
	 * Creates and fills a volume that contains a projection of the system (Number of electrons)
	 * The projection is performed with a gaussian filter
	 * The center of the volume is set over the average center of the macromolecule
	 * The average center is the average position of all atoms (included hydrogens)
	 * @param unit: Voxel size of the volume to be created
	 */
	vlVolume *pdb2map_real(float res,float unit);
	/**
	 * Creates and fills a volume that contains a projection of the system (Number of electrons.3BB2R reduced PDB model. In bfact column). The Volume is not moved.
	 * The projection is performed with a gaussian filter
	 * The center of the volume is set over the average center of the macromolecule
	 * The average center is the average position of all atoms (included hydrogens)
	 * @param unit: Voxel size of the volume to be created
	 */
	vlVolume *pdb2map_real_3BB2R(float res,float unit);
	/**
	 * Creates a volume (centered over the average center of the Macromolecule)
	 * The volume contains a projection of the vdw potential
	 *
	 *@param filename: Name of the file to be created
	 *@param voxel_size: voxel size of the volume to be created
	 *@param opt: indicates what Van Der Waals potential is used (Carbon, Hydrogen, large)
	 *@param opt2: indicates the name format of the atoms (pdb,iupac)
	 *@param opt3: Indicates the force field convention (Rosseta, ICM)
	 *@param probeRad: Radius of the probe atom used to compute the potential
	 *@param probeEmax: Inflexion point in vdw function (recommended: 1.0)
	 *@param withH: Computes the vdw potential using hydrogens too.
	 *@param center_hyd: Computes the center of volume using hydrogens too.
	 */
	vlVolume * createVDW(float voxel_size, int opt, ConventionNames opt2,Convention opt3,float probeRad,float probeEmax ,int withH,bool center_hyd=false,float rm=1.0,float f1=1.0,float f2=1.0);
	vlVolume * createVDW(int start, int end, float voxel_size, int opt, ConventionNames opt2,Convention opt3,float probeRad,float probeEmax ,int withH,bool center_hyd=false,float rm=1.0,float f1=1.0,float f2=1.0);
	/**
	 * Creates a volume file containing a projection of the Van der Waals potential of the Macromolecule
	 * It uses the function createVDW
	 *@param filename: Name of the file to be created
	 *@param voxel_size: voxel size of the volume to be created
	 *@param opt: indicates what Van Der Waals potential is used (Carbon, Hydrogen, large)
	 *@param opt2: indicates the name format of the atoms (pdb,iupac)
	 *@param opt3: Indicates the force field convention (Rosseta, ICM)
	 *@param probeRad: Radius of the probe atom used to compute the potential
	 *@param probeEmax: Inflexion point in vdw function (recommended: 1.0)
	 *@param withH: Computes the vdw potential using hydrogens too.
	 *@param center_hyd: Computes the center of volume using hydrogens too.
	 */
	void writeVDW(char *filename, float voxel_size, int opt, ConventionNames opt2,Convention opt3,float probeRad,float probeEmax,int withH=1,bool center_hyd=false,float rm=1.0,float f1=1.0,float f2=1.0);
	/**
	 * Projects in a volume a mask of the Macromolecule where each atom occupies an sphere with Van der Waals radius
	 *@param voxel_size: voxel size of the volume to be created
	 *@param opt: indicates what Van Der Waals potential is used (Carbon, Hydrogen, large)
	 *@param opt2: indicates the name format of the atoms (pdb,iupac)
	 */
	vlVolume *project_radiusVDW(float voxel_size, int withH,bool center_hyd,float extra_radius=0.0);
	/**
	 * Projects in a volume a hidrophoby contact potential
	 *@param voxel_size: voxel size of the volume to be created
	 *@param opt: indicates what Van Der Waals potential is used (Carbon, Hydrogen, large)
	 *@param opt2: indicates the name format of the atoms (pdb,iupac)
	 *@param ideal,max,sigma: Defines the ChemScore block funtion
	 */
	vlVolume* project_radiusVDW_loop(float voxel_size, int withH,bool center_hyd,float extra_radius=0.0);
	/**
	 * Projects in a volume a hidrophoby contact potential
	 *@param voxel_size: voxel size of the volume to be created
	 *@param opt: indicates what Van Der Waals potential is used (Carbon, Hydrogen, large)
	 *@param opt2: indicates the name format of the atoms (pdb,iupac)
	 *@param ideal,max,sigma: Defines the ChemScore block funtion
	 */
	vlVolume * project_Hydrophoby(float voxel_size,  ConventionNames opt2,Convention opt3,float ideal, float max,float sigma,int withH,bool center_hyd);
	/**
	 * Projects the ASA of a Macromolecule in a volume previously created
	 * ASA should be stored in occupancy value of atoms (Use ASA::localSAS)
	 *@param vol: Volume where the ASA is projected
	 *@param withH: Project hydrogens' ASA?
	 */
	void projectASA(vlVolume *vol,int withH);
	/**
	 * Writes in a file the ASA of a Macromolecule
	 * the volume is centered over the average center of the Macromolecule
	 * It calls the projectASA function
	 * ASA should be stored in occupancy value of atoms (Use ASA::localSAS)
	 *@param vol: Volume where the ASA is projected
	 *@param voxel_size: voxel size of the volume
	 *@param withH: Project hydrogens' ASA?
	 */
	void writeASA_c(char *filename, float voxel_size,int withH);
	/**
	 * Projects the charge of the atoms of a Macromolecule in a volume previously created
	 * Charge should be stored in occupancy value of atoms
	 *
	 *@param vol: Volume where the charge is projected
	 *@param  withH: Project hydrogens' charge?
	 */
	void projectCharge(vlVolume *vol,int withH);
	/**
	 * Writes in a file the charge of a Macromolecule
	 * the volume is centered over the average center of the Macromolecule
	 * It calls the projectASA function
	 * Charge should be stored in occupancy value of atoms
	 *@param vol: Volume where the ASA is projected
	 *@param voxel_size: voxel size of the volume
	 *@param withH: Project hydrogens' ASA?
	 */
	void writeCharge(char *filename, float voxel_size,int withH);


	/**
	 * Creates a volume (centered over the average center of the Macromolecule)
	 * The volume contains a projection of the Electrostatic potential
	 * Charge should be stored in occupancy value of atoms
	 *@param voxel_size: voxel size of the volume to be created
	 *@param threshold: Limitation in the potential value (It should be positive and will limit both positive and negative values)
	 *@param type: Type of interaction that is studied
	 *@param withH: Use hydrogen in the potential computation?
	 */
	vlVolume *createELE(float voxel_size, float threshold, int type, int withH);
	vlVolume *createELE(int start, int end, float voxel_size, float threshold, int type, int withH);

	/**
	 * Idem to createFile but takes charges from PDB file
	 */
	vlVolume *createELE_file(float voxel_size, int type, int withH);


	/**
	 * Creates a file with a volume  (centered over the average center of the Macromolecule)
	 * The volume contains a projection of the Electrostatic potential
	 * Charge should be stored in occupancy value of atoms
	 * It calls writeCharge function
	 *@param filename: Name of the file to be created
	 *@param voxel_size: voxel size of the volume to be created
	 *@param opt2: indicates the name format of the atoms (pdb,iupac)
	 *@param opt3: Indicates the force field convention (Rosseta, ICM)
	 *@param threshold: Limitation in the potential value (It should be positive and will limit both positive and negative values)
	 *@param type: Type of interaction that is studied
	 *@param withH: Use hydrogen in the potential computation?
	 */
	void writeELE(char *filename, float voxel_size, ConventionNames opt2,Convention opt3, float threshold,int type=2,int withH=0);
	void writeELE(char *filename, float voxel_size, int type=2, bool withH=0);

#endif


	//INFORMATION ABOUT THE SYSTEM
	/** Computation of the system geometric center (average positions of all atoms)
	 *@param position: position of the center
	 */
	bool geoCenter(Tcoor position);
	/**
	 * Computation of the system geometric center (ponderated avegarage positions of all atoms)
	 *@param position: position of the center
	 *@param propierty: W (weight) A (Atom number)
	 */
	bool geoProperty(float *position, char propierty);
	/// Computation of the dimensions and the center of the system (center is computed as the middle point between the farthest atoms of the Macromolecule)
	void geoBox();
	///Total electronic density of the system
	float eDensity();
	///Maximun length between a atom of the system to the Origin of coordinates.
	float maxLength();
	///Maximun length between a backbone atom of the system to the Origin of coordinates.
	float maxLength_backbone();
	///Maximum length between two atoms of the system
	float maxDist();
	/**
	 *  Write a Text file with information about the system
	 * @param name: Name of the file to be created
	 */
	void info(FILE *file);
	/**
	 * rmsd between to molecules with same number of atoms
	 * @param mol2: Reference to the second molecule
	 */
	float rmsd(Macromolecule *mol2);
	/*
	 *rmsd between selected residues in two molecules
	 *@param mol2: Reference to the other molecule
	 *@param list: double list width the associated residues (The list must store the positions of associated residues)
	 *@param total: number of common residues
	 */
	float rmsd( Macromolecule * mol2, int **list, int total );
	/**
	 * Combined rmsd between selected residues in two pairs of molecules molecules
	 * @param mol2: Reference to the other molecule (this and mol2 are the first pair)
	 * @param list: double list width the associated residues in first pair(The list must store the positions of associated residues)
	 * @param total: number of common residues in first pair
	 * @param other1: First Macromolecule in second pair
	 * @param other2: Second Macromolecule in second pair
	 * @param list2: double list width the associated residues in second pair(The list must store the positions of associated residues)
	 * @param total2: number of common residues in second pair
	 */
	float rmsd( Macromolecule * mol2, int **list, int total, Macromolecule *other1, Macromolecule *other2,int **list2, int total2);
	/**
	 * rmsd between selected residues in two molecules
	 * @param mol2: Reference to the other molecule
	 * @param list: double list width the associated residues (Each list must store characters indicating valid residues- '0' is not valid. '1' is valid)
	 */
	float rmsd( Macromolecule * mol2, char **list );
	/** Mon made (17/7/2012)
	 * RMSD between to molecules with different number of residues
	 * @param mol2: Reference to the second molecule
	 * @param mask1: Boolean mask containing matching atoms information for molecule 1
	 * @param mask2: Boolean mask containing matching atoms information for molecule 2
	 */
	float rmsd( Macromolecule *mol2, bool *mask1, bool *mask2 );
	/** Mon made (8/5/2008)
	 * Rmsd between two Macromolecules. with same number of atoms The individual squared deviation of each pair is stored in the Bfactor columm of the first (active) Macromolecule
	 *@param mol2: Reference to the other molecule
	 */
	/* the same storing local rmsd in a FIle */
	float rmsd_file( Macromolecule *mol2, bool *mask1, bool *mask2, char *name );

	float rmsd_bf( Macromolecule * mol2 );

	// MON: Weighted RMSD computation
	float wrmsd2(Macromolecule *mol2, double *weights);

	// MON: Computes the Gaussian weights profile from inter-atomic distances between "mol" and "mol2"
	// allocating profile memory if (*p_profile==NULL).
	// Weights are computed according to ec.(7) from: w= exp(-(d^2)/c ) --> c=5A (for big motions)
	// Damm & Carlson. "Gaussian-Weighted RMSD Superposition of Proteins: A Structural
	// Comparison for Flexible Proteins and Predicted Protein Structures".
	// Biophysical Journal. Volume 90, Issue 12, 15 June 2006, Pages 4558-4573
	void gaussian_weight(Macromolecule *mol2, double **p_profile, double c = 5.0);

	// MON: Computes the Gaussian weights profile from inter-atomic distances between "mol" and "mol2"
	// allocating profile memory if (*p_profile==NULL).
	// (17/7/2012) --> Taking into account sequence alignment (at atomic level) for different size molecules.
	//                 *The Weights profile its relative to the first molecule (mol)
	// Weights are computed according to ec.(7) from: w= exp(-(d^2)/c ) --> c=5A (for big motions)
	// Damm & Carlson. "Gaussian-Weighted RMSD Superposition of Proteins: A Structural
	//                   Comparison for Flexible Proteins and Predicted Protein Structures".
	//                   Biophysical Journal. Volume 90, Issue 12, 15 June 2006, Pages 4558-4573
	void gaussian_weight(Macromolecule *mol2, double **p_profile, bool *mask1, bool *mask2, double c = 5.0);

	// MON: Converts an alingment mask at residue level (made with read_aln) into one at atomic level (CA)
	bool *maskres2maskatom(bool *mask);

	/**
	 * Compute minimum rmsd between two molecules with same number of atoms
	 *@param mol2: Reference to the other molecule
	 */
	float minRmsd(Macromolecule *mol2);
	/** compute minimum rmsd between two molecules with same number of atoms
	 *  It saves rotational matrix to move mol2 over mol1 in order to obtain minimum rmsd
	 *  New "minRmsd" implementation (Mon 1/6/2010) from:
	 *  (It fixes some problems related to highly symmetric structures)
	 * @param mol2: The other molecule
	 * @param mdest: Rotational matrix
	 */
	float minRmsd(Macromolecule *mol2, float mdest[4][4]);
	/** compute minimum rmsd between two molecules with different number of atoms
	 *  It saves rotational matrix to move mol2 over mol1 in order to obtain minimum rmsd
	 *  New "minRmsd" implementation (Mon 1/6/2010) from: (see implementation)
	 *  (It fixes some problems related to highly symmetric structures)
	 * Overloading input to take into account alignment masks at atomic level (17/7/2012)
	 * @param mol2: The other molecule
	 * @param mdest: Rotational matrix
	 * @mask1: Alignment mask at atomic level for macromolecule 1
	 * @mask2: Alignment mask at atomic level for macromolecule 2
	 */
	float minRmsd( Macromolecule * mol2, float mdest[4][4], bool *mask1, bool *mask2 );
	/** compute weighted minimum rmsd between two molecules with same number of atoms
	 *  It saves rotational matrix to move mol2 over mol1 in order to obtain minimum rmsd
	 *  Each residue pair is weighted to identify its importance in the rmsd
	 *  New "minWRmsd" implementation (Mon 1/6/2010) from:
	 *  (It fixes some problems related to highly symmetric structures)
	 * @param mol2: The other molecule
	 * @param mdest: Rotational matrix
	 * @param profile: Residue level weights (One per residue. See Damm & Carlson)
	 */
	float minWRmsd(Macromolecule *mol2, float mdest[4][4], double *profile);
	/** Mon made (17/7/2012)
	 * compute weighted minimum rmsd between two molecules with different number of atoms
	 *  It saves rotational matrix to move mol2 over mol1 in order to obtain minimum rmsd
	 *  Each residue pair is weighted to identify its importance in the rmsd
	 *  New "minWRmsd" implementation (Mon 1/6/2010) from:
	 *  (It fixes some problems related to highly symmetric structures)
	 * @param mol2: The other molecule
	 * @param mdest: Rotational matrix
	 * @param profile: Residue level weights (One per residue. See Damm & Carlson)
	 * @param mask1: Boolean mask containing matching atoms information for molecule 1
	 * @param mask2: Boolean mask containing matching atoms information for molecule 2
	 */
	float minWRmsd( Macromolecule * mol2, float mdest[4][4], double *profile, bool *mask1, bool *mask2 );
	/**
	 * Finds whether two macromolecule collide.
	 * Accessibility must be stored in Macromolecules' occupancy columm
	 * @param mol: Second Macromolecule
	 * @param limit: minimum distant between atoms to avoid collision
	 * @acc: Minimum accessibility value to consider an atom in the surface
	 */
	bool collision(Macromolecule *mol,float limit=1,float acc=0);
	/**Counts the number of clashes of atoms between two structures.
	 *(hydrogens are not taken into account)
	 *@param mol: The second Macromolecule
	 *@param limit: minimum distant between atoms to avoid collision
	 */
	int count_clashes(Macromolecule *mol,float limit=3);
	/**
	 * Finds minimum distance between two Macromolecules and the closer pair of atoms
	 *
	 * @param mol: Second Macromolecule
	 * @param index1: Pointer of the atom of first Macromolecule closer to the second macromolecule
	 * @param index2: Pointer of the atom of second Macromolecule closer to the first macromolecule
	 */
	float minDistance(Macromolecule *mol,int *index1,int *index2);
	/**
	 * Computes pair-to-pair distance matrix between the atoms of the Macromolecule
	 *@param dmatrix: Bidimensional distance matrix
	 */
	void distanceMatrix(double **dmatrix);
	/**
	 * Save list with atom coordinates of macromolecule
	 * @param coord: 3D position list
	 */
	void coordMatrix(float **coord);
	/**
	 * Set macromolecule coordinates from single row coordinates
	 * @param coord: 3D position list
	 */
	void coordMatrixSet(float *coord);
	/**
	 * Save list with atom property in a pointer
	 */
	void getPtrAtomProperty(float **coord);
	/**
	 * Save list with a kind of formfactor in a pointer
	 */
	void getPtrAtomPropertyBFS(float **coord);

	/* new one with Tobi & Soap */
	int ppimatrices(float **cxyz, int **tpatom, int **spatom, int **pfirst, int **cas, int **nres);


	int pdbmatrices(float **coord, int **nres, int **pfirst, int **cas, int **tpatom);

	/**
	 *  Exchange the bfactor from an external table
	 *@param table: table that stores the new bfactor values
	 *Mon's (8/4/2010): if fragment=true "table" has num_fragment elements
	 *                  otherwise "table" will have num_atoms elements (default)
	 */

	void exchange_Pdbfact(double *table, bool fragment=false);
	/**
	 * Obtains the bfactor of all atoms that are stored in a table. This table must be initialized outside of the function
	 *@param table: table where bfactor are stored
	 */
	void get_Pdbfact(double *table);
	/// gets number of atoms in pedestrian way
	int get_num_atoms();
	/// gets number of fragments in pedestrian way
	int get_num_fragments();
	/**
	 * Creates a new system with a selection of elements of the original system
	 * An element will be selected if it matchs any of the Conditions
	 * The old and new system will share the elements selected
	 *@param cond: Collection of Conditions to be applied
	 *@param positive: If true, the elements selected must fit the conditions, otherwise the selected elements will be the ones which do not.
	 *@param include_hetero: If true, to add all Hetero-atoms without any condition check
	 */
	Macromolecule * select( Conditions * cond, bool positive=true, bool include_hetero=false);
	/**
	 * Creates a new system with a selection of elements of the original system
	 * An element will be selected if it matchs any of the Conditions
	 * This function creates copies of the atoms, the new system
	 * does not share the atoms of the current system
	 *@param cond: Collection of Conditions to be applied
	 *@param positive: If true, the elements selected must fit the conditions, otherwise the selected elements will be the ones which do not.
	 */
	Macromolecule * select_cpy( Conditions * cond, bool positive=true );
	/**
	 * Creates a new system with a selection of elements of the original system
	 * An element will be selected if it if it matchs any of the AA names from **aaName
	 * The old and new system will share the elements selected
	 *@param cond: Collection of Conditions to be applied
	 *@param positive: If true, the elements selected must fit the conditions, otherwise the selected elements will be the ones which do not.
	 *@param aaName: Array with pointers to AA 3 letter names (NULL ended)
	 */
	Macromolecule * select_aaName(char **aaName);
	/**
	 *Determines the interaction surface between two macromolecules
	 *@param mol2: The other macromolecule
	 *@param  distance: amstrongs distance between two atoms interacting
	 *@param numContacts: Number of contacs between the macromolecules
	 *@return Double list with pairs of contacts (number of residues in contact)
	 */
	int** getInterface( Macromolecule * mol2, float distance, int *numContacts, char *** , bool withHET=false );
	/**
	 * Returns a  list with the corresponding residues in another macromolecule for a set of atoms in this macromolecule
	 *@param mol2: The other macromolecule
	 *@param list: List of selected residues in this macromolecule
	 *@param distance: distance in amstrong between atoms for considering that two residues could be the same
	 *@param total: number of selected residues
	 *@return list with position of corresponding residues
	 */
	int* relative(Macromolecule *mol2, int *list,float distance, int total);
	///Returns a string list with the names of all the aminoacids that conforms the system
	char * getResInfo();
	/**
	 * Saves backbone coordinates in a vector
	 */
	void coordBackbone(vector<vector<float> > & backbone);
	/**
	 * Saves Alpha Carbons coordinates in a vector
	 */
	void loadCAs(vector < vector<float> > & CAs);
	/**
	 * Saves Backbone coordinates of the first residue in a vector ????
	 */
	void loadNCACFirstAA( vector < vector<float> > & cas );
	/**
	 * Return Macromolecule residue sequence
	 *@param seq: Char list of characters representing the residues in the macromolecule
	 */
	char* get_sequence();

	/** OBSOLETE
	 * Erase all the hydrogens of the Macromolecule;
	 * @param self: Change the macromolecule itself (true) or return a modified copy (false)
	 */
	//Macromolecule *delete_hyd(bool self=true);

	///Check whether occupancy column has asa values
	bool check_asa();
	///Set secondary structure in residues (only for proteins)
	char* secondary_structure(bool output_char=false);
	///Read secondary structure
	char * get_ss();

	///ROTAMER OPERATIONS
	/**
	 *Computes all dihedral angles
	 *@param List of dihedral angles calculated
	 */
	void all_dihedrals(float **);
	/** Write file with dihedral angles. Format AA PHI PSI OMEGA ROT_INDEX
	 * @param List of dihedral angles
	 */
	void show_all_dihedrals(float *, char *);
	/**
	 * Write file with dihedral angles. Format AA PHI PSI OMEGA ROT_INDEX
	 * @param List of dihedral angles
	 * POSSIBLY DEPRECATED
	 */
	void show_all_dihedrals_mon(float *, char *);
	/**
	 * Saves Backbone coordinates of first residue of each segment in a vector
	 */
	void CalculateFirstNCACofSegments( vector < vector<float> > & NCACsPDB);

	// It will make an inter-atomic distance profile between two macromolecules
	// (It allocates memory! When *profile = NULL) (It also returns the RMSD)
	// "mask" and "mask2" are the atom-wise masks to compute inter atomic distances (if != NULL, only CA atoms will be considered)
	double dist_profile( Macromolecule *mol2, double **profile, bool *mask=NULL, bool *mask2=NULL );

	/**
	 * Projects in an existent volume the surface of the macromolecule
	 *
	 *@param vol: volume for projection. It must be created.
	 *@param PR: Dilate/contraction value. Allows smoothing effect. (best value: 4)
	 *@param inner: expand surface (false) shrink surface (true)
	 */





#ifdef VOLUME_INCLUDED
	///SURFACE
	/**
	 * Creates a volume with the surface of the macromolecule
	 *
	 *@param unit: Voxel size
	 *@param PR: Dilate/contraction value. Allows smoothing effect. (best value: 4)
	 *@param inner: expand surface (false) shrink surface (true)
	 */
	vlVolume* projectSurface( float unit, float PR, bool inner=false );
	/**
	 * Projects in an existent volume the surface of the macromolecule
	 *
	 *@param vol: volume for projection. It must be created.
	 *@param PR: Dilate/contraction value. Allows smoothing effect. (best value: 4)
	 *@param inner: expand surface (false) shrink surface (true)
	 */
	void projectSurface( vlVolume *vol, float PR, bool inner=false );

#endif

};
#endif
