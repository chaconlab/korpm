/*
 * io.h
 *
 *  Created on: Oct 15, 2013
 *      Author: mon
 */

#ifndef IO_H_
#define IO_H_


// MON: Parses a Clustal's ".aln" file to determine residue correspondence between two sequences
// type = 1, 2 or 3 (match *, *: or *:. ,respectively)
void read_aln(char *file, bool **p_mask1, int *nres1, bool **p_mask2, int *nres2, int *nres, int type=2);


#endif /* IO_H_ */
