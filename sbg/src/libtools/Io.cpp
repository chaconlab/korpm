/*
 * io.cpp
 *
 *  Created on: Oct 15, 2013
 *      Author: mon
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// MON: Parses a Clustal's ".aln" file to determine residue correspondence between two sequences
// type = 1, 2 or 3 (match *, *: or *:. ,respectively)
void read_aln(char *file, bool **p_mask1, int *nres1, bool **p_mask2, int *nres2, int *nres, int type)
{
	bool debug = false;
	int max = 200;
	int maxseq = 60;
	char line[max]; // readed line buffer
	char seq1[max],seq2[max],homo[max]; // sequence buffers
	int cont,nline,i,length,index,length1,index1,length2,index2,size1,size2;
//	char *fmt = "%*16c%60c";
//	char *fmt = "%*s %60c"; // for seq1 and seq2
	char fmt[20]; // for seq1 and seq2
	char fmt2[20]; // for homology
	sprintf(fmt,"%%*s %%%dc",maxseq);

	bool *mask1 = (bool *)malloc(sizeof(bool)); // this outputs mask1
	bool *mask2 = (bool *)malloc(sizeof(bool)); // this outputs mask2
	*nres = 0; // number of matched residues
	*nres1 = 0;
	*nres2 = 0;

	char *check1; // this will store the first sequence 1st sequence
	char *check2; // this outputs mask2

	if(debug)
	{
		check1 = (char *)malloc(sizeof(char)); // this will store the first sequence 1st sequence
		check2 = (char *)malloc(sizeof(char)); // this outputs mask2
	}

	FILE *f_file;
	if( !(f_file = fopen(file,"r")) )
	{
		printf("Msg(read_aln): Reading clustal's alignment file: %s, Failed!\n",file);
		exit(2);
	}

	cont = 0; // "alignment chunk" line counter
	nline = 0; // absolute index of the line
	while( fgets(line,max,f_file) != NULL ) // Reads the whole file (line by line)
	{
		if(nline >=3) // ignoring the first 3 lines
		{
			if(cont%4 == 0) // each "alignment chunk" repeats each 4 lines
				cont = 0; // reset chunk line counter

			if(cont == 0) // 1st sequence
			{
				sscanf(line,fmt,seq1); // get just sequence
//				seq1[60]='\0'; // mandatory

				// Counting 1st sequence elements
				length1 = 0;
				index1 = 0;
				while( seq1[index1] != ' ' && seq1[index1] != '\n' && index1 < 60 ) // screens full line (valid chunk only)
				{
					if(seq1[index1] != '-')
						length1++; // counting "seq" size form 1st sequence (to avoid "homology" mistakes later)
					index1++;
				}
				*nres1 += length1; // updating 1st sequence length

				// Capping 1st sequence
				seq1[index1]='\0'; // mandatory

				size1 = strlen(line)-strlen(seq1);
				if(debug)
					printf("name size 1= %d  (line= %lu) (seq1= %lu)\n",size1,strlen(line),strlen(seq1));

				if(debug)
				{
					printf("seq1: %s\n",seq1);
					printf("*nres1= %d\n",*nres1);
				}

				// Once sequence length is known, masks memory can be allocated...
				if( *nres1 > 0 && !(mask1 = (bool *) realloc(mask1,*nres1 * sizeof(bool))) )
				{
					fprintf(stderr,"Msg(read_ali): Memory allocation failure (%d)!\nForcing exit!\n",*nres1);
					exit(1);
				}

				// Checking sequence ...
				if(debug)
				{
					printf("Mon1: ");
					if( !(check1 = (char *) realloc(check1,*nres1 * sizeof(char))) )
					{
						fprintf(stderr,"Msg(read_ali): Memory allocation failure!\nForcing exit!\n");
						exit(1);
					}
					length = 0;
					index = 0;
					while( seq1[index] != ' ' && seq1[index] != '\n' && seq1[index] != '\0' && index < 60 ) // screens full line (valid chunk only)
					{
						if(seq1[index] != '-')
						{
							check1[*nres1-length1+length] = seq1[index];
							printf("%c",seq1[index]);
							length++; // counting "seq" size form 1st sequence (to avoid "homology" mistakes later)
						}
						else
							printf("x");
						index++;
					}
					printf("\n");
				}
			}

			if(cont == 1) // 2nd sequence
			{
				sscanf(line,fmt,seq2); // get just sequence
//				seq2[60]='\0'; // mandatory
//				size2 = strlen(line)-strlen(seq2);
//				if(debug)
//					printf("name size 2= %d  (line= %d) (seq2= %d)\n",size2,strlen(line),strlen(seq2));

				// Counting 2nd sequence elements
				length2 = 0;
				index2 = 0;
				while( seq2[index2] != ' ' && seq2[index2] != '\n' && index2 < 60 ) // screens full line (valid chunk only)
				{
					if(seq2[index2] != '-')
						length2++; // counting "seq" size form 2nd sequence (to avoid "homology" mistakes later)
					index2++;
				}
				*nres2 += length2; // updating 2nd sequence length

				// Capping 2nd sequence
				seq2[index2]='\0'; // mandatory

				size2 = strlen(line)-strlen(seq2);
				if(debug)
					printf("name size 2= %d  (line= %lu) (seq2= %lu)\n",size2,strlen(line),strlen(seq2));

				// Some integrity checking...
				if(size1 != size2)
				{
					fprintf(stderr,"Msg(read_aln): Sorry, \".aln\" file corruption detected!\nForcing exit!\n");
					exit(2);
				}

				sprintf(fmt2,"%%*%dc%%%dc",size2-1,(int)strlen(seq2));
				if(debug)
					fprintf(stderr,"fmt2= %s\n",fmt2);

				if(debug)
					printf("seq2: %s\n",seq2);

				if(debug)
					printf("*nres2= %d\n",*nres2);

				// Once sequence length is known, masks memory can be allocated...
				if( *nres2 > 0 && !(mask2 = (bool *) realloc(mask2,*nres2 * sizeof(bool))) )
				{
					fprintf(stderr,"Msg(read_ali): Memory allocation failure (%d)!\nForcing exit!\n",*nres2);
					exit(1);
				}

				// Checking sequence ...
				if(debug)
				{
					printf("Mon2: ");
					if( !(check2 = (char *) realloc(check2,*nres2 * sizeof(char))) )
					{
						fprintf(stderr,"Msg(read_ali): Memory allocation failure!\nForcing exit!\n");
						exit(1);
					}
					length = 0;
					index = 0;
					while( seq2[index] != ' ' && seq2[index] != '\n' && seq2[index] != '\0' && index < 60 ) // screens full line (valid chunk only)
					{
						if(seq2[index] != '-')
						{
							check2[*nres2-length2+length] = seq2[index];
							printf("%c",seq2[index]);
							length++; // counting "seq" size form 1st sequence (to avoid "homology" mistakes later)
						}
						else
							printf("x");
						index++;
					}
					printf("\n");
				}

				// Some integrity checking...
				if(index1 != index2)
				{
					fprintf(stderr,"Msg(read_aln): Sorry, \".aln\" file corruption detected!\nForcing exit!\n");
					exit(2);
				}
			}

			if(cont == 2) // Homology line
			{
				sscanf(line,fmt2,homo); // get just homology

				// Capping homology
				homo[index1] = '\0'; // mandatory

				if(debug)
					printf("homo: %s\n",homo);

				// Parsing 1st sequence into a boolean mask
				length = 0;
				index = 0;
				while( seq1[index] != ' ' && seq1[index] != '\n' && seq1[index] != '\0' && index < 60 ) // screens full line (valid chunk only)
				{
					if(seq1[index] != '-')
					{
//						if(homo[index] == '*')
						if( (type==1 && homo[index] == '*') ||
							(type==2 && (homo[index] == '*' || homo[index] == ':')) ||
							(type==3 && (homo[index] == '*' || homo[index] == ':' || homo[index] == '.')) )
						{
							mask1[*nres1-length1+length] = true; // residue match
							(*nres)++; // counting the number of matched residues
						}
						else
							mask1[*nres1-length1+length] = false; // residue mismatch
						length++; // counting "seq" size form 1st sequence (to avoid "homology" mistakes later)
					}
					index++;
				}

				// Parsing 2nd sequence into a boolean mask
				length = 0;
				index = 0;
				while( seq2[index] != ' ' && seq2[index] != '\n' && seq2[index] != '\0' && index < 60 ) // screens full line (valid chunk only)
				{
					if(seq2[index] != '-')
					{
//						if(homo[index] == '*')
						if( (type==1 && homo[index] == '*') ||
							(type==2 && (homo[index] == '*' || homo[index] == ':')) ||
							(type==3 && (homo[index] == '*' || homo[index] == ':' || homo[index] == '.')) )
							mask2[*nres2-length2+length] = true; // residue match
						else
							mask2[*nres2-length2+length] = false; // residue mismatch
						length++; // counting "seq" size form 1st sequence (to avoid "homology" mistakes later)
					}
					index++;
				}

				if(debug)
					printf("\n\n");
			}
			cont++; // counts the number of line within each "alignment chunk"
		}

		nline++; // number of line absolute index
	}

	// Checking that sequences and their respective masks match!
	if(debug)
	{
		printf("Checking whether sequences and their respective masks do match:\n");

		// 1st sequences are parsed into a boolean mask
		printf("check1(%d): ",*nres1);
		for(i=0; i< *nres1; i++)
		{
			if(mask1[i])
				printf("%c",check1[i]);
			else
				printf("%c",'-');
		}

		// 2nd sequences are parsed into a boolean mask
		printf("\ncheck2(%d): ",*nres2);
		for(i=0; i< *nres2; i++)
		{
			if(mask2[i])
				printf("%c",check2[i]);
			else
				printf("%c",'-');
		}
		free(check1);
		free(check2);
		printf("\nEND\n");
	}

	*p_mask1 = mask1; // this outputs mask1
	*p_mask2 = mask2; // this outputs mask2
	fclose(f_file);
}




