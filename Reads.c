/*
 * Copyright (c) <2008 - 2009>, University of Washington, Simon Fraser University
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this list
 * of conditions and the following disclaimer.
 * - Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or other
 *   materials provided with the distribution.
 * - Neither the name of the <ORGANIZATION> nor the names of its contributors may be
 *   used to endorse or promote products derived from this software without specific
 *   prior written permission.
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * Author         : Faraz Hach
 * Email          : fhach AT cs DOT sfu
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <zlib.h>
#include "Common.h"
#include "Reads.h"


//FILE *_r_fp1;
//FILE *_r_fp2;
//gzFile *_r_gzfp1;
//gzFile *_r_gzfp2;
//gzFile O_r_gzfp1;
//gzFile O_r_gzfp2;
Read *_r_seq;
int _r_seqCnt;
int *_r_samplingLocs;

/**********************************************/
//char *(*readSeq)(char *);
//char *(*readSeq)(char *);

char *(*readSeq)(void *, char *);

/*############################################*/
char *readline_GZ(gzFile * fp, char *seq )
{
	return gzgets(fp, seq, SEQ_MAX_LENGTH);
}
/*############################################*/

/*############################################*/
char *readline_TXT(FILE * fp, char *seq )
{
	return fgets(seq, SEQ_MAX_LENGTH,fp);
}
/*############################################*/


///int tokBUF(char* buffer, char* read_buf, int *BUF_POS){
int scanBUF(char* buffer, char* read_buf, long unsigned int *BUF_POS){
	char * tok;	
	int l;
	tok = strtok(&buffer[*BUF_POS],"\n\r");
	if (tok==NULL){
		return 0;
	}else{
		strcpy(read_buf,&buffer[*BUF_POS]);
		l = strlen(&buffer[*BUF_POS]);
		*BUF_POS+=l+1;
		return 1;
	}
}

int _scanBUF(char* buffer, char * read_buf,long unsigned int *BUF_POS){
	int ret, read_in;
	//printf("->%d"%())
	ret = sscanf(&buffer[*BUF_POS],"%1024[^\n]%n",read_buf,&read_in);
	*BUF_POS+=read_in+1;
	return ret;
}

int readAllReads(char *fileName1,
						char *fileName2,
						int compressed,
						unsigned char *fastq,
						unsigned char pairedEnd,
						Read **seqList,
						unsigned int *seqListSize)
{
	double startTime=getTime();

	char seq1[SEQ_MAX_LENGTH];
	char rseq1[SEQ_MAX_LENGTH];
	char name1[SEQ_MAX_LENGTH];
	char qual1[SEQ_MAX_LENGTH];
	char seq2[SEQ_MAX_LENGTH];
	char rseq2[SEQ_MAX_LENGTH];
	char name2[SEQ_MAX_LENGTH];
	char qual2[SEQ_MAX_LENGTH];

	char dummy[SEQ_MAX_LENGTH];
	char ch;
	int err1, err2;
	int nCnt;
	int discarded = 0;
	int seqCnt = 0;
	int maxCnt = 0;
	int i;
	Read *list = NULL;
	// new vars
	
	char * BUFF_1 = NULL;
	char * BUFF_2 = NULL;
	unsigned long int BUFF_1_pos = 0;		
	unsigned long int BUFF_2_pos = 0;		
	
	char READ_BUFFER[SEQ_MAX_LENGTH];
	size_t read_len;
	size_t curr_pos;

	void * _r_fp1;
	void * _r_fp2;

	if (!compressed)
	{
		_r_fp1 = fileOpen( fileName1, "r");

		if (_r_fp1 == NULL) return 0;

		if ( pairedEnd && fileName2 != NULL ){
			_r_fp2 = fileOpen ( fileName2, "r" );
			if (_r_fp2 == NULL)	return 0;
		}
		else{
			_r_fp2 = _r_fp1;
		}
		readSeq = &readline_TXT;
	}
	else{
		_r_fp1 = fileOpenGZ (fileName1, "r");
		if (_r_fp1 == NULL){
			return 0;
		}

		if ( pairedEnd && fileName2 != NULL ){
			_r_fp2 = fileOpenGZ ( fileName2, "r" );
			if (_r_fp2 == NULL)	return 0;
		}
		else{
			_r_fp2 = _r_fp1;
		}
		readSeq = &readline_GZ;
	}
	
	//READ INTO 1 or 2 buffers
		
	///READ IN read 1 into the buffer
	curr_pos=0;

	while (readSeq(_r_fp1,READ_BUFFER)){
		read_len = strlen(READ_BUFFER);
		BUFF_1=(char*)(realloc(BUFF_1,curr_pos+read_len));
		//strcpy(&BUFF_1[curr_pos],READ_BUFFER);
		memcpy(&BUFF_1[curr_pos],READ_BUFFER,read_len);
		curr_pos+=read_len;
		maxCnt++;	
		//printf("%s",READ_BUFFER);
	}
	
	if (pairedEnd){
		curr_pos=0;
		while (readSeq(_r_fp2,READ_BUFFER)){
			read_len = strlen(READ_BUFFER);
			BUFF_2=(char*)(realloc(BUFF_2,curr_pos+read_len));
			//strcpy(&BUFF_2[curr_pos],READ_BUFFER);
			memcpy(&BUFF_1[curr_pos],READ_BUFFER,read_len);
			curr_pos+=read_len;
		}
	}	
	//printf("%s",BUFF_1);
	//exit(1);
	printf("read in complete\n");
	printf("%d lines\n",maxCnt);
	//printf("%s",BUFF_1);
	if (BUFF_1[0] == '>')
		*fastq = 0;
	else
		*fastq = 1;

	// Counting the number of lines in the file
	//while (readSeq(dummy)) maxCnt++;
	
	if (!compressed){
		fclose(_r_fp1);
		if (pairedEnd) fclose(_r_fp2);
	}else{
		gzclose(_r_fp1);
		if (pairedEnd) gzclose(_r_fp2);
	}
	
	///AFTER HERE, no changes except reads in from stream
	// Calculating the Maximum # of sequences
	if (*fastq)
	{
		maxCnt /= 4;
	}
	else
	{
		maxCnt /= 2;
	}

	if (pairedEnd && fileName2 != NULL )
		maxCnt *= 2;

	list = getMem(sizeof(Read)*maxCnt);

	
	//while( readSeq(name1) )
	while(scanBUF(BUFF_1,name1,&BUFF_1_pos)==1)
	{
		err1 = 0;
		err2 = 0;
		//readSeq(seq1);
		scanBUF(BUFF_1,seq1,&BUFF_1_pos);
		name1[strlen(name1)-1] = '\0';
		for (i=0; i<strlen(name1);i++)
		{
			if (name1[i] == ' ')
			{
				name1[i] = '\0';
				break;
			}

		}

		if ( *fastq )
		{
			scanBUF(BUFF_1,dummy,&BUFF_1_pos);
			scanBUF(BUFF_1,qual1,&BUFF_1_pos);
			//readSeq(dummy);
			//readSeq(qual1);
			qual1[strlen(qual1)-1] = '\0';
		}
		else
		{
			sprintf(qual1, "*");
		}


		// Cropping
		if (cropSize > 0)
		{
			seq1[cropSize] = '\0';
			if ( *fastq )
				qual1[cropSize] = '\0';
		}


		nCnt = 0;
		for (i=0; i<strlen(seq1); i++)
		{
			seq1[i] = toupper (seq1[i]);
			if (seq1[i] == 'N')
			{
				nCnt++;
			}
			else if (isspace(seq1[i]))
			{

				seq1[i] = '\0';
				break;
			}
		}

		if (nCnt > errThreshold)
		{
			err1 = 1;
		}

		// Reading the second seq of pair-ends
		if (pairedEnd)
		{
			scanBUF(BUFF_2,name2,&BUFF_2_pos);
			scanBUF(BUFF_2,seq2,&BUFF_2_pos);
			//readSeq(name2);
			//readSeq(seq2);
			name2[strlen(name2)-1] = '\0';
			for (i=0; i<strlen(name2);i++)
			{
				if (name2[i] == ' ')
				{
					name2[i] = '\0';
					break;
				}
			}

			if ( *fastq )
			{
				//readSeq(dummy);
				//readSeq(qual2);
				scanBUF(BUFF_2,dummy,&BUFF_2_pos);
				scanBUF(BUFF_2,qual2,&BUFF_2_pos);

				qual2[strlen(qual2)-1] = '\0';
			}
			else
			{
				sprintf(qual2, "*");
			}


			// Cropping
			if (cropSize > 0)
			{
				seq2[cropSize] = '\0';
				if ( *fastq )
					qual2[cropSize] = '\0';
			}

			nCnt = 0;
			for (i=0; i<strlen(seq2); i++)
			{
				seq2[i] = toupper (seq2[i]);
				if (seq2[i] == 'N')
				{
					nCnt++;

				}
				else if (isspace(seq2[i]))
				{
					seq2[i] = '\0';
				}
			}
			if (nCnt > errThreshold)
			{
				err2 = 1;
			}
		}

		if (!pairedEnd && !err1)
		{

			int _mtmp = strlen(seq1);
			list[seqCnt].hits = getMem (1+3*_mtmp+3+strlen(name1)+1);
			list[seqCnt].seq = list[seqCnt].hits + 1;
			list[seqCnt].rseq = list[seqCnt].seq + _mtmp+1;
			list[seqCnt].qual = list[seqCnt].rseq + _mtmp+1;
			list[seqCnt].name = list[seqCnt].qual + _mtmp+1;

			reverseComplete(seq1, rseq1, _mtmp);
			int i;

			list[seqCnt].hits[0] = 0;

			for (i=0; i<=_mtmp; i++)
			{
				list[seqCnt].seq[i] = seq1[i];
				list[seqCnt].rseq[i] = rseq1[i] ;
				list[seqCnt].qual[i] = qual1[i];
			}
			list[seqCnt].rseq[_mtmp]=list[seqCnt].qual[_mtmp]='\0';
			sprintf(list[seqCnt].name,"%s%c", ((char*)name1)+1,'\0');
			seqCnt++;
		}
		else if (pairedEnd && !err1 && !err2)
		{
			// Naming Conventions X/1, X/2 OR X
			int tmplen = strlen(name1);
			if (strcmp(name1, name2) != 0)
			{
				tmplen = strlen(name1)-2;
			}
		
			//first seq
			int _mtmp = strlen(seq1);
			list[seqCnt].hits = getMem (1+3*_mtmp+3+tmplen+1);
			list[seqCnt].seq = list[seqCnt].hits + 1;
			list[seqCnt].rseq = list[seqCnt].seq + _mtmp+1;
			list[seqCnt].qual = list[seqCnt].rseq + _mtmp+1;
			list[seqCnt].name = list[seqCnt].qual + _mtmp+1;

			reverseComplete(seq1, rseq1, _mtmp);
			int i;

			list[seqCnt].hits[0] = 0;

			for (i=0; i<=_mtmp; i++)
			{
				list[seqCnt].seq[i] = seq1[i];
				list[seqCnt].rseq[i] = rseq1[i] ;
				list[seqCnt].qual[i] = qual1[i];
			}

			name1[tmplen]='\0';
			list[seqCnt].rseq[_mtmp]=list[seqCnt].qual[_mtmp]='\0';
			sprintf(list[seqCnt].name,"%s%c", ((char*)name1)+1,'\0');
			
			seqCnt++;

			//second seq
			list[seqCnt].hits = getMem (1+3*_mtmp+3+tmplen+1);
			list[seqCnt].seq = list[seqCnt].hits + 1;
			list[seqCnt].rseq = list[seqCnt].seq + _mtmp+1;
			list[seqCnt].qual = list[seqCnt].rseq + _mtmp+1;
			list[seqCnt].name = list[seqCnt].qual + _mtmp+1;

			reverseComplete(seq2, rseq2, _mtmp);

			list[seqCnt].hits[0] = 0;

			for (i=0; i<=_mtmp; i++)
			{
				list[seqCnt].seq[i] = seq2[i];
				list[seqCnt].rseq[i] = rseq2[i];
				list[seqCnt].qual[i] = qual2[i];
			}

			name2[tmplen]='\0';
			list[seqCnt].rseq[_mtmp]=list[seqCnt].qual[_mtmp]='\0';
			sprintf(list[seqCnt].name,"%s%c", ((char*)name2)+1,'\0');

			seqCnt++;
		}
		else
		{
			discarded++;
		}
	}

	if (seqCnt > 0)
	{
		QUAL_LENGTH = SEQ_LENGTH = strlen(list[0].seq);
		if (! *fastq)
		{
			QUAL_LENGTH = 1;
		}
		//fprintf(stderr, "%d %d\n", SEQ_LENGTH, QUAL_LENGTH);
	}
	else
	{
		fprintf(stdout, "ERR: No reads can be found for mapping\n");
		return 0;
	}


	if (pairedEnd)
	{
//		seqCnt /= 2;
	}


	*seqList = list;
	*seqListSize = seqCnt;

	_r_seq = list;
	_r_seqCnt = seqCnt;
	

	free(BUFF_1);
	if (pairedEnd) free(BUFF_2);
	
	fprintf(stdout, "%d sequences are read in %0.2f. (%d discarded) [Mem:%0.2f M]\n", seqCnt, (getTime()-startTime), discarded, getMemUsage());
	//totalLoadingTime+=getTime()-startTime;
	/*	
	fprintf(stdout,"HERE\n");
	int j;
	for (j=0;j<maxCnt;j++){
		fprintf(stdout,"%s\n",list[j].seq);
	}
	exit(1);
	*/	
	return 1;
}

void loadSamplingLocations(int **samplingLocs, int * samplingLocsSize)
{
	int i;
	int samLocsSize = errThreshold + 1;
	int *samLocs = getMem(sizeof(int)*samLocsSize);

	for (i=0; i<samLocsSize; i++)
	{
		samLocs[i] = (SEQ_LENGTH / samLocsSize) *i;
		if ( samLocs[i] + WINDOW_SIZE > SEQ_LENGTH)
			samLocs[i] = SEQ_LENGTH - WINDOW_SIZE;
	}

	*samplingLocs = samLocs;
	*samplingLocsSize = samLocsSize;
	_r_samplingLocs = samLocs;
}

void finalizeReads(char *fileName)
{
	FILE *fp1=NULL;

	if (fileName != NULL)
	{
		fp1 = fileOpen(fileName, "w");
	}
	if (pairedEndMode)
		_r_seqCnt /=2;

	int i=0;
	for (i = 0; i < _r_seqCnt; i++)
	{
		if (pairedEndMode)
		{
			if (_r_seq[2*i].hits[0] == 0 &&  strcmp(_r_seq[2*i].qual,"*")!=0)
			{
				fprintf(fp1,"@%s/1\n%s\n+\n%s\n@%s/2\n%s\n+\n%s\n", _r_seq[i*2].name, _r_seq[i*2].seq, _r_seq[i*2].qual, _r_seq[i*2].name, _r_seq[i*2+1].seq, _r_seq[i*2+1].qual);
			}
			else if (_r_seq[2*i].hits[0] == 0)
			{
				fprintf(fp1, ">%s/1\n%s\n>%s/2\n%s\n", _r_seq[i*2].name, _r_seq[i*2].seq, _r_seq[i*2].name, _r_seq[i*2+1].seq);
			}
		}
		else
		{
			if (_r_seq[i].hits[0] == 0 && strcmp(_r_seq[i].qual, "*")!=0)
			{
				fprintf(fp1,"@%s\n%s\n+\n%s\n", _r_seq[i].name, _r_seq[i].seq, _r_seq[i].qual);
			}
			else if (_r_seq[i].hits[0] == 0)
			{
				fprintf(fp1,">%s\n%s\n", _r_seq[i].name, _r_seq[i].seq);
			}
		}
	}

	fclose(fp1);
	if (pairedEndMode)
		_r_seqCnt *= 2;

	for (i = 0; i < _r_seqCnt; i++)
	{
		freeMem(_r_seq[i].hits,0);
	}


	freeMem(_r_seq,0);
	freeMem(_r_samplingLocs,0);
}
