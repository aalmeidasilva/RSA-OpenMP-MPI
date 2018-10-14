#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <gmp.h>
#include <sys/time.h>
#include <mpi.h>

#ifdef _OPENMP
	#include <omp.h>
#endif

#include "almeidamacros.h"

#define BITSTRENGTH  1024               /* size of modulus (n) in bits */
#define PRIMESIZE    (BITSTRENGTH / 2)  /* size of the primes p and q  */

#define CYPHER_BLOCK_SIZE 100
#define CYPHER_WORK_SIZE 10*CYPHER_BLOCK_SIZE
#define CYPHER_LINE_LEN 310

#define ROOT 0

/* Declare global variables */
mpz_t d,e,n;
mpz_t M,c;

int myRank, numProcs;
char mpiName[MPI_MAX_PROCESSOR_NAME];

/* Declare time-related variables */
struct timeval tv1,tv2;
struct timeval tvdiff;
struct timezone tz;

/* Declare core routines */
void RSA_generateKeys();
int  RSA_checkKeys();
void RSA_encrypt(char *infile);
void RSA_decrypt(char *infile);

/* Initialization related routines */
void initializeGMP();
void clearGMP();
void initializeRandom();

/* Timing routine */
void timediff(struct timeval*,struct timeval*,struct timeval*);

/* Helper routines */
inline char *process(char*);
inline void encrypt(char*, FILE*);
char *readCypher(char *infile);

// MPI Functions
inline void initMPI(int argc, char **argv);

/* Main subroutine */
int main(int argc, char **argv)
{

	if(argc < 3){
		verbose("Usage: ./rsa encrypt/decrypt input_file");
		exit(EXIT_FAILURE);
	}	

	// MPI Stuff
	initMPI(argc, argv);

	/* Initialize the GMP integers first */
	initializeGMP();

	if(myRank == ROOT)
		verbose("Number of available OpenMP threads: %d", omp_get_max_threads());

	// Checking keys
	if(!RSA_checkKeys()){
		if(myRank == ROOT){
			verbose("Creating new RSA Key Files...\n");
			RSA_generateKeys();
		}
	}

	// Encryting
	if(!strncmp(argv[1], "encrypt", 7)){
		if(myRank == ROOT){ // Encryption is not parallel, so it executes only on root process
			verbose("Encrypting file %s", argv[2]);
			RSA_encrypt(argv[2]);
		}
	}
	// Decrypting
	else if(!strncmp(argv[1], "decrypt", 7)){
		if(myRank == ROOT) verbose("Decrypting file %s", argv[2]);
		RSA_decrypt(argv[2]);
	}
	// Invalid operation, aborting
	else{
		verbose("Invalid operation");
		verbose("Usage: ./rsa encrypt/decrypt input_file");
		exit(EXIT_FAILURE);
	}
	
	/* Clear the GMP integers */
	clearGMP();
	MPI_Finalize();

	return 0;
}

/* Initialize all the GMP integers once and for all */
void initializeGMP()
{
	mpz_init(d);
	mpz_init(e);
	mpz_init(n);

	mpz_init(M);
	mpz_init(c);
}

/* Clean up the GMP integers */
void clearGMP()
{
	mpz_clear(d);
	mpz_clear(e);
	mpz_clear(n);

	mpz_clear(M);
	mpz_clear(c);
}

/* This initializes the random number generator */
void initializeRandom()
{
	/* sleep for one second (avoid calls in the same second) */
	// Holy cow this is dumb. Anyway...
	sleep(1);

	/* Set seed for rand() by system time() ... */
	unsigned int time_elapsed;
	time((time_t*)&time_elapsed);
	srand(time_elapsed);
}

/* This function calculates and returns the time
 * difference between two timeval struct
 */
void timediff(struct timeval* a,struct timeval* b,struct timeval* result)
{
	(result)->tv_sec  = (a)->tv_sec  - (b)->tv_sec;
	(result)->tv_usec = (a)->tv_usec - (b)->tv_usec;

	if((result)->tv_usec < 0){
		--(result)->tv_sec;
		(result)->tv_usec += 1000000;
	}
}

/* This function checks whether the keys exist 
 * in the file ~/.rsaprivate and ~/.rsapublic
 */
int RSA_checkKeys()
{
	char publicFile[100];
	char privateFile[100];

	strcpy(publicFile,getenv("HOME"));
	strcpy(privateFile,getenv("HOME"));

	strcat(publicFile,"/.rsapublic");
	strcat(privateFile,"/.rsaprivate");

	FILE* fpublic  = fopen(publicFile,"r");
	FILE* fprivate = fopen(privateFile,"r");

	if((!fpublic) || (!fprivate)){
		/* Key files do not exist */
		return 0;
	}

/*
	verbose("Using RSA Key Files:");
	verbose("Public Key File : %s",publicFile);
	verbose("Private Key File: %s\n",privateFile);
*/
	char d_str[1000];
	char e_str[100];
	char n_str[1000];

	/* Get keys */
	fscanf(fpublic,"%s\n",e_str);
	fscanf(fpublic,"%s\n",n_str);

	fscanf(fprivate,"%s\n",d_str);

	mpz_set_str(d,d_str,10);
	mpz_set_str(e,e_str,10);
	mpz_set_str(n,n_str,10);

	fclose(fpublic);
	fclose(fprivate);

	return 1;
}



void RSA_generateKeys()
{
	/* initialize random seed */
	initializeRandom();

	/* first, record the start time */
	if(gettimeofday(&tv1,&tz)!=0)
		printf("\nWarning: could not gettimeofday()!");

	mpz_t p,q;
	mpz_init(p);
	mpz_init(q);

	char* p_str = new char[PRIMESIZE+1];
	char* q_str = new char[PRIMESIZE+1];

	p_str[0] = '1';
	q_str[0] = '1';

	for(int i=1;i<PRIMESIZE;i++)
		p_str[i] = (int)(2.0*rand()/(RAND_MAX+1.0)) + 48;

	for(int i=1;i<PRIMESIZE;i++)
		q_str[i] = (int)(2.0*rand()/(RAND_MAX+1.0)) + 48;

	p_str[PRIMESIZE] = '\0';
	q_str[PRIMESIZE] = '\0';

	mpz_set_str(p,p_str,2);
	mpz_set_str(q,q_str,2);

	mpz_nextprime(p,p);
	mpz_nextprime(q,q);

	mpz_get_str(p_str,10,p);
	mpz_get_str(q_str,10,q);

	verbose("Random Prime 'p' = %s",p_str);
	verbose("Random Prime 'q' = %s",q_str);

	char n_str[1000];

	mpz_t x;
	mpz_init(x);

	mpz_mul(n,p,q);

	mpz_get_str(n_str,10,n);
	verbose("n = %s",n_str);

	mpz_t p_minus_1,q_minus_1;
	mpz_init(p_minus_1);
	mpz_init(q_minus_1);

	mpz_sub_ui(p_minus_1,p,(unsigned long int)1);
	mpz_sub_ui(q_minus_1,q,(unsigned long int)1);

	mpz_mul(x,p_minus_1,q_minus_1);

	mpz_t gcd;
	mpz_init(gcd);

	unsigned long int e_int = 65537;
	while(true){
		mpz_gcd_ui(gcd,x,e_int);

		if(mpz_cmp_ui(gcd,(unsigned long int)1)==0)
			break;

		/* try the next odd integer... */
		e_int += 2;
	}
	mpz_set_ui(e,e_int);

	char d_str[1000];
	if(mpz_invert(d,e,x)==0){
		verbose("Warning: Could not find multiplicative inverse.");
		verbose("Trying again...");
		RSA_generateKeys();

	}
	mpz_get_str(d_str,10,d);

	printf("\n");

	verbose("Public Keys (e,n):");
	verbose("Value of 'e': %ld",e_int);
	verbose("Value of 'n': %s ",n_str);

	printf("\n");

	verbose("Private Key:");
	verbose("Value of 'd': %s",d_str);

	/* get finish time of key generation */
	if(gettimeofday(&tv2,&tz)!=0)
		printf("\nWarning: could not gettimeofday()!");

	timediff(&tv2,&tv1,&tvdiff);
	verbose("Key Generation took (including I/O) %ld seconds and %d microsecondsi\n", tvdiff.tv_sec, tvdiff.tv_usec);

	/* Write values to file $HOME/.rsapublic and $HOME/.rsaprivate */
	char publicFile[100];
	char privateFile[100];

	strcpy(publicFile,getenv("HOME"));
	strcpy(privateFile,getenv("HOME"));

	strcat(publicFile,"/.rsapublic");
	strcat(privateFile,"/.rsaprivate");

	FILE* fpublic  = fopen(publicFile,"w");
	FILE* fprivate = fopen(privateFile,"w");

	if((!fpublic) || (!fprivate)){
		fprintf(stderr,"FATAL: Could not write to RSA Key Files!");
		exit(EXIT_FAILURE);
	}

	/* Write ~/.rsapublic */
	fprintf(fpublic,"%ld\n",e_int);
	fprintf(fpublic,"%s\n",n_str);

	/* Write ~/.rsaprivate */
	fprintf(fprivate,"%s\n",d_str);

	fclose(fpublic);
	fclose(fprivate);

	verbose("RSA Key files stored:");
	verbose("Public  Key file: %s",publicFile);
	verbose("Private Key file: %s\n",privateFile);

	/* clean up the gmp mess */
	mpz_clear(p);
	mpz_clear(q);
	mpz_clear(x);
	mpz_clear(p_minus_1);
	mpz_clear(q_minus_1);
	mpz_clear(gcd);
}

/* The RSA Encryption routine */
void RSA_encrypt(char *infile)
{
	char pubkeyfile[200];   /* file containing public key */
	char outfile[200];    /* filename to decrypt */
	FILE *fin,*fout;    /* file pointers */
	FILE *fpublic;        /* file pointer to public keyfile */

	int i;            /* string index */
	char chread;        /* character read */
	char stread[CYPHER_BLOCK_SIZE + 1];    /* string read */

	// Opening input file
	fin = fopen(infile,"r");
	if(!fin){
		fprintf(stderr,"FATAL : Could not open %s for reading",infile);
		return;
	}

	// Creating output file name
	sprintf(outfile, "%s.encrypted", infile);
	verbose("Output file: %s", outfile);

	// Opening output file
	fout = fopen(outfile,"w");
	if(!fout){
		fprintf(stderr,"FATAL : Could not open %s for writing",outfile);
		return;
	}

	/* Get time before encryption */
	if(gettimeofday(&tv1,&tz)!=0)
		fprintf(stderr,"\nWARNING : could not gettimeofday() !");

	i = 0;
	chread = 'a';
	do{

		chread = fgetc(fin);

		if(chread == EOF){
			if(i != 0){
				stread[i] = '\0';
				encrypt(stread,fout);
			}
		}
		else if(i == CYPHER_BLOCK_SIZE-1){
			stread[i] = chread;
			stread[i+1] = '\0';
			encrypt(stread,fout);
			i = 0;

		}
		else{
			stread[i] = chread;
			i++;
		}

	} while(chread != EOF);

	/* Get time after encryption */
	if(gettimeofday(&tv2,&tz)!=0)
		printf("\nWarning : could not gettimeofday() !");

	timediff(&tv2,&tv1,&tvdiff);
	verbose("Encryption took %ld seconds and %ld microseconds", tvdiff.tv_sec, tvdiff.tv_usec);

	fclose(fin);
	fclose(fout);
}

/* This function actually does the encrypting of each message */
inline void encrypt(char* msg,FILE* fout){

	unsigned int i;
	char tmps[4];
	char* intmsg = new char[strlen(msg)*3 + 1];

	/* Here, (mpz_t) M is the messsage in gmp integer 
	 * and (mpz_t) c is the cipher in gmp integer */

	char ciphertext[CYPHER_LINE_LEN];

	strcpy(intmsg,"");

	for(i = 0; i < strlen(msg) ; i++){
		/* print it in a 3 character wide format */
		sprintf(tmps,"%03d", (int)msg[i]);
		strcat(intmsg,tmps);
	}

	mpz_set_str(M,intmsg,10);
	delete [] intmsg;

	/* c = M^e(mod n) */
	mpz_powm(c,M,e,n);

	/* get the string representation of the cipher */
	mpz_get_str(ciphertext,10,c);

	/* write the ciphertext to the output file */
	fprintf(fout,"%s\n",ciphertext);
}

char *readCypher(char *infile, int *lines){

	// Opening file		
	FILE *fp = fopen(infile, "r");
	if(!fp){
		perror("readCypher(): Could not open file for reading");
		exit(EXIT_FAILURE);
	}

	// Getting file size
	fseek(fp, 0, SEEK_END); // seek to end of file
	int fsize = ftell(fp); // get current file pointer
	fseek(fp, 0, SEEK_SET); // seek back to beginning of file

	// Allocating Cypher Text memory
	char *cypherText = talloc(char, fsize * 1.02); // Padding for work space
	
	// Reading file
	*lines = 0;
	while(fscanf(fp, "%s\n", &cypherText[(*lines)*CYPHER_LINE_LEN]) > 0)
		(*lines)++;
	
	// Closing file
	if(fclose(fp) != 0){
		perror("readCypher(): Error closing file");
		exit(EXIT_FAILURE);
	}

	return cypherText;
}	

void RSA_decrypt(char *infile)
{

	// Allocating data per thread, preventing race condition
	char *decrypted = talloc(char, omp_get_max_threads() * CYPHER_LINE_LEN);
	mpz_t *M = talloc(mpz_t, omp_get_max_threads());
	mpz_t *c = talloc(mpz_t, omp_get_max_threads());
	for(int i = 0; i < omp_get_max_threads(); i++){
		mpz_init(M[i]);
		mpz_init(c[i]);
	}

	// Reading input file
	int lines;
	char *cypherText = readCypher(infile, &lines);	

	// Calculating MPI indexes
	int myStart, myEnd;
	myStart = (lines / numProcs) *myRank;
	if(lines % numProcs > myRank){
		myStart += myRank;
		myEnd = myStart + (lines / numProcs) + 1;
	}
	else{
		myStart += lines % numProcs;
		myEnd = myStart + (lines / numProcs);
	}
	int myLines = myEnd - myStart;

	// Allocating processed lines and final text memory
	char **processedLines = talloc(char*, myLines);
	char *processedText = talloc(char, (myLines*CYPHER_BLOCK_SIZE) + 1);

	/* Get time before decryption */
	if(gettimeofday(&tv1,&tz)!=0)
		printf("\nWarning : could not gettimeofday() !");

#pragma omp parallel for	
	for(int i = myStart; i < myEnd; i++){
		int tindex = omp_get_thread_num();
		mpz_set_str(c[tindex], &cypherText[i * CYPHER_LINE_LEN], 10);
		mpz_powm(M[tindex], c[tindex], d, n);
		mpz_get_str(&decrypted[tindex * CYPHER_LINE_LEN], 10, M[tindex]);
		processedLines[i-myStart] = process(&decrypted[tindex * CYPHER_LINE_LEN]);
	}

	/* Get time after decription */
	if(gettimeofday(&tv2,&tz)!=0)
		printf("\nWarning : could not gettimeofday() !");
	timediff(&tv2,&tv1,&tvdiff);

	// Copying processed data from lines to full text
	for(int i = 0; i < myLines; i++)
		memcpy(&processedText[i*CYPHER_BLOCK_SIZE], processedLines[i], CYPHER_BLOCK_SIZE);
	processedText[myLines * CYPHER_BLOCK_SIZE] = 0;

	int *displacements = talloc(int, numProcs);
	int *recvCounts = talloc(int, numProcs);
	char *recvBuffer = talloc(char, lines * CYPHER_BLOCK_SIZE);
		
	if(myRank == ROOT){
		displacements[0] = 0;
		for(int i = 0; i < numProcs; i++){
			myStart = (lines / numProcs) * i;
			if(lines % numProcs > i){
				myStart += i;
				myEnd = myStart + (lines / numProcs) + 1;
			}
			else{
				myStart += i;
				myEnd = myStart + (lines / numProcs);
			}
			recvCounts[i] = CYPHER_BLOCK_SIZE * (myEnd - myStart);
			if(i > 0)
				displacements[i] = recvCounts[i-1] + displacements[i-1];
		}
	}
	
	MPI_Gatherv(processedText, CYPHER_BLOCK_SIZE * myLines, MPI_CHAR,
		    recvBuffer, recvCounts, displacements,
		    MPI_CHAR, ROOT, MPI_COMM_WORLD);

	// Wrinting output file
	if(myRank == ROOT){
		char outfile[100];
		sprintf(outfile, "%s.decrypted", infile);
		FILE *fp = fopen(outfile, "w");
		fprintf(fp, "%s", recvBuffer);
		fclose(fp);
	}

	long local__sec = tvdiff.tv_sec;
	long local_usec = tvdiff.tv_usec;
	long global__sec = 0, global_usec = 0;
	long global__sec_min = 0, global_usec_min = 0;
	long global__sec_max = 0, global_usec_max = 0;

	// Reducing time counters
	verbose("Process %d (blocksize: %d) of %d on %s took: %ld seconds and %ld microseconds", myRank, myLines, numProcs, mpiName, local__sec, local_usec);
	MPI_Reduce(&local__sec, &global__sec, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&local_usec, &global_usec, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&local__sec, &global__sec_min, 1, MPI_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&local_usec, &global_usec_min, 1, MPI_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&local__sec, &global__sec_max, 1, MPI_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&local_usec, &global_usec_max, 1, MPI_LONG, MPI_MAX, 0, MPI_COMM_WORLD);

	if (myRank == ROOT){
		verbose("Run Time Summary:");
		printf("Ave: %3lds, %6ldus\nMin: %3lds, %6ldus\nMax: %3lds, %6ldus\n",
			global__sec / numProcs, global_usec / numProcs, global__sec_min , global_usec_min , global__sec_max, global_usec_max);
	}	

//	printf("\n - Decryption took (including output) %ld seconds and %ld microseconds\n", tvdiff.tv_sec, tvdiff.tv_usec);
}

/* This function shows the decrypted integer 
 * message as an understandable text string 
 */
inline char *process(char* str)
{
	unsigned int i=0;
	int tmpnum;
	char strmod[CYPHER_WORK_SIZE];
	char *output = talloc(char, CYPHER_BLOCK_SIZE);

	if(strlen(str)%3 == 1){
		strcpy(strmod,"00");
		strcat(strmod,str);
	}
	else if(strlen(str)%3 == 2){
		strcpy(strmod,"0");
		strcat(strmod,str);
	}
	else
		strcpy(strmod,str);

	while(i<=strlen(strmod)-3){
		tmpnum = strmod[i] - 48;
		tmpnum = 10*tmpnum + (strmod[i+1] - 48);
		tmpnum = 10*tmpnum + (strmod[i+2] - 48);

		sprintf(&output[i/3], "%c", tmpnum);
		i += 3;
	}
	if(i/3 < CYPHER_BLOCK_SIZE)
		output[i/3] = 0;

	return output;
}

inline void initMPI(int argc, char **argv){
	int mpiNameLen;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Get_processor_name( mpiName, &mpiNameLen);
}
