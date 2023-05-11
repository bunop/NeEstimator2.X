// Modifications:
// July 16, 2012: replace extension "out" by "txt"

// Since releasing version 2.01:
// -----------------------------
// In April 2014:
// * Change PrtPop/PrtMonoLoc so that when printing non-polymorphic loci,
//   include popID (variable for sample ID) if only temporal methods are run.
// * Fix an error in printing jackknife confidence intervals for temporal
//   methods (a typo error which causes the lower bound to print "infinite"
//   when the upper bound is).

#include <time.h>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

//#define INFINITE	(float) 9999999
//#define EPSILON		(float) 0.0000001	// used to compare a number with zero
#define INFINITE	(float) 10E10
#define EPSILON	(float) 10E-10	// used to compare a number with zero
#define MAXDEG	2000000000	// 2 billion, maximum for degree of freedom.
#define MAXLONG		4294967295LU	// maximum for unsigned long.
// (For signed long type in 32-bit machine, which has 4 bytes, or 32 bits,
// the largest integer is 2,147,483,647 = 2^{31}-1.
// For unsigned long, largest is 2^{32}-1 = 4,294,967,295)
// When a number, which should be 0 from a calculation but may become nonzero
// by rounding-off error, is used to compare with zero in a condition
// statement, we may use EPSILON instead of zero
#define LEN_LOCUS	10  // max chars for locus names
#define POP_TEMP	20   // chars for population names in temporal method
						// used for printing to output (can be adjusted!)
#define LEN_BLOCK	30
    // LEN_BLOCK = max chars to record when reading a data block from input
	// should be enough to distinguish the names of two populations.
#define GENLEN		10		// maximum length for a genotype
#define LENDIR		250		// max chars for a folder (directory) name
							// (maybe 256 is maximum for Windows?)
#define LENFILE		60		// max chars for a file name
							// FILENAME_MAX is the constant in C
#define PATHFILE	LENDIR+LENFILE

#define FSTAT		1		// for FSTAT format
#define GENPOP		2		// for GENEPOP format
#define FREQUENCY	3		// for frequency format (never used)

#define MINFORM		1		// this should be the minimum of all formats
#define MAXFORM		2		// this should be the maximum of all formats

#define MAXCRIT		10		// maximum number of critical values
#define NCUT_SET	4		// default number of critical values.
#define MAX_SAMP	1000000	// maximum number of samples
#define MAX_POP		1000000	// maximum number of pops
#define MATING  0   // default mating model: 0 for random, 1 for monogamy
#define TABX	0	// default for having tab between columns in tabular-format outputs
// Since the files for Locus data and Burrows coefficients are large,
// especially the latter one, so we limit them:
#define MAXLOCPOP	50	// maximum populations to output in Locus Data file
#define MAXBURRPOP	50	// maximum populations to output in Burrows coef. file

// These are for reading strings
#define WHITESPACE  " \t\f\r\v\n"	// when skipping chars using WHITESPACE
#define CHARSKIP	" ,\t\f\r\v\n"	// or CHARSKIP, allow going new lines
#define PATHCHR		"\\/"			// characters that may represent path,
									// In Windoes, backslash, in Unix, slash
#define STOPCHAR	",\n"
//#define BLANKS		" \t"	// not allow jumping lines for next char
#define BLANKS		" \t\f\r\v"	// not allow jumping lines for next char
#define SPECHR		'*'
#define ENDCHRS		"*\n\f\r\v"	// special char added to new-line char to stop
							// collecting chars for a string
#define XCHRSTOP	" *,\t\f\r\v\n"	// add to CHARSKIP, SPECHR
#define XWHITESTOP	" *\t\f\r\v\n"	// add to WHITESPACE, SPECHR
#define MERGE		0	// when true, jackknife CI includes parameter CI.

#define XFILSUFLD		"xLD.txt"	// append to main file for LD xtra output
#define XFILSUFHET		"xHt.txt"	// append to main file for Het xtra
#define XFILSUFCOAN		"xCn.txt"	// append to main file for Coan xtra
#define XFILSUFTEMP		"xTp.txt"	// append to main file for Temporal xtra
#define EXTENSION	".txt"	// for adding to the prefix of main output file

// The next 3 are default values for LD method, Heterozygote excess method,
// and Nomura method. Set = 1 if Yes, 0 if No.
#define LDACTION	1
#define HETACTION	1
#define COANACTION	1
#define TEMPACTION	1
#define REWEIGH		0	// set = 0 if only weigh sample size when there are
						// missing data (affected methods: LD, Pollak)
						// Can remove this constant and replace REWEIGH
						// by 0 from RunPop0

#define MAXMETHOD	15	// if there are n methods available, this should
						// be set as 2^n-1
#define MAXGENERATION	150	// maximum of generation in temporal method
#define LOCOUTPUT   100		// max # of loci to print details in LocData output
#define LOCBURR		100000	// max # of loci to print in Burrows output
#define BURRSHORT	1	// set = 1 for printing summarized r^2-values for each
						// locus pair, set = 0 for all allele pairs in locus pairs
// The next 2 constants (set = 0 as negation) only take effect
// when Burrows outputted in multi files
#define NONAMEBUR	1	// set = 1 for not using input file name as prefix,
#define NOEXPLAIN	1	// set = 1 for not printing out explanations
						// (only the headers for columns of values are printed)

#define USETMP		1	// set = 1 to use temporary files when possible (At LD)
						// set = 0 then use arrays instead


// for Nomura's method:
#define NONSIBOUT	0	// maximum putative nonsib pairs outputted in outLoc
#define LOCPERLINE	10	// maximum loci per line
#define LOCLIM		100

#define DETAILTEMP	0	// set = 1 if want to have more details in file OUTLOCNAME
#define MAXJACKLD	100000	// maximum number of polymorphic loci in the run
						// to allow calculating jackknife CI in LD
// added in July 2016:
#define MINSAMP	3	// minimum number of individuals for jackknife on samples
#define RESETNE 1	// if set = 0, then r^2, exp(r^2), Ne in LD will NOT be
			// recalculated when there are missing data. The purpose is to
			// produce r^2 on outputs for checking. Normally, set = 1.
// Nov 2016: add this to notify if a Pcrit <= this, it would mean running
// by excluding only singleton alleles.
// Also, in function CritValRead, it will convert Pcrit = 1 into this value,
// so that the user, instead of not knowing which smapll value is
// smaller than this constant, then just enter value "1" instead.
#define PCRITX 10E-8
#define NOSNGL "No S*"

// add in Mar 2015
// -------------------------------------------------------------------------
struct locusMap
{
	int num;
	char name [LEN_LOCUS];
	char chromo[LEN_LOCUS];
};

struct chromosome
{
	char name [LEN_LOCUS];	// name of the chromosome
	int nloci;		// number of loci in the chromosome
	int *locus;		// loci in the chromosome: for i = 0, ..., (nloci-1),
};					// with p = locus[i], then the chromosome contains the
					// (p+1)th locus read in the input file
// -------------------------------------------------------------------------


typedef struct allele *ALLEPTR;
struct allele
{
    int mValue;				// Allele mobility value
    int copy;
    int homozyg;
    float freq;
	float hetx;
    ALLEPTR next;
};

typedef struct fish *FISHPTR;
struct fish
{
    int gene[2];
    FISHPTR next;
};

// for Nomura's method:
typedef struct nonsib *NONSIBPTR;
struct nonsib
{
	int first;
    int second;
//	int simVal;	// this field is not needed
    NONSIBPTR next;
};

// the following is for storing estimate of molecular coancestry
// and the weight of allele frequencies in each locus
typedef struct molecoef *COANPTR;
struct molecoef
{
	int locus;	// this makes the search easier since this structure
	// is used for only loci that are considered (maybe < nloci)
	float fresq;	// for sum of squares of frequencies
    float scoan;
	float diffcoan;
	float weight;	// used for storing weight at locus. The weight
	// at this locus is the product of this field and (1-scoan).
    COANPTR next;
};

// This is for temporal method
typedef struct timefreq *FREQPTR;
struct timefreq
{
    int mValue;				// Allele mobility value
    int *samples;			// number of samples, in different generations
    float *freqs;			// freq. of alleles, in different generations
    FREQPTR next;
};

typedef struct age *AGEPTR;
struct age
{
	float year;
	AGEPTR next;
};


// ------------------------------------------------------------------------
// Function Prototypes: Only list those used by main.
// For the rest, functions called by others should be listed first.
// ------------------------------------------------------------------------
// From GetData
// --------------------------------------------------------------------------
int strcmp0 (char str1[], char str2[]);

int RunDirect (char misFilSuf[]);
int RunMultiFiles (char *mFileName, char mOpt);
int RunMultiCommon (char *mFileName);
int RunOption (char misFilSuf[], char LocSuf[], char BurSuf[],
				char hasOpt, char rem, char *FileOne, char *FileTwo);


int main(int argc, char *argv[])
{

/* We may or may not need to know current directory when the executable
// run, but just in case, try this:
// To get current directory, include the following:
#include <direct.h>
// in Unix, may be:
//#include <unistd.h>
	char currdir[100]="\0";
// Then function
	getcwd (currdir, sizeof(currdir));
// will give current directory
	printf ("curr dir = [%s]\n", currdir);
// If run under workbench,
// argv[0] will give the whole executable name including path name,
// but when run under command line, by typing the name of the program,
// argv[0] is the executable name only (what was typed).

// The Interface should issue a command line including info directive,
// option directive, or the multiple data input files.
// They are preceded by i:, o:, m:, m+:, respectively.
//*/

	char *misFilSuf = "NoDat.txt";
		// added to prefix of output = missing data file name
	char *LocSuf = "Loc.txt";
		// added to prefix of output = Locus Data output file name
	char *BurSuf = "Bur.txt";
		// added to prefix of output = Burrow Coefs output file name
	char *FileOne, *FileTwo;
	int n, p;
	char c;
	char mOpt;
	char hasOpt;
	char rem = 0;

	FileOne = (char *) malloc(PATHFILE * sizeof(char));
	*FileOne = '\0';
	FileTwo = (char *) malloc(PATHFILE * sizeof(char));
	*FileTwo = '\0';
//-----------------------------------------------------------------------
	mOpt = 0;
	if (argc < 2) {
		n = RunDirect (misFilSuf);
		if (n > 1) printf ("*** Number of runs = %d ***\n", n);
	}	// end of run from command line without argument.
	else {
	// run from the command line where arguments were given from the user
	// Each string (besides the name of this program) should start by
	// either 'm', 'm+', 'c', 'i', or 'o'. The next char must be a colon ':';
	// otherwise, the program will stop. The string after those characters
	// should be the name of a directive file for the program to read.

	// The first directive file must be preceded by one of the following
	//	* 'm': for the file that contains list of (multiple) input files
	//	       and some popular options.
	//	* 'm+': for the file that contains list of (multiple) input files
	//	       and more options.
	//  * 'c': for the file that contains common settings followed by list
	//         of input files
	//	* 'i': for the file containing input and output file name.
	// The 'i' file can be followed by a second directive file preceded by
	//	* 'o': This 'o' file is to supplement optiions for 'i'-file.
		n = strlen (argv[1]);
// error messages on screen (no file or not preceded by appropriate chars):
		c = argv[1][0];
		if ((n <= 2) || (c != 'i' && c != 'm' && c!= 'c')) {
			printf ("Illegal argument!\n");
			exit (0);
		};
		if (c == 'i' && argv[1][1] != ':') {
			printf ("Illegal argument!\n");
			exit (0);
		};
		if (c == 'c' && argv[1][1] != ':') {
			printf ("Illegal argument!\n");
			exit (0);
		};
		if (c == 'm') {
			if (argv[1][1] == '+') {
				if (n == 3) {
					printf ("Illegal argument!\n");
					exit (0);
				};
				if (argv[1][2] != ':') {
					printf ("Illegal argument!\n");
					exit (0);
				};
				mOpt = 1;	// run multiple files with more parameters
			} else {
				if (argv[1][1] != ':') {
					printf ("Illegal argument!\n");
					exit (0);
				};
			// things are OK: run multiple files, with less parameters
			};
		};
	// assign FileOne as the name of the first control file to be used
		for (p=0; p<n; p++) {
			if (mOpt == 1) *(FileOne +p) = argv[1][3+p];
			else *(FileOne +p) = argv[1][2+p];
		};
		*(FileOne+n) = '\0';
		if (argv[1][0] == 'm') {
			// "rm" stands for remove, to remove this file
			if (argc > 2) rem = (strcmp0 (argv[2], "rm") == 0) ? 1 : 0;
			n = RunMultiFiles (FileOne, mOpt);

			printf ("\n*** Number of data files = %d ***\n", n);
			if (rem == 1) remove (FileOne);
			exit (0);
		};
		if (argv[1][0] == 'c') {
			// "rm" stands for remove, to remove this file
			if (argc > 2) rem = (strcmp0 (argv[2], "rm") == 0) ? 1 : 0;
			n = RunMultiCommon (FileOne);

			printf ("\n*** Number of data files = %d ***\n", n);
			if (rem == 1) remove (FileOne);
			exit (0);
		};
		if (argv[1][0] == 'i') {
		// first directive file is info on input and output files.
			hasOpt = 0;
			if (argc > 2) {
				hasOpt = (argv[2][0] == 'o' && argv[2][1] == ':')? 1: 0;
				// "rm" stands for remove, to remove those "i", "o" files
				rem = (strcmp (argv[2], "rm") == 0) ? 1 : 0;
				if ((hasOpt == 1) && argc > 3)
					rem = (strcmp (argv[3], "rm") == 0)? 1: 0;
			};
			if (hasOpt == 1) {
			// FileTwo is assigned to be the name of option file:
				n = strlen (argv[2]);
				if (n>2) for (p=0; p<n; p++) *(FileTwo +p) =argv[2][2+p];
				*(FileTwo+n) = '\0';
			};
			// FileOne is the name of info directive file
			// FileTwo is the name of option directive file
			RunOption (misFilSuf, LocSuf, BurSuf, hasOpt, rem,
						FileOne, FileTwo);

		};	// end of "if (argv[1][0] == 'i')"

	};	// end of "if (argc < 2) else .. "
	return 0;

}



//--------------------------------------------------------------------------
// Tools
// -----------------------------

// ------------------------------------------------------------------------

char BinaryDigit (int m, char position)
{
// Suppose m is written as binary, e.g., if m = 9, then m = 1001.
// This function value 0 or 1 as the value of the digit at "position."
// Thus,
// BinaryDigit(9,1)=BinaryDigit(9,3)=1, BinaryDigit(9,2)=BinaryDigit(9,3)=0.
// Algorithm: Write
//      m = a(1)2^0 + a(2)2^1 + ... a(j)2^(j-1) + ... + a(k)2^(k-1),
// where a(j) = 0 or 1 for j <= k, a(k) = 1, and a(j) = 0 for j > k.
// The function will return a(j), where  j = "position." The index k is
// the smallest integer such that 2^k > m.
// Thus, if "position" > k, the return value is 0.
// If position = k, then the return value a(k) is 1, otherwise, write
//    m - 2^(k-1)= a(1)2^0 + a(2)2^1 + ... a(j)2^(j-1) + ... ,
// that is, we replace m by m - 2^(k-1) and then use recursive method.

// This function is used to determine which methods to run.
// With m as input, and parameter "position" as the method, the program
// will read m, and if the value of the function is 1, the method
// associated with "position" will be run.

// m=1: LD (Linkage Disequilibrium) method, (associated with position = 1)
// m=2: Heterozygote method, (associated with position = 2)
// m=4: Molecular Coancestry method (associated with position = 3)
// m=8: Temporal method (associated with position = 4)
// Any sum of those will include the methods represented by the terms.
// For example,
// m=3 = 1+2: both LD and Het methods,
// m=5 = 1+4: both LD and Nomura methods.

// In functions where BinaryDigit is called, input m will be
// assumed to be > 0 and less than 2^(number of methods), which
// is given by constant MAXMETHOD. Thus, BinaryDigit won't be called
// with m being out of this range.
    int k;
	int twopwk, twopwk1;	// twopwk=2^k, twopwk1=2^{k-1}
    if (m <= 0) return 0;	// negative m is seen as 0. Actually, we
			// only need return value when m=0 to make recursive work.
    else {
		twopwk1 = twopwk = 1;
        for (k = 0; twopwk <= m; k++, twopwk1 = twopwk, twopwk *= 2);
        // now twopwk = 2^k > m and twopwk1 = 2^(k-1) <= m
        if (k < position) return 0;	// that is, 2^(position-1) > m
        else {
			if (k == position) return 1;
			return BinaryDigit (m-twopwk1, position);
		};
	};
}


// ------------------------------------------------------------------------
// Compare two strings, case-insensitive; return 0 if they are equal.
// This is case-insensitive vrsion of strcmp function.
int strcmp0 (char str1[], char str2[])
{
	int i, j;
	int m = strlen (str1);
	int n = strlen (str2);
	int c1, c2;
	for (i=0, j=0; i<m && j<n; i++, j++) {
		c1 = tolower(str1[i]); c2 = tolower(str2[j]);
		if (c1 != c2) return (c1-c2);
	};
// Out of the loop, corresponding characters are equal, case insensitive.
// need now to compare the lengths:
	return m-n;
}

// --------------------------------------------------------------------------

int StopSign (char c, char stops[])
{
	int len = strlen(stops);
	int i;
	for (i=0; i < len; i++)
		if (c == stops[i]) return 1;
	return 0;
}


// --------------------------------------------------------------------------

void Reverse (char str[], int start, int stop)
{
// Reverse string str between positions start and stop.
// Suppose str = "abcdefghijk", start=2 (str[2]='c'), stop=9 (str[9]='j'),
// output should be str = "abjihgfedck": the middle "cdefghij" is reversed.
	int k = (stop+1-start)/2;
	// if start = stop, then k = 0: next loop is skipped
	int i;
	char c;
	for (i=0; i<k; i++) {
		c = str[start+i], str[start+i] = str[stop-i], str[stop-i] = c;
	};
}

// --------------------------------------------------------------------------
/* This will not be used:

int SplitName (char *fileName, char *prefix, char *ext,
				int lenPre, int lenExt)
// split fileName into two parts: prefix and extension.
// For prefix, take up to lenPre rightmost chars (no slash "\",
// i.e., no folder name). For ext, take up to lenExt rightmost chars.
// Return the length of the prefix.
// (The extension will be in lower case)
{
	char c;
	int i, j, n, stop;
	int dot = 0;
	i = strlen(fileName);
	// ignore trailing BLANKS:
	for (j=i-1; j>=0 && StopSign(*(fileName+j), BLANKS)==1; j--);
	*(fileName+ j+1)='\0';
	// extract prefix and extension of the file. Extension part is the
	// string after dot '.'. Note that there may be several dots in the
	// file name, only the string after the last dot is the extension,
	// so we need to traverse backward to detect the last one.
	stop = j;	// the index of the last visible character in the file name.
	// count number of dots
	for (n=0, i=0; i<=stop; i++) if (*(fileName + i) == '.') n++;
	// so n is at least 1 if there is a dot in the name
	if (n == 0) dot = 1;	// no extension
	for (n=stop, i=0, j=0; n>=0; n--) {
		c = *(fileName + n);
		if (c == '\\') break;	// no path name is included in prefix
		if (dot == 0 && c!= '.' && i<lenExt) *(ext + i++) = tolower(c);
		if (dot > 0 && j<lenPre) *(prefix + j++) = c;
		if (c == '.') dot++;
	};
	*(ext+i)='\0';
	*(prefix+j)='\0';
	Reverse (ext, 0, i-1);
	Reverse (prefix, 0, j-1);
	return j;
}

*/
// --------------------------------------------------------------------------

int GetPrefix (char *fileName, char *prefix, int lenPre, char path[])
// Get prefix of fileName, take up to lenPre rightmost chars
// starting after char path. If path ='\', that means no folder name
// included in prefix. Return the length of prefix.
{
	char c;
	int i, j, n, stop;
	int dot = 0;
	i = strlen(fileName);
	// ignore trailing BLANKS:
	for (j=i-1; j>=0 && StopSign(*(fileName+j), BLANKS)==1; j--);
	*(fileName+ j+1)='\0';
	// prefix ends at the last dot '.'. Note that there may be several dots
	// in the file name, prefix ends at the last dot.
	stop = j;	// the index of the last visible character in the file name.
	// count number of dots
	for (n=0, i=0; i<=stop; i++) if (*(fileName + i) == '.') n++;
	// so n is at least 1 if there is a dot in the name
	if (n == 0) dot = 1;
	for (n=stop, i=0, j=0; n>=0; n--) {
		c = *(fileName + n);
		// no path name should be included in prefix
		if (StopSign(c, path)==1) break;
//		if (c == path) break;
		if (dot > 0 && j<lenPre) *(prefix + j++) = c;
		if (c == '.') dot++;
	};
// after the previous for loop, prefix is the string read backward
// from the end of fileName, until it reaches character path.
// add in Oct 21, 2011: remove blanks at the end of string prefix,
	for (i=j-1; i>=0 && StopSign(*(fileName+i), BLANKS)==1; i--, j--);

	*(prefix+j)='\0';
	Reverse (prefix, 0, j-1);
	return j;
}

// --------------------------------------------------------------------------

int GetToken (FILE *input, char *token, int maxlen, char skips[],
			  char stops[], int *lastc, int *empty)
// read input, grab the first block of data of "acceptable" characters
// (characters that are not in "stops" list, after skipping over
// leading characters in "skips" list), then
// put it as token, return its length, not to exceed maxlen, i.e.,
// all acceptable characters after this will be ignored. The "lastc"
// serves as the last character in the search, it's not in token.
// After grabbing token, we are still on the line containing it
// unless token is empty.
// The parameter *empty will have value equal to the number of characters
// in "skip" being trimmed from the right before returning "token".
// This value will give information on the number of characters on the
// line, from the last character of "token" to where the reading stops.
{
	char c;
	int i = 0;
	*empty = 0;
	for (; (c=fgetc(input)) != EOF && StopSign(c, skips)==1;);
	if (c != EOF && (StopSign(c, stops)==0)) {
		*token = c;
		// go through all acceptable chars until reaching a char in "stops",
		// but only record up to maxlen.
		for (i=1; (*lastc=c=fgetc(input))!=EOF && (StopSign(c, stops)==0); )
			if (i < maxlen-1) (*(token+ i++)) = c;
		// don't want to go to a new line at the end of this function
		// when we already obtain a nonempty token:
		if (c == '\n') ungetc (c, input);
		// the last char c is in "stops". If it happens that preceding
		// chars are in "skips", i.e., token has trailing characters
		// in "skips", we eliminate them by moving the null terminated
		// character to the last non-skip; increment "empty" so that "empty"
		// will be the number of chars in token being eliminated
		for (; (StopSign (*(token+ i-1), skips)==1); i--, (*empty)++);
		// terminate string here:
		*(token+ i) ='\0';
	} else {
	// if getting nothing (token is empty), we may have gone to next line
	// if either skips or stops contains new line char.
		*token = '\0';
		*lastc = c;
	};
	return i;

}

// --------------------------------------------------------------------------
int Value (char *data)
// convert string of digits in to an integer, only accept digits, no sign.
// If any nondigit is in the string data, return -1.
{
	int i;
	int n=0;
	int k = strlen(data);
	for (i=0; i<k; i++)
		if (!isdigit(*(data+i))) return -1;
	for (i=0; i<k; n=10*n+(*(data+i)-'0'), i++);
	return n;
}


// --------------------------------------------------------------------------
int GetClues (FILE *input, int clueVal[], int nClue, int newline)
{
	// Read a line from input, for a maximum of nClue tokens, to get values
	// from the tokens read, for array clueVal of nonnegative integers.
	// Stop at an invalid one (contains non-digit, e.g., +, -, or a comma)
	// or at SPECHR, or at the end of line.
	// The difference between SPECHR and other non-digit characters is that
	// if a token consists of a valid number immediately followed by
	// SPECHR, which is in XWHITESTOP, then SPECHR is not part of the
	// token returned by GetToken, so it is still a valid integer by the
	// call of Value.
	// Return the number of valid integers, maximum is then nClue.
	// Go to a new line if one of the following occurs:
	// 1. SPECHR is reached (it will not be part of the token read even it
	//    is adjacent on the right of the token on the line)
	// 2. Attempt to read a token will have to pass a line.
	// 3. A token containing non-digit (determined by function Value).
	// 4. Parameter newline was set non-zero at the call.
	// In the first 3 cases, the return value is less than nClue.
	// When the return value is nClue, new line is obtained solely
	// on the value of newline.

	int i, k, c, m, n;
	char *token;
	if (nClue <= 0) return 0;
	token = (char*) malloc(sizeof(char)*10);
	m = newline;
	for (i = 0; i < nClue; ) {
		// use BLANKS so that no attempt to go to next line to get
		// a non-blank char. If the token obtained by GetToken is empty,
		// then either it ends by SPECHR, or cursor is already on next line
//		if (i == 0)	// WHITESPACE has '\n' added to BLANKS. Use it to allow
		// to skip empty lines for the first entry.
//			m = GetToken(input, token, 10, WHITESPACE, XWHITESTOP, &c, &n);
//		else
			m = GetToken(input, token, 10, BLANKS, XWHITESTOP, &c, &n);
//printf("Token = %s, length = %d\n", token, m);
		if (m <= 0) {	// empty token
//			if (c == SPECHR) m = 1;
			if (c != '\n') m = 1;
			break;
		};
		if ((k = Value(token)) < 0) break;
		// token contains some non-digit, and m > 0,
		// cursor will be put to the next line by code below
		clueVal[i++] = k;
		m = newline;
		if (c == SPECHR) break;	// stop reading when stopped by SPECHR
		// (When token stopped by a char in XWHITESTOP, unless it is
		// end-of-line char, the cursor is put behind that stopping char,
		// so, another value may be obtained if it is after this SPECHR,
		// therefore, we need to stop right here!
	};
	// note: if *newline > 0, then after the loop, m always > 0 except when
	// token being empty and the cursor already put to next line.
	// If the loop above did not go through its cycle (i.e., there is break),
	// then the cursor will be on next line no matter the value of *newline,
	// as being made sure by the code below.
	// So, if the return value is less than nClue, new line is obtained.
	// If return value = nClue, the last token read has length > 0, the
	// function GetToken does not put cursor to next line
	if (m > 0) for (; (c=fgetc(input)) != EOF && c!='\n';);
	free(token);
	return i;
}


// --------------------------------------------------------------------------
int GetPair (FILE *input, int *low, int *high, int newline)
{
	int pair[2];
	int k;
	// default for low, high are taken as values assigned before calling
	pair[0] = *low; pair[1] = *high;
	k = GetClues (input, pair, 2, newline);
	*low = pair[0];
	*high = pair[1];
	return k;
}

// --------------------------------------------------------------------------
int GetCluesF (FILE *input, float clueVal[], int nClue, int newline, int *last)
{
	// Counterpart of GetClues, to get the array of float values instead.
	// Condition 3 in the comments at GetClues is replaced by:
	// 3. A token which is not a legitimate float value.
	// In this function, a new return value for a parameter is added:
	// Parameter *last = 1 when the reading ends with SPECHR immediately
	// follows the last legitimate value, else *last = 0.

	int i, c, d, m, n;
	float f;
	char *token = (char*) malloc(sizeof(char)*10);
	*last = 0;
	if (nClue <= 0) return 0;
	d = 0;
	m = newline;
	for (i = 0; i < nClue; ) {
//		if (i == 0)	// WHITESPACE has '\n' added to BLANKS. Use it to allow
		// to skip empty lines for the first entry.
//			m = GetToken(input, token, 10, WHITESPACE, XWHITESTOP, &c, &n);
//		else
			m = GetToken(input, token, 10, BLANKS, XWHITESTOP, &c, &n);
		if (m <= 0) {	// empty token, cursor at SPECHR or passed the line
//			if (c == SPECHR) m = 1;
			if (c != '\n') m = 1;
			break;
		};
		if (sscanf(token, "%f", &f) <= 0) break;
		// now, f is a legitimate float value
		d = c;	// c is the character in XWHITESTOP that stops token
		clueVal[i++] = f;
		m = newline;
		if (c == SPECHR) break;	// stop reading when stopped by SPECHR
		// (When token stopped by a char in XWHITESTOP, unless it is
		// end-of-line char, the cursor is put behind that stopping char,
		// so, another value may be obtained if it is after this SPECHR,
		// therefore, we need to stop right here!
	};
	if (d == SPECHR) *last = 1;
	if (m > 0) for (; (c=fgetc(input)) != EOF && c!='\n';);
	free(token);
	return i;
}


// --------------------------------------------------------------------------
int GetPairI (FILE *input, int *low, int *high, int newline)
{
	// counterpart of GetPair, but allow negative values for low, high.
	float pair[2];
	int k;
	int c;
	pair[0] = (float) *low; pair[1] = (float) *high;
	k = GetCluesF (input, pair, 2, newline, &c);
	*low = (int) pair[0];
	*high = (int) pair[1];
	return k;
}


// --------------------------------------------------------------------------
int GetInt (FILE *input, int *value, int newline)
{
	// similar to GetPairI, but for attempting one integer only.
	int i, k, c, m, n;
	char *token = (char*) malloc(sizeof(char)*10);
	k = 0;
	m = GetToken(input, token, 10, BLANKS, XWHITESTOP, &c, &n);
// WHITESPACE has '\n' added to BLANKS. Use it to allow to skip empty lines.
//	m = GetToken(input, token, 10, WHITESPACE, XWHITESTOP, &c, &n);
	if (m <= 0) {
//		if (c == SPECHR) m = 1;
		if (c != '\n') m = 1;
	} else if (sscanf(token, "%d", &i) > 0) {
		*value = i;
		m = newline;	// so that cursor still on the line if newline = 0
		k = 1;
	};
	if (m > 0) for (; (c=fgetc(input)) != EOF && c!='\n';);
	free(token);
	return k;
}

// --------------------------------------------------------------------------

int GetRanges (FILE *inpFile, int ranges[], int size, int maxVal, char *byRange)
{
// Read ranges in pairs of nonnegative integers, to put in array "ranges"
// whose maximum size is "size" (which should be an even number), then
// return the number of valid ranges. Each pair read is a tentatvive range,
// and those tentative ranges are allowed to overlap.
// Use function GetClues to get a maximum of "size" numbers on the line.
// Let n be the true number of entries obtained by GetClues.
// If n = 0: no number given, so no limit.
// If n > 0 (at least one number read) and the first one is 0: no limit
// If n = 1: one number only. If it is 0, no limit as said above.
// If it is positive, then assume one range, from 1 to that number.

// If n is odd > 1, the last entry (without pairing) will be ignored.
// If n is even, we have nRanges = n/2 pairs.
// If there is an illegitimate pair, then the all entries after that pair
// will be ignored. (A legitimate pair should contain two nonzero numbers
// i, j with i <= j.) So, if the first pair is illegitimate, we assume no
// range is entered, hence we will take all possible values as "no limit"

// For no limit, we set nRanges = 1, ranges[0] = 1, ranges[1] = maxVal
	int i, k, m, n, nRanges;
	for (i=0; i<size; i++) ranges[i] = 0;
	n = GetClues(inpFile, ranges, size, 1);	// the last 1 is to put cursor
	nRanges = n/2;		// to the next line after reading entries for ranges
	*byRange = 0;
	if (ranges[0] == 0) {	// (n=0) => *ranges=0 (default), so this covers n=0
		nRanges = 1;
		ranges[0] = 1;
		ranges[1] = maxVal;
		return nRanges;
	} else if (n == 1) {	// the case (*ranges=0) was excluded
		nRanges = 1;		// This is the case where there is only one entry,
		ranges[1] = ranges[0];	// so it is assumed that there is one range,
		ranges[0] = 1;			// the lower bound is 1, the upper bound is
		*byRange = 1;			// that entry.
		return nRanges;
	};
	// Now, n at least 2 and the first entry is > 0. If the second entry is
	// smaller than the first, then the first pair is an illegitimate pair,
	// we assume no limit as in (*ranges == 0):
	if (ranges[1] < ranges[0]) {
		nRanges = 1;
		ranges[0] = 1;
		ranges[1] = maxVal;
		return nRanges;
	};
	// since loci will be limited by ranges specified here, set *byRange = 1:
	*byRange = 1;
	// The first pair is good. Next, if there is an illegitimate pairing
	// after the first, then will ignore all from that pair.
	for (k = 1; k < nRanges; k++) {
		if (ranges[2*k] > ranges[2*k+1] || ranges[2*k] == 0)
			break;
	};
	nRanges = k;
// The following lines of code are for perfection.
// The calling function RunMultiCommon can do OK without these.
//* -------------------------------------------------------------------
	// Now combine pairings that either overlap or adjacent.
	// Starting from the first range (k=0), look to see any subsequent ranges
	// with common contents, or low-end of one range bigger than high-end
	// of the other by exactly 1. Then combine them and set that subsequent
	// range as empty by reassign the endpoints both to be 0. As we do this,
	// we turn some nonempty ranges to empty. As we move along from k = 0,
	// we will encounter such empty ranges and count them, using variable n.
	// At the end of the next "for" loop, all ranges are disjoint and there
	// are no empty ranges between them. They may not be in ascending order.
	for (k=0, n=0; k+n < nRanges; k++) {
		if (ranges[2*k] == 0) {
		// if this Rk, the(k+1)st range, is empty, reassign it to be the next
		// nonempty range, say S, and assign that range S to be empty.
		// (In other words, swap them.)
			n++;	// since this range Rk is empty, increment n
			for (i=k+1; i < nRanges; i++) if (ranges[2*i] > 0) break;
			if (i < nRanges) {
				ranges[2*k] = ranges[2*i];
				ranges[2*k+1] = ranges[2*i+1];
				ranges[2*i] = 0;
				ranges[2*i+1] = 0;
			};
		};
		// check all ranges after k for possible merging
		for (i = k+1; i < nRanges; i++) {
			if (ranges[2*i] == 0) continue;
			// When low-end of range Ri is in Rk, or is equal to high-end
			// of Ri plus 1, we combine Rk with Ri to obtain new Rk, and
			// set Ri to be empty.
			if (ranges[2*i] >= ranges[2*k] && ranges[2*i] <= ranges[2*k+1]+1)
			{
				if (ranges[2*i+1] > ranges[2*k+1])
					ranges[2*k+1] = ranges[2*i+1];
				ranges[2*i] = 0;
				ranges[2*i+1] = 0;
				continue;
			};
			// Similar to the above, with roles being swapped.
			if (ranges[2*i] < ranges[2*k] && ranges[2*i+1] + 1 >= ranges[2*k])
			{
				ranges[2*k] = ranges[2*i];
				if (ranges[2*i+1] > ranges[2*k+1])
					ranges[2*k+1] = ranges[2*i+1];
				ranges[2*i] = 0;
				ranges[2*i+1] = 0;
				continue;
			};
		};	// merge for range index k is done
	};
	nRanges -= n;
// Now, arrange so that the ranges in ascending order, nicer if outputted
	for (k=0; k < nRanges; k++) {
		m = ranges[2*k];
		n = ranges[2*k + 1];
		for (i = k+1; i < nRanges; i++) {
			if (ranges[2*i] < m) {	// swap Rk and Ri
				ranges[2*k] = ranges[2*i];
				ranges[2*k + 1] = ranges[2*i + 1];
				ranges[2*i] = m;
				ranges[2*i + 1] = n;
				m = ranges[2*k];
				n = ranges[2*k + 1];
			};
		};
	};
// testing:
//printf ("Ranges: ");
//for (k=0; k<nRanges; k++) printf (" [%d, %d] ", ranges[2*k], ranges[2*k+1]);
//printf ("\n");
//*/ ------------------------------------------------------------------
	return nRanges;
}


// --------------------------------------------------------------------------
int SetMethod (int m, char *mLD, char *mHet, char *mNomura, char *mTemporal)
{
	// all methods to run if the input m is at least the max.
	// No method if m <= 0.
	int mCount = 0;
	*mLD = 0;
	*mHet = 0;
	*mNomura = 0;
	*mTemporal = 0;
	if (m >= MAXMETHOD) {
		*mLD = 1;
		*mHet = 1;
		*mNomura = 1;
		*mTemporal = 1;
	} else {
		if (m > 0) {						// Digits of m in binary
			*mLD = BinaryDigit (m, 1);		//		at position 1
			*mHet = BinaryDigit (m, 2);		//		at position 2
			*mNomura = BinaryDigit (m, 3);	//		at position 3
			*mTemporal = BinaryDigit (m, 4);	//	at position 4
		};
	};
	if (*mLD == 1) mCount++;
	if (*mHet == 1) mCount++;
	if (*mNomura == 1) mCount++;
	if (*mTemporal == 1) mCount++;
	return mCount;
}
// --------------------------------------------------------------------------
// added in Nov 2014:
//
// Modified in Dec 2016
// Name for Burrows file according to spec. value cutoff by string "-S",
void GetBurrName(char *outBurrName, int popRead, float cutoff)
// create outBurrName string:
// outBurrName = prefix + "Pop" + popRead + "Bur" + cutoff + ".txt"
// Example: if prefix is empty, popRead = 5, cutoff = 0.02, then
// outBurrName = "Pop5Bur02.txt"
// (Starting after the dot, take cutoff value up to significant digits.)
{
	int i, j;
	char *ptr;
	int len;
	i = popRead;
	char *cutoffstr  = (char *) malloc(20 * sizeof(char));
// Dec 2016: with special cutoff for dropping singleton, use "-S" as suffix
	if (cutoff > 0 && cutoff <= PCRITX) {
		sprintf (cutoffstr, "%s", "-S");
		*(cutoffstr+2) = '\0';
	} else {
		sprintf (cutoffstr, "%f", cutoff);
		ptr = strchr(cutoffstr, '.');
	// remove characters in cutoffstr up to the "dot"
		if (ptr != NULL) cutoffstr = (ptr+1);
		len = strlen(cutoffstr);
		for (j=0, i=len-1; i >=0; i--) {
	// starting from rightmost, j = number of zeros until seeing a non-zero
			if (*(cutoffstr+i) == '0') j++;
			else break;
		}
		len -= j;
		if (len == 0) {	// this is when cutoff = 0
			*cutoffstr = '0';
			len = 1;
		}
		*(cutoffstr+len) = '\0';
	}
	len = strlen(outBurrName);
	sprintf (outBurrName+len, "%s%d%s%s%s",
			"Pop", popRead, "Bur", cutoffstr, ".txt");
	free(cutoffstr);
}


// -----------------------  End of Tools.c    -------------------------------

// GetData  --------------
// --------------------------------------------------------------------------

FILE *GetInpFile (char *inpName, char *prefix, int lenPre, char *format)
// read input file name entered on screen, open file of that name,
// then split the name in two parts: prefix and extension.
{
	FILE *input = NULL;
	int i, j, n, c;
	// allow mistyping input file name twice, but if no name is given, quit:
	for (n=0; n<3 ; n++) {
		*inpName = '\0';
		printf ("> Input file name: ");
		for (i=0; ((c=getchar())!=EOF) && (c!='\n') && (c!=SPECHR);
			*(inpName+i)= c, i++);
// added Dec6,11 this line:
		if (c!='\n') for (; getchar() !='\n';);	// finish the line
		// ignore trailing BLANKS:
		for (j=i-1; j>=0 && StopSign(*(inpName+j), BLANKS)==1; j--);
		*(inpName+ j+1)='\0';
		if (*inpName != '\0') {
			if ((input = fopen(inpName, "r")) == NULL) {
				perror(inpName);
				continue;
			} else {
				for (i=0; i<=j; i++) putchar (*(inpName+i));
				putchar (' ');
				break;
			};
		} else return NULL;	// no entry, quit
	};
	GetPrefix (inpName, prefix, lenPre, PATHCHR);
	*format = FSTAT;
	if (j > 3) {
		if (*(inpName+j-3) == '.' && tolower(*(inpName+j-2)) == 'g' &&
			tolower(*(inpName+j-1)) == 'e' && tolower(*(inpName+j)) == 'n')
		*format = GENPOP;
// for now, we do not accept frequency format
//		if (*(inpName+j-3) == '.' && tolower(*(inpName+j-2)) == 'f' &&
//			tolower(*(inpName+j-1)) == 'r' && tolower(*(inpName+j)) == 'q')
//		*format = FREQUENCY;
	};
	return input;
}

// --------------------------------------------------------------------------
FILE *GetOutFile (char *outName, char *prefix, char mLD, char mHet,
				  char mNomura, char mTemporal)
{
	// the length of prefix in the calling function should be at least 6
	// less than the length of outName, since "Ne.txt" will append to it
	FILE *output = NULL;
	int i, c;
	int append = 0;
	*outName = '\0';
	strcat (outName, prefix);
	i = mLD+mHet+mNomura+mTemporal;
	if (i > 1) {
		strcat (outName, "Ne");
	} else {
		if (mLD == 1) strcat (outName, "LD");
		if (mHet == 1) strcat (outName, "Ht");
		if (mNomura == 1) strcat (outName, "Cn");
		if (mTemporal == 1) strcat (outName, "Tp");
	};
	strcat (outName, EXTENSION);
	append = 0;
// the previous lines should be the same as the commented out below
/*
	int n, len = strlen(prefix);
	for (n=0; n<len; n++) outName[n] = *(prefix+n);
	outName[len] = '\0';
	*outName = *strcat(outName, "Ne.txt");
*/
	printf ("\n> Output will be written to file: %s.\n"
			"> If OK, press <Enter>; else, type in output file name.\n"
			"> To append, insert an asterisk (*) before <Enter> key: ", outName);

	for (i = 0; i < PATHFILE && (c = getchar()) != '\n' && c != SPECHR; ) {
		if (c == '\t') c = ' ';	// put this condition to not accept Tab
								// in the name and replace it by a blank
	// skip over leading white spaces
		if (c == ' ' && i == 0) continue;
		*(outName+i) = c;
		i++;
	};
	if (c == SPECHR) append = 1;
	if (c != '\n') for (; getchar() != '\n';);
// clear trailing blanks - actually unneeded
//	for ( ;i>=0, *(outName+ (i-1))== ' '; i--);
	if (i > 0) *(outName+i) = '\0';
//	else { // empty string for name, revert to default
//		*outName = '\0';
//		strcat (outName, prefix);
//		strcat (outName, "Ne.txt");
//	};
//-------------------------------------------------------------------------

	if (append == 1) output = fopen(outName, "a");
	else output = fopen(outName, "w");
	if (output != NULL) {
		printf ("\nOutput will be written to %s", outName);
		if (append == 1) printf (" (append)");
		printf ("\n");
	// accept prefix from outName upto LENFILE chars. This will be used
	// to assign name for an extra output file, shorter one.
	// This will change parameter prefix
//	GetPrefix (outName, prefix, LENFILE-5, '\\');
		GetPrefix (outName, prefix, LENFILE-5, PATHCHR);
	};
	return output;
}


// --------------------------------------------------------------------------

/*
int GetLocUse (FILE *input, int nloci, char *locUse, int nUse, FILE *locFile)
// Get locus names from input, commas will be ignored along with white spaces
// Print up to 6 characters per name, and up to 100 loci on screen.
// If temporary file locFile is not NULL, write all loci in there, up to 10
// characters per name so that when needed, can be pulled out for printing.
// If unneeded, put NULL in place of this parameter at the calling function.
{
	int p, q, k, n, c;
	char locnam[LEN_LOCUS] = "\0";
	if (nUse < nloci) printf ("Number of loci to be used: %d\n", nUse);
	printf ("Locus names - last 6 characters:");
	if (nUse > 100) printf (" (only the last 100 are listed)");
	printf ("\n");
	for (p=0, q=0; p<nloci; p++){
		if (GetToken(input, locnam, LEN_LOCUS,CHARSKIP,CHARSKIP,&c,&n) <= 0)
		{
			printf("\nOnly %d locus names on input file.\n", p);
			return -1;
		};
		// print to console, up to the last 6 characters in locus name,
		// left aligned ("-" sign in format), and up to 10 loci per line:
		// We only print up to 100 loci, the last ones:
		if (*(locUse+p) != 1) continue;
		k = nUse - q;
		n = strlen (locnam) - 6;
		if (n < 0) n = 0;
		if (k <= 100) printf ("%-7.6s", (locnam+n));
		if (locFile != NULL) {
//			n = strlen (locnam) - 10;	// store locname 10 chars
//			if (n < 0) n = 0;
//			fwrite((locnam+n), 10, 1, locFile);
			fwrite(locnam, LEN_LOCUS, 1, locFile);
		};
		if (q == nUse-1 || (k <= 100 && (q+1) % 10 == 0)) printf ("\n");
		q++;
	};
	return 0;
};


// --------------------------------------------------------------------------


void PrtLocUse (FILE *locFile, FILE *output, int nloci, char *locUse,
				int nlocUse, int nPrt)
//
// Print up to 10 characters per name, and up to nPrt loci from temporary
// file locFile to output, which is for Burrow coefs.
{
	int p, q, n;
	char locnam[LEN_LOCUS] = "\0";
	if (locFile == NULL || output == NULL) return;
	fprintf (output, "Locus names are listed after their designated numberings\n");
	fprintf (output, "(Up to 10 rightmost characters are printed");
	if (nPrt < nlocUse)
		fprintf (output, " and only up to %d names are listed", nPrt);
	fprintf (output, ")\n");
	rewind (locFile);
	for (p=0, q=0; p<nloci && q<nPrt; p++) {
		if (*(locUse+p) != 1) continue;
		fread (locnam, LEN_LOCUS, 1, locFile);
		q++;
		n = strlen (locnam) - 10;
		if (n < 0) n = 0;
		// print up to 10 chars per name, and 5 names per line
		fprintf (output, "%5d:%-12.10s", (p+1), (locnam+n));
		if (q == nPrt || (q % 5 == 0)) fprintf (output, "\n");
	};
	fprintf (output, "\n");
};
*/
// --------------------------------------------------------------------------

// replace the old GetLocUse by this, in Mar 2015
// (replace temporary file by array locList of struct locusMap)
int GetLocUsed (input, nloci, locUse, nUse, locList)
	FILE *input;
	int nloci, nUse;
	char *locUse;
	struct locusMap *locList;
// Get locus names from input, commas will be ignored along with white spaces
// Print up to 6 characters per name, and up to 100 loci on screen.
// If locList is not NULL, write all loci in there, up to 10
// characters per name so that when needed, can be pulled out for printing.
// If unneeded, put NULL in place of this parameter at the calling function.
// added in Mar 2015:
// parameter locList to put names in the list
{
	int p, q, k, n, c;
	char locnam[LEN_LOCUS] = "\0";
	if (nUse < nloci) printf ("Number of loci to be used: %d\n", nUse);
	printf ("Locus names - last 6 characters:");
	if (nUse > 100) printf (" (only the last 100 are listed)");
	printf ("\n");
	for (p=0, q=0; p<nloci; p++){
		if (GetToken(input, locnam, LEN_LOCUS,CHARSKIP,CHARSKIP,&c,&n) <= 0)
		{
			printf("\nOnly %d locus names on input file.\n", p);
			return -1;
		}
		// print to console, up to the last 6 characters in locus name,
		// left aligned ("-" sign in format), and up to 10 loci per line:
		// We only print up to 100 loci, the last ones:
		if (*(locUse+p) != 1) continue;
		k = nUse - q;
		n = strlen (locnam) - 6;
		if (n < 0) n = 0;
		if (k <= 100) printf ("%-7.6s", (locnam+n));
		if (q == nUse-1 || (k <= 100 && (q+1) % 10 == 0)) printf ("\n");
		if (locList != NULL) {
			n = strlen(locnam);
			strncpy(locList[q].name, locnam, n);
			(locList[q].name)[n] = '\0';
			locList[q].num = p;
			(locList[q].chromo)[0] = '\0';
		}
		q++;
	}
	return 0;
}


// --------------------------------------------------------------------------
// Mar 2015: replace the old PrtLocUse by this, using array locList instead
// of temporary file to store array of locus.
void PrtLocUsed (struct locusMap *locList, FILE *output, int nloci, char *locUse,
				int nlocUse, int nPrt)
//
// Print up to 10 characters per name, and up to nPrt loci from the list of
// loci locList to output, which is for Burrow coefs.
{
	int p, q, n;
	if (locList == NULL || output == NULL) return;
	fprintf (output, "Locus names are listed after their designated numberings\n");
	fprintf (output, "(Up to 10 rightmost characters are printed");
	if (nPrt < nlocUse)
		fprintf (output, " and only up to %d names are listed", nPrt);
	fprintf (output, ")\n");
	for (p=0, q=0; p<nloci && q<nPrt; p++) {
		if (*(locUse+p) != 1) continue;
		n = strlen (locList[q].name) - 10;
		if (n < 0) n = 0;
		// print up to 10 chars per name, and 5 names per line
		fprintf (output, "%5d:%-12.10s", (p+1), (locList[q].name +n));
		q++;
		if (q == nPrt || (q % 5 == 0)) fprintf (output, "\n");
	}
	fprintf (output, "\n");
}


// --------------------------------------------------------------------------


char GetInfoDat (FILE *input, int *nPop, int *nloci, int *maxMobilVal,
					int *lenM, int maxlen)
// Get nPop, nloci, maxMobilVal, lenM (=number of digits in alleles)
// from an asummed FSTAT format file, commas will be ignored along with white
// spaces, so those input values are allowed to spread over several lines:
// A return value of 0 indicates that this file cannot be a FSTAT format
// (Can return as 0 if any reading of nPop, nloci, etc, is invalid. However,
// we keep reading so that if needed, we can determine which one is invalid.)
{
	char val = 1;
	int c, n;
	char *data = (char*) malloc(sizeof(char)*maxlen);
	GetToken(input, data, maxlen, CHARSKIP, CHARSKIP, &c, &n);
	if ((*nPop = Value(data)) <= 0) val = 0;
	GetToken(input, data, maxlen, CHARSKIP, CHARSKIP, &c, &n);
	if ((*nloci = Value(data)) <= 0) val = 0;
	GetToken(input, data, maxlen, CHARSKIP, CHARSKIP, &c, &n);
	if ((*maxMobilVal = Value(data)) <= 0) val = 0;
	GetToken(input, data, maxlen, CHARSKIP, CHARSKIP, &c, &n);
	if ((*lenM = Value(data)) <= 0) val = 0;
	free (data);
	return val;
}

// --------------------------------------------------------------------------

int ValidGeno (char *data, int gene[2], int lenM)
// return an error code for validity of a genotype, and assign genotype to
// array gene. Parameter lenM is the length of an allele as a character
// string. We view that a data block is good if it contains only digits and
// its length should not exceed the lengths of two alleles. Each allele
// is supposed to contain "lenM" digits, so if the length of data exceeds
// 2(lenM), it is invalid.

// If a data block contains less than 2(lenM) digits, it is least severe,
// considered as containing only zero value.
// The higher the return value of ValidGeno, the severity of the error:
//	* 0: legitimate genotype, no missing
//	* 1: any of the 2 alleles is zero: missing data
//	* 2: data block contains less than 2(lenM) digits,
//	* 3: data block contains more than 2(lenM) digits.
//	* 4: data block contains non digit characters.
{
	int i;
	int k = strlen(data);
	gene[0]=gene[1]=0;
	for (i=0; i<k; i++)
		if (!isdigit(*(data+i))) return 4;
	if (k > 2*lenM) return 3;
	if (k < 2*lenM) return 2;
	else {	// k = 2*lenM
		for (i=0; i<lenM; gene[0]=10*gene[0]+(*(data+i)-'0'),
					gene[1]=10*gene[1]+(*(data+i+lenM)-'0'), i++);
		if (gene[0]<=0 || gene[1]<=0) return 1;
		return 0;
	}
}

// --------------------------------------------------------------------------

int GetSample (FILE *input, int nloci, int *sampData, int lenM, int *samp,
				int maxlen, int *nSampErr, int *currErr,
				char genErr[], int *firstErr, char *locUse)
// Get one sample from input file, return an error code err.
//	* 0: normal,
//	* between nloci and (2*nloci-1): at least 1 missing data
//	* between 2*nloci and (3*nloci-1): at least 1 genotype too short.
//	* between 3*nloci and (4*nloci-1): at least 1 genotype too wide.
//	* >= 4*nloci: at least 1 genotype has nondigit.
//	* -1: end of file encountered.
// Put data from this one sample into array sampData.
{
	char *data;
	int m = 0, p, c, k, mp;
	int err = 0;	// for error code
	if ((data = (char*) malloc(sizeof(char)*maxlen)) == NULL) {
		printf ("Cannot locate memory for genotypes!\n");
		return -2;
	};

	genErr[0] = '\0';
	*firstErr = -1;
// add Sept 2011:
	*currErr = 0;	// to count total errors in this call
	(*samp)++;
	for (p=0; p<nloci; p++){
		*(sampData+ 2*p) = 0; *(sampData+ 2*p+1) = 0;
		k = GetToken(input, data, maxlen, WHITESPACE, WHITESPACE, &c, &mp);
		if (k <= 0) {
//	no more data in input, although not yet going thru nloci loci, so
//	error warnings here, may decide to terminate the program or just return -1
//  and leave the decision at the calling program, along with some error messages.
			if (*currErr == 0) (*nSampErr)++;
			printf ("Data of sample %d end too soon.\n", *samp);
			free (data);
			return -1;
		} else {
		// skip this checking for error in genotype if this locus is not used
			if (*(locUse+p) == 0) continue;
		// some data here, check for validity, and store genotype data
		// into array sampData, and assign error code err
			mp = ValidGeno (data, (sampData+ 2*p), lenM);
			// if m > 0, missing data or some error in genotype block data,
			// the bigger m, the more serious.
			// assign err only the first time an error occurs, so that
			// we know what locus.
			if (mp > 0) {
				// Only count number of errors once for a sample.
				if (*currErr == 0) {
					(*nSampErr)++;
					*firstErr = p;
				};
				(*currErr)++;
				if (mp > m) {
			// The condition mp > m here is to return error code resulted in
			// the most serious error on the data.
			// as a result: (err/nloci) = mp = errcode from ValidGeno
			// and (err%nloci)+1 = (p+1), the locus order.
					err = mp*nloci + p;
				// string genErr is for the most serious genotype error:
					strncpy (genErr, data, GENLEN);
				};
			// Only print to screen unusual errors.
			// Those errors may not effect the program since this sample is
			// read but not used when option limiting samples/pop is in place
//				if (mp == 2)
//				printf ("Too few digits at locus %d, sample %d: [%s]\n", p+1,
//						*samp, data);
				if (mp == 3)
				printf ("Too many digits at locus %d, sample %d: [%s]\n", p+1,
						*samp, data);
				if (mp == 4)
				printf ("Nondigit at locus %d, sample %d: [%s]\n", p+1,
						*samp, data);
			};
			if (mp > m) m = mp;
		};
	};
	free (data);
	return err;

}

// --------------------------------------------------------------------------

int DatPopID (FILE *input, char *popID, int maxlen)
// read pop id: if new, replace *popID. Return -1 if end of file,
// 0 if not new pop, 1 if new pop.
{
// this "for" loop may cause the program to skip a line unnecessary
// if previous reading ends at a new line character. Thus, to prevent
// such case, put back new line character at the reading before this.
// It was done in GetToken. Therefore, this loop ensures that pop name
// must be at the beginning of an input line.
	int i, c, len;
	char *data = (char*) malloc(sizeof(char)*maxlen);
	for (; (c=fgetc(input)) != EOF && c != '\n';);

	if (GetToken(input, data, maxlen, WHITESPACE, WHITESPACE, &c, &i) <= 0)
	{
		// no more population or samples.
		free (data);
		return -1;
	};
	// comparing strings, case insensitive:
	len = strlen(data);
	if (strcmp0(popID, data) != 0){
	// new pop is read, so its name is reassigned
		for (i=0; i<=len; *(popID+i) = *(data+i), i++);
		free (data);
		return 1;	// new pop
	};
	free (data);
	return 0;
}

// --------------------------------------------------------------------------

int GetnLoci (FILE *input, int maxlen, int *lenM)
// This is for finding the number of loci in GENEPOP format file.
// Look for loci and count them until reaching the keyword "pop",
// which signals the beginning of population. The key word "pop"
// should stand alone, i.e, has length 3 without counting BLANKS.
{
	int c, p, k, n;
	char *data;
//	for (; (c=fgetc(input)) == EOF || c !='\n';);
//	if (c == EOF) return -1;

// ----------  Modification of the above 2 commented lines ----------
// Since there was a very special case encountered that the input file
// does not seem to have "end of line" recognizable, the program just keeps
// reading from the start of input forever.
// That input file happens to be in Windows and no "end of line" character
// when the file is opened by Notepad, but shows normal when opened by
// TextPad. If an attempt to print the first "maxlen" characters of the
// file, the first few characters are not seen in the file!
// Thus, we add the maximum 10000 characters here to return error value
	for (k=0; ((c=fgetc(input))==EOF || c!='\n')&& (k<=10000); k++)
	if (c == EOF || k >= 10000) return -1;
// -------------- End of the modification -----------------------------

	data = (char*) malloc(sizeof(char)*maxlen);
	*data = '\0';
	*lenM = 0;
	p = 0;	// count the number of loci - loci are separated
	// by any of characters in CHARSKIP array. Stop when "pop" is seen.
	// Comparing data block read with string "pop", case insensitive
	c = 0;	// keyword "pop" should not be followed by a comma.
// In the "while" loop below:
// data is the string not including any character terminating it, a character
// in the array stated as the 5th parameter in GetToken (which is CHARSKIP).
// If this parameter is replaced by WHITESPACE, then the comma (not in
// WHITESPACE) can be a part of data (which is wrong).
// If data is not "pop", the "while" loop continues. If it is "pop" but ends
// by a comma, then the "while" loop still continues; so key word "pop" must
// not be immediately followed by a comma.
// The reason is that a locus is allowed to have name "pop" as long as it is
// followed immediately by a comma. If the key word is allowed to be followed
// immediately by a comma, then the condition "c == ','" should be dropped,
// and no locus name can be "pop".
	while (strcmp0(data, "pop") != 0 || c == ',') {
		if ((k=GetToken(input, data, maxlen, WHITESPACE, CHARSKIP, &c, &n)) <= 0)
		{
			if (c == EOF) {
				free (data);
				return -1;
			}
		} else p++;
	// data is legitimate. If it is "pop", increment p will be
	// one too many. Eventually we must arrive at "pop", so will
	// need to decrease p outside the loop
	}
// add in Sept 2016, to skip the line containing the word "pop" so that
// the next call GetToken starts reading at the new line:
	for (; (c=fgetc(input)) == EOF || c !='\n';);
// now, looking beyond the word "pop" for the next pop name,
// (Assume that the first individual starts with the pop name it belongs to.)
// pop name must end by a comma (an element in STOPCHAR)
	n = GetToken(input, data, maxlen, WHITESPACE, STOPCHAR, &c, &k);
	if (c != ',') {	// name can be empty.
//	if (n <= 0 || c != ',') {	// this asks for a nonempty name.
		free (data);
		return -1;
	}
	// The next one must be a genotype if the length is an even
	// number, but no checking whether it contains nondigits
	n = GetToken(input, data, maxlen, WHITESPACE, WHITESPACE, &c, &k);
	free (data);
	if (n==0 || n%2 != 0) return -1;	// error: not conformed to
	else {								// FSTAT or GENEPOP format
		*lenM = n/2;
		return p-1;
	}

}

// --------------------------------------------------------------------------

int GenPopID (FILE *input, char key[], char *popID, int maxlen)
// look for the key word to signal next pop. Return -1 if end of file,
// 0 if not new pop, 1 if new pop.
// Assuming the next data block is either a pop name (current pop)
// or key word "key" (which is string "pop" in GENEPOP format)
// signaling a new pop is to begin.
// After this function, the next reading will be the first genotype
{
	int i, k, c;
	char *data;
// The next "for" loop may cause the program to skip a line wrongly
// if previous reading ends at a new line character. Thus, to prevent
// such case, put back new line character at the previous reading
// before calling this. It was done in GetToken, called before this function.
// Therefore, this loop ensures that the key word "pop" must be at the
// beginning of an input line without causing that skipping line problem.
	for (; (c=fgetc(input)) != EOF && c != '\n';);
	data = (char*) malloc(sizeof(char)*maxlen);
// Change in Sept 2016: want "key" word to match the first string ending by
// a blank, not by a comma (STOPCHAR does not have blank, only end-of-line).
// Then if a match is found, go to the next line for name of the population,
// which is ended by a comma. Since GetToken will not go to new line,
// will need a lind of code to go to new line (added at the beginning of
// the next "if"). This change allows the data file to have strings after the
// "key" word starting new pop (the "key" word is "pop", case insensitive)
// The changes are also applied in function GetnLoci.
//	if ((k = GetToken(input, data, maxlen, WHITESPACE, STOPCHAR, &c, &i)) <= 0)
	if ((k = GetToken(input, data, maxlen, WHITESPACE, WHITESPACE, &c, &i)) <= 0)
	{	// no more population or samples.
		free (data);
		return -1;
	}
// see if data block "data" is the same as keyword "key," case insensitive
//	if (strcmp0(key, data) == 0 && c != ',') {
	if (strcmp0(key, data) == 0) {
// modify in Sep 2016: finish the line after getting "key" word since GetToken
// does not go to next line:
		for (; fgetc(input) != '\n';);
// The next block is new pop, so collect the next pop name,
		k = GetToken(input, data, maxlen, WHITESPACE, STOPCHAR, &c, &i);
//		if (k <= 0) {	// this condition on k rejects any empty pop name.
		// Using this condition on k, return value -1 will cause the
		// program end by RunPop function as it is the end of file.
		// Thus, it may be more informative to have a different return
		// value when k=0 and c=',' and let RunPop inform the error.

		// If c = ',', we may still have data=empty, so excluding only
		// c != ',' will allow pop name empty.
		if (c != ',') {
			free (data);;	//
			return -1;
		} else {	// this data block represents pop name for a new one.
			for (i=0; i<k; *(popID+i) = *(data+i), i++);
			*(popID+k) = '\0';
			return 1;	// new pop
		}
	}
// Added in Apr 2017: when the key word not found, then this is another
// sample of the current pop. The first token just read is stopped by a
// blank. So, if this token is not the key, it should be a sample id
// (the first sample id is used as pop name and was read right after key word
// is found)). There is a comma to end a sample id.
// There could be a case when a sample id is ended by a blank,
// then followed by a comma; in such case, the cursor is positioned before
// the comma, and the comma will be read as a genotype later, which is wrong.
// This "else" is designed to address that problem.
	else if (data[k-1] != ',') {
			for (; (c=fgetc(input)) != '\n' && c != ',';);
	}
	free (data);
	return 0;
}

// --------------------------------------------------------------------------

void PrtMethod (int nMethod, char mLD, char mHet, char mNomura, char mTemporal)
{
	if (nMethod == 0) return;
	printf ("Method(s):");
	if (mLD > 0) {
		printf (" LD");
		nMethod--;
		if (nMethod > 0) printf (",");
		else printf ("\n");
	};
	if (mHet > 0) {
		printf (" Het-Excess");
		nMethod--;
		if (nMethod > 0) printf (",");
		else printf ("\n");
	};
	if (mNomura > 0) {
		printf (" Molecular Coan.");
		nMethod--;
		if (nMethod > 0) printf (",");
		else printf ("\n");
	};
	if (mTemporal > 0) {
		printf (" Temporal");
		nMethod--;
		if (nMethod > 0) printf (",");
		else printf ("\n");
	};
}


// --------------------------------------------------------------------------

int MethodRead (char *mLD, char *mHet, char *mNomura, char *mTemporal,
				int *nGeneration, float timeline[])
{
	int m, n;
	int mCount = 0;	// count the number of methods
	*mLD = 0;
	*mHet = 0;
	*mNomura = 0;
	*mTemporal = 0;
	*nGeneration = 0;
	printf ("\nWhich method(s) to run?\n");
	printf ("  1 = Linkage Disequlibrium\n");
	printf ("  2 = Heterozygote Excess\n");
	printf ("  4 = Molecular Coancestry\n");
	printf ("  8 = Temporal\n");
	printf ("For multiple methods, enter their sum\n");
	printf ("(for all methods, enter %d or larger): ", MAXMETHOD);
	scanf ("%d", &m);
	// need to clear up all remaining characters on this
	// input line, up to and including the end of line character
	// so that they are not part of the next input on screen
	for (; getchar() !='\n';);	// finish the line
	if ((mCount = SetMethod (m, mLD, mHet, mNomura, mTemporal)) <= 0)
		return mCount;
	m = mCount;
	PrtMethod (m, *mLD, *mHet, *mNomura, *mTemporal);
	if (*mTemporal == 0) return mCount;
// from now, take care of inputs for temporal method:
	printf ("For Temporal Method, enter number of samples (at least 2, max = %d): ",
			MAXGENERATION);
	scanf ("%d", nGeneration);
//for (; getchar() != '\n';);
	for (; getchar() !='\n';);	// finish the line
	if (*nGeneration > MAXGENERATION) {
		*nGeneration = MAXGENERATION;
		printf
			("Number of Samples per population is allowed up to %d only.\n", MAXGENERATION);
	};
	if (*nGeneration <= 1) {
		*mTemporal = 0;
		printf ("Temporal method is aborted!\n");
		// decrease mCount and return
		return --mCount;
	};
	n = 0; // count number of mistakes, limit to MAXGENERATION.
	for (m=0; m<*nGeneration; m++) {
		if (m == 0) {
			printf
			("\nEnter generations for %d samples (must be nonnegative, ascending)\n",
				*nGeneration);
			printf
			("Enter -1 (negative 1) to start over, -2 to abort temporal method.\n");
		};
		timeline[m] = 0;
		printf ("* for sample %d: ", m+1);
		scanf ("%f", timeline+m);
		if (timeline[m] == -2) {
			*mTemporal = 0;
			for (; getchar() !='\n';);	// finish the line
			printf ("Temporal method is aborted!\n");
			return --mCount;
		} else {
			if (timeline[m] == -1) {
				m = -1;	// start over
				n = 0;	// reset mistakes
				continue;
			} else if (timeline[m] <0) {
				m--;
				n++;
				printf ("Must be nonnegative!");
				if (n+1==MAXGENERATION) printf (" Last try");
				printf ("\n");
				if (n >= MAXGENERATION) {
					*mTemporal = 0;
					printf ("You reach %d mistakes, Temporal method is aborted!\n",
							MAXGENERATION);
					return --mCount;
				} else continue;
			};
		};

		if (m > 0) {
			if (timeline[m]-timeline[m-1] <= 0) {
				m--;
				n++;
				printf ("Must be in ascending order!");
				if (n+1==MAXGENERATION) printf (" Last try");
				printf ("\n");
				if (n >= MAXGENERATION) {
					*mTemporal = 0;
					printf ("You reach %d mistakes, Temporal method is aborted!\n",
							MAXGENERATION);
					return --mCount;
				} else continue;
			};
		};
	};
	for (; getchar() !='\n';);	// finish the line
	return mCount;
}
// --------------------------------------------------------------------------

int Ordering (float arr[], int size, char asc, char strict) {
// Order array arr in ascending if asc = 1, and descending if asc != 1.
// If strict = 1, return the true number of distinct elements of the array.
// If strict is not 1, keep repeated values (still in desired order), then
// return the number of arranged elements of the array, up to the first
// occurence of the last one in the desired order. For example, if the array
// is 3 1 7 2 1 7, asc=1, strict=0, then the outcome will be: 1 1 2 3 7 7.
// However, the return value is 5, counting from the "1" up to the first "7".
// If the return value n is less than size, then all elements starting from
// position n will be assigned the same value.

// In calling this function here, strict is set = 1 always.
// (So, part of this function is not applied here, but we keep it for
// maybe application elsewhere.)

// This is for array of small sizes, so not designed for time-efficiency.
// For array of n elements, there are about n^2 comparisons (O(n^2)).
	int i, j, n, repeat;
	float left, right, nextVal;
	float *stored;
	if (size <= 0) return 0;
	j = 1;
	right = arr[0];
	left = arr[0];
	for (i=1; i<size; i++) {
		if ((asc==1 && right<=arr[i]) || (asc!=1 && right>=arr[i])) {
		// true when either:
		// * ascending  is required and "right" <= arr[i], or
		// * descending is required and "right" >= arr[i]
			right = arr[i];
			if (right != arr[i] || strict != 1) j++;
		};
		if ((asc==1 && left>arr[i]) || (asc!=1 && left<arr[i])) left=arr[i];

	};
	// In the array, "left" and "right" will be:
	//    * the smallest and largest in ascending case; or
	//    * the largest and smallest in descending case.
	// If j = size, the "if" condition inside the loop happens every time,
	// meaning that the array elements are in desired order.
	if (j == size) return size;	// good order, no repeated values.
	stored = (float*) malloc(sizeof(float)*size);
	n = 0;		// to hold the number of different elements in the array
	repeat = 0;
	for (; ;) {
		for (i = 0; i <= repeat; i++) *(stored + (n++)) = left;
		if (left == right) break;
		// look for the next value of arr that should be on the immediate
		// right of the values recorded in "stored"
		nextVal = right;
		for (i = 0; i < size; i++) {
			// skip the i where arr[i] is already recorded in "stored"
			if ((asc==1 && left>=arr[i]) || (asc!=1 && left<=arr[i]))
				continue;
			if ((asc==1 && arr[i]<=nextVal) || (asc!=1 && nextVal<=arr[i]))
			{
				if (nextVal != arr[i]) repeat = 0;
				else if (strict != 1) repeat++;
				nextVal = arr[i];
			};
		};
		left = nextVal;
		// for the last one, and if strict order not enforced, then the
		// assigned "repeat" below is the extra copies of this value in
		// the array. It's OK to just assign repeat = 0 here, since the
		// assignment in the "for" loop at the end will do the same thing.
		if (left == right && strict != 1) repeat = size-n-1;
	};
	// if n < size, assign "right" to the rest
	for (i = 0; i < size; i++) {
		if (i < n) arr[i] = *(stored+i);
		else arr[i] = right;
	};
	free(stored);
	return n;
}

// --------------------------------------------------------------------------
// Modified Dec 2016:
// For Pcrit >= 0.5, only accept the value exactly equal 1, then
// change it to PCRITX to drop only singleton alleles. Therefore, the user
// can use critical value 1 for this purpose, or can enter a small value
// less or equal PCRITX directly!
int CritValRead (FILE *input, int maxCrit, float critVal[], int *readVal)
{
// Read critical values (Pcrits). Return the number of Pcrits.
// If < zero, something wrong in input.
// The last parameter, when > 0, will indicate Pcrits were read,
// so that we can advance the number of lines read. This information is
// not crucial, however.
// This function will read at most 2 lines, Examples: suppose 2 lines are:
// 3
// 0.1  0.05  0.02
// The first line is for the number of positive Pcrits; if > 0, then the
// second line is for the number of positive Pcrits to be read.
// If there is a negative inserted or an asterisk next to a value, then
// the list stops there. Otherwise, Pcrit = 0 will be added to the list.
// In the above, Pcrits at the end will be: 0.1, 0.05, 0.02, and 0.
// The return value is 4.
// >>> 	If the second line is
// 		0.1  0.05  0.02* (asterisk next to 0.02) or  0.1  0.05  0.02 -1,
// 		Then Pcrits will be 0.1, 0.05, 0.02 (no zero). Return value = 3.
// >>>  If the second line is
// 		0.1*  0.05  0.02 (asterisk next to 0.1) or  0.1  -1  0.05  0.02,
// 		Then Pcrits will be 0.1 only. Return value = 1.
// >>> 	If the second line is
// 		* 0.1  0.05  0.02 (asterisk in front) or  -1  0.1  0.05  0.02,
// 		Then Pcrits will be 0 only. Return value = 1.
// -----------------------------------------------------------------------
	int i, j, m, n;
	int signal;
//	char *token;
	char noZero = 0;	// used to indicate if 0 is included in critVal.
						// Default = 0: 0 is included in the array critVal
	critVal[0] = 0;
	*readVal = 0;
// read n = number of Pcrits, and then put cursor to a new line:
	if (GetInt (input, &n, 1) <= 0) {
		printf("The number of critical values is not given\n");
		return -1;	// no number on the line
	}
	if (n < 0) return -1;	// consider negative for number of Pcrits
							// as an input error; 0 is accepted!
	// n is the number of positive Pcrits to be read.
	// When added by 0, the maximum allowed is maxCrit.
	if (n >= maxCrit) n = maxCrit-1;

// only read freq. values if the n read is positive
	*readVal = n;
//	if (n == 0) return n;
	if (n == 0) {
		critVal[0] = 0;
		return 1;
	}
	// read one extra number for a negative that signals the end of Pcrits,
	// so try to get up to (n+1) values: The return value of parameter signal
	// will play the same role as a negative: if signal = 1 (happens when a
	// special character is next to a critical value read), Pcrits will stop
	// there (in most cases, we will have no zero in the list of Pcrits).
	// The 1 in the parameter list is to always put the cursor to a new line.
	m = GetCluesF(input, critVal, n+1, 1, &signal);	// maximum for m is n+1
	if (m == 0) {	// no number read, so take only Pcrit = 0 and return
		critVal[0] = 0;
		return 1;
	}
	// the stop signal only effected if it is put within n values.
	if (signal != 0 && m <= n) noZero = 1;
	// now looking for a stop signal as a negative value:
	for (i = 0; i < m; i++) if (critVal[i] < 0) break;
	if (i < m) noZero = 1;	// there is a negative read
	else if (m == n+1) i--;	// i=m=n+1, only take n Pcrits read!
	// In the case i = 0, the first one is negative, take only Pcrit = 0
	if (i == 0) {
		critVal[0] = 0;
		return 1;
	}
	// Reorder of the first i elements (before the first negative one if any)
	// to be strictly descending if necessary:
	m = Ordering(critVal, i, 0, 1);	// m is the number of distinct elements.
	// If there is a zero read, it is the mth element if no negative read!
	// Now, get rid of all values >= 0.5:

// Dec 2016: add a line to this "for" loop and then set special value for
// notifying that only single alleles to be dropped. The special value is
// set to replace the value that is exactly equal 1
	char speCrit = 0;
	for (j=0; j<m; j++) {
		if (critVal[j] == 1.0) speCrit = 1;
		if (critVal[j] < 0.5) break;
	}
	m -= j;	// m now is the number of values < 0.5
	for (i=0; i<m; i++) critVal[i] = critVal[i+j];
// Dec 2016: add special value if there was a signal (by "1" value entered)
// Note: Later, we determine if a Pcrit is the special Pcrit by checking if
// it is > 0 and <= PCRITX, so we set it = PCRITX/10, instead of = PCRITX
// to make sure the checking is correct when there is a possibility of
// rounding error (in fact, a setting = PCRITX results in "false" at the
// statement "if critVal[.] <= PCRITX" when we set = PCRITX)
	if (speCrit == 1) critVal[m++] = PCRITX/10;
	// assign the rest to 0 (Pcrit list will consist of the first m):
	for (i=m; i<maxCrit; i++) critVal[i] = 0;
	if (noZero == 0 && critVal[m-1] > 0) m++;	// include 0 in Pcrit list
// for chccking only, cover later:
/*
FILE *check = fopen ("checkNeCrit.txt", "w");
fprintf(check,"\nNumber of Freq read = %d:\n", m);
for (n=0; n<m; n++)fprintf (check, "%e,  ", *(critVal+n));
fprintf (check , "\n\n");
fclose(check);
//*/
	return m;

}


// --------------------------------------------------------------------------
int GeneratnRead (FILE *info, int *nGeneration,
					float timeline[], int maxGen, int *census)
// Read generations until a string not representing a number is
// encountered. Count the total of valid numbers as nGeneration.
// Thus, the input for generations should end with a nondigit character.
// If the values are legitimate, return the number of values read, minus 1.
{

//	float age;
	int c, i, k, m, n;
	int last;
	float *copy;
	if (info == NULL) return 0;
	// the first number of timeline actually is not "timeline", but the
	// census size and we still record it as timeline[0] for convenience
	m = GetCluesF (info, timeline, maxGen, 1, &last);
	n = m-1;
	if (n <= 1) return n;
	if (last != 0) {	// no strictly ascending required
		// for generations, only require that they are distinct.
		copy = (float*) malloc(sizeof(float)*n);
		for (i=0; i<n; i++) *(copy+i) = timeline[i+1];
		k = Ordering(copy, n, 1, 1);	// order copy in strict ascending
		free (copy);	// no longer needed
		if (k < n) {
			printf ("Error: Generations are not distinct:\n    ");
			for (i=1; i<m; i++) printf ("%6.1f", timeline[i]); printf ("\n");
			return -1;
		} else {
			*nGeneration = n;
		// timeline[0] is for census size. If it is 0, plan II; else plan I
			*census = (int) timeline[0];
			return n;
		};
	};
	// Now, strict ascending of generations is required.
	// make sure the times are in ascending order, starting from index 1
	for (c=1; c<n; c++)	// from index 1, timeline is for time
		if (timeline[c] >= timeline[c+1]) break;
	// if m is at least 3 (one for census size, and at least 2 generations),
	// then the loop executes at least once. In such case, the "for" loop
	// executes to the end if and only if the times read are in ascending
	// order. If the loop executes fully, then c = n
	// Only get new value for nGeneration if at least 2 ages are read:
	if (c < n) {
		// not in ascending order, return -1 so that no more time line
		// will be read for the rest of the current input file.
		// The return value -1 is special for this case so that error
		// messages could be posted elsewhere if desired.
		printf ("Error: Generations not strictly increasing:\n    ");
		for (i=1; i<m; i++) printf ("%6.1f", timeline[i]); printf ("\n");
		return -1;
	} else {
		*nGeneration = n;
		// timeline[0] is for census size. If it is 0, plan II; else plan I
		*census = (int) timeline[0];
		return n;
	};
}


// --------------------------------------------------------------------------

FILE *Prompt (char *inpName, char *prefix, int lenPre, int *nPop,
				int *nloci, int *maxMobilVal, int *lenM, char *format,
				char *mLD, char *mHet, char *mNomura, char *mTemporal,
				int *nGeneration, float timeline[])
// user interact, asking for input file, then determine if it is
// of FSTAT or GENEPOP format. Position pointer to the file at the
// beginning of the list of locus names, return pointer to input file.
// Also return prefix of input file for use in creating default output.
{
	FILE *input = NULL;
	int m, n;
	input = GetInpFile(inpName, prefix, lenPre, format);

	if (input == NULL) {
		printf ("No input file is given, Program aborted!\n");
		exit (EXIT_FAILURE);
	};
	m = strlen(prefix);

// read DAT file:
	if (*format == FSTAT) {
		printf ("(FSTAT format)\n");
		if (GetInfoDat (input, nPop, nloci, maxMobilVal, lenM, LEN_BLOCK) == 0)
		{
			printf ("Top lines of input file indicate this is not FSTAT format,\n"
				"Assuming now it is of GENEPOP format.\n");
			*format = GENPOP;
			rewind (input);
		};
	};
// format may change value at previous "if"
	if (*format == GENPOP) {	// determine nloci
		printf ("(GENEPOP format)\n");
		*nloci = GetnLoci (input, LEN_BLOCK, lenM);
		if (*nloci <= 0) {
			printf("\nError in input file, program aborted.\n");
			exit (EXIT_FAILURE);
		} else {
			rewind (input);
			for (; (n=fgetc(input)) !='\n';);	// position at first locus
		};
	// should assign maxMobilVal to make sure it doesn't get garbage.
	// Since the length is lenM, assign max possible value:
		for (m=1, n=1; n<=*lenM; n++) m *= 10;
		*maxMobilVal = --m;
	};
	printf ("Number of loci = %d, %d-digit alleles\n", *nloci, *lenM);
	if (MethodRead (mLD, mHet, mNomura, mTemporal, nGeneration, timeline) == 0)
	{
		printf ("No method is given!\n");
		exit (0);
	};
	return input;
}


// --------------------------------------------------------------------------
void GetXoutName (char *xOutName, char *outName, int lenFile,
				  char suffix[], char path[])
{
	char prefix [PATHFILE] = "\0";

// when path = character '\', GetPrefix get file name without path
// when path = null character '\0', get whole name
	int m;
	m = lenFile - strlen (suffix);
	GetPrefix (outName, prefix, m, path);
	*xOutName = '\0';
	strcat (xOutName, strcat (prefix, suffix));
}


// --------------------------------------------------------------------------
// The rest are for dealing with files that list all parameters to run the
// program. A few of those parameters can be entered from the console when
// running the program, but it will be cumbersome for the user to reply to
// all on the console. Thus, the files are needed and can be run directly
// from the console with a command line referring to them in "main" function.

// --------------------------------------------------------------------------

void ErrMsg (char *inpName, char *msg, int lineNum)
{
	printf ("Error in file %s at or about line %d:\n*** %s ***\n\n",
		inpName, lineNum, msg);
}
// --------------------------------------------------------------------------
int FindMethod(FILE *info, char *infoName, int *line,
				char *mLD, char *mHet, char *mNomura, char *mTemporal,
				int *tempClue)
{
//	int c;
	int m, n;

	n = 0;
	m = 0;
	GetPair(info, &n, &m, 1);
	*tempClue = m;
	(*line)++;
	// now determine which methods to run
	if (SetMethod (n, mLD, mHet, mNomura, mTemporal) <= 0)
	{
		ErrMsg (infoName, "No method is given!", *line);
		return -1;	// the main program will exit
	};
	return 0;
}

// --------------------------------------------------------------------------


int LociDropped (FILE *optFile, char *locArr, int nloci, int *linedone,
				 int maxline, char mode, char *byRange)
// Read loci to be dropped:
// Return number of loci dropped nlocSkip, -1 if input error.
//
// If mode = 0, assign *(locArr+p) = 0 if locus (p+1) is to be dropped.
// On input, locUse will take place of locArr.
//
// If mode = 1, assign locArr[i] = p > 0 if locus p is dropped, and set
// locArr[j] = 0 for j >= number of dropped alleles nlocSkip. Then in the
// calling function, will need to loop over the first nlocSkip of locArr
// to assign locUse[p-1] = 0 whenever locArr[i] = p > 0. This mode is used
// when we don't know the number of loci yet (used in RunMultiFiles).
// (mode = 1 when no range allowed)

// If maxline > 0, only read a maximum of "maxline" lines to find loci
// dropped.

// return *byRange = 1 if loci to be used are given by specifying ranges
// The first part is to read data on a line to determine number of
// loci to be dropped or ranges of loci to be used.
//

{
	int nlocSkip, c, i, j, p, n, line;
//	int maxloc = nloci;
	char *token;
	int m, num;
	int low, high;
	int *included;	// to mark loci
	*byRange = 0;	// default: loci dropped will be listed

	c = 0;	// set = 0 so that GetPair still let cursor on the current line
	// unless it reaches SPECHR or passed the line because token = empty.
	// If that happens, then c is set to 1 in GetPair.
	// Always set low, high before calling GetPair:
	low = 0; high = 0;
	num = GetPair(optFile, &low, &high, c);	// num = 0, 1, or 2: numbers read
	(*linedone)++;
	// if k < 2, cursor will be on the next line.
	if (num == 0) return -1;	// no number was read
	// num = 1: only one number on the line, should be #loci dropped
	// num = 2: two numbers on the line, but if mode = 1: only take loci
	// dropped, no ranges of loci (so the second number "high" is irrelevant)
	// Since num = 2, cursor is still on the current line by GetPair,
	// need to move to the next line (linedone already incremented)
	if (num == 2 && mode == 1) for (; (c=fgetc(optFile)) != EOF && c!='\n';);
	if (low <= 0) return 0;
	// read loci dropped starting from the next line if only one number read:
	if (num == 1 || mode == 1) {
		nlocSkip = low;	// which is > 0
		if (nlocSkip > nloci) nlocSkip = nloci;
	// loci to be deleted might be separated by commas, so use GetToken
	// instead of fscanf (which will treat comma as illegitimate).
	// The numbering of a locus cannot be more than 10 digits, unless
	// it is preceded by a lot of zeros, so token of size 10 is enough:
		token = (char*) malloc(sizeof(char)*10);
		line = 0;
		if (maxline <= 0) maxline = -1;	// set = -1 so that in the next loop,
	// line will never equal to maxline, then the loop stops by value i only.
		m = 1;
		for (i = 0, j = 0; i < nlocSkip && line != maxline; ) {
		// use CHARSKIP to allow going next line for loci dropped,
		// note that CHARSKIP contains white space and comma, XCHRSTOP
		// have SPECHR added.
		// Example: line contains (enclosed in quotes): "10,2"
		// Then this loop will get 2 tokens: 10, 2.
		// Comma is a part of XCHRSTOP, so the first token will be "10".
		// If comma is not a part, then the first token will be "10,2".
		// If comma is not a part, and the line is "10, 2" (a blank between
		// comma and 2), then will get 2 tokens 10, 2.
		// So, by having comma as a stopping character, we allow commas
		// between numbers without having blanks between.
			m = GetToken(optFile, token, 10, CHARSKIP, XCHRSTOP, &c, &n);
			if (m == 0) {
			// about reassigning m, see comment at CritValRead
				if (c == SPECHR) m = 1;
				break;
			};
			if (c == '\n') {
				(*linedone)++;
				line++;
			};
		// locus p is to be removed:
			p = Value(token);
		// increment j only when this p is valid for locus numbering
			if (p <= nloci && p-1 >= 0) {
				if (mode == 0) *(locArr+ p-1) = 0;
				else *(locArr+ i) = p;
				j++;
			};
			if (p >= 0) i++;	// increment i, as long as the number read is
			// a valid number, not necessarily representing a locus.
			else break;	// no more reading if illegal entry (no number)
			if (i == nlocSkip) break;
		};
		// j is the true number of loci being skipped:
		nlocSkip = j;
		// only finish the line if token has positive length
		if (m > 0) {
			// if (c == '\n') then *linedone is already incremented.
			// if not, we finish the line, so need to increment *linedone
			if (c != '\n') (*linedone)++;
			for (; (c=fgetc(optFile)) != EOF && c!='\n';);
		};
		free(token);
		return nlocSkip;
	};
	// from this, k = 2 and mode = 0 (range allowed)
	// loci are given by ranges, need to read ranges:
	if (high < low) {
		// first range has error: no loci dropped, but line still at current
		// because k = 2, so need to finish the line
		for (; (c=fgetc(optFile)) != EOF && c!='\n';);
		return 0;
	};
	// although it is OK to work directly with parameter locArr (resetting
	// all values to 0), but to be safe, allocate another array as to not
	// disturb the previously assigned values of locArr when it enters this
	// function, just in case some of its values should be kept when the
	// program is modified later.
	included = (int*) malloc(sizeof(int)*nloci);
	for (i = 0; i < nloci; i++) *(included + i) = 0;
	*byRange = 1;	// loci are determined by range
	c = 0;	// reassign c, for parameter in GetPair
	// Entering this loop, we have 0 < low <= high
	while (num == 2) {
	// any pairs (low, high) read will be seen as a range to be added.
	// If low > high, the range is empty. If low > nloci, the range is
	// also empty. Any range is limited to loci from 1 to nloci.
	// Allow ranges to overlap.
		for (i = low; i <= high && i<= nloci; i++) *(included+i-1) = 1;
	// read next range
		num = GetPair(optFile, &low, &high, c);
	};
	nlocSkip = 0;
	for (p = 0; p < nloci; p++) if (*(included+p) == 0) nlocSkip++;
	// make sure all ranges combined are not empty: if so, must be an error
	if (nlocSkip < nloci) {
		for (p = 0; p < nloci; p++) *(locArr+p) = *(included+p);
	} else nlocSkip = 0;	// since cannot have all loci dropped, reset
	free(included);
	return nlocSkip;
}

// --------------------------------------------------------------------------

AGEPTR AgeSeq (FILE *infofile, int *nSeq, int *nPlan)
{
	int n, count;
	float timeline[MAXGENERATION+1];
	AGEPTR head, newptr, preptr;
	int i, istop;
	int census;
	char plan1 = 0, plan2 = 0;
	int size = sizeof(struct age);
	int numGen = 0;
	head = NULL;
	istop = 2;
	// a generetion set includes a census size (0 if plan II), and at least
	// 2 generations; so to be a generation set, return value of GeneratnRead
	// should be at least 2, i.e. the number of entries on the line >= 3.
	count = 0;
	while (istop > 1) {
		for (i = 0; i <= MAXGENERATION; timeline [i++] = 0);
		// the value of n in the next call is not used here
		istop = GeneratnRead (infofile, &n, timeline, MAXGENERATION+1, &census);
		if (istop > 1) {	// valid (new) generation batch
			if (census > 0) plan1 = 1;
			else plan2 = 1;
		}
		// So, istop <= 1 when there are less than 2 values before a non-digit
		// character or the end of the line, then the loop is done after this.
		// Put all values in timeline as a sequence, which will be retrieved
		// by GetGeneration and recognized if it is a valid generation set.
		for (i=0; i<=istop; i++) {
			// Every number read is put as a node
			if ((newptr = (AGEPTR) malloc(size)) != NULL) {
				newptr->year = timeline[i];
				newptr->next = NULL;
				if (numGen == 0) {
					head = newptr;
				} else {
					preptr->next = newptr;
				};
				preptr = newptr;
				numGen++;
			};
		};
		if (numGen == 0) break;	// when numGen still 0, preptr is undefined,
		// which will cause error at the next code. This (numGen = 0) happens
		// only when istop < 0 at the first call in this loop, then the
		// return value of this function will be NULL.
		// Add (-1) at the end of generation sequence to signal it is stopped
		if ((newptr = (AGEPTR) malloc(size)) != NULL) {
			newptr->year = -1.0F;
			newptr->next = NULL;
			preptr->next = newptr;
			preptr = newptr;
		};
		if (istop > 1) count++;
	};
	*nSeq = count;
	*nPlan = 2*plan1+plan2;	// so that if nPlan = 1: plan II only,
			// nPlan =2: plan I only, nPlan = 3: both.
			// We use binary to determine from this return value *nPlan,
			// how many plans to be run
	return head;
}

// --------------------------------------------------------------------------

void GetGeneration (AGEPTR *ageSeq, float timeline[], int *nGeneration,
					int *census)
{
// ageSeq is a series of generation batches; the first one in each batch
// is for census size.
// Each batch represents generations of samples, the first one reserved for
// census size (0 if plan II), and ends by (-1) value
// This function reads one batch of generation, store in variable census and
// array timeline (starting from index 0 as normal) for generations
	int i;
	float yr;
//	char more = 1;
	AGEPTR temp, ptr;
	AGEPTR ageptr = *ageSeq;
	if (ageptr == NULL) return;
	// First, count census + how many generations, using variable i:
	i = 0;
	while (ageptr != NULL) {
		yr = (ageptr->year);
		ageptr = ageptr->next;
		if (yr == -1) break;	// end of a batch of generations.
		i++;
	};
	if (i > 2) { // to be valid, should have at least two generations.
	// variable i = (number of generations) + 1 (for the census pop.size).
	// Only valid generations are counted, otherwise, previous timeline
	// should be unchanged (so that the rest follows this last valid one)
		*nGeneration = i-1;
		i = 0;
		ageptr = *ageSeq;
		while (ageptr != NULL) {
			yr = (ageptr->year);
			ageptr = ageptr->next;
			if (yr > -1) {
				if (i == 0) *census = (int) yr;
				else timeline[i-1] = yr;
			} else break;	// break when at the end of a batch of generations.
			i++;
		};
	};
	// Now, free all nodes that have been traversed, i.e., from top down
	// preceding the node that has not been touched, which is ageptr.
	ptr = *ageSeq;
	while (ptr != ageptr) {
		temp = ptr;
		ptr = ptr->next;
		temp->next = NULL;
		free(temp);
	};
	*ageSeq = ageptr;
}

// -------------------------------------------------------------------------
// add this in Apr 2015:
FILE *GetInp (char *inpFolder, char *inpName) {
	char *inpFile;
	FILE *input;
//	printf("File for chromosomes/loci: %s\n", inpName);
	inpFile = (char*) malloc(sizeof(char)*(PATHFILE));
	*inpFile = '\0';
// statement: inpFile = strcat(inpFolder, inpName) will put pointer
// inpFile to the same as inpFolder, that may cause memory reserved
// for inpFile inaccessible. In fact, if use that statement, inpFolder
// string is extended to include string inpName, and a call to free
// both inpFolder and inpFile will cause the program to crash.
	inpFile = strcat(inpFile, inpFolder);
	inpFile = strcat(inpFile, inpName);
	input = fopen(inpFile, "r");
	free (inpFile);
	return input;

}

// -------------------------------------------------------------------------

// --------------------------------------------------------------------------

FILE *InfoDirective (char *mLD, char *mHet, char *mNomura, char *mTemporal,
				char *infoName, char *format, int *nCrit, float critVal[],
				char *matingMod, char *inpFolder, char *inpName,
				char *outFolder, char *outName, int *nPop, int *nloci,
				int *maxMobilVal, int *lenM, FILE *infofile, char *append,
				AGEPTR *ageSeq, int *nSeq, int *tempClue, int *nPlan)
{

	FILE *input = NULL, *output = NULL;
	int i, c, f, m, n, len;
	char *outFile, *prefix;
	int line = 0;
	*nSeq = 0;
	*matingMod = 0;
	*nPop = MAX_POP;	// default value, will be reassigned if FSTAT format

	// read methods to run, represented by an integer m.
	// The algorithm to determine the methods from the value of m is
	// based on base 2 (binary). Suppose m = 9. Then in binary expansion,
	// 9 = 01001
	// 1:  LD (Linkage Disequilibrium) method, position 1
	// 2:  Heterozygote method, position 2
	// 4:  Nomura method position 3
	// 8:  Temporal method, position 4
	// Any sum of those will include the methods represented by the terms.
	// For example,
	// 3 = 1+2: both LD and Het methods,
	// 5 = 1+4: both LD and Nomura methods.

	*append = 0;	// output is assumed to be overwritten
	if (infofile == NULL) return NULL;
// line 1:
	if (FindMethod (infofile, infoName, &line, mLD, mHet, mNomura,
		mTemporal, tempClue) == -1) {
		fclose (infofile);
		return NULL;
	};

// line 2:
// folder containing input, can be empty (= current folder, the folder
// where the interface executable is running)
	line++;
	len = GetToken (infofile, inpFolder, LENDIR, BLANKS, ENDCHRS, &c, &n);
	if (len > 0 || c == SPECHR)	// still on the line containing inpFolder.
		for (; (c=fgetc(infofile)) != EOF && c!='\n';)
			if (c == EOF) {
				fclose (infofile);
				ErrMsg (infoName, "END OF FILE too soon", line);
				return NULL;
			};
	// input file name, can NOT be empty.
// line 3:
	line++;
	if (GetToken (infofile, inpName, LENFILE, BLANKS, ENDCHRS, &c, &n) <= 0)
	{
		fclose (infofile);
		ErrMsg (infoName, "Fail to obtain input file name", line);
		return NULL;
	}
// now see if can open input file:
	if ((input = GetInp(inpFolder, inpName)) == NULL) {
		printf ("Input file [%s] not found in directory %s\n",
				inpName, inpFolder);
		fclose (infofile);
		return NULL;
	}

/* // block out this, use function GetInp inserted in Apr 2015
	inpFile = (char*) malloc(sizeof(char)*(PATHFILE));
	*inpFile = '\0';
	inpFile = strcat(inpFile, inpFolder);
	inpFile = strcat(inpFile, inpName);
// statement: inpFile = strcat(inpFolder, inpName) will put pointer
// inpFile to the same as inpFolder, that may cause memory reserved
// for inpFile inaccessible. In fact, if use that statement, inpFolder
// string is extended to include string inpName, and a call to free
// both inpFolder and inpFile will cause the program to crash.
	if ((input = fopen(inpFile, "r")) == NULL) {
		printf ("Input file [%s] not found\n", inpFile);
		fclose (infofile);
		free (inpFolder);
		free (inpFile);
		return NULL;
	};
	// no longer need the names:
	free (inpFile);
//*/
	// continue reading this infoName file file to obtain info:
	for (; (c=fgetc(infofile)) != EOF && c!='\n';)
		if (c == EOF) {
			fclose (infofile);
			ErrMsg (infoName, "END OF FILE too soon", line);
			return NULL;
		};
// line 4:
	line++;
	f = 0;
// read input format indicator
	if (GetInt (infofile, &f, 1) <= 0) {
		fclose (infofile);
		ErrMsg (infoName, "No format indicator for input file given", line);
		return NULL;
	};
	if (f < MINFORM || f > MAXFORM) {
		fclose (infofile);
		ErrMsg (infoName, "Illegal format indicator for input file", line);
		return NULL;
	};
	*format = f;

// line 5:
	line++;
// folder containing output, can be empty (= current folder)
// n = # of trailing blanks (in BLANKS)  before reaching char in ENDCHRS
	len = GetToken (infofile, outFolder, LENDIR, BLANKS, ENDCHRS, &c, &n);
	if (len > 0 || c == SPECHR)
		for (; (c=fgetc(infofile)) != EOF && c!='\n';)
			if (c == EOF) {
				fclose (infofile);
				ErrMsg (infoName, "END OF FILE too soon", line);
				return NULL;
			};
// line 6:
	line++;
// output file name, can NOT be empty.
// n = # of trailing blanks (in BLANKS)  before reaching char in ENDCHRS
	if (GetToken (infofile, outName, LENFILE, BLANKS, ENDCHRS, &c, &n) <= 0)
	{
		fclose (infofile);
		ErrMsg (infoName, "No OUTPUT file name", line);
		return NULL;
	};
// append output file when no trailing blank between the last character of
// output name to SPECHR (asterisk), i.e., SPECHR immediately after the name
	if (n == 0 && c == SPECHR) *append = 1;
	for (; (c=fgetc(infofile)) != EOF && c!='\n';);	// finish the line
// that's enough, consider the rest for critical values and for max samples
// per pop as optional, so no error messages are given. If anything wrong,
// default values take place.

// lines 7, and possibly 8:
// number of critical values and the values
	line++; line++;
	if ((n = CritValRead (infofile, MAXCRIT, critVal, &c)) <= 0) {
	// the c at the last parameter is no use here
		// n = -1 if error at reading number of crit. value (line 7),
		// otherwise, n is at least 1.
		// The value of "line" now is 8.
		*nCrit = 0;
		fclose (infofile);
		ErrMsg (infoName, "ERROR on Number of Critical Value", line);
		return NULL;
	} else *nCrit = n;
	if (n > 0) line++;
// line 9, or 8 if number of critical values is 0.

// mating model: 0 for random, 1 for monogamy
	m = 0;
	if (GetInt (infofile, &m, 1) <= 0) {
		fclose (infofile);
		ErrMsg (infoName, "ERROR on Entry for Mating Model", line);
		return NULL;
	};
	line++;
	*matingMod = (m != 0)? 1: 0;

// ----------------------------------------------------------------------
// now to read generations for temporal method:
	if (*mTemporal == 1) {
		*ageSeq = AgeSeq (infofile, nSeq, nPlan);
		if (*ageSeq == NULL || (*ageSeq)->year == -1 || *nPlan == 0) {
			printf ("No temporal method: invalid generation set!\n");
			*mTemporal = 0;
		};
	};
// test
/*
	AGEPTR ageptr;
	printf("Number of sequences = %d\n", len);
	ageptr = *ageSeq;
	while (ageptr != NULL) {
		if ((ageptr->year) > -1) {
			printf ("%6.2f", ageptr->year);
		} else printf ("%6.1f\n", ageptr->year);
		ageptr = ageptr->next;
	};
*/
// ----------------------------------------------------------------------
// Now need to look into input file whose name was given by this infofile
// to see if format is alright, and also check if output can be opened.
// First, read input file to obtain necessary parameter values
	if (f==FSTAT) {	// FSTAT format
		if (GetInfoDat (input, nPop, nloci, maxMobilVal, lenM, LEN_BLOCK) == 0)
		{
			printf ("Error in (FSTAT format) input file \"%s\"\n", inpName);
			return NULL;
		};
	};
	if (f==GENPOP) {	// GENEPOP format
		if ((*nloci = GetnLoci (input, LEN_BLOCK, lenM)) <= 0) {
			printf ("Error in (GENEPOP format) input file \"%s\"\n", inpName);
			return NULL;
		} else {
			rewind (input);
			for (; (c=fgetc(input)) !='\n';);	// position at first locus
		};
		printf ("Number of loci = %d, %d-digit alleles\n", *nloci, *lenM);
	// should assign maxMobilVal to make sure it doesn't get garbage.
	// Since the length is lenM, assign max possible value:
		for (*maxMobilVal=1, i=1; i<=*lenM; i++) *maxMobilVal *= 10;
	};
	// also make sure output file can be opened later.
	outFile = (char*) malloc(sizeof(char)*(PATHFILE));
	*outFile = '\0';
	// concatenate outFolder and outName to get outFile.
	// Since outFolder is a return value, don't want it to change by
	// the call of function strcat and also by the same reason given
	// at concatenating inpFolder with inpName; so, assign
	// outFile = outFolder first:
	outFile = strcat (outFile, outFolder);
	prefix = (char*) malloc(sizeof(char)*LENFILE);
	*prefix = '\0';
	if (*outName != '\0') outFile = strcat (outFile, outName);
	else {
		GetPrefix (inpName, prefix, LENFILE-6, PATHCHR);
		outName = strcat (outName, strcat(prefix, "Ne"));
		outName = strcat (outName, EXTENSION);
		outFile = strcat (outFile, outName);
	};
	// open but not destroy its content if the file exists, use "a"
	if ((output=fopen(outFile, "a")) == NULL) {
		fclose (input);
		fclose (output);
		free (outFile);
		free (prefix);
		printf ("Cannot open file \"%s\" for output.\n", outFile);
		return NULL;
	} else fclose (output);	// done, output file is OK to open
	free (outFile);
	free (prefix);

	return input;
}

// -------------------------------------------------------------------------
// April 2015
// This function is to read the last line of control file for OptDirective,
// to obtain the necesary option for dealing with chromosomes.
// The last line of the control file consists of two entries:
// The first one is a number, chroGrp, either 1 or 2, to indicate if pairs
// of loci are to be taken within a chromosome (chroGrp = 1), or
// must be in different chromosomes (chroGrp = 2).
// The second entry is the file name that lists pairs (chromosome, locus)
// on each line. This file is supposed to be in the same folder as the
// genotype input file (inpFolder).
// The return value of the function is the file being opened, NULL if failed.
// -------------------------------------------------------------------------
int ChroInfo (FILE *infoFile, char *inpName) {
	int p, c, n;
// get first string as an integer, then stay on the line (last parameter = 0)
// p is the chroGrp for using locus pairs in the same or distinct chromosomes
	if (GetInt (infoFile, &p, 0) <= 0) return 0;
	if (p != 1 && p != 2) return 0;
	*inpName = '\0';
	if (GetToken (infoFile, inpName, LENFILE, BLANKS, ENDCHRS, &c, &n) <= 0)
	{
		return 0;
	}
	return p;

}

// -------------------------------------------------------------------------
// April 2015
// This function is to read the last line of control file for OptDirective,
// to obtain the necesary option for dealing with chromosomes.
// The last line of the control file consists of two entries:
// The first one is a number, chroGrp, either 1 or 2, to indicate if pairs
// of loci are to be taken within a chromosome (chroGrp = 1), or
// must be in different chromosomes (chroGrp = 2).
// The second entry is the file name that lists pairs (chromosome, locus)
// on each line. This file is supposed to be in the same folder as the
// genotype input file (inpFolder).
// The return value of the function is the file being opened, NULL if failed.
// -------------------------------------------------------------------------
FILE *ChromoInp (FILE *infoFile, char *inpFolder, int *chroGrp) {
	char *inpFile;
	char *inpName;
	FILE *input;
	int p, c, n;
	*chroGrp = 0;
// get first string as an integer, then stay on the line (last parameter = 0)
// p is the chroGrp for using locus pairs in the same or distinct chromosomes
	if (GetInt (infoFile, &p, 0) <= 0) return NULL;
	if (p != 1 && p != 2) return NULL;
	inpName = (char*) malloc(sizeof(char)*LENFILE);
	*inpName = '\0';
	if (GetToken (infoFile, inpName, LENFILE, BLANKS, ENDCHRS, &c, &n) <= 0)
	{
		free (inpName);
		return NULL;
	}
	printf("File for chromosomes/loci: %s\n", inpName);
	inpFile = (char*) malloc(sizeof(char)*(PATHFILE));
	*inpFile = '\0';
	inpFile = strcat(inpFile, inpFolder);
	inpFile = strcat(inpFile, inpName);
	input = fopen(inpFile, "r");
	// no longer need the names:
	free (inpName);
	free (inpFile);
	*chroGrp = p;
	return input;

}


// --------------------------------------------------------------------------

int OptDirective (char *optName, char *xOutLD, char *xOutHet, char *xOutCoan,
				char *xOutTemp, int *maxSamp,
				int *minPop, int *maxPop, int nPop,
				int *popLoc1, int *popLoc2,
				int *popBurr1, int *popBurr2, int *topBCrit, char *misDat,
				char *param, char *jacknife,
				int nloci, char *locUse, int *nLocDel, int *tempxClue,
				char *byRange, int *topXCrit, char *tabX,
// add Jan 2015/ Apr 2015:
//				char *inpFolder, FILE *chroInp, int *chroGrp,
				char *inpFolder, char *chrofileName, int *chroGrp,
				char *sepBurOut, char *moreCol, char *BurAlePair)

// popLoc1, popLoc2, popBurr1, popBurr2 for the ranges of populations to have
// Freq. Data output (to outLoc file) and Burrows coeffs (for outBurr file).
// popLoc1 < 0 indicates all pop to have freq output, which will be
// interpreted in RunOption, regardless of the value of popLoc2 (if entered).
// If there are two popLoc1, popLoc2 entered, it will indicate a range
// of populations from popLoc1 to popLoc2 to have freq output.
// If there is only one value nonnegative popLoc1 entered, it will indicate
// populations up to popLoc1 to have freq output (so if = 0, no freq.output).
// Same rule for Burrows.

// add topCrit for limiting to the top critical values in x-tra files.

// Return value: indicate how many data lines have been read from file
// with name optName, while collecting parameters from this file.
{
// add d, r in 2015 for division of 2 integers, to get sepBurOut, BurAlePair
	int d, i, r;
	int c, m, n;
	int linedone = -1;
	int *xClues;
	FILE *optFile;
	if ((optFile = fopen (optName, "r")) == NULL) return -1;

// Added in Mar/Apr 2015
	char *chroFileName = (char*) malloc(sizeof(char)*LENFILE);
	*chroFileName = '\0';
	*chroGrp = 0;

	linedone++;
// Default values:
	*popLoc1 = 0;
	*popLoc2 = 0;
	*popBurr1 = 0;
	*popBurr2 = 0;
	*topBCrit = MAXCRIT;
	*misDat = 1;
	*param = 1;
	*jacknife = 1;
	*maxSamp = 0;
	*minPop = 1;
	*maxPop = 0;
	*nLocDel = 0;
// add jan 2015:
	*sepBurOut = 0;	// no separating Burrows outputs
	*moreCol = 0;	// =0: less columns for summarizing Burrows table
					// =1: more columns
	*BurAlePair = 0;
// If a number of numeric data to be read in input file can be varied,
// then scanf may direct cursor to the next line if its list contains more
// than the ones in input, unless some non-number is on the way to stop it.
// For example, if a maximum of 3 numbers are to be read on a line of input
// file where some are optional, then there can be less than 3 entered in
// input (the ones not entered will have default values).
// Therefore, to prevent going to next line which is reserved for data of
// different info, we use function GetClues instead
// -------------------------------------------------------------------
// extra output and its option: (modified July 2013 for Tab-delimiter)
	xClues = (int*) malloc(sizeof(int)*4);
// default values for methods to have extra outputs, which temporal,
// and which critical values:
	*xClues = 0; *(xClues+1) = 0; *(xClues+2) = MAXCRIT; *(xClues+3) = TABX;
	m = GetClues(optFile, xClues, 4, 1);
	n = *xClues; *tempxClue = *(xClues+1); *topXCrit = *(xClues+2);
	// July 2013:
	*tabX = (*(xClues+3) == 0)? 0: 1;
	free(xClues);
	SetMethod (n, xOutLD, xOutHet, xOutCoan, xOutTemp);
	if (m <= 0) {	// no number read!
		fclose (optFile);	// no more reading, the rest: default values
		return linedone;
	};
	linedone++;

// max samples per pop:
	if (GetInt (optFile, maxSamp, 1) <= 0) {
		fclose (optFile);
		return linedone;
	};
	linedone++;

// # pops to have locus data output:
	m = GetPairI(optFile, popLoc1, popLoc2, 1);
// if popLoc1 <0, all pops will have freq. output, which will be adjusted
// at the calling function (RunOption). If popLoc1 = 0, popLoc2 is ignored,
// no Freq. output!
	if (m <= 0) {
		fclose (optFile);
		return linedone;
	};
	if (m == 1 && *popLoc1 > 0) {	// only one number on this line,
	// then Freq. output from population 1 to that entry (popLoc1).
		*popLoc2 = *popLoc1;	// Reassign popLoc1 and popLoc2 so that
		*popLoc1 = 1;			// populations run from popLoc1 to popLoc2.
	};
	// If the above condition is false, either there is a second number
	// reserved for popLoc2 or popLoc1 < 0.
	// If there is popLoc2 read, and if popLoc2 < popLoc1: no Freq output.
	// If popLoc1 < 0, no matter popLoc2 is read or not, the calling
	// function RunOption will make all pops to have freq output.
	linedone++;

// # pops to have Burrows coefs. output: could have up to 3 entries here.
// Get the first entry (which can be negative for unrestriction),
// then for the next two, use GetPair.
// The last 0 in GetInt is to stay on the line if a legitimate integer read
	if (GetInt (optFile, popBurr1, 0) <= 0) {
		fclose (optFile);
		return linedone;
	};
// change in Jan 2015: to be able to read entry on separating Burrows files
//	m = GetPair(infofile, popBurr2, topBCrit, 1);
	m = GetPair(optFile, popBurr2, topBCrit, 0);
	if (*popBurr1 > 0) {
		if (m == 0) {	// 1 entry only, no new values for popBurr2, topBCrit
			*popBurr2 = *popBurr1;
			*popBurr1 = 1;
// add Jan 2015 for this if (m == 2):
		} else if (m == 2) {	// have both popBurr2, topBCrit on the line,
		// now look for clue on separating Burrows output
			if ((n = GetPair (optFile, &i, &c, 1)) > 0) {
				d = i/2;
				r = i%2;	// remainder when dividing i by 2, so r = 0 or 1.
				*sepBurOut = (r != 0)? 1: 0;
				*BurAlePair = (d != 0)? 1: 0;
				if (n > 1) *moreCol = (c != 0)? 1: 0;
			}
		}
	} else if (*popBurr1 < 0) {	// no limit on population, the 2nd number is
// added Jan 2015 for condition (m>1): last entry for separating Burrows outputs
		if (m > 1) {
			i = *topBCrit;	// value read as topCrit is actually
			d = i/2;		// for sepBurOut, Buralepair
			r = i%2;
			*sepBurOut = (r != 0)? 1: 0;
			*BurAlePair = (d != 0)? 1: 0;
			if (GetInt (optFile, &c, 1) > 0) *moreCol = (c != 0)? 1: 0;
		}
		if (m > 0) *topBCrit = *popBurr2;	// the number of highest Pcrits
	}
	linedone++;
// run parametric CIs or not
	if (GetInt (optFile, &n, 1) <= 0) {
		fclose (optFile);
		return linedone;
	};
	linedone++;
	*param = (n != 0)? 1: 0;
// run jacknife CIs or not
	if (GetInt (optFile, &n, 1) <= 0) {
		fclose (optFile);
		return linedone;
	};
	linedone++;
	*jacknife = (n != 0)? 1: 0;
// ----------------------------------------------------------------------
// This line is for population range
	m = 0;
	c = GetPair(optFile, &m, &n, 1);	// c = 0: no number or negative for
	if (m > nPop) m = nPop;			// the first entry, then m still = 0
	if (m > 0) {	// No declaration of error if c=0. If m > 0, then c >=1
		*maxPop = m;
		if (c == 2) {	// then n is obtained a value from GetPair
			if (n >= m) {
				*minPop = m;	// m and n, so the range is from m to n.
				*maxPop = n;
			};	// else, n < m: the second entry is erroneous, maxPop is not
				// reassigned, so the range is from minPop=1 to maxPop=m
		};	// else, only 1 number: population range is from *minPop=1 to m
	} else	// if the first entry is <= 0, assuming all populations are used
		*maxPop = nPop;


	linedone++;
	if ((n = LociDropped (optFile,locUse,nloci,&linedone,0,0,byRange)) == -1)
	{
		fclose (optFile);
		return linedone;
	};
	*nLocDel = n;
	linedone++;
// Add in 2015: have missing data file or not:
	if (GetInt (optFile, &n, 1) <= 0) return linedone;
	linedone++;
	*misDat = (n != 0)? 1: 0;

// Add in Apr 2015 for dealing with chromosomes:
	*chroGrp = ChroInfo (optFile, chrofileName);
// FILE should not be a return value as a parameter!
//	chroInp = ChromoInp (optFile, inpFolder, chroGrp);

// Close the file:
	fclose (optFile);
	return linedone;	// maximum data lines have been read
}

// -------------------------End of GetData.c    -----------------------------

// --------------------------------------------------------------------------
// PutGene:

//------------------------------------------------------------------
/*
 When reading input file of genotypes, assumed to have "nloci" loci,
 the program will store data of a population in 2 arrays of pointers,
 each has size nloci.
 Suppose the population has N samples.
	* Each element of the first array is a pointer to the list of genotypes
	of all samples at the same locus.
	* Each element of the second array points to a list of mobilities.

*/
//------------------------------------------------------------------
FISHPTR MakeFish (int allele[])
{
	FISHPTR newptr;

	int size = sizeof(struct fish);
	if ((newptr = (FISHPTR) malloc(size)) != NULL)
	{
		newptr->gene[0] = allele[0];
		newptr->gene[1] = allele[1];
		newptr->next = NULL;
	}
	return newptr;
};


//------------------------------------------------------------------

char AddFishWide(FISHPTR *fishHead, FISHPTR *fishTail, int nloci,
				int sample[], char *locUse, char quit)
{

// add one sample to each of nloci list of genotypes
// (each list corresponds to one locus).
// On return, each of fishHead points to the top of the list, fishTail
// points to the last.
// The order of elements in the lists is the order of samples read,
// so always add new sample to the last one, using fishTail to locate.
	int allele[2];
	int p;
	FISHPTR newfish;
	if (quit == 0) return 1;
	for (p=0; p<nloci; p++) {	// only accept samples with full data
		if (sample[2*p] <= 0 || sample[2*p+1] <=0) {
			allele[0] = 0;
			allele[1] = 0;
		} else {
			allele[0] = sample[2*p];
			allele[1] = sample[2*p+1];
		};
		if ((newfish = MakeFish (allele)) == NULL) return 0;
		if (*(fishHead+p) == NULL) {
			*(fishHead+p) = newfish;
		} else (*(fishTail+p))->next = newfish;
		*(fishTail+p) = newfish;
	}
	return 1;

}


//------------------------------------------------------------------
void RemoveFish(FISHPTR *fishList, int nloci)
{

	int p;
	for (p=0; p<nloci; p++) {
		FISHPTR curr = *(fishList+p), temp;
		for (; curr != NULL; curr = temp) {
			temp = curr->next;
			curr->next = NULL;
			free (curr);
		};
	}
// Since fishList is allocated at the beginning and only reinitialized
// after each population, this should be commented out.
//	free(fishList);

}

//------------------------------------------------------------------



ALLEPTR MakeAlle (int mvalue)
// This create a new allele node, return NULL if failed,
{
	ALLEPTR newptr;
	if ((newptr = (ALLEPTR) malloc(sizeof(struct allele))) != NULL)
	{
		newptr->mValue = mvalue;
		newptr->freq = 0;
		newptr->hetx = 0;
		newptr->copy = 1;
		newptr->homozyg = 0;
		newptr->next = NULL;
	}
	return newptr;
}

//------------------------------------------------------------------

ALLEPTR AddAlle (int p, ALLEPTR head, int alleleK, ALLEPTR *ptrPos,
				 int *nMobil, int *errcode)
{
// (counter part: AddMobil)
// This routine adds one allele to the list of alleles at one locus.
// Return the pointer to the head of the list,
// ptrPos points to where the new one is added (if this mobility alleleK
// is new) or points to the current one having value = alleleK.
// This is needed since we want to know where it is added, so that
// homozygosity (determined when we have both alleles) can be entered.
// The list of alleles is in ascending order by its mValue's.
// nMobil is incremented by 1 if new allele mobility is added.

    ALLEPTR ptr1, ptr2, newnode;


	if (head == NULL) {
		if ((newnode = MakeAlle (alleleK)) == NULL) {
			printf  ("Out of memory for adding allele at locus %d!\n", p+1);
			*errcode = -1;
			return head;
//			exit (EXIT_FAILURE);
		};
		head = newnode;
		*ptrPos = head;
		(*nMobil)++;
	} else {
		if (alleleK < head->mValue) {	// add to front
			if ((newnode = MakeAlle (alleleK)) == NULL) {
				printf ("Out of memory for adding allele at locus %d!\n", p+1);
				*errcode = -1;
				return head;
//				exit (EXIT_FAILURE);
			};
			newnode->next = head;
			head = newnode;
			*ptrPos = head;
			(*nMobil)++;
		} else {
			ptr1 = head;
			ptr2 = ptr1->next;
// search for the last node ptr1 such that: ptr1->mValue <= alleleK
			while (ptr2 != NULL) {
				if (ptr2->mValue <= alleleK) {
					ptr1 = ptr2;	// then ptr1->mValue <= alleleK
					ptr2 = ptr2->next;
				} else {
					break;
				}
			}	// end of search

// now, ptr1->mValue <= alleleK, either ptr1 is at the end of the list,
// or alleleK is between mValues of ptr1 and ptr2.
			if (ptr1->mValue == alleleK) {
				(ptr1->copy)++;
				*ptrPos = ptr1;
			} else {
				if ((newnode = MakeAlle (alleleK)) == NULL) {
					printf ("Out of memory for adding allele at locus %d!\n", p+1);
					*errcode = -1;
					return head;
//					exit (EXIT_FAILURE);
				};
				newnode->next = ptr2;
				ptr1->next = newnode;
				*ptrPos = newnode;
				(*nMobil)++;
			}
		}	// end of "if (alleleK < head->mValue) ... else "
	}
	return head;
}

//------------------------------------------------------------------

ALLEPTR AddGeno (int p, ALLEPTR head, int gene[], ALLEPTR *ptrPos1,
				ALLEPTR *ptrPos2, int *nMobil, int *missptr,
				int maxMobilVal, int *errcode)
{
// (counter part: AddGene)
// This routine calls AddAlle twice to add genotype to the mobility list
// at a locus. Only add both alleles or none. If one is zero or exceeds
// maxMobilVal, considered as missing data
	ALLEPTR temp;
	*errcode = 0;
	if ((gene[0]<=0) || (gene[1]<=0) ||
		(gene[0]>maxMobilVal) || (gene[1]>maxMobilVal)) {
		(*missptr)++;
		return head;
	};
	temp = AddAlle (p, head, gene[0], ptrPos1, nMobil, errcode);
	if (*errcode != 0) return head;	// no change
	temp = AddAlle (p, temp, gene[1], ptrPos2, nMobil, errcode);
	if (*errcode != 0) return head;	// no change
	if (gene[0] == gene[1])
// ptrPos1, ptrPos2 point to the same address: change one will change both.
		((*ptrPos2)->homozyg)++;

	return temp;
}

//------------------------------------------------------------------


int AddAlleWide (ALLEPTR *alleList, int nloci, int sample[],
				   int *nMobil, int *missptr, int maxMobilVal,
				   int popRead, int samp)
{
// (counter part: AddMobilWide)
// This routine put genotype data for one sample across all loci
// in (nloci) alleList.

	ALLEPTR *ptrPos1, *ptrPos2;
	ALLEPTR head;
	int p;
	int gene[2];
	int errcode = 0;
	if ((ptrPos1= (ALLEPTR*) malloc(sizeof(ALLEPTR)*nloci)) == NULL) {
		errcode = -1;
	};
	if ((ptrPos2 = (ALLEPTR*) malloc(sizeof(ALLEPTR)*nloci)) == NULL) {
		free (ptrPos1);
		errcode = -1;
	};
	if (errcode == -1) {
		printf ("Out of memory at population %d, sample %d when adding allele!\n",
				popRead, samp);
		return errcode;
	};

	for (p=0; p<nloci; p++) {
		*(ptrPos1+p) = NULL; *(ptrPos2+p) = NULL;
		gene[0] = sample[2*p];
		gene[1] = sample[2*p+1];

		// *(missptr+p) represents missing data at locus (p+1)
		head = *(alleList+p);
		*(alleList+p) = AddGeno (p, head, gene, (ptrPos1+p), (ptrPos2+p),
								nMobil+p, missptr+p,
								maxMobilVal, &errcode);
		if (errcode != 0) break;

	}	// end of "for (p=0; p<nloci; p++)"

// Deallocate memory for local variables ptrPos1, ptrPos2, missing:
	free (ptrPos1); free (ptrPos2);
	return errcode;
}

//------------------------------------------------------------------


void RemoveAlle(ALLEPTR *alleList, int nloci)
{

	int p;
	for (p=0; p<nloci; p++) {
		ALLEPTR curr = *(alleList+p), temp;
		for (; curr != NULL; curr = temp) {
			temp = curr->next;
			free (curr);
		};
    }
// Since alleList is allocated at the beginning and only reinitialized
// after each population, this should be commented out.
//	free(alleList);

}

//-------------------------   End of PutGene.c   -----------------------


//----------------------------------------------------------------------

// CalFreq
//----------------------------------------------------------------------

/*
The calculations here on the Linkage Disequilibrium method
(Burrows coefficients, pairwise analysis of loci and the estimate of the
effective population based on coefficient R2) are based on the paper:

	A bias correction for estimates of effective population size
	based on linkage disequilibrium at unlinked gene loci

	by Robin Waples, Conservation Genetics (2006) 7:167-184

The calculation of effective breeders by Heterozygote Excess is based on:

	On the Potential for Estimating the Effective Number of Breeders
	From Heterozygo-Excess in Progeny
	by Pudovkin et al, Genetics 144:383-387 (1996)

	Nb_HetEx: A program to Estimate the Effective Number of Breeders
	by Zhdanova & Pudovkin, J. of Heredity 2008:99(6).

*/
// ------------------------------------------------------------------------


// this may not be used
int NumLocTake (char *locUse, int nloci)
{
	int p, q;
	for (p=0, q=0; p<nloci; p++)
		if (*(locUse+p) == 1) q++;
	return q;
}

//-------------------------------------------------------------------------
// Janked from Print
// --------------------------------------------------------------------------

void PrtLines (FILE *output, int ndash, char dash)
{
	int n;
	if (output == NULL) return;
	for (n=0; n<ndash; n++) fprintf (output, "%c", dash);
	fprintf(output, "\n");
	fflush (output);

}
// --------------------------------------------------------------------------
// called by Loci_Eligible
// add Jan 2015: sepBurOut, moreCol
// Modified in Dec 2016:
// Modify headline for spec cutoff (not using "Lowest Freq...")
void WriteLoci (FILE *outFile, int nloci, char *okLoc, float cutoff,
				int locChk, int *nMobil, int *nInd, int nLocOK,
				char more, char burr, int mLoc, char sepBurOut, char moreCol)
// nLocOK is the number of loci that satisfies the "cutoff" requirement,
// which is the true number of loci that will be used, given as parameter.
// that is, loci where all frequencies are > cutoff.
// added Sept 2011: mLoc, for maximum loci printed
// (burr = 1 when this is for printing to Burrow file, else for Freq. File)
{
	int i, p, m;
	long totInd = 0, totAlle = 0;
	if (outFile == NULL || (more == 0)) return;
	// change in 2015 for Burrows files:
	if (burr == 1 && sepBurOut == 1 && NOEXPLAIN == 1) return;
//	for (p=0, q=0; p<nloci; p++) {
//		if (*(okLoc+p) == 0) continue;
//		q++;
//	};	// q becomes the number of loci being considered,
	if (burr == 1) {	// outFile is Burrows file, either sepBurOut = 0 or
// change in Nov 2014: 	// NOEXPLAIN = 0, or both
		m = 68;
		if (moreCol == 1) m += 43;
		for (i=0; i<m; i++) fprintf (outFile, "*");
		fprintf (outFile, "\n");
	}
// Dec 2016
	if (cutoff > 0 && cutoff <= PCRITX)
		fprintf (outFile, "\nExclude Singleton Alleles at each locus\n");
	else
		fprintf (outFile, "\nLowest Allele Frequency Used =%9.5f\n", cutoff);
		//                   123456789012345678901234567890123456789
//	for (i=0; i<39; i++) fprintf (outFile, "="); fprintf (outFile, "\n");
	PrtLines (outFile, 39, '=');

// print # alleles and # ind. of alleles at each considered locus.
	if (burr != 1) {
//		fprintf (outFile, "\n   Locus   N.Alleles   N. Ind\n");
		fprintf (outFile, "   Locus   N.Alleles   N. Ind\n");
		for (p=0, m=0; p<nloci; p++) {
			if (*(okLoc+p) == 0) continue;
			m++;
			totAlle += *(nMobil+p);
			totInd += *(nInd+p);
//			if (m > mLoc) break;
			if (m <= mLoc)
			fprintf (outFile, "%7d%10d%11d\n", p+1, *(nMobil+p), *(nInd+p));
		};
		if (nLocOK > mLoc)
			fprintf (outFile, "(Only the first %d loci are printed)\n", mLoc);
//		else fprintf (outFile, "\n");
// add SUM in Aug 2012
//		for (i=0; i<39; i++) fprintf (outFile, "-");
//		fprintf (outFile, "\n   SUM:%10lu%11lu\n", totAlle, totInd);
//		for (i=0; i<39; i++) fprintf (outFile, "-");fprintf (outFile, "\n");
		PrtLines (outFile, 39, '-');
		fprintf (outFile, "   SUM:%10lu%11lu\n", totAlle, totInd);
		PrtLines (outFile, 39, '-');
/*
	} else {
		fprintf (outFile, "\nLoci being shown:");
		for (p=0, m=0; p<nloci; p++) {
			if (*(okLoc+p) == 0) continue;
			m++;
			if (m > mLoc) break;
// print up to 12 loci per line
			if (m>1 && m%12==1) fprintf (outFile, "\n%17s", " ");
			fprintf (outFile, "%6d", p+1);
		};
//*/
	};
	// change in Nov 2014
	if (burr != 1 || NOEXPLAIN != 1) {
		fprintf (outFile, "Total number of loci to be used =%6d\n", nLocOK);
		fprintf (outFile, "#Loci rejected by required freq.=%6d\n", locChk - nLocOK);
	}
	fflush (outFile);
}
// --------------------------------------------------------------------------

//
// The next 7 having the word "Het" belong to HetExcess module
// -----------------------------------------------------------
float HetNb (float d)
// calculate effective breeder based on Pudovkin paper, formula (4),
// which is Nb = 1/(2D) + 1/2(D+1) = (2D+1)/(2D(D+1)).
// This function is always decreasing, going from INF to -INF in (-1, 0)
// and from INF to 0 in (0, INF); it is zero at -0.5 (so it is positive
// for d=-1 to d=-0.5, and negative for d=-0.5 to d=0).
{
	float d1, d2;
// add this line in Aug 2012, i.e., no longer having negative value for Nb
	if (d <= EPSILON) return INFINITE;
	d2 = d+d;
	d1 = d+1;
	return (d2+1)/(d2*d1);
}
// --------------------------------------------------------------------------
void HetAverage (int locPoly, long mTotal, float *dSum, float *dXSum,
				 float *dXWSum, float *NbUw, float *NbH, float *NbHW,
				 float dWtotal, float *NeWt, float *hSamp,
// add Jan 2012:
				 long nInd, float totWeigSq, float wSumhetSq, float *stdErr)
// locPoly: number of polymorphic loci
// mTotal: total number of alleles in all polymorphic loci
// dWtotal: total weight in all loci
// The two below are added Jan 2012, used to output Standard Error stdErr
// totWeigSq: sum of weight^2
// wSumhetSq: sum of weighted (het-ex)^2

// On input:
// dSum: sum of all het excesses for all alleles in all loci
// dXSum: sum across loci of all averages of het excesses at each locus
// dXWSum: similar to dXSum, but on weighted averages.
// Then on output, those 3 become average across loci, which then are used
// to calculate  corresponding effective breeders:
// (except Nb from "dSum" which will be calculated in PrtHetSum)
// NbUw: on unweigted average of het.excess across loci,
// NeWt: on weigted average of het.excess across loci, for main output,

// NbH: harmonic mean on Nb across loci
// NbHW: weighted harmonic mean on Nb across loci
{
// add Jan 2012
//	float dMoment = 0;
	float dWmeanSq = 0;
	float ratio1 = 0, ratio2 = 0;
	float r = (float) nInd;
	if (r <= 0) r = (float) mTotal - locPoly;	// this is for the case that
	// a zero is in place of nInd when this function is called. This happens
	// when we first evaluate het-excess at all alleles in Calcul_Freq where
	// the number of independent alleles are not computed.
	// The calling of this function in Calcul_Freq (and its counter part)
	// is mainly to output to auxiliary file for locus data.
	*stdErr = 0;
	if (locPoly > 0) {
		*dXSum /= locPoly;
		// add Mar 2012: hSamp input as the sum of inverses of samp. sizes
		*hSamp = locPoly/(*hSamp);	// becoming harmonic mean of samp. sizes
		*dSum /= mTotal;
	};
	if (fabs(*dXSum) > EPSILON) *NbUw = HetNb(*dXSum);
	else *NbUw = INFINITE;

	if (dWtotal > 0) {
		*dXWSum /= dWtotal;
		dWmeanSq = (*dXWSum) * (*dXWSum);
		wSumhetSq /= dWtotal;
		ratio1 = (wSumhetSq - dWmeanSq)/r;
//		ratio1 = (wSumhetSq - dWmeanSq)/(dWtotal*r);
		if (ratio1 < 0) ratio1 = 0;
		r = (dWtotal)*(dWtotal);	// square of the total weight.
		// watch out: if number of loci and samples are big, r could overflow
		ratio2 = r / (r - totWeigSq);
		*stdErr = (float) sqrt(ratio1 * ratio2);
	};
//	if (fabs(*dXWSum) > EPSILON) *NeWt = HetNb(*dXWSum);
//	else *NeWt = INFINITE;
	*NeWt = HetNb(*dXWSum);

	if (fabs(*NbH) > EPSILON) *NbH = ((float) locPoly)/(*NbH);
	else *NbH = INFINITE;
	if (fabs(*NbHW) > EPSILON) *NbHW = dWtotal/(*NbHW);
	else *NbHW = INFINITE;
}


// --------------------------------------------------------------------------
void HetForeword (FILE *outLoc) {

	fprintf (outLoc,
			"\nSummarizing Table for Effective Breeders by "
			"Heterozygote Excess across loci\n");
	fprintf (outLoc, "#smp  = number of individuals having data,\n#alle = "
			"number of alleles (all low frequency alleles are combined into one)\n");
	fprintf (outLoc,
			"Dm  = Mean of Het-Ex,      "
			"Wt  = Sum[allele weights] = (#alle - 1)*sqrt(#smp)\n");
	fprintf (outLoc,
			"Dm2 = Mean of (Het-Ex)^2,  Wt2 = Sum[(alle. weight)^2] = "
			"(Wt^2)/(#alle)\n");
	fprintf (outLoc, "N[eb] = INF if Dm <= 0.\n");

}

// --------------------------------------------------------------------------
// add in Aug 2012:
void PrtDetailHet (FILE *outLoc, float dWAve, long totInd,
				float totWeight, float wSumhetSq, float totWeigSq)
{
	float r, s, t;
	fprintf (outLoc, "\n* I   = Total of Ind. alleles       = %10lu", totInd);
	fprintf (outLoc, "\n* D   = Weighted mean of Dm         = %10.6f", dWAve);
	t = dWAve*dWAve;
	fprintf (outLoc, "\n* D^2 = (Square of D)               = %10.6f", t);
	wSumhetSq /= totWeight;
	fprintf (outLoc, "\n* D2  = Weighted mean of Dm2        = %10.6f", wSumhetSq);
	r = (wSumhetSq - t)/((float) totInd);
	t = totWeight*totWeight;
	fprintf (outLoc, "\n* W^2 = (Sum of Wt)^2               = %10.2f", t);
	fprintf (outLoc, "\n* W2  = (Sum Wt2)                   = %10.2f", totWeigSq);
	s = t/(t-totWeigSq);
	fprintf (outLoc, "\n* R   = (D2 - D^2)/I                = %10.6f", r);
	fprintf (outLoc, "\n* S   = W^2/(W^2 - W2)              = %10.6f", s);
	t = r*s;
	fprintf (outLoc, "\n* The product P = RxS               = %10.6f", t);
	t = (float) sqrt(t);
	fprintf (outLoc, "\n* Square root of P                  = %10.6f\n", t);
}

// --------------------------------------------------------------------------
void PrtHetSum (FILE *outLoc, float NbH, float NbHW, float NeWt,
				float NbUw, float dSum, float dWAve, float stdErr)
{
	if (NeWt < INFINITE)
		fprintf (outLoc, "\n* Standard Error of D               = %10.6f", stdErr);
	fprintf (outLoc, "\n* N[eb] from weighted mean of Dm    = ");
	if (NeWt < INFINITE) {
		fprintf (outLoc, "%10.1f\n", NeWt);
	} else {
		fprintf (outLoc, "%10s\n", "INFINITE");
	};

	fprintf (outLoc, "\nFor Information only:");
	fprintf (outLoc, "\nCalculated over polymorphic loci:\n");
	fprintf (outLoc, "* Unweighted Harmonic Mean of N[eb] = ");
	if (NbH < INFINITE) fprintf (outLoc, "%10.1f\n", NbH);
	else fprintf (outLoc, "%10s\n", "INFINITE");
	fprintf (outLoc, "* Weighted Harmonic Mean of N[eb]   = ");
	if (NbHW < INFINITE) fprintf (outLoc, "%10.1f\n", NbHW);
	else fprintf (outLoc, "%10s\n", "INFINITE");
	fprintf (outLoc, "\n");
	fprintf (outLoc, "* N[eb] from unweighted mean of Dm  = ");
	if (NeWt < INFINITE) fprintf (outLoc, "%10.1f\n", NbUw);
	else fprintf (outLoc, "%10s\n", "INFINITE");
	fprintf (outLoc, "* N[eb] from mean of all Het-Ex D's = ");
	if (fabs(dSum) >= EPSILON)
		fprintf (outLoc, "%10.1f\n", HetNb(dSum));
	else fprintf (outLoc, "%10s\n", "INFINITE");
	fprintf (outLoc, "\n");
}


// --------------------------------------------------------------------------
float HetExp (float freq, int count)
// write q = freq, frequency of the current allele at current locus.
// *hetExp = (float) 2*q*(1 - q)*(1 + r). Let n = count, number of parents.
// Then 2n is the number of genes in either male or female parents assuming
// there are equal males and females. In Pudovkin's paper (1996),
// formula (1): *hetExp = 2*q*(1-q)*(1+ 1/(2n)),
//             r = 1/2n,
// and in Nb_HetEx: A Program ... , J. of Heredity 2008:99(6),
// *hetExp = 2*q*(1-q)* 2n/(2n-1) = 2*q*(1-q)*(1+ 1/(2n-1)),
//             r = 1/(2n-1).
// This is probably comes from the expected value of the genes descended from
// parents (when males may not equal females) is 2n-1.
// Here, we follow the 2008-paper.
{
	float r, f;
	r = (float) 1.0/((float)2.0*count - 1);
	f = (float) 2*freq*(1 - freq)*(1 + r);
	return f;
}

// --------------------------------------------------------------------------

float HetExAtLoc (int nA, float *hExp, float *hObs, float freq,
				int hetero, int count)
// Returns the het. excess for an allele at a particular locus.
// Input:
// nA is the number of alleles (to be accepted) at current locus.
// freq: frequency of the allele being considered, at current locus.
// hetero: number of heterozygotes having an allele being considered
// count: number of samples having data
// Output: also return values for
// hExp: expected freq. of heterozygote,
// hObs: observed freq. of heterozygotes,
// and update totAlle: total alleles considered so far.
{
//	float r
	float d = 0;
	if (nA > 1) {
		// polymorphic locus
		*hObs = ((float) hetero)/count;
/* // These lines are replaced by function HetExp
		r = (float) 1.0/((float)2.0*count - 1);
		*hExp = (float) 2*freq*(1 - freq)*(1 + r);
//*/
		*hExp = HetExp (freq, count);
// d is the D-value given in Pudovkin, p. 384, column 2.
		d = (*hObs - *hExp)/(*hExp);
		if (fabs(d) < EPSILON) d = 0.0;
	} else {
		*hObs = 0;
		*hExp = 0;
	};
	return d;
}

// ------------------------------------------------------------------

void HetSumUp (float *hetSumAll, float *hetSumAve, float *hetWSumAve,
				float *hetLoc, float *wt, float *totWeight, float *NbHet,
				float *sumHrmonic, float *wsumHrmonic, int *polyLoc,
				long *mTotal, int nA, int count, float *hSamp,
				float *totWeigSq, float *wDsq, float *wSumhetSq)
{
// This is called after all alleles at a locus have been through,
// then parameters are either the measures for that locus or the
// cumulative ones.

// hetLoc: hetero. excess at current locus. On input, it is the sum of het.
// 	       excess of all alleles at the locus. On exit, it is the average.
// wDsq: (het-excess)^2 at current locus. On input, it is the sum of square
// 	      of het-excess of all alleles at the locus. On exit, it is the average.
// wt: weight of this locus
// These two are mainly for calculating the weight:
// nA: number of alleles considered in this locus. If there are j alleles
//     to be dropped, and the total alleles is m, then nA = m-j; so if
//     only one allele is left for this locus, the weight = 0.
// count: number of samples having data

// Cumulative:
// These two are for final het. excess used to calculate effective breeder
// totWeight: total of weights for loci up to this one.
// hetWSumAve: sum of (weighted ave. of het. excess per locus) up to this.

// add in Jan 2012:
// ---------------
// Formula (3) for the standard error SE in Pudovkin et al. paper
//       Sampling Properties ... Conserv. Genet. (2010)11,
// is given as:
//			SE = sqrt{[(S - Y^2/W)/(W*N)] * [W^2/(W^2 - V)]},
// where (i for alleles, j for loci, w_{ij} for weight, d_{ij} for het-ex)
//		S = Sum{i}Sum{j} w_{ij}d_{ij}^2			(weighted sum of (het-ex)^2)
//		Y = Sum{i}Sum{j} w_{ij}d_{ij}			(weighted sum of (het-ex))
//		W = Sum{i}Sum{j} w_{ij}					(total weights)
//		V = Sum{i}Sum{j} w_{ij}^2				(Sum of (weight)^2)
//		N = Sum{j}(n_j-1)						(n_j = # alleles at locus j)
//
// Let
//		N_j = number of individuals sampled at locus j
//		W_j = (n_j-1)*sqrt{N_j} = sum of weight w_{ij} of all alleles at locus j
//		W_{j,2} = W_j^2/n_j = sum of weight w_{ij}^2 of all alleles at locus j
//		V = Sum{j} W_{j,2}
//		D_j = unweighted average of d_{ij} over all alleles i at locus j
//		D_{j,2} = unweighted average of d_{ij}^2 over all alleles i at locus j
//		D   = weighted average of D_j over all loci j
//			= Y/W
//		D_2 = weighted average of D_{j,2} over all loci j
//			= S/W
// Then
//			SE = sqrt{[(D_2-D^2)/N}] * [W^2/(W^2 - V)]}
//
// Roles of variables in the code (left side is for the names of variables,
// ------------------------------  right side is for the notations above):
// j stands for the current locus,
//		nA     = n_j (input parameter)
//		hetLoc = D_j (at input, hetLoc = sum of all d_{ij}, then becomes D_j)
//		wDsq   = D_{j,2} (at input, wDsq = sum of all d_{ij}^2)
//		wt     = W_j (output parameter)
//		wt2    = W_{j,2} (local variable)
//		totWeight = Sum of W_p, p <= j,
//					then at j being the final locus, this will be W
//		totWeigSq = Sum of W_{p,2}, p <= j
//					then at j being the final locus, this will be V
//		hetWSumAve= Sum of W_pD_p, p <= j.
//					then at j being the final locus, this will be Y
//					(In here, hetWSumAve incremented by wt*hetLoc = W_j*D_j)
//		wSumhetSq = Sum of W_pD_{p,2}, p <= j,
//					then at j being the final locus, this will be S
//					(In here, wSumhetSq incremented by wt*wDsq = W_j*D_{j,2})
// -------------------------------------------------------------------------

// The following are for info purposes, not for the final Eff. Breeder output:
// hetSumAll: het. excess of all alleles up to this locus.
// hetSumAve: sum of (averages of het. excess per locus) up to this locus.
// NbHet: Effective breeder based on het. excess at this locus.
// sumHrmonic = sum of inverses of NbHet up to this locus
// wsumHrmonic = sum of weighted inverses of NbHet up to this locus
// polyLoc: number of loci having at least 2 alleles (weight > 0)
// mTotal: total number of alleles used in all loci up to this locus.


	float r, wt2, x;

	if (fabs(*hetLoc) < EPSILON) *hetLoc = 0.0;
	*hetSumAll += *hetLoc;
// bring those down: no need to divive hetLoc by nA if nA = 1
//	*hetLoc /= nA;	// hetLoc now is the average of het. excess at this locus
//	*hetSumAve += *hetLoc;	// cumulative sum of averages.
// On input, nA is the number of alleles used in this locus,
// On output, mTotal is the total alleles used up to this locus. Similar for
// totWeight, hetWSumAve (for calculation of weighted overall D-value).
// The others, hetSumAll, hetSumAve are similar, outputted in auxiliary file
// Actually, because we lump together all low freq alleles, nA always > 1
	if (nA > 1) {
		(*polyLoc)++;
		*hetLoc /= nA;	// hetLoc now is the average of het. excess at this locus
// add in Jan 2012:
		*wDsq /= nA;	// wDsq now is the average of (het. excess)^2 at this locus
// bring this out of the "if"
//		*mTotal += nA;
// add Mar 2012:
		*hSamp += 1.0F/((float) count);
	};
	*mTotal += nA;
	*hetSumAve += *hetLoc;	// cumulative sum of averages.
// add Jan 2012 for sum of het^2:
//	*wSumhetSq += *wDsq;	// cumulative sum of averages of (het-excess)^2.
	x = (float) nA;
	r = (float) sqrt((float) count) * (x-1);
// make sure we only consider loci having at least 2 alleles, never negative
// by rouding-off errror:
	if (r < 0) r = 0;
// The weight formula "wt" is extracted from the paper
// "Nb_HetEx: A Program ..., J. of Heredity, 2008:99(6)
// (The weight of each allele is r/x as stated in the paper)
	*wt = r;	// there are nA = x alleles, so the weight for locus is r.
	*totWeight += r;	// This will represent (Sum{i}Sum{j} w_{ij}) in
						// Pudovkin paper, formula (3) on p. 761
	*hetWSumAve += (*hetLoc)* r;
// add in Jan 2012 for the weighted sum of (het-excess)^2:
	*wSumhetSq += (*wDsq)* r;
// add in Jan 2012 for the sum of weight^2:
	wt2 = r*(r/x);	// this is the sum of weight^2 of all alleles at this locus.
	*totWeigSq += wt2;

	*NbHet = HetNb (*hetLoc);	// eff. breeder at this locus
	if (*NbHet < INFINITE) {
		*sumHrmonic += 1/(*NbHet);		// sum of inverses of NbHet
		*wsumHrmonic += (*wt/(*NbHet));	// sum of weighted inv. of NbHet
	};

}

// ------------------------------------------------------------------
// This belongs to LDmethod module:
// -------------------------------
float ExpR2Samp(float harmonic)
{

    if (harmonic == 0) return 0;
    if (harmonic >= 30)
		return (float) (1.0/harmonic + 3.19/(harmonic*harmonic));
    else
        return (float) (0.0018 + 0.907/harmonic + 4.44/(harmonic*harmonic));
}

// ---------------------------------------------------------------------------
// add in Nov 2012 for printing tabular form Frequencies
// modify in July 2013 for priting locus names
// Apr 2015: replace parameter FILE *locFile by struct locusMap *locList

void PrtLocFreq (ALLEPTR *alleList, int nloci, int nfish, int nMobil[],
				char *locUse, int *missptr, FILE *outLoc, char moreDat,
				int lenM, struct locusMap *locList)
// This routine go through all allele mobilitiy nodes in all "nloci" loci,
// calculate frequencies of these mobilities, put the values in the field
// "freq" of each node.
// Also output the minimum and maximum frequency values at each locus,
// calculate harmonic mean of populations across loci, based on samples
// having data, and use it to calculate r^2-sample value.
// If file outLoc not null, print outputs from this routine to an
// auxiliary output file.
{

	int nLocDat = nloci;	// nLocDat is the number of loci having data.
	int p, k, m, count;
//	int ndashes;

	int *alleNum, *alleDes;
	float q;
	float *freq;
	char locName[LEN_LOCUS];	// size for loc. names in locList
	char *alle;		// reserved string for prining alleles, will have size 5,
					// but Mac does not seem to like declaration alle[5]
	ALLEPTR curr, temp;
	int i, nAlle;
// ---------------------------------------------------------------------
	if (outLoc == NULL || moreDat == 0) return;
// This part is to count the number of different allele mobilities in
// population, put them in an array.
// Assuming that there are no more than 1000 different alleles
// (GENPOP format allows up to 3 digits for an allele), so we state
// the maximum as 1000. If allowing unlimited alleles, then should change
// to pointers (list) for counting.
	alle = (char*) malloc(sizeof(char)*5);
	for (m=0; m<5; m++) alle[m] = ' ';
	alleNum = (int*) malloc(sizeof(int)*(1000));
	for (k=0; k < 1000; *(alleNum + k) = 0, k++);
	nAlle = 0;	// count number of different alleles
	for (p=0; p<nloci; p++) {
		if (*(locUse+p) == 0) {
			nLocDat--;	// count the number of loci
			continue;
		}
// Go thru all mobility values at locus (p+1):
		curr = *(alleList+p);
		for ( ; curr != NULL; curr = temp) {
			temp = curr->next;
			m = curr->mValue;	// m always > 0
		// if m is an allele, then at index m, set alleNum = 1
		// (otherwise, alleNum still = 0), so when we go thru all
		// indices of alleNum, the ones with value 1 are allele mobilities
		// For convenience, we ignore index 0: *alleNum = 0 unchanged.
			if (*(alleNum+m) == 0) {
				*(alleNum+m) = 1;
				nAlle++;
			}
			count = nfish - (*(missptr+p));
			q = (float) (curr->copy)/(2*count);
			curr->freq = q;
		}
	}
	alleDes = (int*) malloc(sizeof(int)*(nAlle+1)); // we ignore index 0,
	for (i=0, m=1; m < 1000; m++) {			// so both i. m start from 1
		if (*(alleNum+m) != 0) {	// if non-zero, it is 1: there is
			i++;					// allele m, which is the ith allele
			*(alleNum+m) = i;	// so now if m is the value of an allele,
			// by looking at value at index m of alleNum, we know its order.
			*(alleDes+i) = m;
		}
	}
	freq = (float*) malloc(sizeof(float)*(nAlle+1)); // we ignore index 0,
	fprintf (outLoc, "Number of loci listed = %d\t(Column 3 is for the number of alleles)\n\n", nLocDat);

	if (locList != NULL) {
		fprintf (outLoc, "(Up to 10 righmost characters for locus names)\n");
		fprintf (outLoc, "Locus [#:Name]  \t Size \tAlleles:");
	} else fprintf (outLoc, "Locus           \t Size \tAlleles:");
	fflush(outLoc);
	for (i=1; i <= nAlle; i++) {
		// alle is string of 5 characters for writing allele up to its length
		sprintf(alle, "%5d", *(alleDes+i));
		// then fill in 0 in front of the number to make the length = lenM
		for (m=0; m<5; m++) if (alle[m] == ' ') alle[m] = '0';
		for (m=0; m<(5-lenM); m++) if (alle[m] == '0') alle[m] = ' ';
		fprintf (outLoc, "\t");
		for (m=0; m<5; m++) fprintf (outLoc, "%c", alle[m]);
		fprintf (outLoc, "   ");
// Mac does not like printing format below, so need to split up as above
//		fprintf (outLoc, "\t%s   ", alle);
//*/ The next code for sttraightforward printing alleles as numbers
//		fprintf (outLoc, "\t%5d   ", *(alleDes+i));
	}
	free(alle);
	fprintf (outLoc, "\n----------------\t------\t--------");
	for (i=1; i <= nAlle; i++) fprintf (outLoc, "\t--------");
	fprintf (outLoc, "\n");
	fflush(outLoc);
	k = 0;	// index for locList
	for (p=0; p<nloci; p++) {
		if (*(locUse+p) == 0) continue;
		if (locList != NULL) strcpy (locName, locList[k++].name);
		for (i=1; i<= nAlle; i++) *(freq+i) = 0;
		// will store frequencies at locus (p+1) into array freq.
		// At index i, it is frequency of this locus for the ith allele of
		// the population
		count = nfish - (*(missptr+p));
		curr = *(alleList+p);
		for ( ; curr != NULL; curr = temp) {
			temp = curr->next;	// allele m, as ordered in the population,
			m = curr->mValue;	// is the ith allele, i = alleNum at index m
			*(freq + (*(alleNum+m))) = (curr->freq);
		}
		i = *(nMobil+p);
		if (locList != NULL) fprintf (outLoc, "%5d:%-10s", (p+1), locName);
		else fprintf (outLoc, "%5d           ", (p+1));
		fprintf (outLoc, "\t%6d\t%7d ", count, i);
		for (i=1; i<= nAlle; i++) {
			if (*(freq+i) < EPSILON) fprintf (outLoc, "\t%8s", "0");
			else fprintf (outLoc, "\t%8.6f", *(freq+i));
		}
		fprintf (outLoc, "\n");
	}
	fprintf (outLoc, "----------------\t------\t--------");
	for (i=1; i <= nAlle; i++) fprintf (outLoc, "\t--------");
	fprintf (outLoc, "\n");

	fflush(outLoc);
	free (alleNum);
	free (alleDes);
	free (freq);

}

// --------------------------------------------------------------------------

// Apr 2015: replace parameter File *locFile by struct locusMap *locList

void Loc_Freq (ALLEPTR *alleList, int nloci, int nfish, float *NeWt,
				int nMobil[], int *missptr, char *locUse, float *minFreq,
				float *maxFreq, FILE *outLoc, char *outLocName,
				char moreDat, int popRead, char mHet, int lenM,
				struct locusMap *locList)
// This routine go through all allele mobilitiy nodes in all "nloci" loci,
// calculate frequencies of these mobilities, put the values in the field
// "freq" of each node.
// Also output the minimum and maximum frequency values at each locus,
// calculate harmonic mean of populations across loci, based on samples
// having data, and use it to calculate r^2-sample value.
// If file outLoc not null, print outputs from this routine to an
// auxiliary output file.
{
	int nLocDat = nloci;	// nLocDat is the number of loci having data.
	int p, k, count;
	float rmin, rmax;
// add few more in Aug 2012;
	float r, s;

	float q;
	float *freq;
	int *mobilVal;
// added Sept 2011:
	int mLoc = 0;
	char locPause = 0;

	ALLEPTR curr, temp;
//	float expR2;
	float harmonic = 0.0;

// To deal with heterozygote excess method --------------------------
// Add some other variables
	long mTotal = 0;	// mTotal = total alleles whose freq < 1.
	int locPoly = 0;	// #loci having more than 1 allele
//	int outNLoc = 100;	// maximum number of loci to output in outLoc file
	float *hObs, *hExp, *Nebp, NbAtp;
	float *dX, dXp, dSum, dXSum, dXWSum;
	float NbH, NbUw, NbHW;	// NbH = harmonic mean of Neb, NbHW for weighted one,
	// NbUw = eff. breeders based on average all Ds of all alleles across loci
	float hx;	// to hold the value of HetExAtLoc for an allele at a locus
	float *dW, dWp, dWtotal = 0.0; 	// weights for the D factor
// add Jan 2012
	float totWeigSq = 0;
	float wDsq, wSumhetSq, stdErr;
	float hSamp = 0;
	wSumhetSq = 0;
	// in the following, (het-ex/loc) = mean of het-ex for alleles at a locus
	NbH = 0;		// for harmonic mean of Nb across loci
	NbUw = 0;		// for Nb based on unweighted mean of (het-ex/loc)
	NbHW = 0;
	dSum = 0.0;		// for mean of all heterozygote excess
	dXSum = 0.0;	// for mean of (het-ex/loc)
	dXWSum = 0.0;	// for weighted mean of (het-ex/loc)
// ---------------------------------------------------------------------

	if (outLoc != NULL && moreDat == 1) {
		printf ("   Allele frequencies are being written to file %s.\n",
				outLocName);
		fprintf (outLoc, "\n\nPOPULATION %6d\t(Sample Size = %d)\n", popRead, nfish);
		for (k=0; k<17; k++) fprintf (outLoc, "*");
		fprintf (outLoc, "\n\n");

// To deal with heterozygote excess method --------------------------
// only allocate when auxiliary file outLoc is needed
		if (mHet == 1) {
			Nebp = (float*) malloc(sizeof(float)*(nloci));
			// for weights at loci:
			dW = (float*) malloc(sizeof(float)*(nloci));
			// for sum of unweighted het-excess:
			dX = (float*) malloc(sizeof(float)*(nloci));
			for (p=0; p<nloci; p++) {
				*(Nebp+p) = 0;
				*(dW+p) = 0;
				*(dX+p) = 0;
			};
		};
// ---------------------------------------------------------------------

	};
// add in Nov 2012:
	PrtLocFreq (alleList, nloci, nfish, nMobil, locUse, missptr,
				outLoc, moreDat, lenM, locList);

	for (p=0; p<nloci; p++) {
		if (*(locUse+p) == 0) {
			nLocDat--;
			continue;
		};

		// number of samples having data at locus (p+1):
		count = nfish - (*(missptr+p));
// test, and print to auxiliary file if needed:
// (added Sept 2011, locPause)
		if (outLoc != NULL && (moreDat == 1 && locPause == 0)) fprintf (outLoc,
			"\nLocus %d, individuals having data = %d\n", p+1, count);
		rmin = 1; rmax = 0;
		if (count == 0) {	// no data at this locus, then *(nMobil+p) = 0.
			nLocDat--;
		} else {
		// these two arrays are used to collect data to write to an auxiliary
		// file for allele mobilities, frequencies, harmonic mean, expected R2
			freq = (float*) malloc(sizeof(float)*(*(nMobil+p)));
			mobilVal = (int*) malloc(sizeof(int)*(*(nMobil+p)));

			harmonic += (float) (1.0/count);
			curr = *(alleList+p);

// To deal with heterozygote excess method --------------------------
// The d variable is the "D" at locus (p+1), as expressed in formula
// after formula (4) in Pudovkin & al. paper.
			if (mHet == 1) {
				dXp = 0;		// D factor at locus (p+1)
				NbAtp = 0;	// eff, breeeder at locus (p+1)
				hObs = (float*) malloc(sizeof(float)*(*(nMobil+p)));
				hExp = (float*) malloc(sizeof(float)*(*(nMobil+p)));
			// add Jan 2012:
				wDsq = 0;
			};
// ------------------------------------------------------------------

// Go thru all mobility values at locus (p+1):
			for (k=0; curr != NULL; curr = temp, k++) {
				temp = curr->next;
				*(mobilVal+k) = curr->mValue;
				q = (float) (curr->copy)/(2*count);
				curr->freq = q;
				if (q < rmin) rmin = q;
				if (q > rmax) rmax = q;
				*(freq+k) = q;

// To deal with heterozygote excess method --------------------------
// Frequency q of this allele (=mValue) is always < 1 except when
// this allele is unique at locus (p+1), or *(nMobil+p) = 1
				if (mHet == 1) {
				// dXp: total het. excess at this locus (p+1)
					hx = HetExAtLoc (*(nMobil+p), hExp+k, hObs+k, q,
					(curr->copy) - 2*(curr->homozyg), count);
					curr->hetx = hx;
					dXp += hx;
				// add Jan 2012:
					wDsq += (hx*hx);
				};
// ------------------------------------------------------------------

			};	// end of "(for k = 0; ..."

// To deal with heterozygote excess method --------------------------
// this block is within "else" of "if count == 0", so *(nMobil+p) > 0:
			if (mHet == 1) {
				HetSumUp (&dSum, &dXSum, &dXWSum, &dXp, &dWp, &dWtotal,
						&NbAtp, &NbH, &NbHW, &locPoly, &mTotal,
						*(nMobil+p), count, &hSamp, &totWeigSq,
						&wDsq, &wSumhetSq);
				if (outLoc != NULL && moreDat == 1)
				{
					*(dX+p) = dXp;
					*(dW+p) = dWp;
					*(Nebp+p) = NbAtp;
				};
			};

// ------------------------------------------------------------------

			// Output mobilities and their frequencies at locus (p+1)
			// as needed here (added Sept 2011, locPause):
			if (outLoc != NULL && (moreDat == 1) && locPause == 0) {
				fprintf (outLoc, "\tAlleles:    ");
				for (k=0; k<*(nMobil+p); k++)
//					fprintf (outLoc, "%10d", *(mobilVal+k));
					fprintf (outLoc, "%9d", *(mobilVal+k));
				fprintf (outLoc, "\n\tFrequencies:   ");
				for (k=0; k<*(nMobil+p); k++)
//					fprintf (outLoc, "%10.5f", *(freq+k));
					fprintf (outLoc, "%9.5f", *(freq+k));
				fprintf (outLoc, "\n");

// To deal with heterozygote excess method --------------------------
				if (mHet == 1 && (*(nMobil+p) > 1)) {
					fprintf (outLoc,  "\tExpected Het:  ");
					for (k=0; k<*(nMobil+p); k++)
//						fprintf (outLoc, "%10.5f", *(hExp+k));
						fprintf (outLoc, "%9.5f", *(hExp+k));
					fprintf (outLoc, "\n");
					fprintf (outLoc, "\tObserved Het:  ");
					for (k=0; k<*(nMobil+p); k++)
//						fprintf (outLoc, "%10.5f", *(hObs+k));
						fprintf (outLoc, "%9.5f", *(hObs+k));
					fprintf (outLoc, "\n");
//* For more details
					fprintf (outLoc, "\td=(Ob-Exp)/Exp:");
					dXp = 0;
// add Aug 2012
// add Average of d^2, and weight of locus (reuse variable q)
					wDsq = 0;
					r = 0;	// for overall Het-Obs
					s = 0;	// for overall Het-Exp
					for (k=0; k<*(nMobil+p); k++) {
						hx = ((*(hObs+k))-(*(hExp+k)))/(*(hExp+k));
						dXp += hx;
						wDsq += (hx*hx);
						fprintf (outLoc, "%9.5f", hx);
						r += (*(hObs+k));
						s += (*(hExp+k));
					};
					r /= 2;
					s /= 2;
					fprintf (outLoc, "\n");
					k = (*(nMobil+p));
					dXp /= ((float) k);
					wDsq /= ((float) k);
					q =(float) sqrt((float) count) * (k-1);
					fprintf (outLoc,
					"Mean d =%9.5f,     Mean d^2 =%8.5f, Weight =%8.2f\n", dXp, wDsq, q);
					q = (r-s)/s;
					fprintf (outLoc, "Overall Het (= Sum/2): Obs. (O) =%8.5f, "
								"Exp.(E)=%8.5f, (O-E)/E=%8.5f\n", r, s, q);
//*/
				};
// ------------------------------------------------------------------
				fflush (outLoc);
			}

			// release memories for those, to reclaim at next locus.
			free (mobilVal);
			free (freq);

// To deal with heterozygote excess method --------------------------
			if (mHet == 1) {
				free (hObs);
				free (hExp);
			};
// ------------------------------------------------------------------

		}	// end of "if (count == 0) { ... } else {"
		// just in case there are no data at all, rmin=1 has not been changed,
		// we don't want minimum value of frequencies to be 1.
//		if (rmin >= 1) rmin = 0;
		if (rmin >= 1) rmin = rmax;
		*(maxFreq+p) = rmax;
		*(minFreq+p) = rmin;
// test, and write to auxiliary file if needed (added locPause Sept 2011):
		if (outLoc != NULL && (moreDat == 1 && locPause == 0)) {
			fprintf(outLoc, "\nMin and Max Freq at locus %d:%10.5f,%10.5f\n",
					p+1, rmin, rmax);
			fflush (outLoc);
		};
// added Sept 2011:
		mLoc++;
		if (mLoc >= LOCOUTPUT) locPause = 1;

	}	// end of for (p=0; ...)
	if (harmonic > 0) harmonic = (float) nLocDat/harmonic;
//	expR2 = ExpR2Samp(harmonic);

// To deal with heterozygote excess method --------------------------
// Calculate effective breeders in various ways from:
//	* Average of all Ds in all loci,(dSum)
//	* Unweighted average of all averages Ds in loci (dXSum -> NbUw)
//	* Weighted average (weighted by locus) of all averages Ds in loci.
//	  (dXWSum -> *NeWt, a parameter of this function).
// and by taking
//	* Unweighted harmonic mean of all eff. breeders across loci (NbH)
//	* Weighted harmonic mean of all eff. breeders across loci (NbHW)
	if (mHet == 1) {
		HetAverage (locPoly, mTotal, &dSum, &dXSum, &dXWSum, &NbUw,
					&NbH, &NbHW, dWtotal, NeWt, &hSamp,
					0, totWeigSq, wSumhetSq, &stdErr);
	};
// ------------------------------------------------------------------

// test, and write to auxiliary file if needed:
	if (outLoc != NULL && (moreDat == 1)) {
		if (nLocDat > LOCOUTPUT)
			fprintf(outLoc, "\nOnly the first %d loci are listed", LOCOUTPUT);
		fprintf(outLoc, "\nTotal loci considered = %d\n", nLocDat);
		fprintf(outLoc,
			"Single-locus Harmonic Mean Sample Size  =%10.2f\n", harmonic);
// don't need this "Expected R^2 Sample" value in outLoc file anymore:
//		fprintf(outLoc, "Expected R^2 Sample based on this value =%14.6f\n", expR2);
		fprintf(outLoc, "\n");

// To deal with heterozygote excess method --------------------------
// print to auxiliary file data from heterozygote excess method
// Since the cut-off frequency 0+ is added (no restriction on frequencies)
// in the list of frequency cut-off values, the print out is generated
// from that case, so this print out becomes redundant and commented out:
//		if (mHet == 1) HetResults (outLoc, nloci, locUse, dX, dW,
//						Nebp, NbH, NbHW, *NeWt, NbUw, dSum);
// ------------------------------------------------------------------

		fflush (outLoc);
	};

// To deal with heterozygote excess method --------------------------
	if (mHet == 1 && outLoc != NULL && moreDat == 1) {
		free (Nebp);
		free (dW);
		free (dX);
	};
// ------------------------------------------------------------------
}


// --------------------------------------------------------------------------
// This function is used only for Het and LD methods

// Add Jan 2015: sepBurOut, moreCol
// Modified Dec 2016 for handling spec. cutoff
// Add 2 parameters: int nfish, int *missptr, 2 variables:  cutoff0, nSamp
// Reset cutoff at each locus in the "for" loop running across loci.

int Loci_Eligible (int nfish, int *missptr, float cutoff, ALLEPTR *alleList,
					int nloci, int *nMobil, float *minFreq, float *maxFreq,
					char *okLoc, int *lastOK,
					char *locUse, FILE *outLoc, FILE *outBurr,
					char moreDat, char moreBurr,
					char sepBurOut, char moreCol)
{
// This routine evaluates which locus should be accepted.
// Only choose loci that are not monomorphic, in the sense that their
// frequencies not so close to 1. A locus is accepted if all its
// frequencies of its alleles are at most 1 - cutoff, and has at least one
// allele having positive freq. >= cutoff, where cutoff is an input
// parameter, e.g. cutoff= 0.01.
// At this routine, we don't enforce yet accepted criteria for an allele.

// On input, locUse is an array of nloci boolean values representing loci
// that are to be considered (regardless of frequencies of its alleles).
// If locUse(p) = 0, locus (p+1) is to be skipped, 1 otherwise.
// On return, okLoc is an array of nloci values, okLoc(p) = 0 if locus
// (p+1) is not accepted. These values depend on cutoff value and locUse.
// okLoc(p) = 0 could mean that at locus (p+1):
//     * the locus is to be ignored by the user, delivered by locUse,
//     * all alleles have freq < cutoff
//     * one allele is dominant, that is, its freq > 1-cutoff, so the locus
//       is considered as monomorphic
//     * the locus is monomorphic (this falls into previous case if cutoff>0)
// okLoc(p) is the total number of alleles satisfying the requirement.
// Local variable nInd(p) represents the number of independent alleles.
// Suppose there are m alleles in the locus.
//     * if all have frequencies at least cutoff, then okLoc(p) = m,
//       but nInd(p) = m-1.
//     * if there is some allele having frequency < cutoff (case cutoff>0),
//       then okLoc(p) = nInd(p) = m.
// Thus, if okLoc(p)=1, it means that
//     * if no allele is dropped, the locus has two alleles,
//     * after dropping some allele, there is only one accepted allele left
//       that is not a dominant one, and the locus is not monomorphic
//       even there is only one allele to be used.
//       For a locus that has a dominant allele or is monomorphic, okLoc = 0.

// The number of independent of alleles in a locus, nInp(p), is outputted
// to an auxiliary file if needed.

// The value lastOK (when added 1 as usual) is the last accepted locus.
// It will be used when doing pairwise analysis, we want to stop the first
// locus before reaching the last one. This value directs the first locus
// not to reach it (so the second, which must differ from the first, can do).
// This is not crucial, however.

	int nLocOK;
	int p, q;
	ALLEPTR curr;
// Dec 2016: add 2 variables:
	float cutoff0;
	int nSamp;
	float plim1 = 1.0F;
	float plim2 = 1.0F -cutoff;
	int *nInd = (int*) malloc(sizeof(int)*nloci);

	*lastOK = -1;
	// assign okLoc = locUse, value of locUse is 1 by default, but may be
	// altered by the user. okLoc may take value 0 later when there is no
	// accepted allele in corresponding locus.
	for (p=0; p<nloci; p++) *(okLoc+p) = *(locUse+p);
	nLocOK = 0;
	for (p=0, q=0; p<nloci; p++) {
	// q becomes the total number of loci being considered

// Changed in Dec 2016, reset cutoff to drop only singleton alleles
// (All referenced cutoff in this loop (3 of them) are changed to cutoff0.)
		cutoff0 = cutoff;
		if (cutoff0 > 0 && cutoff0 <= PCRITX) {
			nSamp = nfish - *(missptr+p);
			// set cutoff0 to be slightly bigger than 1/(2*nSamp),
			// but smaller than 1/nSamp, so that if an allele has only one
			// copy, it will be dropped (having 2 copies will be accepted)
			// [1/(2*nSamp) < 1/(2*nSamp-1) <= 1/nSamp, equality iff nSamp=1]
			cutoff0 = 1.0;
			if (nSamp > 0) cutoff0 = 1.0/(2*nSamp-1);
		// either plim2 was set for all p at its declaration (before this
		// "for" loop), or reset at every p through this "if" statement
			plim2 = 1.0F - cutoff0;
		}

		if (*(okLoc+p) == 0) continue;
		q++;
		*(nInd+p) = *(nMobil+p);
//		if ((*(maxFreq+p) > 0) && (*(maxFreq+p) <= 1-cutoff)){
// having two conditions *(maxFreq+p) > 0 and >= cutoff) to prevent the case
// cutoff = 0, since we want only alleles having frequencies >= cutoff, so if
// cutoff = 0, we actually want all alleles. But if set condition >= cutoff
// only, then loci having no data (maxFreq then = 0) will be accepted!
// If cutoff = 0, then set condition (maxFreq < 1) to eliminate monomorphic
// locus.
//		if ((*(maxFreq+p) > 0) && (*(maxFreq+p) >= cutoff) &&
//								(*(maxFreq+p) <= 1-cutoff))
		if ((*(maxFreq+p) > 0) && (*(maxFreq+p) >= cutoff0) &&
						(*(maxFreq+p) <= plim2) && (*(maxFreq+p) < plim1))

		{	// at least one allele having freq >= cutoff and the locus is not
			// over represented by one allele. However, the locus may still
			// have only one accepted allele
			nLocOK++;
			*lastOK = p;
		// if all freq >= cutoff, N. indep alleles is decreased by 1.
		// Otherwise, it is decreased by the number of alleles having
		// freq. < cutoff.
		// The number of ind. alleles (of the accepted ones) is the number of
		// alleles that when their frequencies are known, then we know
		// frequencies  of all. Thus, when all freq. >= cutoff, all alleles
		// are accepted, the sum of their frequencies is 1, so we only
		// need to know frequencies of all except one.
		//
		// On the other hand, if we have to drop some allele, then the sum of
		// freq of the remaining alleles becomes unknown, so the number of
		// independent alleles is the number of accepted alleles left, since
		// the freq. of each one cannot be found based on the rest.
			if (*(minFreq+p) >= cutoff0) {
				(*(nInd+p))--;
				*(okLoc+p) = (*(nMobil+p)>0)? 1 : 0;
			} else {	// there is some allele to be dropped
				curr = *(alleList+p);
				for (; curr != NULL; curr = curr->next)
					if (curr->freq < cutoff0) (*(nInd+p))--;
				// since maxFreq >= cutoff, at least one allele accepted,
				// so after this loop, *(nInd+p) must still be > 0.
					*(okLoc+p) = (*(nInd+p)>0)? 1 : 0;
			}
		} else {
			*(nInd+p) = 0;
			*(okLoc+p) = 0;
		}
	}
	// write to auxiliary file for locus data as needed
	// parameters = 0 or 1 is for either Loc File or Burrow file:
	WriteLoci (outLoc, nloci, okLoc, cutoff, q, nMobil, nInd,
				nLocOK, moreDat, 0, LOCOUTPUT, sepBurOut, moreCol);
	WriteLoci (outBurr, nloci, okLoc, cutoff, q, nMobil, nInd,
				nLocOK, moreBurr, 1, LOCBURR, sepBurOut, moreCol);
	free (nInd);
	return nLocOK;

}
// --------------------------------------------------------------------------




// The rest are for dealing with Jackknife process


//------------------------------------------------------------------------
long JackKnifeInd(float mean, float variance)
{
	long iBig;
	float phi;
	if (mean == 0) return 1;
    phi = variance/(mean*mean);
    if (phi <= EPSILON) iBig = MAXDEG;
	else iBig = (long) floor(2.0/phi + 0.5);	// round off 2/phi
	if (iBig == 0) iBig = 1;
//printf("Get Degree of Freedom = %10lu\n", iBig);
    return iBig;
}

//------------------------------------------------------------------------


float GetChi (float z, long degfree)
{
    float a, sqrta;
	float d = (float) degfree;
    a = 2.0F/(9.0F*d);
    sqrta = (float) sqrt(a);
    return ((float)pow((1.0F - a + z*sqrta), 3));

}

//-------------------------------------------------------------------------

void Confid95(long degfree, float fmean, float *lowlim, float *uplim)
{
	float xhi, xlo;
	long n = degfree-1;
	float high[100], low[100];
	high[0]=0.001F;    high[1]=0.05F;    high[2]=0.22F;    high[3]=0.48F;
	high[4]=0.83F;     high[5]=1.24F;    high[6]=1.69F;    high[7]=2.18F;
	high[8]=2.70F;     high[9]=3.25F;    high[10]=3.82F;   high[11]=4.40F;
	high[12]=5.01F;    high[13]=5.63F;   high[14]=6.27F;   high[15]=6.91F;
	high[16]=7.56F;    high[17]=8.23F;   high[18]=8.91F;   high[19]=9.59F;
	high[20]=10.28F;   high[21]=10.98F;  high[22]=11.69F;  high[23]=12.40F;
	high[24]=13.12F;   high[25]=13.84F;  high[26]=14.57F;  high[27]=15.31F;
	high[28]=16.05F;   high[29]=16.79F;  high[30]=17.55F;  high[31]=18.32F;
	high[32]=19.08F;   high[33]=19.85F;  high[34]=20.61F;  high[35]=21.37F;
	high[36]=22.14F;   high[37]=22.90F;  high[38]=23.67F;  high[39]=24.43F;
	high[40]=25.22F;   high[41]=26.02F;  high[42]=26.81F;  high[43]=27.60F;
	high[44]=28.40F;   high[45]=29.19F;  high[46]=29.98F;  high[47]=30.77F;
	high[48]=31.57F;   high[49]=32.36F;  high[50]=33.17F;  high[51]=33.98F;
	high[52]=34.80F;   high[53]=35.61F;  high[54]=36.42F;  high[55]=37.23F;
	high[56]=38.04F;   high[57]=38.86F;  high[58]=39.67F;  high[59]=40.48F;
	high[60]=41.31F;   high[61]=42.14F;  high[62]=42.96F;  high[63]=43.79F;
	high[64]=44.62F;   high[65]=45.45F;  high[66]=46.28F;  high[67]=47.10F;
	high[68]=47.93F;   high[69]=48.76F;  high[70]=49.60F;  high[71]=50.44F;
	high[72]=51.2648F; high[73]=52.12F;  high[74]=52.96F;  high[75]=53.79F;
	high[76]=54.63F;   high[77]=55.47F;  high[78]=56.31F;  high[79]=57.15F;
	high[80]=58.00F;   high[81]=58.85F;  high[82]=59.70F;  high[83]=60.55F;
	high[84]=61.40F;   high[85]=62.25F;  high[86]=63.10F;  high[87]=63.95F;
	high[88]=64.80F;   high[89]=65.65F;  high[90]=66.51F;  high[91]=67.36F;
	high[92]=68.22F;   high[93]=69.08F;  high[94]=69.94F;  high[95]=70.79F;
	high[96]=71.65F;   high[97]=72.51F;  high[98]=73.36F;  high[99]=74.22F;

	low[0]=5.02F;      low[1]=7.38F;     low[2]=9.35F;     low[3]=11.14F;
	low[4]=12.83F;     low[5]=14.45F;    low[6]=16.01F;    low[7]=17.53F;
	low[8]=19.02F;     low[9]=20.48F;    low[10]=21.92F;   low[11]=23.34F;
	low[12]=24.74F;    low[13]=26.12F;   low[14]=27.49F;   low[15]=28.85F;
	low[16]=30.19F;    low[17]=31.53F;   low[18]=32.85F;   low[19]=34.17F;
	low[20]=35.48F;    low[21]=36.78F;   low[22]=38.08F;   low[23]=39.36F;
	low[24]=40.65F;    low[25]=41.92F;   low[26]=43.19F;   low[27]=44.46F;
	low[28]=45.72F;    low[29]=46.98F;   low[30]=48.22F;   low[31]=49.45F;
	low[32]=50.69F;    low[33]=51.92F;   low[34]=53.16F;   low[35]=54.40F;
	low[36]=55.63F;    low[37]=56.87F;   low[38]=58.10F;   low[39]=59.34F;
	low[40]=60.55F;    low[41]=61.76F;   low[42]=62.96F;   low[43]=64.17F;
	low[44]=65.38F;    low[45]=66.59F;   low[46]=67.80F;   low[47]=69.00F;
	low[48]=70.21F;    low[49]=71.42F;   low[50]=72.61F;   low[51]=73.80F;
	low[52]=74.98F;    low[53]=76.17F;   low[54]=77.36F;   low[55]=78.55F;
	low[56]=79.74F;    low[57]=80.92F;   low[58]=82.11F;   low[59]=83.30F;
	low[60]=84.47F;    low[61]=85.64F;   low[62]=86.82F;   low[63]=87.99F;
	low[64]=89.16F;    low[65]=90.33F;   low[66]=91.50F;   low[67]=92.68F;
	low[68]=93.85F;    low[69]=95.02F;   low[70]=96.18F;   low[71]=97.34F;
	low[72]=98.5162F;  low[73]=99.66F;   low[74]=100.83F;  low[75]=101.99F;
	low[76]=103.15F;   low[77]=104.31F;  low[78]=105.47F;  low[79]=106.63F;
	low[80]=107.78F;   low[81]=108.93F;  low[82]=110.08F;  low[83]=111.23F;
	low[84]=112.39F;   low[85]=113.54F;  low[86]=114.69F;  low[87]=115.84F;
	low[88]=116.99F;   low[89]=118.14F;  low[90]=119.28F;  low[91]=120.42F;
	low[92]=121.57F;   low[93]=122.71F;  low[94]=123.85F;  low[95]=124.99F;
	low[96]=126.13F;   low[97]=127.28F;  low[98]=128.42F;  low[99]=129.56F;

    if (degfree <= 100) {
        xhi = high[n];
        xlo = low[n];
    	*uplim = (degfree*fmean)/xhi;
    	*lowlim = (degfree*fmean)/xlo;
	} else {
        xhi = GetChi(-1.96F, degfree);
        xlo = GetChi(1.96F, degfree);
    	*uplim = fmean/xhi;
    	*lowlim = fmean/xlo;
    };

}

//-------------------------------------------------------------------------

void t_Confid9x(long degfree, float fmean, float stdErr,
				float *lolim, float *hilim, char wide)
// wide > 0 for 95% CI, otherwise, 90% CI
{
	float t;
	long n = degfree-1;
	float hi95[100], hi90[100];	// [n] corresponds to deg. of freedom (n+1)
	float bound;

	hi95[0]=12.706F;  hi95[1]=4.303F;   hi95[2]=3.182F;   hi95[3]=2.776F;
	hi95[4]=2.571F;   hi95[5]=2.447F;   hi95[6]=2.365F;   hi95[7]=2.306F;
	hi95[8]=2.262F;   hi95[9]=2.228F;   hi95[10]=2.201F;  hi95[11]=2.179F;
	hi95[12]=2.160F;  hi95[13]=2.145F;  hi95[14]=2.131F;  hi95[15]=2.120F;
	hi95[16]=2.110F;  hi95[17]=2.101F;  hi95[18]=2.093F;  hi95[19]=2.086F;
	hi95[20]=2.080F;  hi95[21]=2.074F;  hi95[22]=2.069F;  hi95[23]=2.064F;
	hi95[24]=2.060F;  hi95[25]=2.056F;  hi95[26]=2.052F;  hi95[27]=2.048F;
	hi95[28]=2.045F;  hi95[29]=2.042F;  hi95[30]=2.040F;  hi95[31]=2.037F;
	hi95[32]=2.035F;  hi95[33]=2.032F;  hi95[34]=2.030F;  hi95[35]=2.028F;
	hi95[36]=2.026F;  hi95[37]=2.024F;  hi95[38]=2.023F;  hi95[39]=2.021F;
	hi95[40]=2.020F;  hi95[41]=2.018F;  hi95[42]=2.017F;  hi95[43]=2.015F;
	hi95[44]=2.014F;  hi95[45]=2.013F;  hi95[46]=2.012F;  hi95[47]=2.011F;
	hi95[48]=2.010F;  hi95[49]=2.009F;  hi95[50]=2.008F;  hi95[51]=2.007F;
	hi95[52]=2.006F;  hi95[53]=2.005F;  hi95[54]=2.004F;  hi95[55]=2.003F;
	hi95[56]=2.002F;  hi95[57]=2.002F;  hi95[58]=2.001F;  hi95[59]=2.000F;
	hi95[60]=2.000F;  hi95[61]=1.999F;  hi95[62]=1.998F;  hi95[63]=1.998F;
	hi95[64]=1.997F;  hi95[65]=1.997F;  hi95[66]=1.996F;  hi95[67]=1.995F;
	hi95[68]=1.995F;  hi95[69]=1.994F;  hi95[70]=1.994F;  hi95[71]=1.993F;
	hi95[72]=1.993F;  hi95[73]=1.993F;  hi95[74]=1.992F;  hi95[75]=1.992F;
	hi95[76]=1.991F;  hi95[77]=1.991F;  hi95[78]=1.990F;  hi95[79]=1.990F;
	hi95[80]=1.990F;  hi95[81]=1.989F;  hi95[82]=1.989F;  hi95[83]=1.989F;
	hi95[84]=1.988F;  hi95[85]=1.988F;  hi95[86]=1.988F;  hi95[87]=1.987F;
	hi95[88]=1.987F;  hi95[89]=1.987F;  hi95[90]=1.986F;  hi95[91]=1.986F;
	hi95[92]=1.986F;  hi95[93]=1.986F;  hi95[94]=1.985F;  hi95[95]=1.985F;
	hi95[96]=1.985F;  hi95[97]=1.984F;  hi95[98]=1.984F;  hi95[99]=1.984F;

	hi90[0]=6.314F;   hi90[1]=2.920F;   hi90[2]=2.353F;   hi90[3]=2.132F;
	hi90[4]=2.015F;   hi90[5]=1.943F;   hi90[6]=1.895F;   hi90[7]=1.860F;
	hi90[8]=1.833F;   hi90[9]=1.812F;   hi90[10]=1.796F;  hi90[11]=1.782F;
	hi90[12]=1.771F;  hi90[13]=1.761F;  hi90[14]=1.753F;  hi90[15]=1.746F;
	hi90[16]=1.740F;  hi90[17]=1.734F;  hi90[18]=1.729F;  hi90[19]=1.725F;
	hi90[20]=1.721F;  hi90[21]=1.717F;  hi90[22]=1.714F;  hi90[23]=1.711F;
	hi90[24]=1.708F;  hi90[25]=1.706F;  hi90[26]=1.703F;  hi90[27]=1.701F;
	hi90[28]=1.699F;  hi90[29]=1.697F;  hi90[30]=1.696F;  hi90[31]=1.694F;
	hi90[32]=1.692F;  hi90[33]=1.691F;  hi90[34]=1.690F;  hi90[35]=1.688F;
	hi90[36]=1.687F;  hi90[37]=1.686F;  hi90[38]=1.685F;  hi90[39]=1.684F;
	hi90[40]=1.683F;  hi90[41]=1.682F;  hi90[42]=1.681F;  hi90[43]=1.680F;
	hi90[44]=1.679F;  hi90[45]=1.679F;  hi90[46]=1.678F;  hi90[47]=1.677F;
	hi90[48]=1.677F;  hi90[49]=1.676F;  hi90[50]=1.675F;  hi90[51]=1.675F;
	hi90[52]=1.674F;  hi90[53]=1.674F;  hi90[54]=1.673F;  hi90[55]=1.673F;
	hi90[56]=1.672F;  hi90[57]=1.672F;  hi90[58]=1.671F;  hi90[59]=1.671F;
	hi90[60]=1.670F;  hi90[61]=1.670F;  hi90[62]=1.669F;  hi90[63]=1.669F;
	hi90[64]=1.669F;  hi90[65]=1.668F;  hi90[66]=1.668F;  hi90[67]=1.668F;
	hi90[68]=1.667F;  hi90[69]=1.667F;  hi90[70]=1.667F;  hi90[71]=1.666F;
	hi90[72]=1.666F;  hi90[73]=1.666F;  hi90[74]=1.665F;  hi90[75]=1.665F;
	hi90[76]=1.665F;  hi90[77]=1.665F;  hi90[78]=1.664F;  hi90[79]=1.664F;
	hi90[80]=1.664F;  hi90[81]=1.664F;  hi90[82]=1.663F;  hi90[83]=1.663F;
	hi90[84]=1.663F;  hi90[85]=1.663F;  hi90[86]=1.663F;  hi90[87]=1.662F;
	hi90[88]=1.662F;  hi90[89]=1.662F;  hi90[90]=1.662F;  hi90[91]=1.662F;
	hi90[92]=1.661F;  hi90[93]=1.661F;  hi90[94]=1.661F;  hi90[95]=1.661F;
	hi90[96]=1.661F;  hi90[97]=1.661F;  hi90[98]=1.660F;  hi90[99]=1.660F;

	if (wide > 0) { // do 95% CI
		bound = 1.96F;
		if (degfree <= 100) t = hi95[n];
		else t = bound;
	} else {		// do 90% CI
		bound = 1.645F;
		if (degfree <= 100) t = hi90[n];
		else t = bound;
	};
	*lolim = fmean - t*stdErr;
	*hilim = fmean + t*stdErr;

/*
t-values t(0.025):
12.706  4.303  3.182  2.776  2.571  2.447  2.365  2.306  2.262  2.228
 2.201  2.179  2.160  2.145  2.131  2.120  2.110  2.101  2.093  2.086
 2.080  2.074  2.069  2.064  2.060  2.056  2.052  2.048  2.045  2.042
 2.040  2.037  2.035  2.032  2.030  2.028  2.026  2.024  2.023  2.021
 2.020  2.018  2.017  2.015  2.014  2.013  2.012  2.011  2.010  2.009
 2.008  2.007  2.006  2.005  2.004  2.003  2.002  2.002  2.001  2.000
 2.000  1.999  1.998  1.998  1.997  1.997  1.996  1.995  1.995  1.994
 1.994  1.993  1.993  1.993  1.992  1.992  1.991  1.991  1.990  1.990
 1.990  1.989  1.989  1.989  1.988  1.988  1.988  1.987  1.987  1.987
 1.986  1.986  1.986  1.986  1.985  1.985  1.985  1.984  1.984  1.984


t-values t(0.05):
 6.314  2.920  2.353  2.132  2.015  1.943  1.895  1.860  1.833  1.812
 1.796  1.782  1.771  1.761  1.753  1.746  1.740  1.734  1.729  1.725
 1.721  1.717  1.714  1.711  1.708  1.706  1.703  1.701  1.699  1.697
 1.696  1.694  1.692  1.691  1.690  1.688  1.687  1.686  1.685  1.684
 1.683  1.682  1.681  1.680  1.679  1.679  1.678  1.677  1.677  1.676
 1.675  1.675  1.674  1.674  1.673  1.673  1.672  1.672  1.671  1.671
 1.670  1.670  1.669  1.669  1.669  1.668  1.668  1.668  1.667  1.667
 1.667  1.666  1.666  1.666  1.665  1.665  1.665  1.665  1.664  1.664
 1.664  1.664  1.663  1.663  1.663  1.663  1.663  1.662  1.662  1.662
 1.662  1.662  1.661  1.661  1.661  1.661  1.661  1.661  1.660  1.660
*/

}


//-------------------------------------------------------------------------
// From HetExcess module.
//
// ------------------------------------------------------------------

// Following Pudovkin et al. paper 2010
void CI_t_DistHet (long degfree, float dmean, float stdErr,
				  float *lowNe, float *highNe, char mode)
{
	float loLim, hiLim;
	t_Confid9x (degfree, dmean, stdErr, &loLim, &hiLim, mode);
	*lowNe = HetNb (hiLim);
	*highNe = HetNb (loLim);
}


// ------------------------------------------------------------------
float HetxLow (FISHPTR *fishList, int nLowF, float lastHx, int *smAlle,
			   float totF, int nfish)
{
// nlowF: number of alleles with low freq,
// lastHx: the last hx read in the calling function before this is called
// smAllle: list of low freq alleles
// totF: total frequencies of those alleles
// nfish: number of samples having data
// Return value: hx, heterozygote excess, when all those alleles are
// lumped together, considered as only one allele.
	int het, homo, i, j;
	int a0, a1;
	float hExp, hObs, hx;
	FISHPTR curr = *fishList;
	if (nLowF == 0) return 0;
	if (nLowF == 1) return lastHx;
	het = 0;
	homo  = 0;
	for ( ; curr != NULL; curr = curr->next) {
		a0 = curr->gene[0];
		a1 = curr->gene[1];
		for (i=0; i<nLowF; i++) {
			if ((a0 != *(smAlle+i)) && (a1 != *(smAlle+i)))
				continue;
			else break;	// smAlle+i is one of the alleles
		};
		if (i >= nLowF) continue;	// no match
// After the for loop and there is a match at i<nLowF, we know one
// of the two alleles is smAlle[i].
	// smAlle+i is one of the alleles, check the other one
		if (a0 == *(smAlle+i)) {
			for (j=0; j<nLowF; j++) {
				if (a1 != *(smAlle+j))
					continue;
				else break;	// smAlle+j is the other
			};
			if (j >= nLowF) het++;
			else homo++;
		} else { // (a1 == *(smAlle+i))
// At 1st loop of (i=0; ..), when (a0 != *(smAlle+i)) && (a1 != *(smAlle+i)))
// is evaluated, the first condition (a0 != *(smAlle+i)) is evaluated first,
// then (a1 != *(smAlle+i))) is evaluated. Only when the first is OK, the
// second will be evaluated.
// Thus, get out of the loop by "break", we have
// either
// * The first fails (then a0 = smAlle[i]) and the second is unkown,
//   or
// * The first holds (a0 != smAlle[i]), but the second fails (a1 = smAlle[i])
//   However, we can not say right away that this is heterozygote.
//   The reason: index i runs from 0 to (nLowF-1), the second allele a1 may
//   match smAlle[i] then break, but for some j>i, we may have a0=smAlle[j].
//   Thus, the next for loop may start from j=i+1, but as precaution, let j=0
			for (j=0; j<nLowF; j++) {
				if (a0 != *(smAlle+j))
					continue;
				else break;	// smAlle+j is the other
			};
			if (j >= nLowF) het++;
			else homo++;
		};
	};
	hObs = (float) het /((float) nfish);
	hExp = HetExp (totF, nfish);
	hx = (hObs - hExp)/(hExp);
	if (fabs(hx) < EPSILON) hx = 0.0;
	return hx;
}

// ------------------------------------------------------------------



void HetXcess (FISHPTR *fishList, ALLEPTR *alleList, int nloci, int nfish,
				int nMobil[], int *missptr, char *okLoc, FILE *outLoc,
				char moreDat, float cutoff, float *hetWSumAve, float *NeWt,
				long *nIndH, float *hSamp, float *loNe, float *hiNe,
				char param)
{
	float hx, NbHet;	// for each locus
	float *wt, *hetLoc;
	float freq;
	int p, count, k, nA;
// for lumping alleles
	char quit;
	int *smAlle;	// to hold alleles with low freq.
	float totF, hxsm;
	ALLEPTR curr;

	float NbUw;			// Nb based on hetSumAve, (NeWt based on hetWSumAve)
	long mTotal = 0;	// total alleles
	int polyLoc = 0;	// total number of polymorphic loci
	float hetSumAll = 0;	// sum of het. ex. of all alleles
	float hetSumAve = 0;	// unweighted ave. of ave. het ex. across loci
	float totWeight = 0;
	float sumHrmonic = 0;	// unweighted harmonic average of Neb per locus,
	float wsumHrmonic = 0;	// weighted harmonic average of Neb per locus,
	int indAlle = 0;
	int nlowF = 0;
	float x, t;
// add Jan 2012:
	float wDsq, wSumhetSq, stdErr;
	float totWeigSq = 0;
	char outCI = 0;	// set 1 for outputting CI in outLoc file for each locus,
					// 0 otherwise

// temporary add:
/*
float totwHetsq = 0;
float totwHet = 0;
float q, s, u, r1, r2;
//*/
	wSumhetSq = 0;
//	float hetWSumAve = 0;	// weighted ave. of ave. het ex. across loci
	*hetWSumAve = 0;		// weighted ave. of ave. het ex. across loci
	*nIndH = 0;				// total indep. alleles
	*hSamp = 0;	// harmonic mean of sample size among all poly. loci
	hetLoc = (float*) malloc (sizeof(float)*nloci);
	wt = (float*) malloc (sizeof(float)*nloci);

	for (p=0; p<nloci; p++) {
		if (*(okLoc+p) == 0) continue;	// marked by Loci_Eligible to be skipped.
		k = *(nMobil+p);
		if (k < 2) continue;
		smAlle = (int*) malloc (sizeof(int)*k);

		// number of samples having data at locus (p+1):
		count = nfish - (*(missptr+p));
		nA = 0;
		*(hetLoc+p) = 0;
		*(wt+p) = 0;
		indAlle = 0;
		nlowF = 0;
		quit = 1;
		totF = 0;
// add Jan 2012:
		wDsq = 0;
		for (curr = *(alleList+p); curr != NULL; curr = curr->next)
		{
			if (*(nMobil+p) < 2) continue;
			freq = curr->freq;
			if (freq == 0 || freq > 1-cutoff) continue;
			if (freq < cutoff) {
				*(smAlle+nlowF) = curr->mValue;
				nlowF++;		// increment # of low freq. alleles
				totF += freq;	// accumulate their freq.
				hxsm = curr->hetx;
				continue;
			};
			nA++;	// increment the # of alleles used in calculating D
					// this allele has freq >= cutoff and <= 1-cutoff
			hx = curr->hetx;
			*(hetLoc+p) += hx;
		// add Jan 2012:
			wDsq += (hx*hx);
			quit = 0;	// this locus has at least 1 accepted allele
		};	// went through all alleles at locus (p+1)
		if (quit == 1) {
			free (smAlle);
			continue;
		};
		if (nlowF > 0) nA++;	// increase total alleles because all
								// low freq alleles are lumped together
	// change in Jan 2012, set hx to be the value of HetxLow, *(hetLoc+p)
	// is incremented by hx, same result as before, and now, increment
	// wDsq by (hx)^2:
		if (nlowF > 0) {	// fishList is not allocated if no critical value
		// other than 0+. In such case, nlowF = 0, and (fishList+p) is not
		// allocated, so this condition is to avoid referring to (fishList+p)
			hx = HetxLow ((fishList+p), nlowF, hxsm, smAlle, totF, count);
			*(hetLoc+p) += hx;
			wDsq += (hx*hx);
		};
		free (smAlle);
// Put condition nA > 0 just in case this locus has no allele freq accepted,
// but it may be redundant because quit should be 0:
		if (nA > 0) indAlle = nA - 1;
		else indAlle = 0;
		*nIndH += indAlle;

// temporary add
/*
q = (float) nA;
s = (float) count;
u = ((q - 1.0F)/q)*((float)sqrt(s));
s = u*wDsq;
t = u*(*(hetLoc+p));
totwHetsq += s; // will be sum{i}sum{j} of w_{ij} d_{ij}^2 in the paper
totwHet += t;  // will be sum{i}sum{j} of w_{ij} d_{ij} in the paper
//*/
		HetSumUp (&hetSumAll, &hetSumAve, hetWSumAve, (hetLoc+p), (wt+p),
					&totWeight, &NbHet, &sumHrmonic, &wsumHrmonic,
					&polyLoc, &mTotal, nA, count, hSamp, &totWeigSq,
					&wDsq, &wSumhetSq);
// now, *(hetLoc+p) is mean of hets, *wDsq is mean of het^2, at locus (p+1),
// *(wt+p) is the weight of locus (p+1).
		if (outLoc != NULL && (moreDat == 1)) {
			x = *(wt+p);
// add more outputs in Aug 2012:
			t = (x*x)/((float)nA);	// equals to sum of (allele-weight)^2
			// this was calculated in HetSumUp, recall here since we want to
			// add to output for outLoc (can have as parameter in HetSumUp)
			if (polyLoc == 1) {
				PrtLines (outLoc, 79, '=');
				HetForeword (outLoc);
				if (outCI == 1)
					fprintf (outLoc, "\nLocus  #smp. #alle.   Dm     Dm2      Wt     "
						"  Wt2     N[eb]        95%% CI\n");
// headline includes details but no CI:
				else fprintf (outLoc, "\nLocus  #smp. #alle.   Dm     Dm2      Wt      "
					"Dm*Wt   Dm2*Wt      Wt2     N[eb]\n");
			};
			if (outCI == 1) {
			// calculate Standard Error for the average D at this locus
				stdErr = (wDsq - (*(hetLoc+p))*(*(hetLoc+p)))/((float) indAlle);
				if (stdErr < 0) stdErr = 0;
				else stdErr = (float) sqrt(stdErr);
				CI_t_DistHet (indAlle, *(hetLoc+p), stdErr, loNe, hiNe, 1);
				fprintf (outLoc, "%5d%6d%5d%9.4f%8.4f%8.2f%10.1f",
					(p+1), count, nA, *(hetLoc+p), wDsq, x, t);
				if (NbHet < INFINITE) fprintf (outLoc, "%9.1f", NbHet);
				else fprintf (outLoc, "%9s", "INF");
				if ((*loNe > 0) && (*loNe < INFINITE))
					fprintf (outLoc, "%9.1f", *loNe);
				else fprintf (outLoc, "%9s", "INF");
				if ((*hiNe > 0) && (*hiNe < INFINITE))
					fprintf (outLoc, "%9.1f", *hiNe);
				else fprintf (outLoc, "%9s", "INF");
//fprintf (outLoc, "%10.4f", stdErr);
				fprintf (outLoc, "\n");
			} else {
// output includes details but no CI:
				fprintf (outLoc,
					"%5d%6d%5d%9.4f%8.4f%8.2f%10.4f%9.4f%10.1f",
					(p+1), count, nA, *(hetLoc+p), wDsq, x, (*(hetLoc+p))*x,
					wDsq*x, t);
				if (NbHet < INFINITE) fprintf (outLoc, "%9.1f\n", NbHet);
				else fprintf (outLoc, "%9s\n", "INF");
			};

// previous output, no wt2 item:
//			if (outLoc != NULL && (moreDat == 1)) fprintf (outLoc,
//			"%5d%8d%12.4f%11.2f%12.1f\n", (p+1), nA, *(hetLoc+p), *(wt+p), NbHet);
		};
	};
// after calling HetAverage, hetWSumAve becomes the overall D-value,
// NeWt will be the output of Nb in the main output file
	t = *hetWSumAve;	// hold value for *hetWSumAve since it will be changed
	HetAverage (polyLoc, mTotal, &hetSumAll, &hetSumAve,
				hetWSumAve, &NbUw, &sumHrmonic, &wsumHrmonic, totWeight, NeWt,
				hSamp, *nIndH, totWeigSq, wSumhetSq, &stdErr);
// Corresponding parameter names listed in HetAverage
// locPoly = polyLoc, mTotal (same), *dSum = &hetSumAll, *dXSum = &hetSumAve,
// *dXWSum = hetWSumAve (hetWSumAve is a pointer, parameter), *NbUw = &NbUw,
// *NbH = &sumHrmonic, *NbHW = &wsumHrmonic, dWtotal = totWeight,
// *NeWt = NeWt, *hSamp = hSamp, (NeWt, hSamp are pointers, parameters)
// nInd =  *nIndH (nInd is a pointer, parameter),
// the rest have same names.

// for 95% confidence interval:
//	mTotal -= polyLoc;	// represent degree of freedom

// temporary add
/*
q = (float) (*nIndH);
r1 = (totwHetsq - (totwHet*totwHet)/totWeight)/(totWeight * q);
s = totWeight*totWeight;
r2 = s/(s - totWeigSq);
x = (float) sqrt(r1*r2);
//*/
	CI_t_DistHet (*nIndH, *hetWSumAve, stdErr, loNe, hiNe, 1);

// print to console (moved to the calling function RunPop#)
	printf ("     Heterozygote Excess Method\n");
	printf ("       Estimated Neb: ");
	if (*NeWt < INFINITE) printf ("%21.1f\n", *NeWt);
	else printf ("%21s\n", "Infinite");
	if (outLoc != NULL && (moreDat == 1)) {
		PrtLines (outLoc, 79, '-');
		if (outCI == 1)
			fprintf (outLoc, "SUM:%12lu%25.2f%10.1f\n", mTotal,
					totWeight, totWeigSq);
		else fprintf (outLoc, "SUM:%12lu%25.2f%10.4f%9.4f%10.1f\n", mTotal,
					totWeight, t, wSumhetSq, totWeigSq);
//		PrtDetailHet (outLoc, *hetWSumAve, *nIndH, totWeight, wSumhetSq, totWeigSq);
		PrtHetSum (outLoc, sumHrmonic, wsumHrmonic, *NeWt, NbUw, hetSumAll,
					*hetWSumAve, stdErr);
		PrtLines (outLoc, 79, '=');

//fprintf (outLoc, "\nStandard Error directly from formula (3) in the paper =%9.6f", x);
	};
	free (hetLoc);
	free (wt);
}

// --------------------------------------------------------------------------

// ------------------------  End of CalFreq & HetExcess ---------------------

// LDmethod
// --------------------------------------------------------------------------

int Count(int gene[], int m)
{
	int i, j;
	for (i=0, j=0; i<2; i++)
		if (gene[i] == m) j++;
	return j;
}

// --------------------------------------------------------------------------


float LD_Ne (float harmonic, float rPrime, char matingMod, float infinite)
{

    float estNe;

    float x;
    if (rPrime == 0) return (float) infinite;

	if (harmonic >= 30) {
        if (matingMod == 0) {
        // random mating: x = discriminant in quadratic rN^2 - (1/3)N + 0.69 = 0,
        // where r=rPrime, N = estNe. The eq. came from r = 1/(3N) - 0.69/N^2.
            x = (float) (1.0/9.0 - 2.76*rPrime);
            x = (x > 0) ? x : 0;
            estNe = (float) (1.0/3.0 + sqrt(x));
        } else {
            x = (float) (4.0/9.0 - 7.2*rPrime);
            x = (x > 0) ? x : 0;
            estNe = (float) (2.0/3.0 + sqrt(x));
        }
    } else {
        if (matingMod == 0) {
        // random mating: X = discriminant in quadratic rN^2 - .308 N + 0.52 = 0,
        // where r=rPrime, N = estNe. The eq. came from r = .308/N - 0.52/N^2.
            x = (float) (0.094864 - 2.08*rPrime);
            x = (x > 0) ? x : 0;
            estNe = (float) (0.308 + sqrt(x));
        } else {
            x = (float) (0.381924 - 5.24*rPrime);   // (0.618)^2 = 0.381924
            x = (x > 0) ? x : 0;
            estNe = (float) (0.618 + sqrt(x));
        }
    };
	estNe = estNe/(2*rPrime);
	return (estNe > infinite)? infinite : estNe;

}
//-------------------------------------------------------------------------
// added in July-Aug 2016
char JackSamp(int n, double rSmp[], float r, unsigned long long rCount[],
// param opened is temporarily added for checking: print to file checkR2
// char opened,
				float *lowr, float *highr, long *Jdegree)
// On input,
// * r is the interested value obtained from the whole sample set S with
//   n elements.
// * rSmp[k] is assumed to be the interested value obtained from Sk,
//   the sample set where individual (k+1)th is removed, k = 0, ..., n-1.
// * rCount[k] = 0 if rSmp[k] is not calculated for sample set Sk.
//   (When this function is called as part of Jackknife on samples, rCount[k]
//    represents the number of locus pairs that are eligible in the set Sk.)
// On output,
// * lowr and highr are the lower and upper bounds of the 95% confidence
//   interval for r.
// Function returns 0 if no process for confidence interval; 1 otherwise.
// The following algebraic identity is used:
// Suppose A is the average of an array x[] of n elements. Then
// Sum {(x[k] - A)^2} = Sum {x[k]^2} - (1/n)* (Sum {x[k]})^2
{
	int k, nJack;
	double j1, rTot, rSqTot, varJack;
	float rAve;
	float correction = 0.84;

// temporarily added for checking: print to file checkR2
// (The file checkR2.txt was opened, written and closed in Pair_Analysis.)
//FILE *checkR2 = NULL;
//if (opened == 1) checkR2 = fopen ("checkR2.txt", "a");

	for (k = 0, nJack = 0, rTot = 0, rSqTot = 0; k < n; k++) {
		if (rCount[k] > 0) {
			rTot += rSmp[k];
			rSqTot += (rSmp[k]*rSmp[k]);
			nJack++;
		}
	}
	if (nJack <= 0) return 0;
	rAve = rTot/nJack;
// temporarily added for checking: print to file checkR2
/*
if (checkR2 != NULL) {
fprintf(checkR2, "\nFor weighted r^2-values on the last column of the table:\n");
fprintf(checkR2, "* Sum of %d weighted r^2-values = %10.7f\n", nJack, rTot);
fprintf(checkR2, "* Average of %d weighted r^2-values = %10.7f\n", nJack, rAve);
}
//*/

	j1 = 1.0/((float) nJack);
	varJack = (nJack - 1)*j1 * (rSqTot - j1*rTot*rTot);
	correction *= correction;
// temporarily added for checking: print to file checkR2
/*
if (checkR2 != NULL) {
fprintf(checkR2,
"* Variance of r^2 = %11.4e, multiplied by factor (%f) = %11.4e\n",
varJack, correction, varJack*correction);
}
//*/

	varJack *= correction;	// empirical correction
	*Jdegree = JackKnifeInd(rAve, varJack);
//	float phi = varJack/(rAve*rAve);
//	*Jdegree = (long) floor(2.0/phi + 0.5);
//	if (*Jdegree <= 0) *Jdegree = 1;	// as a precaution against round-off

// temporarily added for checking: print to file checkR2
/*
if (checkR2 != NULL) {
fprintf(checkR2, "* Degree of freedom = %lu\n", *Jdegree);
}
//*/

	Confid95(*Jdegree, r, lowr, highr);

// temporarily added for checking: print to file checkR2
/*
if (checkR2 != NULL) {
fprintf(checkR2, "* Confidence interval for r^2: [%e, %e]\n", *lowr, *highr);
fflush(checkR2);
fclose(checkR2);
}
//*/

	return 1;
}

//-------------------------------------------------------------------------
// modified for Jackknife on samples - 2016
// Replace LDConfidInt (which was for Jackknife on loci, with parameter drift
// as an option for having CI with r2-drift, before translated to CI for Ne).
int LDConfidInt95 (float harmonic, int nfish, float wExpR2, float rB2WAve,
				double nIndSum, double r2WRemSmp[], unsigned long long rCount[],
// param opened is temporarily added for checking: print to file checkR2
//char opened,
				char modify, float *confidL, float *confidH,
				long *Jdegree, float infinite,
				char mating, int mode, char moreBurr, FILE *outBurr)
// mode = 0: parameter, 1: Jackknife.
{
	float lowR2, hiR2;
    float lowNe, hiNe;
	long indR2, indRdrift;
	float lowR2drift, hiR2drift;

	*Jdegree = 0;
	if (mode != 0) {	// do Jackknife
		if (JackSamp(nfish, r2WRemSmp, rB2WAve, rCount,
// param opened is temporarily added for checking: print to file checkR2
//opened,
		&lowR2, &hiR2, Jdegree) == 0)
		{
			printf("*** Jackknife on samples is not possible.\n");
			return 1;
		}
	} else {
		indR2 = (long) nIndSum;
		Confid95 (indR2, rB2WAve, &lowR2, &hiR2);
	}
	if (outBurr != NULL && moreBurr== 1 && nIndSum > 0 &&
	// add in Nov 2014: conditions on NOEXPLAIN (effected whether
	// Burrows output in multiple files or not)
		(NOEXPLAIN != 1)) {
		fprintf (outBurr, "\n");
		if (mode == 0) {
			fprintf (outBurr, "# Parametric CI for r^2:      ");
			fprintf (outBurr, "%10.6f%12.6f", lowR2, hiR2);
			lowNe = LD_Ne (harmonic, hiR2 - wExpR2, mating, infinite);
			hiNe = LD_Ne (harmonic, lowR2 - wExpR2, mating, infinite);
			fprintf (outBurr, "   >>> CI for Ne:");
			if (lowNe < 0 || lowNe > infinite)
				fprintf (outBurr, "%11s", "infinite");
			else fprintf (outBurr, "%11.1f", lowNe);
			if (hiNe < 0 || hiNe > infinite)
				fprintf (outBurr, "%11s\n", "infinite");
			else fprintf (outBurr, "%11.1f\n", hiNe);
		} else {
			fprintf (outBurr, "# Jackknife CI for r^2:       ");
			fprintf (outBurr, "%10.6f%12.6f", lowR2, hiR2);
			lowNe = LD_Ne (harmonic, hiR2 - wExpR2, mating, infinite);
			hiNe = LD_Ne (harmonic, lowR2 - wExpR2, mating, infinite);
			fprintf (outBurr, "   >>> CI for Ne:");
			if (lowNe < 0 || lowNe > infinite)
				fprintf (outBurr, "%11s", "infinite");
			else fprintf (outBurr, "%11.1f", lowNe);
			if (hiNe < 0 || hiNe > infinite)
				fprintf (outBurr, "%11s\n", "infinite");
			else fprintf (outBurr, "%11.1f\n", hiNe);
		}
		fflush(outBurr);
	}
	lowR2drift = lowR2 - wExpR2;
	hiR2drift = hiR2 - wExpR2;
	lowNe = LD_Ne (harmonic, hiR2drift, mating, infinite);
	hiNe = LD_Ne (harmonic, lowR2drift, mating, infinite);
	if (hiNe > infinite || hiNe <= 0) hiNe = infinite;
	if (modify == 1) {	// lower confidL to lowNe if it is bigger,
		// and raise confidH to hiNe if it is smaller
		if (lowNe < *confidL) *confidL = lowNe;
		if (hiNe > *confidH) *confidH = hiNe;
	} else {	// new values for confidL and confidH
		*confidL = lowNe;
		*confidH = hiNe;
	}
	return 0;
}

// --------------------------------------------------------------------------

//int NeAdjustedTmp (FILE *rAveTemp, FILE *weighFile,
// Remove parameter weighFile
int NeAdjustedTmp (FILE *rAveTemp,
// param opened is temporarily added for checking: print to file checkR2
//char opened,
				long nBurrVal, float harmonic, char matingMod,
				float infinite, float *adjNe, float *r2driftAve,
				float *totW, float *totR2, float *totRdrift, float *expR2,
				float *rBurrAve)
// Calculate Ne based on arrays r-values, number of Ind alleles, sample size,
// Parameter "adjNe" at input is the tentative Ne, output: new value for Ne
// If <= 0 or too big: no adjusting.
// Reassign the weights in the array pairWt
// Return 0 if no change, 1 otherwise
{
	long ind;
	float indAlle;
	float weight, nsamp, a, b, r2, r2drift;
	double bigW, bigR2, bigRdrift, r2ExpW;
	long count=0;
	a = (*adjNe)*3;
	// to avoid overflow, quit if a is 0 or too big
	if (a >= infinite || a <= 0) return 0;
//	if (rAveTemp == NULL || weighFile == NULL) return 0;
	if (rAveTemp == NULL) return 0;


// temporarily added for checking: print to file checkR2
// (The file checkR2.txt was opened, written and closed in Pair_Analysis.)
//FILE *checkR2 = NULL;
//if (opened == 1) checkR2 = fopen ("checkR2.txt", "a");

	r2ExpW = 0;
	bigW = 0;
	bigR2 = 0;
	bigRdrift = 0;
//	rewind (rAveTemp);	// this is done in the calling function
	printf ("     Initial estimate of Ne: %12.1f\n", *adjNe);
	for (ind = 0; ind < nBurrVal; ind++) {
		// retrieve values: start with the product of ind. alleles,
		// then sample having data at pair, Rdrift-value.
		// The final weight at pair will be adjusted
		indAlle = 0;
		fread (&indAlle, sizeof(float), 1, rAveTemp);
		count++;
		if (indAlle < 0.5F) break;
		fread (&nsamp, sizeof(float), 1, rAveTemp);
		fread (&r2, sizeof(float), 1, rAveTemp);
//		fwrite (&r2, sizeof(float), 1, weighFile);
		fread (&r2drift, sizeof(float), 1, rAveTemp);
//		fwrite (&r2drift, sizeof(float), 1, weighFile);
// NOTE: the weight may be adjusted here before multiplying with (nsamp)^2,
//
		weight = indAlle*(nsamp*nsamp);
// the condition is actually unnecessary; just put here if we want to
// proceed without exiting by this condition
		if (a < infinite && a > 0) {
			b = a + nsamp;
			b *= b;
			weight /= b;
		}
		bigR2 += (r2*weight);
		bigW += weight;
		bigRdrift += (r2drift * weight);

// if try to update the weight recorded in rAveTemp, it seems that the pointer
// will go to the end of rAveTemp when this is run in Mac.
// Since we will not use rAveTemp, instead, replace by weighFile, so need
// to read one more float to position the pointer to the next data set
// (Here, use b as a dummy variable.)
//this line is not good:	fwrite (&weight, sizeof(float), 1, rAveTemp);
		fread (&b, sizeof(float), 1, rAveTemp);
//		fwrite (&weight, sizeof(float), 1, weighFile);
		r2ExpW += (ExpR2Samp (nsamp) * weight);
	}
	*totR2 = (float) bigR2;
	*totW = (float) bigW;
	*totRdrift = (float) bigRdrift;
	bigR2 /= bigW;
	bigRdrift /= bigW;
	r2ExpW /= bigW;
	*r2driftAve = (float) bigRdrift;		// weighted average of r2-drift.
	*rBurrAve = (float) bigR2;				// weighted average of r2.
	*expR2 = (float) r2ExpW;
	*adjNe = LD_Ne(harmonic, *r2driftAve, matingMod, infinite);
	printf ("     Final estimate of Ne: %14.1f\n", *adjNe);

// temporarily added for checking: print to file checkR2
/*
if (checkR2 != NULL) {
fprintf(checkR2, "\n");
fprintf(checkR2, "On the whole sample set, after reweighting because of missing data:\n");
fprintf(checkR2, "* Weighted average r^2 = %10.7f\n", *rBurrAve);
fprintf(checkR2, "* Weighted average of expected r^2-sample = %10.7f\n", *expR2);
fprintf(checkR2, "* Weighted average r^2-drift = %10.7f\n", *r2driftAve);
fflush(checkR2);
fclose(checkR2);
}
//*/

	return 1;
}


// --------------------------------------------------------------------------

int NeAdjustedArr (float *pairWt, float *rB2, float *rBdrift, float *prodInd,
				float *sampCount, long nBurrVal, float harmonic, char matingMod,
				float infinite, float *adjNe, float *r2driftAve,
				float *totW, float *totR2, float *totRdrift, float *expR2,
				float *rBurrAve)
// Calculate Ne based on arrays r-values, number of Ind alleles, sample size,
// Parameter "adjNe" at input is the tentative Ne, output: new value for Ne
// If <= 0 or too big: no adjusting.
// Reassign the weights in the array pairWt
// Return 0 if no change, 1 otherwise
{
	long ind;
	float weight, nsamp, a, b, r2, r2drift, r2ExpW;
	a = (*adjNe)*3;
	// to avoid overflow, quit if a is 0 or too big
	if (a >= infinite || a <= 0) return 0;
	*totRdrift = 0;
	*totR2 = 0;
	*totW = 0;
	r2ExpW = 0;

	printf ("     Initial estimate of Ne: %12.1f\n", *adjNe);
	for (ind = 0; ind < nBurrVal; ind++) {
		// retrieve values: start with the product of ind. alleles,
		// then sample having data at pair, Rdrift-value.
		// The final weight at pair will be adjusted
		weight = *(prodInd+ind);
		// the product of indep. alleles should be a whole number, which was
		// changed to float, so if < 0.5, it means that it is zero,
		// therefore all values are read.
		// (the number of values nBurrVal can be an overestimate)
		if (weight < 0.5F) break;
		nsamp = *(sampCount+ind);
		r2 = *(rB2+ind);
		r2drift = *(rBdrift+ind);
// NOTE: the weight may be adjusted here before multiplying with (nsamp)^2,
//
		weight *= (nsamp*nsamp);
		b = a + nsamp;
		b *= b;
		weight /= b;
		*totR2 += (r2*weight);
		*totRdrift += (r2drift * weight);
		*totW += weight;
		*(pairWt+ind) = weight;
		r2ExpW += (ExpR2Samp (nsamp) * weight);
	}
	r2drift = (*totRdrift)/(*totW);		// weighted average of r2-drift.
	r2 = (*totR2)/(*totW);				// weighted average of r2.
	*r2driftAve = r2drift;
	*expR2 = r2ExpW/(*totW);			// weighted average of expR2-sample
	*rBurrAve = r2;
//	*rBurrAve = (*expR2) + r2drift;	// weighted average of r^2
	*adjNe = LD_Ne(harmonic, r2drift, matingMod, infinite);
	printf ("     Final estimate of Ne: %14.1f\n", *adjNe);
	return 1;
}

// --------------------------------------------------------------------------

// ---------------------------------------------------------------------------

// Jackknife on loci is blocked out:
/* --------------------------------
// for Jackknife on loci: JackKnifeInd and LDJackKnifeInd
//------------------------------------------------------------------------

// Adjusted Variance in Jackknife
void LDJackKnifeInd(float r2Ave, float r2driftAve, float *pairWt, float *rB2,
					float *rBdrift,
//					float *varJR2, float *varJRdrift,
					unsigned long long nBurrVal, float totW, float totR2,
					float totRdrift, char tmpUsed,
					FILE *weighFile, long *indR2, long *indDrift)
// Adjusted Variance in Jackknife
{
//
// Doing Jackknife:
// * for weighted r: find R(k) = weighted mean of weighted r after deleting r
//   at record k, then find mean R(:) of all R(k). Then take the sum of
//   all (R(k) - R(:))^2, multiply by (N-1)/N.
//   Let wk = weight for individual rk, from array pairWt.
//   Let rBar = weighted mean of all rk.
//   Assume W = sum of all weights wk, k = 1, ..., N. (as paramater totW)
//   Then
//   rBar = sum of all wk*rk, divided by W, is the overall weighted mean
//          (can be introduced as an input parameter, here, it is r2driftAve).
//   R(k) = sum of wj*rj for all j /= k, divided by sum(wj, j/=k) (= W - wk)
//        = (W*rBar - wk*rk)/(W - wk).
//   Suppose R(k) is weighted by sum(wj, j/=k) = W - wk.
//   Then the weighted average of R(k) is
//   Rw(:) = Sum{k}[(W-wk)*R(k)] / Sum{k}(W-wk)
//         = Sum{k}[W*rBar - wk*rk] / (N-1)*W = [N*W*rBar - rBar*W] / (N-1)*W
//         = rBar.
//   R(k) - Rw(:) = (W*rBar - wk*rk)/(W - wk) - rBar
//                = [wk /(W - wk)]* (rBar - rk)
//   Now, take the sum over k of [R(k) - Rw(:)]^2, then multiply by (N-1)/N
//   to get the Adjusted Variance.
//
	unsigned long long k;
	float varR2, varRdrift, factor, wk, r2, rdrift, r2W, rdriftW;
	*indR2 = 0;
	*indDrift = 0;
	varR2 = 0;
	varRdrift = 0;
	if (tmpUsed != 0) rewind (weighFile);	// the weighFile records pairs of
						// r-drift and its weight for each pair of loci when
						// Ne is recalculated with revised weights.

	for (k=0; k<nBurrVal; k++) {
		if (tmpUsed != 0) {
			fread (&r2, sizeof(float), 1, weighFile);
			fread (&rdrift, sizeof(float), 1, weighFile);
			fread (&wk, sizeof(float), 1, weighFile);
		} else {
			wk = *(pairWt+k);
			r2 = *(rB2+k);
			rdrift = *(rBdrift+k);
		}
		factor = wk /(totW-wk);
		r2W = (r2-r2Ave)*factor;
		rdriftW = (rdrift-r2driftAve)*factor;
		varRdrift += (rdriftW*rdriftW);
		varR2 += (r2W*r2W);
	}
    if (nBurrVal > 0) {
		varR2 = ((nBurrVal - 1)*varR2) / ((float) nBurrVal);
		varRdrift = ((nBurrVal - 1)*varRdrift) / ((float) nBurrVal);
	}
//	*varJR2 = varR2; *varJRdrift = varRdrift;
	*indR2 = JackKnifeInd(r2Ave, varR2);
	*indDrift = JackKnifeInd(r2driftAve, varRdrift);
}
// ---------------------------------------------------------------------------

int LDConfidInt (char drift, unsigned long long nBurrAve, float harmonic,
				float wExpR2, float rB2WAve, float r2driftAve, double nIndSum,
				char modify, float *confidL, float *confidH, float infinite,
				char mating, int mode, float *pairWt, float *rB2, float *rBdrift,
				char tmpUsed, FILE *weighFile, float totW,
				float totR2, float totRdrift, char moreBurr, FILE *outBurr)
// mode = 0: parameter, 1: Jackknife.
{
	float lowR2drift, hiR2drift;
	float lowR2, hiR2;
    float lowNe, hiNe;
	long indR2, indRdrift;
//	float varR2, varRdrift;

	indR2 = (long) nIndSum;
	indRdrift = (long) nIndSum;
	if (mode != 0) {	// do Jackknife
		LDJackKnifeInd(rB2WAve, r2driftAve, pairWt, rB2, rBdrift,
					nBurrAve, totW, totR2, totRdrift,
					tmpUsed, weighFile, &indR2, &indRdrift);
	}
	Confid95 (indRdrift, r2driftAve, &lowR2drift, &hiR2drift);
	Confid95 (indR2, rB2WAve, &lowR2, &hiR2);
	if (outBurr != NULL && moreBurr== 1 && nBurrAve > 0 &&
	// add in Nov 2014: conditions on NOEXPLAIN (effected whether
	// Burrows output in multiple files or not)
		(NOEXPLAIN != 1)) {
		fprintf (outBurr, "\n");
		if (mode == 0) {
			fprintf (outBurr, "# Parametric CI for r^2:      ");
			fprintf (outBurr, "%10.6f%12.6f", lowR2, hiR2);
			lowNe = LD_Ne (harmonic, hiR2 - wExpR2, mating, infinite);
			hiNe = LD_Ne (harmonic, lowR2 - wExpR2, mating, infinite);
			fprintf (outBurr, "   >>> CI for Ne:");
			if (lowNe < 0 || lowNe > infinite)
				fprintf (outBurr, "%11s", "infinite");
			else fprintf (outBurr, "%11.1f", lowNe);
			if (hiNe < 0 || hiNe > infinite)
				fprintf (outBurr, "%11s\n", "infinite");
			else fprintf (outBurr, "%11.1f\n", hiNe);
			fprintf (outBurr, "# Parametric CI for r^2-drift:");
			fprintf (outBurr, "%10.6f%12.6f", lowR2drift, hiR2drift);
			lowNe = LD_Ne (harmonic, hiR2drift, mating, infinite);
			hiNe = LD_Ne (harmonic, lowR2drift, mating, infinite);
			fprintf (outBurr, "   >>> CI for Ne:");
			if (lowNe < 0 || lowNe > infinite)
				fprintf (outBurr, "%11s", "infinite");
			else fprintf (outBurr, "%11.1f", lowNe);
			if (hiNe < 0 || hiNe > infinite)
				fprintf (outBurr, "%11s\n", "infinite");
			else fprintf (outBurr, "%11.1f\n", hiNe);
		} else {
//			fprintf (outBurr, "Jackknife Variance for r^2 = %12.4e\n", varR2);
			fprintf (outBurr, "# Jackknife CI for r^2:       ");
			fprintf (outBurr, "%10.6f%12.6f", lowR2, hiR2);
			lowNe = LD_Ne (harmonic, hiR2 - wExpR2, mating, infinite);
			hiNe = LD_Ne (harmonic, lowR2 - wExpR2, mating, infinite);
			fprintf (outBurr, "   >>> CI for Ne:");
			if (lowNe < 0 || lowNe > infinite)
				fprintf (outBurr, "%11s", "infinite");
			else fprintf (outBurr, "%11.1f", lowNe);
			if (hiNe < 0 || hiNe > infinite)
				fprintf (outBurr, "%11s\n", "infinite");
			else fprintf (outBurr, "%11.1f\n", hiNe);
//			fprintf (outBurr, "Jackknife Variance for r^2-drift = %12.4e\n", varRdrift);
			fprintf (outBurr, "# Jackknife CI for r^2-drift: ");
			fprintf (outBurr, "%10.6f%12.6f", lowR2drift, hiR2drift);
			lowNe = LD_Ne (harmonic, hiR2drift, mating, infinite);
			hiNe = LD_Ne (harmonic, lowR2drift, mating, infinite);
			fprintf (outBurr, "   >>> CI for Ne:");
			if (lowNe < 0 || lowNe > infinite)
				fprintf (outBurr, "%11s", "infinite");
			else fprintf (outBurr, "%11.1f", lowNe);
			if (hiNe < 0 || hiNe > infinite)
				fprintf (outBurr, "%11s\n", "infinite");
			else fprintf (outBurr, "%11.1f\n", hiNe);
		}
		fflush(outBurr);
	}
	if (drift != 1) {	// use r^2-CIs to generate CI for Ne
		lowR2drift = lowR2 - wExpR2;
		hiR2drift = hiR2 - wExpR2;
	}
	lowNe = LD_Ne (harmonic, hiR2drift, mating, infinite);
	hiNe = LD_Ne (harmonic, lowR2drift, mating, infinite);
	if (hiNe > infinite || hiNe <= 0) hiNe = infinite;
	if (modify == 1) {	// lower confidL to lowNe if it is bigger,
		// and raise confidH to hiNe if it is smaller
		if (lowNe < *confidL) *confidL = lowNe;
		if (hiNe > *confidH) *confidH = hiNe;
	} else {	// new values for confidL and confidH
		*confidL = lowNe;
		*confidH = hiNe;
	}
	return 0;
}
//*/

// ----- End of common functions for use -------

// --------------------------------------------------------------------------


// ---------------------------------------------------------------------------
// Rewritten in March 2016, so this only needed to be called once for a pair
// of loci (p1, p2).
// Return at pair of loci (p1, p2):
// * Number of samples *nSamp among "nfish" samples having data at both loci.
// * Number of independent alleles nInd1, nInd2, at loci p1, p2 with freq.
//   calculated against samples having data at both loci.
// * Values for arrays mValp1, freqp1, homop1, where nEff (nEff1, nEff2) will
//   be their effective size (do not access beyond the first nEff elements),
//   mValp1[0] < mValp1[1] < ... < mValp1[nEff-1] are mobility values,
//   freqp1[0], ..., freqp1[nEff-1] are frequencies corresponding to those
//   mobilities, and homop1 are for freq of homozygotes.
// (When nEff1 or nEff2 = 0, its product = 0, then Burrows_Calcul will be
//  skipped right after calling this function.)
// * Similarly for mValp2, freqp2, homop2.
// (The arrays are allocated once for all pairs of loci to save running time)
// On input, nMp1, nMp2 are the number of nodes in allep1, allelp2,
// the total number of allele mobilities at loci p1, p2.
// On input, missing = 1 when the input file has missing data

// Dec 2016:
// ---------
// Modify to deal with the case that only singleton alleles are dropped
// when cutoff takes a special value.

// Add parameter *cutoffRev to calculate the real cutoff value when
// cutoff parameter is the special value for dropping singlton alleles.
// In this function, alleles at the pairs are accepted by this reassigned
// cutoff value, and the parameter *cutoffRev will be used by
// Burrows_Calcul to proceed the jackknife algorithm. This parameter
// is only used in Burrows_Calcul, and its value is varied by locus pairs.

void IndAlle2 (int **p1Gen, int** p2Gen, int *noDatFish,
				FISHPTR fishp1, FISHPTR fishp2, float cutoff, int nfish,
				ALLEPTR allep1, ALLEPTR allep2, int nMp1, int nMp2,
				int *nEff1, int *mValp1, int *nEff2, int *mValp2, float *nSamp,
				float *freqp1, float *homop1, float *freqp2, float *homop2,
				int *nInd1, int *nInd2, char missing,
    			// add fmin in Nov 2014 for lowest freq.
				float *fminp1, float *fminp2, float *cutoffRev)
{

    // add in Nov 2014 for lowest freq.
    *fminp1 = 1.0;
    *fminp2 = 1.0;
    FISHPTR fish1, fish2;
	int misdat = 0;
	int i, k, m, n;
	int totAlle = 2*nfish;
	int nzero, ndrop, mhomo, mcount;
	float x;
	int noDat, noDat1, noDat2;
    ALLEPTR curr;
	for (k = 0; k < nfish; k++) {
		noDatFish[k] = 0;	// assuming no data missing at sample (k+1)
	}
// In the comments, locus p means the (p+1)th locus
// --------------------------------------------------------------------------
	int homo1, homo2;	// count homoz. in samples missing data at one locus
	fish1 = fishp1; fish2 = fishp2;
	for (k = 0, homo1 = 0, homo2 = 0; (fish1 != NULL) && (fish2 != NULL);
				fish1 = fish1->next, fish2 = fish2->next, k++)
	{
		if (missing != 0) {	// there are missing data in input, so check here
			noDat1 = (fish1->gene[0] == 0)? 1: 0;
			noDat2 = (fish2->gene[0] == 0)? 2: 0;
			noDat = noDat1 + noDat2;	// when = 3: no data at both loci
			noDatFish[k] = noDat;	// = 1 or 2: missing data at one locus
			if (noDat > 0) {
				misdat++;
				if ((noDat==1) && (fish2->gene[0]==fish2->gene[1])) homo2++;
				if ((noDat==2) && (fish1->gene[0]==fish1->gene[1])) homo1++;
			}
		}
		for (i = 0; i < 2; i++) {
			p1Gen[k][i] = fish1->gene[i];
			p2Gen[k][i] = fish2->gene[i];
		}
	}
	*nSamp = (float) nfish - misdat;
// --------------------------------------------------------------------------
// Dec 2016: to deal with value of cutoff for dropping singleton alleles only
	if (*nSamp > 0 && cutoff > 0 &&  cutoff <= PCRITX) {
		// set cutoff to drop only singleton
		// let n = *nSamp. Set x, (1/2n) < 1/x <= 1/n, so that only singleton
		// allele (freq = 1/2n) is rejected (if not singleton, freq >= 1/n).
		// when one sample having data is removed from the set, the remaining
		// has (n-1). Singleton alleles will have freq = 1/(2(n-1)), others
		// have freq at least 1/(n-1). Thus, we wish to have
		// 1/2(n-1) < 1/x <= 1/(n-1). Combining, x will be set so that
		// 				(n-1) < n <= x < 2(n-1) < 2n.
		// Only the two middle inequalities are needed. They cannot be
		// satisfied concurrently when n <= 2. Obvious when n = 1.
		// When n = 2, no x such that 2 = n <= x < 2(n-1) = 2
		// However,
		// * if n = 1, the sample is either a homozygote or both alleles are
		//   singletons. We can set cutoff > 1/2 so that *nEff1 or *nEff2 = 0
		//   (*nEff1, *nEff2 are # of alleles >= cutoff and <= 1-cutoff)
		// * if n = 2, the cutoff value that works for the whole set (n=2)
		//   is in the interval (1/4, 1/2], then it will allow singleton
		//   when one sample is taken out (one sample left). However, Burrows
		//   correlation r will be 0 then, so it does not effect much.
		//   Hence we can set cuttoff = 1 when n = 1 and 1/3 when n = 2,
		//   corresponding to x = 1 and x = 3, so set x = 2n - 1.
		// When n > 2, we have n <= 2n - 5/2 < 2(n-1).
		//
		x = 2*(*nSamp) - 1;
		if (*nSamp > 2) x -= 1.5;
		cutoff = 1.0/x;
	}
	*cutoffRev = cutoff;
// cutoff is used through the rest, cutoffRev is assigned by this function
// to be used in the calling function Burrows_Calcul.
// --------------------------------------------------------------------------

	// Index starts at 0, so loci p1, p2 mean loci (p1+1), (p2+1) by counting
	// We have 2 arrays p1Gen, p2Gen for loci p1, p2. Each element of an array
	// is an array of 2 genes from the sample corresponding to that element.
	// Array noDatFish is to show which element of the two arrays has no data
	// If there is no missing data between two loci,
	if (misdat == 0) {
		curr = allep1;
		for (n=0, ndrop=0; curr != NULL; curr=curr->next)
		{
        	m = curr->mValue;
			x = curr->freq;
			mhomo = curr->homozyg;
			if (x < cutoff) ndrop++;	// if cutoff = 0, ndrop stays unchanged
										// (also, there is no x = 0)
        // just in case cutoff = 0, we don't want to include x = 1,
		// so the condition x < 1 is needed
			else if (x < 1 && x <= 1 - cutoff) {
        // if there is no such x, then n stays 0.
        // There is at most one x > 1 - cutoff (provided cutoff < 0.5).
        // If there is one allele with > 1 - cutoff, then all others have
        // x < cutoff; therefore
        // nMp1 - 1 = (number of alleles at locus p1 with freq. < cutoff);
        // or nMp1 - 1 = ndrop.
        // Suppose there is such x violating the condition stated:
        //     * If cutoff = 0, then x = 1, ndrop still = 0.
        //     * If cutoff > 0, x could be 1 or just merely > 1 - cutoff.
        //       In case x < 1, then ndrop > 0; otherwise (x = 1), ndrop = 0
        // After this loop, nInd1 = nMp1 - ndrop. If ndrop=0, nInd1 = nMp1-1.
        // Since this ELSE IF is checked only when x >= cutoff, n will
        // then stays 0 if there is such dominant allele A.
        // n will be the number of alleles between cutoff and (1 - cutoff)
        // to be used in the pairwise comparison for locus pair p1, p2.
        // At the end, n > 0: no dominant allele
				*(mValp1+n) = m;
				*(freqp1+n) = x;
				*(homop1+n) = (float) mhomo/(*nSamp);
				// add in Nov 2014:
				if (x < *fminp1) *fminp1 = x;
				n++;
			}
		}
		*nEff1 = n;
    	if (ndrop > 0)	// ndrop = # alleles having freq < cutoff
        	nMp1 -= ndrop;
    	else					// all freq >= cutoff
    	// Number of freq needed to know is one less, because the sum = 1
        	nMp1--;
	// add in Nov 2014: The next line might be redundant, just for assurance
		if (nMp1 == 0) *fminp1 = 0;
		if (ndrop > 0 && nMp1 == 1) *fminp1 = 1-(*fminp1);

    	if (n == 0) nMp1 = 0;	// if there is a dominant allele (freq > 1-cutoff),
                            // then all others have frequencies either = 0
							// or < cutoff, so n = 0, but nMp1 calculated = 1
		*nInd1 = nMp1;
		// Now, do the same for allele list at locus p2:
		curr = allep2;
		for (n=0, ndrop=0; curr != NULL; curr=curr->next)
		{
        	m = curr->mValue;
			x = curr->freq;
			mhomo = curr->homozyg;
			if (x < cutoff) ndrop++;
			else if (x < 1 && x <= 1 - cutoff) {
				*(mValp2+n) = m;
				*(freqp2+n) = x;
				*(homop2+n) = (float) mhomo/(*nSamp);
				if (x < *fminp2) *fminp2 = x;
				n++;
			}
		}
		*nEff2 = n;
    	if (ndrop > 0)	// ndrop = # alleles having freq < cutoff
        	nMp2 -= ndrop;
    	else nMp2--;
		if (nMp2 == 0) *fminp2 = 0;
		if (ndrop > 0 && nMp2 == 1) *fminp2 = 1-(*fminp2);
    	if (n == 0) nMp2 = 0;
    	*nInd2 = nMp2;
		return;
	}
	// Now for the case that there are missing data at either locus p1 or p2
	totAlle -= 2*misdat;	// total number of alleles at each locus, in
							// samples having data at both loci.
	int mLeft, aLeft;
	// mLeft: number of alleles not counted at locus p1
	// aLeft: number of alleles left from allep1 list
	// Recalculate frequencies at locus p1 based on samp having data at both
    curr = allep1;
	for (n = 0, nzero = 0, ndrop = 0, mLeft = nMp1, aLeft = totAlle;
				(aLeft > 0) && (curr != NULL); curr = curr->next)
	{
        m = curr->mValue;
		mhomo = curr->homozyg;
		if (mLeft == 1) {
			mhomo -= homo1;
			mcount = aLeft;
			aLeft = 0;
			mLeft--;
		} else {
			mcount = curr->copy;
			for (i = 0; i < nfish; i++)	// [i] = (i+1)th samp, (p) = locus p
			{
				if (noDatFish[i] == 2) {	// [i] has data at (p1), not (p2)
					k = Count(p1Gen[i], m); // k = # allele m at [i], (p1)
					mcount -= k;	// against samp having data at both (p)
					if (k == 2) {
						mhomo--;
						homo1--;	// homo1 is the number of homo at (p1) in
					}				// samples that have missing data at (p2)
				}
			}
			aLeft -= mcount;
			mLeft--;
// When aLeft = 0, it is possible that there are still allele mobilities
// appeared only in samples having data at locus (p1) but no data at locus
// (p2), i.e., mLeft may be nonzero. In such case, those allele mobilities
// have zero frequency in set of samples having data at both loci.
 			if (aLeft == 0 && mLeft > 0) nzero += mLeft;
		}
        // x is freq of allele m calculated against nSamp fish having
        // data at both loci p1, p2:
        x = (*nSamp > 0)? (float) mcount/(2*(*nSamp)): 0;
        if (x == 0) nzero++;	// nzero is to take care of the next condition
        else {
			if (x < cutoff) ndrop++;	// if cutoff = 0, ndrop stays unchanged,
										// so nzero is needed
        // just in case cutoff = 0, we don't want to include x = 1,
		// so the condition x<1 is needed
			else if (x < 1 && x <= 1 - cutoff) {
        // if there is no such x, then n stays 0.
        // There is at most one > 1 - cutoff (provided cutoff < 0.5).
        // Thus, if there is one such allele, say A (its x(A) > 1 - cutoff),
        // then all others have x < cutoff; therefore nMp1 - 1
        // is the number of alleles at locus p1+1 having recalculated
        // freq x = 0 or < cutoff; so, nMp1 - 1 = ndrop + nzero.
        // If that x(A) < 1, then ndrop > 0; otherwise ndrop = 0.
        // Immediately after this loop, nInd1 = nMp1 - nzero.
        // Then x(A) < 1 => ndrop > 0 => nInd1 = nMp1 - nzero - ndrop = 1,
        // x(A) = 1 => ndrop = 0 => nInd1 = nMp1 - nzero - 1 = 0.
        // Since this ELSE IF is checked only when x >= cutoff, n will
        // then stays 0 if there is such dominant allele A.
        // n will be the number of alleles between cutoff and (1 - cutoff)
        // to be used in the pairwise comparison for locus pair p1, p2.
        // At the end, n > 0: no dominant allele
				*(mValp1+n) = m;
				*(freqp1+n) = x;
				*(homop1+n) = (float) mhomo/(*nSamp);
				// add in Nov 2014:
				if (x < *fminp1) *fminp1 = x;
				n++;
			}
		}
	}
	*nEff1 = n;
   	nMp1 -= nzero;	// exclude alleles that are absent when only samples
					// having data at both loci p1+1, p2+1 are counted
   	if (ndrop > 0)			// ndrop = # alleles having freq < cutoff
       	nMp1 -= ndrop;
   	else					// all freq >= cutoff
   	// Number of freq needed to know is one less, because the sum = 1
       	nMp1--;
	if (nMp1 == 0) *fminp1 = 0;
	// The next line might be redundant, just for assurance
	if (ndrop > 0 && nMp1 == 1) *fminp1 = 1-(*fminp1);
   	if (n == 0) nMp1 = 0;
   	*nInd1 = nMp1;
    // done with alleles at locus p1

	// Now do the same for locus p2
    curr = allep2;

	for (n = 0, nzero = 0, ndrop = 0, mLeft = nMp2, aLeft = totAlle;
				(aLeft > 0) && (curr != NULL); curr = curr->next)
	{
        m = curr->mValue;
		mhomo = curr->homozyg;
		if (mLeft == 1) {
			mhomo -= homo2;
			mcount = aLeft;
			aLeft = 0;
			mLeft--;
		} else {
			mcount = curr->copy;
			for (i = 0; i < nfish; i++)	// [i] = (i+1)th samp, (p) = locus p
			{
				if (noDatFish[i] == 1) {	// [i] has data at (p2), not (p1)
					k = Count(p2Gen[i], m); // k = # allele m at [i], (p2)
					mcount -= k;	// against samp having data at both (p)
					if (k == 2) {
						mhomo--;
						homo2--;	// homo2 is the number of homo at (p2) in
					}				// samples that have missing data at (p1)
				}
			}
			mLeft--;
			aLeft -= mcount;
// When aLeft = 0, it is possible that there are still allele mobilities
// appeared only in samples having data at locus (p2) but no data at locus
// (p1), i.e., mLeft may be nonzero. In such case, those allele mobilities
// have zero frequency in set of samples having data at both loci.
			if (aLeft == 0 && mLeft > 0) nzero += mLeft;
		}

        // x is freq of allele m calculated against nSamp fish having
        // data at both loci p1, p2:
        x = (*nSamp > 0)? (float) mcount/(2*(*nSamp)): 0;
        if (x == 0) nzero++;	// nzero is to take care of the next condition
        else {
			if (x < cutoff) ndrop++;
			else if (x < 1 && x <= 1 - cutoff) {
				*(mValp2+n) = m;
				*(freqp2+n) = x;
				*(homop2+n) = (float) mhomo/(*nSamp);
				if (x < *fminp2) *fminp2 = x;
				n++;
			}
		}
	}
	*nEff2 = n;
   	nMp2 -= nzero;
   	if (ndrop > 0) nMp2 -= ndrop;
   	else nMp2--;
	if (nMp2 == 0) *fminp2 = 0;
	// The next line might be redundant, just for assurance
	if (ndrop > 0 && nMp2 == 1) *fminp2 = 1-(*fminp2);
   	if (n == 0) nMp2 = 0;
   	*nInd2 = nMp2;
	// done with alleles at locus p2
}

// ---------------------------------------------------------------------------
void AlleInSamp (int nfish, int m, int **pGen, int *noDatFish, char *countm)
{
	int k;
	for (k = 0; k< nfish; k++) {
		if (noDatFish[k] > 0) countm[k] = 0;
		else countm[k] = Count(pGen[k], m);
	}
}


// --------------------------------------------------------------------------

void Burrows_Delta (float f1, float f2, float x, float y,
					float nSamp, int nfish,
					float *dBur, float *rBur, float *rBur2,
					float *pSum, char *countm1, char *countm2)
{
// Inputs to the function:
// countm1, countm2 are two arrays holding the number of alleles, say m1, m2,
// at two loci; e.g., countm1[i] is the number of allele m1 at sample[i] at
// locus p1.
// f1, f2 are frequencies of m1, m2,
// Let f = frequency of m, and h = frequency of homozygote mm,
// then z = f(1-f) + h - f^2 = f(1-2f) + h. Here, x and y on input should
// be z, calculated at locus 1 and locus 2, where f replaced by f1, f2,
// and h replaced by h1, h2, respectively.
// (h1 is frequency of homozygotes [m1,m1] at locus 1,
// h2 is frequency of homozygotes [m2,m2] at locus 2.)
// Those h1, h2 are not in input list, replaced by x and y instead
// (to avoid repeating the expression for x or y when one h is not changed).
// Must have x, y > 0 in input (otherwise, Burrows coeffs are assigned 0).
// nSamp is the number of individuals having data at both loci,

	int i, countM;
	*dBur = 0;

    // countM is incremented if there are alleles m1 at locus p1, allele m2
    // at locus p2. The number of m1 at locus p1 is countm1, of m2 at locus
    // p2 is countm2. countM is incremented by:
    //    * 1 if heterozygotes at both loci,
    //    * 2 if hetero at one and homo at the other
    //    * 4 if homo at both.
	for (countM = 0, i = 0; i < nfish; i++)
			countM += (countm1[i]*countm2[i]);
	*pSum = (float) countM;
	if (nSamp > 0) *dBur = *pSum /((float) 2.0*nSamp) - 2.0*f1*f2;
	// (unbias) adjusting factor: nSamp/(nSamp-1)
	if (nSamp > 1) *dBur *= (nSamp/(nSamp -(float) 1.0));
//	Note: The condition x,y > 0 must hold before this function is called
	*rBur = (*dBur)/sqrt(x*y);
   	*rBur2 = (*rBur)*(*rBur);
	// although absolute value of rBur is at most 1 if not for the factor
	// nSamp/(nSamp-1) applied to "dBur" above, so we bring back to 1
	if (*rBur2 > 1.0) *rBur2 = 1.0;

}


// --------------------------------------------------------------------------
// This function is for determining if an allele is eligible after one
// sample (having data) is removed
// On input:
// remv = number of the allele (of interest) counted in the removed sample,
// f = frequency of that allele before the sample being removed
// On output,
// *fx = frequency of the interested allele after the sample is removed,
// Return 1 when allele is rejected by freq. restriction cutoff,
// return 0 otherwise.
char Rejected (float cutoff, float nSamp, float f, int remv, float *fx)
{
	float val = 2*nSamp;
// totAlle = total number of alleles after one sample is removed,
	float totAlle = val - 2;
	if (totAlle < 0.5) return 1;	// => totAlle = 0: no sample!
	val *= f;
	val -= remv;	// number of interested alleles are left;
					// theoretically, val is a whole number
	*fx = val/totAlle;
	if (cutoff == 0) {
		if (val < 0.5) *fx = 0;
		else if (val > totAlle - 0.5) *fx = 1;
		if (*fx == 0 || *fx == 1) return 1;
	} else {
		if (*fx < cutoff || *fx > (1.0 - cutoff)) return 1;
	}
	return 0;
}

// --------------------------------------------------------------------------
// On input:
// count 1, count2 are the number of allele (of interest) exist in the
// removed sample (maximum = 2)
// homo1, homo2 are frequencies of homozygotes at the two loci of interest
// before the sample is removed
// frac = 1/(nSamp-1) to be calculated before calling this function, since
// this will be called repeatedly (to save a repeated arithmetics)
// On output, for the values at two loci of interest:
// *var = fx(1-2fx) + h, where h is the frequency of homozygote.
// Return value = 0: Burrows correlation r^2 is not determined by fraction
// whose dnominator is the product (var1)*(var2)

char r2Default (float nSamp, float frac, float f1x, float homo1,
				int count1, float f2x, float homo2, int count2,
				float *var1, float *var2, float epsilon)
{
	float z, h1x, h2x;
	char rSet = 0;
	z = homo1*(nSamp);
	if (count1 == 2) z -= 1.0;
	h1x = frac*z;	// freq. of homozygote after the sample is removed
	z = homo2*(nSamp);
	if (count2 == 2) z -= 1.0;
	h2x = frac*z;
	*var1 = f1x - 2*f1x*f1x + h1x;
	*var2 = f2x - 2*f2x*f2x + h2x;
	rSet = (*var1 < epsilon || *var2 < epsilon)? 0: 1;
	return rSet;
}

// ---------------------------------------------------------------------------
// Modified Dec 2016:
// Add a local variable cutoffRev to be parameter of IndAlle2, and reset
// cutoff = cutoffRev. When Burrows_Calcul is called again for another
// locus pair, cutoff as parameter of Burrows_Calcul is still the old one.
void Burrows_Calcul (float cutoff, ALLEPTR allep1, ALLEPTR allep2,
                FISHPTR popLoc1, FISHPTR popLoc2, int p1, int p2,
                int nMp1, int nMp2, int nfish, float *nSamp,
                int *nInd1, int *nInd2, int *nMpairs, float *rB,
                int currPop, int *missptr, FILE *outBurr,
                char *outBurrName, char moreBurr, char BurrPause,
                float *expR2, char weighsmp,
// add in Nov 2014/ Jan 2015/ Mar 2015:
                char sepBurOut, char moreCol, char BurAlePair,
// add Mar 2016:
// unidented variables are for checking, remove them together with those in
// LDRunPairs, Pair_Analysis, LDMethod, etc., later:
// for checking with checkR2
//char *opened,
                char jack, int **p1Gen, int **p2Gen, int *noDatFish,
                char *countm1, char *countm2,
                int *mValp1, float *freqp1, float *homop1,
                int *mValp2, float *freqp2, float *homop2,
                float *r2AtPairX,
// temporarily add for checking with checkR2:
//double *r2JackTot,
                float *JweighPair,
                unsigned long long *r2Count, float epsilon)

// calculate Burrows coefficients for a pair of loci (p1, p2).
// nMp1, nMp2: number of alleles at loci p1, p2
// added in Sept 2011 parameter BurrPause
// weighsmp to denote if the population has missing data
// epsilon is a small number, used to check some rational number should be 0

// r2AtPairX is array of r^2 calculated at this locus pair, for each sample
// set where one sample is removed from S. Name such sample set as S# or Sk
// (when sample k is removed); there are "nfish" such sets.
// So r2AtPairX[k] = average of r^2 over all allele pairs at locus pair
// (p1, p2) for sample set Sk.
//
// r2JackTot contains the sum of r^2 up to this locus pair (p1, p2), for
// each sample set S#.
//
// JweighPair contains product of independent alleles at each S#;
// this is used for calculating weighted average of r^2 in each S#.
// r2Count is the number of r^2 evaluated for each sample set S#. Some sample
// set S# may not have r^2 calculated (ineligible for Burrows calculations)
{

	int i, j, k, nEff1, nEff2;
	int countM, m1, m2;
	char writeBur;
	float f1, f2, x, xy;
	float varp1;
	float dBur, rBur, rBur2;
	float rMean, dBurMean, r2Mean;   // for averages of rBur and dBur, r2Mean
	// *rB = r2Mean is the output for average of rBur2 (square of rBur)
	// Added in Nov 2014 for min freq at loci p1, p2:
	float fminp1, fminp2;
// add cutoffRev in Dec 2016, for reassigning cutoff value:
	float cutoffRev;

	IndAlle2 (p1Gen, p2Gen, noDatFish, popLoc1, popLoc2, cutoff, nfish,
			allep1, allep2, nMp1, nMp2, &nEff1, mValp1, &nEff2, mValp2,
			nSamp, freqp1, homop1, freqp2, homop2, nInd1, nInd2,
			weighsmp, &fminp1, &fminp2, &cutoffRev);
// Dec 2016: the rest use cutoff, so we don't want to go to change them,
// just set cutoff to be this reassignment when cutoff is a special value
// for dropping singletons (old cutoff is still unchanged when this exits):
	cutoff = cutoffRev;

// nInd1, nInd2 are the number of alleles having freq > cutoff
// when measured against samples having data at both loci p1, p2.
// However, if there are no freq <= cutoff, nIndJ will be the number
// of alleles at locus J, less 1 (since the sum of freq = 1, only nIndJ
// values of freq are needed in this case, the one left is then known).
// They are defined as numbers of independent alleles at loci p1, p2
// taken against samples having data at both lci.
// nEff1 and nEff2 are the number of alleles in loci p1, p2 whose
// freqs are between cutoff and < 1. Allele having frequency (1-cutoff) and
// above is not counted. Those alleles (counted nEff1, nEff2) are used in
// calculating Burrows coeffs. We have nInd1 <= nEff1, nInd2 <= nEff2.

// nInd1 < nEff1: no allele being dropped (no allele having freq < cutoff)
// Similarly for nInd2 < nEff2.

// Added in Jan 2016:
	float rowSum;
	char rSkip1, rSkip2, rSkip, dSkip;
// Add in March 2016:
	float pSum, pSumx, rBur2x, f1x, f2x, var1, var2;
	char reject;
	char rSet;
	float frac, frac2, dBurx;

/* **********************************************************************
   The following can be proved mathematically.
   (dBur, rBur are Burrows disequilibrium and correlation, rBur2 = (rBur)^2.)
   (1) nInd1 < nEff1 and nInd2 < nEff2 (no allele dropped).
        Then dBurMean = 0.
    a. If nEff1 = 2 (=> nInd1 = 1), there are only 2 alleles A1, A2
       at locus p1, then the Burrows coefficients dBur, rBur for
       the two pairs of loci (A1,B), (A2,B) are equal in absolute value
       and of opposite signs for any allele B at the other locus.
       Similarly for nEff2 = 2.
    b. If both nEff1 = 2 and nEff2 = 2 (each locus has 2 alleles),
       then average of rBur2 (square of rBur), which is one of returned
       value *rB, is any rBur2 at one of the 4 pairs.
   (2) nInd1 < nEff1 and nInd2 = nEff2, no allele dropped at locus
    p1 but some are dropped at locus p2.
        Then dBurMean = 0 (averaged over all eligible locus pairs).
        If also nEff1 = 2, then we have a similar result as in (1)a.
        It follows that rMean = 0.
        Similar results are for nInd1 = nEff1 and nInd2 < nEff2.

*/
// **********************************************************************

// nMpairs is the number of allele pairs at locus pair (p1, p2) at which
// Burrows coefficients are calculated; nInd1*nInd2 is number of ind. alleles
	*nMpairs = nEff1 * nEff2;

	*rB = 0;
	*expR2 = ExpR2Samp(*nSamp);
// added in Sept 2011 to prevent division by zero:
	if (*nMpairs <= 0) return;

// write to Burrows file, add BurrPause in Sept 2011, BurAlePair in Apr 2015:
	writeBur = (outBurr != NULL && moreBurr == 1 && BurrPause == 0)? 1: 0;
	if (writeBur == 1  && BurAlePair == 1)
	{
		if (sepBurOut == 0) {
			fprintf (outBurr, "\n      Pop.    Loc._Pairs   Allele_Pairs    P1"
			"    P2    Burrows->D       r         r^2\n");
		fprintf (outBurr, "   ");
        PrtLines (outBurr, 85, '-');
//  	} else {
//          fprintf (outBurr, "\nLoc._Pairs   Allele_Pairs    P1"
//          "    P2    Burrows->D       r         r^2\n");
		}
		fflush (outBurr);
	}

	rMean= 0.0; dBurMean= 0.0, r2Mean = 0.0;
	// rSkip# = 1 if exactly 2 alleles at the locus, and no allele dropped.
	// Either rSkip1 = 1 or rSkip2 = 1 will imply that the sum of
	// coefficients rBur across all pairs of loci is zero
	rSkip1 = ((*nInd1<nEff1) && (nEff1==2))? 1: 0;
	rSkip2 = ((*nInd2<nEff2) && (nEff2==2))? 1: 0;
	rSkip = rSkip1 + rSkip2;
	// rSkip > 0 implies that the sum of rBur will be zero.
	dSkip = ((*nInd1<nEff1) || (*nInd2<nEff2))? 1: 0;
	// dSkip > 0 implies that the sum of dBur will be zero;
	// rSkip > 0 implies that dSkip > 0.
	// Thus, if rSkip > 0, sum of dBur and sum of rBur will be zero,
	// no need to keep track their sums for their averages dBurMean, rMean

	frac = 1.0/(*nSamp - 1); // *nSamp > 1 is known
	frac2 = frac/2.0;

	dBurMean = 0;
	rMean = 0;
	if (rSkip == 2) {    // each locus has 2 alleles, none is dropped
		m1 = *mValp1;
		f1 = *freqp1;
		varp1 =  f1 * (float) (1.0 - 2*f1) + (*homop1);
		m2 = *mValp2;
		f2 = *freqp2;
		float t =  f2 * (float) (1.0 - 2*f2) + (*homop2);
		if ((varp1 < epsilon) || (t < epsilon)) { // => varp1 = 0 or t = 0
			dBur = 0;   // varp1 = 0 or t = 0: heterozygote throughout with
			rBur = 0;   // m1 at loc p1, or with m2 at loc p2
			rBur2 = 0;
// for jackknife on samples:
// sj ----- sj ----- sj ----- sj ----- sj ----- sj ----- sj ----- sj -----
        // for any sample removed, same thing applies to the remaining:
        // r^2 = 0 for each sample set with one removed,
        // only the number of r^2 evaluated is increased
			if (jack != 0) {
				for (k = 0; k< nfish; k++) {
					r2Count[k]++;
					r2AtPairX[k] = 0;
					JweighPair[k] = 1;
				}
			}
// ej ----- ej ----- ej ----- ej ----- ej ----- ej ----- ej ----- ej -----
		} else {    // in this case, we have *nSamp >= 2 (otherwise, only one
    // sample having data at 2 loci, implying it is heterozygote (varp1=t=0)
			AlleInSamp (nfish, m1, p1Gen, noDatFish, countm1);
			AlleInSamp (nfish, m2, p2Gen, noDatFish, countm2);
			Burrows_Delta (f1, f2, varp1, t, *nSamp, nfish,
				 &dBur, &rBur, &rBur2, &pSum, countm1, countm2);
	// The rest of this "else" are for jackknife on Samples
// sj ----- sj ----- sj ----- sj ----- sj ----- sj ----- sj ----- sj -----
			if (jack != 0) {
				for (k = 0; k < nfish; k++) {
        // for each sample (k+1) removed, let's call the rest as Sk
        // This loop calculates r^2 at each sample set Sk
					if (noDatFish[k] > 0) {		// sample (k+1)th has no data
						r2AtPairX[k] = rBur2;	// same r^2 when removing it
// temporarily add for checking with checkR2:
//r2JackTot[k] += rBur2;
						r2Count[k]++;
						JweighPair[k] = 1;
					} else {
						reject = Rejected(cutoff, *nSamp,f1,countm1[k], &f1x)
							+ Rejected (cutoff, *nSamp,f2, countm2[k], &f2x);
				// if reject = 0: both alleles are accepted
						if (reject == 0) {
							rSet = r2Default (*nSamp, frac, f1x, *homop1,
									countm1[k], f2x, *homop2, countm2[k],
									&var1, &var2, epsilon);
				// rSet = 0 when var1 or var2 = 0, then r^2 = 0 on Sk
							if (rSet != 0) {
								dBurx = pSum - countm1[k]*countm2[k];
								dBurx *= frac2;
								dBurx -= (2*f1x*f2x);   // Burrows for Sk
								if (*nSamp > 2.5)   // i.e., is at least 3
						// adjust Burrows disequilibrium by a factor
								dBurx = dBurx *((*nSamp-1)/(*nSamp-2));
								rBur2x = (dBurx*dBurx)/(var1*var2);
								if (rBur2x > 1.0) rBur2x = 1.0; // pull back!
								r2AtPairX[k] = rBur2x;
// temporarily add for checking with checkR2:
//r2JackTot [k] += rBur2x;    // r^2 for Sk
							}
							r2Count[k]++;
							JweighPair[k] = 1;
						} else {	// this pair of loci is rejected in Sk
							JweighPair[k] = 0;
							r2AtPairX[k] = 0;
						} // end of "if (reject == 0) ... else"
				// if reject != 0, no r^2 on Sk, r2Count[k] not increased
					}   // end of "if (noDatFish[k] > 0) ... else ..."
				}   // end of "for (k = 0; k< nfish; k++)"
			}   // end of "if (jack != 0)"
// ej ----- ej ----- ej ----- ej ----- ej ----- ej ----- ej ----- ej -----
		}   // end of "if ((varp1 < epsilon) || (t < epsilon)) ... else ..."
		*rB = rBur2;
// write to auxiliary file here
		if (writeBur == 1) {
// This next line is for testing purpose only, comment out later
//fprintf (outBurr, "Short cut\n");
			if (BurAlePair == 1)
			{
				for(i = 0; i < 2; i++) {
					if (i == 1) {   // Note: at i = 1, the loop j just ends
				// with i = 0, j = 1. The dBur, rBur change signs by the
				// code in next loop, which are then the values at allele
				// pair corresponding to i = 1,j = 0 (starting next loop)
						m1 = *(mValp1+1); m2 = *mValp2;
						f1 = *(freqp1+1); f2 = *freqp2;
					}
					for (j = 0; j < 2; j++) {
						if (j == 1) {
							m2 = *(mValp2+1); f2 = *(freqp2+1);
							dBur = -dBur; rBur = -rBur;
						}
						if (sepBurOut == 0) fprintf (outBurr,
							"%9d%8d%6d%8d%6d  %7.3f%7.3f%11.6f%12.6f%12.6f\n",
							currPop, p1+1, p2+1, m1, m2, f1, f2, dBur,
							rBur, rBur2);
						else fprintf (outBurr,
							"%3d%6d%8d%6d  %7.3f%7.3f%11.6f%12.6f%12.6f\n",
							p1+1, p2+1, m1, m2, f1, f2, dBur, rBur, rBur2);
					}
				}
				if (sepBurOut == 0) {
					fprintf (outBurr, "   ");
        			PrtLines (outBurr, 85, '-');
					fprintf (outBurr, "   Number of Allele Pairs:%8d,      "
							"Means:  %15.3e%12.3e%12.3e\n",
							*nMpairs, dBurMean, rMean, *rB);
// added, but can be removed (Aug 2016) ----------------------------------
					float w1 = 1.0;
					float w2 = (*nSamp)*(*nSamp);
					if (weighsmp > 0)
						fprintf (outBurr, "%49cIndp. = (1, 1), Size =%5.0f,"
									" Wt:%7.0f\n", ' ', *nSamp, w1*w2);
					else
						fprintf (outBurr, "%78cWt:%7.0f\n", ' ', w1);
// -----------------------------------------------------------------------
				}
			} else {
				if (moreCol == 0)
					fprintf (outBurr, "%6d %6d %7.4f %7.4f %8d %14.5e %14.5e\n",
							p1+1, p2+1, fminp1, fminp2, (int) *nSamp, *rB,
							(*rB)-(*expR2));
				else fprintf (outBurr,  "%6d %6d %7.4f %7.4f %5d%5d "
							"%7d%7d %13.4e %13.4e %13.4e %13.4e\n",
							p1+1, p2+1, fminp1, fminp2, *nInd1, *nInd2,
							*nMpairs, (int) *nSamp, dBurMean, rMean,
							*rB, (*rB)-(*expR2));
			}
			fflush (outBurr);
		} // end of "if (writeBur == 1)"
		return;
	}

// Added in Jan 2016:
//	float *y = (float*) malloc(sizeof(float)*nEff2); // changed y to array
	float *varp2 = (float*) malloc(sizeof(float)*nEff2);
	float *colSum = (float*) malloc(sizeof(float)* nEff2);
	float *rRow = (float*) malloc(sizeof(float)*nEff2);
	float *r2Row = (float*) malloc(sizeof(float)*nEff2);

// for jackknife on Samples
// In the comments, we use the symbol S for the whole set of samples,
// Sk for the sample set S with one sample removed (assuming sample (k+1)th)
// sj ----- sj ----- sj ----- sj ----- sj ----- sj ----- sj ----- sj -----
	int c;
	float *f1xAt = (float*) malloc(sizeof(float)*nfish);
	float *f2xAt = (float*) malloc(sizeof(float)*nfish);
// to keep track of r^2 at all Sk, for a pair of alleles
	float *r2xAt = (float*) malloc(sizeof(float)*nfish);
// to count number of eligible alleles at each Sk
	int *m1Acc = (int*) malloc(sizeof(int)*nfish);
	int *m2Acc = (int*) malloc(sizeof(int)*nfish);
// to determine if current allele pairs are rejected at Sk
	char *m1Rej = (char*) malloc(sizeof(char)*nfish);
	char *m2Rej = (char*) malloc(sizeof(char)*nfish);
	char gotr2x;
// for the sum of frequencies of eligible alleles at Sk
	float *f1xSum = (float*) malloc(sizeof(float)*nfish);
	float *f2xSum = (float*) malloc(sizeof(float)*nfish);

	for (k = 0; k < nfish; k++) {
		f1xAt[k] = 0;
		f2xAt[k] = 0;
		r2xAt[k] = 0;
		m1Acc[k] = 0;
		m2Acc[k] = 0;
		m1Rej[k] = 0;
		m2Rej[k] = 0;
		f1xSum[k] = 0;
		f2xSum[k] = 0;
		r2AtPairX[k] = 0;
	}
// ej ----- ej ----- ej ----- ej ----- ej ----- ej ----- ej ----- ej -----
	for (j = 0; j < nEff2; j++) {
		f2 = *(freqp2+j);   // frequency of allele (j+1)th at locus p2.
		varp2[j] = f2 * (float) (1.0 - 2*f2) + *(homop2+j); // for variance
                        // related to allele (j+1)th at locus p2
		colSum[j] = 0;              // related to allele (j+1)th at locus p2
		rRow[j] = 0;    // colSum, rRow, r2Row for holding values calculated
		r2Row[j] = 0;
// for jackknife on Samples:
// sj ----- sj ----- sj ----- sj ----- sj ----- sj ----- sj ----- sj -----
// Determine m2Acc[k] = number of alleles at locus p2 that are accepted at Sk
		if (jack != 0) {
			if (varp2[j] < epsilon) {   // locus p2 is het. with allele (j+1)th,
	// so allele (j+1)th is eligible at locus p2 (freq = 0.5) in all Sk.
				for (k = 0; k < nfish; k++) m2Acc[k]++;
			} else {
				m2 = *(mValp2+j);   // m2 is allele (j+1)th
				for (k = 0; k < nfish; k++) {
					if (noDatFish[k] > 0) m2Acc[k]++;
					else {
						c = Count(p2Gen[k], m2);
						if (Rejected(cutoff, *nSamp, f2, c, (f2xAt+k)) == 0) {
							m2Acc[k]++;
							f2xSum[k] += f2xAt[k];
						}
					}
				}
			}
		}
// ej ----- ej ----- ej ----- ej ----- ej ----- ej ----- ej ----- ej -----
	}

// In the next loop, for each i, rRow[.] is a vector of rBur-values corresp.
// to i, at various j (in loop of j). It contains rBur(i,j) for various j.
// The values in this array are only needed for storing at i = 0 in the case
// 1 = nInd1 < nEff1 = 2.  When going to the next i = 1 in the loop, those
// values in rRow (stored rBur values at i = 0 for various j) are reused.
// In such case, then for each j, rBur(1,j) = -rBur(0,j). Thus, at each j,
// just set rBur(1,j) = -rRow[j].
// Simlarly for r2Row (storing values of rBur2).

// colSum[j] stores the sum of dBur-values at column j, up to row i = nInd1-1.
// This will be used at row i = nInd1, which can only happen if nInd1 < nEff1.
// When this happens (at i = nInd1), the dBur at this i, for any column j,
// is the negative of the sum of all previous dBur values under the same j.
// That sum was stored as colSum[j].
//
// On the other hand, rowSum is the sum of dBur at row i, up to j = (nInd2-1).
// This will be used at each i, and only at j = nInd2 (this can only happen
// if nInd2 < nEff2). The value of dBur at (i, nInd2) is negative of the sum
// of all dBur at the same i, of all j < nInd2. That sum is rowSum.
	for (i = 0; i < nEff1; i++) {
		m1 = *(mValp1+i);
		f1 = *(freqp1+i);       // frequency of allele m1 at locus p1.
		varp1 =  f1 * (float) (1.0 - 2*f1) + *(homop1+i);   // represent variance
                                                // related to allele m1
		if (varp1 < epsilon) {  // this means variance related to m1 is zero;
		// so, all samples are heterozygotes with one allele being m1.
			dBur = 0; rBur = 0; rBur2 = 0;
// write to auxiliary file here
			if (writeBur == 1 && BurAlePair == 1)
			{
				for (j = 0; j < nEff2; j++) {
					m2 = *(mValp2+j);
					f2 = *(freqp2+j);   // frequency of allele m2 at locus p2.
					if (sepBurOut == 0) fprintf (outBurr,
						"%9d%8d%6d%8d%6d  %7.3f%7.3f%11.6f%12.6f%12.6f\n",
						currPop, p1+1, p2+1, m1, m2, f1, f2, dBur, rBur, rBur2);
					else fprintf (outBurr,
						"%3d%6d%8d%6d  %7.3f%7.3f%11.6f%12.6f%12.6f\n",
						p1+1, p2+1, m1, m2, f1, f2, dBur, rBur, rBur2);
					fflush (outBurr);
				}
			}
        // The next "for" loop is for jackknife on Samples.
// sj ----- sj ----- sj ----- sj ----- sj ----- sj ----- sj ----- sj -----
        // The same thing will happen when a sample is removed (whether
        // that sample has data at both loci or not): all samples in Sk
        // are heterozygotes with one allele being m1 (allele m1 in Sk
        // has freq 1/2); so r^2 = 0 for each Sk when m1 is paired
        // with any eligible allele at locus p2).

        // r2AtPairX[k] = (sum of r^2 over all allele pairs in Sk) unchanged
			if (jack != 0) {
				for (k = 0; k < nfish; k++) {
//          r2xAt[k] = 0;   // r^2 set = 0 for any pairing with m1
				// this allele m1 is eligible in all Sk:
					m1Rej[k] = 0;
					m1Acc[k]++;
				}
			}
// ej ----- ej ----- ej ----- ej ----- ej ----- ej ----- ej ----- ej -----


		} else {    // of "(varp1 < epsilon)" => varp1 (related to m1) > 0.
		// obtain countm1[k] = number of copies of allele m1 at each sample k
			AlleInSamp (nfish, m1, p1Gen, noDatFish, countm1);
// for jackknite on Sample
// sj ----- sj ----- sj ----- sj ----- sj ----- sj ----- sj ----- sj -----
// must check if this allele m1 is existed and accepted in each Sk
			if (jack != 0) {
				for (k = 0; k < nfish; k++) {
					if (noDatFish[k] > 0) {
						m1Rej[k] = 0;
						m1Acc[k]++;
					} else {
						m1Rej[k] = Rejected(cutoff, *nSamp, f1, countm1[k],
										(f1xAt+k));
						if (m1Rej[k] == 0) {
							m1Acc[k]++;
							f1xSum[k] += f1xAt[k];
						}
					}
				}
			}
// As i goes through the "for" loop, m1Acc[k] is the total number of eligible
// alleles in the sample set Sk (sample set after a sample k is deleted).
// The number of alleles accepted in Sk was found before this "for" loop.
// Therefore, the number of pairs accepted in Sk is m1Acc[k] * m2Acc[k]
// after this "for i = " loop is through.
// ej ----- ej ----- ej ----- ej ----- ej ----- ej ----- ej ----- ej -----
			rowSum = 0;
			for (j = 0; j < nEff2; j++) {
				m2 = *(mValp2+j);
				f2 = *(freqp2+j);   // frequency of allele m2 at locus p2.
				gotr2x = 0;
				if (varp2[j] < epsilon) {   // varp2[j] = 0, variance at m2
			// all samples are heterozygotes at locus p2 with one allele = m2
					dBur = 0; rBur = 0; rBur2 = 0;
// sj ----- sj ----- sj ----- sj ----- sj ----- sj ----- sj ----- sj -----
				// for jackknife on samples:
				// For all k, r^2 at Sk = 0,  r2AtPairX[k] is unchanged
				// This allele m2 is accepted in all Sk
					if (jack != 0) {
						for (k = 0; k < nfish; k++) m2Rej[k] = 0;
					}
// ej ----- ej ----- ej ----- ej ----- ej ----- ej ----- ej ----- ej -----

				} else {    // both varp1 and varp2[j] > 0 (rel. to m1, m2)
		// obtain countm2[k] = number of copies of allele m2 at each sample k
					AlleInSamp (nfish, m2, p2Gen, noDatFish, countm2);
// sj ----- sj ----- sj ----- sj ----- sj ----- sj ----- sj ----- sj -----
				// jackknife on samples: check eligibility of m2 in each Sk
					if (jack != 0) {
						for (k = 0; k < nfish; k++) {
							if (noDatFish[k] > 0) m2Rej[k] = 0;
							else m2Rej[k] = Rejected(cutoff, *nSamp, f2,
											countm2[k], (f2xAt+k));
						}
					}
// ej ----- ej ----- ej ----- ej ----- ej ----- ej ----- ej ----- ej -----

// Always have nInd1 <= nEff1. Then nInd1 < nEff1 <=> no allele dropped at
// locus p1. The next "if" statement is to deal with this case.
// First, dBur for the last i (which is nInd1), at any j, will be the negative
// of the sum of all dBur at previous i (of the same j), which is colSum[j].
// Moreover, if also nEff1 = 2 (or nInd1 = 1), then rBur at i = nInd1 = 1 and
// at any j will be the negative of rBur at i = 0 and the same j (= rRow[j])
					if (i == *nInd1) {  // only if *nInd1<nEff1, i=nEff1-1
						dBur = - colSum[j];
						if (nEff1 == 2) {
							rBur = - rRow[j];
							rBur2 = r2Row[j];
						} else {
							xy = varp1 * varp2[j];  // this must be > 0
							rBur = dBur/sqrt(xy);
							rBur2 = rBur * rBur;
							if (rBur2 > 1.0) rBur2 = 1.0;
						}
// for Jackknife on samples:
// sj ----- sj ----- sj ----- sj ----- sj ----- sj ----- sj ----- sj -----
// At allele m1 of locus p1, we know if m1 is rejected or eligible under Sk
// for each sample (k+1)th removed, from the value m1Rej[k] evaluated earlier
// at the start of "else" of "if (varp1 < epsilon)". When m2 is not rejected
// under Sk, r^2 will be calculated.
// When allele pair (m1, m2) is eligible under Sk, and i = *nInd1 = 1, then
// the previous m1 (corresponding to i=0) is also eligible at locus p1,
// we actually have r^2 at this pair (m1, m2) being the same as r^2 at
// (prev m1, m2). However, to be able to use this previous r^2-value, we will
// need to store it in a 2-dimensional array (one for m2, the other for k,
// to store r^2 with the first allele being fixed).
// Without such array, we will need to calculate r^2 for each Sk
// (when "nEff1 == 2"), based on the Burrows disequilibrium dBur for the
// whole sample set S, instead of grabbing previous values in each Sk.
// So, we will need to use the following value as in the case "else" above:
						if (jack != 0) {
							if (j == *nInd2 && j == 1) {
								gotr2x = 1;
								for (k = 0; k < nfish; k++) {
// Note: In this "j=1" case, when allele m1, m2 are eligible in Sk, then
// allele m1 ane previous m2 (corresponding to j = 0) will also eligible
// in Sk, so r2xAt[k] is the available r2-value for previous pair
									if (m1Rej[k] == 0 && m2Rej[k] == 0) {
										r2AtPairX[k] += r2xAt[k];
									}
								}
							} else //  of "(j == *nInd2 && j == 1)"
								pSum = roundf((dBur + 2*f1*f2)*2*(*nSamp));
						}	// end of "if (jack != 0)"
// ej ----- ej ----- ej ----- ej ----- ej ----- ej ----- ej ----- ej -----
					} else {    // of "if (i == *nInd1)", i.e., i < *nInd1.
					// If *nInd1 = nEff1, then case "else" always happens.
						if (j == *nInd2) {  // can happen if nInd2 < nEff2
							dBur = - rowSum;
							if (j == 1) {
								rBur = -rBur; // prev. rBur is for (i, 0)
// for Jackknife on samples:
// sj ----- sj ----- sj ----- sj ----- sj ----- sj ----- sj ----- sj -----
// At each k, either the allele pair corresponding to (i, j=1) is rejected
// or r^2 at this allele pair is the same as for previous pair (i, j=0).
// If pair at (i, j=1) is eligible, then so is pair at (i, j=0).
								if (jack != 0) {
									gotr2x = 1;
									for (k = 0; k < nfish; k++) {
// Note: In this "j=1" case, when allele m1, m2 are eligible in Sk, then
// allele m1 ane previous m2 (corresponding to j = 0) will also eligible
// in Sk, so r2xAt[k] is the available r2-value for previous pair
										if (m1Rej[k] == 0 && m2Rej[k] == 0) {
											r2AtPairX[k] += r2xAt[k];
										}
									}
								}
// ej ----- ej ----- ej ----- ej ----- ej ----- ej ----- ej ----- ej -----
							} else { // j = nInd2 > 1: need to find rBur
								xy = varp1 * varp2[j];
								rBur = dBur/sqrt(xy);
								rBur2 = rBur * rBur;
								if (rBur2 > 1.0) rBur2 = 1.0;
// for Jackknife on samples:
// sj ----- sj ----- sj ----- sj ----- sj ----- sj ----- sj ----- sj -----
								pSum = roundf((dBur + 2*f1*f2)*2*(*nSamp));
// ej ----- ej ----- ej ----- ej ----- ej ----- ej ----- ej ----- ej -----
							}
						} else {	// of "if (j == *nInd2)"
        // Burrows coeff. need to be calculated at pair (i, j), for any j.
        // (countm1, countm2 were determined)
							Burrows_Delta (f1, f2, varp1, varp2[j], *nSamp, nfish,
								&dBur, &rBur, &rBur2, &pSum, countm1, countm2);
							rowSum += dBur; // sum of dBur across j, for this i
						}
						colSum[j] += dBur;  // colSum[j] is sum of entries of
                                // dBur-matrix at column j, up to row i
						rRow[j] = rBur; // rRow, r2Row are reassigned at each i;
						r2Row[j] = rBur2;// only values for the last i are kept!
					}   // end of "if (i == *nInd1) ... else ..."
// for Jackknife on samples:
// sj ----- sj ----- sj ----- sj ----- sj ----- sj ----- sj ----- sj -----
// After the condition set for i (= nInd1 or not), except for the case
// j = *nInd2 = 1 where r^2 in each Sk takes previous r^2-value (at j=0),
// the rest (unconditional on i, and for j not in the case mentioned) need
// value pSum. This pSum was assigned in case i = *nInd1. In case i /= *nInd1,
// it was set at j = *nInd2 /= 1 through assignment, and at j /= *nInd2
// through the call Burrows_Delta.
// The exceptional case, where r^2 in each Sk takes previous r^2-value,
// is signaled by gotr2x = 1 (otherwise, it is still 0).
					if (jack != 0) {
						if (gotr2x == 0) {
// From Burrows disequilibrium at (m1, m2) under whole sample set S, which is
// dBur = pSum - 2(f1*f2), we can derive Burrows disequilibrium at (m1,m2)
// under sample set Sk. The calculations are based on "pSum". This pSum is
// obtained through the call Burrows_Delta where dBur was calculated, or
// straight from dBur if dBur was deduced from previous values.
							for (k = 0; k < nfish; k++) {
								if (m1Rej[k] == 0 && m2Rej[k] == 0) {
                    // rBur2x will be r^2 for Sk
									if (noDatFish[k] > 0) { // missing data at k
										r2xAt[k] = rBur2;
										r2AtPairX [k] += r2xAt[k];
									} else {
										rSet = r2Default (*nSamp, frac, f1xAt[k],
											*(homop1+i), countm1[k], f2xAt[k],
											*(homop2+j), countm2[k],
											&var1, &var2, epsilon);
										if (rSet != 0) {
											pSumx = pSum - countm1[k]*countm2[k];
											dBurx = frac2*pSumx -
													2*f1xAt[k]*f2xAt[k];
											if (*nSamp > 2.5) // i.e. at least 3
							// adjust Burrows disequilib. by unbiased factor
												dBurx *= (*nSamp-1)/(*nSamp-2);
											r2xAt[k] = (dBurx*dBurx)/(var1*var2);
											if (r2xAt[k] > 1.0) r2xAt[k] = 1.0;
											r2AtPairX [k] += r2xAt[k];   // r^2 for Sk
										} else {    // if rset=0, r^2 = 0
											r2xAt[k] = 0;
										}
									}

								}   // end of "if (m1Rej[k]==0 && m2Rej[k]==0)"
							}   // end of "for (k = 0; k < nfish; k++)"
						}   // end of "if (gotr2x == 0)"

					}   // end of "if (jack != 0)"
// ej ----- ej ----- ej ----- ej ----- ej ----- ej ----- ej ----- ej -----
				}   // end of "if (varp2[j] < epsilon) ...  else ..."
				r2Mean += rBur2;
				if (rSkip == 0) rMean += rBur;  // rMean = 0 if rSkip > 0
				if (dSkip == 0) dBurMean += dBur;
// write to auxiliary file here
				if (writeBur == 1 && BurAlePair == 1)
				{
					if (sepBurOut == 0) fprintf (outBurr,
							"%9d%8d%6d%8d%6d  %7.3f%7.3f%11.6f%12.6f%12.6f\n",
						currPop, p1+1, p2+1, m1, m2, f1, f2, dBur, rBur, rBur2);
					else fprintf (outBurr,
							"%3d%6d%8d%6d  %7.3f%7.3f%11.6f%12.6f%12.6f\n",
							p1+1, p2+1, m1, m2, f1, f2, dBur, rBur, rBur2);
					fflush (outBurr);
				}
			}   // end of "for (j = 0; j < nEff2; j++)"
		// This is the end of "else": varp1 is not 0; so r2, dBur can be
		// nonzero when j runs in the "for" loop
		}   // end of "if (varp1 < epsilon) ... else ..."
	}    // end of "for (i = 0; i < nEff1; i++)"

// for Jackknife on Samples:
// sj ----- sj ----- sj ----- sj ----- sj ----- sj ----- sj ----- sj -----
	if (jack != 0) { // wrap up information on Sk
		x = 1.0 - epsilon;  // basically, if y <= 1 and y > x, then y = 1.
		for (k = 0; k < nfish; k++) {
			i = m1Acc[k] * m2Acc[k];
			if (i > 0) {
				r2AtPairX[k] /= i;   // mean of r^2 over all elig. alle pairs
// temporarily add for checking with checkR2:
// r2JackTot[k] += r2AtPairX[k];
				r2Count[k]++;
			}
		// now for product of ind. alleles at each Sk:
			if (noDatFish[k] > 0) JweighPair[k] = (*nInd1)*(*nInd2);
			else {
				m1 = m1Acc[k];  // number of eligible alleles at locus p1
				m2 = m2Acc[k];  // number of eligible alleles at locus p2
			// if Sum(eligible alleles) = 1, #(ind. alleles) is one less
				if (f1xSum[k] > x) m1--;
				if (f2xSum[k] > x) m2--;
				JweighPair[k] = m1*m2;
			}
		}
	}
// Now, r2AtPairX[k] is the average of r^2 taken over all allele pairs for
// sample set Sk at locus pair (p1, p2); so it represents r^2 at this pair.
// ej ----- ej ----- ej ----- ej ----- ej ----- ej ----- ej ----- ej -----
	*rB = r2Mean/(*nMpairs); // Burrows coeff. averaging over all allele
                            // pairs for this locus pair (p1, p2).
	rMean = rMean/(*nMpairs);
	dBurMean = dBurMean/(*nMpairs);
	// write to auxiliary file here (added BurrPause in Sept 2011):
	if (writeBur == 1)
	{
		if (BurAlePair == 1) {
			if (sepBurOut == 0) {
				fprintf (outBurr, "   ");
				for (j=0; j<85; j++) fprintf (outBurr, "-");
				fprintf (outBurr, "\n");
				fprintf (outBurr, "   Number of Allele Pairs:%8d,"
						"      Means:  %15.3e%12.3e%12.3e\n",
						*nMpairs, dBurMean, rMean, *rB);
// added, but can be removed (Aug 2016) ----------------------------------
				float w1 = (*nInd1)*(*nInd2);
				float w2 = (*nSamp)*(*nSamp);
				if (weighsmp > 0)
					fprintf (outBurr, "%47cIndp. = (%2d, %2d), Size =%5.0f,"
							" Wt:%7.0f\n", ' ', *nInd1, *nInd2, *nSamp, w1*w2);
				else
					fprintf (outBurr, "%78cWt:%7.0f\n", ' ', w1);
// -----------------------------------------------------------------------
			}
		} else {
			if (moreCol == 0)
				fprintf (outBurr, "%6d%7d%8.4f%8.4f %8d %14.5e %14.5e\n",
				p1+1, p2+1, fminp1, fminp2, (int) *nSamp, *rB, (*rB)-(*expR2));
			else fprintf (outBurr,
				"%6d%7d%8.4f%8.4f%6d%5d%8d%7d%14.4e%14.4e%14.4e%14.4e\n",
				p1+1, p2+1, fminp1, fminp2, *nInd1, *nInd2, *nMpairs,
				(int) *nSamp, dBurMean, rMean, *rB, (*rB)-(*expR2));
		}
		fflush (outBurr);
	};

	free (colSum);
	free (varp2);
	free (rRow);
	free (r2Row);
// for jackknife:
	free (f1xAt);
	free (f2xAt);
	free (r2xAt);
	free (m1Acc);
	free (m2Acc);
	free (m1Rej);
	free (m2Rej);
	free (f1xSum);
	free (f2xSum);
}
//-------------------------------------------------------------------------

void AddBurrVal (int nInd1, int nInd2, float rB, float nSamp,		// in
					float expR2, char weighsmp, int locSkip,		// in
					unsigned long long nLocPairs,					// in
// the next 5 are arrays
					float *rB2, float *rBdrift, float *prodInd,		// in-out
					float *sampCount, float *pairWt,				// in-out
					FILE *rAveTemp,									// in-out
// because of rouding off errors, use double to calculate (since Apr 2012):
					double *totInd,	// in-out, total ind. alle.
					double *wMeanSamp, // in-out, total wt. inv. of samp.size
//					float *rBAveTotal,
					double *rWeight, double *bigExpR2,				// in-out
					double *bigRprime, double *bigR)				// in-out
{
	float rdrift;
	float reNSamp, nIndtot, weight, rBweight;

	nIndtot = (float) nInd1 * nInd2;
	weight = nIndtot;
	// adjust weight by sample size: (added Jul 18 2010?)
	// (In previous LD, no weighsmp condition is needed)
	// weighsmp > 0 when there are missing data; otherwise nSamp as
	// sample having data should be a constant across loci
	if (weighsmp > 0) weight *= (nSamp*nSamp);
	reNSamp = 0.0;
	if (nSamp > 0) reNSamp = (float) 1.0/nSamp;
	*wMeanSamp += nIndtot*reNSamp;
	*totInd += nIndtot;
	*rWeight += weight;
	rBweight = rB * weight;
	*bigR += rBweight;
	*bigExpR2 += (expR2 * weight);
	rdrift = rB - expR2;
	if (rAveTemp != NULL) {
//		*rBAveTotal += rB;
		fwrite (&nIndtot, sizeof(float), 1, rAveTemp);
		fwrite (&nSamp, sizeof(float), 1, rAveTemp);
		fwrite (&rB, sizeof(float), 1, rAveTemp);
		fwrite (&rdrift, sizeof(float), 1, rAveTemp);
		fwrite (&weight, sizeof(float), 1, rAveTemp);
	} else {	// because temporary file cannot be opened
		*(prodInd + (nLocPairs)) = nIndtot;
		*(sampCount + (nLocPairs)) = nSamp;
		*(rB2 + (nLocPairs)) = rB;
		*(rBdrift + (nLocPairs)) = rdrift;
		*(pairWt + (nLocPairs)) = weight;
	}

	rdrift *= weight;		// put in weight
	*bigRprime += rdrift;	// total weighted sum of rBdrift
}
//-------------------------------------------------------------------------
void JackWeight (char weighsmp, float nSamp, int nfish, int *noDatFish,
				float *r2AtPairX, double *r2WRemSmp, float *JweighPair,
				double *JweightTot)
// on input, JweighPair[k] is assumed to be the product of ind. alleles
// at a locus pair, for a collection of sample sets Sk (sample set S minus
// the kth element, for k = 0, ..., nfish-1)
// nSamp is the number of samples in S having data at both loci
// r2AtPairX[k] is r^2-value at the locus pair, for each Sk
// noDatFish[k] > 0 when the kth sample has missing data at locus pair
// JweighPair[k] is holding the product of ind. alleles at the locus pair

// On output:
// JweighPair[k] is multiplied by square of #samples having data when needed;
// however, this will be reset by the calculation of r^2 at next locus pair.
// r2WRemSmp[k] is the sum of (r^2-value)x(weight), up to this locus pair
// JweightTot[k] is the sum of the weights up to this locus pair, for each Sk
// When this function is finished for all locus pairs, then
// r2WRemSmp[k]/JweightTot[k] is the weighted average of r^2 for Sk.
// After this function is called, r2WRemSmp[k] is reassigned as
// r2WRemSmp[k]/JweightTot[k] when needed.
{
	int k;
	if (weighsmp > 0) {
		float w0 = nSamp*nSamp;
		float w1 = (nSamp-1)*(nSamp-1);
		for (k = 0; k < nfish; k++) {
			if (noDatFish[k] > 0) JweighPair[k] *= w0;
			else JweighPair[k] *= w1;
		}
	}
	for (k = 0; k < nfish; k++) {
		r2WRemSmp[k] += (r2AtPairX[k] * JweighPair[k]);
		JweightTot[k] += JweighPair[k];
	}

}

//-------------------------------------------------------------------------
// return the number of r^2-values in parameter list,
// The function returns the number of locus pairs calculated in LD method

unsigned long long LDRunPairs (
					float cutoff, ALLEPTR *alleList, int currPop,	// in
					int nfish, FISHPTR *fishHead, int *nMobil,		// in
					int *missptr, int lastOK, char *okLoc,			// in
					FILE *outBurr,									// in-out
					char moreBurr, char *outBurrName,				// in
// the next 5 are arrays for storing value/locus pair in case rAveTemp = null
					float *rB2, float *rBdrift, float *prodInd,		// in-out
					float *sampCount, float *pairWt,				// in-out
					char weighsmp, int locSkip,						// in
					FILE *rAveTemp,									// in-out
// because of rouding off errors, use double to calculate (Apr 2012):
					double *totInd,	// in-out, total ind. alle.
					double *wMeanSamp, // in-out, total wt. inv. of samp.size
//					float *rBAveTotal,
// currently, tot. weight *rWeight=*totInd, since ind.alle. is used as weight
					double *rWeight, double *bigExpR2,				// in-out
					double *bigRprime, double *bigR,				// in-out
// *nPairPtr = total pairs to be listed in Burrows file
					unsigned long long *nPairPtr,					// in-out
// *npairTot = total loc. pairs, *npairSkip = number of loc. pairs skipped
					unsigned long long *npairTot, long *npairSkip,	// in-out
					unsigned long long prompt,	// in
					// (to inform the user after "prompt" pairs calculated)
					char sepBurOut, char moreCol, char BurAlePair,	// in
// add in Mar 2016
					char jack, int **p1Gen, int **p2Gen, int *noDatFish,
					char *countm1, char *countm2,
					int *mValp1, float *freqp1, float *homop1,
					int *mValp2, float *freqp2, float *homop2, float *r2AtPairX,
// temporarily add for checking:
//double *r2JackTot,
					unsigned long long *r2Count,
					double *r2WRemSmp, float *JweighPair, double *JweightTot,
// temporarily add for checking:
//double *r2Ave,
//char *opened,

					float epsilon)
{
	char BurrPause = 0;
	int p1, p2;
	int nInd1, nInd2;
	float expR2, rB, nSamp;
	int nMpairs;
// comment out the next line sinse those are used in
//	float rdrift, reNSamp, nIndtot, weight, rBweight;

	unsigned long long pairval;

    unsigned long long nLocPairs = 0;	// number of locus pairs, return value

	ALLEPTR allep1, allep2;
	FISHPTR popLoc1, popLoc2;
	// to inform the user after "prompt" pairs calculated;
	pairval = prompt;

// temporarily add for checking
//*r2Ave = 0;

// In the next "for" loop, pick a locus in the ascending order, another locus
// from the set of loci at the order after the first, then go through all
// allele in the mobility lists corresponding to this pair of loci. These pairs
// of loci are in the set of accepted pairs given by array okLoc, which was
// determined by function Loci_Eligible.
// Then calculate Burrows coefficients by function Burrows_Calcul.
	for (p1=0; (p1<lastOK); p1++) {
		if (*(okLoc+p1) == 0) continue;
		allep1 = *(alleList+p1);
		popLoc1 = *(fishHead+p1);
		for (p2=p1+1; (p2<=lastOK); p2++) {
			if (*(okLoc+p2) == 0) continue;	// locus (p2+1) is skipped.
			(*npairTot)++;
			allep2 = *(alleList+p2);
			popLoc2 = *(fishHead+p2);
// changed in Nov 30, 11:
			if (p1 - locSkip >= LOCBURR || p2 - locSkip >= LOCBURR) {
				BurrPause = 1;
			} else {
				BurrPause = 0;
				(*nPairPtr)++;
			}
			Burrows_Calcul (cutoff, allep1, allep2, popLoc1, popLoc2,
					p1, p2, *(nMobil+p1), *(nMobil+p2), nfish, &nSamp,
					&nInd1, &nInd2, &nMpairs, &rB, currPop, missptr,
					outBurr, outBurrName, moreBurr, BurrPause, &expR2,
					weighsmp, sepBurOut, moreCol, BurAlePair,
// temporarily added for checking:
//opened,
					jack, p1Gen, p2Gen, noDatFish, countm1, countm2,
					mValp1, freqp1, homop1, mValp2, freqp2, homop2,
					r2AtPairX,
// temporarily add for checking:
//r2JackTot,
					JweighPair, r2Count, epsilon);
			if (nMpairs <= 0) {
				(*npairSkip)++;
				continue;
			}

// temporarily add to check
//*r2Ave += rB;
			AddBurrVal (nInd1, nInd2, rB, nSamp, expR2, weighsmp, locSkip,
					nLocPairs, rB2, rBdrift, prodInd, sampCount, pairWt,
					rAveTemp, totInd, wMeanSamp, rWeight, bigExpR2,
					bigRprime, bigR);
			if (jack != 0) JackWeight (weighsmp, nSamp, nfish, noDatFish,
							r2AtPairX, r2WRemSmp, JweighPair, JweightTot);
// add this prompt to inform the user the progress:
			if ((nLocPairs) == pairval) {
				printf ("%18llu done, at loc. pair (%d, %d)\n", pairval, p1+1, p2+1);
				pairval +=prompt;
			}
            // count the number of locus pairs nLocPairs.
            // It's important that increment happens after assignments above,
            // to avoid going over the limit of the declared sizes of arrays
            // prodInd, etc.
			nLocPairs++;
		}
	}
	return nLocPairs++;
}
//-------------------------------------------------------------------------
// Version of LDRunPairs, but locus pairs are taken across chromosomes

unsigned long long LDTwoChromo (
				float cutoff, ALLEPTR *alleList, int currPop,	// 3 in
				int nfish, FISHPTR *fishHead, int *nMobil,		// 3 in
				int *missptr, int lastOK, char *okLoc,			// 3 in
				FILE *outBurr,									// 1 in-out
				char moreBurr, char *outBurrName,				// 2 in
// the next 5 are arrays for storing value/locus pair in case rAveTemp = null
				float *rB2, float *rBdrift, float *prodInd,		// 3 in-out
				float *sampCount, float *pairWt,				// 2 in-out
				char weighsmp, int locSkip,						// 2 in
				FILE *rAveTemp,									// 1 in-out
// because of rouding off errors, use double to calculate (Apr 2012):
				double *totInd,	// total ind. alle.				// 1 in-out
				double *wMeanSamp, // total wt. inv. of samp.size  1 in-out,
//				float *rBAveTotal,
// currently, tot. weight *rWeight=*totInd, since ind.alle. is used as weight
				double *rWeight, double *bigExpR2,				// 2 in-out
				double *bigRprime, double *bigR,				// 2 in-out
// *nPairPtr = total pairs to be listed in Burrows file
				unsigned long long *nPairPtr,					// 1 in-out
// *npairTot = total loc. pairs, *npairSkip = number of loc. pairs skipped
				unsigned long long *npairTot, long *npairSkip,	// 2 in-out
				unsigned long long prompt,						// 1 in
					// (to inform the user after "prompt" pairs calculated)
				char sepBurOut, char moreCol, char BurAlePair,	// 3 in
// add 2 in-parameters in Apr 2015:
				struct chromosome *chromoList, int nChromo,		// 2 in
// add in March 2016:
				char jack, int **p1Gen, int **p2Gen, int *noDatFish,
				char *countm1, char *countm2,
				int *mValp1, float *freqp1, float *homop1,
				int *mValp2, float *freqp2, float *homop2, float *r2AtPairX,
// temporarily add for checking:
//double *r2JackTot,
				unsigned long long *r2Count,
				double *r2WRemSmp, float *JweighPair, double *JweightTot,
// temporarily add for checking:
//double *r2Ave,
//char *opened,
				 float epsilon)
{

	char BurrPause = 0;
	int p1, p2, nInd1, nInd2;
	int m, k1, k2, n;
	int pair12;
// add in Feb 2012
	float rdrift, expR2;
	float rB, nSamp;
	int nMpairs;
	float reNSamp, nIndtot, weight, rBweight;

	unsigned long long pairval;

    unsigned long long nLocPairs = 0;	// number of locus pairs, return value

	ALLEPTR allep1, allep2;
	FISHPTR popLoc1, popLoc2;

// temporarily add for checking
//*r2Ave = 0;

	// to inform the user after "prompt" pairs calculated;
	pairval = prompt;

// In the next "for" loop, pick a locus in the ascending order, another locus
// from the set of loci at the order after the first, then go through all
// allele in the mobility lists corresponding to this pair of loci. These pairs
// of loci are in the set of accepted pairs given by array okLoc, which was
// determined by function Loci_Eligible.
// Then calculate Burrows coefficients by function Burrows_Calcul.

// If we switch the order of the two "for" loops:
//		* "for (n = 0; n > m, m < nChromo; n++)"
//		* "for (k1 = 0; k1 < chromoList[m].nloci; k1++)"
// then the run will slightly faster, since we don't need to reassign
// variable p1 every time that n (index for second chromoList) is incremented
// That way, we would pair each locus in the first chromosome, chromoList[m],
// with every locus outside of chromoList[m], for each instance of m.
// However, doing as below allows us to get summaries (if desired) of pairing
// two chromosomes.
	for (m = 0; m < (nChromo-1); m++) {
		for (n = m+1; n < nChromo; n++) {
			// pairing loci in chromoList[m] and chromoList[n] --------------
			// For future modification:
			// We can calculate Burrows coefficient (weighted or not) average
			// over all pairs taken across chromoList[m] and chromoList[n],
			// leaving the overall average calculated at the calling function.
			// For example, pair12 is the number of pairs collected across
			// the two chromosomes.
			// For now, no attempt is made for others, pair12 is not used yet.
			pair12 = 0;
			for (k1 = 0; k1 < chromoList[m].nloci; k1++) {
				p1 = (chromoList[m].locus)[k1]; // "p" is increasing with "k"
				if (p1 > lastOK) break;	// quit when it passes last accepted one
				if (*(okLoc+p1) == 0) continue;
				allep1 = *(alleList+p1);
				popLoc1 = *(fishHead+p1);
				for (k2 = 0; k2 < chromoList[n].nloci; k2++) {
					p2 = (chromoList[n].locus)[k2];
					if (*(okLoc+p2) == 0) continue;
					if (p2 > lastOK) break;
					(*npairTot)++;
					allep2 = *(alleList+p2);
					popLoc2 = *(fishHead+p2);
					if (p1 - locSkip >= LOCBURR || p2 - locSkip >= LOCBURR) {
						BurrPause = 1;
					} else {
						BurrPause = 0;
						(*nPairPtr)++;
					}
					Burrows_Calcul (cutoff, allep1, allep2, popLoc1, popLoc2,
						p1, p2, *(nMobil+p1), *(nMobil+p2), nfish, &nSamp,
						&nInd1, &nInd2, &nMpairs, &rB, currPop, missptr,
						outBurr, outBurrName, moreBurr, BurrPause, &expR2,
						weighsmp, sepBurOut, moreCol, BurAlePair,
// temporarily added for checking:
//opened,
						jack, p1Gen, p2Gen, noDatFish, countm1, countm2,
						mValp1, freqp1, homop1, mValp2, freqp2, homop2,
						r2AtPairX,
// temporarily add for checking:
//r2JackTot,
						JweighPair, r2Count, epsilon);
					if (nMpairs <= 0) {
						(*npairSkip)++;
						continue;
					}

// temporarily add for checking
//*r2Ave += rB;

					AddBurrVal (nInd1, nInd2, rB, nSamp, expR2, weighsmp,
							locSkip, nLocPairs, rB2, rBdrift, prodInd,
							sampCount, pairWt, rAveTemp, totInd, wMeanSamp,
							rWeight, bigExpR2, bigRprime, bigR);
					if (jack != 0)
						JackWeight (weighsmp, nSamp, nfish, noDatFish,
							r2AtPairX, r2WRemSmp, JweighPair, JweightTot);

				// add this prompt to inform the user the progress:
					if ((nLocPairs) == pairval) {
						printf ("%18llu done, at loc. pair (%d, %d)\n",
								pairval, p1+1, p2+1);
						pairval +=prompt;
					}
            	// count the number of locus pairs nLocPairs.
            	// Increment happens after assignments above,
            	// to avoid going over the limit of the declared sizes
            	// of arrays prodInd, etc.
					nLocPairs++;
					pair12++;	// number of locus pairs from chromo[m], [n].
				}	// end of "for (k2 = 0; k2 < chromoList[n].nloci; k2++)"
			// thru all loc. pairs (p1, p2), with p2 in chromoList[n]
			}	// end of "for (k1 = 0; k1 < chromoList[m].nloci; k1++)"
		// thru all (p1, p2), p1 in chromoList[m], p2 in chromoList[n], n fixed
		}	// end of "for (n = 1; n > m, n < nChromo; n++)"
		// thru all (p1, p2), p1 in chromoList[m], p2 in chromoList[n], all n
	}	// end of "for (m = 0; m < (nChromo-1); m++)"
	return nLocPairs++;
}


//-------------------------------------------------------------------------
// Version of LDRunPairs, but locus pairs are taken within each chromosome

unsigned long long LDOneChromo (
				float cutoff, ALLEPTR *alleList, int currPop,	// 3 in
				int nfish, FISHPTR *fishHead, int *nMobil,		// 3 in
				int *missptr, int lastOK, char *okLoc,			// 3 in
				FILE *outBurr,									// 1 in-out
				char moreBurr, char *outBurrName,				// 2 in
// the next 5 are arrays for storing value/locus pair in case rAveTemp = null
				float *rB2, float *rBdrift, float *prodInd,		// 3 in-out
				float *sampCount, float *pairWt,				// 2 in-out
				char weighsmp, int locSkip,						// 2 in
				FILE *rAveTemp,									// 1 in-out
// because of rouding off errors, use double to calculate (Apr 2012):
				double *totInd,	// total ind. alle.				// 1 in-out
				double *wMeanSamp, // total wt. inv. of samp.size  1 in-out
//				float *rBAveTotal,
// currently, tot. weight *rWeight=*totInd, since ind.alle. is used as weight
				double *rWeight, double *bigExpR2,				// 2 in-out
				double *bigRprime, double *bigR,				// 2 in-out
// *nPairPtr = total pairs to be listed in Burrows file
				unsigned long long *nPairPtr,					// 1 in-out
// *npairTot = total loc. pairs, *npairSkip = number of loc. pairs skipped
				unsigned long long *npairTot, long *npairSkip,	// 2 in-out
				unsigned long long prompt,						// 1 in
					// (to inform the user after "prompt" pairs calculated)
				char sepBurOut, char moreCol, char BurAlePair,	// 3 in
// add 2 in-parameters in Apr 2015:
				struct chromosome *chromoList, int nChromo,		// 2 in-out
				char jack, int **p1Gen, int **p2Gen, int *noDatFish,
				char *countm1, char *countm2, int *mValp1, float *freqp1,
				float *homop1, int *mValp2, float *freqp2, float *homop2,
				float *r2AtPairX,
// temporarily add for checking:
//double *r2JackTot,
				unsigned long long *r2Count,
				double *r2WRemSmp, float *JweighPair, double *JweightTot,
// temporarily add for checking:
//double *r2Ave,
//char *opened,
				float epsilon)
{

	char BurrPause = 0;
	int p1, p2, nInd1, nInd2;
	int m, k1, k2;
	int pair12;
// add in Feb 2012
	float rdrift, expR2;
	float rB, nSamp;
	int nMpairs;
	float reNSamp, nIndtot, weight, rBweight;

	unsigned long long pairval;

    unsigned long long nLocPairs = 0;	// number of locus pairs, return value

	ALLEPTR allep1, allep2;
	FISHPTR popLoc1, popLoc2;
	// to inform the user after "prompt" pairs calculated;
	pairval = prompt;

// temporarily add for checking
//*r2Ave = 0;


// In the next "for" loop, pick two loci within a chromosome, then go
// through all allele in the mobility lists corresponding to this pair
// of loci. These pairs of loci are in the set of accepted pairs given
// by array okLoc, which was determined by function Loci_Eligible.
// Then calculate Burrows coefficients by function Burrows_Calcul.

	for (m = 0; m < nChromo; m++) {
		// pairing loci in chromoList[m] --------------------------------
		// For future modification:
		// We can calculate Burrows coefficient (weighted or not) average
		// over all pairs taken within chromoList[m], outputted to outBurr,
		// leaving the overall average calculated at the calling function.
		// For example, pair12 is the total pairs collected in chromoList[m],
		// For now, no attempt is made for others, pair12 is not used yet.
		pair12 = 0;
		for (k1 = 0; k1 < (chromoList[m].nloci - 1); k1++) {
			p1 = (chromoList[m].locus)[k1]; // "p" is increasing with "k"
			if (p1 > lastOK) break;	// quit when it passes last accepted one
			if (*(okLoc+p1) == 0) continue;
			allep1 = *(alleList+p1);
			popLoc1 = *(fishHead+p1);
			for (k2 = k1+1; k2 < chromoList[m].nloci; k2++) {
				p2 = (chromoList[m].locus)[k2];
				if (*(okLoc+p2) == 0) continue;
				if (p2 > lastOK) break;
				(*npairTot)++;
				allep2 = *(alleList+p2);
				popLoc2 = *(fishHead+p2);
				if (p1 - locSkip >= LOCBURR || p2 - locSkip >= LOCBURR) {
					BurrPause = 1;
				} else {
					BurrPause = 0;
					(*nPairPtr)++;
				}
				Burrows_Calcul (cutoff, allep1, allep2, popLoc1, popLoc2,
						p1, p2, *(nMobil+p1), *(nMobil+p2), nfish, &nSamp,
						&nInd1, &nInd2, &nMpairs, &rB, currPop, missptr,
						outBurr, outBurrName, moreBurr, BurrPause, &expR2,
						weighsmp, sepBurOut, moreCol, BurAlePair,
// temporarily added for checking:
//opened,
						jack, p1Gen, p2Gen, noDatFish, countm1, countm2,
						mValp1, freqp1, homop1, mValp2, freqp2, homop2,
						r2AtPairX,
// temporarily add for checking:
//r2JackTot,
						JweighPair, r2Count, epsilon);
				if (nMpairs <= 0) {
					(*npairSkip)++;
					continue;
				}

// temporarily add to check
//*r2Ave += rB;

				AddBurrVal (nInd1, nInd2, rB, nSamp, expR2, weighsmp,
							locSkip, nLocPairs, rB2, rBdrift, prodInd,
							sampCount, pairWt, rAveTemp, totInd, wMeanSamp,
							rWeight, bigExpR2, bigRprime, bigR);
				if (jack != 0)
					JackWeight (weighsmp, nSamp, nfish, noDatFish,
							r2AtPairX, r2WRemSmp, JweighPair, JweightTot);

				// add this prompt to inform the user the progress:
				if ((nLocPairs) == pairval) {
					printf ("%18llu done, at loc. pair (%d, %d)\n",
								pairval, p1+1, p2+1);
					pairval +=prompt;
				}
            	// count the number of locus pairs nLocPairs.
            	// Increment happens after assignments above,
            	// to avoid going over the limit of the declared sizes
            	// of arrays prodInd, etc.
				nLocPairs++;
				pair12++;	// number of locus pairs from chromoList[m]
			}	// end of "for (k2 = k1+1; k2 < chromoList[n].nloci; k2++)"
			// thru all loc. pairs (p1, p2), with p2 > p1 in chromoList[n]
		}	// end of "for (k1 = 0; k1 < chromoList[m].nloci; k1++)"
		// thru all (p1, p2), p1 < p2 in chromoList[m], m fixed
	}	// end of "for (m = 0; m < (nChromo-1); m++)"
	return nLocPairs++;
}


// --------------------------------------------------------------------------

void Pair_Analysis (float cutoff, ALLEPTR *alleList, int currPop,
						int nfish, FISHPTR *fishHead, int *nMobil,
						int *missptr, int lastOK, char *okLoc,
						double *nIndSum, float *rB2WAve, float *wHarmonic,
						float *wExpR2, FILE *outBurr,
//						FILE *outLoc, char moreDat,
						char moreBurr, char *outBurrName,
						unsigned long long *nBurrAve, float *rB2, float *rBdrift,
						float *prodInd, float *sampCount, float *pairWt,
						float *r2driftAve, float *totWeight, float *bigR2,
						float *bigRdrift, char weighsmp, FILE *rAveTemp,
// add in jan 2015
						char sepBurOut, char moreCol, char BurAlePair,
// temporarily added for checking:
//char *opened,
// add in Apr 2015
						struct chromosome *chromoList, int nChromo, char chroGrp,
// add in Jan, ..., 2016
						char jack, int *mValp1, float *freqp1, float *homop1,
						int *mValp2, float *freqp2, float *homop2,
						double *r2WRemSmp, unsigned long long *r2Count)
// rBdrift stores r2-drift for all locus pairs
// prodInd stores product of ind. alleles at locus pairs
// sampCount stores sample sizes for all locus pairs

{
// change to double (Apr 2012)
    double bigR = 0.0;	// total of rB taken over all locus pairs, weighted
						// by independent of allele pairs.
//    float rBAveTotal = 0.0;	// sum of r2-averaged values to be stored,
							// then retrieved later for Jackknife
	int i, p;
	float sign;
// because of rouding off errors, use these to calculate (Apr 2012):
	double wMeanSamp, totInd;	// weighted harmonic mean samp size, ind. alleles
	double rWeight, bigExpR2, bigRprime;
// add time to inform user:
//	time_t rawtime;
//	char info[15], *timeptr;
// added Sep.- Nov 2011 the following variables:
	int locSkip = 0;
//	char BurrPause = 0;
	long npairSkip = 0;
	unsigned long long prompt;
	unsigned long long nLocPairs = 0;
	unsigned long long nPairPtr = 0;
	unsigned long long npairTot = 0;
	unsigned long long maxpairs;

// added Mar 2016 -----------------------------------------------------
// The following arrays are to store info at a pair of loci, say (p1, p2).
// p1Gen, p2Gen are for genotypes of samples at two loci p1, p2
// noDatFish is to indicate which sample has no data, or both have data
// For sample (i+1)th:
// noDatFish[i] = 0: both loci have data,
//              = 1: no data at locus p1,
//              = 2: no data at locus p2,
//              = 3: both loci have no data.
// countm1[i] = number of copies of a particular allele m1 at locus p1,
// countm2[i] = number of copies of a particular allele m2 at locus p2,
// Allocating them here, not at Burrows_Calcul (indirectly called by this
// function) to decrease execution time, since Burrows_Calcul is called
// repeatedly (at each pair of loci). Memory allocations consume time!
	int **p1Gen = (int**) malloc(sizeof(int*)*nfish);
	int **p2Gen = (int**) malloc(sizeof(int*)*nfish);
	int *noDatFish = (int*) malloc(sizeof(int)*nfish);
	char *countm1 = (char*) malloc(sizeof(char)*nfish);
	char *countm2 = (char*) malloc(sizeof(char)*nfish);

// those are actually needed to declare here:
	float *r2AtPairX = (float*) malloc(sizeof(float)*nfish);
// temporarily add for checking:
//double *r2JackTot = (double*) malloc(sizeof(double)*nfish);

	float *JweighPair = (float*) malloc(sizeof(float)*nfish);
	double *JweightTot = (double*) malloc(sizeof(double)*nfish);
	for (i = 0; i < nfish; i++) {
		p1Gen[i] = (int*) malloc (2*sizeof(int));
		p2Gen[i] = (int*) malloc (2*sizeof(int));
		r2Count[i] = 0;
		r2AtPairX[i] = 0;
// temporarily add for checking:
//r2JackTot[i] = 0;
		JweighPair[i] = 0;
		JweightTot[i] = 0;
		r2WRemSmp[i] = 0;
	}

	float epsilon = (float) 8*nfish*nfish;
	epsilon = (1/epsilon);
	// for a rational number r (0 <= r) whose denominator <= 4*(nfish)^2,
	// r < epsilon implies r = 0.
// --------------------------------------------------------------------

//	ALLEPTR allep1, allep2;
//	FISHPTR popLoc1, popLoc2;

	*nBurrAve = 0;		// number of r2-averaged values
	*nIndSum = 0;		// Total of independent comparisons
	*rB2WAve = 0;
	*wHarmonic = 0;
	*totWeight = 0;
	*wExpR2 = 0;
	*r2driftAve = 0;
// these are double precision, will be assigned to parameters when done
	wMeanSamp = 0;
	totInd = 0;
	rWeight = 0;
/*
	if (outLoc != NULL && (moreDat == 1)) {
		for (i=0; i<40; i++) fprintf (outLoc, "-");
		fprintf (outLoc, "\n\n");
		fprintf (outLoc, "\nFor Linkage Disequilibrium Method\n");
		for (i=0; i<45; i++) fprintf (outLoc, "-");
		fprintf (outLoc, "\n   Locus pair    N. Ind. pair   #Samp. w/data\n");
		fflush (outLoc);
	}
//*/
// change in Nov 2014-Mar 2015
// Bring labeling from Burrows_Calcul for separate file here, so that only
// label for all pairs of loci
//    if (outBurr != NULL && moreBurr == 1 && BurrPause == 0
    if (outBurr != NULL && moreBurr == 1 && BurAlePair == 1 && sepBurOut == 1)
    {
        	fprintf (outBurr, "Loc._Pairs   Allele_Pairs    P1"
				"    P2    Burrows->D       r         r^2\n");
		fflush (outBurr);
    }
//	if (outBurr != NULL && moreBurr == 1 && BURRSHORT == 1) {
	if (outBurr != NULL && moreBurr == 1 && BurAlePair == 0) {
//		if (sepBurOut == 1 && NOEXPLAIN != 1) {
// Use "or" so that explanations always given when only one file outputted
		if (sepBurOut != 1 || NOEXPLAIN != 1) {
			fprintf (outBurr, "\n> LowP1/LowP2: Lowest allele freq. at Loc1/Loc2 if more than one allele used,\n"
							  "               = (1 - q) if only one allele is considered, whose freq. = q\n");
		}
		if (moreCol == 0)
		{
			if (sepBurOut != 1 || NOEXPLAIN != 1) fprintf (outBurr, "\n");
        	fprintf (outBurr, "  Loc1   Loc2   LowP1   LowP2  Samp.Size    Mean_r^2     r^2_drift\n");
//			if (NOEXPLAIN != 1) {
			if (sepBurOut != 1 || NOEXPLAIN != 1) {
				for (i=0; i<68; i++) fprintf (outBurr, "-");
				fprintf (outBurr, "\n");
			}
		} else {
//			if (sepBurOut == 1 && NOEXPLAIN != 1) {
			if (sepBurOut != 1 || NOEXPLAIN != 1) {
				fprintf (outBurr, "> Ind1/Ind2: Number of independent alleles in Loc1/Loc2\n");
				fprintf (outBurr, "> #Pairs: Number of allele pairs considered in (Loc1, Loc2)\n\n");
			}
			fprintf (outBurr, "  Loc1   Loc2   LowP1   LowP2  Ind1 Ind2  #Pairs  Samples   ");
			fprintf (outBurr, "Mean_D        Mean_r        Mean_r^2     r^2_drift\n");
//			if (NOEXPLAIN != 1) {
			if (sepBurOut != 1 || NOEXPLAIN != 1) {
				for (i=0; i<111; i++) fprintf (outBurr, "-");
				fprintf (outBurr, "\n");
			}
		}
		fflush (outBurr);
	}

	locSkip = 0;	// count number of the first LOCBURR loci,
					// not counting loci skipped. This is used for printing
					// pairs of loci using the first LOCBURR loci.
	i = LOCBURR - 1;
	for (p=0; (p<lastOK); p++) {
		if (*(okLoc+p) == 0) {
			locSkip++;
			i++;
		}
		if (p >= i) break;
	}
//printf ("last Locus = %d, loci skipped = %d\n", lastOK+1, locSkip);
// add Apr 2012 for putting info to the console
	prompt = (unsigned long long) 1000000;	// inform the user after prompt pairs calculated
	maxpairs = (unsigned long long) (lastOK-locSkip);
	maxpairs *= (maxpairs+1);
	maxpairs /= 2;
	if (maxpairs > prompt) {
		printf ("     Calculating r^2");
		if (chroGrp == 0 || nChromo <= 1)
			printf (" (at most %llu values)", maxpairs);
		printf (":\n");
	}
// In the next "for" loop, pick a locus in the ascending order, another locus
// from the set of loci at the order after the first, then go through all
// allele in the mobility lists corresponding to this pair of loci. These pairs
// of loci are in the set of accepted pairs given by array okLoc, which was
// determined by function Loci_Eligible.
// Then calculate Burrows coefficients by function Burrows_Calcul. The outputs
// from this function are then averaged over all pairs.

// add Feb 2012
	*bigRdrift = 0;
	*bigR2 = 0;
	bigExpR2 = 0;
	bigRprime = 0;
// next few lines are for checking, comment out later
/*
double r2Ave;
FILE *checkR2;	// the file will be named "checkR2.txt"
char toChk = 1;	// if this = 1, will check
int popChk = 1;	// only check for one population
//if (currPop != popChk) toChk = 0;
if (currPop != popChk || jack == 0) toChk = 0;
if (toChk != 0) {
if (*opened == 0) {
checkR2 = fopen ("checkR2.txt", "w");
*opened = 1;
}
else checkR2 = fopen ("checkR2.txt", "a");
}
//*/

//	info[15] = '\0';
// Added Apr 2015
// If there is only one chromosome (nChromo = 1), then all pairs are taken
	if (chroGrp > 0 && nChromo > 1) {
		if (chroGrp == 1) {
			printf ("       Loci are paired within each chromosome\n");
			nLocPairs = LDOneChromo (cutoff, alleList, currPop, nfish,	//4
							fishHead, nMobil, missptr, lastOK, okLoc,	//5
							outBurr, moreBurr, outBurrName, rB2,		//4
							rBdrift, prodInd, sampCount, pairWt,		//4
							weighsmp, locSkip, rAveTemp, &totInd,		//4
							&wMeanSamp, &rWeight, &bigExpR2, &bigRprime,//4
							&bigR, &nPairPtr, &npairTot, &npairSkip,	//4
							prompt, sepBurOut, moreCol, BurAlePair,		//4
							chromoList, nChromo,						//2
							jack, p1Gen, p2Gen, noDatFish, countm1,
							countm2, mValp1, freqp1, homop1, mValp2,
							freqp2, homop2, r2AtPairX,
// temporarily add for checking:
//r2JackTot,
							r2Count,
							r2WRemSmp, JweighPair, JweightTot,
// temporarily add for checking:
//&r2Ave, opened,
							epsilon);
		} else {
			printf ("       Loci are paired across chromosomes\n");
			nLocPairs = LDTwoChromo (cutoff, alleList, currPop, nfish,
							fishHead, nMobil, missptr, lastOK, okLoc,
							outBurr, moreBurr, outBurrName, rB2,
							rBdrift, prodInd, sampCount, pairWt,
							weighsmp, locSkip, rAveTemp, &totInd,
							&wMeanSamp, &rWeight, &bigExpR2, &bigRprime,
							&bigR, &nPairPtr, &npairTot, &npairSkip,
							prompt,	sepBurOut, moreCol, BurAlePair,
							chromoList, nChromo,
							jack, p1Gen, p2Gen, noDatFish, countm1,
							countm2, mValp1, freqp1, homop1, mValp2,
							freqp2, homop2, r2AtPairX,
// temporarily add for checking:
//r2JackTot,
							r2Count,
							r2WRemSmp, JweighPair, JweightTot,
// temporarily add for checking:
//&r2Ave, opened,
							epsilon);
		}
	} else {
		nLocPairs = LDRunPairs (cutoff, alleList, currPop, nfish, fishHead,	//5
						nMobil, missptr, lastOK, okLoc, outBurr, moreBurr,	//6
						outBurrName, rB2, rBdrift, prodInd, sampCount,		//5
						pairWt, weighsmp, locSkip, rAveTemp, &totInd,		//5
						&wMeanSamp, &rWeight, &bigExpR2, &bigRprime, &bigR,	//5
						&nPairPtr, &npairTot, &npairSkip, prompt,			//4
						sepBurOut, moreCol, BurAlePair,
						jack, p1Gen, p2Gen, noDatFish, countm1, countm2,
						mValp1, freqp1, homop1, mValp2, freqp2, homop2,
						r2AtPairX,
// temporarily add for checking:
//r2JackTot,
						r2Count,
						r2WRemSmp, JweighPair, JweightTot,
// temporarily add for checking:
//&r2Ave, opened,
						epsilon);
	}

// --------------------------------------------------------------------
// These are for checking r^2-calculations for sample sets minus one.
// (Put r^2-outputs in a file.)
// Add or delete a slash "/" at the line before the beginning of code
// to cover or open the code block (1 slash: cover; 2 slashes: open):
/*

if (toChk != 0) {
if (nLocPairs > 0) r2Ave /= (nLocPairs);
fprintf (checkR2,
"\nPopulation %d, #samples = %d, Pcrit =%6.3f, up to locus %d\n"
"              Number of Pairs =%9llu,  (unweighted mean) r^2 = %9.6f\n\n"
"Table for unweighted and weighted r^2 when one individual is removed\n\n",
currPop, nfish, cutoff, lastOK+1, nLocPairs, r2Ave);
fprintf(checkR2, "Remv  N. of r^2    Total Wt      r^2-Sum     Wt-r^2-Sum    r^2Ave    Wt-r^2Ave\n");
}
//*/ //----------------------------------------------------------------

// added Mar 2016 ----------------
	for (i = 0; i < nfish; i++) {
		free (p1Gen[i]);
		free (p2Gen[i]);
// for checking, comment out later:
/*
if (toChk != 0) fprintf(checkR2,
"%3d%11llu%13.0f%13.3f%15.3f", // 2nd format: "eleven ell ell u"
i+1, r2Count[i], JweightTot[i], r2JackTot[i], r2WRemSmp[i]);
//*/
		if (JweightTot[i] != 0) {
// temporarily add for checking:
//r2JackTot[i] /= (float)r2Count[i];
			r2WRemSmp[i] /= JweightTot[i];
		}
// for checking, comment out later:
//if (toChk != 0) fprintf(checkR2, "%11.6f%12.7f\n", r2JackTot[i], r2WRemSmp[i]);
	}

	free (p1Gen);
	free (p2Gen);
	free (noDatFish);
	free (countm1);
	free (countm2);
// temporarily add for checking:
//free (r2JackTot);
	free (r2AtPairX);
	free (JweighPair);
	free (JweightTot);

// -------------------------------
	if (nLocPairs == 0) return;
	*nBurrAve = nLocPairs;
	*nIndSum = totInd;
	*totWeight = (float) rWeight;
	if (*nIndSum > 0) {
		*bigR2 = (float) bigR;
		*bigRdrift = (float) bigRprime;
// this is used for calculating Ne from r^2
		bigR /= rWeight;
// this is used for calculating Ne from r^2-(exp. r^2): changed in Feb 2012
		bigRprime /= rWeight;
		if (wMeanSamp > 0) wMeanSamp = totInd/wMeanSamp;
// change in Feb 2012: calculate expected R^2-sample as weighted average
// of all R^2 in pairs:
		bigExpR2 /= rWeight;

		*rB2WAve = (float) bigR;
		*wExpR2 = (float) bigExpR2;
		*r2driftAve = (float) bigRprime;
		*wHarmonic = (float) wMeanSamp;
	}


// next 3 lines for checking, comment out later:
/*
if (toChk != 0) {
fprintf(checkR2, "\nr^2 Weighted for whole sample = %10.6f\n", *rB2WAve);
fflush (checkR2);
fclose (checkR2);
}
//*/

//    *wExpR2 = ExpR2Samp(*wHarmonic);
	if (maxpairs > prompt)
		printf("     Actual number of r^2-values evaluated = %llu\n", *nBurrAve);
    // write the total weight, the sum of ave. r2-values of allele pairs,
	// unweighted and weighted.
	// These will be the last records for those temporary files, which will be
	// used for doing Jackknife.
	// write 0 when done with rAveTemp to mark the end of this file
	sign = 0;
	if (rAveTemp != NULL && (*nBurrAve) > 0) {
		fwrite (&sign, sizeof(float), 1, rAveTemp);
//		fwrite (&rBAveTotal, sizeof(float), 1, rAveTemp);
//		fwrite (&bigR, sizeof(float), 1, rAveTemp);
	}
/*
	// write to auxiliary file (for locus data) as needed
    if (outLoc != NULL && (moreDat == 1)) {
		for (i=0; i<45; i++) fprintf (outLoc, "-");
		fprintf (outLoc,"\n\n");
		if (nLocPairs > nPairPtr) fprintf (outLoc,
				"Only %lu accepted locus pairs are listed, up to locus %d\n",
				nPairPtr, LOCOUTPUT+locSkip);
		fprintf (outLoc,"Total locus pairs investigated =%13lu\n", npairTot);
		fprintf (outLoc,"    * Number of pairs rejected =%13lu\n", npairSkip);
		fprintf (outLoc,"    * Number of pairs accepted =%13lu\n", nLocPairs);
		fprintf (outLoc,
		"\nWeighted Harmonic Mean Sample Size =%9.2f\n(by Ind. Alleles)\n", *wHarmonic);
//		fprintf (outLoc, "Weighted Mean of Expected R^2 Sample =    %11.6f\n\n", *wExpR2);
		for (i=0; i<45; i++) fprintf (outLoc, "*"); fprintf (outLoc, "\n");
		fflush (outLoc);
	}
//*/
	if (outBurr != NULL && moreBurr == 1 &&
	// change in Nov 2014: add condition
		(sepBurOut != 1)) {
		if (nLocPairs > nPairPtr) fprintf (outBurr,
				"\nOnly %llu accepted locus pairs are listed, up to locus %d",
				nPairPtr, LOCBURR+locSkip);
		fprintf (outBurr,"\nTotal locus pairs investigated =%13llu\n", npairTot);
		if ((*nBurrAve) == 0) {
			fflush (outBurr);
			return;
		}
		fprintf (outBurr,"    * Number of pairs rejected =%13lu\n", npairSkip);
		fprintf (outBurr,"    * Number of pairs accepted =%13llu\n", nLocPairs);
		fprintf (outBurr,
				"\nWeighted (by Ind. Alleles) Harmonic Mean Sample Size =%11.2f\n",
				*wHarmonic);
		fprintf (outBurr,
			"Expected R^2-sample calculated from this sample size = %10.6f\n",
				ExpR2Samp(*wHarmonic));
		fprintf (outBurr, "\n# Weighted Mean of r^2 =%22.6f\n", *rB2WAve);
		fprintf (outBurr, "# Weighted Mean of Exp. r^2 Sample =%10.6f\n", *wExpR2);
		fprintf (outBurr, "# Weighted Mean of r^2-drift =%16.6f  (%11.3e), ",
				*r2driftAve, *r2driftAve);
	// Note: this last output line does not end with a new line, since a value
	// for Ne will be added after the call for this function, under the same
	// conditions on outBurr here
		fflush (outBurr);
	}

}


// --------------------------------------------------------------------------

float LDmethod (float cutoff, ALLEPTR *alleList, int popRead, int samp,
				FISHPTR *fishHead, int *nMobil, int *missptr, int lastOK,
				char *okLoc, double *nIndSum, float *rB2WAve, float *r2driftAve,
				float *wHarmonic, float *wExpR2, FILE *outBurr, FILE *outLoc,
				char moreDat, char moreBurr, char *outBurrName, char mating,
				float infinite, char param, char jacknife, char *jackOK,
				float *confJacklow, float *confJackhi, long *Jdegree,
				float *confParalow, float *confParahi, char weighsmp,
				int *memOut, int icount,
// added in Jan/Mar 2015
				char sepBurOut, char moreCol, char BurAlePair,
// temporarily added for checking with checkR2:
// char *opened,
// add in Apr 2015
				struct chromosome *chromoList, int nChromo, char chroGrp)
{
	FILE *rAveTemp = NULL;
	FILE *weighFile = NULL;
	char tmpUsed = USETMP;
// add Feb 2012:
	int j, k;
	float *rB2=NULL, *rBdrift=NULL, *sampCount=NULL, *prodInd=NULL, *pairWt=NULL;
//	float *rBdrift, *sampCount, *prodInd, *pairWt;
//	float r2driftAve;
	float totW, totR2, totRdrift;
	float r2, r2drift, weight, dummy;
	unsigned long long nB;
	unsigned long long nBurrAve = 0;	// to count the total number of r^2-values
						// which will be used for jacknife.
	char modify;
	char drift = 0;		// if = 1, use r2-drift for CI; otherwise, use r2.
	float estNe = 0;
//	float iniNe = 0;
	*memOut = 0;

// add in Mar 2016:
	int maxNAlle = 0;
	for (k = 0; k < lastOK; k++) {
		if (*(nMobil+k) > maxNAlle) maxNAlle = *(nMobil+k);
	}
// crash sometimes without this increase:
	maxNAlle++;
	int *mValp1 = (int*) malloc(sizeof(int)*maxNAlle);
	float *freqp1 = (float*) malloc(sizeof(float)*maxNAlle);
	float *homop1 = (float*) malloc(sizeof(float)*maxNAlle);
	int *mValp2 = (int*) malloc(sizeof(int)*maxNAlle);
	float *freqp2 = (float*) malloc(sizeof(float)*maxNAlle);
	float *homop2 = (float*) malloc(sizeof(float)*maxNAlle);

// For Jackknife on samples:
// Sk represents the sample set S with the (k+1)th individual removed
// r2Count[k]: count the number of eligible pairs of loci in Sk
	unsigned long long *r2Count =
		(unsigned long long*) malloc(sizeof(unsigned long long)*samp);
// r2WRemSmp[k] is the weighted r2 for sample set Sk
	double *r2WRemSmp = (double*) malloc(sizeof(double)*samp);
// Initialize those 2 arrays will be done in Pair_Analysis

// Give an early estimate how many r^2-values to considered to set the size of arrays
	for (j=0; (j<lastOK); j++) {
		if (*(okLoc+j) == 0) continue;	// marked by Loci_Eligible to be skipped.
		for (k=j+1; (k<=lastOK); k++) {
			if (*(okLoc+k) == 0) continue;	// locus (k+1) is skipped.
			nBurrAve++;
		}
	}
// print to console LD is started:
//	printf ("     Linkage Disequilibrium Method\n");
	if (tmpUsed == 1) {
		if (((rAveTemp = tmpfile()) == NULL) || ((weighFile = tmpfile()) == NULL))
		{
			printf ("     The System does not allow creating temporary file. RAM is used\n");
			tmpUsed = 0;
		}
	}

	if (tmpUsed == 0) {
// for storing rB-values Ind. alleles and sample sizes of pairs of loci.
// array pairWt to eliminate some calculations
		if ((rB2 = (float*) calloc (nBurrAve, sizeof(float))) == NULL ||
			(rBdrift = (float*) calloc (nBurrAve, sizeof(float))) == NULL ||
			(prodInd = (float*) calloc (nBurrAve, sizeof(float))) == NULL ||
			(sampCount = (float*) calloc (nBurrAve, sizeof(float))) == NULL ||
			(pairWt = (float*) calloc (nBurrAve, sizeof(float))) == NULL)
		{
			printf ("Out of memory for doing LD method at c = %5.3f!\n", cutoff);
			if (rB2 != NULL) free (rB2);
			if (rBdrift != NULL) free (rBdrift);
			if (prodInd != NULL) free (prodInd);
			if (sampCount != NULL) free (sampCount);
			*nIndSum = 0;
			*rB2WAve = 0;
			*wHarmonic = 0;
			*memOut = 1;
			*jackOK = 0;
			return 0;
		}
		for (nB=0; nB<nBurrAve; nB++) {
			*(rB2) = 0;
			*(rBdrift+nB) = 0;
			*(prodInd+nB) = 0;
			*(sampCount+nB) = 0;
			*(pairWt+nB) = 0;
		}
	}
// don't do jackknife when not needed
	if (*jackOK == 0) jacknife = 0;
	nBurrAve = 0;	// reset this, which will be calculated correctly
					// in the next function
	Pair_Analysis (cutoff, alleList, popRead, samp, fishHead, nMobil,
						missptr, lastOK, okLoc, nIndSum, rB2WAve,
						wHarmonic, wExpR2, outBurr, // outLoc, moreDat,
						moreBurr, outBurrName, &nBurrAve, rB2, rBdrift,
						prodInd, sampCount, pairWt, r2driftAve,
						&totW, &totR2, &totRdrift, weighsmp, rAveTemp,
						sepBurOut, moreCol, BurAlePair,
// temporarily add for checking with checkR2:
//opened,
// add in Apr 2015
						chromoList, nChromo, chroGrp,
// add in Mar 2016:
						jacknife, mValp1, freqp1, homop1,
						mValp2, freqp2, homop2, r2WRemSmp, r2Count);

	free (mValp1);
	free (freqp1);
	free (homop1);
	free (mValp2);
	free (freqp2);
	free (homop2);

// this calculates Ne from (r^2 - (exp r^2-sample)), changed in Feb 2012
	estNe = LD_Ne (*wHarmonic, *r2driftAve, mating, infinite);
//	iniNe = estNe;
	if (outBurr != NULL && moreBurr == 1 && (sepBurOut != 1))
		fprintf (outBurr, "        Ne =%10.1f\n", estNe);
	j = 0;
	// weighsmp > 0 when there are missing data
//	if (weighsmp > 0) {
	if (weighsmp > 0 && RESETNE != 0) {
// print to console initial estimate:
// recalculate with adjusted weights based on estNe above
		if (tmpUsed == 1) {
			rewind (rAveTemp);
// Remove parameter weighFile
//			j = NeAdjustedTmp (rAveTemp, weighFile,
			j = NeAdjustedTmp (rAveTemp,
// temporarily add for checking with checkR2
//*opened,

				nBurrAve, *wHarmonic, mating, infinite, &estNe, r2driftAve,
				&totW, &totR2, &totRdrift, wExpR2, rB2WAve);
		} else {
			j = NeAdjustedArr (pairWt, rB2, rBdrift, prodInd, sampCount, nBurrAve,
				*wHarmonic, mating, infinite, &estNe, r2driftAve,
				&totW, &totR2, &totRdrift, wExpR2, rB2WAve);
		}
//		j = NeRevisedTmp (tmpUsed, pairWt, rB2, rBdrift, prodInd, sampCount,
//					rAveTemp, weighFile,
//					nBurrAve, *wHarmonic, mating, infinite, &estNe,
//					&r2driftAve, &totW, &totR, wExpR2, rB2WAve);
//		if (j > 0 && icount == 0) {
		// info for Ne estimate is printed in the calling function above
//		if (j > 0) {
//			printf ("     Initial estimate of Ne: %12.1f\n", iniNe);
//			printf ("     Final estimate of Ne: %14.1f\n", estNe);
//		}
	}
	if (j == 0) {	// there is no attempt to reweight
		// icount = 0 when the program does not run with multiple files
//		if (icount == 0)
		printf ("       Estimate of Ne: %20.1f\n", estNe);

// this part is blocked out since jackknife on loci is no longer used:
/*
		if (rAveTemp != NULL) {
			rewind (rAveTemp);
			for (nB=0; nB<nBurrAve; nB++) {
				fread (&dummy, sizeof(float), 1, rAveTemp);
				fread (&dummy, sizeof(float), 1, rAveTemp);
				fread (&r2, sizeof(float), 1, rAveTemp);
				fread (&r2drift, sizeof(float), 1, rAveTemp);
				fread (&weight, sizeof(float), 1, rAveTemp);
				fwrite (&r2, sizeof(float), 1, weighFile);
				fwrite (&r2drift, sizeof(float), 1, weighFile);
				fwrite (&weight, sizeof(float), 1, weighFile);
			}
		}
//*/
	}
	if (outBurr != NULL && moreBurr == 1 && j != 0 &&
	// change in Nov 2014: add conditions on NOEXPLAIN
//		(NOEXPLAIN != 1)) {
		(sepBurOut != 1)) {
		if (weighsmp > 0) fprintf (outBurr,
			"\nWeights on locus pairs are revised on the initial estimate Ne\n");
		fprintf (outBurr, "# Weighted Mean of r^2 =%22.6f\n", *rB2WAve);
		fprintf (outBurr, "# Weighted Mean of Exp. r^2 Sample =%10.6f\n", *wExpR2);
		fprintf (outBurr,
			"# Weighted Mean of r^2-drift =%16.6f  (%11.3e), Revised Ne =%10.1f\n",
				*r2driftAve, *r2driftAve, estNe);
	}

	modify = 0;
	if (param==1)
	{
		LDConfidInt95 (*wHarmonic, samp, *wExpR2, *rB2WAve,
				*nIndSum, r2WRemSmp, r2Count,
// param opened is temporarily added for checking: print to file checkR2
//*opened,
				modify, confParalow, confParahi, Jdegree, infinite,
				mating, 0, moreBurr, outBurr);
/* Comment out this call for Jackknife on loci:
		LDConfidInt (drift, nBurrAve, *wHarmonic, *wExpR2, *rB2WAve, *r2driftAve,
					*nIndSum, modify, confParalow, confParahi, infinite,
					mating, 0, pairWt, rB2, rBdrift, tmpUsed, weighFile,
					totW, totR2, totRdrift, moreBurr, outBurr);
//*/
// print to console when not running multiple files::
		if (icount == 0) {
			printf ("     Parameter CI: ");
			if (*confParalow < 0 || *confParalow >= infinite)
				printf("%15s", "infinite");
			else printf("%15.1f", *confParalow);
			if (*confParahi < 0 || *confParahi >= infinite)
				printf("%16s\n", "infinite");
			else printf("%16.1f\n", *confParahi);
		}
	}
	*confJacklow = *confParalow;
	*confJackhi = *confParahi;
	modify = (param==1 && MERGE==1)? 1: 0;

	if (jacknife==1)
	{
		LDConfidInt95 (*wHarmonic, samp, *wExpR2, *rB2WAve,
				*nIndSum, r2WRemSmp, r2Count,
// param opened is temporarily added for checking: print to file checkR2
//*opened,
				modify, confJacklow, confJackhi, Jdegree, infinite,
				mating, 1, moreBurr, outBurr);
/* Comment out this call for Jackknife on loci:
		LDConfidInt (drift, nBurrAve, *wHarmonic, *wExpR2, *rB2WAve, *r2driftAve,
					*nIndSum, modify, confJacklow, confJackhi, infinite,
					mating, 1, pairWt, rB2, rBdrift, tmpUsed, weighFile,
					totW, totR2, totRdrift, moreBurr, outBurr);
//*/
// print to console when not running multiple files::
		if (icount == 0) {
			printf ("     Jackknife CI: ");
			if (*confJacklow < 0 || *confJacklow >= infinite)
				printf("%15s", "infinite");
			else printf("%15.1f", *confJacklow);
			if (*confJackhi < 0 || *confJackhi >= infinite)
				printf("%16s\n", "infinite");
			else printf("%16.1f\n", *confJackhi);
		}
	}
	free (r2Count);
	free (r2WRemSmp);
	if (outBurr != NULL && moreBurr == 1) fprintf (outBurr, "\n");
	if (tmpUsed == 0) {
		free (rB2);	// ok to free, but not survived in pop2
		free (rBdrift);	// ok to free, but not survived in pop2
		free (prodInd);	// crash when free
		free (sampCount);	// ok to free together with rBdrift
		free (pairWt);	// not survived to the end of pop 1 with others, but go alone to pop2
	}
	if (rAveTemp != NULL) fclose (rAveTemp);
	if (weighFile != NULL) fclose (weighFile);
	return estNe;
}


//---------------------------------------------------------------------------


//*************************************************************************
// Nomura
//
// For Molecular Coancestry Method, from the paper:
// "Estimation of Effective Number of Breeders From Molecular Coancestry
// of Single Cohort Sample," Tetsuro Nomura, Evolutionary. Appl. March 2008.

// In the code below, parameters fishList, mobilList, alleList are used.
// The mobilList and alleList are similar, only one is actually used
// in the main program depending on the number of loci "nloci" given in input.

// fishList is allocated as nloci addresses. Address (fishList+p) points to
// the list of samples at locus (p+1); each sample is represented as a node
// having field as the sample genotype. If "node" begins from the top of
// the list pointed to by (fishList+p), the assignment node = node->next
// will allow us to access genotypes at locus (p+1) of the whole population.

// Similarly, (mobilList+p) and (alleList+p) point to the list of alleles
// at locus (p+1). Each node in either list has fields: "mValue" for allele
// and "freq" for frequency of the corresponding allele

// nMobil[] is the array of the number of alleles at loci: nMobil[p] = number
// of alleles at locus (p+1). If it is 0, no sample has data at locus (p+1).



//------------------------------------------------------------------
void RemSampLoc(FISHPTR *fishList)
{
	FISHPTR curr = *fishList, temp;
	for (; curr != NULL; curr = temp) {
		temp = curr->next;
		curr->next = NULL;
		free (curr);
	}
	free (fishList);

}

//------------------------------------------------------------------

NONSIBPTR MakeNonsib (int samp1, int samp2)
{
	NONSIBPTR newptr=NULL;

	int size = sizeof(struct nonsib);
	if ((newptr = (NONSIBPTR) malloc(size)) != NULL)
	{
		newptr->first = samp1;
		newptr->second = samp2;
//		newptr->simVal = sim;
		newptr->next = NULL;
	}
	return newptr;
}
//------------------------------------------------------------------
void RemNonSib(NONSIBPTR *nonsibList)
{

	NONSIBPTR temp, curr = *nonsibList;
	*nonsibList = NULL;
	for (; curr != NULL; curr = temp) {
		temp = curr->next;
		curr->next = NULL;
		free (curr);
	};
	free (nonsibList);
}
//------------------------------------------------------------------
void RemFirstNode(NONSIBPTR *nonsibList, NONSIBPTR nodeRem)
// remove the fist node in nonsibList, assumed to be node nodeRem,
// that is, nonsibList initially points to node nodeRem.
{

	*nonsibList = nodeRem->next;
	nodeRem->next = NULL;
	free (nodeRem);
}
//------------------------------------------------------------------
void RemMextNode(NONSIBPTR nonsibNode, NONSIBPTR nodeRem)
// remove node nodeRem, assumed to be the node right after nonsibNode
{

	nonsibNode->next = nodeRem->next;
	nodeRem->next = NULL;
	free (nodeRem);
}
//------------------------------------------------------------------

COANPTR MakeCoan (int p, float s, float f, float wp, float freq2)
{
	COANPTR newptr;

	int size = sizeof(struct molecoef);
	if ((newptr = (COANPTR) malloc(size)) != NULL)
	{
		newptr->locus = p;
		newptr->scoan = s;
		newptr->fresq = freq2;
		newptr->diffcoan = f;
		newptr->weight = wp;
		newptr->next = NULL;
	}
	return newptr;
};
//------------------------------------------------------------------
void RemCoan(COANPTR *coanList)
{

	COANPTR temp, curr = *coanList;
	for (; curr != NULL; curr = temp) {
		temp = curr->next;
		curr->next = NULL;
		free (curr);
	};
	free (coanList);
}
//------------------------------------------------------------------

int SimilarInd (int xId[2], int yId[2], int defval, char *hasdat)
// Calculate molecular coancestry fxy (molecular similarity index)
// between two individuals x, y, using formula (4) in the paper,
// However, we don't divide by 4, but return value as integer.
// The return value is from 0 to 4, or the default value defval.

{
    int i, j;
    int k = 0;
    *hasdat = 1;
	// If one has no data, return default value defval:
	if (xId[0] == 0 || yId[0] == 0) {
		*hasdat = 0;
		return defval;
	};
    for (i=0; i<2; i++)
        for (j=0; j<2; j++)
            if (xId[i] == yId[j]) k++;
    return k;
}

//------------------------------------------------------------------


char AddSampRow(FISHPTR *sampTop, FISHPTR *sampTail, int sample[])
{

	FISHPTR newfish;
	if ((newfish = MakeFish (sample)) == NULL) return 0;
//printf ("sample = %d %d\n", newfish->gene[0],newfish->gene[1]);
	if (*sampTop == NULL) {
		*sampTop = newfish;
	} else (*sampTail)->next = newfish;
	*sampTail = newfish;
	return 1;

}

//------------------------------------------------------------------


/*
char AddSampCount(int *nfish, FISHPTR *sampTop, FISHPTR *sampTail, int sample[])
{

	FISHPTR newfish;
	if ((newfish = MakeFish (sample)) == NULL) return 0;
	if (*nfish == 0) {
		*sampTop = newfish;
	} else (*sampTail)->next = newfish;
	*sampTail = newfish;
	*nfish +=1;
	return 1;

}
//------------------------------------------------------------------
*/


float PutativeNonSib (NONSIBPTR *nonsibList, NONSIBPTR *nonsibTail,
					  int i, int *jmin, int p, int *npairs, float *ctotal,
					  FISHPTR *fishList, int nMobil[], int nloci,
					  char *okLoc, int nSamp, char *gotNoSib,
					  int *sibNodes, char *errcode,
					  FILE *outLoc, char moreDat, char detail)
//					  int *sibNodes, char *errcode)
// At sample i (i=0, ..., nSamp-1), pick another sample so that the pair
// can be assumed a putative nonsib pair. Then add a new node to nonsibList
// to represent this new putative nonsib. This node will be the (i+1)th node
// added to the list, field "first" of this node is i, field "second"
// represents the sample that sample i is paired with.
// At each locus q, starting position at the first sample, go thru all
// samples at that locus, indexing those samples by j=0, ..., (nSamp-1).
// Only choose j != i and such that the pair {i, j} is not identical to
// any pairs {m, n} in the NONSIBPTR nonsibList.

// Strategy for choosing j:
// Upon entering this function, nonsibList consists of at most (i-1) nodes
// (we say "at most" since some nodes will be deleted along the way), and
// for a node whose "first field" = m, its field "second" = n, and we have
// pair {m,n} where m != n. Each such pair should NOT be the same as pair
// {i,j} (order does not count, e.g., {0,2} and {2,0} are identical).
// To see which j should be excluded, we note that i > all such m, so need
// to see if i is one of the n's, i.e., check to see if i is in the field
// "second" of any node in nonsibList.
// * If it is NOT, then j can be any except i.
// * If it is, identify those nodes m whose "second" fields = i. Then j != i
// and j cannot be any of those m.
//
// For each such pair (i, j), calculate the average of similarity indices
// across all loci q !=p, cf. formula (8) in the paper. Then among those j,
// choose the j=jmin that results in the smallest index value.
// Then add the ith node to nonsibList with its "second" field = jmin.
// The pair (i,jmin) is assumed to be a putative nonsib pair at locus (p+1).

// To be more efficient in searching in nonsibList, we adopt the strategy:
// * Only add node to nonsibList for the pair (i, jmin) if jmin > i
// * Delete all nodes whose "second" field <= i at the end of this function,
// since those nodes are not needed in the search for the one "second" with
// sample (i+1) and beyond. Thus, we should delete the node whose second
// field is i (in this function), since it is irrelevant in the next i.

// Along the way, we also calculate the sum of molecular coancestries for all
// pairs (i,j): j>i at this locus (p+1), as a return value *ctotal, together
// with the number of pairs *npairs.
// The return value for the function is (the coancestry for the newly
// putative nonsib pair), of which sample i is one. Set parameter gotNoSib = 1
// if such pair is found.



// Return 0 if memory is eshausted or there is no other locus besides (p+1).
//------------------------------------------------------------------
{
	int q, j, k, n;
	int genei[2], genej[2];
	int fij;
	float f, fmin;
	float sp;			// s-value at locus p, will be return value
	float tolerance;	// difference between average of similarity indices
						// of two samples should be at least this.
	int nLocUsed = nloci;
	char hasdat, skip = 0;
	int jcount;
//	int nfish = 0;
	int maxSibs = NONSIBOUT;
	FISHPTR *nextSamp;
	FISHPTR sampi;
	NONSIBPTR node, prevNode;
	FISHPTR *sampTop, *sampTail;
	sampTop = (FISHPTR*) malloc(sizeof(FISHPTR));
	sampTail = (FISHPTR*) malloc(sizeof(FISHPTR));
	*sampTop = NULL;
	*sampTail = NULL;
	if (p >= LOCOUTPUT) maxSibs = 0; // only print nonsibs if locus <= p

// default values:
	*errcode = 0;
	sp = 0;
	*npairs = 0;
	*ctotal = 0;
	*gotNoSib = 0;
	*jmin = i;	// jmin will be sample pairing with i to be nonsib pair
	tolerance = (float) 0.5/nloci;
	tolerance *= tolerance;	// in effect, tolerance = 1/(4*nloci^2)
	for (q=0; q<nloci; q++) if (*(okLoc+q)==0 || *(nMobil+q)<=1) nLocUsed--;
	if (nLocUsed < 2 || nSamp < 2) return sp;
// create list of genotypes of sample i across all polymorphic loci
	for (q=0; q<nloci; q++) {
		if (*(okLoc+q)==0 || *(nMobil+q)<=1) continue;
		sampi = *(fishList +q);
		for (k=0; k<=i; k++) {
			for (n=0; n<2; n++) genei[n] = sampi->gene[n];
			if (k==i) break;
			sampi = sampi->next;
		};
//		for (k=0; k<i; k++, sampi = sampi->next);
//		for (n=0; n<2; n++) genei[n] = sampi->gene[n];
		// only consider sample i if it has data at locus (p+1):
		if (q==p && genei[0]==0) {
			RemSampLoc (sampTop);
			return sp;
		};
		if ((AddSampRow(sampTop, sampTail, genei) == 0)) {
//		if ((AddSampCount(&nfish, sampTop, sampTail, genei) == 0)) {
			RemSampLoc (sampTop);
			*errcode = 1;
			return sp;
		};
	};

// First, create pointers to sample lists at all loci to travel the lists
// Parameter i is the (i+1)th sample.
// Variable j is the (j+1)th sample in every list to be traveled.
// The j considered must not be i, and must be such that the pair (i,j)
// (unordered) is not one of the putative nonsib pairs chosen before.
// Those pairs are stored in nonsibList. Each node in nonsibList contains
// two fields: field "first" registers the "i" and field "second" registers
// the chosen "j" in previous calls of this function. The list is inserted
// a node representing the newly putative nonsib, whose "first" field is i
// (the parameter) at the end of this function. This function will be called
// for each sample i in the ascending order, so nonsibList was in ascending
// order by the values of fields "first".

// For each j, we search nonsibList to see if the pair (i, j) is already in
// the list. As we search, some in nonsibList are deleted since they cannot
// match at later searches in here or in subsequent calls to this function.
	if ((nextSamp = (FISHPTR*) malloc(sizeof(FISHPTR)*nloci)) != NULL) {
	// position those newly created pointers to point to all sample lists,
	// these are for samples j as in loop "for (j=0; j<nSamp; j++)" below.
		for (q=0; q<nloci; q++) *(nextSamp+q) = *(fishList+q);
		fmin = (float) 5.0;	// represent min of (f values),
		// where f is the value associated with sample pair (i, j), j is
		// an "eligible" sample paired with sample i.
		// For each eligible sample j, we have pair (i,j), We calculate the
		// value f for this pair, which is the average of similarity
		// coancestry of pair (i,j) across loci (q+1) different from (p+1).
		// Since 4 is max value for each similarity index [we temporarily do
		// not divide by 4 as in Nomura's paper, formula (4), since we want
		// function SimilarInd returns integer value], fmin is reassigned
		// after the first f obtained.
		for (j=0; j<nSamp; j++) {
			skip = 0;
		// see if this j is eligible: first, it must not be i
		// Also exclude j that has no data at locus (p+1), sample j at
		// locus (p+1) is pointed to by (nextSamp+p)
			if ((j == i) || ((*(nextSamp+p))->gene[0] == 0)) skip = 1;

		// Now make sure pair (i, j) is not a putative nonsib pair previously
		// chosen from the calls of this function at previous "i". The pairs
		// are stored as nodes in nonsibList.

		// Since the current i here is greater than previous "i", the only
		// case that (i,j) matches to a previously chosen pair is that the
		// j here is one of previous "i" where the chosen "j" (corresponding
		// to that "i") is the current i.
		// Thus, we only need to search for match when j in this "for" loop
		// is less than i, since a match requires j to equal previous "i"
			if (j < i && skip == 0) {
				node = *nonsibList;	// if i = 0, the list is empty
				prevNode = *nonsibList;

		// we look for a node from top down whose field "second" is i.
		// Once such node is found, see if its field "first" is j. If it is,
		// this j is skipped, and this node will be deleted since it won't
		// match the next j and also won't match the i in the next call of
		// this function (the i in the next call is larger than this i).
		// Since we delete such node, then at the next j, if another node
		// has field "second" = i, its "first" field must be at least j,
		// (If it is less, it was matched with previous j and supposed to
		// be deleted!) so if it doesn't match j, it must be bigger, and
		// all nodes down the list should have bigger "first" field!
		// Thus, when a node whose field "second" is i but field "first"
		// is not j, we can quit searching down the list.
				for (k=0; node != NULL; k++) {
					if (node->second == i) // j should != "first" field
					{
						if (node->first == j) {
							skip = 1;
					// this j must be skipped. Also the node is no longer
					// needed for next j, so delete it. Since we will delete
					// such node, then the search for (node->second == i)
					// at (next j) should have (node->first >= next j).
					// At the end (when j reaches i-1), all nodes whose
					// "second" = i are deleted.

// -- FOR CHECKING the ALGORITHM -------------------------------------------
							if (outLoc!=NULL && moreDat==1 && i<maxSibs
								&& detail==1) fprintf(outLoc,
									"%9sRemove from reference: (%d,%d)\n",
									" ", (node->first)+1, (node->second)+1);
// -------------------------------------------------------------------------

							// this resetting of *nonsibTail is crucial,
							// since *nonsibTail is used to add node
							if (*nonsibTail == node) *nonsibTail = prevNode;
							if (k == 0)	// match at the first node
								RemFirstNode(nonsibList, node);
							else
								RemMextNode(prevNode, node);
							(*sibNodes)--;
						};
					// here, (node->first > j) and later nodes with
					// (node->second == i) yield the same thing, so don't
					// need to search more
						break;
					};
					prevNode = node;
					node = node->next;

				};	// end of "for (k=0; node != NULL; k++) "

			};	// end of "if (j < i && skip == 0)"

			if (skip == 0) {	// this includes both cases j<i and j>i.
			// If j>i, skip is still 0 (at this time) if sample j has data at
			// locus (p+1). The condition (skip=0) is for a pair (i, j) to be
			// included in the search of finding a pair that has the smallest
			// value for f = sum of similarity indices for pair (i,j) across
			// all loci differing from locus (p+1), all are assumed to be
			// polymorphic. If not or if a locus is not considered
			// (i.e. (okLoc+q)=0), that locus is skipped as shown in
			// the for loop below
				f = 0;
				jcount = 0;	// number of loci !=(p+1) that (i,j) has data
				sampi = *sampTop;
				// The next "for" loop is to
				//	* sum up all coancestry indices across all loci != (p+1)
				//	  of the pair (i,j), then average them when exit,
				//	* grab coancestry index of the pair (i,j) when j>i, at
				//	  locus (p+1), as part of the sum of coancestry indices
				//	  for all pairs (i,j), j>i at this locus (p+1).
				// A little bit awkward here: every time j changes, sampi
				// must travel the whole list.

// -- FOR CHECKING the ALGORITHM -------------------------------------------
				if (outLoc!=NULL && moreDat==1 && i<maxSibs && detail==1)
					fprintf(outLoc,	"(%2d,%5d)  ", i+1, j+1);
// -------------------------------------------------------------------------

				for (q=0; q<nloci; q++) {
					// samples across loci were built for polymorphic loci
					// that are in consideration:
					if (*(okLoc+q)==0 || *(nMobil+q)<=1) continue;
					for (n=0; n<2; n++) {
					// genotype of sample i at locus (q+1).
						genei[n] = sampi->gene[n];
					// genotype of sample j at locus (q+1).
						genej[n] = (*(nextSamp+q))->gene[n];
					};
					fij = SimilarInd(genei, genej, 0, &hasdat);
				// position sampi at the next eligible locus
					sampi = sampi->next;
					if (q==p && j>i) {	// happens at most once in "q" loop
				// accumulate coancestries for pairs (i,j), j>i, locus (p+1)
				// npairs is the number of pairs having data at locus (p+1)
				// Since sample j and i have data at locus (p+1), npairs is
				// incremented
						*ctotal += fij;
						*npairs += hasdat;
						continue;	// goto next locus, locus (q+1) != (p+1)
					};
					if (q!=p) {

// -- FOR CHECKING the ALGORITHM -------------------------------------------
						if (outLoc!=NULL && moreDat==1 && i<maxSibs
							&& detail==1) fprintf(outLoc,	"%6d", fij);
// -------------------------------------------------------------------------

				// accumulate coancestries at loci different from (p+1)
				// sample j may or may not have data at locus (q+1), so
				// jcount may not increase
						f += fij;
						jcount += hasdat;
					};
				};	// end of "for (q=0; q<nloci; q++)"

// -- FOR CHECKING the ALGORITHM -------------------------------------------
				if (outLoc!=NULL && moreDat==1 && i<maxSibs && detail==1)
					fprintf(outLoc,	"\n");
// -------------------------------------------------------------------------

				// compare f and fmin to reassign jmin:
				// using "<" is to take the first j that yields min value,
				// using "<=" is to take the last j that yields min value.
				// (f is valid only when pair (i,j) has data at some locus)
				if (jcount>0) {
					f /= (float) jcount;
				// Compare f just calculated and fmin in storage, use
				// tolerance as a guard against rounding error.
				// The use of tolerance is to prevent the case that f = fmin
				// but the conversion or rounding error causes f < fmin.
				// The value of tolerance should be large enough so that
				// it is larger than the rounding error, and small enough
				// so that it does not affect the comparison between f, fmin.
				//
				// Explanation:
				// -----------
				// Suppose f, fmin are exact, no rounding error, and suppose
				// fr, fminr are values of f, fmin stored in the machine.
				// Let r be the maximum rounding error. Then
				// |f - fr| <= r, |fmin - fminr| <=r.
				// Suppose the difference between f values is more than 2c.
				// (Here, we have f calculated by taking the sum of similarity
				// indices as integers across nloci-1 at most, then take the
				// average, so the difference between any two is less than
				// 1/n^2, where n = nloci; so we can take 2c = 1/n^2 or less.)
				// Therefore, if f < fmin, then f + c < fmin - c; and since
				// the rounding error r should be far less than c/4, we have
				// fr + c/2 <= f+2r+ c/2 < f+c < fmin-c < fminr+r-c < fminr.
				// The value of tolerance here is at most c/2, and when we
				// add f+tolerance, we have rounding error r in f and then
				// another r when adding; that's why we use f+2r+ c/2.
					if (f+tolerance < fmin) {
						*jmin = j;
						fmin = f;
					};
				};	// end of "if (jcount>0)"
			};	// end of "if (skip == 0)"
			// advance sample, to match index j
			for (q=0; q<nloci; q++) {
				if (*(okLoc+q)==0 || *(nMobil+q)<=1) continue;
				*(nextSamp+q) = (*(nextSamp+q))->next;
			};
// this caused err:
//			for (q=0; q<nloci; q++) *(nextSamp+q) = (*(nextSamp+q))->next;
		};	// end of "for (j=0; j<nSamp; j++)", we went thru all eligible j,

// -- FOR CHECKING the ALGORITHM -------------------------------------------
		if (outLoc!=NULL && moreDat==1 && i!=*jmin  && i<maxSibs
			&& detail==1)
		{
			fprintf(outLoc, "%9sChosen pair", " ");
			if (*jmin > i) fprintf(outLoc, ", added to reference\n");
		};
// The pair will be inserted in the function that calls this function
// -------------------------------------------------------------------------

		// After going thru all eligible samples j != i:
		// * In general, got putative nonsib pair (i, jmin) with minimum f.
		//   At the same time, accumulated coan. "ctotal" for pairs (i,j),
		//   j>i, at locus (p+1)
		// * if no eligible j, jmin=i. This could happen if i=nSamp-1 in the
		//   case all samples i<nSamp-1 paired with (nSamp-1) as nonsibs.
		// * if there are eligible j, but the pairs (i,j) have missing data
		//   at every locus besides (p+1), they are not counted, so jmin=i.
		// There will be no new putative nonsib pair when jmin=i.
		// * npairs = 0 when no j > i such that pair (i,j) has data. This
		//   could happen even if jmin != i, i.e., a nonsib pair is found.
		//   (trivial case: when i = nSamp-1, npair = 0, but we likely can
		//   find a jmin < i)

		// Note : SimiarInd gives values as integer, should be divided by 4
		// to get the coancestry index
		*ctotal = *ctotal/4;
		if (*jmin == i) {	// no nonsib pair found.
		// only ctotal may be collected - npairs=0 or not.
			*ctotal = *ctotal/4;
			RemSampLoc (sampTop);
			free(nextSamp);
			return sp;	// which is 0
		};
		*gotNoSib = 1;	// signaling that a nonsib pair is found.
		// point to sample jmin at locus p+1, then evaluate similarity index
		// for the pair (i, jmin) at that locus for the return values.
		// if jmin > i, this value was actually evaluated before, at the time
		// the jmin was not known to pair with i as a nonsib pair.
		*(nextSamp+p) = *(fishList+p);
		for (j=0; j<nSamp; j++) {
			if (j== *jmin)
				for (n=0; n<2; n++) genej[n] = (*(nextSamp+p))->gene[n];
			if (j==i)
				for (n=0; n<2; n++) genei[n] = (*(nextSamp+p))->gene[n];
			if (j>= *jmin && j>=i) break;
			*(nextSamp+p) = (*(nextSamp+p))->next;
		};
/*
		for (j=0; j<jmin; j++) *(nextSamp+p) = (*(nextSamp+p))->next;
		sampi = *(fishList+p);
		for (k=0; k<i; k++) sampi = sampi->next;
		for (n=0; n<2; n++) {
			genej[n] = (*(nextSamp+q))->gene[n];
			genei[n] = sampi->gene[n];
		};
*/
		fij = SimilarInd(genei, genej, 0, &hasdat);
		sp = (float) fij / (float) 4.0;
		// only add putative nonsib node to nonsibList if jmin > i:
		if (*jmin > i) {
			if ((node = MakeNonsib (i, *jmin)) == NULL) {
				*errcode = 1;
				free(nextSamp);
				RemSampLoc (sampTop);
				*sampTail = NULL;
				free (sampTail);
				return sp;
			};
			if (*nonsibList == NULL) {
				*nonsibList = node;
			} else {
				(*nonsibTail)->next = node;
			};
			*nonsibTail = node;
			(*sibNodes)++;
		};
		free (nextSamp);
		*sampTail = NULL;
		RemSampLoc (sampTop);
		free (sampTail);
		return sp;
	} else {	// memory cannot be allocated for nloci pointers
		*sampTail = NULL;
		free (sampTail);
		RemSampLoc (sampTop);
		*errcode = 1;
		return sp;	// still 0
	};
}



//------------------------------------------------------------------

float CoanDiff (FISHPTR *fishList, int nMobil[], int p, int nloci,
			int nSamp, char *okLoc, float *sp, FILE *outLoc, char moreDat,
// add missptr, hSamp in Mar 2012 to calculate harmonic mean of sample size
			int count, float *hSamp, int *polyLoc)

// Estimate the average molecular coancestry s^(p) for locus (p+1) over nSamp
// pairs of putative nonsibs (i, j), where i = 0, 1, ..., nSamp-1.
// The pairs are to be distinct, i.e. no two pairs are identical.
// (Pairs (1,3) and (3,1) are identical, but (1,3) and (8,3) are distinct.)
// Those pairs are created in function PutativeNonSib, which evaluates
// their coancestry values at locus (p+1).
//
// Calculate the average molecular coancestry "sp" of all pairs at locus
// (p+1) based on putative nonsib pairs.
// The total molecular "totcoan" of all pairs (i,j), j>i, is averaged from
// ctotal(i), which is calculated from PutativeNonSib as the sum of all
// coancestries  of pairs (i,j), i fixed, j>i.
// The function then returns the difference between the two..
{
	int i, jmin;
	NONSIBPTR *nonsibList, *nonsibTail;
	NONSIBPTR node;
	float ctotal, totcoan;
	float totpairs;
	char gotNoSib, errcode;
	int npairs;
	int nputSibs;
	int size;
	int sibNodes = 0;
	int maxSibs = NONSIBOUT;

// To turn off detailed how putative nonsib pairs are chosen to outLoc file,
// set *detail = 0. To turn on, set it = 1:
	char detail = 0;
	if (p >= LOCOUTPUT) maxSibs = 0;	// only print nonsibs up to LOCOUTPUT loci
	// no matter what detai is set, only affect if maxSibs > 0
	if (maxSibs <= 0) detail = 0;

	size = sizeof(NONSIBPTR);
	nonsibList = (NONSIBPTR*) malloc(size);
	nonsibTail = (NONSIBPTR*) malloc(size);
	*nonsibList = NULL; *nonsibTail = NULL;
// default values:
	*sp = 0;
	totcoan = 0;
	totpairs = 0;
	nputSibs = 0;
	if (*(okLoc+p)==0 || *(nMobil+p)<=1) return 0;
// add in Mar 2012 to calculate harmonic mean of sample size:
	(*polyLoc)++;
	if (count > 0) *hSamp += 1.0F/((float) count);
	if (outLoc != NULL && moreDat == 1 && maxSibs > 0) {
		if (detail == 1) fprintf (outLoc,
			"Molecular coancestries of pairs (i,j) listed below should be "
			"divided by 4\n\n");
		fprintf (outLoc, "Locus %d:", p+1);
		if (detail == 1) fprintf (outLoc,
			"\nSample Pair    4*Molecular coancestry at other polymorphic loci\n");
	};
	for (i=0; i<nSamp; i++)
	{
		*sp += PutativeNonSib (nonsibList, nonsibTail, i, &jmin, p, &npairs,
							&ctotal, fishList, nMobil, nloci, okLoc,
							nSamp, &gotNoSib, &sibNodes, &errcode,
							outLoc, moreDat, detail);
		totcoan += ctotal;
		totpairs += (float) npairs;
		nputSibs += gotNoSib;
		if (outLoc != NULL && moreDat == 1 && i < maxSibs)
			fprintf (outLoc, "  (%d,%d)", i+1, jmin+1);
// this "if" can be removed (used to check the list of referenced nonsibs):
		if (outLoc != NULL && moreDat == 1 && i < maxSibs && detail==1) {
			fprintf (outLoc, ". Reference list =");
			for (node = *nonsibList; node != NULL; node=node->next)
				fprintf (outLoc, " (%d,%d)", node->first+1,node->second+1);
			fprintf (outLoc, "\n");
		};

	};
	if (nputSibs > 0) *sp /= (float) nputSibs;
	if (totpairs > 0) totcoan /= (float) totpairs;
	if (outLoc != NULL && moreDat == 1 && maxSibs > 0) {
		fprintf (outLoc, "\n\n  [n0 = %d,    s^ = %12.8f,    fm = %12.8f]\n\n",
				nputSibs, *sp, totcoan);
		fflush (outLoc);
	};
	*nonsibTail = NULL;
	free (nonsibTail);
// actually, nonsibList is already empty since its nodes are deleted
// at the calls to function PutativeNonSib, but we still need to remove
// memory for the pointer to nonsibList, just as doing with nonsibTail.
	RemNonSib (nonsibList);
	return (totcoan - *sp);

}

//------------------------------------------------------------------
void PrintCoan (FILE *outLoc, char moreDat, COANPTR *coanList)
{
// print molecular coancestry and related value at each locus
	COANPTR node;
	float fm, coan, f2, wt;
	int p;
	if (outLoc == NULL || moreDat == 0) return;
	fprintf (outLoc,
		"  Locus   sum(freq^2)      s^           fm         weight\n");
	for (node = *coanList; (node != NULL); node =node->next)
	{
		wt = node->weight;
		p = (node->locus) + 1;
		coan = (node->scoan);
		f2 = node->fresq;
		fm = (node->diffcoan) + (node->scoan);
		if (wt == 0) fprintf (outLoc, "%6d*", p);
		else fprintf (outLoc, "%6d ", p);
		fprintf (outLoc, "%12.6f %12.6f %12.6f %12.6f\n", f2, coan, fm, wt);
		fflush (outLoc);
	};
}




//------------------------------------------------------------------


//------------------------------------------------------------------
float NbCoan (float f1)
{
	float Nb, f2;
	f2 = 2*f1;
	if (fabs(f2) < EPSILON) Nb = INFINITE;
	else Nb = (float) 1.0/f2;
// add this to make sure no negative:
	if (Nb < 0) Nb = INFINITE;
	return Nb;
}



//------------------------------------------------------------------


//------------------------------------------------------------------
// CoanList consists of all information to do Jackknife on loci:
// locus: p
// scoan: sp, ave. of coan coefs from putative nonsib pairs
// freqsq: sum of squares of frequencies in locus
// diffcoan : (Ave. of coan coef. of all sample pairs) - sp
// weight: wp, weight of locus p
void ConfidCoan(COANPTR *coanList, float *loNe, float *hiNe)
// The first part of this function is to calculate the adjusted variance
// with Jackknife ethod; this part is parallel to the one for temporal,

// After the first part, the second part is to call JackKnifeInd to get
// degree of freedom, then use Chi-square approximation to get
// CI for f1
{
	int n;
	float fval, fpw, wp, sp, diff, fBar, totW, varF;
	long iBig;
	float loLim, hiLim;
	COANPTR node;
	*loNe = INFINITE, *hiNe = INFINITE;
// The first for loop is to find the weighted average.
// Actually, in the calling function, the function PutCoanInd will be called
// first, and in that function, both the overall f and the weight are
// calculated (total weight is totW here); so we may skip this for loop
// if we insert more parameters in the list.
// However, the calculation is not much long, keep this for loop
// as a double check if necessary. Note that this "for" loop will count
// the number n of loci eligible for jackknife. If this loop is deleted,
// then the "n" should be in the next "for" loop (with incremented)
	for (n=0,fBar=0,totW=0, node = *coanList; node != NULL; node = node->next)
	{
		wp = node->weight;
		if (wp <= 0) continue;	// then 1-sp=0
		n++;
		totW += wp;
		sp = node->scoan;
		diff = node->diffcoan;
		fval = diff/(1-sp);
		fBar += fval*wp;
	};
	// should have at least 2 loci having weights
	if (n <= 1) return;
	fBar = fBar/totW;
//printf ("fBar = %10.7f\n",fBar);
	varF = 0;
	for (node = *coanList; node != NULL; node = node->next)
	{
		wp = node->weight;
		if (wp <= 0) continue;	// then 1-sp=0
		sp = node->scoan;
		diff = node->diffcoan;
		fval = diff/(1-sp);
		fpw = (wp/(totW-wp))*(fBar-fval);
		varF += fpw*fpw;
	};
    varF = ((n - 1)*varF) / ((float) n);
	iBig = JackKnifeInd(fBar, varF);
	Confid95 (iBig, fBar, &loLim, &hiLim);
//printf ("iBig = %10d, fBar = %12.4f, varF = %12.4e, loLim = %12.4e, hiLim = %12.4e\n",
//iBig, fBar, varF, loLim, hiLim);

	*loNe = NbCoan (hiLim);
	*hiNe = NbCoan (loLim);
}
//------------------------------------------------------------------
// CoanList consists of all information to do Jackknife on loci:
// locus: p
// scoan: sp, ave. of coan coefs from putative nonsib pairs
// freqsq: sum of squares of frequencies in locus
// diffcoan : (Ave. of coan coef. of all sample pairs) - sp
// weight: wp, weight of locus p

void CoanConfid (COANPTR *coanList, float *loNe, float *hiNe)
// The difference between this and the above function is: this does not
// weight the Jackknife blocks. (Suppose there are N items, Pick
// (N-1) items to find weighted average. This does not put a weight to
// this average value.)

// The first part of this function is to calculate the adjusted variance
// with Jackknife method; this part is parallel to the one for temporal,

// After the first part, the second part is to call JackKnifeInd to get
// degree of freedom, then use Chi-square approximation to get
// CI for f1
{
	int n;
	float fval, wp, sp, diff, fBar, totW, totf, b, nf, fvalj, fAve, varF;
	long iBig;
	float loLim, hiLim;
	COANPTR node;
	*loNe = INFINITE, *hiNe = INFINITE;
// get total weights and total weighted fs
	for (n=0,totf=0,totW=0, node=*coanList; node != NULL; node = node->next)
	{
		wp = node->weight;
		if (wp <= 0) continue;	// then 1-sp=0
		n++;
		totW += wp;
		sp = node->scoan;
		diff = node->diffcoan;
		fval = diff/(1-sp);
		totf += fval*wp;
	};
	// should have at least 2 loci having weights
	if (n <= 1) return;
	fAve = totf/totW;
	for (n=0, fBar=0, varF=0, node=*coanList; node != NULL; node = node->next)
	{
		wp = node->weight;
		if (wp <= 0) continue;	// then 1-sp=0
		n++;
		sp = node->scoan;		// these 4 lines from here are for f with weight
		diff = node->diffcoan;
		fval = diff/(1-sp);
		fval *= wp;				// fval now is f with weight at this node
		fvalj = (totf - fval)/(totW - wp);	// weighted ave. excluding the node
		nf = 1.0F/((float) n);
		b = fvalj - fBar;
		fBar += nf*b;
		if (n > 1)
//				varF += nf*(n-1)*(b*b);
// the following formula is the same, it uses updated fBar
				varF += b*(fvalj - fBar);
	};
	nf *= ((float) n - 1);
	varF *= nf;
	iBig = JackKnifeInd(fBar, varF);
//	Confid95 (iBig, fBar, &loLim, &hiLim);
// Choose fAve, the weighted average when all nodes are used. It may differ
// from fBar above.
	Confid95 (iBig, fAve, &loLim, &hiLim);
	*loNe = NbCoan (hiLim);
	*hiNe = NbCoan (loLim);
}

//------------------------------------------------------------------




//------------------------------------------------------------------

float WeightAtLoc0 (int p, ALLEPTR *alleList, float sp, float *freq2)
{
// Suppose there are *(nMobil+p) alleles at locus (p+1), say
// a(1), ... a(m). Let
//		r(p) = sum {j=1,...,m}(a(j)^2),
//		f(p) = r(p) * (1 - r(p)).
// This f(p) and sp, the average coancestry indices of putative nonsib pairs,
// are used to evaluate the weight (cf. 2nd column of page 464 in Nomura's).
// The weight w(p) at locus p is given by:
//	   w(p) = v(p)*(1-sp), where v(p) = (1-sp) / f(p).

	float r, fp, vp;
	ALLEPTR curr;
	r = 0;
	for (curr = *(alleList+p); (curr != NULL); curr =curr->next)
		r += (curr->freq) * (curr->freq);
	*freq2 = r;
	fp = r*(1-r);
	// reassign r as 1-sp:
	r = (float) 1.0 - sp;
	if (fabs(r) < EPSILON) r = 0;
// assign weght for monomorphic locus:
//	if (fp == 0) vp = 1;
	if (fp == 0) vp = 0;
	else vp = r/fp;
	return vp;
}
//------------------------------------------------------------------

int PutCoanInd0 (FISHPTR *fishList, ALLEPTR *alleList, int nMobil[],
				int nloci, int nSamp, char *okLoc, COANPTR *coanList,
				float *f1, FILE *outLoc, char moreDat,
// add missptr, hSamp in Mar 2012 to calculate harmonic mean of sample size
				int *missptr, float *hSamp)
// Estimate the average molecular coancestry per locus, for all loci,
// and put the values in coanList, which is created as going thru all loci.
// Return -1 if the list cannot be completed, 0 if OK. However, the list
// is considered as extra, not crucial.
// The main purpose is to evaluate f1, used to find effective breeder.
// f1 = (1/W)*[sum aross loci of (vp*(fp - sp))], where
// W = sum of (weight wp at locus (p+1)),
// vp is calculated in WeightAtLoc,
// sp is calculated from CoanDiff, and (fp-sp) is the return value of
// that function.
{
	int p, count, polyLoc;
	float fdiff, sp, vp, wp, freq2;
// add time to inform user:
//	time_t rawtime;
//	char info [15], *timeptr;
	long prompt, scount, locprt;
//	int i=0;

	float totW = 0;
	COANPTR node, *coanTail;
	int errcode = 0;
	coanTail = (COANPTR*) malloc(sizeof(COANPTR));
	*coanTail = NULL;
	*f1 = 0;
	*hSamp = 0;
	polyLoc = 0;
// for informing the progress on console:
	prompt = 1000000;//000;
	scount = nloci*nSamp;
	locprt = 0;
//	info[15] = '\0';
	for (p=0; p<nloci; p++) {
		if (*(okLoc+p)==0 || *(nMobil+p)==0) continue;
		count = nSamp - (*(missptr+p));
// put out information to the console, so the user can see how much progress:

		if (locprt >= prompt) {
			printf ("       Coan-Coeffs done up to locus %d\n", p+1);
//			time ( &rawtime );
//			timeptr = ctime (&rawtime);
//			for (i = 0; i < 15; info[i] = *(timeptr+4+i), i++);
//			printf (", %s\n", info);
			locprt = 0;
		};

// estimate average molecular coancestry at locus (p+1):
		fdiff = CoanDiff (fishList, nMobil, p, nloci, nSamp, okLoc, &sp,
						outLoc, moreDat, count, hSamp, &polyLoc);
		vp = WeightAtLoc0 (p, alleList, sp, &freq2);
		wp = ((float)1.0-sp)*vp;
		totW += wp;
		*f1 += vp*fdiff;
		if ((node = MakeCoan (p, sp, fdiff, wp, freq2)) == NULL) {
			errcode += 1;
		} else {
		// add to coanList this new node
			if (*coanList == NULL) *coanList = node;
			else (*coanTail)->next = node;
			*coanTail = node;
		};
// for putting out information to the console:
		locprt += scount;
	};
	if (*hSamp > 0) *hSamp = ((float) polyLoc)/(*hSamp);
	if (outLoc != NULL && moreDat == 1) fprintf (outLoc, "\n");
	*f1 /= totW;
	return errcode;	// number of times that a node cannot be added.

}
//------------------------------------------------------------------

float CoanMethod (FISHPTR *fishList, ALLEPTR *alleList, int nMobil[],
				int nloci, int nSamp, char *okLoc, float *f1, FILE *outLoc,
				char moreDat, float *loNbCoan, float *hiNbCoan, char jack,
				int *missptr, float *hSamp)
{
	COANPTR *coanList;
	float Nb;
	coanList = (COANPTR*) malloc(sizeof(COANPTR));
	*coanList = NULL;
	if (outLoc != NULL && moreDat == 1) {
		PrtLines (outLoc, 15, '-');
		fprintf (outLoc, "\nNOMURA's METHOD\n");
		PrtLines (outLoc, 15, '-');
		fprintf (outLoc, "At each locus,\n");
		fprintf (outLoc,
			"* s^ = molecular coancestry average of putative nonsib pairs\n");
		fprintf (outLoc,
			"* fm = molecular coancestry average of all pairs (i,j), i < j\n");
		if (NONSIBOUT > 0) {
			fprintf (outLoc, "* n0 = total number of putative nonsib pairs\n\n");
			fprintf (outLoc,
				"Putative nonsib pairs (i,j) are displayed below");
			if (nSamp > NONSIBOUT) fprintf (outLoc, ", up to i = %d", NONSIBOUT);
			fprintf (outLoc, "\n\n");
		};
	};
	printf ("     Molecular Coancestry Method\n");	// next code may run slow!
	if (PutCoanInd0 (fishList, alleList, nMobil, nloci, nSamp, okLoc,
					coanList, f1, outLoc, moreDat, missptr, hSamp) == 0)
	{
		if (jack == 1) CoanConfid(coanList, loNbCoan, hiNbCoan);
		if (outLoc != NULL && moreDat == 1)
		{
			PrtLines (outLoc, 58, '-');
			PrintCoan (outLoc, moreDat, coanList);
			PrtLines (outLoc, 58, '-');
			fprintf (outLoc, "\n");
		};
		RemCoan (coanList);
	};
//	f2 = 2*(*f1);
//	if (fabs(f2) < EPSILON) Nb = INFINITE;
//	else Nb = (float) 1.0/f2;
	Nb = NbCoan (*f1);
	if (outLoc != NULL && moreDat == 1) {
		fprintf(outLoc, "   f1^ = %10.6f    ", *f1);
		if (Nb < INFINITE) fprintf(outLoc, "Neb^ = 1/[2(f1^)] =%12.2f\n", Nb);
		else fprintf(outLoc, "Neb^ = INFINITE\n");
		PrtLines (outLoc, 58, '-');
		fprintf (outLoc, "\n\n");
	};
	printf ("       Estimated Neb^: ");
	if (Nb < INFINITE) printf ("%20.1f\n", Nb);
	else printf ("%20s\n", "Infinite");
	return Nb;
}


//------------------------------------------------------------------

// Temporal module
//------------------------------------------------------------------

FREQPTR MakeFreq (int mvalue, int samp, float freq,
				  int nGeneration, int generation)
// This create a new Freq node, return NULL if failed,
// The new node is filled by data from the list of alleles across loci
// that were stored when reading input, whose frequencies were also
// calculated. Those data are parameters samp and freq.
// The node is assumed to be created at generation "generation."
// (Each block of populations consists of nGeneration populations:
// Generation 0 up to Generation "nGeneration-1")
// A NULL value indicates that memory runs out.
{
	int i, size;
	FREQPTR newptr;
	size = sizeof(struct timefreq);
	if ((newptr = (FREQPTR) malloc(size)) != NULL)
	{
		newptr->mValue = mvalue;
		if ((newptr->samples = (int*) malloc(sizeof(int)*nGeneration))
				== NULL) return NULL;
		if ((newptr->freqs = (float*) malloc(sizeof(float)*nGeneration))
				== NULL) return NULL;
		for (i=0; i<nGeneration; i++){
			*(newptr->samples +i) = 0;
			*(newptr->freqs +i) = 0;
		};
		*(newptr->freqs +generation) = freq;
		*(newptr->samples +generation) = samp;
		newptr->next = NULL;
	}
	return newptr;
};

//------------------------------------------------------------------

FREQPTR AddFreq (FREQPTR head, int alleleK, int samp, float freq,
				int nGeneration, int generation, int *errcode)
{
// This routine adds one allele to the list of alleles designed for
// temporal method (at one locus). Once it is done for one generation,
// then for subsequence generations, it is rarely added, except in the case
// there is an allele in the next generation that is not in previous one.

// Return the pointer to the head of the list,
// The list of alleles is in ascending order by its mValue's.
// Note that the list of alleles is supposed to be known

    FREQPTR ptr1, ptr2, newnode;

	*errcode = 0;
	if (head == NULL) {
		if ((newnode = MakeFreq (alleleK, samp, freq, nGeneration, generation))
				== NULL)
		{;
			printf  ("Out of memory for storing allele in temporal method!\n");
			*errcode = -1;
			return head;	// no change
		};
		head = newnode;
	} else {
		if (alleleK < head->mValue) {	// add to front
			if ((newnode = MakeFreq (alleleK, samp, freq, nGeneration, generation))
				== NULL)
			{;
				printf ("Out of memory for storing allele in temporal method!\n");
				*errcode = -1;
				return head;	// no change
			};
			newnode->next = head;
			head = newnode;
		} else {
			ptr1 = head;
			ptr2 = ptr1->next;
// search for the last node ptr1 such that: ptr1->mValue <= alleleK
			while (ptr2 != NULL) {
				if (ptr2->mValue <= alleleK) {
					ptr1 = ptr2;	// then ptr1->mValue <= alleleK
					ptr2 = ptr2->next;
				} else {
					break;
				}
			}	// end of search

// now, ptr1->mValue <= alleleK, either ptr1 is at the end of the list,
// or alleleK is between mValues of ptr1 and ptr2.
//
			if (ptr1->mValue == alleleK) {
			// this is common when this function is called in the next
			// generation: alleleK is already appeared in previous generation,
			// so only add its frequency and samples having this allele.
			// This generation is assumed to be different from the generation
			// that this allele was added to the list.
				*(ptr1->freqs +generation) = freq;
				*(ptr1->samples +generation) = samp;
			} else {
				if ((newnode = MakeFreq (alleleK, samp, freq, nGeneration, generation))
					== NULL) {;
					printf ("Out of memory for storing allele in temporal method!\n");
					*errcode = -1;
					return head;	// no change
				};
				newnode->next = ptr2;
				ptr1->next = newnode;
			}
		}	// end of "if (alleleK < head->mValue) ... else "
	}
	return head;
}

//------------------------------------------------------------------


void AddFreqWide (FREQPTR *freqList, ALLEPTR *alleList, int nloci,
				  int nfish, int *missptr, char *locUse, int nGeneration,
				  int generation, int *errcode, char weighsmp)
// Create list of frequencies from AlleList for one generation.
// *(alleList+p) is the list of alleles at locus (p+1)
{

	ALLEPTR ptr1;
	FREQPTR head;
	int alleleK;
	int p, count;
	float freq;
	*errcode = 0;
	for (p=0; p<nloci; p++) {
		// "count" is the number of samples having data at locus (p+1).
		// Note that this will come in every node of the list at locus (p+1),
		// so it is quite a redundant and waste of storage (but it is
		// convenient, otherwise, we will need to create a separate list)
		// Note that when a node is added, the field samples are filled
		// with zeros for all generations except at the current one,
		// so it will need to be adjusted when the list is complete.
		count = nfish;
		if (weighsmp > 0) count -= (*(missptr+p));
		for (ptr1 = *(alleList+p); ptr1 != NULL; ptr1 = ptr1->next)
		{
			alleleK = ptr1->mValue;
			freq = ptr1->freq;
			head = *(freqList+p);
			*(freqList+p) = AddFreq (head, alleleK, count, freq,
								nGeneration, generation, errcode);
			if (*errcode != 0) return;
		};
	};
}
//------------------------------------------------------------------



void FreqAdjnPrt (FILE *output, char moreDat, FREQPTR *freqList, int nloci,
				  int nGeneration, char *locUse)
{
	int p, i, j, samp;
	float freq;
	FREQPTR ptr1;
	// Adjust number of samples (as mentioned in function AddFreqWide),
	for (i=0; i<nGeneration; i++)
	{
		for (p=0; p<nloci; p++)
		{
			samp = 0;
			for (ptr1=*(freqList+p); ptr1!=NULL; ptr1=ptr1->next) {
				samp = *(ptr1->samples +i);
				if (samp > 0) break;
			};
			for (ptr1=*(freqList+p); ptr1!=NULL; ptr1=ptr1->next)
				*(ptr1->samples +i) = samp;
		};
	};

	if (output == NULL || moreDat == 0) return;
	fprintf (output, "\nTEMPORAL METHOD: From %d samples\n", nGeneration);
	for (p=0; p<31; p++) fprintf (output, "=");
	fprintf (output, "\n");
	fprintf (output, "%18s\n", "Having");
//	fprintf (output, "%17s\n", "data");
	for (p=0; p<nloci; p++) {
		if (*(locUse+p)==0) continue;
		if (p == 0)
			fprintf (output, "Locus %d      Data\n", p+1);
		else
			fprintf (output, "Locus %d\n", p+1);
		fprintf (output, "%-15s", "Alleles:");
//		fprintf (output, "Locus %d\n%-15s", p+1, "Alleles:");
	// go over all allele mobility values for locus (p+1):
		samp = 0;
		for (ptr1 = *(freqList+p); ptr1 != NULL; ptr1 = ptr1->next)
			fprintf (output, "%8d", ptr1->mValue);
		fprintf (output, "\n");
		for (i=0; i<nGeneration; i++)
		{
			fprintf (output, "Sample %3d:", i+1);
			for (j=0, ptr1 = *(freqList+p); ptr1!=NULL; ptr1=ptr1->next, j++)
			{
				samp = *(ptr1->samples +i);
				freq = *(ptr1->freqs +i);
				if (j==0) fprintf (output, "%6d", samp);
				fprintf (output, "%8.4f", freq);
			};
			fprintf (output, "\n");
		};
		fprintf (output, "\n");
	};
}


//------------------------------------------------------------------

void RemoveFreq(FREQPTR *freqList, int nloci)
{

	int p;
	FREQPTR curr, temp;
	for (p=0; p<nloci; p++) {
		curr = *(freqList+p);
		for (; curr != NULL; curr = temp) {
			temp = curr->next;
			free ((curr->samples));
			free ((curr->freqs));
			free (curr);
		};
		*(freqList+p) = NULL;
    }

}
//------------------------------------------------------------------

float FprimeS (float fs, float invhmonic, float invcount2, int census)
// while invhmonic is the inverse of harmonic mean of two samples at a locus,
// invcount2 is the inverse of 2*(size of sample 2).
{
	float extra, x, y;
// fsprime is given by (13) p. 931 in Jorde/Ryman paper
//	x = fs*(1.0F - invhmonic*4) - invhmonic;
//	y = (1.0F + fs/4)*(1.0F - invcount2);
	extra = 0;
	// extra term is for Plan I (census size is given > 0)
	if (census > 0) extra = 1.0F/((float) census);
	x = fs*(1.0F - invhmonic/4 + extra/4) - invhmonic + extra;
	y = (1.0F + fs/4)*(1.0F - invcount2);
	return x/y;

}

//------------------------------------------------------------------


float NeFromFp (float fprime, float timegap)
// Calculate Nb = Eff number of breeders, for Jorde/Ryman method
// of sample size across loci and generations at time t1=time1, t2=time2,
// fprime-value is calculated from the method
// Nb = |t1 - t2|/(2*fprime), where t1=time1, t2=time2.
{
	if (fabs(fprime) <= EPSILON) return INFINITE;
	else return timegap/(2*fprime);
}

//------------------------------------------------------------------

float FvalRev (float fpIni, float *fprime, float *hSampMean, float *nAlle,
				 float aux, int nloci, float timegap, int excluded)
// from arrays fprime, hSampMean (harmonic mean of sample sizes), nAlle
// (containing number of independent alleles) at loci, and a constant aux
// (tentative Ne calculated based on all loci),
// we recalculate fprimeOverAll as weighted mean with new weights. The new
// weight at a locus is (nAlle*hSampMean^2)/(hSampMean*timegap + aux)^2
// The parameter "excluded" is the locus to be excluded from calculation:
// this is useful when doing jackknife.
{
	float a;
	int p;
	float fp, fpMean, hm, topw, bottomw, w, totW;
// Algorithm to find the weighted mean of f(1), ..., f(n) with weights
// w(1), ..., w(n): (here, n is the number of loci that are considered)
// Start with W0 = 0, F0 = 0,
// W(i) = W(i-1) + w(i)
// F(i) = F(i-1) +(w(i)/W(i))*(f(i)-F(i-1)).
// Then F(i) will be the weighted mean of all f(j), j <= i;
// in particular, F(n) is the weighted mean of all f(j).
// It should be noted that W0 must be initialized = 0, but F0 can be anything
// In fact, at i=1,
//		W1 = w1, F1 = F0 + (w1/W1)(f1 - F0) = F0 - (f1 - F0) = f1,
// so initial value F0 doesn't matter.

	a = 2*aux;
	if (a >= sqrt(INFINITE) || a <= 0) return fpIni;
	totW = 0;
	fpMean = 0;
	for (p=0; p<nloci; p++) {
		if (p == excluded) continue;
		hm = *(hSampMean+p);	// harmonic mean of sample size, locus (p+1)
		fp = *(fprime+p);		// fprime at locus (p+1)
		if (hm <= 0) continue;	// this locus is not considered
		// calculate weight w at locus (p+1):
		topw = ((float) *(nAlle+p))*hm*hm;	// (ind. alleles)*(samp size)^2
		bottomw = (hm*timegap)+a;
		bottomw *= bottomw;
		w = topw/bottomw;
		totW += w;
		fp -= fpMean;
		fpMean += (w/totW)*fp;
	};
	return fpMean;	// the revised fprime mean
}


//------------------------------------------------------------------

float VarTemJR (float totFsTop, float *ftop,
					float totFsBottom, float *fbottom,
					float invhSum, float *invhmonic,
					float inv2Sum, float *inv2count,
					int nloci, float *fprimeBar, int census)

// This is for Jorde/Ryan method.
// totFsTop: sum of all numerators in the ratio representing Fs - all loci,
// ftop: array of numerators for Fs at each locus.
// Similarly for topFsBottom, fbottom.
// invhSum: sum of inverses of harmonic sizes of two samples across loci,
// invhmonic: array of harmonic sample size for each locus.
// Similarly, inv2Sum: sum across loci of inverses of sizes of sample 2,
// inv2count: array of inverses of sizes of sizes of sample 2 at loci.

// Algorithm to find variance S of the mean:
// Given x1, ..., xn, the n sample values obtained; in this case, fsprime(p)
// from all loci except locus p, for jackknife process.
// Let M(p) be the mean of fsprime(j) over all j=1,...,p; so M = M(n) will be
// the mean of all fsprime(p) (p=1,...,n).
// Let S(p) be the sum of (fsprime(j) - M(p))^2, over all j=1,...,p; so
// S = S(n) = sum of (fsprime(j) - M)^2, then SE = sqrt{((n-1)/n)* S}.
// We can obtain M(p), S(p) in one loop (p=1,..., n) as follows.
// Starting with M(0) = 0, S(1)= 0.
// (Note that not all loci are eligible, so p will be the pth eligible locus,
// not necessarily locus p.)
// M(p) = M(p-1) + [fsprime(p) - M(p-1)]/p
// S(p) = S(p-1) + [(p-1)/p]*[fsprime(p) - M(p-1)]^2
//      = S(p-1) + [fsprime(p) - M(p-1)]*[fsprime(p) - M(p)]
// Some algebraic manupilations are needed to see that M(p), S(p)
// can be derived this way.
// In the code, fBar stands for M(p), and varF for S(p).
//
// Suppose the n sample values  x1, ..., xn are weighted w1, ..., wn.
// The induction process starts with the weight W(0) = 0, M(0) = 0, S(1) = 0,
// W(p) = W(p-1) + wp
// M(p) = M(p-1) + (wp/W(p-1))*(xp - M(p-1))
// S(p) = S(p-1) + [(wp*W(p-1)) / W(p)]*(xp - M(p-1))^2.
// Then, at each step,
//      W(p) = w1 + ... + wp,
//      M(p) = (1/W(p))*(w1x1 + ... + wpxp),
//      S(p) = w1(x1 - M(1))^2 + ... + wp(xp - M(p))^2.
//
// If wp = 1 for all p, we come back to unweighted case.
// At final step, M(n) is the weighted average, S(n) = sum of (x(j) - M(n))^2
// The weighted version is not implemented here.
{
//	float loLim, hiLim, stdErr;
//	float infinite = INFINITE;
	int p, locElig, locCount;
	float nf, diff, top, bottom, totTop, totBottom, totinvh, totinv2;
	float invh, inv2;
	float varF, fBar, fminusp, fsprimep;
//	float temp, fBar=0, weigh = 0;
// count total eligible loci:
	for (p=0, locElig=0; p<nloci; p++) {
		bottom = *(fbottom+p);
		if (bottom <= 0) continue;	// this locus is not considered
		locElig++;
	};
	if (locElig <= 1) return 0;

	for (p=0, locCount=0, fBar=0, varF=0; p<nloci; p++) {
		bottom = *(fbottom+p);
		top = *(ftop+p);
// When locus p is ineligible, the denominator in f is 0.
		if (bottom <= 0) continue;	// this locus is not considered
		locCount++;
		totBottom = totFsBottom - bottom;
		totTop = totFsTop - top;
		totinvh = invhSum - (*(invhmonic+p));
		totinv2 = inv2Sum - (*(inv2count+p));
		invh = totinvh / (locElig - 1);
		inv2 = totinv2 / (locElig - 1);
		fminusp = totTop/totBottom;	// this is Fs without locus p
		fsprimep = FprimeS (fminusp, invh, inv2, census);
		nf = 1.0F/locCount;
		diff = fsprimep - fBar; // diff. between current f and previous mean
		fBar += nf*diff;
		if (locCount > 1)
//			varF += nf*(locCount-1)*(diff*diff);
// the following formula is the same, it uses updated fBar
			varF += diff*(fsprimep - fBar);
	};
	*fprimeBar = fBar;
// out of the loop, nf is still 1/locCount
	nf *= ((float)locCount-1);
	varF *= nf;
	return varF;
}

//------------------------------------------------------------------

void ConfidTempoJR (float totFsTop, float *ftop, float totFsBottom,
					float *fbottom, float invhSum, float *invhmonic,
					float inv2Sum, float *inv2count, int nloci,
					float *lowNe, float *highNe, float timegap,
					float infinite, float fprimeMean, int census, char mode)
// if mode = 0: use Chi-square, else: normal
// fprimeMean is the fprime calculated from all loci.
// Calculate confidence interval following Jorde/Ryman method, using Jackknife
// for finding variance of the mean.
// For normal approximation, find standard error s. Then the lower and upper
// limits of fprime, loLim and hiLim are given by fprimeBar -/+ 1.96s.
{
	float fprimeBar, varF, stdErr;
	float loLim, hiLim;
	long iBig;

	fprimeBar = fprimeMean;	// fprimeBar will be reassigned in VarTemJR
	varF = VarTemJR(totFsTop, ftop, totFsBottom, fbottom, invhSum,
					invhmonic, inv2Sum, inv2count, nloci, &fprimeBar, census);
//printf ("varF = %10.7f\n",varF);
//printf ("fprimeBar = %10.7f\n",fprimeBar);
//printf ("fprimeMean = %10.7f\n",fprimeMean);

	if (mode == 0) {
// Use chi-square distribution to have 95%CI:
	// choose fprimeMean, the value when all loci are used, or
	// fprimeBar: the mean from doing Jackknife (mean among values calculated
	// by dropping one locus).
	// In this algo, the two may be different.
		iBig = JackKnifeInd(fprimeMean, varF);
		Confid95 (iBig, fprimeMean, &loLim, &hiLim);
	} else {
// Use normal distribution to have 95%CI:
		stdErr = (float) sqrt(varF);
		loLim = fprimeBar - 1.96F * stdErr;
		hiLim = fprimeBar + 1.96F * stdErr;
	};
	// loLim and hiLim already have census factored in thru fprimeBar, which
	// is calculated in VarTemJR.
	*lowNe = NeFromFp (hiLim, timegap);
	*highNe = NeFromFp (loLim, timegap);
	if (*highNe <= 0 || *highNe > infinite) *highNe = infinite;
}

//------------------------------------------------------------------

void CIParamF (float f, float invSize, float timegap, long totA,
				   float *lowNe, float *highNe, float infinite, float extra)
// For parameter CI, using f, to get lower and upper bound for f.
// Then subtracted by invSize, which is 1/S, where S is sample size,
// to get lower and upper bound for fprime.
{
	float loLim, hiLim, adjust;
	Confid95 (totA, f, &loLim, &hiLim);	// CI for f
// print to console for checking
//printf ("PARAMETER CI:  f =%9.6f,  Deg. Free =%6lu,  f-CI: %9.6f,%10.6f\n",
//f, totA, loLim, hiLim);
	adjust = extra - invSize; // on input, extra=0 if plan 2, 1/census if plan 1
	loLim += adjust;	// becomes lower bound for fprime
	hiLim += adjust;	// becomes upper bound for fprime
	*lowNe = NeFromFp (hiLim, timegap);
	*highNe = NeFromFp (loLim, timegap);
//printf ("\tSubtracted by 1/S = %6.4f, f-CI implies f'-CI: %9.6f  %9.6f\n",
//invSize, loLim, hiLim);
//printf ("\tCI via f for Ne: %13.6f  %15.6f\n", *lowNe, *highNe);
	if (*highNe <= 0 || *highNe > infinite) *highNe = infinite;
}

//------------------------------------------------------------------

void CIParamFprime (float fprime, float timegap, long totA,
				   float *lowNe, float *highNe, float infinite)
// For parameter CI, using fprime.
{
	float loLim, hiLim;
	Confid95 (totA, fprime, &loLim, &hiLim);
	*lowNe = NeFromFp (hiLim, timegap);
	*highNe = NeFromFp (loLim, timegap);
	if (*highNe <= 0 || *highNe > infinite) *highNe = infinite;
// print to console for checking
//printf ("PARAMETER CI: f' = %9.6f, Deg. Free =%6lu, f'-CI: %9.6f  %9.6f\n",
//fprime, totA, loLim, hiLim);
//printf ("\tCI via f' for Ne: %12.6f  %15.6f\n\n", *lowNe, *highNe);
}

//------------------------------------------------------------------


long JacKnifeTemp (float *fArr, float *hSampMean, float *nAlle,
					int nloci, float fWtot, float totWeight,
					float timegap, char reWeigh)

//
{
	int p, locCount;
	float s, s2, wp, nf, diff, Neraw;
	float wtMinusp, varF, fBar, fMinusp;
	long iBig;
//	float invS = 0;
	// when reWeigh = 1, there is a need for reweight, then *invSize
	// will be recalculated as weighted mean of 1/hSampMean.
	// The weights used for having *invSize will be the ones for
	// having overall f, with full number of loci.
	// (overall f is not necessarily fBar calculated below)
	for (p=0, locCount=0, fBar=0, varF=0; p<nloci; p++) {
		wp = *(nAlle+p);	// Pollak: ind alleles, Nei/Tajima: total alleles
		if (wp <= 0) continue;	// this locus is not considered
		locCount++;
		if (reWeigh == 1) {	// initial weight includes samp size s^2
		// The wp here was calculated this way in TemporalNeEst, when
		// there are missing data, which results in having reWeigh = 1.
		// Therefore, totWeight will aways be the sum of all weights
		// (Pollak only)
			s = *(hSampMean+p);
			s2 = s*s;
			wp *= s2;	// weight at this locus
		};
		wtMinusp = totWeight - wp;	// total weights except at locus (p+1)
		// Since fWtot, as input parameter, represents total weighted
		// sum of fArr across loci, so subtracted by weighted fArr at
		// locus (p+1) will be weighted sum of all loci except locus (p+1),
		// then divide by total weights except locus (p+1) to get fMinusp:
		fMinusp = (fWtot - (*(fArr+p)) * wp)/wtMinusp;
		// recalculate fMinusp with new weights based on tentative Neraw
		if (reWeigh == 1) {	// fMinusp revised without locus (p+1):
			Neraw = NeFromFp (fMinusp, timegap);
			fMinusp = FvalRev (fMinusp, fArr, hSampMean, nAlle,
									Neraw, nloci, timegap, p);
		};

		nf = 1.0F/((float)locCount);
		diff = fMinusp - fBar;
		fBar += nf*diff;
		if (locCount > 1)
//			varF += nf*(locCount-1)*(diff*diff);
// the following formula is the same, it uses updated fBar
			varF += diff*(fMinusp - fBar);
	};
// out of the loop, nf is still 1/locCount, update to have (locCount-1)/locCount
	nf *= ((float)locCount-1);
	varF *= nf;
//printf ("fBar = %15.6f, variance = %15.8f\n", fBar, varF);
    iBig = JackKnifeInd(fBar, varF);
//printf ("locCount = %10d, iBig = %15d\n", locCount,  iBig);
	return iBig;
}


//------------------------------------------------------------------

void JackTempkc (float fMean, float *fArr, float invSize, float *hSampMean,
				   float *nAlle, int nloci, float fWtot, float totWeight,
				   float timegap, float *lowNe, float *highNe,
				   float infinite, char method)
// Jackknife for Pollak, Nei/Tajima using F-measure
// method = 1 is for Pollak with reweighting in effect
// hSampMean is array (index=locus) of harmonic mean of sample sizes
// of two samples.
// invSize is 1/S, where S is the weighted harmonic mean of sample sizes
// across loci
{
	float loLim, hiLim;
	long iBig = JacKnifeTemp (fArr, hSampMean, nAlle, nloci, fWtot,
								totWeight, timegap, method);
	Confid95 (iBig, fMean, &loLim, &hiLim);
// print to console for checking
//printf ("JACKKNIFE CI:\nUsing f, Deg. Free =%6lu, Mean f =%9.6f, f-CI:"
//"%9.6f,%10.6f\n", iBig, fMean, loLim, hiLim);
	loLim -= invSize;
	hiLim -= invSize;
//printf ("\tSubtracted 1/S =%7.4f, f'-CI: %9.6f%11.6f\n",invSize,loLim,hiLim);
	*lowNe = NeFromFp (hiLim, timegap);
	*highNe = NeFromFp (loLim, timegap);
//printf ("\tCI via f for Ne: %13.6f%17.6f\n", *lowNe, *highNe);
	if (*highNe <= 0 || *highNe > infinite) *highNe = infinite;
}


//------------------------------------------------------------------

void JackTempkcFp (float fprimeMean, float *fprime, float *hSampMean,
				   float *nAlle, int nloci, float fkprimeWtot,
				   float totWeight, float timegap, float *lowNe,
				   float *highNe, float infinite, char method)
// Jackknife for Pollak, Nei/Tajima using Fprime
// method = 1 is for Pollak with reweighting in effect
// hSampMean is array (index=locus) of harmonic mean of sample sizes
// of two samples.
{
	float loLim, hiLim;
	long iBig = JacKnifeTemp (fprime, hSampMean, nAlle, nloci,
								fkprimeWtot, totWeight, timegap, method);
	Confid95 (iBig, fprimeMean, &loLim, &hiLim);
	*lowNe = NeFromFp (hiLim, timegap);
	*highNe = NeFromFp (loLim, timegap);
// print to console for checking
//printf ("Using f', Deg. Freedom = %6lu,\n\tMean f' = %9.6f,\t "
//"f'-CI: %9.6f  %9.6f\n", iBig, fprimeMean, loLim, hiLim);
//printf ("\tCI via f' for Ne: %12.6f  %15.6f\n", *lowNe, *highNe);
	if (*highNe <= 0 || *highNe > infinite) *highNe = infinite;
}


//------------------------------------------------------------------
// Modified Dec 2016:
// At the "for" loop at Critical values, skip when the value is PCRITX

void TemporalNeEst (FILE *outLoc, char moreDat, FREQPTR *freqList, int nloci,
				  char *locUse, int g1, int g2, int nCrit, float critVal[],
				  long *nTotAlle, long *nIndAlle,
				  float *Hkmean, float *Hcmean, float *Hsmean,
				  float *fkmean, float *fcmean, float *fsmean,
				  float *fkprimeMean, float *fcprimeMean, float *fsprimeAll,
				  float *Nek, float *Nec, float *Nes,
				  float *loNek, float *hiNek, float *loNec, float *hiNec,
				  float *loNes, float *hiNes,
				  float *jloNek, float *jhiNek, float *jloNec, float *jhiNec,
				  float *jloNes, float *jhiNes, char param, char jack,
				  float timeline[], int census,
//				  float timeline[], int nGeneration,
				  char tempk, char tempc, char temps, float infinite,
				  char weighsmp)
{
// References:
// * R.Waples's "A Generalized Approach for Estimating Effective
//   Population Size From Temporal Changes in Allele Frequency" formulas
//   (8), (9), in Genetics 121: 379-391, February 1989.
//   This is used as reference for weighFk, weighFc, unweigFk, unweigFc.
// * Jordes & Ryman's "Unbiased Estimator for Genetic Drift and Effective
//   Population Size" formulas (3)-(5), (9), (13)
//   in Genetics 177: 927-935, October 2007.
//   This is used as reference for unweigFs, weighFs.
// Estimate Effective Population based on Fk, Fc, Fs between two generations
// g1 and g2.

// List of parameters: nTotAlle, nIndAlle, ... up to hiNes correspond to
// critical values critVal[nCrit]

	FREQPTR ptr1;
	float freq1, freq2, pbar, pmean, totsmf1, totsmf2, totsm, totsmbar, x,xx;
// topFs, bottomFs are sums of numerators, denominators use in calculating Fs
	float *topFs, *bottomFs;
	float totFsTop, totFsBottom, zz;
// add Jan 2012, overall fsprime:
// This could be in the list of param. (then declared in the calling func)
//	float *fsprimeAll;
//	float totnormbar;	// this might be deleted later
	int p, n, count1, count2, totcount, smAlle, locCount, alle, nA;
	float *indAlle, *nAlle;
// add Jan 2012 for the weight at loci in Pollak f-values. The weights will
// be reassigned after the first tentative Ne is calculated
//	float *locWtk;

// fk, fc, fs for F-values, hLoc for harmonic mean in two generations,
// nkLoc, ncLoc, nsLoc are effective breeders corresponding to F-values.
	float *fk, *fc, *fs, *hLoc, *nkLoc, *ncLoc, *nsLoc;
// add Feb 2012 for fprime in Pollak and Nei/Tajima:
// wfkprime is fkprime*weight, where the weights at loci are not just
// ind. alleles, but multiply by square of sample size.
	float *fkprime, fkprimeWtot, fkWtot, totwtK;
	float *fcprime, fcprimeWtot, fcWtot, totwtC;
//	float *wfkprime;
	float *ncMean, *nkMean, *nsMean;
	float hsampk, hsamp, harmonic, harmonick;
	float unweigFk, unweigFc, unweigFs;
	float weighFk, weighFc, weighFs;
	long totAlle, totInd;
	float coefc, coefk, sumfc, sumfk;// fc, fk at each allele, and sum at loc
	float bigsumfc, bigsumfk;		// sum all fc, fk in all loci
	float crit;
// The next 2 lines are for Jorde/Ryan method, added Dec 9, 11:
	float invhMean, inv2Mean, invh, inv2, invhSum, inv2Sum;
	float *invhmonic, *invcount2, *fsprime;
// variable ss is used for a choice of printing the last decimal place
// without rouding it up; comment out if that option is not to be used
//	char ss[20];
	float timegap;
// Add for Plan I:
	float extra = 0;	// for temporal Polak and Nei/Tajima
	if (census > 0) extra = 1.0F/((float) census);
	x = timeline[g1];
	xx = timeline[g2];
	timegap = (x > xx)? x-xx: xx-x;
	if (g1 + g2 == 1) {
		printf ("\nTemporal Method ... ");
		if (census > 0) printf ("Plan I, Census Size = %d\n", census);
		else printf ("Plan II\n");
	};
// hLoc holds the harmonic means of individuals having data in 2 generations
	hLoc = (float*) malloc (sizeof(float)*nloci);
// indAlle holds the number of alleles minus 1
	indAlle = (float*) malloc (sizeof(float)*nloci);
	nAlle = (float*) malloc (sizeof(float)*nloci);
// bring these 2 out from "if (temps==1)"
	invhmonic = (float*) malloc (sizeof(float)*nloci);
	invcount2 = (float*) malloc (sizeof(float)*nloci);
// to lessen memory allocation, only allocate when necessary:
	if (DETAILTEMP == 1) {	// Since details are given in Locdata out
		tempk = 1; tempc = 1; temps = 1;
	};
	if (tempk == 1) {
		fk = (float*) malloc (sizeof(float)*nloci);
// added Feb 2012:
		fkprime = (float*) malloc (sizeof(float)*nloci);
//		wfkprime = (float*) malloc (sizeof(float)*nloci);
		nkLoc = (float*) malloc (sizeof(float)*nloci);
		nkMean = (float*) malloc (sizeof(float)*nCrit);
		for (n = 0; n < nCrit; n++) {
			*(nkMean+n) = 0;
		};
// add Jan 2012:
//		locWtk = (float*) malloc (sizeof(float)*nloci);
	};
	if (tempc == 1) {
		fc = (float*) malloc (sizeof(float)*nloci);
// added Feb 2012:
		fcprime = (float*) malloc (sizeof(float)*nloci);
		ncLoc = (float*) malloc (sizeof(float)*nloci);
		ncMean = (float*) malloc (sizeof(float)*nCrit);
		for (n = 0; n < nCrit; n++) {
			*(ncMean+n) = 0;
		};
	};
	if (temps == 1) {
		fs = (float*) malloc (sizeof(float)*nloci);
// added Dec 9, 11 - Jan 2012:
		fsprime = (float*) malloc (sizeof(float)*nloci);

		topFs = (float*) malloc (sizeof(float)*nloci);
		bottomFs = (float*) malloc (sizeof(float)*nloci);
		nsLoc = (float*) malloc (sizeof(float)*nloci);
		nsMean = (float*) malloc (sizeof(float)*nCrit);
// add Jan 2012, overall fsprime:
// This could be in the list of param. (then declared in the calling funct.)
//		fsprimeAll = (float*) malloc (sizeof(float)*nCrit);
//		for (n = 0; n < nCrit; n++) {
//			*(fsprimeAll+n) = 0;
//		};
	};
	for (n = 0; n < nCrit; n++) {
		crit = critVal[n];

// Dec 2016: skip the special value reserved for dropping singleton alleles:
		if (crit > 0 && crit <= PCRITX) continue;

		locCount = 0;
		unweigFk = 0, unweigFc = 0, unweigFs = 0;
		weighFk = 0, weighFc = 0, weighFs = 0;
		totAlle = 0, totInd = 0;
		harmonic = 0, harmonick = 0;
		totFsTop = 0, totFsBottom = 0;
// add Jan 2012:
		inv2Sum = 0;
		invhSum = 0;
		bigsumfc = 0, bigsumfk = 0;
		totwtK = 0;
		fkprimeWtot = 0;
		fkWtot = 0;
		totwtC = 0;
		fcprimeWtot = 0;
		fcWtot = 0;
		if (outLoc != NULL && moreDat == 1) {
			if (n == 0) {
				fprintf (outLoc, "\n");
				for (p=0; p<45; p++) fprintf (outLoc, "-");
				fprintf (outLoc,
						"\nSamples %d and %d:   Generations%5.1f and %5.1f\n",
						g1+1, g2+1, timeline[g1], timeline[g2]);
				if (census > 0) fprintf (outLoc,
						"        (Plan I,   Census Size = %d)\n", census);
				else  fprintf (outLoc, "   (Plan II)\n");
				for (p=0; p<45; p++) fprintf (outLoc, "-");
				fprintf (outLoc, "\n");
			};
			if (crit > 0) fprintf (outLoc,
				"\nWith lowest frequency set at %5.3f:\n", crit);
			else fprintf (outLoc, "\nWithout restriction on frequencies:\n");
			for (p=0; p<35; p++) fprintf (outLoc, "-");
			fprintf (outLoc, "\n");
		};
		for (p=0; p<nloci; p++) {
			sumfk = 0;
			sumfc = 0;
			*(hLoc+p) = 0;
			*(indAlle+p) = 0;
			*(nAlle+p) = 0;	// # alleles whose freq.>=crit ("normal" alleles)
			nA = 0;
			if (tempk == 1) {
				*(fk+p) = 0;
				*(nkLoc+p) = 0;
			// add Feb 2012
				*(fkprime+p) = 0;
//				*(wfkprime+p) = 0;
// add Jan 2012, initialize weight at each locus:
//				*(locWtk+p) = 0;
			};
			if (tempc == 1) {
				*(fc+p) = 0;
			// add Feb 2012
				*(fcprime+p) = 0;
				*(ncLoc+p) = 0;
			};
// bring these 2 out from "if (temps==1)" so it can be used on others
			*(invhmonic+p) = 0;
			*(invcount2+p) = 0;
			if (temps == 1) {
				*(fs+p) = 0;
// added Dec 9, 11:
				*(fsprime+p) = 0;

				*(topFs+p) = 0;
				*(bottomFs+p) = 0;
				*(nsLoc+p) = 0;
			};
			if (*(locUse+p) == 0) continue;
			if (*(freqList+p) == NULL) continue;
			smAlle = 0;	// # alleles whose freq. < crit ("small" alleles)
			totsmf1 = 0, totsmf2 = 0, totsm = 0, totsmbar = 0;
//			totnormbar = 0;	// this may be deleted (not used at all)
// This is for details: controlled by constant DETAILTEMP
//			if (outLoc != NULL && MoreDat == 1 && DETAILTEMP == 1) {
//				fprintf (outLoc, "\nLocus %d\n", p+1);
//				fprintf (outLoc, "Allele    x       y    (x-y)^2   ");
//				fprintf (outLoc, "z=(x+y)/2   z(1-z)       Fc          Fk\n");
//			};

			for (ptr1 = *(freqList+p), count1 = 2*(ptr1->samples[g1]),
				count2 = 2*(ptr1->samples[g2]), totcount = count1 + count2;
				ptr1 != NULL; ptr1 = ptr1->next)
			{
				alle = ptr1->mValue;
				freq1 = ptr1->freqs[g1];
				freq2 = ptr1->freqs[g2];
				if (freq1 == 0 && freq2 == 0) continue;
			// this allele is in neither generation. Note that the allele
			// list consists of all alleles that appear in any generations.
			// There may be some allele in just one generation only.

			// Note: count1, count2 are 2*the number of samples having data
			// at locus (p+1): they are unchanged from allele to allele.
			// Thus, we can assign only once, for each locus, and the
			// assignment is enclosed in "for" statement.
//				count1 = 2*(ptr1->samples[g1]);
//				count2 = 2*(ptr1->samples[g2]);
//				totcount = count1 + count2;
// pbar is the weighted mean of freq1 and freq2
// pbar is used to compare with crit value, to determine if the allele
// should be dropped
				pbar = (freq1*count1 + freq2*count2)/((float)totcount);
				pmean = (freq1+freq2)/2;
				if (pbar < crit) {	// if pbar = 0, both freq1, freq2 = 0
					smAlle++;
					totsmf1 += freq1;
					totsmf2 += freq2;
					totsmbar += pbar;
					totsm += pmean;
				} else {
					// need to have condition pbar==1 because crit can be 0
					if (pbar>(1-crit) || pbar==1) break; // almost mono,
				// If pbar = 1, the locus is monomorphic in both generations
				// for the same allele.
				// In case pbar>(1-crit), all other alleles have pbar < crit,
				// so nAlle is still 0. This locus will be skipped.
				// Note: if this locus is mono. at both generations, but
				// for different alleles, then pbar for each is < 1, so the
				// locus may still be in play.
//					(*(nAlle+p))++;
					nA++;
					x = freq1 - freq2;
					xx = x*x;
/*
					*(fk+p) += xx/pbar;	// a general term in Waples's (9)
					*(fc+p) += xx/(pbar - freq1*freq2); // the right side is
				// a general term in Waples's (8) or Jorde/Ryman's (1).
				// Its denominator is always > 0 since pbar < 1 in this case.
				// (pbar is between freq1, freq2. If it equals to one of them,
				// then it must also equal to both. Since both freq1, freq2
				// are <= 1, its product freq1*freq2 <= both, so the product
				// is <= pbar. It cannot be equal to pbar, otherwise, one
				// of freq1 or freq2 must be 1, say freq1 =1. Then pbar=freq2,
				// which then implies pbar=freq1 also, hence pbar=1!
					topFs += xx;
					bottomFs += pbar*((float)1.0 - pbar);
			// this variable totnormbar does not seem to be used at all:
					totnormbar += pbar;
*/
// use weighted mean pmean instead of unweighted mean pbar:
					coefk = xx/pmean;			// Waples' (9)
					coefc = xx/(pmean - freq1*freq2);
					sumfk += coefk;
					sumfc += coefc;
							// Waples' (8) or Jorde & Ryman's (1)
//					if (tempk == 1) *(fk+p) += coefk;
//					if (tempc == 1) *(fc+p) += coefc;

					zz = pmean*((float)1.0 - pmean);
			// values top, bottom are for fs (F-values in  Jorde & Ryman).
					if (temps == 1) {
						*(topFs+p) += xx;
						*(bottomFs+p) += zz;
					};
// This is for details: controlled by DETAILTEMP
					if (outLoc != NULL && moreDat == 1 && DETAILTEMP == 1) {
						if (nA == 1) {
							fprintf (outLoc, "\nLocus %d\n", p+1);
							fprintf (outLoc, "Allele    x       y    (x-y)^2   ");
							fprintf (outLoc,
									"z=(x+y)/2   z(1-z)       Fc          Fk\n");
						};
						fprintf (outLoc,
						"%4d%9.4f%8.4f%9.5f%10.4f%11.5f%12.6f%12.6f\n",
						alle, freq1, freq2, xx, pmean, zz, coefc, coefk);
					};
				};
			}; // gone thru all alleles at locus (p+1): for (ptr1=..."
		// when nAlle = 0, it is either an almost monomorphic locus,
		// or every allele is of frequency < crit.
//			if (*(nAlle+p) == 0) continue;	// this locus is skipped: fk,fc=0.
			if (nA == 0) continue;	// this locus p is skipped: fk,fc=0.
		// From now, there is at least one normal allele <= 1-crit
			locCount++;
			if (smAlle > 0) {	// for all alleles whose frequencies
			// across two generations being small than "crit", we consider
			// as one "normal" allele:
			// totsmf1, totsmf2 are frequencies of all of those,
			// totsmbar is the frequency of this combined allele
			// as it is considered across two generations
//				(*(nAlle+p))++;	// so this should be at least 2.
				nA++;	// so this should be at least 2.
				x = totsmf1 - totsmf2;
				xx = x*x;
/*
				*(fk+p) += xx/totsmbar;
				*(fc+p) += xx/(totsmbar - totsmf1*totsmf2);
				topFs += xx;
				bottomFs += totsmbar*((float)1.0 - totsmbar);
*/
			// the denominator cannot be zero since totsmbar < 1
			// as in the case explained before.
				coefk = xx/totsm;
				coefc = xx/(totsm - totsmf1*totsmf2);
				sumfk += coefk;
				sumfc += coefc;
				zz = totsm*((float)1.0 - totsm);
//				if (tempk == 1) *(fk+p) += coefk;
				if (temps == 1) {
					*(topFs+p) += xx;
					*(bottomFs+p) += zz;
				};
// This is for details: controlled by DETAILTEMP
				if (outLoc != NULL && moreDat == 1 && DETAILTEMP == 1)
					fprintf (outLoc,
						"%4s%9.4f%8.4f%9.5f%10.4f%11.5f%12.6f%12.6f\n",
						"rare", totsmf1, totsmf2, xx, totsm, zz, coefc, coefk);
			};	// end of "if (smAlle > 0)"
			if (tempk == 1) *(fk+p) = sumfk;
			if (tempc == 1) *(fc+p) = sumfc;
// This is for details: controlled by DETAILTEMP
			if (outLoc != NULL && moreDat == 1 && DETAILTEMP == 1) {
			// these are for outLocting to auxiliary file only:
				bigsumfk += sumfk;	// note that when DETAILTEMP = 1, temps =1
				bigsumfc += sumfc;	// too, so topFs, bottomFs are allocated.
				fprintf (outLoc,
						"%4s%9s%8s%9.5f%10s%11.5f%12.6f%12.6f\n",
						"SUM:", ">>>>>>", ">>>>>>", *(topFs+p), ">>>>>>",
						*(bottomFs+p), sumfc, sumfk);
			};
			*(nAlle+p) = (float) nA;
			*(indAlle+p) = *(nAlle+p) - 1;
			inv2 = (float) 1.0/count2;
			invh = (float) 1.0/count1 + inv2;
			if (temps == 1) {
				totFsTop += *(topFs+p);
				totFsBottom += *(bottomFs+p);
			};
// modified Dec 9, 11, to hold values for:
// * 1/count2, which is 1/2*(individuals having data at sample 2),
// * 1/N - where N = harmonic mean of individuals in two samples having data
// (Prev.: *(hLoc+p) is in place of invhmonic then changed to its inverse)
			*(invcount2+p) = inv2;
			*(invhmonic+p) = invh;
			inv2Sum += inv2;
			invhSum += invh;
			// then weighted by ind. alleles:
			hsampk = (*(indAlle+p)) * invh ;
			hsamp = (*(nAlle+p)) * invh ;
			*(hLoc+p) = (float) 1.0/invh;		// *(hLoc+p) now is
												// the harmonic mean
			harmonic += hsamp;
			harmonick += hsampk;
			if (temps == 1) {
				if ((*(bottomFs+p)) > 0) *(fs+p) = (*(topFs+p))/(*(bottomFs+p));
				*(fsprime+p) = FprimeS (*(fs+p), invh, inv2, census);
				// Eff. breeder based on Fs at this locus
				*(nsLoc+p) = NeFromFp (*(fsprime+p), timegap);
				x = *(nAlle+p);
				*(nsMean+n) += x / (*(nsLoc+p));
// these 2: weighFs, unweigFs, might not be needed:
				weighFs += (*(fs+p)) * x;
				unweigFs += *(fs+p);
			};
			if (tempk == 1) {
			// at this step, (fk+p) is the sum of all f-values at locus (p+1)
				x = *(indAlle+p);
				if (x > 0) {
					weighFk += *(fk+p);	// sum of all f-values at all alleles
			// with K-1=indAlle=x, s=hLoc (sample size at this locus), then
			// weight for fkprime is taken to be (K-1)*s^2.
			// If sample sizes are the same across loci (weighsmp = 0), don't
			// need to include them in the weights (since all s are the same)
					xx = *(hLoc+p);	// xx represents the weight at this locus
					if (weighsmp == 0) xx = x;	// (K-1)
					else xx *= (x*xx);			// (K-1)*s^2
					totwtK += xx;
			// dividing by x, (fk+p) will be the Pollak f-value for locus p+1:
					*(fk+p) /= x;
					*(fkprime+p) = *(fk+p) - invh; // fprime-value for loc p+1
				// for plan I:
					*(fkprime+p) += extra;	// = 1/census, at the beginning
				// Eff. breeder based on Fk at this locus
					*(nkLoc+p) = NeFromFp (*(fkprime+p), timegap);
					*(nkMean+n) += xx / (*(nkLoc+p));
					fkprimeWtot += (*(fkprime+p))*xx;	// sum(fprime*weight)
					fkWtot += (*(fk+p))*xx;	// sum(fk*weight)
				} else {
					*(fk+p) = 0; *(fkprime+p) = 0;
				};
				unweigFk += *(fk+p);
			};
			if (tempc == 1) {
				weighFc += *(fc+p);
				x = *(nAlle+p);
				xx = x;		// xx is to denote weight at locus (p+1)
				totwtC += xx;
				*(fc+p) /= x;	// Nei/Tajima f-value
				*(fcprime+p) = *(fc+p) - invh; // fprime-value
				// for plan I:
				*(fcprime+p) += extra;	// = 1/census, given at the beginning
				fcprimeWtot += (*(fcprime+p))*xx;
				fcWtot += (*(fc+p))*xx;	// sum(fc*weight)
				unweigFc += *(fc+p);
				// Eff. breeder based on Fc at this locus
				*(ncLoc+p) = NeFromFp (*(fcprime+p), timegap);
				*(ncMean+n) += xx / (*(ncLoc+p));
			};
			totInd += (int) (*(indAlle+p));	// total ind. alleles
			totAlle += (int) (*(nAlle+p));	// total all alleles
		};
		// ------------------------------------------------------------------
		// end of "for (p=0; ..." gone thru all loci, now calculate means
		if (harmonic > 0) harmonic = (float) totAlle/harmonic;
		if (harmonick > 0) harmonick = (float) totInd/harmonick;
		*(Hkmean+n) = harmonick;
		*(Hcmean+n) = harmonic;
		*(Hsmean+n) = harmonic;
		if (locCount > 0) {	// these 3 may not be used
			if (tempk == 1) unweigFk = unweigFk/((float) locCount);
			if (tempc == 1) unweigFc = unweigFc/((float) locCount);
			if (temps == 1) unweigFs = unweigFs/((float) locCount);
		// but this is used, added Jan 2012:
		// This is the inverse of unweighted harmonic mean, taken across loci
		// of all harmonic means of two samples taken at each locus.
		// That is, 1/invhMean is the unweighted harmonic mean across loci.
			invhMean = invhSum / locCount;
		// Similarly for harmonic mean across loci of sample 2.
			inv2Mean = inv2Sum / locCount;
		};
	// weighFk is the weighted mean of (fk+p) across loci, where
	// the weight of each locus is the number of its independent alleles
		if (totInd > 0) {
			if (tempk == 1) {
				*(fkmean+n) = weighFk / ((float) totInd);
//if (fabs(*(nkMean+n)) > EPSILON)
				*(nkMean+n) = totwtK / (*(nkMean+n));
				*(fkprimeMean+n) = fkprimeWtot/totwtK;	// w. mean of fkprime
				// with weights = (K-1)s^2 or (K-1) as described above.
//				x = *(fkmean+n) - invhMean;	// = weighted mean of fkprime if
									// the weight is number of ind. alleles
				*(Nek+n) = NeFromFp (*(fkprimeMean+n), timegap);
// add Feb 2012, recalculate Ne:
				if (weighsmp > 0) {
//x = 1.0F/(*(Hkmean+n));
//if (outLoc != NULL && moreDat == 1)
//fprintf(outLoc, "Before reweigh, fk = %10.8f, fkprime = %10.8f, fk - fkprime"
//" = %10.8f, 1/(Harmonic Size) = %12.8f\n",
//*(fkmean+n), *(fkprimeMean+n), *(fkmean+n) - *(fkprimeMean+n), x);
				// recalculate Fprime for final Ne
					*(fkprimeMean+n) = FvalRev (*(fkprimeMean+n), fkprime, hLoc,
								indAlle, *(Nek+n), nloci, timegap, nloci);
				// recalculate Fk for outLoc
					*(fkmean+n) = FvalRev (*(fkmean+n), fk, hLoc,
								indAlle, *(Nek+n), nloci, timegap, nloci);
// recalculate Harmonic Sample Size for outLoc
//x = FvalRev (x, invhmonic, hLoc, indAlle, *(Nek+n), nloci, timegap, nloci);
					*(Nek+n) = NeFromFp (*(fkprimeMean+n), timegap);
//if (outLoc != NULL && moreDat == 1)
//fprintf(outLoc, "After reweigh,  fk = %10.8f, fkprime = %10.8f, fk - fkprime"
//" = %10.8f, 1/(Harmonic Size) = %12.8f\n",
//*(fkmean+n), *(fkprimeMean+n), *(fkmean+n) - *(fkprimeMean+n), x);
				};
			};
		};	// end of "if (totInd > 0)"
	// weighFc is the weighted mean of (fc+p) across loci, where
	// the weight of each locus is the number of its alleles
		if (totAlle > 0) {
			if (tempc == 1) {
				*(fcmean+n) = weighFc / ((float) totAlle);
				*(ncMean+n) = totwtC / (*(ncMean+n));
				*(fcprimeMean+n) = fcprimeWtot/totwtC;	// w. mean of fcprime
//				x = *(fcmean+n) - invhMean;	// = weighted mean of fcprime
				*(Nec+n) = NeFromFp (*(fcprimeMean+n), timegap);
			};
			// this is what Jorde/Ryan suggested for overall value Fs:
			if (temps == 1) {
				*(fsmean+n) = totFsTop / totFsBottom;
				*(nsMean+n) = (float) totAlle / (*(nsMean+n));
// add Jan 2012:
				*(fsprimeAll+n) = FprimeS (*(fsmean+n), invhMean, inv2Mean, census);
				*(Nes+n) = NeFromFp (*(fsprimeAll+n), timegap);
			};
		};	// end of "if (totAlle > 0)"
		*(nTotAlle+n) = totAlle;
		*(nIndAlle+n) = totInd;
	// now we can print fk, fc at each locus before releasing the memory
		if (outLoc != NULL && moreDat == 1) {
/*
			if (n == 0) {
				fprintf (outLoc, "\n");
				for (p=0; p<45; p++) fprintf (outLoc, "-");
				fprintf (outLoc,
						"\nSamples %d and %d:   Generations%5.1f and %5.1f\n",
						g1+1, g2+1, timeline[g1], timeline[g2]);
				for (p=0; p<45; p++) fprintf (outLoc, "-");
				fprintf (outLoc, "\n");
			};
			if (crit > 0) fprintf (outLoc,
				"\nWith lowest frequency set at %5.3f:\n", crit);
			else fprintf (outLoc, "\nWithout restriction on frequencies:\n");
			for (p=0; p<35; p++) fprintf (outLoc, "-");
			fprintf (outLoc, "\n");
*/
			if (DETAILTEMP == 1) {
				fprintf (outLoc, "\nSUMMARIZE%20.4f%21.4f%13.5f%12.5f\n",
						totFsTop, totFsBottom, bigsumfc, bigsumfk);
				fprintf (outLoc, "%23sTotal Alleles & Ind. Alleles:", " ");
				fprintf (outLoc, "%11d%12d\n", (int)(totAlle+0.5F),
								(int)(totInd+0.5F));
//				fprintf (outLoc, "\nSUMMARIZE%21.5f%21.5f%12s%12.6f\n",
//						totFsTop, totFsBottom, "Fs/All =", *(fsmean+n));
			};
			fprintf (outLoc, "\nLocus   H.Mean");
			if (tempk == 1) fprintf (outLoc, "%9s%12s", "Fk ", "Nb_k");
			if (tempc == 1) fprintf (outLoc, "%9s%12s", "Fc ", "Nb_c");
			if (temps == 1) fprintf (outLoc, "%9s%10s%12s", "Fs ", "Fs'", "Nb_s");
			fprintf (outLoc, "\n");
			for (p=0; p<nloci; p++) {
				if (*(locUse+p) == 0) continue;
				// put an asterisk to locus if it is dropped by freq required
				if (*(nAlle+p)==0)
					fprintf (outLoc, "%5d*%8.1f", (p+1), *(hLoc+p));
				else
					fprintf (outLoc, "%5d%9.1f", (p+1), *(hLoc+p));
				if (tempk == 1) fprintf (outLoc, "%11.6f%10.1f",
										*(fk+p), *(nkLoc+p));
				if (tempc == 1) fprintf (outLoc, "%11.6f%10.1f",
										*(fc+p), *(ncLoc+p));
				if (temps == 1) fprintf (outLoc, "%11.6f%10.6f%10.1f",
										*(fs+p), *(fsprime+p), *(nsLoc+p));
				fprintf (outLoc, "\n");
			};
//			fprintf (outLoc, "\nMean\n");
/* Comment out the printing of unweighted means:
			fprintf (outLoc, "* Unweighted: ");
			if (tempk == 1) fprintf (outLoc, "%-12.6f", unweigFk);
		    if (tempc == 1) fprintf (outLoc, "%-12.6f", unweigFc);
		    if (temps == 1) fprintf (outLoc, "%-12.6f", unweigFs);
		    fprintf (outLoc, "\n");
*/
			if (tempk+tempc+temps > 0) {
				fprintf (outLoc, "\nWeighted Mean:");
				if (tempk == 1) fprintf (outLoc, "%11.6f%10.1f",
							*(fkmean+n), *(nkMean+n));
				if (tempc == 1) fprintf (outLoc, "%11.6f%10.1f",
							*(fcmean+n), *(ncMean+n));
// block out this:
//				if (temps == 1) fprintf (outLoc, "%21s%10.1f",
//							" ", *(nsMean+n));
				fprintf (outLoc, "\n");
			};
			if (temps == 1) {
				fprintf (outLoc, "Fs overall:%14.6f\n", *(fsmean+n));
				fprintf (outLoc, "Fs' overall:%13.6f\n", *(fsprimeAll+n));
			};
//			if (temps==1) fprintf(outLoc,"Fs overall:%14.6f\n",*(fsmean+n));
		};	// end of "if outLoc != NULL", outLoc is for the Locus data
		// For parametric confidence intervals, use Chi-square distribution.
		// Take totInd = number of Indep. alleles as degree of freedom
		// Find CI either for F, or for F', then translated to CI for Ne.
		// Use CIParamF to find CI for F, CIParamFprime to find CI for F'.
		// (In Jorde/Ryman method, use CIParamFprime.)
		// For jackknife CI, either jackknife on values of F or F'.
		// On values of F, use JackTempkc; on values of F', use JackTempkcFp.
		// Jorde/Ryman method is treated differently.)
		if (tempk==1) {
// print to console for checking
//printf ("\nPollak CIs\n");
			x = 1.0F/(*(Hkmean+n));
			if (param == 1) {
				CIParamF (*(fkmean+n), x, timegap, totInd, (loNek+n),
							(hiNek+n), infinite, extra);
//				CIParamFprime (*(fkprimeMean+n), timegap, totInd, (loNek+n),
//							(hiNek+n), infinite);
			};
			if (jack == 1) {
				if (weighsmp > 0)
					x = FvalRev (x, invhmonic, hLoc, indAlle, *(Nek+n),
								nloci, timegap, nloci);
				JackTempkc (*(fkmean+n), fk, x, hLoc, indAlle, nloci,
					fkWtot, totwtK, timegap, (jloNek+n), (jhiNek+n),
					infinite, weighsmp);	// last = 1 if reweight in Pollak
//				JackTempkcFp (*(fkprimeMean+n), fkprime, hLoc, indAlle, nloci,
//					fkprimeWtot, totwtK, timegap, (jloNek+n), (jhiNek+n),
//					infinite, weighsmp);	// last = 1 if reweight in Pollak
			};
		};
		if (tempc==1) {
// print to console for checking
//printf ("\nNei/Tajima CIs\n");
			x = 1.0F/(*(Hcmean+n));
			if (param == 1) {
				CIParamF (*(fcmean+n), x, timegap, totInd, (loNec+n),
							(hiNec+n), infinite, extra);
//				CIParamFprime (*(fcprimeMean+n), timegap, totInd, (loNec+n),
//							(hiNec+n), infinite);
			};
			if (jack == 1) {
				JackTempkc (*(fcmean+n), fc, x, hLoc, nAlle, nloci,
					fcWtot, totwtC, timegap, (jloNec+n), (jhiNec+n),
					infinite, 0);
//				JackTempkcFp (*(fcprimeMean+n), fcprime, hLoc, nAlle, nloci,
//					fcprimeWtot, totwtC, timegap, (jloNec+n), (jhiNec+n),
//					infinite, 0);	// last: 0 = not Pollak
			};
		};
		if (temps==1) {
// print to console for checking
//printf ("\nJorde/Ryman CIs\n");
			if (param == 1) {
				CIParamFprime (*(fsprimeAll+n), timegap, totInd, (loNes+n),
							(hiNes+n), infinite);
			};
			if (jack == 1) {
				ConfidTempoJR (totFsTop, topFs, totFsBottom, bottomFs,
								invhSum, invhmonic, inv2Sum, invcount2,
								nloci, (jloNes+n), (jhiNes+n),
								timegap, infinite, *(fsprimeAll+n), census, 1);
		// the last constant, 0: uses chi-square, 1: normal (as in Jorde/Ryan)
			};
		};

		if (outLoc != NULL && moreDat == 1) {
			fprintf (outLoc, "Effective Pop:");
			if (tempk == 1) fprintf (outLoc, "%21.1f", *(Nek+n));
			if (tempc == 1) fprintf (outLoc, "%21.1f", *(Nec+n));
			if (temps == 1) fprintf (outLoc, "%21.1f", *(Nes+n));
			fprintf (outLoc, "\n");
		};
	};
	free (invhmonic);
	free (invcount2);
	if (tempk == 1) {
		free (fk);
		free (fkprime);
		free (nkLoc);
		free (nkMean);
	};
	if (tempc == 1) {
		free (fc);
		free (fcprime);
		free (ncLoc);
		free (ncMean);
	};
	if (temps == 1) {
		free (fs);
		free (fsprime);
//		free (fsprimeAll);
		free (topFs);
		free (bottomFs);
		free (nsLoc);
		free (nsMean);
	};
	free (hLoc);
	free (indAlle);
	free (nAlle);

}

//------------------------------------------------------------------




//-------------------------------------------------------------------------
// Print
// --------------------------------------------------------------------------

/*
void PrtLines (FILE *output, int ndash, char dash)
{
	int n;
	if (output == NULL) return;
	for (n=0; n<ndash; n++) fprintf (output, "%c", dash);
	fprintf(output, "\n");
	fflush (output);

}
//*/
// --------------------------------------------------------------------------

int PrtMisDat(FILE *missDat, int m, int hiErr, int samp,
					int noGen, char genErr[], int firstErr)
// return nonzero if serious error in genotype data, which is determined
// by parameter m, whose values are obtained from function GetSample.
// Leave the decision whether to abandon the rest of input file or not
// to the calling of this.
// Note: firstErr is the first locus that has missing data,
// hiErr is the locus that has the most serious error for the sample.
// Both parameters are indices of locus array which is 0-based,
// so need to add 1 for the output.
{
	if (noGen <= 0) return 0;
	if (missDat != NULL) {
		if (hiErr > firstErr)
			fprintf (missDat,  " %7d %8d,%7d   %10s   %11d\n",
					samp, firstErr+1, hiErr+1, genErr, noGen);
		else
			fprintf (missDat, " %7d %12d       %10s   %11d\n",
					samp, hiErr+1, genErr, noGen);
		fflush (missDat);
	};
	if (m == 1) return 0;
	else return (m-2);
}
//--------------------------------------------------------------------------


void PrtMisLabel (FILE *missDat, int popRead, char popID[])
{
	if (missDat ==  NULL ) return;
	fprintf (missDat, "Population %d [%s]\n", popRead, popID);
	PrtLines (missDat, 59, '-');
	fprintf (missDat, "Individual       Locus         Genotype     "
						"Number of Loci\n%41c with missing data\n", ' ');
	fflush (missDat);
}

// --------------------------------------------------------------------------

void PrtSumMisDat (FILE *missDat, int popRead, int nErr, char newID[], int next)
{
	if (missDat == NULL) return;
	PrtLines (missDat, 59, '-');
	fprintf (missDat,
		"Total missing data for population%5d: %12d\n\n", popRead, nErr);

	if (next != -1) PrtMisLabel (missDat, popRead+1, newID);
	fflush (missDat);
}

//--------------------------------------------------------------------------

FILE *PrtMisHead (char missFileName[], char inpName[], int popRead, char newID[])
{
	FILE *missDat;
	if ((missDat=fopen (missFileName, "w")) == NULL) return NULL;
	fprintf (missDat, "Missing data from input file %s.\n\n", inpName);
	fprintf (missDat, "Possible four types of missing data at a locus:\n");
	fprintf (missDat, "\t1. Genotype contains only zeros or partially scored.\n");
	fprintf (missDat, "\t2. Genotype has less digits than normal one.\n");
	fprintf (missDat, "\t3. Genotype has more digits than normal one.\n");
	fprintf (missDat, "\t4. Genotype contains non-digit character.\n");
	fprintf (missDat, "Types 3 and 4 stop the program.\n\n");
	fprintf (missDat,
		"In the table, each row is for an individual with missing data\n"
		"(a) If column 'Locus' has only one number, then it is the first\n"
		"    locus with missing data and also of highest missing data type.\n");
	fprintf (missDat,
		"(b) If column 'Locus' has 2 numbers, then the first number is\n"
		"    the first locus with data missing, and the second number is\n"
		"    the first locus that has highest missing data type.\n");
	fprintf (missDat,
		"(c) Genotype column contains the genotype of the locus in case (a)\n"
		"    or the second locus in case (b).\n\n");
	PrtMisLabel (missDat, popRead, newID);
	return missDat;
}

//--------------------------------------------------------------------------


int PrtError (FILE *output, FILE *missDat, int nloci, int nSampErr,
			int popRead, int samp, char popID[], int err, int noGen,
			char genErr[], int firstErr)
{
	int p, m;
	int errCode;
	if (err == 0) return 0;
	if (err == -1) {
		if (missDat != NULL) {
			fprintf (missDat, "Population %d: Sample %d ends too soon.\n",
					popRead, samp);
			fprintf (output, "\nPopulation %d: Sample %d ends too soon.\n",
					popRead, samp);
		};
		errCode = 3;
		return errCode;
	};	// now, err > 0, it is of the form m*nloci + p, p=0,...,nloci-1.
	m = err/nloci; p = err % nloci;
	errCode = PrtMisDat(missDat, m, p, samp, noGen, genErr, firstErr);
	if (errCode != 0) {	// quit running with this input file.
	// since we quit, should output error messages to output file
		if (errCode == 1)
			fprintf (output, "\nFatal error: At locus %d, "
				"Sample %d (population %d [%s]) has too many characters"
				" for a genotype.\n", p+1, samp, popRead, popID);
		else	// errCode = 2
			fprintf (output, "\nFatal error: At locus %d, "
				"Sample %d (population %d [%s]) has non-digit character"
				" for a genotype.\n", p+1, samp, popRead, popID);
		fflush (output);
	};
	return errCode;
}

// --------------------------------------------------------------------------

void PrtVersion (FILE *output) {
//		fprintf (output, "Draft output from NeEstimator v.2 (Beta)\n");
		fprintf (output, "Output from NeEstimator v.2\n");
}

// --------------------------------------------------------------------------

void PrtHeader (FILE *output, char append, char *inpName, int icount, int outype)
// outype = 0 when this is the short output; otherwise, the main output
{
	time_t rawtime;
	if (output == NULL) return;
	// print line of "=" corresponding to 4 critial values = 26+12*4 = 74
	if (append == 1) {
		fprintf (output, "\n");
		if (outype == 0) PrtLines (output, 77, '=');
		else PrtLines (output, 74, '=');
	} else {
		PrtVersion (output);
//		time ( &rawtime );
//		fprintf (output, "Starting time: %s", ctime (&rawtime));
	};
	time ( &rawtime );
	fprintf (output, "Starting time: %s", ctime (&rawtime));
	if (outype != 0) printf ("Starting time: %s", ctime (&rawtime));

	fprintf (output, "Input File");
	if (icount > 0) fprintf (output, " #%d", icount);
	fprintf (output, ": \"%s\"\n\n", inpName);
	fflush (output);
}
// --------------------------------------------------------------------------

int PrtLimitUse (FILE *output, char *locUse, int nloci, char byRange,
				int popStart, int popEnd, int nPop, int maxSamp, char term[])
// return the number of loci to be used, print to output info on
// populations, samples, and loci, which are limited
// Currently, when this is called, parameter nPop = popEnd, so there will
// be no printing for population range. Can change later if necessary.

{
	int m, n, p;
	int locSt, k, num;
	n = 0;
	if (output == NULL) return n;
	locSt = nloci;
	// m will be the number of loci being dropped, n is the number used.
	m = 0;
	for (p=0; p<nloci; p++) {
		if (*(locUse+p) == 0) m++;
		else n++;
	};
	if (popEnd < nPop) {	// populations to run are limited
		if (popStart == 1) {
			if (popEnd == 1) fprintf (output, "Only run for %s 1\n", term);
			else fprintf (output, "Run up to %s %d \n", term, popEnd);
		} else {
			if (popStart < popEnd) fprintf (output,
				"Limit to %ss from %d to %d \n", term, popStart, popEnd);
			else fprintf (output, "Only run for %s %d\n", term, popEnd);
		};
	} else if (popStart > 1) {
		fprintf (output, "Run from %s %d \n", term, popStart);
	};

	if (maxSamp < MAX_SAMP) fprintf (output,
			"Up to %d individuals are processed per %s.\n", maxSamp, term);
	fprintf (output, "Number of Loci = %d\n", nloci);
	if (m > 0) { // there are loci being dropped
		fprintf (output, "Number of loci being dropped: %d\n", m);
		if (n == 0) {
			fflush (output);
			return 0;
		};
		if (byRange == 1) {
			fprintf (output, "Loci in Use: ");
			locSt = 0;
			k = 0;	// # loci accepted, starting from locSt, left endpt of a range
			num = 0;	// number of ranges
			for (p=0; p<nloci; p++) {
				if (*(locUse+p) == 0) {
					if (k > 0) {	// k>0 when loci accepted were traveled
					// and this is the first locus that is rejected
						if (num > 0 && num % 10 == 0) fprintf (output, "\t\n");
						else if (num > 0) fprintf (output, ", ");
						if (k==1)
							fprintf (output, " %d", locSt+1);
						else
							fprintf (output, " %d - %d", locSt+1, locSt+k);
						num++;
					};
					k = 0;	// reset k for the number of loci in the next range
					locSt++;
				} else { // locus (p+1) is in the range
					if (k == 0) locSt = p;	// the first locus accepted
					k++;
			// Fix in June 3, 2013, so that if the last range including the
			// last locus, the range will be printed.
					// if this p happens to be the last locus, it means that
					// we are at the last range.
					if (p == nloci-1) {
						if (num > 0 && num % 10 == 0) fprintf (output, "\t\n");
						else if (num > 0) fprintf (output, ", ");
						if (k==1)
							fprintf (output, " %d", locSt+1);
						else
							fprintf (output, " %d - %d", locSt+1, locSt+k);
					};
				}
			};
//			if (num % 10 != 1)
			fprintf (output, "\n");
		} else {	// list loci dropped;
			fprintf (output, "Loci dropped:");
			num = 0;	// number of loci being dropped
			for (p=0; p<nloci; p++) {
				if (*(locUse+p) == 0) {
					if (num > 0 && num % 12 == 0) fprintf (output, "\n%13s"," ");
					else if (num > 0) fprintf (output, ", ");
					fprintf (output, " %5d", (p+1));
					num++;
				};
			};
//			if (num % 12 != 1)
			fprintf (output, "\n");
		};
	};
	fprintf (output, "\n");
	fflush (output);
	return n;
}


// --------------------------------------------------------------------------
// add in Mar 2013 for printing monomorphic loci
// add in apr 2014 to print popID if only temporal method (single = 0)
void PrtMonoLoc (FILE *output, int nloci, int nMobil[], char *locUse,
				char *popID, char single)
{
	int p, n, i, line;
	int perLine = 10;
	if (output == NULL) return;
	i = 0;
	n = 0;		// count number of monomorphic
	line = 0;	// count number of lines printed
	for (p=0; p<nloci; p++) {
		if (*(locUse+p) == 0) continue;
		if (*(nMobil+p) <= 1) {
			if (n == 0) {
				if (single == 0) fprintf (output, "Sample [%s]\n", popID);
				fprintf (output, "Non-polymorphic loci:");
			}
			if (n > 0) fprintf (output, ",");
			// the first line permits up to (perLine - 2) loci
			// (just in case perLine < 2, we use >= instead of ==)
			if ((i == perLine) || (line == 0 && i >= perLine-2)) {
			// indent by 7 blanks
				fprintf (output, "\n%7s", " ");
				i = 1;
				line++;
			} else i++;
			fprintf (output, "%6d", p+1);
			n++;
		};
	};
	if (n > 0) {
		fprintf(output, "\nTotal non-polymorphic = %d\n", n);
		fprintf(output, "\n");
		fflush(output);
	};
}
// --------------------------------------------------------------------------
// Jan 2017: add parameter specP
void PrtPop (FILE *output, int popRead, char *popID, int samp, char mLD,
				char mHet,	char mNomura, char mating,
				int nloci, int nMobil[], char *locUse, char specP)
{
	int n;
	if (output == NULL) return;
	char single = (mLD+mHet+mNomura > 0)? 1: 0;
	PrtMonoLoc (output, nloci, nMobil, locUse, popID, single);
	if (mLD > 0 && popRead == 1) {
		if (mating == 0) fprintf(output, "LD mating model: Random\n");
		else fprintf(output, "LD mating model: Monogamy\n");
		if (specP != 0) fprintf (output, "\n(Symbol \"%s\" in Frequency "
			"means that NO Singleton Alleles are accepted.)\n", NOSNGL);
	};
	if (single > 0) {
		fprintf (output, "\nPopulation%6d [%s]  (Number of Individuals = %d)\n",
			popRead, popID, samp);
		for (n=0; n<16; n++) fprintf (output, "*");
		fprintf (output, "\n");
	};
	fflush (output);
}

// --------------------------------------------------------------------------
void PrtFreq (FILE *output, char mLD, float *critVal, int nCrit,
				char spec1, char spec2)
// Dec 2016:
// called by PrtTemporal, RunPop0 (when either Het or LD is run)
// Since the special value for dropping singleton alleles is not used in
// Het or Temporal, we should add parameter mLD
{
	int n;
	int m = 26+12*nCrit;
	if (output == NULL) return;
	// print lines consisting of m = 26+12*nCrit characters spec1:
	PrtLines (output, m, spec1);
//	PrtDashes (output, nCrit, spec1, 0);
	fprintf (output, "Lowest Allele Frequency Used");
	for (n=0; n<nCrit; n++) {
		if (*(critVal+n) > 0 && *(critVal+n) <= PCRITX) {
			if (mLD > 0) fprintf (output, "%10s  ", NOSNGL);
			continue;
		}
		if (*(critVal+n) > 0) fprintf (output, "%10.3f  ", *(critVal+n));
		else fprintf (output, "%9s", "0+");
	};
	fprintf (output , "\n");
	// print lines consisting of m = 26+12*nCrit characters spec2:
	PrtLines (output, m, spec2);
//	PrtDashes (output, nCrit, spec2, 0);
	fflush (output);
}

// --------------------------------------------------------------------------

void PrtLDResults (FILE *output, int nCrit, float *wHarmonic,
				   double *nIndSum, float *rB2WAve, float *wExpR2,
				   float *estNe, float infinite, char bigInd)
{
	int n;
	float indMax = (float) MAXLONG;
	unsigned long indPrt;
	if (output == NULL) return;
	fprintf (output, "\nLINKAGE DISEQUILIBRIUM METHOD\n\n");
	fprintf (output, "Harmonic Mean Sample Size =");
	if (bigInd == 0) {
		fprintf (output, "%11.1f", *wHarmonic);
		for (n=1; n<nCrit; n++) fprintf (output, "%12.1f", *(wHarmonic+n));
	} else {
		fprintf (output, "%12.1f", *wHarmonic);
		for (n=1; n<nCrit; n++) fprintf (output, "%14.1f", *(wHarmonic+n));
	};
	fprintf (output, "\nIndependent Comparisons =");
	if (bigInd == 0) {
		if (indMax <= *nIndSum) indPrt = MAXLONG;
		else indPrt = (unsigned long) *nIndSum;
	// format right below is: "eleven ell u"
		fprintf (output, "%11lu", indPrt);
		for (n=1; n<nCrit; n++) {
			if (indMax <= *(nIndSum+n)) indPrt = MAXLONG;
			else indPrt = (unsigned long) *(nIndSum+n);
			fprintf (output, "%12lu", indPrt);
		};
	} else {
		if (indMax <= *nIndSum) indPrt = MAXLONG;
		else indPrt = (unsigned long) *nIndSum;
		fprintf (output, "%13lu", indPrt);
		for (n=1; n<nCrit; n++) {
			if (indMax <= *(nIndSum+n)) indPrt = MAXLONG;
			else indPrt = (unsigned long) *(nIndSum+n);
			fprintf (output, "%14lu", indPrt);
		};
	};
	fprintf (output, "\nOverAll r^2 =");
	if (bigInd == 0) {
		fprintf (output, "%25.6f", *rB2WAve);
		for (n=1; n<nCrit; n++) fprintf (output, "%12.6f", *(rB2WAve+n));
	} else {
		fprintf (output, "%27.6f", *rB2WAve);
		for (n=1; n<nCrit; n++) fprintf (output, "%14.6f", *(rB2WAve+n));
	};
	fprintf (output, "\nExpected r^2 Sample =");
	if (bigInd == 0) {
		fprintf (output, "%17.6f", *wExpR2);
		for (n=1; n<nCrit; n++) fprintf (output, "%12.6f", *(wExpR2+n));
	} else {
		fprintf (output, "%19.6f", *wExpR2);
		for (n=1; n<nCrit; n++) fprintf (output, "%14.6f", *(wExpR2+n));
	};
	fprintf (output, "\nEstimated Ne^ =");
	if (bigInd == 0) {
		if (*estNe >= infinite || *estNe < 0)
			fprintf (output, "%23s", "Infinite");
		else fprintf (output, "%23.1f", *estNe);
		for (n=1; n<nCrit; n++) {
			if (*(estNe+n) >= infinite || *(estNe+n) < 0)
				fprintf (output, "%12s", "Infinite");
			else fprintf (output, "%12.1f", *(estNe+n));
		};
	} else {
		if (*estNe >= infinite || *estNe < 0)
			fprintf (output, "%25s", "Infinite");
		else fprintf (output, "%25.1f", *estNe);
		for (n=1; n<nCrit; n++) {
			if (*(estNe+n) >= infinite || *(estNe+n) < 0)
				fprintf (output, "%14s", "Infinite");
			else fprintf (output, "%14.1f", *(estNe+n));
		};
	};
/*
	fprintf (output, "Harmonic Mean Sample Size = %10.1f", *wHarmonic);
	for (n=1; n<nCrit; n++) fprintf (output, "%12.1f", *(wHarmonic+n));
	fprintf (output, "\n");
	fprintf (output, "Independent Comparisons = %10lu", (long) *nIndSum);
	for (n=1; n<nCrit; n++) fprintf (output, "%12lu", (long) *(nIndSum+n));
	fprintf (output, "\n");
	fprintf (output, "OverAll r^2 = %24.5f", *rB2WAve);
	for (n=1; n<nCrit; n++) fprintf (output, "%12.5f", *(rB2WAve+n));
	fprintf (output, "\n");
	fprintf (output, "Expected r^2 Sample = %16.5f", *wExpR2);
	for (n=1; n<nCrit; n++) fprintf (output, "%12.5f", *(wExpR2+n));
	fprintf (output, "\n");
	if (*estNe >= infinite || *estNe < 0)
		fprintf (output, "Estimated Ne^ = %22s", "Infinite");
	else
		fprintf (output, "Estimated Ne^ = %22.1f", *estNe);
	for (n=1; n<nCrit; n++) {
		if (*(estNe+n) >= infinite || *(estNe+n) < 0)
			fprintf (output, "%12s", "Infinite");
		else
			fprintf (output, "%12.1f", *(estNe+n));
	};
*/
	fprintf (output, "\n\n");
	fflush (output);
}

// --------------------------------------------------------------------------

void PrtLDConfid (FILE *output, int nCrit, float confidLow[],
					float confidHi[], float infinite, int mode, char *header,
					char *jackOK, char bigInd)

{
	int n, k;
	if (output == NULL) return;
	if (*header >= 1) {	// print header the first time this function is called.
		fprintf (output, "95%% CIs for Ne^\n");
	};
	if (mode == 0)
		fprintf (output, "* Parametric              ");
	else {
	// in case no jackknife is calculated because *jackOK = 0:
		for (n=0, k=0; n<nCrit; n++) if (*(jackOK+n) == 1) k++;
		if (k > 0) {
//			fprintf (output, "* JackKnife on Loci       ");
			fprintf (output, "* JackKnife on Samples    ");
//                            123456789012345678901234567890
		} else {
			fprintf (output,
				"* CIs by Jackknife are skipped when number of polymorphic loci > %d", MAXJACKLD);
			fprintf (output, "\n\n");
			fflush (output);
			// reassign header, so that the next call, no "95% ..." printed
			*header = 0;
			return;
		};
	};
	k = 0;
	if (bigInd == 0) {
		for (n=0; n<nCrit; n++) {
			if (mode == 0 || *(jackOK+n) == 1) {
				if (confidLow[n] < infinite && confidLow[n] >= 0)
					fprintf (output, "%12.1f", confidLow[n]);
				else fprintf (output, "%12s", "Infinite");
			} else {
				fprintf (output, "%12s", "SKIPPED");
//				fprintf (output, "%12s", "UNAVAILABLE");
				k++;
			};
		};
		fprintf (output, "\n%26s", " ");
		for (n=0; n<nCrit; n++) {
			if (mode == 0 || *(jackOK+n) == 1) {
				if (confidHi[n] < infinite && confidHi[n] > 0)
					fprintf (output, "%12.1f", confidHi[n]);
				else fprintf (output, "%12s", "Infinite");
			};
		};
	} else {
		for (n=0; n<nCrit; n++) {
			if (mode == 0 || *(jackOK+n) == 1) {
				if (confidLow[n] < infinite && confidLow[n] >= 0)
					fprintf (output, "%14.1f", confidLow[n]);
				else fprintf (output, "%14s", "Infinite");
			} else {
				fprintf (output, "%14s", "SKIPPED");
//				fprintf (output, "%14s", "UNAVAILABLE");
				k++;
			};
		};
		fprintf (output, "\n%26s", " ");
		for (n=0; n<nCrit; n++) {
			if (mode == 0 || *(jackOK+n) == 1) {
				if (confidHi[n] < infinite && confidHi[n] > 0)
					fprintf (output, "%14.1f", confidHi[n]);
				else fprintf (output, "%14s", "Infinite");
			};
		};
	};
	if (k > 0) fprintf (output,
		"\n\n  CIs by Jackknife are skipped when number of polymorphic loci > %d", MAXJACKLD);
/*
	for (n=0; n<nCrit; n++) {
		if (mode == 0 || *(jackOK+n) == 1) {
			if (confidLow[n] < infinite && confidLow[n] >= 0)
				fprintf (output, "%12.1f", confidLow[n]);
			else fprintf (output, "%12s", "Infinite");
		} else fprintf (output, "%12s", "UNAVAILABLE");
	};
	fprintf (output, "\n%26s", " ");
	for (n=0; n<nCrit; n++) {
		if (mode == 0 || *(jackOK+n) == 1) {
			if (confidHi[n] < infinite && confidHi[n] >= 0)
				fprintf (output, "%12.1f", confidHi[n]);
			else fprintf (output, "%12s", "Infinite");
		};
	};
*/
	fprintf (output, "\n\n");
	fflush (output);
// reassign header, so that the next call, there is no "95% ..." printed
	*header = 0;

}

// --------------------------------------------------------------------------

void PrtPair (FILE *output, int num, char *id, int nChar, char pair) {
// If (pair=1), print to output a pair consisting of a number "num" and
// a string "id", separated by a colon, for a total length = nChar.
// Take the rightmost chars of "id" if it is too long, or fill in
// with blanks if it is too short.
// Print nChar blanks if pair /= 1.
	int i, j, k, n;
	k = strlen(id);
	for (j=0, n=num; n > 0; n /= 10, j++);	// j = # digits in num
	// Let i be the maximum characters in id that can be fit
	i = nChar - (j+1);
	if (pair == 1) {
		fprintf (output, "%d:", num);
		for (j=k-i; j < k || j < nChar; j++) {
			if (j >= 0 && j < k) fprintf (output, "%c", id[j]);
			else if (j < i && j >= k) fprintf (output, "%c", ' ');
		};
	} else for (j=0; j < nChar; j++) fprintf (output, "%c", ' ');
}

// --------------------------------------------------------------------------
// Jan 2017: adjust printing when spec Pcrit is present
void PrtLDxFile (char *inpName, FILE *xOutput, int samp, float *wHarmonic,
				  int popRead, int popStart, char *popID, float *critVal,
				  int nCrit, double *nIndSum, float *rB2WAve, float *wExpR2,
				  float *estNe, char param, char jacknife, float infinite,
				  float *confParalow, float *confParahi, float *confJacklow,
				  float *confJackhi, long *Jdegree, char *jackOK,
				  char mating, int topCrit,
// add May 2013 to be able to print input file name and loci as columns,
// and not to print header when necessary.
// count represents input file numbering, common = 1 if this xOutput is for
// all input files, with only one header to be printed.
// In this, count is the number of input file so far.
				  int nloci, int count, char common)
// print to extra ouput (shorter form) file.
{
	int i, k, m, n, lenInp;
	int stCrit, critOut;
	char pair;
	int dashes = 77;
	float indMax = (float) MAXLONG;
	unsigned long indPrt;
	if (xOutput == NULL) return;
// Jan 2017:
	char spec = 0;
// don't print jackknife if none is OK
// If this block is in effect, no header for jackknife CI if the first
// population doesn't have jackknife calculated (see PrtLDHeader). As a
// result, for later population, jackknife may be calculated and printed here,
// but without header saying what it is. So, we comment out
//	for (n=0, k=0; n<nCrit; n++) {
//		if (*(jackOK+n) == 1) k++;
//	}
//	if (k == 0) jacknife = 0;

// These are for testing, comment out when done
//topCrit = 2;
//param = 0;
//jacknife = 0;

// Critical values are indexed from big values to small, i.e., the first one,
// index=0, is the largest value. Parameter topCrit indicates that only the
// "topCrit" largest values will be outputted here.
// Numnber of critical values unchanged for output if topCrit < 0 or >= nCrit;
// otherwise, take the top "topCrit" critical values if topCrit > 0.
// when topCrit =0, take only critical value 0; in such case, index for
// critical value takes only one value, which is (nCrit - 1)
	stCrit = 0;
	if (topCrit > 0 && topCrit < nCrit) nCrit = topCrit;
	else if (topCrit == 0) stCrit = nCrit - 1;
	critOut = nCrit - stCrit;
	if (critOut > 1) dashes += 8;	// if more than one crit. values,
									// then they will be printed in each row
	if (common == 1) dashes += 27;	// more columns when this file is common
	k = param + jacknife;
	dashes += 20*k;
	if (jacknife == 1) dashes += 4;
	lenInp = 19;	// max width for printing Input name and number.
	if (popRead == popStart) {	// print header at the first pop read only
		if (common == 0 || count == 1) {	// print header
//			fprintf (xOutput, "\nLINKAGE DISEQUILIBRIUM METHOD, Mating Model: ");
			fprintf (xOutput, "\nMating Model: ");
			if (mating == 0) fprintf (xOutput, "Random");
			else fprintf (xOutput, "Monogamy");
			fprintf (xOutput, "\n\n");
			if (critOut > 1) fprintf (xOutput,
			"Lowest allele frequencies used, ordered in each population:\n");
			else fprintf (xOutput, "Lowest allele frequency used:");
			for (i=stCrit; i<nCrit; i++) {
// Jan 2017
				if (critVal[i] > 0 && critVal[i] <= PCRITX) {
					fprintf (xOutput, "%10s", NOSNGL);
					spec = 1;
				} else
					fprintf (xOutput, "%10.4f", critVal[i]);
			}
			fprintf (xOutput, "\n");
			if (spec != 0) fprintf (xOutput,
				"(\"%s\": No Singleton Allele is accepted.)\n", NOSNGL);
			if (common != 0) fprintf (xOutput,
				"Input Names are shown up to %d righmost characters.\n",
				lenInp - 2);
			fprintf (xOutput,
				"Up to 17 righmost characters can be shown for population names.\n");
			PrtLines (xOutput, dashes, '-');
			// if output is common for all input files, have columns for
			// input file name and #loci        123456789012345678901234567
			if (common != 0) fprintf (xOutput, "Input File Number   #Loci  ");

// if want the m at the end of this if (popRead == ...) to be increased
// by some j, add that j to the number after the first % in the format below
//			if (critOut > 1)  fprintf (xOutput, "Population Number%2s Samp"
			if (critOut > 1)  fprintf (xOutput, "Population #%2s Samp"
				"%2sCrit.%2sWeighted%6s#Indep.    r^2%5sExp(r^2)%7sNe^%9s",
				" ", " ", " ", " ", " ", " ", " ");
//			else fprintf (xOutput, "Population Number%2s Samp"
			else fprintf (xOutput, "Population #%2s Samp"
				"%2sWeighted%6s#Indep.   r^2%5sExp(r^2)%7sNe^%9s",
				" ", " ", " ", " ", " ", " ");
			if (k == 2)
				fprintf (xOutput, "%8sCIs for Ne^", " ");
			else if (k == 1)
				fprintf (xOutput, "CIs for Ne^");
			fprintf (xOutput, "\n");
			// if output is common for all input files, have columns for
			// input file name and #loci        123456789012345678901234567
			if (common != 0) fprintf (xOutput, "then :Name                 ");
//			if (common != 0) fprintf (xOutput, "followed by :Name          ");
// if want the m at the end of this if (popRead == ...) to be increased
// by some j, add that j to the number after the first % in the format below
//			if (critOut > 1) fprintf (xOutput, "followed by :Name%2s Size  Value%2s"
			if (critOut > 1) fprintf (xOutput, "then :Name   %2sSize  Value%2s"
				"H. Mean %6sAlleles%12sSample%18s", " ", " ", " ", " ", " ");
//			else fprintf (xOutput, "followed by :Name%2s Size%2s"
			else fprintf (xOutput, "then by :Name%2sSize%2s"
				"H. Mean %6sAlleles%12sSample%18s", " ", " ", " ", " ", " ");
			if (k == 2)
//				fprintf (xOutput, "   Parametric       Jackknife Loci");
				fprintf (xOutput, "  Parametric       Jackknife Samp  (Eff.df)");
			else if (param == 1)
				fprintf (xOutput, "  Parametric");
			else if (jacknife == 1)
//				fprintf (xOutput, "  Jacknife Loci");
				fprintf (xOutput, "Jacknife Samp  (Eff.df)");
			fprintf (xOutput, "\n");

			for (i=0; i<dashes; i++) fprintf (xOutput, "-");
			fprintf (xOutput, "\n");
		}
	}
//	m = 17;	// 17 is the length of string "Population Number". m will be the
	m = 12;	// 12 is the length of string "Population #". m will be the
			// max width for printing Population name and number.
	for (n=stCrit; n<nCrit; n++) {
		if (common != 0) {	// print common values for the whole input file,
		// the order and name of file, and number of loci
			pair = (n==stCrit && popRead == popStart)?1:0;
			// insert a blank line before printing output for this input file
//			if (pair == 1 && count > 1) fprintf (xOutput, "\n");
			PrtPair (xOutput, count, inpName, lenInp, pair);
//			if (pair == 1) fprintf (xOutput, "%6d  ", nloci);
//			else  fprintf (xOutput, "%8c", ' ');
		// comment out the previous 2 lines, which are for printing input
		// name and loci on the same first line of data.
		// The following are for putting name, loci# in a separate line
			if (pair == 1) {
				fprintf (xOutput, "%6d\n", nloci);
				PrtLines (xOutput, lenInp+6, '-');
				PrtPair (xOutput, count, inpName, lenInp, 0);
			}
			fprintf (xOutput, "%8c", ' ');
		}
		pair = (n==stCrit)?1:0;
		PrtPair (xOutput, popRead, popID, m, pair);
		// print Ind. Alleles, r^2, exp(r^2), Ne:
		if (indMax <= *(nIndSum+n)) indPrt = MAXLONG;
		else indPrt = (unsigned long) *(nIndSum+n);
		if (critOut > 1) {	// critival values are printed
			if (n == stCrit) fprintf (xOutput, "%6d", samp);
			else fprintf (xOutput, "%6s", " ");
// Jan 2017: deal with spec Pcrit:
			if (critVal[n] > 0 && critVal[n] <= PCRITX)
				fprintf (xOutput, "%8s", NOSNGL);
			else fprintf (xOutput, "%8.4f", critVal[n]);
			fprintf (xOutput, "%9.1f%12lu%10.6f%10.6f",
					*(wHarmonic+n), indPrt, *(rB2WAve+n), *(wExpR2+n));
		} else {	// no need to print critical value (only one)
//			fprintf (xOutput, "%6d%9.1f%14lu%10.6f%10.6f",
			fprintf (xOutput, "%6d%9.1f%12lu%10.6f%10.6f",
					samp, *(wHarmonic+n), indPrt, *(rB2WAve+n), *(wExpR2+n));
		}
		if (*(estNe+n) < infinite) fprintf (xOutput, "%11.1f", *(estNe+n));
		else fprintf (xOutput, "%11s", "Infinite");
		if (param == 1) {
			if (confParalow[n] < infinite && confParalow[n] >= 0)
				fprintf (xOutput, "%10.1f", confParalow[n]);
			else fprintf (xOutput, "%10s", "Infinite");
			if (confParahi[n] < infinite && confParahi[n] >= 0)
				fprintf (xOutput, "%10.1f", confParahi[n]);
			else fprintf (xOutput, "%10s", "Infinite");
		}
		if (jacknife == 1) {
			if (*(jackOK+n) == 1) {
				if (confJacklow[n] < infinite && confJacklow[n] >= 0)
					fprintf (xOutput, "%10.1f", confJacklow[n]);
				else fprintf (xOutput, "%10s", "Infinite");
				if (confJackhi[n] < infinite && confJackhi[n] >= 0)
					fprintf (xOutput, "%10.1f", confJackhi[n]);
				else fprintf (xOutput, "%10s", "Infinite");
// add in Aug 2016:
				fprintf (xOutput, "%10lu", Jdegree[n]);
			} else {
				fprintf (xOutput, "%10s", "*");
				fprintf (xOutput, "%10s", "*");
			}
		}
		fprintf (xOutput, "\n");
	}
	fflush (xOutput);
}

// --------------------------------------------------------------------------
// Dec 2016:
// Add parameter critVal to avoid printing when PCrit takes spec value for
// dropping singletons. Note that the printing of Het method is under the
// same headlines as for LD, so when Pcrit is spec. value, printing for
// Het under that Pcrit should be blanks.
void PrtHetNe (FILE *output, float *hetWAve, float *Ne,
			   float *loNe, float *hiNe, float *hSamp,
			   char param, int nCrit, float critVal[],
			   long *nIndH, float infinite)
{
	int i;
	if (output==NULL || nCrit <= 0) return;

	fprintf (output, "\nHETEROZYGOTE EXCESS METHOD\n\n");
	fprintf (output, "%-27s", "Harmonic Mean Sample Size =");
	for (i=0; i<nCrit; i++) {
// Dec 2016:
		if (critVal[i] > 0 && critVal[i] <= PCRITX) {
			fprintf (output, "%12s", " ");
			continue;
		}
		fprintf (output, "%11.1f ", *(hSamp+i));
	}
	fprintf (output, "\n");
	fprintf (output, "%-26s", "Independent Alleles =");
	for (i=0; i<nCrit; i++) {
// Dec 2016:
		if (critVal[i] > 0 && critVal[i] <= PCRITX) {
			fprintf (output, "%12s", " ");
			continue;
		}
		fprintf (output, "%10lu  ", *(nIndH+i));
	}
	fprintf (output, "\n");

	fprintf (output, "%-26s", "Weighted Mean D =");
	for (i=0; i<nCrit; i++) {
// Dec 2016:
		if (critVal[i] > 0 && critVal[i] <= PCRITX) {
			fprintf (output, "%12s", " ");
			continue;
		}
		fprintf (output, "%12.5f", *(hetWAve+i));
	}
// end of the for loop, n is now the number of CIs corresponding to
// positive values of D's (which are *(hetWAve+i) for i = 1, ..., nCrit).
// Actually, this is used in the part commented out below.
// That part of code is for the case that we only print CIs
// if the D-value hetWAve is positive.
	fprintf (output, "\n");
	fprintf (output, "%-26s", "Estimated Neb^  =");
	for (i=0; i<nCrit; i++) {
// Dec 2016:
		if (critVal[i] > 0 && critVal[i] <= PCRITX) {
			fprintf (output, "%12s", " ");
			continue;
		}
		if (*(Ne+i) <= 0 || *(Ne+i) == infinite)
			fprintf (output, "%12s", "Infinite");
//			fprintf (output, "%12s", "UNDEFINED");
		else
			fprintf (output, "%12.1f", *(Ne+i));
	}
	fprintf (output, "\n\n");
//	if (param+jack == 0) return;
	if (param == 0) return;
	fprintf (output, "95%% CIs for Nb:\n");

// temporarily block printing CI:
//	fprintf (output, "%26s", " ");
//	for (i=0; i<nCrit; i++) fprintf (output, "%12s", "PENDING");
//	fprintf (output, "\n");
//	return;


	if (param == 1) {
		fprintf (output, "%-26s", "* Parametric ");
		for (i=0; i<nCrit; i++) {
// Dec 2016:
			if (critVal[i] > 0 && critVal[i] <= PCRITX) {
				fprintf (output, "%12s", " ");
				continue;
			}
			if (*(loNe+i) > 0 && *(loNe+i) < infinite)
				fprintf (output, "%12.1f", *(loNe+i));
			else
				fprintf (output, "%12s", "Infinite");
		};
		fprintf (output, "\n");
		fprintf (output, "%26s", " ");
		for (i=0; i<nCrit; i++) {
// Dec 2016:
			if (critVal[i] > 0 && critVal[i] <= PCRITX) {
				fprintf (output, "%12s", " ");
				continue;
			}
			if (*(hiNe+i) > 0 && *(hiNe+i) < infinite)
				fprintf (output, "%12.1f", *(hiNe+i));
			else
				fprintf (output, "%12s", "Infinite");
		};
		fprintf (output, "\n");
	};
// block jackknife:
/*
	if (jack == 1) {
		fprintf (output, "%-26s", "* Jackknife on Loci");
		for (i=0; i<nCrit; i++) {
			if (*(jloNe+i) > 0)
				fprintf (output, "%12.1f", *(jloNe+i));
			else
				fprintf (output, "%12s", "Infinite");
		};
		fprintf (output, "\n");
		fprintf (output, "%26s", " ");
		for (i=0; i<nCrit; i++) {
			if (*(jhiNe+i) > 0 && *(jhiNe+i) < infinite)
				fprintf (output, "%12.1f", *(jhiNe+i));
			else
				fprintf (output, "%12s", "Infinite");
		};
		fprintf (output, "\n");
	};
*/
}
// --------------------------------------------------------------------------
//
// Function CritEndLine sets the starting and ending index for array critVal,
// for outputting results corresponding to critVal[] (Pcrits) in xFiles.
// This function is called by PrtHetxFile, PrtTempxFile,
// and by PrtLeading (called by Prt**Header's used in Prt**TabFile's),
// so needs to be on top of PrtHetxFile, PrtTempxFile, and PrtLeading
// --------------------------------------------------------------------------
// There is a spec Pcrit that is only applied to LD method.
// The param topCrit represents the top Pcrits to  output in these xFiles,
// running from the highest down. If those Pcrits include this spec Pcrit,
// then the actual number of Pcrits to output will be one less if the
// xFile is not for LD.
// Practical outputs of this function: stCrit = starting index,
// critOut = number of Pcrits to output,
// ending index = return value of the function, minus 1
// Those will be used in the "for" loop at the function calling this function
// --------------------------------------------------------------------------
// When this function is called from non-LD methods, a "for" loop for running
// through necessary Pcrits will be run after calling this function:
// If this spec Pcrit is present, then the index for this Pcrit is either
// the next-to-last index in the Pcrit list if the list contains 0, or the
// last index if the list does not contain 0. In the loop, if index n is the
// index for this spec Pcrit, we just need to skip it. The only problem is
// when the loop has only one index which is the index for this spec Pcrit.
// This happens when
// (a) topCrit = 1, and the spec Pcrit is the top (so 0 is the only remaining,
//     since there must be another Pcrit). In this instance, the Pcrit
//     list has only 2 members, the spec Pcrit corresponds to index n = 0,
//     the Pcrit=0 corresponds to index n = 1.
//     We want output results for Pcrit = 0 whose index is 1, so the loop will
//     run through all Pcrits (then skip at index = 0 for this spec Pcrit).
// (b) topCrit = 0, so it was then reassigned so that the only Pcrit outputted
//     is the smallest Pcrit, and when there is no Pcrit=0, the smallest Pcrit
//     turns out to be the spec Pcrit, corresponding to the last index. Then
//     we need to output results for the Pcrit preceding this spec Pcrit.

int CritEndLine (float *critVal, int nCrit, int topCrit, int *stCrit,
			int *critOut, char forLD)
{
	int i;
	*stCrit = 0;	// only reassigned if topCrit = 0, then also *critOut = 1
	if (topCrit < nCrit) {
		if (topCrit > 0) {
		// When all three: topCrit = 1, critVal[0] = 0 or the spec Pcrit,
		// and forLD = 0: no reassignment for nCrit. If critVal[0] = 0, then
		// only one Pcrit, so nCrit = 1, no spec PCrit. Otherwise, nCrit = 2.
		// Then "for" loop for PCrits at the calling function has 2 indices:
		// starting at stCrit=0 and ending at index=1, critVal[0] is skipped.
		// When one of them is false, nCrit is reassigned below, the "for"
		// loop may or may not contain spec Pcrit. The spec Pcrit will be
		// skipped at that "for" loop if forLD = 0.
		// (a) forLD != 0: this function is called by LD, spec Pcrit is OK.
		// (b) topCrit > 1, reassigning nCrit = topCrit implies the "for"
		//     loop has at least 2 indices; so Pcrits for output exist if
		//     the spec Pcrit is skipped.
		// (c) critVal[0] > PCRITX: critVal[0] is not the spec Pcrit.
		//     Then it is accepted at Het or Temporal.
		//     (If critVal[0] = 0: only one Pcrit = 0.)
		// So we are assured that the "for" loop will go through with at least
		// one accepted Pcrit for the method.
			if (topCrit > 1 || critVal[0] > PCRITX || forLD != 0)
				nCrit = topCrit;
		} else {	// i.e., topCrit = 0
			if (critVal[nCrit-1] > 0 && critVal[nCrit-1] <= PCRITX && forLD == 0)
			// No Pcrit=0, the smallest Pcrit is the spec Pcrit, method is
			// not LD, there must be another Pcrit preceding the spec Pcrit,
			// which is then the smallest Pcrit for that method.
				*stCrit = nCrit - 2;
			else	// the smallest Pcrit is not the spec Pcrit or method is LD,
			// so the smallest Pcrit is accepted for the method that calls this
				*stCrit = nCrit - 1;
		}
	}
	*critOut = nCrit - *stCrit;
	// We have *critOut = 2 in the "if" case at the case "topCrit = 0" above,
	// then it will be reassigned below so that it is 1.
	if (forLD !=0) return nCrit;
	// the next loop is for non-LD method, which adjusts *critOut
	for (i=*stCrit; i<nCrit; i++) {
		if (critVal[i] > 0 && critVal[i] <= PCRITX) {
			(*critOut)--;
			break;
		}
	}
	// If *critOut > 1, then we must have *stCrit = 0
	return nCrit;
}

// --------------------------------------------------------------------------
// Modified Dec 2016:
// Recalculate stCrit, nCrit and critOut when the spec Pcrit is present
// in the PCrit array. Then adjust the "for" loop.

void PrtHetxFile (char *inpName, FILE *xOutput, int popRead, int popStart,
				  char *popID, float *critVal, int nCrit, long *indAlleH,
				  float *hetD, float *estHetN, char param, float infinite,
				  float *loHetNe, float *hiHetNe, int samp, float *hSamp,
				  int topCrit,
// add May 2013 to be able to print input file name and loci as columns,
// and not to print header when necessary.
// count represents input file numbering, common = 1 if this xOutput is for
// all input files, with only one header to be printed.
// In this, count is the number of input so far.
				  int nloci, int count, char common)
// print to extra ouput (shorter form) file.
{
	int i, m, n, lenInp;
	int stCrit, critOut;
	char pair;
	int dashes = 62;
	if (xOutput == NULL) return;
//topCrit = 0;
//param = 0;
	lenInp = 19;	// max width for printing Input name and number.

// Critical values are indexed from big values to small, i.e., the first one,
// index = 0, is the largest value. Parameter topCrit indicates that only the
// "topCrit" largest values will be outputted here.
//
// Number of critical values unchanged for output if topCrit < 0 or >= nCrit;
// otherwise, take the top "topCrit" critical values accepted for this method
// if topCrit > 0.
// In the "for" loop for Pcrits, index n runs from 0 to nCrit-1. Therefore,
// according to the previous sentence, nCrit will be reasigned = topCrit when
// topCrit > 0 and topCrit < nCrit (no change if topCrit < 0 or >= nCrit).
//
// When topCrit = 0, take only the smallest Pcrit; in such case, the loop
// over critical values contains only the index for this Pcrit. Hence, the
// starting index for n will be (nCrit-1), same as the ending index.
// For this, we use stCrit as the starting index, which is nCrit-1 in this
// case, and 0 for previous case.

// Jan 2017: Use function CritEndLine to set index limits
	nCrit = CritEndLine (critVal, nCrit, topCrit, &stCrit, &critOut, 0);

	if (critOut > 1) dashes += 7;	// more chars for printing crit. values
	if (param == 1) dashes += 20;
	if (common == 1) dashes += 27;	// more columns when this file is common
	if (popRead == popStart) {	// print header at the first pop read only
		if (common == 0 || count == 1) {	// print header
			fprintf (xOutput, "\nHETEROZYGOTE-EXCESS METHOD\n\n");
			if (critOut > 1) fprintf (xOutput, "Lowest allele frequencies used, "
								"ordered in each population:\n");
			else fprintf (xOutput, "Lowest allele frequency used:");
			for (i=stCrit; i<nCrit; i++) {
// Modify Jan 2017: skip when critVal[i] is the spec Pcrit
				if (critVal[i] > 0 && critVal[i] <= PCRITX) continue;
				fprintf (xOutput, "%10.4f", critVal[i]);
			}
			fprintf (xOutput, "\n");
			if (common != 0) fprintf (xOutput,
				"Input Names are shown up to %d righmost characters.\n",
				lenInp - 2);
			fprintf (xOutput,
				"Up to 17 righmost characters can be shown for population names.\n");
			PrtLines (xOutput, dashes, '-');
			// if output is common for all input files, have columns for
			// input file name and #loci        123456789012345678901234567
			if (common != 0) fprintf (xOutput, "Input File Number   #Loci  ");

// if want the m at the end of this if (popRead == ...) to be increased
// by some j, add that j to the number after the first % in the format below
			if (critOut > 1) fprintf (xOutput, "Population Number%2sSamp  Crit."
				"  Harmonic   #Indep.     D%10sNe%10s", " ", " ", " ");
			else fprintf (xOutput, "Population Number%2sSamp  "
				"Harmonic   #Indep.     D%10sNe%10s", " ", " ", " ");

			if (param == 1) fprintf (xOutput, "CIs for Ne");
			fprintf (xOutput, "\n");
			// if output is common for all input files, have columns for
			// input file name and #loci        123456789012345678901234567
			if (common != 0) fprintf (xOutput, "followed by :Name          ");
// if want the m at the end of this if (popRead == ...) to be increased
// by some j, add that j to the number after the first % in the format below
			if (critOut > 1) fprintf (xOutput, "followed by :Name%2sSize  Value  "
							"Mean Size%2sAlleles%26s", " ", " ", " ");
			else fprintf (xOutput, "followed by :Name%2sSize  "
							"Mean Size%2sAlleles%26s", " ", " ", " ");
			if (param == 1) fprintf (xOutput, "  Parametric");
			fprintf (xOutput, "\n");

			for (i=0; i<dashes; i++) fprintf (xOutput, "-");
			fprintf (xOutput, "\n");
		};
	};
	m = 17;	// 17 is the length of "Population Number". m will be the
			// max width for printing Population name and number.
// Jan 2017: add stCrit0 as the truly starting index since Pcrit at index
// stCrit may be the spec Pcrit, which is skipped.
// Then stCrit inside the loop is replaced by stCrit0.
	int stCrit0 = stCrit;
	for (n=stCrit; n<nCrit; n++) {
// Modify Jan 2017: skip when critVal[i] is the spec Pcrit
		if (critVal[n] > 0 && critVal[n] <= PCRITX) {
			if (n == stCrit) stCrit0++;
			continue;
		}
		if (common != 0) {	// print common values for the whole input file,
		// the order and name of file, and number of loci
			pair = (n==stCrit0 && popRead == popStart)?1:0;
			// insert a blank line before printing output for this input file
//			if (pair == 1 && count > 1) fprintf (xOutput, "\n");
			PrtPair (xOutput, count, inpName, lenInp, pair);
//			if (pair == 1) fprintf (xOutput, "%6d  ", nloci);
//			else  fprintf (xOutput, "%8c", ' ');
		// comment out the previous 2 lines, which are for printing input
		// name and loci on the same first line of data.
		// The following are for putting name, loci# in a separate line
			if (pair == 1) {
				fprintf (xOutput, "%6d\n", nloci);
				PrtLines (xOutput, lenInp+6, '-');
				PrtPair (xOutput, count, inpName, lenInp, 0);
			};
			fprintf (xOutput, "%8c", ' ');
		};
		pair = (n==stCrit0)?1:0;
		PrtPair (xOutput, popRead, popID, m, pair);
		// print samp size, harmonic samp size, Ind. Alleles, D measure, Ne:
		if (critOut > 1) {
			if (n==stCrit0) fprintf (xOutput, "%6d%8.4f%8.1f%9lu%10.5f", samp,
							critVal[n], *(hSamp+n), *(indAlleH+n), *(hetD+n));
			else fprintf (xOutput, "%6s%8.4f%8.1f%9lu%10.5f", " ",
							critVal[n], *(hSamp+n), *(indAlleH+n), *(hetD+n));
		} else fprintf (xOutput, "%6d%9.1f%9lu%10.5f", samp, *(hSamp+n),
				*(indAlleH+n), *(hetD+n));
		if (*(estHetN+n) >= infinite || *(estHetN+n) < 0)
			fprintf (xOutput, "%11s", "Infinite");
		else fprintf (xOutput, "%11.1f", *(estHetN+n));
		if (param == 1) {
			if (loHetNe[n] < infinite)
				fprintf (xOutput, "%10.1f", loHetNe[n]);
			else fprintf (xOutput, "%10s", "Infinite");
			if (hiHetNe[n] < infinite)
				fprintf (xOutput, "%10.1f", hiHetNe[n]);
			else fprintf (xOutput, "%10s", "Infinite");
		};
		fprintf (xOutput, "\n");
	};
	fflush (xOutput);
}


// --------------------------------------------------------------------------

void PrtNomuraNe (FILE *output, float f1, float Ne, int nCrit, int lastCrit,
				  float loNbCoan, float hiNbCoan, char jack, float hSamp)
// Add May 1, 2013 parameter lastCrit:
// On input, lastCrit = 0 when there are methods run before Nomura, and the
// last critical value is 0; otherwise that value was assigned non-zero (= 1).
// If this lastCrit is non-zero, then 0 is not labeled as a critical value,
// so we should emphasize outputs for Coancestry method are for crit value 0!
{
	int i, n;
	// this m is the same as in PrtFreq
	int m = 26+12*nCrit;
	if (output==NULL) return;
	n = (nCrit-1) * 12;	// # of blanks to skip
	if (lastCrit != 0) PrtLines (output, m, '-');
	fprintf (output, "\nMOLECULAR COANCESTRY METHOD");
	if (lastCrit != 0) fprintf (output, "\n(No frequency restriction)");
	fprintf (output, "\n\n");
	fprintf (output, "Harmonic Mean Sample Size =");
	for (i=0; i<n; i++) fprintf (output, " ");
	fprintf (output, "%11.1f\n", hSamp);
	fprintf (output, "OverAll f1^    = ");
	for (i=0; i<n; i++) fprintf (output, " ");
	fprintf (output, "%21.5f\n", f1);
	fprintf (output, "Estimated Neb^ = ");
	for (i=0; i<n; i++) fprintf (output, " ");
	if (Ne < 0 || Ne == INFINITE)
		fprintf (output, "%21s\n\n", "Infinite");
//		fprintf (output, "%21s\n\n", "UNDEFINED");
	else
		fprintf (output, "%21.1f\n\n", Ne);
	if (jack == 1) {
		fprintf (output, "95%% CIs for Ne^\n* Jackknife on Loci");

// temporarily block printing CI:
//		for (i=0; i<n; i++) fprintf (output, " ");
//		fprintf (output, "%19s\n", "PENDING");
//		return;

// the string "* Jackknife on Loci" occupies 19 characters,
// the last character for the result of the first critical value
// is at column 38: 39-19=19 should be the length reserved for CIs outputs
		for (i=0; i<n; i++) fprintf (output, " ");
		if (loNbCoan >= 0 && loNbCoan < INFINITE)
			fprintf (output, "%19.1f\n", loNbCoan);
		else
			fprintf (output, "%19s\n", "Infinite");
// add 19 to n, to replace the length of the string "* Jackknife on Loci"
		for (i=0; i<n+19; i++) fprintf (output, " ");
		if (hiNbCoan >= 0 && hiNbCoan < INFINITE)
			fprintf (output, "%19.1f\n", hiNbCoan);
		else
			fprintf (output, "%19s\n", "Infinite");

//		fprintf (output, "\n");
	};
// print to console (moved to the calling function RunPop#):
//	printf ("     Molecular Coancestry Method\n");
//	printf ("       Estimated Neb^:%21.1f\n", Ne);
}

// --------------------------------------------------------------------------

void PrtCoanxFile (char *inpName, FILE *xOutput, int popRead, int popStart,
				  char *popID, float f1, float coanNeb, char jacknife,
				  float infinite, float loNbCoan, float hiNbCoan, int samp,
				  float hSamCoan,
// add May 2013 to be able to print input file name and loci as columns,
// and not to print header when necessary.
// count represents input file numbering, common = 1 if this xOutput is for
// all input files, with only one header to be printed.
// In this, count is the number of input so far.
				  int nloci, int count, char common)
// print to extra ouput (shorter form) file.
{
	int i, m, lenInp;
	int dashes = 53;
	char pair;
	if (xOutput == NULL) return;
	lenInp = 19;	// max width for printing Input name and number.
	if (jacknife == 1) dashes += 20;
	if (common == 1) dashes += 27;	// more columns when this file is common
	if (popRead == popStart) {	// print header at the first pop read only
		if (common == 0 || count == 1) {	// print header
			fprintf (xOutput, "\nMOLECULAR COANCESTRY METHOD\n\n");
			if (common != 0) fprintf (xOutput,
				"Input Names are shown up to %d righmost characters.\n",
				lenInp - 2);
			fprintf (xOutput,
				"Up to 17 righmost characters can be shown for population names.\n");
			PrtLines (xOutput, dashes, '-');
			// if output is common for all input files, have columns for
			// input file name and #loci        123456789012345678901234567
			if (common != 0) fprintf (xOutput, "Input File Number   #Loci  ");

// if want the m at the end of this if (popRead == ...) to be increased
// by some j, add that j to the number after the first % in the format below
			fprintf (xOutput, "Population Number%2sSamp  Harmonic     f1^  "
							"%6sNeb^ %8s", " ", " ", " ");
			if (jacknife == 1) fprintf (xOutput, "CIs for Ne");
			fprintf (xOutput, "\n");
			// if output is common for all input files, have columns for
			// input file name and #loci        123456789012345678901234567
			if (common != 0) fprintf (xOutput, "followed by :Name          ");
// if want the m at the end of this if (popRead == ...) to be increased
// by some j, add that j to the number after the first % in the format below
			fprintf (xOutput, "followed by :Name%2sSize    Mean %26s",
					" ", " ");
			if (jacknife == 1) fprintf (xOutput, "    Jackknife");
			fprintf (xOutput, "\n");

			for (i=0; i<dashes; i++) fprintf (xOutput, "-");
			fprintf (xOutput, "\n");
		};
	};
	m = 17;	// 17 is the length of "Population Number". m will be the
			// max width for printing Population name and number.

	if (common != 0) {	// print common values for the whole input file,
		// the order and name of file, and number of loci
		pair = (popRead == popStart)?1:0;
		// insert a blank line before printing output for this input file
//		if (pair == 1 && count > 1) fprintf (xOutput, "\n");
		PrtPair (xOutput, count, inpName, lenInp, pair);
//		if (pair == 1) fprintf (xOutput, "%6d  ", nloci);
//		else  fprintf (xOutput, "%8c", ' ');
		// comment out the previous 2 lines, which are for printing input
		// name and loci on the same first line of data.
		// The following are for putting name, loci# in a separate line
		if (pair == 1) {
			fprintf (xOutput, "%6d\n", nloci);
			PrtLines (xOutput, lenInp+6, '-');
			PrtPair (xOutput, count, inpName, lenInp, 0);
		};
		fprintf (xOutput, "%8c", ' ');
	};
	pair = 1;
	PrtPair (xOutput, popRead, popID, m, pair);
	fprintf (xOutput, "%6d%9.1f%10.5f", samp, hSamCoan, f1);
	if (coanNeb >= infinite || coanNeb < 0)
		fprintf (xOutput, "%11s", "Infinite");
	else fprintf (xOutput, "%11.1f", coanNeb);
	if (jacknife == 1) {
		if (loNbCoan < infinite)
			fprintf (xOutput, "%10.1f", loNbCoan);
		else fprintf (xOutput, "%10s", "Infinite");
		if (hiNbCoan < infinite)
			fprintf (xOutput, "%10.1f", hiNbCoan);
		else fprintf (xOutput, "%10s", "Infinite");
	};
	fprintf (xOutput, "\n");
	fflush (xOutput);
}

// --------------------------------------------------------------------------
// Modified Dec 2016: add parameter critVal

void PrtTempVal (FILE *output, int nCrit, float critVal[],
				float *Hmean, float *fmean, float *fprimeMean,
				float *Ne, float *loNe, float *hiNe,
				float *jloNe, float *jhiNe, char param, char jack,
				char label1[], char label2[], char allTemp, float infinite)
// Dec 2016: Since Temporal method will skip for spec value of cutoff
// (which reject only singleton alleles), the part corresponding to that
// value is unassigned, and will be skipped in the printing process

{
	int n;
	if (output == NULL) return;
// Dec 2016:
	char skip;
	int nCrit0 = nCrit;
	for (n=0; n<nCrit; n++) {
		if (critVal[n] > 0 && critVal[n] <= PCRITX) {
			skip = 1;
			break;
		}
	}
	// nCrit0 is the real number of Pcrits, used only for counting dashes
	nCrit0 -= skip;

	if (allTemp > 1) {
		fprintf (output, "   ");
// print lines of 26+12*nCrit characters '-'
// Dec 2016: change nCrit to nCrit0:
		PrtLines (output, 26+12*nCrit0-3, '-');
//		PrtDashes (output, nCrit, '-', 3);
	};

	if (allTemp > 1) fprintf (output, "%s\n", label1);
	fprintf (output, "%s", "   Harmonic Mean Samp Size =");
	for (n=0; n<nCrit; n++) {
// Dec 2016
		if (critVal[n] > 0 && critVal[n] <= PCRITX) {
//			fprintf (output, "%12s", " ");
			continue;
		}
//		fprintf (output, "%10.1f  ", *(Hmean+n));
// want no blank at the end, so replacing the above code by:
		if (n == 0) fprintf (output, "%10.1f", *(Hmean+n));
		else fprintf (output, "%12.1f", *(Hmean+n));
	}
	fprintf (output, "\n");
//	fprintf (output, "%s", "   Independent Alleles =");
//	for (n=0; n<nCrit; n++) fprintf (output, "%12lu", *(nAlle+n));
//	fprintf (output, "\n");
	fprintf (output, "%19s =     ", label2);
	for (n=0; n<nCrit; n++) {
// Dec 2016
		if (critVal[n] > 0 && critVal[n] <= PCRITX) {
//			fprintf (output, "%12s", " ");
			continue;
		}
		fprintf (output, "%12.5f", *(fmean+n));
	}
	fprintf (output, "\n");
	fprintf (output, "%19s =     ", "F'");
	for (n=0; n<nCrit; n++) {
// Dec 2016
		if (critVal[n] > 0 && critVal[n] <= PCRITX) {
//			fprintf (output, "%12s", " ");
			continue;
		}
		fprintf (output, "%12.5f", *(fprimeMean+n));
	}
	fprintf (output, "\n");
	fprintf (output, "%21s     ", "* Ne =");
	for (n=0; n<nCrit; n++) {
// Dec 2016
		if (critVal[n] > 0 && critVal[n] <= PCRITX) {
//			fprintf (output, "%12s", " ");
			continue;
		}
		if (*(Ne+n) < infinite)
			fprintf (output, "%12.1f", *(Ne+n));
		else
			fprintf (output, "%12s", "Infinite");
	};
	fprintf (output, "\n");
	if (param+jack == 0) return;
	fprintf (output, "\n     95%% CIs for Ne:\n");
	if (param == 1) {
		fprintf (output, "%-26s", "     * Parametric ");
		for (n=0; n<nCrit; n++) {
// Dec 2016
			if (critVal[n] > 0 && critVal[n] <= PCRITX) {
//				fprintf (output, "%12s", " ");
				continue;
			}
			if (*(loNe+n) < infinite)
				fprintf (output, "%12.1f", *(loNe+n));
			else
				fprintf (output, "%12s", "Infinite");
		};
		fprintf (output, "\n%26s", " ");
		for (n=0; n<nCrit; n++) {
// Dec 2016
			if (critVal[n] > 0 && critVal[n] <= PCRITX) {
//				fprintf (output, "%12s", " ");
				continue;
			}
			if (*(hiNe+n) < infinite)
				fprintf (output, "%12.1f", *(hiNe+n));
			else
				fprintf (output, "%12s", "Infinite");
		};
		fprintf (output, "\n");
	};
	if (jack == 1) {
		fprintf (output, "%-26s", "     * Jackknife on Loci");
		for (n=0; n<nCrit; n++) {
// Dec 2016
			if (critVal[n] > 0 && critVal[n] <= PCRITX) {
//				fprintf (output, "%12s", " ");
				continue;
			}
			if (*(jloNe+n) < infinite)
				fprintf (output, "%12.1f", *(jloNe+n));
			else
				fprintf (output, "%12s", "Infinite");
		};
		fprintf (output, "\n%26s", " ");
		for (n=0; n<nCrit; n++) {
// Dec 2016
			if (critVal[n] > 0 && critVal[n] <= PCRITX) {
//				fprintf (output, "%12s", " ");
				continue;
			}
			if (*(jhiNe+n) < infinite)
				fprintf (output, "%12.1f", *(jhiNe+n));
			else
				fprintf (output, "%12s", "Infinite");
		};
		fprintf (output, "\n");
	};
// The commented out part is for printing CI only when Ne not negative
/*
	m = 0;	// count how many CIs corresponding to critical values
	for (n=0; n<nCrit; n++) {
		if (*(Ne+n) < 0)
			fprintf (output, "%12s", "NONE");
		else {
			fprintf (output, "%12.1f", *(loNe+n));
			m++;
		};
	};
	fprintf (output, "\n");
	if (m > 0) {
		fprintf (output, "%26s", " ");
		for (n=0; n<nCrit; n++) {
			if (*(Ne+n) < 0)
				fprintf (output, "%12s", " ");
			else {
				if (*(hiNe+n) < INFINITE)
					fprintf (output, "%12.1f", *(hiNe+n));
				else
					fprintf (output, "%12s", "Infinite");
			};
		};
		fprintf (output, "\n");
	};
*/
}

// --------------------------------------------------------------------------
// Modified Dec 2016
// At the "for" loop at Critical values, skip when the value is PCRITX
// As PrtFreq is modified with added parameter mLD, "0" will be in place
// to exclude the spec Pcrit when that function is called here
void PrtTemporal (FILE *output, int nCrit, float critVal[],
					int g1, int g2, int nGeneration, float timeline[],
					long *nTotAlle, long *nIndAlle,
					float *Hkmean, float *Hcmean, float *Hsmean,
					float *fkmean, float *fcmean, float *fsmean,
					float *fkprimeMean, float *fcprimeMean, float *fsprimeAll,
					float *NeTempk, float *NeTempc, float *NeTemps,
					float *loNek, float *loNec, float *loNes,
					float *hiNek, float *hiNec, float *hiNes,
					float *jloNek, float *jloNec, float *jloNes,
					float *jhiNek, float *jhiNec, float *jhiNes,
					char param, char jack,
					float infinite, char tempk, char tempc, char temps,
					char *popID[], int *popSize, int census)
{
	int n, nc, m;
	char allTemp = tempk + tempc + temps;
	if (output == NULL) return;
	if (allTemp == 0) return;
	if (g1 + g2 == 1) {
//		printf ("\nTemporal Method ...\n");
		fprintf (output, "\nPopulation: ");
		for (m=0; m<nGeneration; ) {
			fprintf (output, "%s", popID [m]);
			if (++m < nGeneration) fprintf (output, "/");
			else fprintf (output, ",\t%d Samples,\tSample Sizes: ", nGeneration);
//			else fprintf (output, ", \tSample Sizes: ");
		};
		for (m=0; m<nGeneration; ) {
			fprintf (output, "%d", *(popSize+m));
			if (++m < nGeneration) fprintf (output, "/");
			else fprintf (output, "\n");
		};
		m = 0;
		fprintf (output, "\nTEMPORAL METHOD");
		if (allTemp > 1) fprintf (output, "S");
		fprintf (output, " (");
		if (tempk == 1) {
			fprintf (output, "Pollak");
			m++;
		};
		if (tempc ==1) {
			if (m > 0) fprintf (output, ", ");
			fprintf (output, "Nei/Tajima");
			m++;
		};
		if (temps ==1) {
			if (m > 0) fprintf (output, ", ");
			fprintf (output, "Jorde/Ryman");
		};
		fprintf (output, "), ");
		if (census > 0) {
			fprintf (output, "Plan I,");
			if (allTemp > 1) fprintf (output, "\n%16s", " ");
			fprintf (output, " Census Population Size = %d.", census);
		} else fprintf (output, "Plan II.");
		fprintf (output, "\n");
//		fprintf (output, " From %d samples\n", nGeneration);
// Dec 2016: the 0 is for ignoring spec falue for excluding singleton alleles
		PrtFreq (output, 0, critVal, nCrit, '-', '=');
	};
	// 11 characters will be used for writing sample # and sample ID
	// including the two square brackets.
	// Sample # uses up to 2 characters, and sample ID uses up to 7.
	// We want to right align the combination, that is, fill in the left
	// by blanks so that the whole will occupy 11 characters. Since we
	// have to deal with both together without having any blanks between
	// them, we should find the correct number of blanks to fill.
	nc = 11 - 4;
	m = nc - strlen (popID[g1]);		// this is the number of blanks to be
			// added or the number of chars to be taken out from popID[g1]
	fprintf (output, "Samples ");
	for (n = m; n > 0; n--) fprintf (output, "%c", ' ');
	m = (m > 0)? 0: -m;
	n = nc - strlen (popID[g2]);
	n = (n > 0)? 0: -n;
	fprintf (output, "%2d[%s] & %d[%s]\n", g1+1, popID[g1]+m,
										g2+1, popID[g2]+n);
	// Note: the word Generations has 4 chars more than Samples.
	// Thus, number 7 in the first format below should be the same as nc.
	fprintf (output, "Generations %7.1f & %-7.1f\n", timeline[g1], timeline[g2]);
	fprintf (output, "\n%s", "   Independent Alleles =");

	// Dec 2016:
	char skip = 0;
	for (n=0; n<nCrit; n++) {
		if (critVal[n] > 0 && critVal[n] <= PCRITX) {
			skip = 1;
			continue;
		}
		fprintf (output, "%12lu", *(nIndAlle+n));
	}
	fprintf (output, "\n");
	if (tempk == 1) {
		PrtTempVal (output, nCrit, critVal, Hkmean, fkmean, fkprimeMean,
					NeTempk, loNek, hiNek, jloNek, jhiNek, param, jack,
					"  (Pollak)", "Fk", allTemp, infinite);
	};
	if (tempc == 1) {
		PrtTempVal (output, nCrit, critVal, Hcmean, fcmean, fcprimeMean,
					NeTempc, loNec, hiNec, jloNec, jhiNec, param, jack,
					"  (Nei/Tajima)", "Fc", allTemp, infinite);
	};
	if (temps == 1) {
		PrtTempVal (output, nCrit, critVal, Hsmean, fsmean, fsprimeAll,
					NeTemps, loNes, hiNes, jloNes, jhiNes, param, jack,
					"  (Jorde/Ryman)", "Fs", allTemp, infinite);
	};
// Dec 2016:
	nCrit -= skip;
	// print lines of 26+12*nCrit characters '='
	PrtLines (output, 26+12*nCrit, '=');
//	PrtDashes (output, nCrit, '=', 0);

}

// --------------------------------------------------------------------------
// Dec 2016: critOut is stated as parameter, instead of local variable,
// the function is called from PrtTempxFile, where critOut is modified.

void PrtTempValx (FILE *xOutput, int nCrit, float *critVal, long *nAlle,
				float *Hmean, float *fmean, float *fprimeMean,
				float *Ne, float *loNe, float *hiNe,
				float *jloNe, float *jhiNe, char param, char jack,
				char label[], char popPair[], float time1, float time2,
				int nSkip1, int nSkip2, char nTemp, float infinite,
				int stCrit, int critOut, int nPlan, int census,
				char common, int k, char *called)
// when this is called, parameter called will be set to 1, to notify
// that this function is called at least once. It will be used
// in tandem with parameter common. So, this "called" should be set
// to 0 before calling this function the first time
{
	int i, n;
	if (xOutput == NULL || nTemp == 0) return;
//	char method[3] = "";
//	for (i=0; i<2; i++) method[i]=label[i];
//	method[2] ='\0';
	char *method = (char*) malloc(sizeof(char)*3);
	*method = '\0';
	strncat(method, label, 2);
	for (n=stCrit; n<nCrit; n++) {
	// Dec 2016:
		if (*(critVal+n) > 0 && *(critVal+n) <= PCRITX) continue;
		// print k blanks if common = 1:
		if (common == 1 && (*called == 1)) PrtPair (xOutput, 0, "", k, 0);
		if (nPlan > 1) {	// has only plan I, or has both plans in input
			if (census > 0) fprintf (xOutput, "%6d   ", census);
			else fprintf (xOutput, "%9s", " ");
		};
		for (i = nSkip1; i > 0; i--) fprintf (xOutput, "%c", ' ');
		fprintf (xOutput, "%s", popPair);
		for (i=0; i<nSkip2; i++) fprintf (xOutput, "%c", ' ');
		fprintf (xOutput, "%6.1f & %-6.1f", time1, time2);
		if (nTemp > 1) fprintf (xOutput, "%4s ", method);
		if (critOut > 1) fprintf (xOutput, "%7.4f ", *(critVal+n));
		fprintf (xOutput, "%9.1f  ", *(Hmean+n));
		fprintf (xOutput, "%8lu", *(nAlle+n));
		fprintf (xOutput, "%9.5f", *(fmean+n));
		fprintf (xOutput, "%10.5f", *(fprimeMean+n));
		if (*(Ne+n) < infinite)
			fprintf (xOutput, "%9.1f", *(Ne+n));
		else
			fprintf (xOutput, "%9s", "Infinite");
		if (param == 1) {
			if (*(loNe+n) < infinite)
				fprintf (xOutput, "%10.1f", *(loNe+n));
			else
				fprintf (xOutput, "%10s", "Infinite");
			if (*(hiNe+n) < infinite)
				fprintf (xOutput, "%9.1f", *(hiNe+n));
			else
				fprintf (xOutput, "%9s", "Infinite");
		};
		if (jack == 1) {
			if (*(jloNe+n) < infinite)
				fprintf (xOutput, "%10.1f", *(jloNe+n));
			else
				fprintf (xOutput, "%10s", "Infinite");
			if (*(jhiNe+n) < infinite)
				fprintf (xOutput, "%9.1f", *(jhiNe+n));
			else
				fprintf (xOutput, "%9s", "Infinite");
		};
		fprintf (xOutput, "\n");
		*called = 1;
	};
	free (method);
}


// --------------------------------------------------------------------------
// Modified Dec 2016:
// Recalculate stCrit, nCrit and critOut when the spec Pcrit is present
// in the PCrit array. Then adjust the "for" loop.
void PrtTempxFile (FILE *xOutput, int popRead, int popStart, char newSet,
					char lastSet, float *critVal, int nCrit, int g1, int g2,
					int nGeneration, float timeline[], long *nIndAlle,
					float *Hkmean, float *Hcmean, float *Hsmean,
					float *fkmean, float *fcmean, float *fsmean,
					float *fkprimeMean, float *fcprimeMean, float *fsprimeAll,
					float *NeTempk, float *NeTempc, float *NeTemps,
					float *loNek, float *loNec, float *loNes,
					float *hiNek, float *hiNec, float *hiNes,
					float *jloNek, float *jloNec, float *jloNes,
					float *jhiNek, float *jhiNec, float *jhiNes,
					char param, char jack,
					float infinite, char tempk, char tempc, char temps,
					char *popID[], int topCrit, int nPlan, int census,
// add May 2013 to be able to print input file name and loci as columns,
// and not to print header when necessary.
// count represents input file numbering, common = 1 if this xOutput is for
// all input files, with only one header to be printed.
// In this, count is the number of input so far.
// and smPair is the number of sample pairs so far.
					char *inpName, int nloci, int count, int smPair, char common)
// print to extra ouput (shorter form) file.
{
	int nCI;
	int stCrit, critOut;
	int dashes = 90;
	int i, k, m, n1, n2, lenInp;
	int nc = 10;	// number of characters to print popID, and set size
					// of popPair to be 2*(nc+2) to accomodate 2 pops,
					// together with string "<->".
	char nTemp;
	char newInp, called;
	char *popPair, *method;
	nTemp = tempk + tempc + temps;
	if (nTemp == 0) return;
	if (xOutput == NULL) return;
//	char popPair[24] = "\0";
	popPair = (char*) malloc(sizeof(char)*24);
	*popPair = '\0';

//	char method[45] = "\0";
	method = (char*) malloc(sizeof(char)*45);
	*method = '\0';
// Uncomment for testing only
//topCrit = 2;
//tempk = 0;
//tempc = 0;
//temps = 0;
//param = 0;
//jack = 0;
	nTemp = tempk + tempc + temps;
	if (nTemp == 0) return;
	if (xOutput == NULL) return;
	if (nPlan > 1) dashes +=9;	// has only plan I, or has both plans in input
	stCrit = 0;

// The next 3 lines are replaced by CritEndLine, as done in PrtHetxFile
//	if (topCrit > 0 && topCrit < nCrit) nCrit = topCrit;
//	else if (topCrit == 0) stCrit = nCrit - 1;
//	critOut = nCrit - stCrit;

// Jan 2017: Use function CritEndLine to set index limits
	nCrit = CritEndLine (critVal, nCrit, topCrit, &stCrit, &critOut, 0);

	if (critOut > 1) dashes += 8;	// more chars for printing crit. values

	nCI = param + jack;
	if (nTemp == 1) dashes -= 5;
	dashes += nCI*19;
	if (common == 1) dashes += 27;	// add 27 chars for input name and #loci
	lenInp = 19;	// #chars for input name and order, when common = 1
	newInp = (popRead == popStart-1+nGeneration)? 1: 0;
	if ((newInp == 1 && g1+g2 == 1) && (common == 0 || smPair == 1)) {
	// print header
		if (nTemp > 1) {	// have more than one temporal method
			if (nPlan == 1)		// plan II
				fprintf (xOutput,
					"\nOutput for %d TEMPORAL METHODS, Plan II:\t", nTemp);
			else if (nPlan == 2)	// plan I
				fprintf (xOutput,
						"\nOutput for %d TEMPORAL METHODS, Plan I:\t", nTemp);
			else	// then nPlan = 3, have both plans in input
				fprintf (xOutput, "\nOutput for %d TEMPORAL METHODS:\t", nTemp);
		} else {
			if (nPlan == 1)
				fprintf (xOutput, "\nOutput for TEMPORAL METHOD, Plan II");
			else if (nPlan == 2)	// plan I
				fprintf (xOutput, "\nOutput for TEMPORAL METHOD, Plan I");
			else // may have more than one plan
				fprintf (xOutput, "\nOutput for TEMPORAL METHOD");
		}
		if (tempk == 1) {
			strcat(method, "Pk (Pollak)");
			if (nTemp > 1) strcat(method, ", ");
		};
		if (tempc == 1) {
			strcat(method, "NT (Nei/Tajima)");
			if (nTemp > 2 || (nTemp == 2 && tempk == 0)) strcat(method, ", ");
		};
		if (temps == 1) {
			strcat(method, "JR (Jorde/Ryman)");
		};
		if (nTemp > 1) fprintf (xOutput, "%s\n", method);
		else fprintf (xOutput, "%s\n", strchr(method, ' '));
	// else: print method starting from the first blank
		if (nPlan > 2) fprintf (xOutput,
			"Census size N is listed for plan I. If absent, it is Plan II.\n");
		fprintf (xOutput, "\n");
		if (critOut > 1) {
			fprintf (xOutput, "Lowest allele frequencies used: ");
			for (i=0; i<nCrit; i++) {
// Jan 2017: skip when critVal[i] is the spec Pcrit
				if (critVal[i] > 0 && critVal[i] <= PCRITX) continue;
				fprintf (xOutput, "%10.4f", critVal[i]);
			}
			fprintf (xOutput, "\n\n");
			fprintf (xOutput, "For each sample pair of one population, ");
			if (nTemp > 1) fprintf (xOutput, "and for each method, ");
			fprintf (xOutput, "outputs are\n");
			fprintf (xOutput, "in consecutive lines corresponding "
						"to frequencies in the order above.\n");
			fprintf (xOutput, "Consecutive pairs are separated by a blank line.\n");
			fprintf (xOutput, "Consecutive populations are separated by 2 blank lines.\n");
		} else {
			fprintf (xOutput, "Lowest allele frequency used: %8.4f\n", critVal[stCrit]);
			fprintf (xOutput, "Each sample pair of one population outputs one line");
			if (nTemp > 1) fprintf (xOutput, " for each method");
			fprintf (xOutput, ".\n");
			fprintf (xOutput, "Consecutive populations are separated by one blank line.\n");
		};
		if (common != 0) fprintf (xOutput,
				"Input Names are shown up to %d righmost characters.\n",
			lenInp - 2);
		fprintf (xOutput,
			"Up to 10 righmost characters are printed per sample name in Sample Pair.\n");

		PrtLines (xOutput, dashes, '-');
		// if output is common for all input files, have columns for
		// input file name and #loci        123456789012345678901234567
		if (common != 0) fprintf (xOutput, "Input File Number   #Loci  ");

		if (nPlan > 1) fprintf (xOutput, " Census  ");
		// Distinguish when this output is for one or more than one temporal
		// methods. If more than one, need to have a column specifying methods.
		if (critOut > 1) {
			if (nTemp > 1) {
				fprintf (xOutput, "%4sSample Pair IDs%6sGenerations  Method Crit.  Harmonic   #Indep."
								"%5sF%8sF'%8sNe", " ", " ", " ", " ", " ");
				if (nCI > 0) {
					if (nCI == 1) fprintf (xOutput, "%9sCIs for Ne", " ");
					else fprintf (xOutput, "%19sCIs for Ne", " ");
				};
				fprintf (xOutput, "\n");
				// if output is common for all input files, have columns for
				// input file name and #loci        123456789012345678901234567
				if (common != 0) fprintf (xOutput, "followed by :Name          ");
				if (nPlan > 1) fprintf (xOutput, " Size N  ");
				fprintf (xOutput, "(last 10 chars each ID)%22sValue  Mean Size%2sAlleles%26s",
									" ", " ", " ");
			} else {
				fprintf (xOutput, "%4sSample Pair IDs%6sGenerations   Crit.  Harmonic    #Indep."
								"%5sF%8sF'%8sNe", " ", " ", " ", " ", " ");
				if (nCI > 0) {
					if (nCI == 1) fprintf (xOutput, "%9sCIs for Ne", " ");
					else fprintf (xOutput, "%18sCIs for Ne", " ");
				};
				fprintf (xOutput, "\n");
				// if output is common for all input files, have columns for
				// input file name and #loci        123456789012345678901234567
				if (common != 0) fprintf (xOutput, "followed by :Name          ");
				if (nPlan > 1) fprintf (xOutput, " Size N  ");
				fprintf (xOutput, "(last 10 chars each ID)%16sValue  Mean Size%3sAlleles%26s",
									" ", " ", " ");
			};
		} else {
			if (nTemp > 1) {
				fprintf (xOutput, "%4sSample Pair IDs%6sGenerations  Method Harmonic   #Indep."
									"%5sF%8sF'%8sNe", " ", " ", " ", " ", " ");
				if (nCI > 0) {
					if (nCI == 1) fprintf (xOutput, "%9sCIs for Ne", " ");
					else fprintf (xOutput, "%19sCIs for Ne", " ");
				};
				fprintf (xOutput, "\n");
				// if output is common for all input files, have columns for
				// input file name and #loci        123456789012345678901234567
				if (common != 0) fprintf (xOutput, "followed by :Name          ");
				if (nPlan > 1) fprintf (xOutput, " Size N  ");
				fprintf (xOutput, "(last 10 chars each ID)%22sMean Size%2sAlleles%26s",
									" ", " ", " ");
			} else {
				fprintf (xOutput, "%4sSample Pair IDs%6sGenerations   Harmonic    #Indep."
								"%6sF%8sF'%8sNe", " ", " ", " ", " ", " ");
				if (nCI > 0) {
					if (nCI == 1) fprintf (xOutput, "%8sCIs for Ne", " ");
					else fprintf (xOutput, "%18sCIs for Ne", " ");
				};
				fprintf (xOutput, "\n");
				// if output is common for all input files, have columns for
				// input file name and #loci        123456789012345678901234567
				if (common != 0) fprintf (xOutput, "followed by :Name          ");
				if (nPlan > 1) fprintf (xOutput, " Size N  ");
				fprintf (xOutput, "(last 10 chars each ID)%16sMean Size%3sAlleles%26s",
									" ", " ", " ");
			};
		};
		if (param == 1) fprintf (xOutput, "%9sParametric", " ");
		if (jack == 1) fprintf (xOutput, "%10sJackKnife", " ");
		fprintf (xOutput, "\n");
		for (i=0; i<dashes; i++) fprintf (xOutput, "-");
		fprintf (xOutput, "\n");
	};	// end printing header
	free (method);
	fflush (xOutput);
// now, print the values:

	// This part is similar to the one in PrtTemporal.
	// 10 characters will be used for writing sample ID
	// We want to have string "<->" between two IDs. Since we
	// have to deal with both together without having any blanks between
	// them, we should find the correct number of blanks to fill.
	n1 = nc - strlen (popID[g1]);		// this is the number of blanks to be
			// added or the number of chars to be taken out from popID[g1]
	m = (n1 > 0)? 0: -n1;
	n2 = nc - strlen (popID[g2]);
	k = (n2 > 0)? 0: -n2;
//	for (i = n1; i > 0; i--) strcat(popPair, " ");
	strcat(popPair, popID[g1]+m);
	strcat(popPair, "<->");
	strcat(popPair, popID[g2]+k);
// for test only:
//	fprintf (xOutput, "[%s], ", popID[g1]+m);
//	fprintf (xOutput, "[%s], ", popID[g2]+k);
//	fprintf (xOutput, "%s\n", popPair);

	called = 0;
	if (common != 0) {	// print common values for the whole input file,
		// the order and name of file, and number of loci
		// insert a blank line before printing output for this input file
//		if (newInp == 1 && count > 1) fprintf (xOutput, "\n");
		PrtPair (xOutput, count, inpName, lenInp, newInp);
//		if (newInp == 1) fprintf (xOutput, "%6d  ", nloci);
//		else  fprintf (xOutput, "%8c", ' ');
		// comment out the previous 2 lines, which are for printing input
		// name and loci on the same first line of data.
		// The following are for putting name, loci# in a separate line
		if (newInp == 1) {
			fprintf (xOutput, "%6d\n", nloci);
			PrtLines (xOutput, lenInp+6, '-');
			PrtPair (xOutput, count, inpName, lenInp, 0);
		};
		fprintf (xOutput, "%8c", ' ');
	};
	lenInp += 8;	// add 8 to include chars reserved for printing #loci
	if (tempk == 1) PrtTempValx (xOutput, nCrit, critVal, nIndAlle,
				Hkmean, fkmean, fkprimeMean, NeTempk, loNek, hiNek,
				jloNek, jhiNek, param, jack, "Pk (Pollak)", popPair,
				timeline[g1], timeline[g2], n1, n2, nTemp, infinite,
				stCrit, critOut, nPlan, census, common, lenInp, &called);
	if (tempc == 1) PrtTempValx (xOutput, nCrit, critVal, nIndAlle,
				Hcmean, fcmean, fcprimeMean, NeTempc, loNec, hiNec,
				jloNec, jhiNec, param, jack, "NT (Nei/Tajima)", popPair,
				timeline[g1], timeline[g2], n1, n2, nTemp, infinite,
				stCrit, critOut, nPlan, census, common, lenInp, &called);
	if (temps == 1) PrtTempValx (xOutput, nCrit, critVal, nIndAlle,
				Hsmean, fsmean, fsprimeAll, NeTemps, loNes, hiNes,
				jloNes, jhiNes, param, jack, "JR (Jorde/Ryman)", popPair,
				timeline[g1], timeline[g2], n1, n2, nTemp, infinite,
				stCrit, critOut, nPlan, census, common, lenInp, &called);
	free (popPair);
// separate by a blank line if more than one critical value:
//	if (nCrit > 1) fprintf (xOutput, "\n");
// Jan 2017: change nCrit to critOut:
	if (critOut > 1) fprintf (xOutput, "\n");
// separate by one more blank line at the end of each population:
	if (g1+g2 == 2*nGeneration-3) fprintf (xOutput, "\n");
	fflush (xOutput);
}
// --------------------------------------------------------------------------
void PrtTempPop (FILE *output, int generation, int nGeneration,
				int nPoptemp, int popRun, int nPairTmp, float timeline[])
{
	int n, g1;
	if (output == NULL) return;
	fprintf(output, "\nNumber of populations = %d, with %d samples\n",
					nPoptemp, popRun);
	fprintf(output, "Number of sample pairs analyzed = %d\n", nPairTmp);
	n = 42;	// for number of dashes to draw for the last line
	if (generation < nGeneration-1) {
		n += 12;
		fprintf(output, "Last population has sample(s) taken only at generation");
		if (generation > 0) fprintf(output, "s");
		for (g1=0; g1<= generation; g1++) {
			fprintf(output, " %5.1f", timeline[g1]);
			n += 6;
		};
		fprintf(output, "\n");
	};
	PrtLines (output, n, '-');
	fflush (output);
}

// --------------------------------------------------------------------------
// Add PrtLeading in July 2013 for printing to extra file with Tab delimiter,
// which is called by PrtLDHeader, PrtCnHeader, PrtHtHeader, PrtTpHeader.
// These 4 functions are called by Prt**TabFile's.
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
// Jan 2017:
// Have *stCrit as a parameter (not local), and return value (not "void")
// for the reassigned value of nCrit, since those will be used in
// Prt**TabFile's, after the call Prt**Header that calls PrtLeading.
// The added return values for PrtLeading and stCrit will be registered as
// parameters to be added in Prt**Header's, so that we don't have to redo
// the code of function PrtLeading

int PrtLeading (FILE *xOutput, float *critVal, int nCrit, int topCrit,
				int *stCrit, int *critOut, char method[], char keyName[],
				int lenInp, int lenName, char common, char skip, char forLD)
{
	int i;
	fprintf (xOutput, "\n%s\n\n", method);
	char specP = 0;

	if (nCrit > 0) {
		// the next 4 lines are replaced by CritEndLine (stCrit was local)
//		stCrit = 0;
//		if (topCrit > 0 && topCrit < nCrit) nCrit = topCrit;
//		else if (topCrit == 0) stCrit = nCrit - 1;
//		*critOut = nCrit - stCrit;
// Jan 2017: Use function CritEndLine to set index limits
		nCrit = CritEndLine (critVal, nCrit, topCrit, stCrit, critOut, forLD);
		if (*critOut > 1) {	// then *stCrit must be 0
			fprintf (xOutput, "Lowest allele frequencies used: ");
			for (i=*stCrit; i<nCrit; i++) {
// Jan 2017:
				if (critVal[i] > 0 && critVal[i] <= PCRITX) {
					if (forLD != 0) {
						fprintf (xOutput, "   \"%s\"", NOSNGL);
						specP = 1;
					}
					continue;
				}
				fprintf (xOutput, "%10.4f", critVal[i]);
			}
			fprintf (xOutput, "\n");
		} else {
// Jan 2017: cover the printing below to handle spec Pcrit
//			fprintf (xOutput, "Lowest allele frequency used: %8.4f\n",
//							critVal[*stCrit]);
			fprintf (xOutput, "Lowest allele frequency used: ");
			if (forLD != 0 && critVal[i] > 0 && critVal[i] <= PCRITX) {
				fprintf (xOutput, "   \"%s\"\n", NOSNGL);
				specP = 1;
			} else fprintf (xOutput, "%8.4f\n", critVal[*stCrit]);
		}
		if (specP != 0) fprintf (xOutput,
			"(\"%s\": No Singleton Allele is accepted)\n", NOSNGL);
	}	// nCrit = 0 when this is called by Coan method

	if (common != 0) fprintf (xOutput, "A maximum of %d"
			" rightmost characters can be shown for Input name.\n", lenInp);
	fprintf (xOutput, "Up to %d rightmost characters are shown for %s.\n",
//	                   1234567890123456789012345678901234567890123456
						lenName, keyName);
	i = 46 + strlen(keyName);
	PrtLines (xOutput, i, '-');
	if (skip != 0) fprintf (xOutput, "\n");
	return nCrit;
}

// --------------------------------------------------------------------------
// Jan 2017: add param stCrit, which is a param in PrtLeading, and change
// var critOut to be parameter. Function will return nCrit (instead of void),
// which was reassigned by the return value of PrtLeading.
// These will be used later in the calling function PrtLDTabFile.

int PrtLDHeader (FILE *xOutput, float *critVal, int nCrit, int topCrit,
				int *stCrit, int *critOut, char param, char jack, int lenInp,
				int lenName, char mating, char common, long *Jdegree)
{
//	char nTemp, nCI;
//	int i, k, m;
//	int critOut = 1;
	*critOut = 1;
	char *method = (char*) malloc(sizeof(char)*200);
	*method = '\0';
	sprintf (method, "%s", "LINKAGE DISEQUILIBRIUM METHOD, Mating Model: ");
	if (mating == 0) strcat (method, "Random");
	else strcat (method, "Monogamy");
// Jan 2017: add param stCrit, and the last parameter = 1 for LD:
	nCrit = PrtLeading (xOutput, critVal, nCrit, topCrit, stCrit, critOut,
				method, "Population name", lenInp-2, lenName-2, common, 0, 1);
	// if output is common for all input files, have columns for
	// input file name and #loci        12345678901234567  12345
	int i;
	char c = ' ';
	if (common != 0) {	// 1234567890123456789
		fprintf (xOutput,
			"%19c\t\%6c\t%19c\t%9c\t\%6c\t%13c\t%12c\t%9c\t%10c\t%9c",
			c,c,c,c,c,c,c,c,c,c);
// first 19 = length of "Input File [#:Name]", then 6 = length of " #Loci",
// next 19 = length of "Population [#:Name]", 9 = length of "Samp Size",
// 6 = length of "PCrit.", 13 = "Weighted Mean", 12 = "Ind. Alleles",
// 9 = "    r^2  ", 10 = "  Exp(r^2)", 9 = "     Ne  ".
// The spaces are reserved for those headlines before headlines for CI.
// When the output files are for only one input file (common = 0),
// the first two reserved spaces (19 blanks, \t, 6 blanks, \t) are gone
// The length of all these headlines should match the lengths
// for outputs of these values to be put in PrtLDTabFile
		if (param == 1) fprintf (xOutput, "\t    Parametric CI\t");
		if (jack == 1) fprintf (xOutput, "\t    Jackknife CI");
		fprintf (xOutput, "\n");
		fprintf (xOutput, "Input File [#:Name]\t");	// this has length 19
		for (i=0; i<lenInp-19; i++) fprintf (xOutput, " ");
		fprintf (xOutput, " #Loci\t");
	} else {
		fprintf (xOutput, "%19c\t%9c\t\%6c\t%13c\t%12c\t%9c\t%10c\t%9c",
			c,c,c,c,c,c,c,c);
		if (param == 1) fprintf (xOutput, "\t    Parametric CI\t");
		if (jack == 1) fprintf (xOutput, "\t    Jackknife CI");
		fprintf (xOutput, "\n");
	}
	fprintf (xOutput, "Population [#:Name]\tSamp Size");
//                     1234567890123456789  123456789
	if (*critOut > 1) fprintf (xOutput, "\tPCrit.");
//                                        123456
	fprintf (xOutput,
		"\tWeighted Mean\tInd. Alleles\t    r^2  \t  Exp(r^2)\t     Ne  ");
//	       1234567890123  123456789012  123456789  1234567890  123456789
// Headlines "Low, High" occupy 10 spaces to match with 10 spaces for values
// to be printed in PrtLDTabFile
	if (param == 1) fprintf (xOutput, "\t       Low\t      High");
//                                       1234567890  1234567890
	if (jack == 1) fprintf (xOutput, "\t       Low\t      High\t  (Eff.df)");
//                                      1234567890  1234567890  1234567890
	fprintf (xOutput, "\n");
	return nCrit;
}

// --------------------------------------------------------------------------
// Jan 2017:
// Add stCrit, critOut into the call PrtLDHeadeer, and remove the code for
// calculating these before the "for" loop (n runs from stCrit to nCrit).

void PrtLDTabFile (char *inpName, FILE *xOutput, int samp, float *wHarmonic,
				  int popRead, int popStart, char *popID, float *critVal,
				  int nCrit, double *nIndSum, float *rB2WAve, float *wExpR2,
				  float *estNe, char param, char jack, float infinite,
				  float *confParalow, float *confParahi, float *confJacklow,
				  float *confJackhi, char *jackOK, char mating, int topCrit,
// add May-July 2013 to be able to print input file name and loci as columns,
// and not to print header when necessary.
// count represents input file numbering, common = 1 if this xOutput is for
// all input files, with only one header to be printed.
// In this, count is the number of input file so far.
				  int nloci, int count, char common, long *Jdegree)
// print to extra ouput (shorter form) file.
{
	int i, n, lenInp, lenName, lenSkip;
	int stCrit, critOut;
	char pair;
	float indMax = (float) MAXLONG;
	unsigned long indPrt;
	char *skipStr;
	if (xOutput == NULL) return;
// don't print jackknife if not OK
// If this block is in effect, no header for jackknife CI if the first
// population doesn't have jackknife calculated (see PrtLDHeader). As a
// result, for later population, jackknife may be calculated and printed here,
// but without header saying what it is. So, we comment out
//	for (n=0, i=0; n<nCrit; n++) {
//		if (*(jackOK+n) == 1) i++;
//	};
//	if (i == 0) jack = 0;
//	stCrit = 0;

	skipStr = (char*) malloc(sizeof(char)*100);
	*skipStr = '\0';
	lenInp = 19;	// max width for printing Input name and number.
	lenName = 19;	// 19 is the length of "Population [#:Name]".
					// lenName= max width for printing Population + number.
	if ((popRead == popStart) && (common == 0 || count == 1))
		nCrit = PrtLDHeader (xOutput, critVal, nCrit, topCrit, &stCrit,
							&critOut, param, jack, lenInp, lenName,
							mating, common, Jdegree);
// Jan 2017:
// stCrit, nCrit, critOut need to be modified when spec Pcrit is present
	else
		nCrit = CritEndLine (critVal, nCrit, topCrit, &stCrit, &critOut, 1);
	lenSkip = 0;
	if (common != 0) {	// print order and name of input file,
						// and number of loci
		pair = (popRead == popStart)?1:0;
		PrtPair (xOutput, count, inpName, lenInp, pair);
		if (pair == 1) fprintf (xOutput, "\t%6d\t", nloci);
		else fprintf (xOutput, "\t%6c\t", ' ');
		for (i=0; i<lenInp; i++) *(skipStr+i) = ' ';
		*(skipStr+lenInp) = '\t';
		lenSkip = lenInp + 1;
		for (i=0; i<6; i++) *(skipStr+lenSkip+i) = ' ';
		*(skipStr+lenSkip+6) = '\t';
		lenSkip += 7;
		*(skipStr+lenSkip) = '\0';
	};
// now, print the values:

// Jan 2017:
// Skip next 4 lines: stCrit, nCrit, critOut were done in calling PrtLDHeader
//	stCrit = 0;
//	if (topCrit > 0 && topCrit < nCrit) nCrit = topCrit;
//	else if (topCrit == 0) stCrit = nCrit - 1;
//	critOut = nCrit - stCrit;

	for (n=stCrit; n<nCrit; n++) {
		pair = (n==stCrit)?1:0;	// (if pair ever = 0, then critOut > 1)
		if (common !=0 && pair == 0) fprintf (xOutput, "%s", skipStr);
		PrtPair (xOutput, popRead, popID, lenName, pair);
		// print Ind. Alleles, r^2, exp(r^2), Ne:
		if (indMax <= *(nIndSum+n)) indPrt = MAXLONG;
		else indPrt = (unsigned long) *(nIndSum+n);
		if (critOut > 1) {
		// 9 spaces for Samp. size (printed for the first PCrit: 7+2 blanks),
		// 6 spaces for PCrit, 13 (10+3 blanks)for Harmonic mean,
		// 12 (10 + 2) for Ind. Alleles, 9 for r^2, 10 for Exp(r^2).
			if (n==stCrit) fprintf (xOutput, "\t%7d  ", samp);
			else fprintf (xOutput, "\t%9s", " ");
// Jan 2017, adjust printing when spec Pcrit is present:
			if (critVal[n] > 0 && critVal[n] <= PCRITX)
				fprintf (xOutput, "\t%6s", NOSNGL);
			else fprintf (xOutput, "\t%6.4f", critVal[n]);
			fprintf (xOutput, "\t%10.1f   \t%10lu  \t%9.6f\t%10.6f\t",
						*(wHarmonic+n), indPrt, *(rB2WAve+n), *(wExpR2+n));
		} else fprintf (xOutput,
					"\t%7d  \t%10.1f   \t%10lu  \t%9.6f\t%10.6f\t",
						samp, *(wHarmonic+n), indPrt,
						*(rB2WAve+n), *(wExpR2+n));
		// 9 spaces for Ne, then 10 for Low, High in CI
		if (*(estNe+n) < infinite) fprintf (xOutput, "%9.1f", *(estNe+n));
		else fprintf (xOutput, "%9s", "Infinite");
		if (param == 1) {
			if (confParalow[n] < infinite && confParalow[n] >= 0)
				fprintf (xOutput, "\t%10.1f", confParalow[n]);
			else fprintf (xOutput, "\t%10s", "Infinite");
			if (confParahi[n] < infinite && confParahi[n] >= 0)
				fprintf (xOutput, "\t%10.1f", confParahi[n]);
			else fprintf (xOutput, "\t%10s", "Infinite");
		};
		if (jack == 1) {
			if (*(jackOK+n) == 1) {
				if (confJacklow[n] < infinite && confJacklow[n] >= 0)
					fprintf (xOutput, "\t%10.1f", confJacklow[n]);
				else fprintf (xOutput, "\t%10s", "Infinite");
				if (confJackhi[n] < infinite && confJackhi[n] >= 0)
					fprintf (xOutput, "\t%10.1f", confJackhi[n]);
				else fprintf (xOutput, "\t%10s", "Infinite");
// add in Aug 2016:
				fprintf (xOutput, "\t%10lu", Jdegree[n]);
			} else {
				fprintf (xOutput, "\t%10s", "     *  ");
				fprintf (xOutput, "\t%10s", "   *    ");
			};
		};
		fprintf (xOutput, "\n");
	};
	free (skipStr);
	fflush (xOutput);
}

// --------------------------------------------------------------------------
// Jan 2017: add param stCrit, which is a param in PrtLeading, and change
// var critOut to be parameter. Function will return nCrit (instead of void),
// which was reassigned by the return value of PrtLeading.
// These will be used later in the calling function PrtHtTabFile.

int PrtHtHeader (FILE *xOutput, float *critVal, int nCrit, int topCrit,
					int *stCrit, int *critOut, char param, int lenInp,
					int lenName, char common)
{
//	int critOut;
// Jan 2017: add param stCrit, and the last parameter = 0 for non-LD:
	nCrit = PrtLeading (xOutput, critVal, nCrit, topCrit, stCrit, critOut,
				"HETEROZYGOTE EXCESS METHOD", "Population name",
				lenInp-2, lenName-2, common, 1, 0);
	// if output is common for all input files, have columns for
	// input file name and #loci        12345678901234567  12345
	if (common != 0) {	// 1234567890123456789
		int i;
		fprintf (xOutput, "Input File [#:Name]\t");	// this has length 19
		for (i=0; i<lenInp-19; i++) fprintf (xOutput, " ");
		fprintf (xOutput, " #Loci\t");
	}
	fprintf (xOutput, "Population [#:Name]\tSamp Size");
//                     1234567890123456789  123456789
	if (*critOut > 1) fprintf (xOutput, "\tPCrit.");
//                                         123456
	fprintf (xOutput, "\tHarmonic Mean\tInd. Alleles\t    D   \t    Ne  ");
//	                     1234567890123  123456789012  12345678  12345678
	if (param == 1) fprintf (xOutput, "\t   Parametric CI for Ne");
//                                       12345678901234567890123
	fprintf (xOutput, "\n");
	return nCrit;

}

//------------------------------------------------------------------

// Modified Jan 2017:
// Recalculate stCrit, nCrit and critOut when the spec Pcrit is present
// in the PCrit array, from PrtHtHeader. Then adjust the "for" loop.
void PrtHetTabFile (char *inpName, FILE *xOutput, int popRead, int popStart,
				  char *popID, float *critVal, int nCrit, long *indAlleH,
				  float *hetD, float *estHetN, char param, float infinite,
				  float *loHetNe, float *hiHetNe, int samp, float *hSamp,
				  int topCrit,
// add May 2013 to be able to print input file name and loci as columns,
// and not to print header when necessary.
// count represents input file numbering, common = 1 if this xOutput is for
// all input files, with only one header to be printed.
// In this, count is the number of input so far.
				  int nloci, int count, char common)
// print to tabular format ouput file.
{
	int stCrit, critOut;
	int i, n, lenSkip, lenInp, lenName;
	char pair;
	char *skipStr;
	if (xOutput == NULL) return;
	skipStr = (char*) malloc(sizeof(char)*100);
	*skipStr = '\0';
	lenInp = 19;	// max width for printing Input name and number.
	lenName = 19;	// 19 is the length of "Population [#:Name]".
					// lenName= max width for printing Population + number.
//	stCrit = 0;
	if ((popRead == popStart) && (common == 0 || count == 1))
		nCrit = PrtHtHeader (xOutput, critVal, nCrit, topCrit, &stCrit,
							&critOut, param, lenInp, lenName, common);
// Jan 2017:
// stCrit, nCrit, critOut need to be modified when spec Pcrit is present
	else
		nCrit = CritEndLine (critVal, nCrit, topCrit, &stCrit, &critOut, 0);

	lenSkip = 0;
	if (common != 0) {	// print order and name of input file,
						// and number of loci
		pair = (popRead == popStart)?1:0;
		PrtPair (xOutput, count, inpName, lenInp, pair);
		if (pair == 1) fprintf (xOutput, "\t%6d\t", nloci);
		else fprintf (xOutput, "\t%6c\t", ' ');
		for (i=0; i<lenInp; i++) *(skipStr+i) = ' ';
		*(skipStr+lenInp) = '\t';
		lenSkip = lenInp + 1;
		for (i=0; i<6; i++) *(skipStr+lenSkip+i) = ' ';
		*(skipStr+lenSkip+6) = '\t';
		lenSkip += 7;
		*(skipStr+lenSkip) = '\0';
	};
// now, print the values:
// Jan 2017: the caluclation of stCrit, critOut, and reassignment of nCrit
// were done in PrtLeading, that was called by PrtHtHeader, so
// the code for those was removed.

// Jan 2017: add stCrit0 as the truly starting index since Pcrit at index
// stCrit may be the spec Pcrit, which is skipped.
// Then stCrit inside the loop is replaced by stCrit0.
	int stCrit0 = stCrit;
	for (n=stCrit; n<nCrit; n++) {
// Modify Jan 2017: skip when critVal[i] is the spec Pcrit
		if (critVal[n] > 0 && critVal[n] <= PCRITX) {
			if (n == stCrit) stCrit0++;
			continue;
		}
		pair = (n==stCrit0)?1:0;	// (if pair ever = 0, then critOut > 1)
		if (common !=0 && pair == 0) fprintf (xOutput, "%s", skipStr);
		PrtPair (xOutput, popRead, popID, lenName, pair);
		if (critOut > 1) {
			if (n==stCrit0)
				fprintf (xOutput, "\t%7d  \t%6.4f\t%10.1f   \t%10lu  \t%8.5f",
						samp, critVal[n], *(hSamp+n), *(indAlleH+n), *(hetD+n));
			else fprintf (xOutput, "\t%9s\t%6.4f\t%10.1f   \t%10lu  \t%8.5f",
						" ", critVal[n], *(hSamp+n), *(indAlleH+n), *(hetD+n));
		} else fprintf (xOutput, "\t%7d  \t%10.1f   \t%10lu  \t%8.5f",
						samp, *(hSamp+n), *(indAlleH+n), *(hetD+n));
		if (*(estHetN+n) >= infinite || *(estHetN+n) < 0)
			fprintf (xOutput, "\t%8s", "Infinite");
		else fprintf (xOutput, "\t%8.1f", *(estHetN+n));
		if (param == 1) {
			if (loHetNe[n] < infinite)
				fprintf (xOutput, "\t%9.1f", loHetNe[n]);
			else fprintf (xOutput, "\t%9s", "Infinite");
			if (hiHetNe[n] < infinite)
				fprintf (xOutput, "\t%9.1f", hiHetNe[n]);
			else fprintf (xOutput, "\t%9s", "Infinite");
		};
		fprintf (xOutput, "\n");
	};
	free (skipStr);
	fflush (xOutput);
}

//---------------------------------------------------------------------------
void PrtCnHeader (FILE *xOutput, char jack, int lenInp, int lenName, char common)
{
	PrtLeading (xOutput, NULL, 0, 0, NULL, NULL,
				"MOLECULAR COANCESTRY METHOD", "Population name",
				lenInp-2, lenName-2, common, 1, 0);
	// if output is common for all input files, have columns for
	// input file name and #loci        12345678901234567  12345
	if (common != 0) {	// 1234567890123456789
		int i;
		fprintf (xOutput, "Input File [#:Name]\t");	// this has length 19
		for (i=0; i<lenInp-19; i++) fprintf (xOutput, " ");
		fprintf (xOutput, " #Loci\t");
	}
	fprintf (xOutput,
		"Population [#:Name]\tSamp Size\tHarmonic Mean\t   f^1   \t     Neb^ ");
//	     1234567890123456789  123456789  1234567890123  123456789  1234567890
	if (jack == 1) fprintf (xOutput, "\t  Jackknife CI for Neb^");
	fprintf (xOutput, "\n");
}


//---------------------------------------------------------------------------

void PrtCoanTabFile (char *inpName, FILE *xOutput, int popRead, int popStart,
				  char *popID, float f1, float coanNeb, char jack,
				  float infinite, float loNbCoan, float hiNbCoan, int samp,
				  float hSamCoan,
// add May 2013 to be able to print input file name and loci as columns,
// and not to print header when necessary.
// count represents input file numbering, common = 1 if this xOutput is for
// all input files, with only one header to be printed.
// In this, count is the number of input so far.
				  int nloci, int count, char common)
// print to extra ouput (shorter form) file.
{
	int lenInp, lenName;
	char newInp;
	if (xOutput == NULL) return;
	lenInp = 19;	// max width for printing Input name and number.
	lenName = 19;	// 19 is the length of "Population [#:Name]".
					// lenName= max width for printing Population + number.
	if ((popRead == popStart) && (common == 0 || count == 1))
		PrtCnHeader (xOutput, jack, lenInp, lenName, common);

	if (common != 0) {	// print common values for the whole input file,
		// the order and name of file, and number of loci
		newInp = (popRead == popStart)?1:0;
		PrtPair (xOutput, count, inpName, lenInp, newInp);
		if (newInp == 1) fprintf (xOutput, "\t%6d\t", nloci);
		else fprintf (xOutput, "\t%6c\t", ' ');
	};

	PrtPair (xOutput, popRead, popID, lenName, 1);
	fprintf (xOutput, "\t%7d  \t%10.1f   \t%9.5f", samp, hSamCoan, f1);
	if (coanNeb >= infinite || coanNeb < 0)
		fprintf (xOutput, "\t%10s", "Infinite");
	else fprintf (xOutput, "\t%10.1f", coanNeb);
	if (jack == 1) {
		if (loNbCoan < infinite)
			fprintf (xOutput, "\t%10.1f", loNbCoan);
		else fprintf (xOutput, "\t%10s", "Infinite");
		if (hiNbCoan < infinite)
			fprintf (xOutput, "\t%10.1f", hiNbCoan);
		else fprintf (xOutput, "\t%10s", "Infinite");
	};
	fprintf (xOutput, "\n");
	fflush (xOutput);
}

//---------------------------------------------------------------------------
// Revised in Oct 2016
// Jan 2017: add param stCrit, which is a param in PrtLeading, and change
// var critOut to be parameter. Function will return nCrit (instead of void),
// which was reassigned by the return value of PrtLeading.
// These will be used later in the calling function PrtTpTabFile.

int PrtTpHeader (FILE *xOutput, float *critVal, int nCrit, int topCrit,
					int *stCrit, int *critOut,
					int nPlan, char tempk, char tempc, char temps,
					char param, char jack, int lenInp, int lenPair,
					int lenPop, char common, int nGeneration)
// lenInp is the length to print Input Name, lenPop is for population
// lenInp should be at least 19 (currently is).
// lenPair = 10 when called
{
	char nTemp, nCI;
	int i, j, k, m, n;
//	int critOut = 1;
	*critOut = 1;
	char *method = (char*) malloc(sizeof(char)*200);
	*method = '\0';

	nTemp = tempk + tempc + temps;
	nCI = param+jack;
	if (nTemp == 0) return 0;
	if (xOutput == NULL) return 0;

	if (nTemp > 1) {
		if (nPlan == 1)		// plan II
			sprintf (method,
				"Output for %d TEMPORAL METHODS, Plan II:\t", nTemp);
		else if (nPlan == 2)	// plan I
			sprintf (method,
				"Output for %d TEMPORAL METHODS, Plan I:\t", nTemp);
		else	// then nPlan = 3, may have both plans
			sprintf (method, "Output for %d TEMPORAL METHODS:\t", nTemp);
	} else {
		if (nPlan == 1)
			sprintf (method, "Output for TEMPORAL METHOD, Plan II ");
		else if (nPlan == 2)	// plan I
			sprintf (method, "Output for TEMPORAL METHOD, Plan I ");
		else // may have more than one plan
			sprintf (method, "Output for TEMPORAL METHOD ");
	};
	if (nTemp == 1) strcat(method, "(");
	if (tempk == 1) {
		strcat(method, "Pollak");
		if (nTemp > 1) strcat(method, ", ");
	}
	if (tempc == 1) {
		strcat(method, "Nei/Tajima");
		if (nTemp > 2 || (nTemp == 2 && tempk == 0)) strcat(method, ", ");
	}
	if (temps == 1) {
		strcat(method, "Jorde/Ryman");
	}
	if (nTemp == 1) strcat(method, ")");
	if (nPlan > 2) strcat (method,
		"\nCensus size N is listed for plan I. If absent, it is Plan II.");
// Jan 2017:
// nCrit as return value for PrtLeading, to be the return value of this.
// add param stCrit, and the last parameter = 0 for non-LD:
	nCrit = PrtLeading (xOutput, critVal, nCrit, topCrit, stCrit, critOut,
					method, "each name in Sample Pair", lenInp-2, lenPair,
					common, 0, 0);

// 10 is the length of string "Population",
// 22 = length of "Pop. [up to %d samples]"
	m = (common == 0)? 10: 22;
	k = lenPop - m;
	n = 9;	// to be the length of 2 interested headers in each method before CI
//	if (nTemp > 1) {
	if (common !=0) {
		// this will occupy same spaces as headers for input file and #loci
		fprintf (xOutput, "%19s\t", " ");
		// lenInp assigned in PrtTpTabFile, currently = 19, same as header
		for (i=0; i<lenInp-19; i++) fprintf (xOutput, " ");
		fprintf (xOutput, "%6s\t", " ");
	}
		// Write empty spaces for "Population" or "Pop. [up to %d samples]"
	for (i=0; i<m; i++) fprintf (xOutput, " ");
	for (i=0; i<k; i++) fprintf (xOutput, " ");	// do nothing if lenPop<m
	fprintf (xOutput, "\t");
		 // Write spaces for "Census Size N" (length = 13)
	if (nPlan > 1) fprintf (xOutput, "%13s\t", " ");
		// spaces for sample pair name + string "<->":
	for (i=0; i<2*lenPair+3; i++) fprintf (xOutput, " ");
	fprintf (xOutput, "\t");
	// print spaces for "\tGenerations\tPCrit.\tInd. Alleles" or
//	                           12345678901  123456  123456789012
	// "\tGenerations\tInd. Alleles"
//	          12345678901  123456789012
	if (*critOut > 1)
		fprintf (xOutput, "%11s\t%6s\t%12s\t", " ", " ", " ");
	else fprintf (xOutput, "%11s\t%12s\t", " ", " ");


	if (tempk==1) {
		if (nTemp > 1) fprintf (xOutput, "  Pollak\t");	// the length = 8
//	                                      12345678
		else fprintf (xOutput, "%8s\t", " ");
		// there are 3 values for the method before CI:
		// first one occupies 8, the next 2 occupy n spaces:
		for (i = 0; i < 2; i++) {
			for (j = 0; j < n; j++) fprintf (xOutput, " ");
			fprintf(xOutput, "\t");
		}
		if (param == 1) {
			fprintf (xOutput, "Parametric CI\t");
//	                           1234567890123
			for (j = 0; j < n; j++) fprintf (xOutput, " ");
			fprintf(xOutput, "\t");
		}
		if (jack == 1) {
			fprintf (xOutput, "JackKnife CI\t");
//	                           123456789012
			for (j = 0; j < n; j++) fprintf (xOutput, " ");
			fprintf(xOutput, "\t");
		}
	}
	if (tempc==1) {
		if (nTemp > 1) fprintf (xOutput, "Nei/Tajima\t");	// length = 10
//	                                      1234567890
		else fprintf (xOutput, "%10s\t", " ");
		// there are 3 values for the method before CI,
		// first one occupies 10, the next 2 occupy n spaces:
		for (i = 0; i < 2; i++) {
			for (j = 0; j < n; j++) fprintf (xOutput, " ");
			fprintf(xOutput, "\t");
		}
		if (param == 1) {
			fprintf (xOutput, "Parametric CI\t");
			for (j = 0; j < n; j++) fprintf (xOutput, " ");
			fprintf(xOutput, "\t");
		}
		if (jack == 1) {
			fprintf (xOutput, "JackKnife CI\t");
			for (j = 0; j < n; j++) fprintf (xOutput, " ");
			fprintf(xOutput, "\t");
		}
	}
	if (temps==1) {
		if (nTemp > 1) fprintf (xOutput, "Jorde/Ryman\t");	// length = 11
//	                                      12345678901
		else fprintf (xOutput, "%11s\t", " ");
		// there are 3 values for the method before CI:
		// first one occupies 11, the next 2 occupy n spaces:
		for (i = 0; i < 2; i++) {
			for (j = 0; j < n; j++) fprintf (xOutput, " ");
			fprintf(xOutput, "\t");
		}
		if (param == 1) {
			fprintf (xOutput, "Parametric CI\t");
			for (j = 0; j < n; j++) fprintf (xOutput, " ");
			fprintf(xOutput, "\t");
		}
		if (jack == 1) {
			fprintf (xOutput, "JackKnife CI\t");
			for (j = 0; j < n; j++) fprintf (xOutput, " ");
			fprintf(xOutput, "\t");
		}
	}
	fprintf (xOutput, "\n");


	// if output is common for all input files, have columns for
	// input file name and #loci
	if (common != 0) {	// 1234567890123456789
		fprintf (xOutput, "Input File [#:Name]\t");	// this has length 19
		for (i=0; i<lenInp-19; i++) fprintf (xOutput, " ");
		fprintf (xOutput, " #Loci\t");
	}
// now print column headers
	if (common == 0) {	// there can be several Generation sets
		fprintf (xOutput, "Population");
//	                       1234567890
	} else {	// run multiple input with only one generation set
		fprintf (xOutput, "Pop. [up to %d samples]", nGeneration);
//	                       123456789012 3456789012
	};
//	k = lenPop - m;
	for (i=0; i<k; i++) fprintf (xOutput, " ");
	if (nPlan > 1) fprintf (xOutput, "\tCensus Size N");
	//                                  1234567890123
	fprintf (xOutput, "\t");
	n = lenPair - 4;	// 4 = (1/2) of [length("Sample Pair") - length("<->")]
	for (i=0; i<n; i++) fprintf (xOutput, " ");
	fprintf (xOutput, "Sample Pair");
	for (i=0; i<n; i++) fprintf (xOutput, " ");
	if (*critOut > 1) fprintf (xOutput, "\tGenerations\tPCrit.\tInd. Alleles");
	else fprintf (xOutput, "\tGenerations\tInd. Alleles");
	if (tempk==1) {
//		if (nTemp > 1) fprintf (xOutput, "\tPollak Fk\tPollak F'\tPollak Ne");
//	                                        123456789
//		else
		fprintf (xOutput, "\t    Fk  \t    F'   \t    Ne   ");	// n = 9
//	                         12345678  123456789  123456789
		fprintf (xOutput, "\t    Low      ");// 13 = length of "Parametric CI"
//	                         1234567890123
		fprintf (xOutput, "\t   High  ");	// n = 9
//	                         123456789
		fprintf (xOutput, "\t    Low     ");// 12 = length of "Jackknife CI"
//	                         123456789012
		fprintf (xOutput, "\t   High  ");	// n = 9
//	                         123456789
	}
	if (tempc==1) {
//		if (nTemp > 1) fprintf (xOutput, "\tNei/Tajima Fc\tNei/Tajima F'\tNei/Tajima Ne");
//		else
	// 10 is the length of "Nei/Tajima"
		fprintf (xOutput, "\t    Fc    \t    F'   \t    Ne   ");	// n = 9
//	                         1234567890  123456789  123456789
		fprintf (xOutput, "\t    Low      ");// 13 = length of "Parametric CI"
//	                         1234567890123
		fprintf (xOutput, "\t   High  ");	// n = 9
//	                         123456789
		fprintf (xOutput, "\t    Low     ");// 12 = length of "Jackknife CI"
//	                         123456789012
		fprintf (xOutput, "\t   High  ");	// n = 9
//	                         123456789
	}
	if (temps==1) {
//		if (nTemp > 1) fprintf (xOutput, "\tJorde/Ryman Fs\tJorde/Ryman F'\tJorde/Ryman Ne");
//		else
	// 11 is the length of "Jorde/Ryman"
		fprintf (xOutput, "\t    Fs     \t    F'   \t    Ne   ");
//	                         12345678901  123456789  123456789
		fprintf (xOutput, "\t    Low      ");// 13 = length of "Parametric CI"
//	                         1234567890123
		fprintf (xOutput, "\t   High  ");	// n = 9
//	                         123456789
		fprintf (xOutput, "\t    Low     ");// 12 = length of "Jackknife CI"
//	                         123456789012
		fprintf (xOutput, "\t   High  ");	// n = 9
//	                         123456789
	}
	fprintf (xOutput, "\n");
	free (method);
	return nCrit;
}
// --------------------------------------------------------------------------
void PrtTpfVal (FILE *output, float f, float fprime, float Ne,
				int len1, int len2)
// Print the first value f with len1 spaces, the next 2 with len2 spaces
{
	int i;
	fprintf (output, "\t%8.5f", f);
	for (i = 0; i < len1 - 8; i++) fprintf (output, " ");
	fprintf (output, "\t%8.5f", fprime);
	for (i = 0; i < len2 - 8; i++) fprintf (output, " ");
	fprintf (output, "\t%8.5f", Ne);
	for (i = 0; i < len2 - 8; i++) fprintf (output, " ");

}

// --------------------------------------------------------------------------
void PrtTpCI (FILE *output, float loNe, float hiNe, float infinite,
				int len1, int len2)
{
	int i;
	if (loNe < infinite)
		fprintf (output, "\t%8.1f", loNe);
	else
		fprintf (output, "\t%8s", "Infinite");
	for (i = 0; i < len1 - 8; i++) fprintf (output, " ");
	if (hiNe < infinite)
		fprintf (output, "\t%8.1f", hiNe);
	else
		fprintf (output, "\t%8s", "Infinite");
	for (i = 0; i < len2 - 8; i++) fprintf (output, " ");
}


// --------------------------------------------------------------------------

int GetLenTpPop (int nGeneration, char *popID[], char *abbrev, int *lenTpPop,
				char newInp)
// determine the length (*lenTpPop) of string to be used for population name
// in temporal, the name is a combination of sample names.
// Return value k could be bigger than *lenTpPop.
// *lenTpPop will be updated to k if k is bigger and this is at the start of
// a new input file (k remains bigger if at the middle of input file).
// At the calling function:
// If k < *lenTpPop, then from k to *lenTpPop, population name to be filled
// with blanks. If k > *lenTpPop, some chars will be trimmed.
{
	int k, m, maxLen;
	maxLen = 2*POP_TEMP+5;	// maximum length of popID is POP_TEMP
	k = 0;	// total length if take all sample names
	for (m=0; m<nGeneration; m++) k += strlen(popID [m]);
	// to include the number of slashes separating names:
	k += (nGeneration-1);
	// when newInp != 1, *lenTpPop won't be updated
	if (k > maxLen || (k > (*lenTpPop) && (newInp != 1))) {
	// only take the first and last names of samples:
		k = strlen(popID [0]) + strlen(popID [nGeneration-1]) + 5;
		// (If k still > *lenTpPop, pop name created by these two popID
		// at the calling function will be trimmed down.)
		*abbrev = 1;
	} else *abbrev = 0;	// k <= maxLen and (k <= *lenTpPop or newInp = 1).
	// Since *lenTpPop is always <= maxLen, this else implies that (as
	// *lenTpPop being updated by the following code) we'll have the return
	// value <= *lenTpPop. So, when *abbrev = 0; the population name
	// taken by all sample names will be fit to the length *lenTpPop.
	if (k > *lenTpPop && newInp == 1) *lenTpPop = k;
	return k;
}

//---------------------------------------------------------------------------
// This is called at PrtTpTabFile, after PrtTpHeader being called there.
// Jan 2017:
// add critOut as parameter, instead of local, which was given by PrtTpHeader
void PrtTpValues (FILE *xOutput, char newInp, int nGeneration,
					float timeline[], char *popID[], int g1, int g2,
					int stCrit, int nCrit, int critOut, float *critVal,
					int nPlan, int census, int lenTpPop, long *nIndAlle,
					char tempk, char tempc, char temps,
					float *fkmean, float *fcmean, float *fsmean,
					float *fkprimeMean, float *fcprimeMean, float *fsprimeAll,
					float *NeTempk, float *NeTempc, float *NeTemps,
					float *loNek, float *loNec, float *loNes,
					float *hiNek, float *hiNec, float *hiNes,
					float *jloNek, float *jloNec, float *jloNes,
					float *jhiNek, float *jhiNec, float *jhiNes,
					char param, char jack, float infinite,
					int lenGot, char abbrev, int lenPair,
					char *inpName, int lenInp, int nloci, int count,
					char common)
// g1,g2 are indices for generations
// count represents input file numbering, common = 1 if this xOutput is for
// all input files, with only one header to be printed.
// In this, count is the number of input so far.
// and smPair is the number of sample pairs so far.
// lenTpPop is the length for priting string representing populations.
// It should be set about 10 (length of the word "Population")before this
// function is called the first time, per input file
// lenPair is the maximum number of chars per sample name to print sample pair
{
// print pops
	char *popList;
	int i, k, m, n, n1, n2, lenSkip;
//	int critOut;
	// set size of popPair to be 2*(lenPair+2) to accomodate 2 pops,
	// together with string "<->".
	char *popPair = (char*) malloc(sizeof(char)*(2*(lenPair+2)));
	char *skipStr = (char*) malloc(sizeof(char)*200);
	*popPair = '\0';
	*skipStr = '\0';

	lenSkip = 0;
	if (common != 0) {	// print order and name of input file, and #loci
		// if the last parameter is 0, PrtPair just prints "count" blanks
		newInp = ((newInp == 1) && (g1+g2 == 1))? 1: 0;
		PrtPair (xOutput, count, inpName, lenInp, newInp);
		if (newInp == 1) fprintf (xOutput, "\t%6d\t", nloci);
		else fprintf (xOutput, "\t%6c\t", ' ');
		for (i=0; i<lenInp; i++) *(skipStr+i) = ' ';
		*(skipStr+lenInp) = '\t';
		lenSkip = lenInp + 1;
		for (i=0; i<6; i++) *(skipStr+lenSkip+i) = ' ';
		*(skipStr+lenSkip+6) = '\t';
		lenSkip += 7;
		*(skipStr+lenSkip) = '\0';
	};

	popList = (char*) malloc(sizeof(char)*(lenGot));
	*popList = '\0';
	if (g1+g2 == 1) {	// print population name by combining all sample names
		*popList = '\0';
		if (abbrev == 1) {	// only take the first and last names of samples:
			strcat (popList, popID [0]);
			strcat (popList, "/.../");
			strcat (popList, popID [nGeneration-1]);
		} else {
			for (m=0; m<nGeneration; ) {
				strcat (popList, popID [m]);
			// Choose a character to separate sample names (of the same pop):
//				if (++m < nGeneration) strcat (popList, "/");
//				if (++m < nGeneration) strcat (popList, "&");
				if (++m < nGeneration) strcat (popList, "+");
			};
		};
		k = lenTpPop - lenGot;
		// lenGot could be bigger than lenTpPop, i.e., k < 0.
		// In such case, trim popList to match the length (lenTpPop)
		if (k < 0) {	// take the first m chars, then put 3 dots,
		// then for the rest, take the last chars of popList:
			m = lenTpPop / 2;
			for (i=0; i<3; i++) *(popList+m+i) = '.';
			m += 3;
			n = lenTpPop - m;
			for (i=0; i<n; i++) *(popList+m+i) = *(popList+m+i-k);
			*(popList + lenTpPop) = '\0';
		};
		fprintf (xOutput, "%s", popList);
		for (m=0; m<k; m++) fprintf (xOutput, " ");
	} else for (i=0; i<lenTpPop; i++) fprintf (xOutput, " ");

	for (i=0; i<lenTpPop; i++) *(skipStr+(lenSkip+i)) = ' ';
	lenSkip += lenTpPop;
	*(skipStr+lenSkip) = '\0';

	if (nPlan > 1) {	// 13 chars for census size ("Census Size N")
	// When common = 1, only one plan for temporal methods, we only need
	// to print census size for each population at the first sample pair.
	// When common = 0, it is possible that both plans are run in the same
	// input file; so we need to print census size at every sample pair,
	// to distinguish with no census at plan II
		if (census > 0) {
			if (common == 0 || ((common==1) && (g1+g2==1)))
				fprintf (xOutput, "\t%10d   ", census);
			else fprintf (xOutput, "\t%13s", " ");
		} else fprintf (xOutput, "\t%13s", " ");
		*(skipStr+(lenSkip++)) = '\t';
		for (i=0; i<13; i++) *(skipStr+(lenSkip+i)) = ' ';
		lenSkip += 13;
		*(skipStr+lenSkip) = '\0';
	};
// print sample pair names:
	n1 = lenPair - strlen (popID[g1]);		// this is the number of blanks to be
			// added or the number of chars to be taken out from popID[g1]
	m = (n1 > 0)? 0: -n1;
	n2 = lenPair - strlen (popID[g2]);
	k = (n2 > 0)? 0: -n2;
	strcat(popPair, popID[g1]+m);
	strcat(popPair, "<->");
	strcat(popPair, popID[g2]+k);
	fprintf(xOutput, "\t");
	for (i = n1; i > 0; i--) fprintf (xOutput, "%c", ' ');
	fprintf (xOutput, "%s", popPair);
	for (i=0; i<n2; i++) fprintf (xOutput, "%c", ' ');
	n = 2*lenPair + 3;
	*(skipStr+(lenSkip++)) = '\t';
	for (i=0; i<n; i++) *(skipStr+(lenSkip+i)) = ' ';
	lenSkip += n;
	*(skipStr+lenSkip) = '\0';

// generations:
//	fprintf (xOutput, "\t%5.1f & %-5.1f", timeline[g1], timeline[g2]);
// This occupies the same as "Generations": the '&' is at 't' letter
	fprintf (xOutput, "\t%5.1f &%4.1f", timeline[g1], timeline[g2]);
	*(skipStr+(lenSkip++)) = '\t';
	for (i=0; i<lenPair; i++) *(skipStr+(lenSkip+i)) = ' ';
	lenSkip += lenPair;
	*(skipStr+lenSkip) = '\0';

	// with critical values
	// Jan 2017: cover the next line, critOut already in parameter
	// critOut = nCrit - stCrit;

	for (n=stCrit; n<nCrit; n++) {
// Jan 2017: no spec pcrit in Tp
		if (critVal[n] > 0 && critVal[n] <= PCRITX) continue;
		if (n > stCrit) fprintf (xOutput, "%s", skipStr); // 6=length("PCrit.")
		if (critOut > 1) fprintf (xOutput, "\t%6.4f", *(critVal+n));
		// independent alleles:
		fprintf (xOutput, "\t%10lu  ", *(nIndAlle+n));
		// reserve m spaces for printing those f, fprime, Ne, following a Tab
		m = 9;	// length for printing the last 2 values for each method,
				// same length as their headers (see PrtTpHeader)
		if (tempk==1) {
			PrtTpfVal (xOutput,
						*(fkmean+n), *(fkprimeMean+n), *(NeTempk+n), 8, m);
			if (param==1) PrtTpCI (xOutput, *(loNek+n), *(hiNek+n), infinite, 13, 9);
			if (jack==1) PrtTpCI (xOutput, *(jloNek+n), *(jhiNek+n), infinite, 12, 9);
		}
		if (tempc==1) {
			PrtTpfVal (xOutput,
						*(fcmean+n), *(fcprimeMean+n), *(NeTempc+n), 10, m);
			if (param==1) PrtTpCI (xOutput, *(loNec+n), *(hiNec+n), infinite, 13, 9);
			if (jack==1) PrtTpCI (xOutput, *(jloNec+n), *(jhiNec+n), infinite, 12, 9);
		}
		if (temps==1) {
			PrtTpfVal (xOutput,
						*(fsmean+n), *(fsprimeAll+n), *(NeTemps+n), 11, m);
			if (param==1) PrtTpCI (xOutput, *(loNes+n), *(hiNes+n), infinite, 13, 9);
			if (jack==1) PrtTpCI (xOutput, *(jloNes+n), *(jhiNes+n), infinite, 12, 9);
		}
		// ------------------------------------------------------------------
		fprintf (xOutput, "\n");
	};	// end of for n=stCrit ...
	free (popList);
	free(popPair);
	free(skipStr);

}

// --------------------------------------------------------------------------
// Jan 2017: modified according to the modifications of PrtTpHeader
void PrtTpTabFile (FILE *xOutput, int popRead, int popStart, char newSet,
					char lastSet, float *critVal, int nCrit, int g1, int g2,
					int generation, int nGeneration, float timeline[], long *nIndAlle,
					float *Hkmean, float *Hcmean, float *Hsmean,
					float *fkmean, float *fcmean, float *fsmean,
					float *fkprimeMean, float *fcprimeMean, float *fsprimeAll,
					float *NeTempk, float *NeTempc, float *NeTemps,
					float *loNek, float *loNec, float *loNes,
					float *hiNek, float *hiNec, float *hiNes,
					float *jloNek, float *jloNec, float *jloNes,
					float *jhiNek, float *jhiNec, float *jhiNes,
					char param, char jack,
					float infinite, char tempk, char tempc, char temps,
					char *popID[], int topCrit, int nPlan, int census,
// add May 2013 to be able to print input file name and loci as columns,
// and not to print header when necessary.
// count represents input file numbering, common = 1 if this xOutput is for
// all input files, with only one header to be printed.
// In this, count is the number of input so far.
// and smPair is the number of sample pairs so far.
					char *inpName, int nloci, int count, int smPair,
					char common, int *lenTpPop)
// print to extra ouput (shorter form) file.
{
	char newInp;
	// Jan 2017: add critOut as param in PrtTpHeader, used by PrtTpValues
	int stCrit, critOut;
	char abbrev;
	int lenGot;
	int lenPair = 10;	// maximum #chars to print sample name for sample pair
	int lenInp = 19;	// maximum number of chars to print input file name
	if (xOutput == NULL) return;
	newInp = (popRead == popStart-1+generation)? 1: 0;
	if (g1+g2 == 1)
	// At the calling of this function, g1 < g2 are indices running from 0
	// to at most (nGeneration-1) [may not have enough samples left to reach].
	// The condition (g1+g2 == 1) means that we are at the first sample pair.
	// (First sample is at generation corresponding to g1, 2nd sample for g2.)
		lenGot = GetLenTpPop (generation, popID, &abbrev, lenTpPop, newInp);
	// when newInp = 1, this function is first called with g1+g2 = 1, so
	// lenTpPop might be reassigned by GetLenTpPop. PrtTpHeader will use this
	// lenTpPop value.

	// Parameter common = 1 when running multiple files and this output file
	// is for all, When common = 0, this output file is for each input file
	// whether the program is run for a single or multiple input files.
	// Therefore, header should be written only once for all input files when
	// common = 1, and written once for each input file when common = 0.
	// Parameter smPair represents the total number of sample pairs so far,
	// in running single or multiple input files (for both common = 0 and 1).
	// The condition smPair = 1 implies that newInp = 1 and g1+g2 = 1.
	// (The first sample happens at the beginning of the first input file.)
	// The condition [newInp = 1 and g1+g2 = 1] means a first pair at a new
	// input file (which may not be the first input file), so header is
	// written only if common = 0.
	//
	// The commented out conditions here represent the above discussion, and
	// are the same as the ones in effect.
	// if ((smPair == 1) || (common == 0 && newInp == 1 && g1+g2 == 1))
// Jan 017: add stCrit, critOut in the param list for PrtTpHeader
//	stCrit = 0;
	if ((newInp == 1 && g1+g2 == 1) && (common == 0 || smPair == 1))
		nCrit = PrtTpHeader (xOutput, critVal, nCrit, topCrit,
							&stCrit, &critOut,
					nPlan, tempk, tempc, temps, param, jack, lenInp,
					lenPair, *lenTpPop, common, nGeneration);
// Jan 2017:
// stCrit, nCrit, critOut need to be modified when spec Pcrit is present
	else
		nCrit = CritEndLine (critVal, nCrit, topCrit, &stCrit, &critOut, 0);
	if (nCrit == 0) return;	// something wrong
// now, print the values:
// Cover 3 lines since those are already calculated from PrtTpHeader
//	stCrit = 0;
//	if (topCrit > 0 && topCrit < nCrit) nCrit = topCrit;
//	else if (topCrit == 0) stCrit = nCrit - 1;

	// in PrtTpValues, parameter lenGot is used only when (g1+g2 == 1)
	PrtTpValues (xOutput, newInp, generation, timeline, popID, g1, g2,
				stCrit, nCrit, critOut, critVal, nPlan, census, *lenTpPop,
				nIndAlle, tempk, tempc, temps, fkmean, fcmean, fsmean,
				fkprimeMean, fcprimeMean, fsprimeAll, NeTempk, NeTempc,
				NeTemps, loNek, loNec, loNes, hiNek, hiNec, hiNes,
				jloNek, jloNec, jloNes, jhiNek, jhiNec, jhiNes, param, jack,
				infinite, lenGot, abbrev, lenPair,
				inpName, lenInp, nloci, count, common);

	fflush (xOutput);
}


// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
// RunPop
// --------------------------------------------------------------------------


int RunPop0 (int icount, char *inpName, FILE *input, char append, FILE *output,
			char *outFolder, struct locusMap *locList,	// added in Nov 2014
			FILE *outLoc, char *outLocName, FILE *outBurr,
			char *outBurrName, FILE *shOutputLD, FILE *shOutputHet,
			FILE *shOutputCoan, FILE *shOutputTemp, int popLoc1, int popLoc2,
			int popBurr1, int popBurr2, int topBCrit, int popStart,
			int popEnd, int maxSamp, int nloci, int lenM, int maxMobilVal,
			int nCrit, float critVal[], char format, char param,
			char jacknife, char *locUse, char mating, char *missFileName,
			float infinite, int lenBlock,
			char mLD, char mHet, char mNomura, char mTemporal,
			int nGeneration, float timeline[], AGEPTR *ageSeq,
			int nSeq, char getAge, int tempClue, int tempxClue, char byRange,
			int topCrit, int nPlan, int census, int *totPop, int *totPairTmp,
			char common, char tabX,
// add parameters in Apr 2015:
			char sepBurOut, char moreCol, char BurAlePair,
			struct chromosome *chromoList, int nChromo, int chroGrp)
// Return values
// * 0: things are OK, everything else is error.
// * 1,2: serious error in genotype data: either nondigits are present or
//		too many digits.
// * 3: Data for a sample don't cover all loci.
// * 4: no input file, or no method at all.
// * 5: all loci are dropped as designated by option, no loci to run
// * -1: OUT OF MEMORY.
// If nonzero is returned, the running stops at the population
// where the error occurs.
// The first parameter icount is used only for printing Input #,
// to distinguish if this function is called in doing multiple files
{

	time_t rawtime;
// to store all mobility values and their occurences
	ALLEPTR *alleList;
// to store the number of different mobility values,
// and missing data information:
	int *nMobil;
// *(missptr+p) = number of samples missing data at locus (p+1)
	int *missptr;
// to store all samples:
	FISHPTR *fishHead = NULL, *fishTail = NULL;
// Variables for calling Loc_Freq:
// to store min and max values of frequencies at each locus:
	float *minFreq, *maxFreq;
// to determine among loci not deleted from input, which one will
// satisfy criteria for consideration
	char *okLoc;
	int nLocOK;
// calculation of exp r2, estimated Ne,  harmonic mean,
// r^2 from Burrows coeffs, indep. of alleles
	float *wExpR2, *estNe, *wHarmonic, *rB2WAve, *r2Drift;
	char bigInd = 0;	// to notify if # ind alle. too big, to adjust output
	int *sampData;
// confidence intervals:
	float *confJacklow, *confJackhi, *confParalow, *confParahi;
// pop names:
	char *popID, *newID;
	char *jackOK;
	char weighsmp;
	// added in june 2016 for jackknife on sample, to notice if the number
	// of individuals in a population sample is less than MINSAMP,
	// no jackknife on samples. The condition for having jackknife on samples
	// depends also on the population since different populations may have
	// different number of individuals; so jackknife is ok for one but not
	// ok for others
	char jSamp;

	char moreDat, moreBurr, moreBurr0;	// to determine which population
		// that outputs to auxiliary files outLoc, outBurr will stop.
	char header;
	char genErr[GENLEN];
	int firstErr;	// for the first locus with missing data;
//-------------------------------------------------------------------
// add in Nov 2014 for having separate Burrows output files:
	char *outFile = (char *) malloc(PATHFILE * sizeof(char));
	*outFile = '\0';
//-------------------------------------------------------------------

// Coancestry:
	float loNbCoan, hiNbCoan;
//-------------------------------------------------------------------
// To deal with heterozygote excess method
	float hetNeb;
	float *hSamp;
	float *hetD, *estHetN;
//	float *jloHetNe, *jhiHetNe;
	float *loHetNe, *hiHetNe;
	long *indAlleH;
// To deal with Coancestry method:
	float f1, coanNeb;
	float hSamCoan;
//-------------------------------------------------------------------
// To deal with temporal method:
	FREQPTR *freqList;
	long *nTotAlle, *nIndAlle;
	float *NeTempk, *NeTempc, *NeTemps;
	float *loNek, *hiNek, *loNec, *hiNec, *loNes, *hiNes;
	float *jloNek, *jhiNek, *jloNec, *jhiNec, *jloNes, *jhiNes;
	float *Hkmean, *Hcmean, *Hsmean;	// harmonic sample size
	float *fkmean, *fcmean, *fsmean;
	float *fkprimeMean, *fcprimeMean, *fsprimeAll;
	int g1, g2;
	int errfreq = 0;
// these will be set acording to input tempClue:
	char tempk = 1, tempc = 1, temps = 1;
// this will be set according to tempxClue and tempk, tempc, temps,
// used to print to extra files
	char tempkx = 1, tempcx = 1, tempsx = 1;
	char *popIDtemp [MAXGENERATION];
	int *popSize;	// for storing sample sizes in generations
	int generation = 0;
//	int readAge = 1;	// > 0: to send notice that age reading happens
						// when 0, no longer reading ages for the rest.
	char newSet = 1;	// to determine if generation set is a new one
	char lastSet= 0;	// to determine if generation set is the last one
	int seq = 0;
	int nPoptemp = 0, nPairTmp = 0;

// add July 2013
	int lenTpPop = 25;	// to accomodate 2 population names, each with 10
						// chars, and a string "/.../"
//-------------------------------------------------------------------
	int popRun = 0;
	int p, m, n, samp, ind;
	int isize, size;	// use for missptr array
	int lastOK, next, popRead, err;
	int nErr;		// count total missing data
	int nSampErr;	// count number of samples having missing data
	int noGen;		// count number of genotype missing in a sample
	int errCode = 0;
	FILE *missDat = NULL;
	int memOut = 0;
	int nLocUsed;	// for upper bound of locus numbering
	char makeFish;	// to determine if fishList is to be created
	if (input == NULL) return 4;	// no input file, don't do a thing
	if (mLD+mHet+mNomura+mTemporal == 0) return 4;
//	if ((nLocUsed = PrtLocInfo (output, locUse, nloci)) == 0) {
	if (mTemporal == 1 && mLD+mHet+mNomura == 0)
		nLocUsed = PrtLimitUse (output, locUse, nloci, byRange,
							popStart, popEnd, MAX_POP, maxSamp, "Sample");
	else nLocUsed = PrtLimitUse (output, locUse, nloci, byRange,
							popStart, popEnd, MAX_POP, maxSamp, "Population");
	if (nLocUsed == 0) {
		fprintf (output, "No loci to run!\n");
		printf ("No loci to run!\n");
		return 5;
	};
	if (nLocUsed == 1) mLD = 0;
// cover this since jackknife is on samples:
//	if (nLocUsed < 3) jacknife = 0;

// Note that array locUse[p] = 1 if locus (p+1) is used, 0 otherwise.
// Add Jan 2013:
// Therefore, we can reassign nlocus to be the largest locus that is used,
// instead of going all over nloci. Note that after the last locus is read,
// the cursor moves to the next line. So, in order for this added code to be
// correct, an individual should only occupy no more than one line of input.
// (This was not assumed in function GetSample but it is followed in
// GENEPOP and FSTAT formats.) When an individual occupies more than one line,
// and the reassigned nloci is not the locus on the last data line of that
// individual, then the next data reading will still read data for that
// individual, not data of the next individual. Such case will cause an error.
// If that possibility can happen, just comment out the next 3 lines of code.

// In case that one individual occupies only one line, these lines of code
// will help in limiting storages. Also, if genotypes at some locus after
// the largest acceptable locus is deleted (intentionally or not), stopping
// before GetSample tries to catch all declared loci is good, because
// GetSample may go to another line if it cannot get all declared loci!
// (If that happens, data will be screwed up!)
	p = nloci;
	while (*(locUse+p-1) == 0) p--;
	nloci = p;

// holds data for one sample:
	if ((sampData = (int*) malloc(sizeof(int)*nloci*2)) == NULL) {;
		printf ("Out of memory for storing sample data!\n");
		return -1;
	};
	size = nloci;
//	size = (size*(size+1))/2;
	if ((missptr = (int*) malloc(sizeof(int)*size)) == NULL) {
//	if ((missptr = (int*) malloc(sizeof(long)*size)) == NULL) {
		free (sampData);
		printf ("Out of memory for missing data array!\n");
		return -1;
	};
	if ((nMobil = (int*) malloc(sizeof(int)*nloci)) == NULL) {
		free (sampData);
		free (missptr);
		printf ("Out of memory for array to count alleles!\n");
		return -1;
	};
	if ((minFreq = (float*) malloc(sizeof(float)*nloci)) == NULL) {
		free (sampData);
		free (missptr);
		free (nMobil);
		printf ("Out of memory for min freq array!\n");
		return -1;
	};
	if ((maxFreq = (float*) malloc(sizeof(float)*nloci)) == NULL) {
		free (sampData);
		free (missptr);
		free (nMobil);
		free (minFreq);
		printf ("Out of memory for max freq array!\n");
		return -1;
	};
	if ((okLoc = (char*) malloc(sizeof(char)*nloci)) == NULL) {
		free (sampData);
		free (missptr);
		free (nMobil);
		free (minFreq);
		free (maxFreq);
		printf ("Out of memory for evaluating loci array!\n");
		return -1;
	};
	if ((alleList = (ALLEPTR*) malloc(sizeof(ALLEPTR)*nloci)) == NULL) {
		free (sampData);
		free (missptr);
		free (nMobil);
		free (minFreq);
		free (maxFreq);
		free (okLoc);
		printf ("Out of memory to reserve alleles!\n");
		return -1;
	};

//----------- for temporal method
	if (mTemporal == 1) {
		if ((freqList = (FREQPTR*) calloc(nloci, sizeof(FREQPTR))) == NULL)
		{
			free (sampData);
			free (missptr);
			free (nMobil);
			free (minFreq);
			free (maxFreq);
			free (okLoc);
			free (alleList);
			printf  ("Out of memory to reserve alleles!\n");
			return -1;
		};
		popSize = (int*) malloc(sizeof(int)*MAXGENERATION);
		for (n=0; n<MAXGENERATION; n++)
			popIDtemp [n] = (char*) malloc(sizeof(char)*POP_TEMP);
		if (getAge != 0) lastSet = 0;
		else lastSet = 1;
// reset tempk, tempc, temps based on values of tempClue. If tempClue
// is between 1 and 6, write it in binary, tempClue = XYZ: X, Y, Z in {0,1}.
// (One of them must be 0 since tempClue < 7).
// Assign: temps = X, tempc = Y, tempk = Z.
		if (tempClue > 0 && tempClue < 7) {
			temps = tempClue/4;
			tempClue -= 4*temps;
			tempc = tempClue/2;
			tempk = tempClue - tempc*2;
		};
	// tempk, tempc, temps are to tell if temporal methods: Pollak(tempk),
	// Nei/Tajima(tempc), or Jorde/Ryman(temps) are to be calculated.
	// tempkx, tempcx, tempsx are for outputting to extra output file.
	// Note that all these are in use only when mTemporal = 1.
	// If tempk = 0, we should not have tempkx = 1, similar for others.
	// First, set values of temp#x the same way as temp#, using tempxClue:
		if (tempxClue > 0 && tempxClue < 7) {
			tempsx = tempxClue/4;
			tempxClue -= 4*tempsx;
			tempcx = tempxClue/2;
			tempkx = tempxClue - tempcx*2;
		};
	// In case that the ones temp#x are zero at those temp# that are 1
	// that is, in principle, no temporal is outputted in extra file,
	// then we should consider this as an error, and should go to default,
	// which is: temp#x = temp#. This case happen only when n = 0:
		n = tempk*tempkx + tempc*tempcx + temps*tempsx;
		if (n == 0) {
			tempkx = tempk; tempcx = tempc; tempsx = temps;
		} else {
		// Note that temp#x should not be 1 if temp# is 0 (if there is
		// no method corresponding to temp#, then should not have output):
			tempkx *= tempk;
			tempcx *= tempc;
			tempsx *= temps;
		}
	};
	// set = 1 to fill in list of fish, 0 if list is unneeded
	makeFish = (mLD > 0 || mNomura > 0 || (mHet > 0 && nCrit > 1))? 1: 0;
	if (makeFish > 0) {
		if ((fishHead = (FISHPTR*) malloc(sizeof(FISHPTR)*nloci)) == NULL) {
//	if ((fishHead = (FISHPTR*) calloc(nloci, sizeof(FISHPTR))) == NULL) {	//7
			free (sampData);
			free (missptr);
			free (nMobil);
			free (minFreq);
			free (maxFreq);
			free (okLoc);
			free (alleList);
			if (mTemporal == 1) free (freqList);
			printf ("Out of memory for sample list!\n");
			return -1;
		};
		if ((fishTail = (FISHPTR*) malloc(sizeof(FISHPTR)*nloci)) == NULL) {
//	if ((fishTail = (FISHPTR*) calloc(nloci, sizeof(FISHPTR))) == NULL) {	//8
			free (sampData);
			free (missptr);
			free (nMobil);
			free (minFreq);
			free (maxFreq);
			free (okLoc);
			free (alleList);
			if (mTemporal == 1) free (freqList);
			free (fishHead);
			printf ("Out of memory for sample list!\n");
			return -1;
		};
	};

//---------------------------------------------------------------
// lenBlock = LEN_BLOCK is small, so those are not likely to fail
	popID = (char*) malloc(sizeof(char)*lenBlock);
	newID = (char*) malloc(sizeof(char)*lenBlock);
	*popID = '\0';
// nCrit is small, those are not likely to fail:
	wExpR2 = (float*) malloc(sizeof(float)*nCrit);
	estNe = (float*) malloc(sizeof(float)*nCrit);
	wHarmonic = (float*) malloc(sizeof(float)*nCrit);
	rB2WAve = (float*) malloc(sizeof(float)*nCrit);
	r2Drift = (float*) malloc(sizeof(float)*nCrit);
	double *nIndSum = (double*) malloc(sizeof(double)*nCrit);
	jackOK = (char*) malloc(sizeof(char)*nCrit);
	for (n=0; n<nCrit; *(jackOK+n)=1, n++);
// for het excess method
	hSamp = (float*) malloc(sizeof(float)*nCrit);
	estHetN = (float*) malloc(sizeof(float)*nCrit);
	hetD = (float*) malloc(sizeof(float)*nCrit);
	loHetNe = (float*) malloc(sizeof(float)*nCrit);
	hiHetNe = (float*) malloc(sizeof(float)*nCrit);
//	jloHetNe = (float*) malloc(sizeof(float)*nCrit);
//	jhiHetNe = (float*) malloc(sizeof(float)*nCrit);
	indAlleH = (long*) malloc(sizeof(long)*nCrit);
// for Temporal method:
	nTotAlle = (long*) malloc(sizeof(long)*nCrit);
	nIndAlle = (long*) malloc(sizeof(long)*nCrit);
	NeTempk = (float*) malloc(sizeof(float)*nCrit);
	NeTempc = (float*) malloc(sizeof(float)*nCrit);
	NeTemps = (float*) malloc(sizeof(float)*nCrit);

	loNek = (float*) malloc(sizeof(float)*nCrit);
	hiNek = (float*) malloc(sizeof(float)*nCrit);
	loNec = (float*) malloc(sizeof(float)*nCrit);
	hiNec = (float*) malloc(sizeof(float)*nCrit);
	loNes = (float*) malloc(sizeof(float)*nCrit);
	hiNes = (float*) malloc(sizeof(float)*nCrit);

	jloNek = (float*) malloc(sizeof(float)*nCrit);
	jhiNek = (float*) malloc(sizeof(float)*nCrit);
	jloNec = (float*) malloc(sizeof(float)*nCrit);
	jhiNec = (float*) malloc(sizeof(float)*nCrit);
	jloNes = (float*) malloc(sizeof(float)*nCrit);
	jhiNes = (float*) malloc(sizeof(float)*nCrit);

	Hkmean = (float*) malloc(sizeof(float)*nCrit);
	Hcmean = (float*) malloc(sizeof(float)*nCrit);
	Hsmean = (float*) malloc(sizeof(float)*nCrit);
	fkmean = (float*) malloc(sizeof(float)*nCrit);
	fcmean = (float*) malloc(sizeof(float)*nCrit);
	fsmean = (float*) malloc(sizeof(float)*nCrit);
	fkprimeMean = (float*) malloc(sizeof(float)*nCrit);
	fcprimeMean = (float*) malloc(sizeof(float)*nCrit);
	fsprimeAll = (float*) malloc(sizeof(float)*nCrit);
// for Jacknife and Parametric confidence intervals
	confJacklow = (float*) malloc(sizeof(float)*nCrit);
	confJackhi = (float*) malloc(sizeof(float)*nCrit);
	confParalow = (float*) malloc(sizeof(float)*nCrit);
	confParahi = (float*) malloc(sizeof(float)*nCrit);
// add in Aug 2016 for storing degree of freedom in Jackknife
	long *Jdegree = (long*) malloc(sizeof(long)*nCrit);
	for (n=0; n<nCrit; *(confJacklow+n)=*(confParalow+n)=-(float)infinite,
					*(confJackhi+n)=*(confParahi+n)=(float)infinite, n++);
// Jan 2017, to deal with spec Pcrit (rejecting singleton alleles)
// This part is added so that we know if this spec Pcrit exists and
// then can print explanation for the symbol used to identify this
// spec Pcrit at the top of the main output (in the calling PrtPop):
	char specP = 0;
	for (n=0; n<nCrit; n++) {
		if (critVal[n] > 0 && critVal[n] <= PCRITX) {
			specP = 1;
			break;
		}
	}

// temporarily add for checking:
//char opened = 0;

	next = 0;
	ind = 0;
	popRead = 0;
	nSampErr = 0;
	nErr = 0;
	for (; next != -1 && popRead <= popEnd; ) {
		strcpy (popID, newID);
		if (format == FSTAT) next = DatPopID (input, newID, lenBlock);
		else next = GenPopID (input, "pop", newID, lenBlock);
		if (next != 0) {
			// either go to next pop (next = 1) or end of file (next = -1),
			// so all data of this pop were read, need to do calculations.
			PrtSumMisDat (missDat, popRead, nErr, newID, next);
			weighsmp = (nSampErr > 0? 1: 0);
			nSampErr = 0;
			nErr = 0;
			if (popRead >= popStart) {
//			if (next == -1 || popRead >= popStart) {
				printf ("-> Total samples = %d", samp);
				if (weighsmp > 0) printf (", with data missing");
				printf ("\n");
				// add in June 2016 for jSamp
				jSamp = (jacknife == 1 && samp >= MINSAMP)? 1: 0;
				PrtPop (output, popRead, popID, samp, mLD, mHet, mNomura, mating,
						nloci, nMobil, locUse, specP);
				if (mHet+mLD > 0) PrtFreq (output, mLD, critVal, nCrit, '-', '-');
				for (p=0; p<nloci; p++) if (*(nMobil+p)==0) *(okLoc+p)=0;
				moreDat = (popRead >= popLoc1 && popRead <= popLoc2)? 1: 0;
				moreBurr0 = (popRead >= popBurr1 && popRead <= popBurr2)? 1: 0;
				moreBurr = moreBurr0;
				Loc_Freq (alleList, nloci, samp, &hetNeb, nMobil, missptr,
						locUse, minFreq, maxFreq, outLoc, outLocName,
						moreDat, popRead, mHet, lenM, locList);
//	change in Nov 2014/Jan 2015 with sepBurOut
				if (outBurr != NULL && moreBurr == 1 && sepBurOut == 0)
					fprintf (outBurr,
					"\nPOPULATION%6d\t(Sample Size = %d)\n", popRead, samp);
				for (n=0; n<nCrit; n++) {// for nCrit frequency cut-off values
				// this loop is for LD and HetExcess methods only
					if (mHet + mLD == 0) break;
// print to console if not running multiple files (in all others, icount=0):
//					if (icount == 0) {
						if (critVal[n] == 0)
							printf ("   * For lowest freq: %5s\n", "0+");
						else
							printf ("   * For lowest freq: %5.3f\n", critVal[n]);
//					};
				// Burrow outputs only to topBCrit highest crit. values
				// If topBCrit = 0: only critical value 0 (at n = nCrit-1)
				// If topBCrit > 0: the crit. values taken are:
				// critVal[0], ..., critVal[m], m = min(topBCrit, nCrit) - 1.
				// If topBCrit < 0: all
					moreBurr = 0;
					if ((topBCrit < 0) || (topBCrit-n > 0) ||
						((topBCrit == 0) && (n == nCrit-1)))
						moreBurr = moreBurr0;
// add in Nov 2014/ Jan 2015:
					if (moreBurr == 1 && sepBurOut == 1) {
						outBurr = NULL;
						*outBurrName = '\0';
						if (NONAMEBUR != 1)
							GetPrefix (inpName, outBurrName, LENFILE-20, PATHCHR);
						GetBurrName(outBurrName, popRead, critVal[n]);
						*outFile = '\0';
						outFile = strcat (outFile, outFolder);
						if ((outBurr=fopen(strcat(outFile, outBurrName),"w"))!=NULL)
						{
							if (NOEXPLAIN != 1) {
								PrtVersion (outBurr);
								fprintf (outBurr, "Input File: %s\n\n", inpName);
								fprintf (outBurr,
									"\nPOPULATION%6d\t(Sample Size = %d)\n",
										popRead, samp);
							}
						} else popBurr2 = 0;	// if cannot open, stop Burrows

					}

					nLocOK = Loci_Eligible (samp, missptr, critVal[n],
							alleList, nloci, nMobil, minFreq, maxFreq,
							okLoc, &lastOK, locUse, outLoc, outBurr,
							moreDat, moreBurr, sepBurOut, moreCol);
					*(jackOK+n) = (nLocOK <= MAXJACKLD)? 1: 0;
					if (mLD == 1) {
						memOut = 0;
						estNe[n] = LDmethod (critVal[n], alleList, popRead,
								samp, fishHead, nMobil, missptr, lastOK, okLoc,
								(nIndSum+n), (rB2WAve+n), (r2Drift+n), (wHarmonic+n),
								(wExpR2+n), outBurr, outLoc, moreDat, moreBurr,
								outBurrName, mating, infinite, param, jSamp,
								(jackOK+n),	// to handle if cannot do jackknife
								(confJacklow+n), (confJackhi+n), (Jdegree+n),
								(confParalow+n), (confParahi+n), weighsmp,
								&memOut, icount, sepBurOut, moreCol, BurAlePair,
// temporarily add for checking with checkR2:
//&opened,
								chromoList, nChromo, chroGrp);
						if (*(nIndSum+n) >= infinite) bigInd = 1;
					// add in Nov 2014/ Jan 2015:
						if (outBurr != NULL && sepBurOut == 1) {
							printf ( "     Burrows coeffs are in file %s.\n",
																outBurrName);
							fclose (outBurr);
							outBurr = NULL;
						}
					};	// end of LD method
				// Dec 2016: no dropping only singletons in Het method:
					if (critVal[n] > 0 && critVal[n] <= PCRITX) continue;
					if (mHet == 1)
						HetXcess (fishHead, alleList, nloci, samp, nMobil, missptr,
							okLoc, outLoc, moreDat, critVal[n], (hetD+n),
							(estHetN+n), (indAlleH+n), (hSamp+n), (loHetNe+n),
							(hiHetNe+n), param);
				};	// end of loop for critical values
// temporarily add for checking:
//if (opened != 0) opened = popRead+1;

				if (mLD == 1 && memOut == 0) {
					PrtLDResults (output, nCrit, wHarmonic, nIndSum,
								rB2WAve, wExpR2, estNe, infinite, bigInd);

// To inform that recalculating Ne when missing data is suppressed:
					if (RESETNE == 0 && weighsmp > 0) fprintf(output,
					"(No attempt to adjust r^2 and Ne for missing data.)\n");

//					header = param+jacknife;
					header = param+jSamp;
					if (param==1)
						PrtLDConfid (output, nCrit, confParalow,
									confParahi, infinite, 0, &header, jackOK, bigInd);
				// if previous call proceeded, header will become false
//					if (jacknife==1)
					if (jSamp==1)
						PrtLDConfid (output, nCrit, confJacklow,
									confJackhi, infinite, 1, &header, jackOK, bigInd);
					if (tabX==0) PrtLDxFile (inpName, shOutputLD, samp, wHarmonic,
							popRead, popStart, popID, critVal, nCrit, nIndSum,
							rB2WAve, wExpR2, estNe, param, jacknife,
							infinite, confParalow, confParahi, confJacklow,
							confJackhi, Jdegree, jackOK, mating, topCrit, nLocUsed, icount, common);
					else PrtLDTabFile (inpName, shOutputLD, samp, wHarmonic,
							popRead, popStart, popID, critVal, nCrit, nIndSum,
							rB2WAve, wExpR2, estNe, param, jacknife,
							infinite, confParalow, confParahi, confJacklow,
							confJackhi, jackOK, mating, topCrit, nLocUsed, icount,
							common, Jdegree);
				};	// end of printing LD results
				// print to main output, Neb from Het Ex. and Coancestry,
				// From Het Ex method, based on lowest freq. req.

				if (mHet == 1) {
					PrtHetNe (output, hetD, estHetN, loHetNe, hiHetNe, hSamp,
								param, nCrit, critVal, indAlleH, infinite);
					if (tabX==0) PrtHetxFile (inpName, shOutputHet, popRead,
								popStart, popID, critVal, nCrit, indAlleH, hetD,
								estHetN, param, infinite, loHetNe, hiHetNe,
								samp, hSamp, topCrit, nLocUsed, icount, common);
					else PrtHetTabFile (inpName, shOutputHet, popRead,
								popStart, popID, critVal, nCrit, indAlleH, hetD,
								estHetN, param, infinite, loHetNe, hiHetNe,
								samp, hSamp, topCrit, nLocUsed, icount, common);
				};

				for (p=0; p<nloci; p++) *(okLoc+p) = *(locUse+p);
				if (mNomura == 1) {
					coanNeb = CoanMethod (fishHead, alleList, nMobil,
								nloci, samp, okLoc, &f1, outLoc, moreDat,
								&loNbCoan, &hiNbCoan, jacknife, missptr, &hSamCoan);
					n = (mHet+mLD > 0)? nCrit: 1;
					m = (mHet+mLD > 0 && critVal[nCrit-1] == 0)? 0: 1;
					PrtNomuraNe (output, f1, coanNeb, n, m, loNbCoan, hiNbCoan,
								jacknife, hSamCoan);
					if (tabX==0) PrtCoanxFile (inpName, shOutputCoan, popRead,
								popStart, popID, f1, coanNeb, jacknife,
								infinite, loNbCoan, hiNbCoan, samp, hSamCoan,
								nLocUsed, icount, common);
					else PrtCoanTabFile (inpName, shOutputCoan, popRead,
								popStart, popID, f1, coanNeb, jacknife,
								infinite, loNbCoan, hiNbCoan, samp, hSamCoan,
								nLocUsed, icount, common);
				};
				n = 26 + 12*nCrit;	// draw a line of asterisks at the end:
				if (bigInd == 1) n += 2*nCrit;
				if (mHet+mLD > 0) PrtLines (output, n, '*');
// show the time after finishing Nomura or LD with big data:
				if (((mLD > 0 && (nloci > 4000 || samp >= 10000)) ||
					(mNomura > 0 && (nloci*samp >= 10000)))
					&& (next != -1) && (popRead < popEnd)) {
					time ( &rawtime );
					fprintf (output, "\nTime: %s\n", ctime (&rawtime));
					fflush (output);
				};
// for temporal method ----------------------------------------------------
				if (mTemporal == 1) {
					if (generation == 0) nPoptemp++;
					if (generation == 0 && getAge != 0 && seq < nSeq) {
					// every time GetGeneration is called, a generation set
					// taken from the list ageSeq is recorded into timeline,
					// then those values are deleted from the list. At the
					// end, the list may not be empty because some generation
					// sets read (by InfoDirective) may not be used (number
					// of population samples to be run is smaller than the
					// number of generations read). However, list ageSeq is
					// constructed only when this program runs with a single
					// input file, so function RunPop0 is called only once
					// then the program ends. Hence we need not worry about
					// memory not deallocated. If needed (when we want to
					// have several generations sets to be read for each
					// input file in running multiple input files), then we
					// can add code to deallocate the list after seq = nSeq.
						newSet = 1;
						seq++;
						GetGeneration (ageSeq, timeline, &nGeneration, &census);
//fprintf (output, "\nSet %d of %d sets, Census size = %d,\tGenerations: ", seq, nSeq, census);
//for (n=0; n<nGeneration; n++) fprintf(output, "%8.1f", timeline[n]);
//fprintf(output, "\n\n");
						if (seq == nSeq) lastSet = 1;
						else lastSet = 0;
					} else {
					// since PrtTempxFile accepts newSet when generation
					// takes value > 0, so must set condition below so that
					// newSet still = 1 by assignment in the if above.
					// Moreover, when generations already used, and no more
					// to be read (the if above fails, no new generations),
					// then set newSet = 0.
						if (generation == 0) newSet = 0;
						if (getAge == 0) {	// no sequence of generations
						// so there is only one generation set, it is not
						// considered new when popRead surpasses the number
						// of populations (samples) corresponding to the number
						// of generations. Note that PrtTempxFile is called
						// after sample (population) for the last generation
						// is read.
//							if (popRead > nGeneration) newSet = 0;
							if ((popRead-popStart+1) > nGeneration) newSet = 0;
							else newSet = 1;
						};
					};

				// Use strcpy: only copy a maximum of POP_TEMP characters
				// from popID to popIDtemp, and choose the last ones, which
				// are more likely to represent the population.
					n = strlen (popID);
					m = (n > POP_TEMP)? n - POP_TEMP: 0;
					strcpy (popIDtemp [generation], popID+m);
					*(popSize + generation) = samp;
					AddFreqWide (freqList, alleList, nloci, samp, missptr,
								locUse, nGeneration, generation, &errfreq, weighsmp);
					if (errfreq != 0) return errfreq;
					if ((generation == nGeneration-1) || (next == -1)
							|| (popRead == popEnd))
					{
						FreqAdjnPrt (outLoc, moreDat, freqList, nloci,
										generation+1, locUse);
						for (g1=0; g1<generation; g1++) {
							for (g2=g1+1; g2<generation+1; g2++) {
								(*totPairTmp)++;
								nPairTmp++;
//								if (g2 == 1 && icount == 0)
//									printf ("\nTemporal Method ...\n");
								TemporalNeEst (outLoc, moreDat, freqList, nloci,
										locUse, g1, g2, nCrit, critVal,
										nTotAlle, nIndAlle,
										Hkmean, Hcmean, Hsmean,
										fkmean, fcmean, fsmean,
										fkprimeMean, fcprimeMean, fsprimeAll,
										NeTempk, NeTempc, NeTemps,
										loNek, hiNek, loNec, hiNec,
										loNes, hiNes, jloNek, jhiNek, jloNec,
										jhiNec, jloNes, jhiNes, param,
										jacknife, timeline, census,
										tempk, tempc, temps, infinite, weighsmp);
								PrtTemporal (output, nCrit, critVal,
										g1, g2, generation+1, timeline,
										nTotAlle, nIndAlle, Hkmean, Hcmean,
										Hsmean, fkmean, fcmean, fsmean,
										fkprimeMean, fcprimeMean, fsprimeAll,
										NeTempk, NeTempc, NeTemps, loNek,
										loNec, loNes, hiNek, hiNec, hiNes,
										jloNek, jloNec, jloNes, jhiNek, jhiNec,
										jhiNes, param, jacknife, infinite,
										tempk, tempc, temps, popIDtemp, popSize, census);
								if (tabX==0) PrtTempxFile (shOutputTemp, popRead, popStart,	//3
										newSet, lastSet, critVal, nCrit, g1, g2,//6
										generation+1, timeline, nIndAlle,	//3
										Hkmean, Hcmean, Hsmean,		//3
										fkmean, fcmean, fsmean,		//3
										fkprimeMean, fcprimeMean, fsprimeAll,//3
										NeTempk, NeTempc, NeTemps,	//3
										loNek, loNec, loNes,		//3
										hiNek, hiNec, hiNes,		//3, 30 so far
										jloNek, jloNec, jloNes,		//3
										jhiNek, jhiNec, jhiNes,		//3
										param, jacknife,			//2
										infinite, tempkx, tempcx, tempsx,	//3
										popIDtemp, topCrit, nPlan, census,	//3, 44 so far
										inpName, nLocUsed, icount, *totPairTmp, common);
								else PrtTpTabFile (shOutputTemp, popRead, popStart,	//3
										newSet, lastSet, critVal, nCrit, g1, g2,//6
										generation+1, nGeneration, timeline, nIndAlle,	//4
										Hkmean, Hcmean, Hsmean,		//3
										fkmean, fcmean, fsmean,		//3
										fkprimeMean, fcprimeMean, fsprimeAll,//3
										NeTempk, NeTempc, NeTemps,	//3
										loNek, loNec, loNes,		//3
										hiNek, hiNec, hiNes,		//3, 30 so far
										jloNek, jloNec, jloNes,		//3
										jhiNek, jhiNec, jhiNes,		//3
										param, jacknife,			//2
										infinite, tempkx, tempcx, tempsx,	//3
										popIDtemp, topCrit, nPlan, census,	//3, 44 so far
										inpName, nLocUsed, icount, *totPairTmp, common, &lenTpPop);
							};
						};
						if ((next==-1 || popRead==popEnd) && (common == 0))
							PrtTempPop (shOutputTemp, generation, nGeneration,
										nPoptemp, popRun+1, nPairTmp, timeline);
						// (note: popRun not incremented yet until all are done)
						RemoveFreq(freqList, nloci);
						generation = 0;
					} else generation++;
				};
// ------------------------------------------------------------------------
//				if (outBurr != NULL && moreBurr0 == 1)
// change in Nov 2014: message on Burrows file if only one file
				if (outBurr != NULL && moreBurr0 == 1 && sepBurOut == 0)
					printf ("Burrows coefficients for LD method are written to file %s.\n",
							outBurrName);

				RemoveAlle (alleList, nloci);
				if (makeFish > 0) RemoveFish (fishHead, nloci);
				popRun++;	// actual number of pops run
				(*totPop)++;
// comment this out if fishHead is not reallocated after each population
//				free (fishTail);
			};
//			if (next == -1) break;	// end of file, done.
// reestablish memory for all arrays
/*
			if ((alleList = (ALLEPTR*) malloc(sizeof(ALLEPTR)*nloci))
				== NULL) {
				printf  ("Out of memory to reserve alleles!\n");
				return -1;
			};
*/
/*
			if ((fishHead = (FISHPTR*) malloc(sizeof(FISHPTR)*nloci))
				== NULL  || (fishTail =
				(FISHPTR*) malloc(sizeof(FISHPTR)*nloci)) == NULL) {
				printf  ("Out of memory for sample list!\n");
				return -1;
			};
//*/
			// reinitialize: next = 1, so it's new pop
			for (isize=0; isize<size; isize++) *(missptr+isize)=0;
			for (p=0; p<nloci; p++){
				*(alleList+p) = NULL;
				*(okLoc+p) = 1;
				*(nMobil+p) = 0;
				*(minFreq+p) = 0;
				*(maxFreq+p) = 0;
				if (makeFish > 0) {
					*(fishHead+p) = NULL;
					*(fishTail+p) = NULL;
				};
			};
		};	// end of "if (next != 0)", that is, finishing calculations
		if (next == 1) {	// we are at the first sample of the next pop
			samp = 0;
			ind = 0;
			popRead++;	// this is then the pop number of the pop just done
			if (popRead > popEnd) break;	// just passed pop "popEnd"
			if (popRead >= popStart)
				printf ("\nPopulation%6d [%s]\n", popRead, newID);
		};
		if (next == -1) break;	// end of file, done.
		// now, next = 0 or 1. If = 1, we start the first sample.
		// nSampErr = the number of samples having data missing or error,
		// use ind for counting number of samples read, samp is used
		// for counting the number of samples used (since some samples
		// read later will be discarded if there is option limiting #sample)
		err = GetSample (input, nloci, sampData, lenM, &ind, lenBlock,
						&nSampErr, &noGen, genErr, &firstErr, locUse);
		samp = ind;
		if (samp > maxSamp) {	// this sample and subsequent ones
			samp = maxSamp;		// are read but not put in the list,
			continue;			// so maxSamp is the total samples used.
		};
		// also, no recording if popRead is not in the range yet
		if (popRead < popStart) continue;
// --------------------------------------------------------------------
// These are for checking if reading is ok. (Put data read in a file.)
// Add or delete a slash "/" at the line before the beginning of code
// to cover or open the code block (1 slash: cover; 2 slashes: open):
// A modification of these lines of code can produce an input file
// in different formats: FSTAT or GENEPOP or something else.
/*
FILE *check;	// the file is named "checkGeno.txt"
if (popRead == 1 && samp==1) check = fopen ("checkGeno.txt", "w");
if (samp==1) fprintf (check, "\nPopulation%6d [%s]\n", popRead, newID);
fprintf(check, "Sample %d\n", samp);
for (p=0; p<nloci; p++){
	fprintf(check, "%2.2d%2.2d ", *(sampData +2*p), *(sampData +2*p+1));
	if ((p+1)%20==0 || p+1==nloci) fprintf(check, "\n");
	fflush (check);
};
//*/ //----------------------------------------------------------------
		nErr += noGen;
		// open missing data file when the first error occurs.
		if (*missFileName !='\0' && nSampErr == 1 && missDat == NULL)
			missDat = PrtMisHead (missFileName, inpName, popRead, newID);
		// quit when sample contains a genotype too wide or has non-digit:
		if ((errCode = PrtError (output, missDat, nloci, nSampErr, popRead,
				samp, popID, err, noGen, genErr, firstErr)) > 0) break;

/*
		{
			if ((missDat=fopen (missFileName, "w")) != NULL)
			{
				fprintf (missDat, "Missing data from input file %s.\n"
				"(Only the first locus with missing data in a sample is listed.)\n\n",
					inpName);
				fprintf (missDat, "Population %d [%s]\n", popRead, newID);
				PrtLines (missDat, 45, '-');
				fprintf (missDat, "   Sample       Locus         Genotype     Total\n");
				fflush (missDat);
			};
		};
		if (err != 0) {
			if (err == -1 && missDat != NULL) {
				fprintf (missDat, "Population %d: Sample %d ends too soon.\n",
						popRead, samp);
				fprintf (output, "\nPopulation %d: Sample %d ends too soon.\n",
						popRead, samp);
				errCode = 3;
				break;
			};	// now, err > 0, it is of the form m*nloci + p, p=0,...,nloci-1.
			m = err/nloci; p = err % nloci;
//			errCode = PrtMisDat(missDat, m, popRead, p, samp, noGen, nSampErr);
			errCode = PrtMisDat(missDat, m, p, samp, noGen, genErr, firstErr);
			if (errCode != 0) {	// quit running with this input file.
			// since we quit, should output error messages to output file
				if (errCode == 1)
					fprintf (output, "\nFatal error: At locus %d, "
						"Sample %d (population %d) has too many characters"
						" for a genotype.\n", p+1, samp, popRead, popID);
				else	// errCode = 2
					fprintf (output, "\nFatal error: At locus %d, "
						"Sample %d (population %d) has non-digit character"
						" for a genotype.\n", p+1, samp, popRead, popID);
				fflush (output);
				break;
			};

		};
*/
		if ((AddAlleWide (alleList, nloci, sampData, nMobil, missptr,
			maxMobilVal, popRead, samp) != 0) ||
			AddFishWide(fishHead, fishTail, nloci, sampData, locUse, makeFish) == 0)
		{
			fprintf (output, "\n\nOut of memory at population %s, sample %d.\n",
					popID, samp);
			errCode = -1;
			fflush (output);
			break;
		};
	};
	if (popRun == 0) {
		fprintf (output, "No population is run!\n");
		printf ("No population is run!\n");
	};
// variables in all methods:
	if (makeFish > 0) {
		free (fishHead);	//1
		free (fishTail);	//2
	};
	free (alleList);	//3
	free (sampData);	//4
	free (missptr);		//5
	free (nMobil);		//6
	free (minFreq);		//7
	free (maxFreq);		//8
	free (okLoc);		//9
	free (popID);		//10
	free (newID);		//11
	free (wExpR2);		//12
// for LD method
	free (estNe);		//13
	free (wHarmonic);	//14
	free (rB2WAve);	//15
	free (r2Drift);	//15
	free (nIndSum);		//16
	free (jackOK);		//17

	free (confJacklow);	//18
	free (confJackhi);	//19
	free (confParalow);	//20
	free (confParahi);	//21
	free (Jdegree);
	if (mTemporal == 1) free (freqList);	//22
// for HetExcess method
	free (hSamp);		//23
	free (estHetN);		//23
	free (hetD);		//24
	free (loHetNe);		//25
	free (hiHetNe);		//26
//	free (jloHetNe);	//27
//	free (jhiHetNe);	//28
	free (indAlleH);	//29
// for temporal method
	free (nTotAlle);	//30
	free (nIndAlle);	//31
	free (NeTempk);		//32
	free (NeTempc);		//33
	free (NeTemps);		//34

	free (loNek);		//35
	free (hiNek);		//36
	free (loNec);		//37
	free (hiNec);		//38
	free (loNes);		//39
	free (hiNes);		//40
	free (jloNek);		//41
	free (jhiNek);		//42
	free (jloNec);		//43
	free (jhiNec);		//44
	free (jloNes);		//45
	free (jhiNes);		//46

	free (Hkmean);		//47
	free (Hcmean);		//48
	free (Hsmean);		//49
	free (fkmean);		//50
	free (fcmean);		//51
	free (fsmean);		//52
	free (fkprimeMean);		//53
	free (fcprimeMean);		//54
	free (fsprimeAll);		//55

	if (mTemporal == 1) {
		free (popSize);
		for (n=0; n<MAXGENERATION; n++)
			free (popIDtemp [n]);		//56
	};

//	if (nErr > 0) fclose (missDat);
	if (missDat != NULL) {
		fclose (missDat);
		printf ("\nMissing data are listed in file %s\n", missFileName);
	};
	return errCode;
}


//------------------------------------------------------------------

void PrintEndTime (FILE *output, time_t rawtime) {
	if (output == NULL) return;
	fprintf (output, "\nEnding time: %s", ctime (&rawtime));
	PrtLines (output, 37, '-');
	fprintf (output, "\n");
	fclose (output);
	output = NULL;
};

//------------------------------------------------------------------

void PrtPopRun (FILE *output, int totPop, int dashes) {
	if (output == NULL) return;
	fprintf (output, "\n");
	PrtLines (output, dashes, '-');
	fprintf(output, "Total number of populations =%8d\n", totPop);
	PrtLines (output, dashes, '-');
};


//------------------------------------------------------------------

void PrtRange (FILE *output, int numList[], int num)
// numList is supposed to be an array of "num" integers in ascending order.
// After adding 1 to each element, print to output as ranges of integers.
// For example, if numList consists of
// 0,1,2,3,7, 9,11, 14, 15, 16.
// Then print 1-4, 8, 10, 12, 15-17.
{
	int k, m, n;
	char buffer[20];
	for (n = 0; n < num; ) {
		for (m = n; m < num; m++) {
			if (numList[m] > numList[n]+1) break;
			else k = m;
		}
		// k is the largest >= n that has numList[k] == numList[n]
		if (k == n) sprintf (buffer, "%d", numList[n]+1);
		else sprintf (buffer, "%d-%d", numList[n]+1, numList[k]+1);
		fprintf (output, "%s", buffer);
		if (m < num) fprintf (output, ", ");
		else fprintf (output, "\n");
		n = m;
	}

}

//------------------------------------------------------------------
void PrtBriefChromo (FILE *output, struct chromosome *chromoList,
				int nChromo, int chroGrp, int unknown)
{
	int j, k, m, n;
	unsigned long long nBurrPair = 0;
	unsigned long nWithin = 0;
	int locSeen = 0;
	if (output == NULL) return;
	if (chromoList == NULL || nChromo <= 1) return;
	int chroRead = nChromo;
	if (unknown > 0) chroRead--;
	fprintf (output,
		"Chromosomes, followed by a colon and the number of loci in Genotype Input File:\n");
	k = 0;
	for (n = 0; n < chroRead; n++) {
		k++;
		m = chromoList[n].nloci;
		locSeen += m;
		nWithin += (m*(m-1))/2;
		j = strlen (chromoList[n].name) - 8;// to have last 8 chars printed
		if (j < 0) j = 0;				// as restricted in the format below
		fprintf (output, "%8s:%6d", (chromoList[n].name + j), m);
		if (n < (chroRead-1)) fprintf (output, ",");
		k = k%5;	// print 5 chromosomes per line
		if (k == 0) fprintf (output, "\n");
	}
	fprintf (output,
		"\nNumber of loci seen in genotype and [chromosomes/loci] files: %d\n",
		locSeen);
	if (unknown > 0) {
		fprintf (output, "Loci not in [chromosomes/loci] file are assigned"
			" default chromosome \"%s\"\n", chromoList[chroRead].name);
		m = chromoList[chroRead].nloci;
		locSeen += m;
	}
	nBurrPair = (unsigned long long) locSeen;
	nBurrPair *= (locSeen - 1);
	nBurrPair /= 2;
	if (chroGrp == 1) {
		fprintf (output,
			"Each pair of loci are taken within a single chromosome.\n");
		printf ("Each pair of loci are taken within a single chromosome.\n");
		nBurrPair = (unsigned long long) nWithin;
	} else {
		fprintf (output,
			"Each pair of loci are taken in distinct chromosomes.\n");
		printf ("Each pair of loci are taken in distinct chromosomes.\n");
		nBurrPair -= nWithin;
	}
	fprintf (output,
			"Maximum number of locus pairs = %llu\n\n", nBurrPair);
	printf ("Maximum number of locus pairs = %llu\n\n", nBurrPair);
}


//------------------------------------------------------------------
void PrtChromo (FILE *output, struct chromosome *chromoList,
				int nChromo, int chroGrp, int unknown)
{
	int m, n;
	unsigned long long nBurrPair = 0;
	unsigned long nWithin = 0;
	int locSeen = 0;
	if (output == NULL) return;
	if (chromoList == NULL || nChromo <= 1) return;
	int chroRead = nChromo;
	if (unknown > 0) chroRead--;
	fprintf (output,
		"Chromosomes and their (numbered) loci in genotype input file:\n");

	for (n = 0; n < chroRead; n++) {
		m = chromoList[n].nloci;
		locSeen += m;
		nWithin += (m*(m-1))/2;
		fprintf (output, "* %s (%d loci):  ", chromoList[n].name, m);
		PrtRange (output, chromoList[n].locus, m);
	}
	fprintf (output,
		"\nNumber of loci seen in [chromosomes/loci] file: %d\n", locSeen);
	if (unknown > 0) {
		fprintf (output, "Loci not in [chromosomes/loci] file are assigned"
			" default chromosome:\n");
		m = chromoList[chroRead].nloci;
		locSeen += m;
		nWithin += (m*(m-1))/2;
		fprintf (output, "* %s (%d loci):  ", chromoList[chroRead].name, m);
		PrtRange (output, chromoList[chroRead].locus, m);
	}
	nBurrPair = (unsigned long long) locSeen;
	nBurrPair *= (locSeen - 1);
	nBurrPair /= 2;
	if (chroGrp == 1) {
		fprintf (output,
			"Each pair of loci are taken within a single chromosome.\n");
		nBurrPair = (unsigned long long) nWithin;
	} else {
		fprintf (output,
			"Each pair of loci are taken in distinct chromosomes.\n");
		nBurrPair -= nWithin;
	}
	fprintf (output,
			"Maximum number of locus pairs = %llu\n\n", nBurrPair);
}


//------------------------------------------------------------------
// add last parameter "common" to have option not to close output files
// when they are needed to write outputs for next input file.
// This parameter should be set = 0 except when running multiple
// input files by calling RunMultiCommon.

int RunPop (int icount, char *inpName, FILE *input, char append,  FILE *output,
			char *outFolder, struct locusMap *locList,	// added in Nov 2014
			FILE *outLoc, char *outLocName, FILE *outBurr,
			char *outBurrName, FILE *shOutputLD, FILE *shOutputHet,
			FILE *shOutputCoan, FILE *shOutputTemp,				// 14 so far
			int popLoc1, int popLoc2, int popBurr1, int popBurr2, int topBCrit,
			int popStart, int popEnd,						// 21 so far
			int maxSamp, int nloci, int lenM, int maxMobilVal, int nCrit,
			float critVal[], char format, char param,		// 29 so far
			char jacknife, char *locUse, char mating, char *missFileName,
			float infinite, int lenBlock,
			char mLD, char mHet, char mNomura, char mTemporal,	// 39 so far
			int nGeneration, float timeline[], AGEPTR *ageSeq, int nSeq,
			char getAge, int tempClue, int tempxClue, char byRange,	//47
			int topCrit, int nPlan, int census,
			int *totPop, int *totPairTmp, char common, char tabX,
//			int nGeneration, float timeline[], FILE *info)
			// add parameters in Apr 2015:
			char sepBurOut, char moreCol, char BurAlePair,
			struct chromosome *chromoList, int nChromo, int chroGrp, int unknown)
{
	int err;
	time_t rawtime;
	char *prefix;
	int i, j, k, m;
	if (output == NULL || input == NULL) return 0;
	PrtHeader (output, append, inpName, icount, 1);
// Add Apr 2015:
	if (mLD == 1) PrtBriefChromo (output, chromoList, nChromo, chroGrp, unknown);
	// The code up to "if (common == 0)" is to remove path names in inpName
	prefix = (char*) malloc (sizeof(char)*PATHFILE);
	k = GetPrefix (inpName, prefix, PATHFILE, PATHCHR);
	m = strlen(inpName) - 1;
	// get the first dot from the right of inpName to determine extension:
	for (i=m; i >= 0; i--) if (*(inpName + i) == '.') break;
	// so if there is a dot in the name, then i >= 0
	if (i >= 0) {
		m -= i;	// m now is the length of extension
		// concatenate the extension including the dot, to the end of prefix
		for (j = 0; j <= m; j++) *(prefix + k + j) = *(inpName + i + j);
		k += (m+1);
	};
	// k now is the length of inpName without path names, reassign inpName:
	for (j = 0; j < k; j++) *(inpName + j) = *(prefix + j);
	*(inpName + k) = '\0';
	free (prefix);
	// the reason of doing the above is to get rid of path characters
	// in the name inpName, which was done by GetPrefix;
	// The reduced name will be printed in column of xtra output if common = 1
	if (common == 0) {
	// Only print limited use of loci, samples, populations when NOT
	// running multiple files with common tabular-format outputs.
		if (shOutputLD != NULL) {
			PrtHeader (shOutputLD, append, inpName, icount, 0);
			PrtLimitUse (shOutputLD, locUse, nloci, byRange, popStart,
					popEnd, MAX_POP, maxSamp, "Population");
		};
		if (shOutputHet != NULL) {
			PrtHeader (shOutputHet, append, inpName, icount, 0);
			PrtLimitUse (shOutputHet, locUse, nloci, byRange, popStart,
					popEnd, MAX_POP, maxSamp, "Population");
		};
		if (shOutputCoan != NULL) {
			PrtHeader (shOutputCoan, append, inpName, icount, 0);
			PrtLimitUse (shOutputCoan, locUse, nloci, byRange, popStart,
					popEnd, MAX_POP, maxSamp, "Population");
		};
		if (shOutputTemp != NULL) {
			PrtHeader (shOutputTemp, append, inpName, icount, 0);
			PrtLimitUse (shOutputTemp, locUse, nloci, byRange, popStart,
					popEnd, MAX_POP, maxSamp, "Sample");
		};
	};

	err = RunPop0 (icount, inpName, input, append, output, outFolder, locList,
				outLoc, outLocName, outBurr, outBurrName,
				shOutputLD, shOutputHet, shOutputCoan, shOutputTemp,
				popLoc1, popLoc2, popBurr1, popBurr2, topBCrit,
				popStart, popEnd, maxSamp, nloci, lenM, maxMobilVal,
				nCrit, critVal, format, param, jacknife, locUse, mating,
				missFileName, infinite, lenBlock, mLD, mHet, mNomura,
				mTemporal, nGeneration, timeline,
				ageSeq, nSeq, getAge, tempClue, tempxClue, byRange,
				topCrit, nPlan, census, totPop, totPairTmp, common, tabX,
//				mTemporal, nGeneration, timeline, info);
	// add parameters in Apr 2015:
				sepBurOut, moreCol, BurAlePair,
				chromoList, nChromo, chroGrp);
	fclose (input);
	if (outLoc != NULL) fclose (outLoc);
	if (outBurr != NULL) fclose (outBurr);
	time ( &rawtime );
	// when there are more to write to output files, then return, not close:
	if (common != 0) return err;
	PrintEndTime (output, rawtime);
	PrintEndTime (shOutputLD, rawtime);
	PrintEndTime (shOutputHet, rawtime);
	PrintEndTime (shOutputCoan, rawtime);
	PrintEndTime (shOutputTemp, rawtime);
	return err;

}


//------------------------------------------------------------------


void SetDefault (int *maxSamp, char *param, char *nonparam, char *mLD,
				char *mHet, char *mNomura, char *mTemporal, int *nCrit,
				float critVal[], int *nGeneration, char *mating,
				float timeline[])
{
	int n;
	*param = 1;
	*nonparam = 1;
	*nCrit = NCUT_SET;
	*mLD = LDACTION;		// run LD method
	*mHet = HETACTION;		// run Het-Excess method
	*mNomura = COANACTION;	// run Coancestry (Nomura) method
	*mating = MATING;
	// for temporal method
//	*mTemporal = TEMPACTION;	// default: run temporal methods
	*mTemporal = 0;			// No temporal
	*maxSamp = MAX_SAMP;	// maximum number of samples per pop
						// that will be processed. This maxSamp
						// may be reassigned from the user input.
	critVal[0] = (float) 0.05;
	critVal[1] = (float) 0.02;
	critVal[2] = (float) 0.01;
	for (n=3; n<MAXCRIT; critVal[n++] = 0);


	// for temporal method
	*nGeneration = 1;	// when = 1: no temporal method, NOT to be changed
	for (n=0; n<MAXGENERATION; timeline [n++] = 0);
}

//------------------------------------------------------------------
int RunMultiFiles (char *mFileName, char mOpt)
// Read file mFileName to run RunPop multiple times, Each run requires
// 3 lines:
//	* first line: input file name
//	* second line: output file name, if left empty, get default output name.
//	  (if default, output will be in the same directory as input)
//	* third line: an "Y" for next run, else to end.
// Return the number of successful runs.
// If input file name is empty, quit.
{
	#define LOCRANGE 100	// max endpoints of locus ranges to be used
	#define MAXLOCI 1000000		// maximum number of loci to be used
	int nRanges;	// Effective number of ranges
	int locRanges[LOCRANGE];	// Endpoints for ranges of loci

	int maxSamp, popStart, nCrit, nGeneration;
	char param, nonparam, mLD, mHet, mNomura, mTemporal, mating;
	float critVal[MAXCRIT];
	float timeline[MAXGENERATION+1];	// add 1 for census size

	char xOutLD = 0;
	char xOutHet = 0;
	char xOutCoan = 0;
	char xOutTemp = 0;
	int popEnd = 0;

	char format;
	char *locUse;
//	char *locDrop;
	char byRange = 0;

	FILE *mInpFile = NULL;
	char append;
	int i, n, count, c, nloci, nPop, maxMobilVal, lenM, k;
	int line;
	int nlocDel = 0;
	int *xClues;
	int tempClue, tempxClue;
	int topCrit = MAXCRIT;
	int census = 0;
	int nPlan = 0;
	int totPop = 0, totPairTmp = 0;
	// add July 2013:
	char tabX;
//	char prefix[LENFILE] = "\0";
//	char prefix[PATHFILE] = "\0";
	FILE *input = NULL;
	FILE *output = NULL;
	FILE *shOutputLD = NULL;
	FILE *shOutputHet = NULL;
	FILE *shOutputCoan = NULL;
	FILE *shOutputTemp = NULL;
	char *prefix = (char*) malloc (sizeof(char)*PATHFILE);
	char *inpName = (char*) malloc (sizeof(char)*PATHFILE);
	char *outName = (char*) malloc (sizeof(char)*PATHFILE);
	char *outNameMore = (char*) malloc (sizeof(char)*PATHFILE);
	mating = 0;
//	locDrop = (char*) malloc(sizeof(char)*100);
	if ((mInpFile = fopen (mFileName, "r")) == NULL) {
		perror (mFileName);
		return 0;
	};
//	SetDefault (&maxSamp, &param, &nonparam, &mLD,
//				&mHet, &mNomura, &mTemporal, &nCrit,
//				critVal, &nGeneration, &mating, timeline);
	*outNameMore = '\0';	// to hold previous output name
	count = 0;
	line = 0;
	nRanges = 1;
	locRanges[0] = 1;
	locRanges[1] = MAXLOCI;
	for ( ; ; ) {
		xOutLD = 0;
		xOutHet = 0;
		xOutCoan = 0;
		xOutTemp = 0;
		tempClue = 0;
		tempxClue = 0;
		nPlan = 0;
		topCrit = MAXCRIT;
		append = 0;	// 0 for overwrite output, 1 for append, assumed overwrite
		shOutputLD = NULL;	// assumed no short output file
		shOutputHet = NULL;	// assumed no short output file
		shOutputCoan = NULL;	// assumed no short output file
		shOutputTemp = NULL;	// assumed no short output file
		nlocDel = 0;	// assumed no loci deleted
		nGeneration = 0;
		SetDefault (&maxSamp, &param, &nonparam, &mLD,
				&mHet, &mNomura, &mTemporal, &nCrit,
				critVal, &nGeneration, &mating, timeline);

		popStart = 1;
		popEnd = MAX_POP;
// insert in Dec 2011, to read methods and critical values first,
// then generation timeline if temporal is desired
		if (FindMethod(mInpFile, mFileName, &line, &mLD, &mHet, &mNomura,
			&mTemporal, &tempClue) < 0) break;
		line++;
		if ((n=CritValRead (mInpFile, MAXCRIT, critVal, &i)) <= 0) {
			ErrMsg (mFileName, "ERROR on Number of Critical Value", line);
			break;
		} else nCrit = n;
		if (i > 0) line++;
		if (mTemporal == 1) {
			line++;
//			if ((n = GeneratnRead (mInpFile, &nGeneration,
//				timeline, MAXGENERATION+1, 1, &census)) <= 1) {
			if ((n = GeneratnRead (mInpFile, &nGeneration,
				timeline, MAXGENERATION+1, &census)) <= 1) {
				ErrMsg (mFileName, "No valid generation timeline for temporal!",
						line);
				mTemporal = 0;
			} else {	// good generation set
			// the first index in timeline is for census size, which is
			// assigned to variable census; so need to move back one index:
			// (actually upto index nGeneration is enough)
				for (k=0; k< nGeneration; k++) timeline[k] = timeline[k+1];
				if (census > 0) nPlan = 2;
				else nPlan = 1;
			};
		};
		if (mOpt == 1) {	// the file contains more options
			line++;	// extra output and its option: ------------------
			// add 1 for tab-delimiter option, July 2013
			xClues = (int*) malloc(sizeof(int)*4);
			// default values:
			xClues[0] = 0;		// sum of methods to have extra outputs
			xClues[1] = 0;		// all temporal methods included
			xClues[2] = MAXCRIT;	// all Pcrits are in extra outputs
			xClues[3] = TABX;		// default values for tab between columns
			c = 1;	// set c = 1 so that GetClues puts cursor to next line
			k = GetClues (mInpFile, xClues, 4, c);	// number of clues read
			if (k <= 0) {
				ErrMsg (mFileName, "At reading clues for tabular-format output!",
						line);
				break;
			} else {
				SetMethod(xClues[0], &xOutLD, &xOutHet, &xOutCoan, &xOutTemp);
				if (mLD == 0) xOutLD = 0;
				if (mHet == 0) xOutHet = 0;
				if (mNomura == 0) xOutCoan = 0;
				if (mTemporal == 0) xOutTemp = 0;
			};
			tempxClue = xClues[1];
			topCrit = xClues[2];
			tabX = (xClues[3] == 0)? 0: 1;
			free(xClues);
/*** Bring these after mating model
			line++;	// max samples per pop:---------------------------
			GetInt (mInpFile, &maxSamp, 1);
			if (maxSamp <= 0) maxSamp = MAX_SAMP;

			line++;	// range of populations to run:
			k = 0;
			i = GetPair(mInpFile, &k, &n, 1);	// i = 0: no number or
							// negative for the first entry, then k still = 0
			if (k > 0) {	// No warning error if i=0. If k > 0, then i >=1
				popEnd = k;
				if (i == 2) {	// then n is obtained a value from GetPair
					if (n >= k) {
						popStart = k;	// k <= n, so range is from k to n.
						popEnd = n;
					};	// else, n < k: the 2nd entry is erroneous, popEnd is
				// not reassigned, the range is from popStart=1 to popEnd=k
				};	// else, only 1 number: pop range is from popStart=1 to k
			};	// else: first entry is <=0, range taken from 1 to MAX_POP
*/
			line++;			// run parametric CIs or not -------------
			n = 1;	// default: have parameter CIs
			GetInt (mInpFile, &n, 1);
			param = (n != 0)? 1: 0;
			nonparam = param;
/* Skip reading for jackknife by assigning nonparam = param on line above
			line++;			// run nonparam CIs or not ---------------
			n = 1;
			GetInt (mInpFile, &n, 1);
			nonparam = (n != 0)? 1: 0;
*/
// Change in Jan 2014: only read mating model if LD is included:
			if (mLD == 1) {
				line++;		// read mating model: 0 for random, 1 for monogamy ---
				n = 0;		// default: random mating
				GetInt (mInpFile, &n, 1);
				mating = (n != 0)? 1: 0;
			};

			line++;	// max samples per pop:---------------------------
			GetInt (mInpFile, &maxSamp, 1);
			if (maxSamp <= 0) maxSamp = MAX_SAMP;

			line++;	// range of populations to run:
			k = 0;
			i = GetPair(mInpFile, &k, &n, 1);	// i = 0: no number or
							// negative for the first entry, then k still = 0
			if (k > 0) {	// No warning error if i=0. If k > 0, then i >=1
				popEnd = k;
				if (i == 2) {	// then n is obtained a value from GetPair
					if (n >= k) {
						popStart = k;	// k <= n, so range is from k to n.
						popEnd = n;
					};	// else, n < k: the 2nd entry is erroneous, popEnd is
				// not reassigned, the range is from popStart=1 to popEnd=k
				};	// else, only 1 number: pop range is from popStart=1 to k
			};	// else: first entry is <=0, range taken from 1 to MAX_POP

// if want to read loci dropped, can do with function LociDropped here:
//			line++;		// read loci dropped, limited to 100 (size of locDrop)
			// (line is incremented in the call LociDropped)
			// In LociDropped, the 5th parameter = 1 to indicate that all
			// loci dropped must be on 1 line, the last 1 is to indicate the
			// array represents the dropped loci: if locDrop [i] = p, locus p
			// is dropped, then locUse[p-1] should be assigned 0.
//			for (i=0; i<100; i++) *(locDrop+ i) = 0;
//			nlocDel = LociDropped (mInpFile, locDrop, 100, &line, 1, 1, &byRange);

			nRanges = GetRanges (mInpFile, locRanges, LOCRANGE, MAXLOCI, &byRange);
		};	// end of reading extra parameters when mOpt = 1
	// get new input file name. If empty, quit
		line++;
		if ((n=GetToken (mInpFile,inpName,PATHFILE,BLANKS,ENDCHRS,&c, &k)) <= 0)
		{
			ErrMsg (mFileName, "At reading input name", line);
			printf("%s\n", inpName);
			break;
		}
		for (; (c=fgetc(mInpFile)) != EOF && c!='\n';);
	// get prefix of input name
	// (minus 7 to reserve for "Ne.out" and "NeTemps.out")
//		GetPrefix (inpName, prefix, PATHFILE-7, '\0');
		GetPrefix (inpName, prefix, PATHFILE-7, "\0");
		format = FSTAT;
		// n is the length of inpName without trailing blanks
		if (n > 4) {
			if (*(inpName+n-4) == '.' && tolower(*(inpName+n-3)) == 'g' &&
				tolower(*(inpName+n-2))=='e' && tolower(*(inpName+n-1))=='n')
			format = GENPOP;
		};	// format will be doubled check later
		line++;
	// get output name, if empty, make output name from the input name.
	// k = # of trailing blanks (in BLANKS)  before reaching char in ENDCHRS
// added Feb 15 2013:
		append = 0;
		if (GetToken (mInpFile, outName, PATHFILE, BLANKS, ENDCHRS, &c, &k) <= 0)
		{
			outName = strcat(strcat(outName, prefix), "Ne");
			outName = strcat(outName, EXTENSION);
			if (c == SPECHR)	// need to go to next line
				for (; (c=fgetc(mInpFile)) != EOF && c!='\n';);
		}
		else {	// need to go to next line
// added Feb 15 2013:
			if (c == SPECHR && k == 0) append = 1;
			for (; (c=fgetc(mInpFile)) != EOF && c!='\n';);
		};
		if ((input = fopen (inpName, "r")) == NULL) {
			printf ("\nERROR:\n");
			perror(inpName);
			if (tolower(fgetc (mInpFile)) != 'y') break;
			for (; (c=fgetc(mInpFile)) != EOF && c!='\n';);
			if (c == EOF) break;
			continue;
		}
	// append if same output name: CASE SENSITIVE (since Unix cares!)
// changed Feb 15 2013:
		if (append == 0) append = (strcmp(outName, outNameMore) == 0)? 1: 0;
		if (append == 1) {
			output = fopen (outName, "a");	// append if the same output
		} else {
			output = fopen (outName, "w");
		}
		*outNameMore = '\0';
		strcat (outNameMore, outName);
		if (output == NULL) {
			printf ("\nCannot open file %s for output\n", outName);
			continue;
		}

		printf ("\n>>> Input %d: [%s], ", count+1, inpName);

// now read input file to obtain necessary parameter values
		// if input file extension is not "gen", it was assumed FSTAT format
		if (format == FSTAT) {
			if (GetInfoDat(input,
						&nPop, &nloci, &maxMobilVal, &lenM, LEN_BLOCK)==0) {
				// this is not FSTAT, now assume it's GENPOP format
				format = GENPOP;
				rewind(input);
			}
		}
		if (format == GENPOP) printf ("GENEPOP format\n");
		if (format == FSTAT)  printf ("FSTAT format\n");
		if (format == GENPOP) {	// GENEPOP format
			if ((nloci = GetnLoci (input, LEN_BLOCK, &lenM)) <= 0) {
				fclose (input);
				fclose (output);
				printf ("Error in input file [%s]\n", inpName);
				if (tolower(fgetc (mInpFile)) != 'y') break;
				for (; (c=fgetc(mInpFile)) != EOF && c!='\n';);
				if (c == EOF) break;
				continue;
			} else {
				rewind (input);
				for (; (c=fgetc(input)) !='\n';);	//position at first locus
			}
	// should assign maxMobilVal to make sure it doesn't get garbage.
	// Since the length is lenM, assign max possible value:
			for (maxMobilVal=1, i=1; i<=lenM; i++) maxMobilVal *= 10;
		}
		n = mLD + mHet + mNomura + mTemporal;
		PrtMethod (n, mLD, mHet, mNomura, mTemporal);
		printf ("Number of loci = %d, %d-digit alleles\n", nloci, lenM);
		locUse = (char*) malloc (sizeof(char)*nloci);


		for (i=0; i<nloci; i++) *(locUse+i) = 0;
		for (i=0; i<nloci; i++) {
			for (k=0; k<nRanges; k++) {
				if ((locRanges[2*k]<=(i+1)) && ((i+1)<= locRanges[2*k+1]))
					*(locUse+i) = 1;
			};
		};
		for (i=0; i<nloci; i++) if (*(locUse+i) != 1) nlocDel++;

		// default: all loci will be used.
//		for (i=0; i<nloci; i++) *(locUse+i) = 1;

//		for (i=0; i<nlocDel; i++) {
//			n = locDrop[i] - 1;
//			if (n >= 0 && n < nloci) *(locUse+ n) = 0;
//		};

		if (GetLocUsed (input, nloci, locUse, nloci-nlocDel, NULL) != 0) {
			fclose (input);
			fclose (output);
			if (tolower(fgetc (mInpFile)) != 'y') break;
			for (; (c=fgetc(mInpFile)) != EOF && c!='\n';);
			free (locUse);
			if (c == EOF) break;
			continue;
		};
		printf ("Output: [%s]", outName);
		if (append == 1) printf (" (Append)"); printf ("\n");
	// outName is no longer used, so assign as shorter output name
	// if so required, also prefix is done after used for input
		*prefix = '\0';
		strcat (prefix, outName);
		if (xOutLD == 1) {
			GetXoutName (outName, prefix, PATHFILE, XFILSUFLD, PATHCHR);
			if (append == 1)
				shOutputLD = fopen (outName, "a");
			else
				shOutputLD = fopen (outName, "w");
			printf ("Tabular-format LD Output: [%s]", outName);
			if (append == 1) printf (" (Append)"); printf ("\n");
		};
		if (xOutHet == 1) {
			GetXoutName (outName, prefix, PATHFILE, XFILSUFHET, PATHCHR);
			if (append == 1)
				shOutputHet = fopen (outName, "a");
			else
				shOutputHet = fopen (outName, "w");
			printf ("Tabular-format Het. Excess Output: [%s]", outName);
			if (append == 1) printf (" (Append)"); printf ("\n");
		};
		if (xOutCoan == 1) {
			GetXoutName (outName, prefix, PATHFILE, XFILSUFCOAN, PATHCHR);
			if (append == 1)
				shOutputCoan = fopen (outName, "a");
			else
				shOutputCoan = fopen (outName, "w");
			printf ("Tabular-format Coancestry Output: [%s]", outName);
			if (append == 1) printf (" (Append)"); printf ("\n");
		};
		if (xOutTemp == 1) {
			GetXoutName (outName, prefix, PATHFILE, XFILSUFTEMP, PATHCHR);
			if (append == 1)
				shOutputTemp = fopen (outName, "a");
			else
				shOutputTemp = fopen (outName, "w");
			printf ("Tabular-format Temporal Output: [%s]", outName);
			if (append == 1) printf (" (Append)"); printf ("\n");
		};
		if (RunPop (++count, inpName, input, append, output,
				NULL,	// this is for Burrows file, so set NULL
				NULL,	// locList
				NULL, "\0", //NULL,		// <-NULL for outLoc file and outLocName
				NULL, "\0", //NULL, 	// <-NULL for outBurr file and outBurrName
				shOutputLD, shOutputHet, shOutputCoan, shOutputTemp,
				0, 0, 0, 0, 0, popStart, popEnd, maxSamp, nloci, lenM,
//				0, 0, 0, 0, 0, 1, popEnd, maxSamp, nloci, lenM,
				// the 0s are for Loc, Bur, the 1 is for population starting number
				maxMobilVal, nCrit, critVal, format, param, nonparam,
				locUse, mating,
				"\0",	// empty for missFileName: no missing data file created
				INFINITE, LEN_BLOCK, mLD, mHet, mNomura, mTemporal,
				nGeneration, timeline, NULL, 0,	// NULL, 0: no list to store
				0, tempClue, tempxClue, byRange, topCrit, nPlan,
				// next 0 is for no getting age from a linked list for timeline
				census, &totPop, &totPairTmp, 0, tabX,
				// for generations, then 0 for "NOT" common
				0, 0, 0, NULL, 0, 0, 0) == 0) // add parameters Apr 2015
			printf("Finish running input %d.\n", count);
		free (locUse);
// these are already closed in RunPop
//		fclose (input);
//		fclose (output);
		if (tolower(fgetc (mInpFile)) != 'y') break;
		line++;
		for (; (c=fgetc(mInpFile)) != EOF && c!='\n';);
		if (c == EOF) break;
	};

	fclose (mInpFile);
	return count;

}


//------------------------------------------------------------------


int RunDirect (char misFilSuf[])
{
	int maxSamp, nCrit, nGeneration;
	char param, nonparam, mLD, mHet, mNomura, mTemporal, mating;
	float critVal[MAXCRIT];
	float timeline[MAXGENERATION];

	char format;
	int nPop, maxMobilVal, lenM;	// lenM = number of digits in alleles
	int popEnd = 0;
//	int popLoc=0, popBurr=0;	// # pops that outLoc, outBurr produce.
						// if nonzero, will produce corresponding file
// nloci: number of loci:
	int n, p, nloci=0;
	int topCrit = MAXCRIT;
	int nRun;
	char done;
	char *locUse;
	char *missFileName, *inpName, *prefix, *outName;
	FILE *input;
	FILE *output;
	int totPop = 0, totPairTmp = 0;

	inpName  = (char *) malloc(LENFILE * sizeof(char));
	*inpName  = '\0';
	prefix  = (char *) malloc(LENFILE * sizeof(char));
	*prefix = '\0';
	outName = (char *) malloc(PATHFILE * sizeof(char));
	*outName = '\0';
	missFileName = (char *) malloc(LENFILE * sizeof(char));
	*missFileName  = '\0';


	SetDefault (&maxSamp, &param, &nonparam, &mLD,
				&mHet, &mNomura, &mTemporal, &nCrit,
				critVal, &nGeneration, &mating, timeline);
	// reserve memory common for all populations

	done = 0;
	nRun = 0;
	while (done == 0) {
		*prefix = '\0';
		*outName = '\0';
		*missFileName  = '\0';
	// argc=1: only the name of the executable is on the command line.
	// The program will run with only one input file.
		if ((input = Prompt (inpName, prefix, LENFILE-9, &nPop, &nloci,
							&maxMobilVal, &lenM, &format,
							&mLD, &mHet, &mNomura, &mTemporal,
							&nGeneration, timeline)) == NULL) {
			perror (inpName);
			return nRun;
		};
//		printf ("Number of loci: %d, %d-digit alleles\n", nloci, lenM);
	// default: all loci will be used.
		locUse = (char*) malloc(sizeof(char)*nloci);
		for (p=0; p<nloci; p++) *(locUse+p) = 1;
		// since nPop is not known in GENEPOP format, assign max value:
		if (format == GENPOP) nPop = MAX_POP;
// popEnd is the number of populations in input file that the
// user wants to run. This should be the default:
		popEnd = nPop;
		// read input, pass locus names and print on screen:
		if (GetLocUsed (input, nloci, locUse, nloci, NULL) != 0)
			exit (EXIT_FAILURE);
		// prefix comes from the call to Prompt above
		// Assign names for auxiliary files based on prefix of input name
		strcat (missFileName, prefix);
		strcat (missFileName, misFilSuf);
//printf ("Missing file name = %s\n", missFileName);

//		output = GetOutFile (outName, prefix);
		output = GetOutFile (outName, prefix, mLD, mHet, mNomura, mTemporal);
		// After GetOutFile, prefix is the prefix of output name.
		// If an option for having tabular-format output files are wanted,
		// then this will be used for naming them.
		if (output == NULL) {
			printf ("\nCannot open file for output. Program aborted!\n");
			exit (EXIT_FAILURE);
		};
		if (RunPop (0, inpName, input, 0, output,
					NULL,	// for Burrows file, set NULL
					NULL,	// locList
					NULL, "\0",		// NULL for outLoc, outLocName,
					NULL, "\0",		// NULL for outBurr, outBurrName,
					NULL, NULL, NULL, NULL,	// NULL for short output files,
					0, 0, 0, 0,	0,		// for popLoc, popBurr,
					1, popEnd, maxSamp,	// the 1 for population starting number
					nloci, lenM, maxMobilVal,
					nCrit, critVal, format, param, nonparam, locUse,
					mating, missFileName, INFINITE, LEN_BLOCK,
					mLD, mHet, mNomura, mTemporal, nGeneration, timeline,
					NULL, 0,	// NULL, 0 for no generation list,
					0, 0, 0, 0, topCrit, 1, 0, &totPop, &totPairTmp, 0, 0,
					// next 0 for no attempt to retrieve
					// next 0 is for default value of tempClue: all or nothing
					// next 0 is for keeping same temporal in xFile as in main
					// next 0 is in place for using range of loci
					// last 1 if for only one plan (which is plan II)
					// next 0 for plan II
					// next-to-last 0 for "Not" common
					// last 0 for no tab in tabular-format output (redundant)
					// Apr 2015: parameters added
					0, 0, 0, NULL, 0, 0, 0)!= 0) return nRun;
		nRun++;
// temporarily exit (i.e., only run one input file, then exit the program):
//		break;

		// ask the user if want to continue with another input file
		printf ("\n> Run another input file? ");
// if do as the next line of code, then lower or upper case will be forgotten
//		if (n=(tolower(getchar())) != 'y') done = 1;
// so we put variable p to hold it.
		p = getchar();
		if ((n=tolower(p)) != 'y') done = 1;
		else {	// will prompt for another data file to run,
			// hence, need to clear up all remaining characters on this
			// input line, up to and including the end of line character
			// so that they are not part of the next input on screen
			for (; getchar() !='\n';);
			// this is for putting what user entered. It is not necessary,
			// only informative when the reading 'y' is from a file served
			// as a batch file:
			putchar (p); printf (": continue with input #%d\n", nRun+1);
		};
		free (locUse);

	};	// end of while (done == 0)
	free (inpName);
	free (prefix);
	free (outName);
	free (missFileName);
	return nRun;
}

// -------------------------------------------------------------------------
// April 2015: add functions RmChromo, GetChromo, and ChromoInp
//--------------------------------------------------------------------------

void RmChromo(struct chromosome *chromoList)
{

	int p;
	int n = sizeof(chromoList);
	for (p=0; p<n; p++) {
		free (chromoList[p].locus);
	}
	free(chromoList);

}

// --------------------------------------------------------------------------

struct chromosome* GetChromo (FILE *chroInp, int nlocUsed,
						struct locusMap *locList, int *nChromo, int *unknown)
// Read file input, which contains chromosome names and their loci.
// Assuming on each line, the first two strings are the names of a chromosome
// and of a locus contained in that chromosome.
// From the list locList containing locus names and locus numberings, compare
// the names read from input file and the names on the list, to determine
// chromosomes, and their loci contained in.
// Put these findings in an array of chromosomes (return value). Each element
// in the array contains the name of the chromosome, the number of loci
// contained in it, and the list of those loci (in term of their numberings).
// Parameter unknown will be the number of loci (in genotype input) that are
// not assigned a chromosome. Those are grouped in an imaginary chromosome.
{
	// This data type is used only here, as a step to create a list containing
	// all necessary info
	typedef struct chro *CHROPTR;
	struct chro
	{
		char name [LEN_LOCUS];	// name of the chromosome
		int nloci;		// the number of loci in this chromosome
		CHROPTR next;
	};
	CHROPTR curr, prev, newptr;
	CHROPTR chroTemp = NULL;
	int num = 0;
	int p, m, c, n;
//	int val;
	struct chroName
	{
		char name [LEN_LOCUS];	// name of the chromosome
	};
	struct chromosome *chromoList;
	struct chroName *chroAtLoc;
	chroAtLoc = (struct chroName*) malloc(sizeof(struct chroName)*nlocUsed);
	for (n=0; n<nlocUsed; n++) chroAtLoc[n].name[0] = '\0';
//  * chroAtLoc[n] will be the name of chromosome containing the nth-locus
//    in the array of locus locList. Recall that locList is an array,
//    where locList[n].name is the name of the nth-locus, locusList[n].num
//    is the numbering of that locus in the genotype input file.
//    The array locList was formed from reading genotype input file, so the
//    chromosome names, stored in locList[n].chromo, were not known yet.
//    When reading lines in "input", grab locus name in "input" then search
//    in locList for that name. When found, then can assign chromosome name
//    for the locus in locList[n].
//  * "chroTemp" will be the linked list of distinct chromosomes from "input".
//    Each element in the list contains a chromosome name and the number of
//    loci in it. The number of elements in the list is known after reading
//	  input (assigned *nChromo). This list will be used to create a more
//	  detailed list, chromoList, to be the return value of the function.

//	  This chromoList will be declared an array (dimension *nChromo), where
//	  each element stores loci by their numberings in the genotype input file.
//	  The lists chroAtLoc, chroTemp, and locList will speed up the process of
//	  determining these loci in each chromosome. Each chromosome and its loci
//	  are stored as fields "name" and "locus" of an element in chromoList.
//    The loci stored in each element of chromoList will be in ascending
//    order according to their orderings in genotype input file.
//	  As an array, instead of a linked list, chromoList is accessed quicker.

//	  The array "chroAtLoc" and the linked list "chroTemp" serve as a bridge to
//	  obtain the detailed list chromoList. They will be eliminated thereafter.
//	  The list locList will still be around.

	*nChromo = 0;
	*unknown = nlocUsed;
	char *chromo = (char*) malloc(sizeof(char)*LEN_LOCUS);
	char *chromo0 = (char*) malloc(sizeof(char)*LEN_LOCUS);
	char *locus = (char*) malloc(sizeof(char)*LEN_LOCUS);
	// this array is to notify which locus is done, so don't have to search
	int *done = (int*) malloc(sizeof(int)*nlocUsed);
	for (n=0; n<nlocUsed; n++) *(done+n) = -1;
	*chromo0 = '\0';
	if (chroInp == NULL) return NULL;
	int len = 0;
	int locSeen = 0;	// to count number of loci in genotype (to be run)
						// that appear in chromsomes/loci file input.

	// The first purpose of the next loop is to assign field chromo of locList
	// to be the first string in each input line. Originally, this field
	// in locList is empty. On the way, a list of chromosome, chroAtLoc,
	// will also be created; this list contains only the names. This list is
	// used to identify the numbering (in genotype input file) of the locus
	// that belongs to the chromosome, the numbering is stored in array "done".

	// The second purpose is to create a linked list "chroTemp", whose elements
	// correspond to distinct chromosomes. The field "name" is assigned the
	// chromosome name, and the field "nloci" is the number of loci in that
	// chromosome.
	for ( ; ; ) {	// each round of this loop is to work with one input line
	// each line contains 1 chromo and 1 locus (comma is allowed to separate)
		// the first string is supposed to be the name of a chromosome
		m = GetToken(chroInp, chromo, LEN_LOCUS, BLANKS, CHARSKIP, &c, &n);
		if (m == 0 || c == '\n') break;	// c = '\n': only 1 string on the line
		// the second string "locus" is supposed to be the name of a locus
		if (GetToken(chroInp, locus, LEN_LOCUS, BLANKS, CHARSKIP, &c, &n) == 0)
			break;
		// finish the line:
		for (; (c=fgetc(chroInp)) == EOF || c !='\n';);
		// the next loop is to find a locus in locList that has the name
		// as the second string, then assign field chromo of that locus to be
		// the first string. Once a locus is assigned this field, it will not
		// be assigned again if later there is an input line where the second
		// string is this locus name. Thus, in the input file, two distinct
		// lines should have distinct second strings.
		// In the process, only chromosomes having loci in genotype input
		// (which are to be used) will be collected (so no chromosome has
		// empty locus list).
		// Measure the lengths of chromosomes so that the "unknown" one
		// can be named appropriately, len will be the maximum length.
		// (The "unknown" chromosome is for collecting the rest of loci in
		// genotype input that do not belong to any chromosome read here.)
		for (n=0; n<nlocUsed; n++) {
			if (done[n] != -1) continue;	// this locus on locList passed
			// curr is the latest node added.
			if (strcmp(locus, locList[n].name) == 0) {	// "locus" on locList
				locSeen++;
				strcpy(locList[n].chromo, chromo);	// register chromosome
				strcpy((chroAtLoc[n].name), chromo);	//  for this locus
				done[n] = locList[n].num;	// locList[n] is done!
				if (strcmp(chromo0, chromo) == 0) {	// same as previous one
					(curr->nloci)++;
					break;	// save time, don't need to search the list
				}
				// check list "chroTemp" to see if this chromo is already there
				curr = chroTemp;
				prev = curr;
				while (curr != NULL) {
					if (strcmp(curr->name, chromo) == 0) break;
					prev = curr;
					curr = curr->next;
				}
				// either curr = NULL: chromo is new, or curr is the node
				// where name is this chromo.
				// If this chromo is new, create another node.
				// If already there, increase number of loci for the found node.
				if (curr == NULL) {	// add this chromo to the list chroTemp
					if (len < strlen(chromo)) len = strlen(chromo);
					newptr = (CHROPTR) malloc(sizeof(struct chro));
					if (newptr == NULL) return 0;	// error, out of memory
					if (prev != NULL) prev->next = newptr; // when list not empty
					strcpy(newptr->name, chromo);
					newptr->nloci = 1;
					newptr->next = NULL;
					curr = newptr;	// for next line of input, if the same chromo
					if (chroTemp == NULL) chroTemp = newptr;	// when list was empty
				// increase number of chromosomes
					num++;
				} else (curr->nloci)++;
				strcpy(chromo0, chromo); // so that chromo0 is one on the list
				break;	// a match has been found, no more search
			}	// end of "if (strcmp(locus, locList[n].name) == 0)"
		}	// end of "for (n=0; n<nlocUsed; n++)", search for
	}	// end of reading input "for ( ; ; )"
	free (chromo0);
	free (chromo);
	free (locus);
// if there are loci (to be used) on genotype input data, but not on this
// input file, then locSeen < nlocUsed, we will collect all remaining loci
// in locList, and put them under an unknown chromosome, named '99 .. 9',
// whose length is the maximum length of all chromosome lengths.
	if (locSeen < nlocUsed) num++;	// num will be the number of chromosomes
	chromoList = (struct chromosome*) malloc(sizeof(struct chromosome)*num);
	if (chromoList != NULL)
	{
		// go from top to bottom of list "chroTemp" to assign names of
		// chromosomes to each element of the array (list) chromoList
		curr = chroTemp;
		n = 0;
		while (curr != NULL) {
			strcpy (chromoList[n].name, curr->name);
			m = curr->nloci;
			chromoList[n].nloci = m;
			chromoList[n].locus = (int*) malloc(sizeof(int)*m);
			c = 0;
			for (p=0; p<nlocUsed; p++) {
				if (strcmp(chroAtLoc[p].name, chromoList[n].name) == 0) {
					(chromoList[n].locus)[c] = done[p];
					c++;
				}
			}
			curr = curr->next;
			n++;
		}
		// now add "unknown" chromosome if necessary, i.e., there are loci
		// that are not assigned chromosome in file input.
		m = nlocUsed - locSeen;
		*unknown = m;
		if (m > 0) {	// n is actually = num-1
			for (c=0; c<len; c++) (chromoList[n].name)[c] = '9';
			// just in case there is a chromosome named '9...9', add an 'X'
			if (len < LEN_LOCUS) (chromoList[n].name)[len] = 'X';
			chromoList[n].locus = (int*) malloc(sizeof(int)*m);
			c = 0;
			for (p=0; p<nlocUsed; p++) {
				if (done[p] == -1) {
					(chromoList[n].locus)[c] = locList[p].num;
					c++;
				}
			}
			chromoList[n].nloci = c;	// = m
		}
	}
	free(chroAtLoc);
	free(done);
	// dispose chroTemp:
	prev = chroTemp;
	for (; prev != NULL; prev = curr) {
		curr = prev->next;
		free (prev);
	}

	*nChromo = num;
	return chromoList;
}


//------------------------------------------------------------------
int RunOption (char misFilSuf[], char LocSuf[], char BurSuf[],
				char hasOpt, char rem, char *FileOne, char *FileTwo)
{

	int maxSamp, nCrit, nGeneration;
	char misDat, param, nonparam, mLD, mHet, mNomura, mTemporal, mating;
	float critVal[MAXCRIT];
	float timeline[MAXGENERATION];

	char append, format;
	char xOutLD, xOutHet, xOutCoan, xOutTemp;
	int nlocUse, nPop, maxMobilVal, lenM;	// lenM = number of digits in alleles
	int popStart = 1;	// run starts from this population.
	int popEnd = MAX_POP;
	// # pops that outLoc, outBurr produce.
	// if nonzero, will produce corresponding file
	// For Burrow, topBCrit is the number of highest crit values for output
	int popLoc1=0, popBurr1=0, topBCrit = MAXCRIT;
	int popLoc2=0, popBurr2=0;
// nloci: number of loci:
	int n, p, nlocDel, nloci=0;
	int topCrit = MAXCRIT;
	char *locUse, byRange;
	char *missFileName, *inpName, *prefix, *outName, *outLocName, *outBurrName;
	char *outFile, *outFile0, *outFolder;
// the file info is for reading info by function InfoDirective.
// Need to keep it open so can continue reading when doing populations:
// it's not closed in function InfoDirective so that function RunPop can
// access this info file.
	FILE *info;
	FILE *input;
	FILE *output;
	FILE *shOutputLD = NULL;
	FILE *shOutputHet = NULL;
	FILE *shOutputCoan = NULL;
	FILE *shOutputTemp = NULL;
	FILE *outLoc=NULL, *outBurr=NULL;
	char mode[2] = "w";
	// add July 2013
	char tabX = 0;
//	FILE *locFile = NULL;


	AGEPTR *ageSeq;
	int nSeq, nPlan, census;
	int tempClue, tempxClue;
	int totPop = 0, totPairTmp = 0;
// Added in Mar/Apr 2015
	struct locusMap *locList = NULL;
	struct chromosome *chromoList = NULL;
	FILE *chroInp = NULL;
	int chroGrp = 0;
	int nChromo = 0;
	int unknown;	// number of loci whose chromosome are unknown
// add Jan 2015/ Apr 2015:
	char *inpFolder;
	char sepBurOut, moreCol, BurAlePair;
	char *chrofileName  = (char *) malloc(LENFILE * sizeof(char));
	*chrofileName  = '\0';

	inpName  = (char *) malloc(LENFILE * sizeof(char));
	*inpName  = '\0';
	outName = (char *) malloc(PATHFILE * sizeof(char));
	*outName = '\0';
	inpFolder = (char *) malloc(LENDIR * sizeof(char));
	*inpFolder = '\0';
	outFolder = (char *) malloc(LENDIR * sizeof(char));
	*outFolder = '\0';
	outFile = (char *) malloc(PATHFILE * sizeof(char));
	*outFile = '\0';
	outFile0 = (char *) malloc(PATHFILE * sizeof(char));
	*outFile0 = '\0';

	// FileOne is the name of info directive file
	if ((info = fopen (FileOne, "r")) == NULL) {
		perror (FileOne);	// inform the file does not exist
		return 0;			// 0 file is run
	};// else
	SetDefault (&maxSamp, &param, &nonparam, &mLD,
				&mHet, &mNomura, &mTemporal, &nCrit,
				critVal, &nGeneration, &mating, timeline);
	ageSeq = (AGEPTR*) malloc(sizeof(AGEPTR));
	*ageSeq = NULL;
	if ((input = InfoDirective (&mLD, &mHet, &mNomura, &mTemporal,
						FileOne, &format, &nCrit, critVal,
						&mating, inpFolder, inpName, outFolder, outName, &nPop,
						&nloci, &maxMobilVal, &lenM,
						info, &append, ageSeq, &nSeq, &tempClue, &nPlan)) == NULL) {
		fclose (info);
		if (rem == 1) {
			remove (FileOne);
			if (hasOpt == 1) remove (FileTwo);
		};
		return 0;
	};
	// close the file before remove
	fclose (info);
	if (rem == 1) remove (FileOne);
// inform to console:
	printf ("Input file: %s -", inpName);
	if (format == FSTAT) printf (" FSTAT format");
	else if (format == GENPOP) printf (" GENEPOP format");
	printf ("\n");
	printf ("Number of loci = %d, %d-digit alleles\n", nloci, lenM);
	nlocDel = 0;
	locUse = (char*) malloc(sizeof(char)*nloci);
	// default: all loci will be used.
	for (p=0; p<nloci; p++) *(locUse+p) = 1;
// in ***Directive function, nPop was assigned default to be MAX_POP
//	if (format == GENPOP) nPop = MAX_POP;
// the next code was also taken care at ***Directive:
//	if (popEnd <= 0 || popEnd > nPop) popEnd = nPop;

	// read input, pass locus names and print on screen:
	n = mLD + mHet + mNomura + mTemporal;
	if (n <= 0) {
		printf ("No method to run\n");
		return 0;
	};
	PrtMethod (n, mLD, mHet, mNomura, mTemporal);

	if (hasOpt == 1) {
		n = OptDirective (FileTwo, &xOutLD, &xOutHet, &xOutCoan, &xOutTemp,
						&maxSamp, &popStart, &popEnd, nPop,
						&popLoc1, &popLoc2, &popBurr1, &popBurr2, &topBCrit,
						&misDat, &param, &nonparam, nloci, locUse, &nlocDel,
						&tempxClue, &byRange, &topCrit, &tabX,
// add Jan 2015/ Apr 2015:
//						inpFolder, chroInp, &chroGrp,
						inpFolder, chrofileName, &chroGrp,
						&sepBurOut, &moreCol, &BurAlePair);
		if (n==-1) {
			perror (FileTwo);
			printf ("Program runs with default options\n");
		} else if (rem == 1) remove (FileTwo);
		if (mLD == 0) {
			xOutLD = 0;
		// since no LD method, no file for Burrow coeff.
			popBurr1 = 0;
			popBurr2 = 0;
		};
		if (mHet == 0) xOutHet = 0;
		if (mNomura == 0) xOutCoan = 0;
		if (mTemporal == 0) xOutTemp = 0;
	};
	// read input, pass locus names and print on screen:
//	locFile = tmpfile();
	nlocUse = nloci-nlocDel;
//	if (GetLocUse (input, nloci, locUse, nlocUse, locFile) != 0) exit(1);

// April 2015: use array instead of temporary file for loci storage
	locList = (struct locusMap *) malloc(nlocUse * sizeof(struct locusMap));
	for (p=0; p < nlocUse; p++) {
		(locList[p].name)[0] = '\0';
		(locList[p].chromo)[0] = '\0';
	}
	if (GetLocUsed (input, nloci, locUse, nlocUse, locList) != 0) {
		printf ("Error when trying to collect locus names\n");
		exit(1);
	}
// Added April 2015:
	// only need to work on file chroInp on chromosomes if chroGrp = 1 or 2,
	// which was determined from function ChromoInp that opened file chroInp
	unknown = nlocUse;
	if (chroGrp == 1 || chroGrp == 2) {
		chroInp = GetInp (inpFolder, chrofileName);
		chromoList = GetChromo (chroInp, nlocUse, locList, &nChromo, &unknown);
		fclose(chroInp);
	}
	free (chrofileName);

// Apr 2015: this block (between 2 dotted lines) was brought here from above
// OptDirective call -------------------------------------------------------
// assign outFile to be output name, including path
	*outFile = '\0';
// open file for output:
	outFile = strcat (outFile, outFolder);
	outFile = strcat (outFile, outName);
	if (append > 0) mode[0] = 'a';
	if ((output=fopen (outFile, mode)) == NULL) {
		printf ("Output file cannot be opened! Program aborted.\n");
		exit (EXIT_FAILURE);
	} else {
		printf ("Outputs are written to file %s", outName);
		if (append > 0) printf (" (append)\n");
		printf ("\n");
	};
// -------------------------------------------------------------------------
	if (xOutLD==1) {
		prefix  = (char *) malloc(LENFILE * sizeof(char));
		*prefix = '\0';
		// prefix becomes the name
		// for short output: used only for printing to console below
		GetXoutName (prefix, outName, LENFILE, XFILSUFLD, PATHCHR);
	// assign outFile0 to be short output name, including path
		GetXoutName (outFile0, outFile, PATHFILE, XFILSUFLD, "\0");
		// FileTwo: full name for short output
		if ((shOutputLD = fopen (outFile0, mode)) != NULL) {
			printf ("Tabular-format LD Output File Name: %s\n", prefix);
		};
		free (prefix);
	};
	if (xOutHet==1) {
		prefix  = (char *) malloc(LENFILE * sizeof(char));
		*prefix = '\0';
		// prefix becomes the name
		// for short output: used only for printing to console below
		GetXoutName (prefix, outName, LENFILE, XFILSUFHET, PATHCHR);
	// assign outFile0 to be short output name, including path
		GetXoutName (outFile0, outFile, PATHFILE, XFILSUFHET, "\0");
		// FileTwo: full name for short output
		if ((shOutputHet = fopen (outFile0, mode)) != NULL) {
			printf ("Tabular-format Het-Excess Output File Name: %s\n", prefix);
		};
		free (prefix);
	};
	if (xOutCoan==1) {
		prefix  = (char *) malloc(LENFILE * sizeof(char));
		*prefix = '\0';
		// prefix becomes the name
		// for short output: used only for printing to console below
		GetXoutName (prefix, outName, LENFILE, XFILSUFCOAN, PATHCHR);
	// assign outFile0 to be short output name, including path
		GetXoutName (outFile0, outFile, PATHFILE, XFILSUFCOAN, "\0");
		// FileTwo: full name for short output
		if ((shOutputCoan = fopen (outFile0, mode)) != NULL) {
			printf ("Tabular-format Coancestry Output File Name: %s\n", prefix);
		};
		free (prefix);
	};
	if (xOutTemp==1) {
		prefix  = (char *) malloc(LENFILE * sizeof(char));
		*prefix = '\0';
		// prefix becomes the name
		// for short output: used only for printing to console below
		GetXoutName (prefix, outName, LENFILE, XFILSUFTEMP, PATHCHR);
	// assign outFile0 to be short output name, including path
		GetXoutName (outFile0, outFile, PATHFILE, XFILSUFTEMP, "\0");
		// FileTwo: full name for short output
		if ((shOutputTemp = fopen (outFile0, mode)) != NULL) {
			printf ("Tabular-format Temporal Output File Name: %s\n", prefix);
		};
		free (prefix);
	};
	if (maxSamp <= 0) maxSamp = MAX_SAMP;
	if (popLoc1 < 0) {	// if the first pop for Freq output is negative,
		popLoc1 = 1;	// assume all pops will have freq output
		popLoc2 = nPop;
	};
	if (popBurr1 < 0) {	// similar to freq.output
		popBurr1 = 1;
		popBurr2 = nPop;
	};
	// Assign auxiliary file names based on inpName without pathnames:
	outLocName = (char *) malloc(LENFILE * sizeof(char));
	*outLocName  = '\0';
	outBurrName = (char *) malloc(LENFILE * sizeof(char));
	*outBurrName  = '\0';
	missFileName = (char *) malloc(LENFILE * sizeof(char));
	*missFileName  = '\0';
	if (misDat != 0)
		GetXoutName (missFileName, inpName, LENFILE, misFilSuf, PATHCHR);
	GetXoutName (outLocName, inpName, LENFILE, LocSuf, PATHCHR);
	GetXoutName (outBurrName, inpName, LENFILE, BurSuf, PATHCHR);

// Since only populations from popStart to popEnd are analyzed, need to adjust:
	if (popLoc1 < popStart) popLoc1 = popStart;
	if (popLoc2 < popLoc1) popLoc2 = 0;
	if (popLoc1 > popEnd) popLoc2 = 0;
	if (popLoc2 > popEnd) popLoc2 = popEnd;
	// so, unless popLoc2 = 0, we have popStart <= popLoc1 <= popLoc2 <= nPop
	if (popBurr1 < popStart) popBurr1 = popStart;
	if (popBurr2 < popBurr1) popBurr2 = 0;
	if (popBurr1 > popEnd) popBurr2 = 0;
	if (popBurr2 > popEnd) popBurr2 = popEnd;
	if (popLoc2 > 0)
// open Freq. data file if the range of populations for Freq. data overlaps
// with the one in analysis (which is from popStart to popEnd <= nPop).
// However, since we assign nPop = MAX_POP in GENPOP format, not true value
// (can do it in InfoDirective if needed), the condition may not be tight
	{
		*outFile = '\0';
		outFile = strcat (outFile, outFolder);
		if ((outLoc = fopen (strcat(outFile, outLocName), "w")) != NULL) {
			PrtVersion (outLoc);
			fprintf (outLoc, "Input File: %s\n\n", inpName);
		} else popLoc2 = 0;
	};
	// change in Nov 2014/ Jan 2015:
	// add condition, so that Burrows file will not be created here, but
	// in RunPop0 for separate Burrows file when sepBurOut = 1
//	if (popBurr2 > 0)
	if (popBurr2 > 0 && sepBurOut == 0)
	{
		*outFile = '\0';
		outFile = strcat (outFile, outFolder);
		if ((outBurr = fopen (strcat(outFile, outBurrName), "w")) != NULL) {
//			if (NOEXPLAIN != 1) {
			if (sepBurOut != 1) {
				PrtVersion (outBurr);
				fprintf (outBurr, "Input File: %s\n\n", inpName);
			// print up to 100 locus names to outBurr
				PrtLocUsed (locList, outBurr, nloci, locUse, nlocUse, 100);
// Added Apr 2015:
				PrtChromo (outBurr, chromoList, nChromo, chroGrp, unknown);
			}
		} else popBurr2 = 0;
	}
//	if (popLoc2+popBurr2 == 0) {
//		fclose(locFile);
//		locFile = NULL;
//	}
	// limit number of pops that can have outputs in those auxiliary
	// files, to prevent the program from producing too big files:
	if ((popLoc2 - popLoc1) >= MAXLOCPOP)
		popLoc2 = popLoc1 + MAXLOCPOP - 1;
	if ((popBurr2 - popBurr1) >= MAXBURRPOP)
		popBurr2 =	popBurr1 + MAXBURRPOP - 1;
	*outFile = '\0';
	outFile = strcat (outFile, outFolder);
	// add outFolder, locList to list of parameters in 2015
	RunPop (0, inpName, input, 0, output, outFolder, locList,
				outLoc, outLocName, outBurr, outBurrName, shOutputLD,
				shOutputHet, shOutputCoan, shOutputTemp,
				popLoc1, popLoc2, popBurr1, popBurr2, topBCrit, popStart, popEnd,
				maxSamp, nloci, lenM, maxMobilVal, nCrit, critVal, format,
//				param, nonparam, locUse, mating, missFileName,
				param, nonparam, locUse, mating, strcat(outFile, missFileName),
				INFINITE, LEN_BLOCK, mLD, mHet, mNomura,
				mTemporal, nGeneration, timeline, ageSeq, nSeq, 1,
				// last 1 is for attempt to create sequence of generations
				// tempClue is for determine which temporal methods for
				// main output, extra output
				tempClue, tempxClue, byRange, topCrit, nPlan, census,
				&totPop, &totPairTmp, 0, tabX,
				// add parameters Apr 2015
				sepBurOut, moreCol, BurAlePair,
				chromoList, nChromo, chroGrp, unknown);

	// close the file before remove
//	fclose (info);
//	if (rem == 1) remove (FileOne);
	free (inpName);
	free (inpFolder);
	free (outName);
	free (outFolder);
	free (outLocName);
	free (outBurrName);
	free (missFileName);
	free (outFile);
	free (outFile0);
// Apr 2015:
	if (locList != NULL) free (locList);
	if (chromoList != NULL) RmChromo (chromoList);
	return 1;
}


//--------------------------------------------------------------------------
FILE *GetOutput (char desc[], char *outName, char append) {
	FILE *output = NULL;
//	printf ("%s: [%s]", desc, outName);
//	if (append == 1) printf (" (Append)");
	printf ("%s", desc);
	if (append == 1) printf (" (Append)");
	printf (": [%s]\n", outName);
	if (append == 1) output = fopen (outName, "a");
	else output = fopen (outName, "w");
	if (output == NULL) printf ("Cannot open file %s\n", outName);
	return output;
}

// --------------------------------------------------------------------------

void PrtLimitCommon (FILE *output, char byRange, int *locRanges, int nRanges,
					int popStart, int popEnd, int maxSamp, char term[])
// If there are limits on samples/pop, restriction on populations or loci,
// then print those to common outputs (this is called in RunMultiCommon).

{
	int k;
	if (output == NULL) return;

	if (popEnd < MAX_POP) {	// populations to run are limited
		if (popStart == 1) {
			if (popEnd == 1) fprintf (output, "Only run for %s 1\n", term);
			else fprintf (output, "Run up to %s %d \n", term, popEnd);
		} else {
			if (popStart < popEnd) fprintf (output,
				"Limit to %ss from %d to %d \n", term, popStart, popEnd);
			else fprintf (output, "Only run for %s %d\n", term, popEnd);
		};
	} else if (popStart > 1) {
		fprintf (output, "Run from %s %d \n", term, popStart);
	};

	if (maxSamp < MAX_SAMP) fprintf (output,
			"Up to %d individuals are processed per %s.\n", maxSamp, term);
	if (byRange == 1) {
		fprintf (output, "Run with Loci in Range");
		if (nRanges > 1) fprintf (output, "s");
		fprintf (output, ": ");
		for (k=0; k<nRanges; k++) {
			if (k > 0) fprintf (output, ", ");
			fprintf (output, " %d - %d", locRanges[2*k], locRanges[2*k+1]);
		};
		fprintf (output, "\n");
	};
	fflush (output);
}



//--------------------------------------------------------------------------
int RunMultiCommon (char *mFileName)
// Run multiple input files having the same options:
// ("applicable" means that the line exists only if needed!)
// 1. Method(s)
// 2. Number of (positive) critical values
// 3. Critical values if applicable (i.e., number on line 2 is positive)
// 4. Plan/Generations for Temporal if applicable.
// 5. Tabular-format output file(s)
// 6. CI or not (parameter and non-parameter)
// 7. Random/Monogamy for LD if applicable
// 8. Max individuals per pop.
// 9. Population range
// 10. Locus ranges to run
// 11. Common output file name
// Then the rest, at each line, are input file names including paths
// if necessary.
//
// Return the number of successful runs.
// If input file name is empty, quit.
{
	#define LOCRANGE 100	// max endpoints of locus ranges to be used
	#define MAXLOCI 1000000		// maximum number of loci to be used

	int nRanges;	// Effective number of ranges
	int locRanges[LOCRANGE];	// Endpoints for ranges of loci
	time_t rawtime;
	int maxSamp, popEnd, popStart, nCrit, nGeneration;
	int totPop = 0, totPairTmp = 0;
	char param, nonparam, mLD, mHet, mNomura, mTemporal, mating;
	float critVal[MAXCRIT];
	float timeline[MAXGENERATION+1];	// add 1 for census size

	char xOutLD = 0;
	char xOutHet = 0;
	char xOutCoan = 0;
	char xOutTemp = 0;
//	int popLoc=0, popBurr=0;	// no outputs for Freq data, Burrows coeffs..

	char format;
	char *locUse;

	FILE *mInpFile = NULL;
	char append;
	int i, n, count, c, nloci, nPop, maxMobilVal, lenM, k;
	int line;
	int nlocUse;
	int tempClue, tempxClue;
	int topCrit = MAXCRIT;
	int census = 0;
	int nPlan = 0;
// add July 2013
	char tabX;
	int *xClues;
	char byRange = 1;
//	char prefix[LENFILE] = "\0";
//	char prefix[PATHFILE] = "\0";
	FILE *input = NULL;
	FILE *output = NULL;
	FILE *shOutputLD = NULL;
	FILE *shOutputHet = NULL;
	FILE *shOutputCoan = NULL;
	FILE *shOutputTemp = NULL;
	char *prefix = (char*) malloc (sizeof(char)*PATHFILE);
	char *inpName = (char*) malloc (sizeof(char)*PATHFILE);
	char *outName = (char*) malloc (sizeof(char)*PATHFILE);
	*prefix = '\0';
	*inpName = '\0';
	*outName = '\0';
//	char prefix[PATHFILE];
//	char inpName[PATHFILE];
//	char outName[PATHFILE];
//	prefix[0] = '\0';
//	inpName[0] = '\0';
//	outName[0] = '\0';
	if ((mInpFile = fopen (mFileName, "r")) == NULL) {
		perror (mFileName);
		return 0;
	};

	count = 0;
	line = 0;
	xOutLD = 0;
	xOutHet = 0;
	xOutCoan = 0;
	xOutTemp = 0;
	tempClue = 0;
	tempxClue = 0;
	topCrit = MAXCRIT;
	append = 0;	// 0 for overwrite output, 1 for append, assumed overwrite
	shOutputLD = NULL;	// assumed no short output file
	shOutputHet = NULL;	// assumed no short output file
	shOutputCoan = NULL;	// assumed no short output file
	shOutputTemp = NULL;	// assumed no short output file
	nGeneration = 0;
	SetDefault (&maxSamp, &param, &nonparam, &mLD,
				&mHet, &mNomura, &mTemporal, &nCrit,
				critVal, &nGeneration, &mating, timeline);

	popStart = 1;
	popEnd = MAX_POP;
// insert in Dec 2011, to read methods and critical values first,
// then generation timeline if temporal is desired
	if (FindMethod(mInpFile, mFileName, &line, &mLD, &mHet, &mNomura,
		&mTemporal, &tempClue) < 0) return 0;
	line++;
	// read critical values
	if ((n=CritValRead (mInpFile, MAXCRIT, critVal, &i)) <= 0) {
		ErrMsg (mFileName, "ERROR on Number of Critical Value", line);
		return 0;
	} else nCrit = n;
	if (i > 0) line++;
	if (mTemporal == 1) {
		line++;
		if ((n = GeneratnRead (mInpFile, &nGeneration,
			timeline, MAXGENERATION+1, &census)) <= 1) {
			ErrMsg (mFileName, "No valid generation timeline for temporal!",
					line);
			mTemporal = 0;
		} else {	// good generation set
		// the first index in timeline is for census size, which is
		// assigned to variable census; so need to move back one index:
		// (actually upto index nGeneration is enough)
			for (k=0; k< nGeneration; k++) timeline[k] = timeline[k+1];
			if (census > 0) nPlan = 2;
			else nPlan = 1;
		};
	};
	line++;
// ------------------------------------------------------------------
// extra output and its option: (modified July 2013 for Tab-delimiter)
	tempxClue = 0;
	xClues = (int*) malloc(sizeof(int)*4);
// default values for methods to have extra outputs, which temporal,
// and which critical values:
	*xClues = 0; *(xClues+1) = 0; *(xClues+2) = MAXCRIT; *(xClues+3) = TABX;
	k = GetClues(mInpFile, xClues, 4, 1);
	n = *xClues;
	if (k <= 0) {	// no number read!
		ErrMsg (mFileName, "At reading clue for tabular-format output!",
					line);
		return 0;
	} else {
		SetMethod (n, &xOutLD, &xOutHet, &xOutCoan, &xOutTemp);
		if (mLD == 0) xOutLD = 0;
		if (mHet == 0) xOutHet = 0;
		if (mNomura == 0) xOutCoan = 0;
		if (mTemporal == 0) xOutTemp = 0;
	};
	tempxClue = *(xClues+1);
	topCrit = (*(xClues+2) < 0)? MAXCRIT: *(xClues+2);
	// July 2013:
	tabX = (*(xClues+3) == 0)? 0: 1;
	free(xClues);

// only one line for the option of having CI (param & jackknife if applicable)
	line++;			// having CIs or not
	n = 1;	// default: have CIs
	GetInt (mInpFile, &n, 1);
	param = (n != 0)? 1: 0;
	nonparam = (n != 0)? 1: 0;
// only read mating model if LD is included in the methods:
	if (mLD == 1) {
		line++;		// read mating model: 0 for random, 1 for monogamy
		n = 0;		// default: random mating
		GetInt (mInpFile, &n, 1);
		mating = (n != 0)? 1: 0;
	};
//-------------------------------------------------------------------


// Comment out if no options for max individual per pop., or range of pops,
// or range of loci.
//*
	line++;	// max samples per pop:
	GetInt (mInpFile, &maxSamp, 1);
	if (maxSamp <= 0) maxSamp = MAX_SAMP;

	line++;	// range of populations to run:
	k = 0;
	i = GetPair(mInpFile, &k, &n, 1);	// i = 0: no number or negative for
										// the first entry, then k still = 0
	if (k > 0) {	// No declaration of error if i=0. If k > 0, then i >=1
		popEnd = k;
		if (i == 2) {	// then n is obtained a value from GetPair
			if (n >= k) {
				popStart = k;		// k <= n, so the range is from k to n.
				popEnd = n;
			};	// else, n < k: the second entry is erroneous, popEnd is not
				// reassigned, the range is from popStart = 1 to popEnd = k
		};	// else, only 1 number: pop. range is from popStart = 1 to k
	};	// else, the first entry is <= 0, assuming range from 1 to MAX_POP
// Now for the raanges of loci
	line++;
	nRanges = GetRanges (mInpFile, locRanges, LOCRANGE, MAXLOCI, &byRange);

/*
	for (i=0; i<LOCRANGE; i++) *(locRanges+i) = 0;
	n = GetClues(mInpFile, locRanges, LOCRANGE, 1);

	// (GetClues assign values >=0 for locRanges, n = number of entries read)
// Expect to read ranges of loci in pairs, so in normal case, n is even.
// If n = 0: no number given, so no limit.
// If n > 0 (at least one number read) and the first one is 0: no limit
// If n = 1: one number only. If it is 0, no limit as said above.
// If it is positive, then assume one range, from 1 to that number.

// If n is odd > 1, the last number will be ignored.
// If n is even, we have nRanges = n/2 pairs. An legitimate pair should
// contain two numbers i, j with i <= j. If illegitimate (i.e., i > j),
// that pair is ignored by reassignining i = 0 = j.
// For no limit, we set nRanges = 1, *locRanges = 1, *(locRanges+1) = MAXLOCI
	nRanges = n/2;
	if (*locRanges == 0) {	// if n = 0, (*locRanges == 0) by default
		nRanges = 1;
		*locRanges = 1;
		*(locRanges+1) = MAXLOCI;
		byRange = 0;
	} else {
		if (n == 1) {	// the case the first one 0 was excluded
			nRanges = 1;
			*(locRanges+1) = *locRanges;
			*locRanges = 1;
		} else {
			for (k = 0; k < nRanges; k++) {
				if (*(locRanges + 2*k) > *(locRanges + 2*k + 1)) {
					*(locRanges + 2*k) = 0;
					*(locRanges + 2*k + 1) = 0;
				};
			};
		};
		byRange = 1;
	};


//*/
//-------------------------------------------------------------------

	line++;
	// read output name, if not given, quit
	append = 0;
	if (GetToken (mInpFile, outName, PATHFILE, BLANKS, ENDCHRS, &c, &k) <= 0)
	{
		ErrMsg (mFileName, "Output File name must be given!", line);
		return 0;
	}
	else {	// need to go to next line, preparing for reading input
// added Feb 15 2013: (if name ends by SPECHR, output files to be appended)
		if (c == SPECHR && k == 0) append = 1;
		for (; (c=fgetc(mInpFile)) != EOF && c!='\n';);
	};

	n = mLD + mHet + mNomura + mTemporal;
	if (n <= 0) {
		printf ("No method to run!\n");
		return 0;
	};
	PrtMethod (n, mLD, mHet, mNomura, mTemporal);
	// now, line is the number of lines read from the beginning up to
	// the name of output file.
	count = 0;	// count legitimate input files listed
	for ( ; ; ) {
		line++;
	// get input file names, and count the legitimate ones
		*inpName = '\0';
		if ((n=GetToken (mInpFile,inpName,PATHFILE,BLANKS,ENDCHRS,&c, &k)) <= 0)
		// string for input name is empty, just quit - no error messages
		{
//			ErrMsg (mFileName, "At reading input name", line);
//			printf("%s\n", inpName);
			break;
		}
// If name is ended by char, then c = '\n', but cursor is still
// on the line of that name, so we still need to go to next line
		for (; (c=fgetc(mInpFile)) != EOF && c!='\n';);
		if ((input = fopen (inpName, "r")) == NULL) {
			printf ("\nERROR in open file %s\n", inpName);
			perror(inpName);
			if (c == EOF) break;
			continue;
		};
		printf ("\n>>> Input %d: [%s], ", count+1, inpName);
		format = FSTAT;
		// n is the length of inpName without trailing blanks
		if (n > 4) {	// assuming GENEPOP if extension is "gen".
			if (*(inpName+n-4) == '.' && tolower(*(inpName+n-3)) == 'g' &&
				tolower(*(inpName+n-2))=='e' && tolower(*(inpName+n-1))=='n')
				format = GENPOP;
		};	// format now will be checked

// now read input file to obtain necessary parameter values
		// if input file extension is not "gen", it was assumed FSTAT format
		if (format == FSTAT) {
			if (GetInfoDat(input,
						&nPop, &nloci, &maxMobilVal, &lenM, LEN_BLOCK)==0) {
				// this is not FSTAT, now assume it's GENPOP format
				format = GENPOP;
				rewind(input);
			}
		}
		if (format == GENPOP) printf ("GENEPOP format\n");
		if (format == FSTAT)  printf ("FSTAT format\n");
		if (format == GENPOP) {	// GENEPOP format
			if ((nloci = GetnLoci (input, LEN_BLOCK, &lenM)) <= 0) {
				fclose (input);
				printf ("Error in input file [%s]\n", inpName);
				continue;
			} else {
				rewind (input);
				for (; (c=fgetc(input)) !='\n';);	//position at first locus
			};
	// should assign maxMobilVal to make sure it doesn't get garbage.
	// Since the length is lenM, assign max possible value:
			for (maxMobilVal=1, i=1; i<=lenM; i++) maxMobilVal *= 10;
		};
		printf ("Number of loci = %d, %d-digit alleles\n", nloci, lenM);
		locUse = (char*) malloc (sizeof(char)*nloci);
	// default: all loci will be used.
//		for (i=0; i<nloci; i++) *(locUse+i) = 1;
	// if ranges of loci are read (as noted before the "for ( ; ; )" loop),
	// then the following loop is needed instead of the previous line of code
		for (i=0; i<nloci; i++) *(locUse+i) = 0;
		nlocUse = nloci;	// assumed no loci deleted
		for (i=0; i<nloci; i++) {
			for (k=0; k<nRanges; k++) {
				if ((locRanges[2*k]<=(i+1)) && ((i+1)<= locRanges[2*k+1]))
					*(locUse+i) = 1;
			};
		};
		for (i=0; i<nloci; i++) if (*(locUse+i) != 1) nlocUse--;
		if (GetLocUsed (input, nloci, locUse, nlocUse, NULL) != 0) {
			printf("This input file is skipped.\n");
			fclose (input);
			continue;
		};

		// only open output files if the first input file looks OK:
		if (count == 0) {
			output = GetOutput("Main Output", outName, append);
			if (output == NULL) return 0;
		// Print headlines in extra output files, including the first
		// input file, the last will be printed after this "for" loop.
		// There will be no headlines printed for those output files
		// at running each input file.

		// string outName will be used repeatedly as shorter output names,
		// also string prefix is done after used for input
			*prefix = '\0';
			strcat (prefix, outName);
			time ( &rawtime );
			if (xOutLD == 1) {
				GetXoutName (outName, prefix, PATHFILE, XFILSUFLD, PATHCHR);
				shOutputLD = GetOutput ("Tabular-format LD Output",
										outName, append);
				if (append == 1) PrtLines (shOutputLD, 60, '-');
				PrtVersion (shOutputLD);
				fprintf (shOutputLD, "Starting time: %s", ctime (&rawtime));
				PrtLimitCommon (shOutputLD, byRange, locRanges, nRanges,
								popStart, popEnd, maxSamp, "Population");
			};
			if (xOutHet == 1) {
				GetXoutName (outName, prefix, PATHFILE, XFILSUFHET, PATHCHR);
				shOutputHet = GetOutput ("Tabular-format Het. Excess Output",
										outName, append);
				if (append == 1) PrtLines (shOutputHet, 60, '-');
				PrtVersion (shOutputHet);
				fprintf (shOutputHet, "Starting time: %s", ctime (&rawtime));
				PrtLimitCommon (shOutputHet, byRange, locRanges, nRanges,
								popStart, popEnd, maxSamp, "Population");
			};
			if (xOutCoan == 1) {
				GetXoutName (outName, prefix, PATHFILE, XFILSUFCOAN, PATHCHR);
				shOutputCoan = GetOutput ("Tabular-format Coancestry Output",
										outName, append);
				if (append == 1) PrtLines (shOutputCoan, 60, '-');
				PrtVersion (shOutputCoan);
				fprintf (shOutputCoan, "Starting time: %s", ctime (&rawtime));
				PrtLimitCommon (shOutputCoan, byRange, locRanges, nRanges,
								popStart, popEnd, maxSamp, "Population");
			};
			if (xOutTemp == 1) {
				GetXoutName (outName, prefix, PATHFILE, XFILSUFTEMP, PATHCHR);
				shOutputTemp = GetOutput ("Tabular-format Temporal Output",
										outName, append);
				if (append == 1) PrtLines (shOutputTemp, 60, '-');
				PrtVersion (shOutputTemp);
				fprintf (shOutputTemp, "Starting time: %s", ctime (&rawtime));
				PrtLimitCommon (shOutputTemp, byRange, locRanges, nRanges,
								popStart, popEnd, maxSamp, "Sample");
			};
		};
//---------------------------------------------------------------------------

		if (RunPop (++count, inpName, input, append, output,
				NULL,	// for outFolder, not needed since no Burrows file
				NULL,	// locList
				NULL, "\0", //NULL,		// <-NULL for outLoc file and outLocName
				NULL, "\0", //NULL, 	// <-NULL for outBurr file and outBurrName
				shOutputLD, shOutputHet, shOutputCoan, shOutputTemp,
				0, 0, 0, 0, 0, popStart, popEnd, maxSamp, nloci, lenM,
				// the 0s are for Loc, Bur
				maxMobilVal, nCrit, critVal, format, param, nonparam,
				locUse, mating,
				"\0",	// empty for missFileName: no missing data file created
				INFINITE, LEN_BLOCK, mLD, mHet, mNomura, mTemporal,
				nGeneration, timeline, NULL, 0,
				0, tempClue, tempxClue, byRange, topCrit, nPlan, census,
				&totPop, &totPairTmp, 1, tabX,
				// NULL, 0: no list for Generations
				// 0 in front of tempClue is for no getting age (no list!)
				// last 1 is for running multiple files with common setting
				 0, 0, 0, NULL, 0, 0, 0) == 0)	// last param added Apr 2015
			printf("Finish running input %d.\n", count);
		free (locUse);
// these are already closed in RunPop
//		fclose (input);
//		fclose (output);
	};
//	free (locRanges);
	time ( &rawtime );
	PrintEndTime (output, rawtime);

	PrtPopRun (shOutputLD, totPop, 37);
	PrintEndTime (shOutputLD, rawtime);
	PrtPopRun (shOutputHet, totPop, 37);
	PrintEndTime (shOutputHet, rawtime);
	PrtPopRun (shOutputCoan, totPop, 37);
	PrintEndTime (shOutputCoan, rawtime);
	if (shOutputTemp != NULL) {
		PrtLines (shOutputTemp, 49, '-');
		fprintf(shOutputTemp,
			"Total number of samples of populations =%9d\n", totPop);
		fprintf(shOutputTemp,
			"Total number of sample pairs analysed  =%9d\n", totPairTmp);
		PrtLines (shOutputTemp, 49, '-');
	};
	PrintEndTime (shOutputTemp, rawtime);
	fclose (mInpFile);
	return count;

}

//------------------------------------------------------------------
