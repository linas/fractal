
/* 
 * FUNCTION:
 * Convert real numbers to Farey numbers, and back
 * Coded in psuedo-OO style for easy conversion to C++
 *
 * HISTORY:
 * Linas Vepstas Januery 16 1994
 * Added TFarey July 1995 -- linas
 */

/* class Farey { */
struct Farey {
   /* public: */
   double real;       /* the number, in float pt rep. */ 

   int intpart;       /* the integer part of the real number */
   double fracpart;   /* the fractional part of the real number  */
   unsigned int num;   /* the numerator of the fractional part */
   unsigned int denom; /* the denominator of the fractional part */

   /* protected: */
   int tinued_frac[32];  /* values of continued fraction */
   int nterms;      /* number of terms in the continued fraction expansion */

   int cutoff;      /* used to avoid instability when converting
                     * floating point numbers */

};

/* function prototypes */
#ifdef ANSI_C

/* public methods */
extern struct Farey *CreateFarey (void);
extern void SetReal (struct Farey *, double);
extern double GetReal (struct Farey *);
extern double GetFarey (struct Farey *);

/* protected methods */
extern double ContinuedFractionTerms (struct Farey *);
extern double ContinuedFractionToEFarey (struct Farey *, double t);
extern double ContinuedFractionToEFareyGap (struct Farey *, double t);
extern double ContinuedFractionToEFraction (struct Farey *, double t);
extern double ContinuedFractionToEReal (struct Farey *, double t);
extern double ContinuedFractionToFarey (struct Farey *);
extern double ContinuedFractionToReal (struct Farey *);
extern double ContinuedFractionToSinFarey (struct Farey *, double t);
extern double ContinuedFractionToTFarey (struct Farey *, double t);
extern double ContinuedFractionToZReal (struct Farey *, double z);
extern void PrintContinuedFraction (struct Farey *);
extern void RatioToContinuedFraction (struct Farey *, int num, int deno);
extern void RealToContinuedFraction (struct Farey *, double);

#else 

/* public methods */
extern struct Farey *CreateFarey ();
extern void SetReal ();
extern double GetReal ();
extern double GetFarey ();

/* protected methods */
extern double ContinuedFractionTerms ();
extern double ContinuedFractionToEFarey ();
extern double ContinuedFractionToEFareyGap ();
extern double ContinuedFractionToEFraction ();
extern double ContinuedFractionToEReal ();
extern double ContinuedFractionToFarey ();
extern double ContinuedFractionToReal ();
extern double ContinuedFractionToSinFarey ();
extern double ContinuedFractionToTFarey ();
extern double ContinuedFractionToZReal ();
extern void PrintContinuedFraction ();
extern void RatioToContinuedFraction ();
extern void RealToContinuedFraction ();

#endif 


/* ---------------------- END OF FILE ------------------------- */
