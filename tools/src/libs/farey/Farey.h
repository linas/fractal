
/* 
 * FUNCTION:
 * Convert real numbers to Farey numbers, and back
 * Coded in psuedo-OO style for easy conversion to C++
 *
 * HISTORY:
 * created Linas Vepstas January 16 1994
 * Added TFarey July 1995 -- linas
 */

#if defined (AIXV3) || defined (AIXV315) || defined (AIXV32)
#define ANSI_C 
#endif

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
extern double ContinuedFractionToCosReal (struct Farey *, double omega);
extern double ContinuedFractionToEFarey (struct Farey *, double t);
extern double ContinuedFractionToEFareyGap (struct Farey *, double t);
extern double ContinuedFractionToEFraction (struct Farey *, double t);
extern double ContinuedFractionToEReal (struct Farey *, double t);
extern double ContinuedFractionToExpReal (struct Farey *, double t);
extern double ContinuedFractionToFarey (struct Farey *);
extern double ContinuedFractionToReal (struct Farey *);
extern double ContinuedFractionToSinFarey (struct Farey *, double t);
extern double ContinuedFractionToSincReal (struct Farey *, double omega);
extern double ContinuedFractionToTFarey (struct Farey *, double t);
extern double ContinuedFractionToXFarey (struct Farey *, double t);
extern double ContinuedFractionToZReal (struct Farey *, double z);

extern double ContinuedFractionToCnReal (struct Farey *, double omega);
extern double ContinuedFractionToSnReal (struct Farey *, double omega);

extern void PrintContinuedFraction (struct Farey *);
extern void RatioToContinuedFraction (struct Farey *, int num, int deno);
extern void RealToContinuedFraction (struct Farey *, double);

extern double CFSum (struct Farey *self, struct Farey *other,
                     double alpha, double beta, double gamma);
extern double CFProd (struct Farey *self, struct Farey *other,
                     double alpha, double beta);
#else 

/* public methods */
extern struct Farey *CreateFarey ();
extern void SetReal ();
extern double GetReal ();
extern double GetFarey ();

/* protected methods */
extern double ContinuedFractionTerms ();
extern double ContinuedFractionToCosReal ();
extern double ContinuedFractionToEFarey ();
extern double ContinuedFractionToEFareyGap ();
extern double ContinuedFractionToEFraction ();
extern double ContinuedFractionToEReal ();
extern double ContinuedFractionToExpReal ();
extern double ContinuedFractionToFarey ();
extern double ContinuedFractionToReal ();
extern double ContinuedFractionToSinFarey ();
extern double ContinuedFractionToSincReal ();
extern double ContinuedFractionToTFarey ();
extern double ContinuedFractionToXFarey ();
extern double ContinuedFractionToZReal ();
extern double ContinuedFractionToCnReal ();
extern double ContinuedFractionToSnReal ();
extern void PrintContinuedFraction ();
extern void RatioToContinuedFraction ();
extern void RealToContinuedFraction ();

extern double CFSum ();
extern double CFProd ();
#endif 


/* ---------------------- END OF FILE ------------------------- */
