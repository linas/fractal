
/* 
 * FUNCTION:
 * Convert real numbers to Farey numbers, and back
 * Coded in psuedo-OO style for easy conversion to C++
 *
 * HISTORY:
 * created Linas Vepstas January 16 1994
 * Added Farey July 1995 -- linas
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
extern double ContinuedFractionToInvZReal (struct Farey *, double z);

extern double ContinuedFractionToCnReal (struct Farey *, double omega);
extern double ContinuedFractionToSnReal (struct Farey *, double omega);

extern double ContinuedFractionToFCnReal (struct Farey *, double omega);
extern double ContinuedFractionToFSnReal (struct Farey *, double omega);

extern double ContinuedFractionToZCnReal (struct Farey *, double omega, double z);
extern double ContinuedFractionToZSnReal (struct Farey *, double omega, double z);

extern void PrintContinuedFraction (struct Farey *);
extern void RatioToContinuedFraction (struct Farey *, int num, int deno);
extern void RealToContinuedFraction (struct Farey *, double);

extern double CFSum (struct Farey *self, struct Farey *other,
                     double alpha, double beta, double gamma);
extern double CFProd (struct Farey *self, struct Farey *other,
                     double alpha, double beta);



extern double Inverse (void *, double (*)(void *, double), double);

extern double InvZReal (double, double);


/* ---------------------- END OF FILE ------------------------- */
