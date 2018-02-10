/* for fitting models */


extern int dqrls_(double *x, int *n, int *p, double *y, int *ny, double *tol, double *b, 
		  double *rsd, double *qty, int *k, int *jpvt, double *qraux, double *work);



/* for varcov estimation*/

/********************************************************************
 **
 ** external declarations for Choleski routines (LINPACK)
 **
 **
 *******************************************************************/

extern int dpofa_(double *x, int *lda, int *n, int *j);
extern int dpodi_(double *x, int *lda, int *n, double *d, int *j);


/********************************************************************
 **
 ** external declarations for Choleski routines (LAPACK)
 **
 **
 *******************************************************************/

extern int dpotrf_(const char *uplo, const int *n, double* a, const int *lda, int *info);
extern int dpotri_(const char *uplo, const int *n, double* a, const int *lda, int *info);


/*****************************************************************
 * svd routine - LINPACK
 *****************************************************************/

extern int dsvdc_(double *x, int *ldx, int *n, int *p, double *s, double *e, double *u, int *ldu,
		 double *v, int *ldv, double *work, int *job, int *info);

/*****************************************************************
 * svd routine - LAPACK
 *****************************************************************/

extern int dgesdd_(const char *jobz,
                      const int *m, const int *n,
                      double *a, const int *lda, double *s,
                      double *u, const int *ldu,
                      double *vt, const int *ldvt,
                      double *work, const int *lwork, int *iwork, int *info);

