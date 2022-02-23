#include <R.h>
#include <Rinternals.h>
#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>



static R_NativePrimitiveArgType eType[] = {
	INTSXP,
  INTSXP,
	INTSXP,
  REALSXP,
};


extern void F77_NAME(mkerrors)(int *nind, int *nsnp, int *genofr, double *eprobfr);


static const R_FortranMethodDef FortranEntries[] = {
	{"mkerrors", (DL_FUNC) &F77_NAME(mkerrors), 4, eType},
  {NULL, NULL, 0, NULL}
};


void attribute_visible R_init_SimGE(DllInfo *info)  // attribute_visible -> error
{
  R_registerRoutines(info,
                     NULL,          // .C
                     NULL,          // .Call
                     FortranEntries, // .Fortran
                     NULL);         // .External
  R_useDynamicSymbols(info, FALSE);
	R_forceSymbols(info, TRUE);  //available from R 3.0.0
}
