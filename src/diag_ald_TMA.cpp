/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "diag_ald_TMA.h"
#include "app.h"
#include "app_ald_TMA.h"
#include "comm_lattice.h"
#include "timer.h"
#include "error.h"
#include "memory.h"

using namespace SPPARKS_NS;


enum{VACANCY,O,Al,AlX,AlXOH,AlOH2,AlOH,OH,OHAlX3,OAlX2,OAlX,OAl,OAlX2H2O,OAlXH2O,OAlXOHH2O,OAlXOH,OAlOH2,OAlOH,QCM,EVENTS,ONE,TWO,THREE};


/* ---------------------------------------------------------------------- */

DiagAldTMA::DiagAldTMA(SPPARKS *spk, int narg, char **arg) : Diag(spk,narg,arg)
{
  if (strcmp(app->style,"ald/TMA") != 0)
    error->all(FLERR,"Diag_style ald requires app_style ald");

  nlist = 0;

  int iarg = iarg_child;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"list") == 0) {
      nlist = narg - iarg - 1;
      list = new char*[nlist];
      int j = 0;
      for (int i = iarg+1; i < narg; i++) {
	int n = strlen(arg[i]) + 1;
	list[j] = new char[n];
	strcpy(list[j],arg[i]);
	j++;
      }
      iarg = narg;
    } else error->all(FLERR,"Illegal diag_style ald command");
  }

  if (nlist == 0) error->all(FLERR,"Illegal diag_style ald command");
  which = new int[nlist];
  index = new int[nlist];
  ivector = new int[nlist];
}

/* ---------------------------------------------------------------------- */

DiagAldTMA::~DiagAldTMA()
{
  for (int i = 0; i < nlist; i++) delete [] list[i];
  delete [] list;
  delete [] which;
  delete [] index;
  delete [] ivector;
}

/* ---------------------------------------------------------------------- */

void DiagAldTMA::init()
{
  appaldTMA = (AppAldTMA *) app;
  
  int none = appaldTMA->none;
  int ntwo = appaldTMA->ntwo;
  int nthree = appaldTMA->nthree;
  for (int i = 0; i < nlist; i++) {
      if (strcmp(list[i],"O") == 0) which[i] = O;
      else if (strcmp(list[i],"Al") == 0) which[i] = Al;
      else if (strcmp(list[i],"OAlOH") == 0) which[i] = OAlOH;
      else if (strcmp(list[i],"OAlOH2") == 0) which[i] = OAlOH2;
      else if (strcmp(list[i],"OAlX") == 0) which[i] = OAlX;
      else if (strcmp(list[i],"OAlX2") == 0) which[i] = OAlX2;
      else if (strcmp(list[i],"OHAlX3") == 0) which[i] = OHAlX3;
      else if (strcmp(list[i],"OAlXOH") == 0) which[i] = OAlXOH;
      else if (strcmp(list[i],"OAl") == 0) which[i] = OAl;
      else if (strcmp(list[i],"VAC") == 0) which[i] = VACANCY;
      else if (strcmp(list[i],"events") == 0) which[i] = EVENTS;
      else if (strcmp(list[i],"OH") == 0) which[i] = OH;
      else if (strcmp(list[i],"OAlX2H2O") == 0) which[i] = OAlX2H2O;
      else if (strcmp(list[i],"OAlXH2O") == 0) which[i] = OAlXH2O;
      else if (strcmp(list[i],"OAlXOHH2O") == 0) which[i] = OAlXOHH2O;
      else if (strcmp(list[i],"QCM") == 0) which[i] = QCM;
      else if (strcmp(list[i],"AlX") == 0) which[i] = AlX;
      else if (strcmp(list[i],"AlXOH") == 0) which[i] = AlXOH;
      else if (strcmp(list[i],"AlOH") == 0) which[i] = AlOH;
      else if (strcmp(list[i],"AlOH2") == 0) which[i] = AlOH2;



    else if (list[i][0] == 's') {
      which[i] = ONE;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > none) 
	error->all(FLERR,"Invalid value setting in diag_style ald");
      index[i] = n - 1;
    } else if (list[i][0] == 'd') {
      which[i] = TWO;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > ntwo) 
	error->all(FLERR,"Invalid value setting in diag_style ald");
      index[i] = n - 1;
    } else if (list[i][0] == 'v') {
      which[i] = THREE;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > nthree) 
	error->all(FLERR,"Invalid value setting in diag_style ald");
      index[i] = n - 1;
    } else error->all(FLERR,"Invalid value setting in diag_style ald");
  }

  siteflag = 1; 

  for (int i = 0; i < nlist; i++) ivector[i] = 0;
}

/* ---------------------------------------------------------------------- */

void DiagAldTMA::compute()
{
  int sites[800],ivalue;
// here as well we have to consider some modification, generally it does not seem so difficult
  if (siteflag) {
    sites[O] = 0; sites[Al] = 0; sites[OAlOH2] = 0; sites[VACANCY] = 0; sites[OAl] = 0; sites[OAlOH2] = 0; sites[AlX] = 0;sites[AlXOH] = 0; sites[AlOH] = 0; sites[AlOH2] = 0;
    sites[OAlOH] = 0; sites[OAlX] = 0; sites[OAlX2] = 0; sites[OHAlX3] = 0;  sites[OAlXOH] = 0;  sites[OAlX2H2O] = 0;  sites[OAlXH2O] = 0;  sites[OAlXOHH2O] = 0;
    int *element = appaldTMA->element;
    int nlocal = appaldTMA->nlocal;
    for (int i = 0; i < nlocal; i++) sites[element[i]]++;
  }

  for (int i = 0; i < nlist; i++) {
    if (which[i] == O) ivalue = sites[O];
    else if (which[i] == Al) ivalue = sites[Al];
    else if (which[i] == VACANCY) ivalue = sites[VACANCY];
    else if (which[i] == OAlX) ivalue = sites[OAlX];
    else if (which[i] == OAlX2) ivalue = sites[OAlX2];
    else if (which[i] == OHAlX3) ivalue = sites[OHAlX3];
    else if (which[i] == OAlOH) ivalue = sites[OAlOH];
    else if (which[i] == OAlOH2) ivalue = sites[OAlOH2];
    else if (which[i] == OH) ivalue = sites[OH];
    else if (which[i] == AlX) ivalue = sites[AlX];
    else if (which[i] == AlXOH) ivalue = sites[AlXOH];
    else if (which[i] == AlOH) ivalue = sites[AlOH];
    else if (which[i] == AlOH2) ivalue = sites[AlOH2];
    else if (which[i] == OAlXOH) ivalue = sites[OAlXOH];
    else if (which[i] == OAl) ivalue = sites[OAl];
    else if (which[i] == OAlX2H2O) ivalue = sites[OAlX2H2O];
    else if (which[i] == OAlXH2O) ivalue = sites[OAlXH2O];
    else if (which[i] == OAlXOHH2O) ivalue = sites[OAlXOHH2O];
    else if (which[i] == EVENTS) ivalue = appaldTMA->nevents;
    else if (which[i] == ONE) ivalue = appaldTMA->scount[index[i]];
    else if (which[i] == TWO) ivalue = appaldTMA->dcount[index[i]];
    else if (which[i] == THREE) ivalue = appaldTMA->vcount[index[i]];
    else if (which[i] == QCM) ivalue = 16*sites[O]+27*sites[Al]+58*sites[OAlX]+73*sites[OAlX2]+89*sites[OHAlX3]+75*sites[OAlXOH]+60*sites[OAlOH]+77*sites[OAlOH2]+17*sites[OH]+91*sites[OAlX2H2O]+76*sites[OAlXH2O]+93*sites[OAlXOHH2O]+43*sites[OAl]+42*sites[AlX]+59*sites[AlXOH]+44*sites[AlOH]+61*sites[AlOH2];

    MPI_Allreduce(&ivalue,&ivector[i],1,MPI_INT,MPI_SUM,world);
  }
}

/* ---------------------------------------------------------------------- */

void DiagAldTMA::stats(char *str)
{
  for (int i = 0; i < nlist; i++) {
    sprintf(str," %6d ",ivector[i]);
    str += strlen(str);
  }
}

/* ---------------------------------------------------------------------- */

void DiagAldTMA::stats_header(char *str)
{
  for (int i = 0; i < nlist; i++) {
    sprintf(str," %6s ",list[i]);
    str += strlen(str);
  }
}
