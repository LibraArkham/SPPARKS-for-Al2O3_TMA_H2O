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


enum{VACANCY,O,Al,OH,AlX3,AlX2,AlX,AlXOH,AlOH2,AlOH};


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
      else if (strcmp(list[i],"AlOH") == 0) which[i] = AlOH;
      else if (strcmp(list[i],"AlOH2") == 0) which[i] = AlOH2;
      else if (strcmp(list[i],"AlX") == 0) which[i] = AlX;
      else if (strcmp(list[i],"AlX2") == 0) which[i] = AlX2;
      else if (strcmp(list[i],"AlX3") == 0) which[i] = AlX3;
      else if (strcmp(list[i],"AlXOH") == 0) which[i] = AlXOH;
      else if (strcmp(list[i],"VAC") == 0) which[i] = VACANCY;
      else if (strcmp(list[i],"events") == 0) which[i] = EVENTS;
      else if (strcmp(list[i],"OH") == 0) which[i] = OH;



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
    sites[O] = 0; sites[Ta] = 0; sites[OTa] = 0;sites[VACANCY] = 0;
    sites[TaX4O] = 0; sites[TaX] = 0; sites[TaO] = 0; sites[TaX5O] = 0;  sites[OH] = 0; sites[TaX5OH] = 0; sites[TaX3O] = 0;
    int *element = appaldTMA->element;
    int nlocal = appaldTMA->nlocal;
    for (int i = 0; i < nlocal; i++) sites[element[i]]++;
  }

  for (int i = 0; i < nlist; i++) {
    if (which[i] == O) ivalue = sites[O];
    else if (which[i] == Al) ivalue = sites[Al];
    else if (which[i] == VACANCY) ivalue = sites[VACANCY];
    else if (which[i] == AlX) ivalue = sites[AlX];
    else if (which[i] == AlX2) ivalue = sites[AlX2];
    else if (which[i] == AlX3) ivalue = sites[AlX3];
    else if (which[i] == AlOH) ivalue = sites[AlOH];
    else if (which[i] == AlOH2) ivalue = sites[AlOH2];
    else if (which[i] == OH) ivalue = sites[OH];
    else if (which[i] == AlXOH) ivalue = sites[AlXOH];
    else if (which[i] == EVENTS) ivalue = appaldTMA->nevents;
    else if (which[i] == ONE) ivalue = appaldTMA->scount[index[i]];
    else if (which[i] == TWO) ivalue = appaldTMA->dcount[index[i]];
    else if (which[i] == THREE) ivalue = appaldTMA->vcount[index[i]];
    else if (which[i] == QCM) ivalue = 16*sites[O]+181*sites[Ta]+197*sites[OTa]+226*sites[TaX]+197*sites[TaO]+377*sites[TaX4O]+422*sites[TaX5O]+17*sites[OH]+423*sites[TaX5OH]+316*sites[TaX3O]+229*sites[O3Ta];

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
