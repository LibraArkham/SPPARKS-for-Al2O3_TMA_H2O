#include "math.h"
#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "app_ald_TMA.h"
#include "solve.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

enum{VACANCY,O,Al,AlX,AlXOH,AlOH2,AlOH,OH,OHAlX3,OAlX2,OAlX,OAl,OAlX2H2O,OAlXH2O,AlXOHH2O,OAlXOH,OAlOH2,OAlOH};


#define DELTAEVENT 100000

/* ---------------------------------------------------------------------- */

AppAldTMA::AppAldTMA(SPPARKS *spk, int narg, char **arg) : 
  AppLattice(spk,narg,arg)
{
  ninteger = 2;
  ndouble = 0;
  delpropensity = 1;
  delevent = 1;
  allow_kmc = 1;
  allow_rejection = 0;
  allow_masking = 0;
  

  create_arrays();

  if (narg != 1) error->all(FLERR,"Illegal app_style command");

  cycle = 0;
  pressureOn = 1;
  hello = 1;
  firsttime = 1;
  esites = NULL;
  echeck = NULL;
  events = NULL;
  maxevent = 0;
  firstevent = NULL;

  // reaction lists

  none = ntwo = nthree = nfour = 0;
  srate = drate = vrate = frate = NULL;
  spropensity = dpropensity = vpropensity = fpropensity =NULL;
  sinput = soutput = NULL;
  dinput = doutput = NULL;
  vinput = voutput = NULL;
  finput = foutput = NULL;
  comneigh = NULL;
  scount = dcount = vcount = fcount = NULL;
  sA = dA = vA = fA = NULL;
  scoord = dcoord = vcoord = fcoord = NULL;
  sexpon = dexpon = vexpon = fexpon = NULL;
  spresson = dpresson = vpresson = fpresson = NULL;
}

/* ---------------------------------------------------------------------- */

AppAldTMA::~AppAldTMA()
{
  delete [] esites;
  delete [] echeck;
  memory->sfree(events);
  memory->sfree(firstevent);
  memory->sfree(srate);
  memory->sfree(drate);
  memory->sfree(vrate);
  memory->sfree(frate);
  memory->sfree(spropensity);
  memory->sfree(dpropensity);
  memory->sfree(vpropensity);
  memory->sfree(fpropensity);
  memory->sfree(sinput);
  memory->sfree(soutput);
  memory->sfree(dinput);
  memory->sfree(doutput);
  memory->sfree(vinput);
  memory->sfree(voutput);
  memory->sfree(finput);
  memory->sfree(foutput);
  memory->sfree(comneigh);
  memory->sfree(scount);
  memory->sfree(dcount);
  memory->sfree(vcount);
  memory->sfree(fcount);
  memory->sfree(sA);
  memory->sfree(dA);
  memory->sfree(vA);
  memory->sfree(fA);
  memory->sfree(scoord);
  memory->sfree(dcoord);
  memory->sfree(vcoord);
  memory->sfree(fcoord);
  memory->sfree(sexpon);
  memory->sfree(dexpon);
  memory->sfree(vexpon);
  memory->sfree(fexpon);
  memory->sfree(spresson);
  memory->sfree(dpresson);
  memory->sfree(vpresson);
  memory->sfree(fpresson);

}

/* ---------------------------------------------------------------------- */

void AppAldTMA::input_app(char *command, int narg, char **arg)
{
  if (strcmp(command,"event") == 0) {
    if (narg < 1) error->all(FLERR,"Illegal event command");
    int rstyle = atoi(arg[0]);
    grow_reactions(rstyle);

    if (rstyle == 1) {
      if (narg != 9) error->all(FLERR,"Illegal event arg command");
//type I
      if (strcmp(arg[1],"O") == 0) sinput[none] = O;
      else if (strcmp(arg[1],"OAlX2") == 0) sinput[none] = OAlX2;
      else if (strcmp(arg[1],"OAlX") == 0) sinput[none] = OAlX;
      else if (strcmp(arg[1],"OAlXOH") == 0) sinput[none] = OAlXOH;
      else if (strcmp(arg[1],"OAlOH2") == 0) sinput[none] = OAlOH2;
      else if (strcmp(arg[1],"OAlOH") == 0) sinput[none] = OAlOH;
      else if (strcmp(arg[1],"OHAlX3") == 0) sinput[none] = OHAlX3;
      else if (strcmp(arg[1],"OAlX2H2O") == 0) sinput[none] = OAlX2H2O;
      else if (strcmp(arg[1],"OAlXH2O") == 0) sinput[none] = OAlXH2O;
      else if (strcmp(arg[1],"AlXOHH2O") == 0) sinput[none] = AlXOHH2O;
      else if (strcmp(arg[1],"OH") == 0) sinput[none] = OH;
      else if (strcmp(arg[1],"OAl") == 0) sinput[none] = OAl;
      else if (strcmp(arg[1],"Al") == 0) sinput[none] = Al;
      else if (strcmp(arg[1],"VAC") == 0) sinput[none] = VACANCY;
      else if (strcmp(arg[1],"AlX") == 0) sinput[none] = AlX;
      else if (strcmp(arg[1],"AlXOH") == 0) sinput[none] = AlXOH;
      else if (strcmp(arg[1],"AlOH2") == 0) sinput[none] = AlOH2;
      else if (strcmp(arg[1],"AlOH") == 0) sinput[none] = AlOH;

      
      else error->all(FLERR,"Illegal event arg1 command");

      if (strcmp(arg[2],"O") == 0) soutput[none] = O;
      else if (strcmp(arg[2],"OAlX2") == 0) soutput[none] = OAlX2;
      else if (strcmp(arg[2],"OAlX") == 0) soutput[none] = OAlX;
      else if (strcmp(arg[2],"OAlXOH") == 0) soutput[none] = OAlXOH;
      else if (strcmp(arg[2],"OAlOH2") == 0) soutput[none] = OAlOH2;
      else if (strcmp(arg[2],"OAlOH") == 0) soutput[none] = OAlOH;
      else if (strcmp(arg[2],"OHAlX3") == 0) soutput[none] = OHAlX3;
      else if (strcmp(arg[2],"OAlX2H2O") == 0) soutput[none] = OAlX2H2O;
      else if (strcmp(arg[2],"OAlXH2O") == 0) soutput[none] = OAlXH2O;
      else if (strcmp(arg[2],"AlXOHH2O") == 0) soutput[none] = AlXOHH2O;
      else if (strcmp(arg[2],"OH") == 0) soutput[none] = OH;
      else if (strcmp(arg[2],"OAl") == 0) soutput[none] = OAl;
      else if (strcmp(arg[2],"Al") == 0) soutput[none] = Al;
      else if (strcmp(arg[2],"VAC") == 0) soutput[none] = VACANCY;
      else if (strcmp(arg[2],"AlX") == 0) soutput[none] = AlX;
      else if (strcmp(arg[2],"AlXOH") == 0) soutput[none] = AlXOH;
      else if (strcmp(arg[2],"AlOH2") == 0) soutput[none] = AlOH2;
      else if (strcmp(arg[2],"AlOH") == 0) soutput[none] = AlOH;

      else error->all(FLERR,"Illegal event command");
      
      sA[none] = atof(arg[3]);
      if (sA[none] == 0.0) error->warning(FLERR,"Illegal coef during reading command");
      sexpon[none] = atoi(arg[4]);
      srate[none] = atof(arg[5]);
      scoord[none] = atoi(arg[6]);
      spresson[none] = atoi(arg[7]);

      none++;
      
//type II 
    } else if (rstyle == 2) {
      if (narg != 11) error->all(FLERR,"Illegal event command");

      if (strcmp(arg[1],"O") == 0) dinput[ntwo][0] = O;
      else if (strcmp(arg[1],"OAlX2") == 0) dinput[ntwo][0]= OAlX2;
      else if (strcmp(arg[1],"OAlX") == 0) dinput[ntwo][0]= OAlX;
      else if (strcmp(arg[1],"OAlXOH") == 0) dinput[ntwo][0]= OAlXOH;
      else if (strcmp(arg[1],"OAlOH2") == 0) dinput[ntwo][0]= OAlOH2;
      else if (strcmp(arg[1],"OAlOH") == 0) dinput[ntwo][0]= OAlOH;
      else if (strcmp(arg[1],"OHAlX3") == 0) dinput[ntwo][0]= OHAlX3;
      else if (strcmp(arg[1],"OAlX2H2O") == 0) dinput[ntwo][0]= OAlX2H2O;
      else if (strcmp(arg[1],"OAlXH2O") == 0) dinput[ntwo][0]= OAlXH2O;
      else if (strcmp(arg[1],"AlXOHH2O") == 0) dinput[ntwo][0]= AlXOHH2O;
      else if (strcmp(arg[1],"OH") == 0) dinput[ntwo][0]= OH;
      else if (strcmp(arg[1],"OAl") == 0) dinput[ntwo][0]= OAl;
      else if (strcmp(arg[1],"Al") == 0) dinput[ntwo][0]= Al;
      else if (strcmp(arg[1],"VAC") == 0) dinput[ntwo][0]= VACANCY;
      else if (strcmp(arg[1],"AlX") == 0) dinput[ntwo][0]= AlX;
      else if (strcmp(arg[1],"AlXOH") == 0) dinput[ntwo][0]= AlXOH;
      else if (strcmp(arg[1],"AlOH") == 0) dinput[ntwo][0]= AlOH;
      else if (strcmp(arg[1],"AlOH2") == 0) dinput[ntwo][0]= AlOH2;

      else error->all(FLERR,"Illegal event command");
      
      if (strcmp(arg[2],"O") == 0) doutput[ntwo][0] = O;
      else if (strcmp(arg[2],"OAlX2") == 0) doutput[ntwo][0]= OAlX2;
      else if (strcmp(arg[2],"OAlX") == 0) doutput[ntwo][0]= OAlX;
      else if (strcmp(arg[2],"OAlXOH") == 0) doutput[ntwo][0]= OAlXOH;
      else if (strcmp(arg[2],"OAlOH2") == 0) doutput[ntwo][0]= OAlOH2;
      else if (strcmp(arg[2],"OAlOH") == 0) doutput[ntwo][0]= OAlOH;
      else if (strcmp(arg[2],"OHAlX3") == 0) doutput[ntwo][0]= OHAlX3;
      else if (strcmp(arg[2],"OAlX2H2O") == 0) doutput[ntwo][0]= OAlX2H2O;
      else if (strcmp(arg[2],"OAlXH2O") == 0) doutput[ntwo][0]= OAlXH2O;
      else if (strcmp(arg[2],"AlXOHH2O") == 0) doutput[ntwo][0]= AlXOHH2O;
      else if (strcmp(arg[2],"OH") == 0) doutput[ntwo][0]= OH;
      else if (strcmp(arg[2],"OAl") == 0) doutput[ntwo][0]= OAl;
      else if (strcmp(arg[2],"Al") == 0) doutput[ntwo][0]= Al;
      else if (strcmp(arg[2],"VAC") == 0) doutput[ntwo][0]= VACANCY;
      else if (strcmp(arg[2],"AlX") == 0) doutput[ntwo][0]= AlX;
      else if (strcmp(arg[2],"AlXOH") == 0) doutput[ntwo][0]= AlXOH;
      else if (strcmp(arg[2],"AlOH") == 0) doutput[ntwo][0]= AlOH;
      else if (strcmp(arg[2],"AlOH2") == 0) doutput[ntwo][0]= AlOH2;

	
      else error->all(FLERR,"Illegal event command2");

      if (strcmp(arg[3],"O") == 0) dinput[ntwo][1] = O;
      else if (strcmp(arg[3],"OAlX2") == 0) dinput[ntwo][1]= OAlX2;
      else if (strcmp(arg[3],"OAlX") == 0) dinput[ntwo][1]= OAlX;
      else if (strcmp(arg[3],"OAlXOH") == 0) dinput[ntwo][1]= OAlXOH;
      else if (strcmp(arg[3],"OAlOH2") == 0) dinput[ntwo][1]= OAlOH2;
      else if (strcmp(arg[3],"OAlOH") == 0) dinput[ntwo][1]= OAlOH;
      else if (strcmp(arg[3],"OHAlX3") == 0) dinput[ntwo][1]= OHAlX3;
      else if (strcmp(arg[3],"OAlX2H2O") == 0) dinput[ntwo][1]= OAlX2H2O;
      else if (strcmp(arg[3],"OAlXH2O") == 0) dinput[ntwo][1]= OAlXH2O;
      else if (strcmp(arg[3],"AlXOHH2O") == 0) dinput[ntwo][1]= AlXOHH2O;
      else if (strcmp(arg[3],"OH") == 0) dinput[ntwo][1]= OH;
      else if (strcmp(arg[3],"OAl") == 0) dinput[ntwo][1]= OAl;
      else if (strcmp(arg[3],"Al") == 0) dinput[ntwo][1]= Al;
      else if (strcmp(arg[3],"VAC") == 0) dinput[ntwo][1]= VACANCY;
      else if (strcmp(arg[3],"AlX") == 0) dinput[ntwo][1]= AlX;
      else if (strcmp(arg[3],"AlXOH") == 0) dinput[ntwo][1]= AlXOH;
      else if (strcmp(arg[3],"AlOH") == 0) dinput[ntwo][1]= AlOH;
      else if (strcmp(arg[3],"AlOH2") == 0) dinput[ntwo][1]= AlOH2;

      else error->all(FLERR,"Illegal event command2");

      if (strcmp(arg[4],"O") == 0) doutput[ntwo][1] = O;
      else if (strcmp(arg[4],"OAlX2") == 0) doutput[ntwo][1]= OAlX2;
      else if (strcmp(arg[4],"OAlX") == 0) doutput[ntwo][1]= OAlX;
      else if (strcmp(arg[4],"OAlXOH") == 0) doutput[ntwo][1]= OAlXOH;
      else if (strcmp(arg[4],"OAlOH2") == 0) doutput[ntwo][1]= OAlOH2;
      else if (strcmp(arg[4],"OAlOH") == 0) doutput[ntwo][1]= OAlOH;
      else if (strcmp(arg[4],"OAlX4") == 0) doutput[ntwo][1]= OHAlX3;
      else if (strcmp(arg[4],"OAlX2H2O") == 0) doutput[ntwo][1]= OAlX2H2O;
      else if (strcmp(arg[4],"OAlXH2O") == 0) doutput[ntwo][1]= OAlXH2O;
      else if (strcmp(arg[4],"AlXOHH2O") == 0) doutput[ntwo][1]= AlXOHH2O;
      else if (strcmp(arg[4],"OH") == 0) doutput[ntwo][1]= OH;
      else if (strcmp(arg[4],"OAl") == 0) doutput[ntwo][1]= OAl;
      else if (strcmp(arg[4],"Al") == 0) doutput[ntwo][1]= Al;
      else if (strcmp(arg[4],"VAC") == 0) doutput[ntwo][1]= VACANCY;
      else if (strcmp(arg[4],"AlX") == 0) doutput[ntwo][1]= AlX;
      else if (strcmp(arg[4],"AlXOH") == 0) doutput[ntwo][1]= AlXOH;
      else if (strcmp(arg[4],"AlOH") == 0) doutput[ntwo][1]= AlOH;
      else if (strcmp(arg[4],"AlOH2") == 0) doutput[ntwo][1]= AlOH2;


      else error->all(FLERR,"Illegal event command2");

      dA[ntwo] = atof(arg[5]);
      dexpon[ntwo] = atoi(arg[6]);
      if (dexpon[ntwo] != 0.0) error->warning(FLERR,"Illegal expon command2");
      drate[ntwo] = atof(arg[7]);
      dcoord[ntwo] = atoi(arg[8]);
      dpresson[ntwo] = atoi(arg[9]);
      ntwo++;
// type III
    }else if (rstyle == 3) {
      if (narg != 11) error->all(FLERR,"Illegal event command31");
 
      if (strcmp(arg[1],"O") == 0) vinput[nthree][0]= O;
      else if (strcmp(arg[1],"OAlX2") == 0) vinput[nthree][0]= OAlX2;
      else if (strcmp(arg[1],"OAlX") == 0) vinput[nthree][0]= OAlX;
      else if (strcmp(arg[1],"OAlXOH") == 0) vinput[nthree][0]= OAlXOH;
      else if (strcmp(arg[1],"OAlOH2") == 0) vinput[nthree][0]= OAlOH2;
      else if (strcmp(arg[1],"OAlOH") == 0) vinput[nthree][0]= OAlOH;
      else if (strcmp(arg[1],"OHAlX3") == 0) vinput[nthree][0]= OHAlX3;
      else if (strcmp(arg[1],"OAlX2H2O") == 0) vinput[nthree][0]= OAlX2H2O;
      else if (strcmp(arg[1],"OAlXH2O") == 0) vinput[nthree][0]= OAlXH2O;
      else if (strcmp(arg[1],"AlXOHH2O") == 0) vinput[nthree][0]= AlXOHH2O;
      else if (strcmp(arg[1],"OH") == 0) vinput[nthree][0]= OH;
      else if (strcmp(arg[1],"OAl") == 0) vinput[nthree][0]= OAl;
      else if (strcmp(arg[1],"Al") == 0) vinput[nthree][0]= Al;
      else if (strcmp(arg[1],"VAC") == 0) vinput[nthree][0]= VACANCY;
      else if (strcmp(arg[1],"AlX") == 0) vinput[nthree][0]= AlX;
      else if (strcmp(arg[1],"AlXOH") == 0) vinput[nthree][0]= AlXOH;
      else if (strcmp(arg[1],"AlOH") == 0) vinput[nthree][0]= AlOH;
      else if (strcmp(arg[1],"AlOH2") == 0) vinput[nthree][0]= AlOH2;

      else error->all(FLERR,"Illegal event command32");


      if (strcmp(arg[2],"O") == 0) voutput[nthree][0]= O;
      else if (strcmp(arg[2],"OAlX2") == 0) voutput[nthree][0]= OAlX2;
      else if (strcmp(arg[2],"OAlX") == 0) voutput[nthree][0]= OAlX;
      else if (strcmp(arg[2],"OAlXOH") == 0) voutput[nthree][0]= OAlXOH;
      else if (strcmp(arg[2],"OAlOH2") == 0) voutput[nthree][0]= OAlOH2;
      else if (strcmp(arg[2],"OAlOH") == 0) voutput[nthree][0]= OAlOH;
      else if (strcmp(arg[2],"OHAlX3") == 0) voutput[nthree][0]= OHAlX3;
      else if (strcmp(arg[2],"OAlX2H2O") == 0) voutput[nthree][0]= OAlX2H2O;
      else if (strcmp(arg[2],"OAlXH2O") == 0) voutput[nthree][0]= OAlXH2O;
      else if (strcmp(arg[2],"AlXOHH2O") == 0) voutput[nthree][0]= AlXOHH2O;
      else if (strcmp(arg[2],"OH") == 0) voutput[nthree][0]= OH;
      else if (strcmp(arg[2],"OAl") == 0) voutput[nthree][0]= OAl;
      else if (strcmp(arg[2],"Al") == 0) voutput[nthree][0]= Al;
      else if (strcmp(arg[2],"VAC") == 0) voutput[nthree][0]= VACANCY;
      else if (strcmp(arg[2],"AlX") == 0) voutput[nthree][0]= AlX;
      else if (strcmp(arg[2],"AlXOH") == 0) voutput[nthree][0]= AlXOH;
      else if (strcmp(arg[2],"AlOH") == 0) voutput[nthree][0]= AlOH;
      else if (strcmp(arg[2],"AlOH2") == 0) voutput[nthree][0]= AlOH2;

      else error->all(FLERR,"Illegal event command33");
      
      if (strcmp(arg[3],"O") == 0) vinput[nthree][1]= O;
      else if (strcmp(arg[3],"OAlX2") == 0) vinput[nthree][1]= OAlX2;
      else if (strcmp(arg[3],"OAlX") == 0) vinput[nthree][1]= OAlX;
      else if (strcmp(arg[3],"OAlXOH") == 0) vinput[nthree][1]= OAlXOH;
      else if (strcmp(arg[3],"OAlOH2") == 0) vinput[nthree][1]= OAlOH2;
      else if (strcmp(arg[3],"OAlOH") == 0) vinput[nthree][1]= OAlOH;
      else if (strcmp(arg[3],"OHAlX3") == 0) vinput[nthree][1]= OHAlX3;
      else if (strcmp(arg[3],"OAlX2H2O") == 0) vinput[nthree][1]= OAlX2H2O;
      else if (strcmp(arg[3],"OAlXH2O") == 0) vinput[nthree][1]= OAlXH2O;
      else if (strcmp(arg[3],"AlXOHH2O") == 0) vinput[nthree][1]= AlXOHH2O;
      else if (strcmp(arg[3],"OH") == 0) vinput[nthree][1]= OH;
      else if (strcmp(arg[3],"OAl") == 0) vinput[nthree][1]= OAl;
      else if (strcmp(arg[3],"Al") == 0) vinput[nthree][1]= Al;
      else if (strcmp(arg[3],"VAC") == 0) vinput[nthree][1]= VACANCY;
      else if (strcmp(arg[3],"AlX") == 0) vinput[nthree][1]= AlX;
      else if (strcmp(arg[3],"AlXOH") == 0) vinput[nthree][1]= AlXOH;
      else if (strcmp(arg[3],"AlOH") == 0) vinput[nthree][1]= AlOH;
      else if (strcmp(arg[3],"AlOH2") == 0) vinput[nthree][1]= AlOH2;

      else error->all(FLERR,"Illegal event command34");
      
      if (strcmp(arg[4],"O") == 0) voutput[nthree][1]= O;
      else if (strcmp(arg[4],"OAlX2") == 0) voutput[nthree][1]= OAlX2;
      else if (strcmp(arg[4],"OAlX") == 0) voutput[nthree][1]= OAlX;
      else if (strcmp(arg[4],"OAlXOH") == 0) voutput[nthree][1]= OAlXOH;
      else if (strcmp(arg[4],"OAlOH2") == 0) voutput[nthree][1]= OAlOH2;
      else if (strcmp(arg[4],"OAlOH") == 0) voutput[nthree][1]= OAlOH;
      else if (strcmp(arg[4],"OAlX4") == 0) voutput[nthree][1]= OHAlX3;
      else if (strcmp(arg[4],"OAlX2H2O") == 0) voutput[nthree][1]= OAlX2H2O;
      else if (strcmp(arg[4],"OAlXH2O") == 0) voutput[nthree][1]= OAlXH2O;
      else if (strcmp(arg[4],"AlXOHH2O") == 0) voutput[nthree][1]= AlXOHH2O;
      else if (strcmp(arg[4],"OH") == 0) voutput[nthree][1]= OH;
      else if (strcmp(arg[4],"OAl") == 0) voutput[nthree][1]= OAl;
      else if (strcmp(arg[4],"Al") == 0) voutput[nthree][1]= Al;
      else if (strcmp(arg[4],"VAC") == 0) voutput[nthree][1]= VACANCY;
      else if (strcmp(arg[4],"AlX") == 0) voutput[nthree][1]= AlX;
      else if (strcmp(arg[4],"AlXOH") == 0) voutput[nthree][1]= AlXOH;
      else if (strcmp(arg[4],"AlOH") == 0) voutput[nthree][1]= AlOH;
      else if (strcmp(arg[4],"AlOH2") == 0) voutput[nthree][1]= AlOH2;

      else error->all(FLERR,"Illegal event command35");

      vA[nthree] = atof(arg[5]);
      vexpon[nthree] = atoi(arg[6]);
      if (vexpon[nthree] != 0.0) error->warning(FLERR,"Illegal vexpon command36");
      vrate[nthree] = atof(arg[7]);
      vcoord[nthree] = atoi(arg[8]);
      vpresson[nthree] = atoi(arg[9]);
      nthree++;

//type 4
    }else if (rstyle == 4) {
      if (narg != 11) error->all(FLERR,"Illegal event command37");
 
      if (strcmp(arg[1],"O") == 0) finput[nfour][0]= O;
      else if (strcmp(arg[1],"OAlX2") == 0) finput[nfour][0]= OAlX2;
      else if (strcmp(arg[1],"OAlX") == 0) finput[nfour][0]= OAlX;
      else if (strcmp(arg[1],"OAlXOH") == 0) finput[nfour][0]= OAlXOH;
      else if (strcmp(arg[1],"OAlOH2") == 0) finput[nfour][0]= OAlOH2;
      else if (strcmp(arg[1],"OAlOH") == 0) finput[nfour][0]= OAlOH;
      else if (strcmp(arg[1],"OHAlX3") == 0) finput[nfour][0]= OHAlX3;
      else if (strcmp(arg[1],"OAlX2H2O") == 0) finput[nfour][0]= OAlX2H2O;
      else if (strcmp(arg[1],"OAlXH2O") == 0) finput[nfour][0]= OAlXH2O;
      else if (strcmp(arg[1],"AlXOHH2O") == 0) finput[nfour][0]= AlXOHH2O;
      else if (strcmp(arg[1],"OH") == 0) finput[nfour][0]= OH;
      else if (strcmp(arg[1],"OAl") == 0) finput[nfour][0]= OAl;
      else if (strcmp(arg[1],"Al") == 0) finput[nfour][0]= Al;
      else if (strcmp(arg[1],"VAC") == 0) finput[nfour][0]= VACANCY;
      else if (strcmp(arg[1],"AlX") == 0) finput[nfour][0]= AlX;
      else if (strcmp(arg[1],"AlXOH") == 0) finput[nfour][0]= AlXOH;
      else if (strcmp(arg[1],"AlOH") == 0) finput[nfour][0]= AlOH;
      else if (strcmp(arg[1],"AlOH2") == 0) finput[nfour][0]= AlOH2;

      else error->all(FLERR,"Illegal event command38");


      if (strcmp(arg[2],"O") == 0) foutput[nfour][0]= O;
      else if (strcmp(arg[2],"OAlX2") == 0) foutput[nfour][0]= OAlX2;
      else if (strcmp(arg[2],"OAlX") == 0) foutput[nfour][0]= OAlX;
      else if (strcmp(arg[2],"OAlXOH") == 0) foutput[nfour][0]= OAlXOH;
      else if (strcmp(arg[2],"OAlOH2") == 0) foutput[nfour][0]= OAlOH2;
      else if (strcmp(arg[2],"OAlOH") == 0) foutput[nfour][0]= OAlOH;
      else if (strcmp(arg[2],"OHAlX3") == 0) foutput[nfour][0]= OHAlX3;
      else if (strcmp(arg[2],"OAlX2H2O") == 0) foutput[nfour][0]= OAlX2H2O;
      else if (strcmp(arg[2],"OAlXH2O") == 0) foutput[nfour][0]= OAlXH2O;
      else if (strcmp(arg[2],"AlXOHH2O") == 0) foutput[nfour][0]= AlXOHH2O;
      else if (strcmp(arg[2],"OH") == 0) foutput[nfour][0]= OH;
      else if (strcmp(arg[2],"OAl") == 0) foutput[nfour][0]= OAl;
      else if (strcmp(arg[2],"Al") == 0) foutput[nfour][0]= Al;
      else if (strcmp(arg[2],"VAC") == 0) foutput[nfour][0]= VACANCY;
      else if (strcmp(arg[2],"AlX") == 0) foutput[nfour][0]= AlX;
      else if (strcmp(arg[2],"AlXOH") == 0) foutput[nfour][0]= AlXOH;
      else if (strcmp(arg[2],"AlOH") == 0) foutput[nfour][0]= AlOH;
      else if (strcmp(arg[2],"AlOH2") == 0) foutput[nfour][0]= AlOH2;

      else error->all(FLERR,"Illegal event command39");
      
      if (strcmp(arg[3],"O") == 0) finput[nfour][1]= O;
      else if (strcmp(arg[3],"OAlX2") == 0) finput[nfour][1]= OAlX2;
      else if (strcmp(arg[3],"OAlX") == 0) finput[nfour][1]= OAlX;
      else if (strcmp(arg[3],"OAlXOH") == 0) finput[nfour][1]= OAlXOH;
      else if (strcmp(arg[3],"OAlOH2") == 0) finput[nfour][1]= OAlOH2;
      else if (strcmp(arg[3],"OAlOH") == 0) finput[nfour][1]= OAlOH;
      else if (strcmp(arg[3],"OHAlX3") == 0) finput[nfour][1]= OHAlX3;
      else if (strcmp(arg[3],"OAlX2H2O") == 0) finput[nfour][1]= OAlX2H2O;
      else if (strcmp(arg[3],"OAlXH2O") == 0) finput[nfour][1]= OAlXH2O;
      else if (strcmp(arg[3],"AlXOHH2O") == 0) finput[nfour][1]= AlXOHH2O;
      else if (strcmp(arg[3],"OH") == 0) finput[nfour][1]= OH;
      else if (strcmp(arg[3],"OAl") == 0) finput[nfour][1]= OAl;
      else if (strcmp(arg[3],"Al") == 0) finput[nfour][1]= Al;
      else if (strcmp(arg[3],"VAC") == 0) finput[nfour][1]= VACANCY;
      else if (strcmp(arg[3],"AlX") == 0) finput[nfour][1]= AlX;
      else if (strcmp(arg[3],"AlXOH") == 0) finput[nfour][1]= AlXOH;
      else if (strcmp(arg[3],"AlOH") == 0) finput[nfour][1]= AlOH;
      else if (strcmp(arg[3],"AlOH2") == 0) finput[nfour][1]= AlOH2;

      else error->all(FLERR,"Illegal event command40");
      
      if (strcmp(arg[4],"O") == 0) foutput[nfour][1]= O;
      else if (strcmp(arg[4],"OAlX2") == 0) foutput[nfour][1]= OAlX2;
      else if (strcmp(arg[4],"OAlX") == 0) foutput[nfour][1]= OAlX;
      else if (strcmp(arg[4],"OAlXOH") == 0) foutput[nfour][1]= OAlXOH;
      else if (strcmp(arg[4],"OAlOH2") == 0) foutput[nfour][1]= OAlOH2;
      else if (strcmp(arg[4],"OAlOH") == 0) foutput[nfour][1]= OAlOH;
      else if (strcmp(arg[4],"OAlX4") == 0) foutput[nfour][1]= OHAlX3;
      else if (strcmp(arg[4],"OAlX2H2O") == 0) foutput[nfour][1]= OAlX2H2O;
      else if (strcmp(arg[4],"OAlXH2O") == 0) foutput[nfour][1]= OAlXH2O;
      else if (strcmp(arg[4],"AlXOHH2O") == 0) foutput[nfour][1]= AlXOHH2O;
      else if (strcmp(arg[4],"OH") == 0) foutput[nfour][1]= OH;
      else if (strcmp(arg[4],"OAl") == 0) foutput[nfour][1]= OAl;
      else if (strcmp(arg[4],"Al") == 0) foutput[nfour][1]= Al;
      else if (strcmp(arg[4],"VAC") == 0) foutput[nfour][1]= VACANCY;
      else if (strcmp(arg[4],"AlX") == 0) foutput[nfour][1]= AlX;
      else if (strcmp(arg[4],"AlXOH") == 0) foutput[nfour][1]= AlXOH;
      else if (strcmp(arg[4],"AlOH") == 0) foutput[nfour][1]= AlOH;
      else if (strcmp(arg[4],"AlOH2") == 0) foutput[nfour][1]= AlOH2;

      else error->all(FLERR,"Illegal event command41");

      fA[nfour] = atof(arg[5]);
      fexpon[nfour] = atoi(arg[6]);
      if (fexpon[nfour] != 0.0) error->warning(FLERR,"Illegal fexpon command42");
      frate[nfour] = atof(arg[7]);
      fcoord[nfour] = atoi(arg[8]);
      fpresson[nfour] = atoi(arg[9]);
      nfour++;
    } else error->all(FLERR,"Illegal event command43");
  } 
  else if (strcmp(command,"pulse_time") == 0) {
    if (narg != 2) error->all(FLERR,"Illegal pulse time");
      T1 = atof(arg[0]);
      T3 = atof(arg[1]);
  }
  else if (strcmp(command,"purge_time") == 0) {
    if (narg != 2) error->all(FLERR,"Illegal purge time");
      T2 = atof(arg[0]);
      T4 = atof(arg[1]);
  }else error->all(FLERR,"Unrecognized command38");
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppAldTMA::grow_app()
{
  element = iarray[0];
  coord = iarray[1];
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppAldTMA::init_app()
{
  if (firsttime) {
    firsttime = 0;

    echeck = new int[nlocal];
    firstevent = (int *) memory->smalloc(nlocal*sizeof(int),"app:firstevent");
    //comneigh was defined to avoid double counting of common neighbor in site_propensity
    comneigh = memory->grow(comneigh,12*maxneigh,2,"app/ald:comneigh");
    // esites must be large enough for 3 sites and their 1st neighbors
    
    esites = (int *) memory->smalloc(12*maxneigh*sizeof(int),"app:esites");
    //esites = new int[12*maxneigh]; 
  }
  // site validity

  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    if (coord[i] < -1 || coord[i] > 8) flag = 1;
    if (element[i] < VACANCY) flag = 1;
  }
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all(FLERR,"One or more sites have invalid values");
}
/* ---------------------------------------------------------------------- */

void AppAldTMA::setup_app()
{
  for (int i = 0; i < nlocal; i++) echeck[i] = 0;
  

  nevents = 0;
  for (int i = 0; i < nlocal; i++) firstevent[i] = -1;
  for (int i = 0; i < maxevent; i++) events[i].next = i+1;
  freeevent = 0;


  if (temperature == 0.0)
    error->all(FLERR,"Temperature cannot be 0.0 for app_ald");
  for (int m = 0; m < none; m++) {
    spropensity[m] = sA[m]*pow(temperature,sexpon[m])*exp(-srate[m]/temperature);
    scount[m] = 0;
  if (spropensity[m] == 0.0) error->warning(FLERR," spropensity cannot be 0.0 for app_ald");
  }
  for (int m = 0; m < ntwo; m++) {
    dpropensity[m] = dA[m]*pow(temperature,dexpon[m])*exp(-drate[m]/temperature);
    dcount[m] = 0;
  if (dpropensity[m] == 0.0) error->warning(FLERR,"dpropensity cannot be 0.0 for app_ald");
  }
  for (int m = 0; m < nthree; m++) {
    vpropensity[m] = vA[m]*pow(temperature,vexpon[m])*exp(-vrate[m]/temperature);
    vcount[m] = 0;
  if (vpropensity[m] == 0.0) error->warning(FLERR,"vpropensity cannot be 0.0 for app_ald");
  }
  for (int m = 0; m < nfour; m++) {
    fpropensity[m] = fA[m]*pow(temperature,fexpon[m])*exp(-frate[m]/temperature);
    fcount[m] = 0;
  if (fpropensity[m] == 0.0) error->warning(FLERR,"fpropensity cannot be 0.0 for app_ald");
  }
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppAldTMA::site_energy(int i)
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   KMC method
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double AppAldTMA::site_propensity(int i)
{
  int j,k,m,l,o,p;


   clear_events(i);
  // std::cout<<"clear_events"<<std::endl;
  double proball = 0.0;


  //type I, check species of sites and consider possible events

  // count_coordO was added here to prevent adsorption of metal precursor
  // on the low coordinate oxygen at sublayer

  for (m = 0; m < none; m++) {
    if (element[i] != sinput[m]) continue;
    int coordi = coord[i] < 0 ? (coord[i] % 10 + 10) : (coord[i] % 10);
    if ((coord[i] == scoord[m] || scoord[m] == -1) && (spresson[m] == pressureOn || spresson[m] == 0) && (coordi < numneigh[i])) {
      // std::cout<<"coord[i]"<<coord[i]<<std::endl;
      // std::cout<<"numneigh[i]"<<numneigh[i]<<std::endl;
	    add_event(i,1,m,spropensity[m],-1,-1,-1);
	    proball += spropensity[m];
    }
    //else if ((element[i]=TaX4O) && (spresson[m] == pressureOn || spresson[m] == 0) && (coordi <= numneigh[i])) {
    //  add_event(i,1,m,spropensity[m],-1,-1,-1);
	  //  proball += spropensity[m];
    //  }
  }
  // std::cout<<"test1"<<std::endl;
  // type II, check species of sites and second neighbor 
  // comneigh variable used to avoid double counting,
  // we consider more than one event between sites, therefore we need 2d comneigh array.

  int nextneib = 1;
  for (int jj = 0; jj < numneigh[i]; jj++) {
    j = neighbor[i][jj];
       for (int kk = 0; kk < numneigh[j]; kk++) {
            k = neighbor[j][kk];
            if (i == k) continue;
            for (m = 0; m < ntwo; m++) {
              int coordi = coord[i] < 0 ? (coord[i] % 10 + 10) : (coord[i] % 10);
              int coordk = coord[k] < 0 ? (coord[k] % 10 + 10) : (coord[k] % 10);
	    if ((element[i] == dinput[m][0] && element[k] == dinput[m][1]) && (dpresson[m] == pressureOn || dpresson[m] == 0) && (coord[i] == dcoord[m] || dcoord[m] == -1) ) {
  	      comevent = 1;
  	      for (int ii = 0; ii < nextneib; ii++) {
  	      if ( comneigh[ii][0] == k && comneigh[ii][1] == dpropensity[m]) comevent = 0;
  	      }
  	      if (comevent){	
                add_event(i,2,m,dpropensity[m],-1,k,-1);
                proball += dpropensity[m];
  	        comneigh[nextneib][0] = k;
                comneigh[nextneib][1] = dpropensity[m];
  	        nextneib++;
  	    } 
          }
        }
      }
    }
    for (m = 0; m < nextneib; m++) comneigh[m][0] = comneigh[m][1] = 0;
    // std::cout<<"test2"<<std::endl;
  //type III, check species of sites and first neighbour 
  for (int jj = 0; jj < numneigh[i]; jj++) {
    j = neighbor[i][jj];
      for (m =0; m < nthree; m++) {
      // std::cout<<"test3!!!"<<std::endl;
      int coordi = coord[i] < 0 ? (coord[i] % 10 + 10) : (coord[i] % 10);
      int coordj = coord[j] < 0 ? (coord[j] % 10 + 10) : (coord[j] % 10);
      if (element[i] == vinput[m][0] && element[j] == vinput[m][1] && (coord[i] == vcoord[m] || vcoord[m] == -1 ) && (vpresson[m] == pressureOn || vpresson[m] == 0) && (coordi < numneigh[i]) && (coordj < numneigh[j])) {
        // std::cout<<"coord[j]:"<<coord[j]<<std::endl;
        // std::cout<<"numneigh[j]:"<<numneigh[j]<<std::endl;
        add_event(i,3,m,vpropensity[m],j,-1,-1);
        proball += vpropensity[m];
      }
    }
  }
  // type IIII, check species of sites and third neighbor 
  for (int jj = 0; jj < numneigh[i]; jj++) {
    j = neighbor[i][jj];
      for (int kk = 0; kk < numneigh[j]; kk++) {
        k = neighbor[j][kk];
         for (int ll = 0; ll < numneigh[k]; ll++) {
        l = neighbor[k][ll];
         for (int oo = 0; oo < numneigh[l]; oo++) {
        o = neighbor[l][oo];
         for (int pp = 0; pp < numneigh[o]; pp++) {
        p = neighbor[o][pp];
          for (int gg = 0; gg < numneigh[p]; gg++) {
            int g = neighbor[p][gg];
            for (m = 0; m < nfour; m++) {
              int coordi = coord[i] < 0 ? (coord[i] % 10 + 10) : (coord[i] % 10);
              int coordg = coord[g] < 0 ? (coord[g] % 10 + 10) : (coord[g] % 10);
              if (element[i] == finput[m][0] && element[g] == finput[m][1] && (coord[i] == fcoord[m] || fcoord[m] == -1 ) && (fpresson[m] == pressureOn || fpresson[m] == 0) && (coordi < numneigh[i]) && (coordg < numneigh[g])) {
                add_event(i,4,m,fpropensity[m],-1,-1,g);
                proball += fpropensity[m];
              }
            }
        }
    }}}}
  }
  // vacant event:5
  add_event(i,5,-1,0.1,-1,-1,-1);
  proball += 0.1;
  // std::cout<<"test3"<<std::endl;
  return proball;
}

/* ----------------------------------------------------------------------
   KMC method
   choose and perform an event for site
------------------------------------------------------------------------- */

void AppAldTMA::site_event(int i, class RandomPark *random)
{
  int j,k,g,m,n,mm,jj;
  //int elcoord = element[i];
  
  double threshhold = random->uniform() * propensity[i2site[i]];
  double proball = 0.0;

  int ievent = firstevent[i];
  while (1) {
    proball += events[ievent].propensity;
    if (proball >= threshhold) break;
    ievent = events[ievent].next;
  }

  int rstyle = events[ievent].style;
  int which = events[ievent].which;
  j = events[ievent].jpartner;
  k = events[ievent].kpartner;
  g = events[ievent].gpartner;
  int elcoord_i = element[i];
  int elcoord_j = element[j];
  int elcoord_k = element[k];
  int elcoord_g = element[g];



  if (rstyle == 1) {
    element[i] = soutput[which];
    scount[which]++;
    } 
  else if (rstyle == 2 && j == -1) {
    element[i] = doutput[which][0];
    element[k] = doutput[which][1];
    dcount[which]++;
    }
  else if (rstyle == 3 && k == -1) {
    element[i] = voutput[which][0];
    element[j] = voutput[which][1];
    vcount[which]++;
    }
  else if (rstyle == 4) {
    element[i] = foutput[which][0];
    element[g] = foutput[which][1];
    fcount[which]++;
  }
  
  else if (rstyle == 5) {}
  else { error->all(FLERR,"Illegal execution event"); }

  //update_coord(elcoord,i,j,k,which);

  // sequence of ALD, 
  // 1 is metal pulse, 3 purge, 2 oxygen pulse.
  if ((cycle+T1) > time ) {pressureOn = 1;}
  else if ((cycle+T1)<= time && time < (cycle+T1+T2)) {pressureOn = 3;}
  else if ((cycle+T1+T2) <= time && time < (cycle+T1+T2+T3)) {pressureOn = 2;}
  else if ((cycle+T1+T2+T3) <= time && time < (cycle+T1+T2+T3+T4)) {pressureOn = 3;}
  else {cycle += T1+T2+T3+T4; }

  

  int nsites = 0;
  int isite = i2site[i];


// mask 
   if (rstyle == 1) {
    if (elcoord_i == OH && element[i] == OHAlX3) {
     put_mask(i);
    }
    else if (elcoord_i == OHAlX3 && element[i] == OH) {
      remove_mask(i);
    }
    else if (elcoord_i == OAlX2H2O && element[i] == OAlXOH) {
      remove_mask(i);
    }
  }


  else if (rstyle == 4) {
    if (elcoord_j == OAlX2 && element[j] == OAlX) {
      remove_mask(j);
    }
    }
 
  count_coord(i);
  propensity[isite] = site_propensity(i);
  esites[nsites++] = isite;
  echeck[isite] = 1;
for (int n = 0; n < numneigh[i]; n++) {
    m = neighbor[i][n];
    isite = i2site[m];
    if (isite >= 0 && echeck[isite] == 0) {
    count_coord(m);
    propensity[isite] = site_propensity(m);
    esites[nsites++] = isite;
    // std::cout<<"nsites:"<<nsites<<std::endl;
    echeck[isite] = 1;
    }
    for (int jj = 0; jj< numneigh[m];jj++) {
    mm = neighbor[m][jj];
    isite = i2site[mm];
        if (isite >= 0 && echeck[isite] == 0) {
            count_coord(mm);
            propensity[isite] = site_propensity(mm);
            esites[nsites++] = isite;
            // std::cout<<"nsites:"<<nsites<<std::endl;
            echeck[isite] = 1;
        }
        for (int ss = 0; ss< numneigh[mm]; ss++) {
            int s = neighbor[mm][ss];
            isite = i2site[s];
            if (isite >= 0 && echeck[isite] == 0) {
              count_coord(s);
              propensity[isite] = site_propensity(s);
              esites[nsites++] = isite;
              // std::cout<<"nsites:"<<nsites<<std::endl;
              echeck[isite] = 1;
            }
            for (int gg = 0; gg< numneigh[s]; gg++) {
              int g = neighbor[s][gg];
              isite = i2site[g];
              if (isite >= 0 && echeck[isite] == 0) {
                count_coord(g);
                propensity[isite] = site_propensity(g);
                esites[nsites++] = isite;
                // std::cout<<"nsites:"<<nsites<<std::endl;
                echeck[isite] = 1;
              }
            }
        }
    }
  }



  solve->update(nsites,esites,propensity);
   // clear echeck array

  for (m = 0; m < nsites; m++)  {echeck[esites[m]] = 0; esites[m]=0;}
  
}

/* ----------------------------------------------------------------------
   clear all events out of list for site I
   add cleared events to free list
------------------------------------------------------------------------- */

void AppAldTMA::clear_events(int i)
{
  int next;
  int index = firstevent[i];
  while (index >= 0) {
    next = events[index].next;
    events[index].next = freeevent;
    freeevent = index;
    nevents--;
    index = next;
  }
  firstevent[i] = -1;
}

/* ----------------------------------------------------------------------
   add an event to list for site I
   event = exchange with site J with probability = propensity
------------------------------------------------------------------------- */

void AppAldTMA::add_event(int i, int rstyle, int which, double propensity,
			  int jpartner, int kpartner, int gpartner)
{
  if (nevents == maxevent) {
    maxevent += DELTAEVENT;
    events = 
      (Event *) memory->srealloc(events,maxevent*sizeof(Event),"app:events");
    for (int m = nevents; m < maxevent; m++) events[m].next = m+1;
    freeevent = nevents;
  }

  int next = events[freeevent].next;

  events[freeevent].style = rstyle;
  events[freeevent].which = which;
  events[freeevent].jpartner = jpartner;
  events[freeevent].kpartner = kpartner;
  events[freeevent].gpartner = gpartner;
  events[freeevent].propensity = propensity;

  if ( propensity == 0 ) error->all(FLERR,"propensity in add_event wrong app ald");
  events[freeevent].next = firstevent[i];
  firstevent[i] = freeevent;
  freeevent = next;
  nevents++;
}

/* ----------------------------------------------------------------------
   grow list of stored reactions for single and double
------------------------------------------------------------------------- */

void AppAldTMA::grow_reactions(int rstyle)
{
  if (rstyle == 1) {
    int n = none + 1;
    srate = (double *) 
      memory->srealloc(srate,n*sizeof(double),"app/ald:srate");
    spropensity = (double *) 
      memory->srealloc(spropensity,n*sizeof(double),"app/ald:spropensity");
    sinput = (int *) 
      memory->srealloc(sinput,n*sizeof(int),"app/ald:sinput");
    soutput = (int *) 
      memory->srealloc(soutput,n*sizeof(int),"app/ald:soutput");
    scount = (int *) 
      memory->srealloc(scount,n*sizeof(int),"app/ald:scount");
    sA = (double *) 
      memory->srealloc(sA,n*sizeof(double),"app/ald:sA");
    sexpon = (int *) 
      memory->srealloc(sexpon,n*sizeof(int),"app/ald:sexpon");
    scoord = (int *) 
      memory->srealloc(scoord,n*sizeof(int),"app/ald:scoord");
    spresson = (int *) 
      memory->srealloc(spresson,n*sizeof(int),"app/ald:spresson");

  } else if (rstyle == 2) {
    int n = ntwo + 1;
    drate = (double *) 
      memory->srealloc(drate,n*sizeof(double),"app/ald:drate");
    dpropensity = (double *) 
      memory->srealloc(dpropensity,n*sizeof(double),"app/ald:dpropensity");
    dinput = memory->grow(dinput,n,2,"app/ald:dinput");
    doutput = memory->grow(doutput,n,2,"app/ald:doutput");
    dcount = (int *) 
      memory->srealloc(dcount,n*sizeof(int),"app/ald:dcount");
    dA = (double *) 
      memory->srealloc(dA,n*sizeof(double),"app/ald:dA");
    dexpon = (int *) 
      memory->srealloc(dexpon,n*sizeof(int),"app/ald:dexpon");
    dcoord = (int *) 
      memory->srealloc(dcoord,n*sizeof(int),"app/ald:dcoord");
    dpresson = (int *) 
      memory->srealloc(dpresson,n*sizeof(int),"app/ald:dpresson");

  } else if (rstyle == 3) {
    int n = nthree + 1;
    vrate = (double *)
      memory->srealloc(vrate,n*sizeof(double),"app/ald:vrate");
    vpropensity = (double *)
      memory->srealloc(vpropensity,n*sizeof(double),"app/ald:vpropensity");
    vinput = memory->grow(vinput,n,2,"app/ald:vinput");
    voutput = memory->grow(voutput,n,2,"app/ald:voutput");
    vcount = (int *)
      memory->srealloc(vcount,n*sizeof(int),"app/ald:vcount");
    vA = (double *)
      memory->srealloc(vA,n*sizeof(double),"app/ald:vA");
    vexpon = (int *)
      memory->srealloc(vexpon,n*sizeof(int),"app/ald:vexpon");
    vcoord = (int *)
      memory->srealloc(vcoord,n*sizeof(int),"app/ald:vcoord");
    vpresson = (int *)
      memory->srealloc(vpresson,n*sizeof(int),"app/ald:vpresson");
  }
  else if (rstyle == 4) {
    int n = nfour + 1;
    frate = (double *)
      memory->srealloc(frate,n*sizeof(double),"app/ald:frate");
    fpropensity = (double *)
      memory->srealloc(fpropensity,n*sizeof(double),"app/ald:fpropensity");
    finput = memory->grow(finput,n,2,"app/ald:finput");
    foutput = memory->grow(foutput,n,2,"app/ald:foutput");
    fcount = (int *)
      memory->srealloc(fcount,n*sizeof(int),"app/ald:fcount");
    fA = (double *)
      memory->srealloc(fA,n*sizeof(double),"app/ald:fA");
    fexpon = (int *)
      memory->srealloc(fexpon,n*sizeof(int),"app/ald:fexpon");
    fcoord = (int *)
      memory->srealloc(fcoord,n*sizeof(int),"app/ald:fcoord");
    fpresson = (int *)
      memory->srealloc(fpresson,n*sizeof(int),"app/ald:fpresson");
  }
}

/* ----------------------------------------------------------------------
   update c.n. for Hf and O, put and remove mask for relative sites
------------------------------------------------------------------------- */
//void AppAldTMA::update_coord(int elcoord, int i, int j, int k, int which)
//{ 
  //type1
  
//    if ((elcoord == O  && element[i] == TaX5O) && (j == -1)) 
//    {//Adsorption of TaX5
//      coord[i]++;
//      //put_mask(i);
//    }
//    else if ((elcoord == TaX5O && element[i] == O ) && (j == -1))
//    {//Desorption of TaX5
//      coord[i]--;
//      //remove_mask(i);
//      count_coordO(i);
//    }
//    else if ((elcoord == Ta  && element[i] == OTa ) && (j == -1))
//    {//Adsorption of oxygen
//      coord[i]++;
//    }
//    else if ((elcoord == OTa && element[i] == Ta ) && (j == -1))
//    {//Desorption of oxygen
//      coord[i]--;
//    }
//    else if ((elcoord == TaX  && element[i] == Ta ) && (j == -1))
//    {//Oxidation
//      coord[i]++;
//    }
//    
//    else if ((elcoord == TaX4O  && element[i] == TaO ) && (j == -1))
//    {//Oxidation
//      coord[i]++;
//   }
  //type2
  /*else if (j!=-1 && k==-1) {
    count_coord(i);
    count_coord(j);
  }*/
  //type3
  
    
//    else if ((elcoord == VACANCY) && (element[i] == TaO) && ( element[j] == O)) {
//     //densification
//        if (element[i] == TaO){
//       	//remove_mask(i, j); // Remove mask from previous O site
//      }
//        count_coord(j, i);
//        if (element[i] == TaO){
//            //put_mask(i); // Put mask on new Zn site
//        }}
//    else if ((elcoord == OTa) && (element[i] == Ta) && (element[j] == O) && ( k == -1))
//  {//oxygen densification
//    count_coord(i,j);
//    count_coordO(j);
//  }

//    else if ((elcoord == Ta) && (element[i] == VACANCY) && (element[j] == TaO))
//  {//TaX4 Reverse densification
    //remove_mask(i, j);
//		count_coord(i, j);
		//put_mask(j);
//  }
  
//    else if ((elcoord == O) && (element[i] == VACANCY) && (element[j] == OTa))
//  {//oxygen Reverse densification
//		count_coord(i,j);
//		coord[j]++;
//  }

//    else if ((elcoord ==TaX5O) && (element[i] == TaX4O) && (element[j] == TaX)) 
//  {//TaX5 Dissociation
  //remove_mask(i);
  //put_mask(i);
//  coord[j]++;
  //put_mask(j);
//  }
  
//}

/* ----------------------------------------------------------------------
   count c.n after densification
------------------------------------------------------------------------- */
void AppAldTMA::count_coord(int i) // i: Oxygen species(does not necessarily hold)
{
 if (0 <= coord[i]) coord[i] = 0;
  else if (-10 < coord[i] && coord[i] < 0) coord[i] = -10;
  else if (-20 < coord[i] && coord[i] < -10) coord[i] = -20;
  else if (-30 < coord[i] && coord[i] < -20) coord[i] = -30;
  else if (-40 < coord[i] && coord[i] < -30) coord[i] = -40;
  else if (-50 < coord[i] && coord[i] < -40) coord[i] = -50;

  for (int s = 0; s < numneigh[i]; s++) {
    int nn = neighbor[i][s];
    if (element[nn] != VACANCY) coord[i]++;
}
}

/* ---------------------------------------------------- ------------------
   count c.n of oxygen before adsorption
------------------------------------------------------------------------- */
void AppAldTMA::count_coordO(int i)
{
    int fullO = 0;
    int emptyO = 0;
    int totalS = 0;

    int isite = i2site[i];
    int nsites = 0;

	for (int m = 0; m < numneigh[i]; m++) {
		int mm = neighbor[i][m];
		for (int s = 0; s < numneigh[mm]; s++) {
			int ss = neighbor[mm][s];
			isite = i2site[ss];
			if (i==ss)  continue;
			if (isite >= 0 && echeck[isite] == 0) {
			  if ( element[ss] >= O && element[ss] <= OAlOH ) {fullO++;}
			  else if (element[ss] == VACANCY) {emptyO++;}
		          esites[nsites++] = isite;
		          echeck[isite] = 1;
		        }
    
		}
	}
   totalS = fullO+emptyO;
   if ( float(fullO) > 4*totalS/5 and coord[i] > -20) {coord[i] += -20;} // decrease the coord of the oxygen site to render it inactive for adsorption
   for (int m = 0; m < nsites; m++)  {echeck[esites[m]] = 0; esites[m]=0;}
}

/* ----------------------------------------------------------------------
   put mask for affected sites
------------------------------------------------------------------------- */
void AppAldTMA::put_mask(int i) {
  int nsites = 0;
  int isite = i2site[i];
  echeck[isite] = 1;
  esites[nsites++] = isite;
  for (int s = 0; s < numneigh[i]; s++) {
    int nn = neighbor[i][s];
    int isite = i2site[nn];
    // std::cout<<"isite:"<<isite<<std::endl;
    // std::cout<<"echeck:"<<echeck[isite]<<std::endl;
    if (isite >= 0 && echeck[isite] == 0) {
        echeck[isite] = 1;
        esites[nsites++] = isite;
        coord[nn] -= 0;
        // std::cout << "put mask " << i << std::endl;
    }
    for (int ss = 0; ss < numneigh[nn]; ss++){
      int nnn = neighbor[nn][ss];
      int isite = i2site[nnn];
      if (isite >= 0 && echeck[isite] == 0) {
        echeck[isite] = 1;
        esites[nsites++] = isite;
        coord[nnn] -= 10;
          }
      for (int sss = 0; sss < numneigh[nnn]; sss++){
        int nnnn = neighbor[nnn][sss];
        int isite = i2site[nnnn];
        if (isite >= 0 && echeck[isite] == 0) {
          echeck[isite] = 1;
          esites[nsites++] = isite;
          coord[nnnn] -= 0;
        }
        for (int ssss = 0; ssss < numneigh[nnnn]; ssss++){
          int nnnnn = neighbor[nnnn][ssss];
          int isite = i2site[nnnnn];
          if (isite >= 0 && echeck[isite] == 0) {
            echeck[isite] = 1;
            esites[nsites++] = isite;
            coord[nnnnn] -= 0;
        }
        }
  }
  }
  }
  // for (int s = 0; s < numneigh[i]; s++) {
  //   int nn = neighbor[i][s];
  //   coord[nn] -= 10;
  // }
  for (int m = 0; m < nsites; m++)  {echeck[esites[m]] = 0; esites[m]=0;}
  // std::cout << "put mask " << i << std::endl;
}

void AppAldTMA::put_mask_2(int i)  {
  int nsites = 0;
  int isite = i2site[i];
  echeck[isite] = 1;
  esites[nsites++] = isite;
  for (int s = 0; s < numneigh[i]; s++) {
    int nn = neighbor[i][s];
    int isite = i2site[nn];
    // std::cout<<"isite:"<<isite<<std::endl;
    // std::cout<<"echeck:"<<echeck[isite]<<std::endl;
    if (isite >= 0 && echeck[isite] == 0) {
        echeck[isite] = 1;
        esites[nsites++] = isite;
        coord[nn] -= 10;
    }
    for (int m = 0; m < nsites; m++)  {echeck[esites[m]] = 0; esites[m]=0;}
  }
}

void AppAldTMA::remove_mask(int i) {
  int nsites = 0;
  int isite = i2site[i];
  echeck[isite] = 1;
  esites[nsites++] = isite;
  for (int s = 0; s < numneigh[i]; s++) {
    int nn = neighbor[i][s];
    int isite = i2site[nn];
    // std::cout<<"isite:"<<isite<<std::endl;
    // std::cout<<"echeck:"<<echeck[isite]<<std::endl;
    if (isite >= 0 && echeck[isite] == 0) {
        echeck[isite] = 1;
        esites[nsites++] = isite;
        coord[nn] += 0;
        // std::cout << "put mask " << i << std::endl;
    }
    for (int ss = 0; ss < numneigh[nn]; ss++){
      int nnn = neighbor[nn][ss];
      int isite = i2site[nnn];
      if (isite >= 0 && echeck[isite] == 0) {
        echeck[isite] = 1;
        esites[nsites++] = isite;
        coord[nnn] += 10;
          }
      for (int sss = 0; sss < numneigh[nnn]; sss++){
        int nnnn = neighbor[nnn][sss];
        int isite = i2site[nnnn];
        if (isite >= 0 && echeck[isite] == 0) {
          echeck[isite] = 1;
          esites[nsites++] = isite;
          coord[nnnn] += 0;
        }
        for (int ssss = 0; ssss < numneigh[nnnn]; ssss++){
          int nnnnn = neighbor[nnnn][ssss];
          int isite = i2site[nnnnn];
          if (isite >= 0 && echeck[isite] == 0) {
            echeck[isite] = 1;
            esites[nsites++] = isite;
            coord[nnnnn] += 0;
        }
        }
  }
  }
  }
  // for (int s = 0; s < numneigh[i]; s++) {
  //   int nn = neighbor[i][s];
  //   coord[nn] += 10;
  // }
  for (int m = 0; m < nsites; m++)  {echeck[esites[m]] = 0; esites[m]=0;}
  // std::cout << "remove mask " << i << std::endl;
}

void AppAldTMA::remove_mask_2(int i) {
  int nsites = 0;
  int isite = i2site[i];
  echeck[isite] = 1;
  esites[nsites++] = isite;
  for (int s = 0; s < numneigh[i]; s++) {
    int nn = neighbor[i][s];
    int isite = i2site[nn];
    if (isite >= 0 && echeck[isite] == 0) {
        echeck[isite] = 1;
        esites[nsites++] = isite;
        coord[nn] += 10;
    }
  }
  for (int m = 0; m < nsites; m++)  {echeck[esites[m]] = 0; esites[m]=0;}
}