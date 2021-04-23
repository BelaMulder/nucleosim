/// MICROTUBULE METHODS
//#include "stdafx.h"
#include "defs.h"
#include "pars.h"


Microtubule::Microtubule(void)
{
   length = inrand()*VG*DT;
    state = GROWING;
    waitingtime =0.0;
}

Microtubule::Microtubule(enum ORIENTATION ornt, Barrier* trgt)
{
  length = inrand()*VG*DT;
   force = 0.0;
   waitingtime = 0.0;
   state = GROWING;
   orient = ornt;
   target = trgt;
}


void Microtubule::grow(double dlg)
       {
           length += dlg;
       }

void Microtubule::shrink(double dls)
       {
           length -= dls;
           if(length < 0.0){state = DORMANT;};
       }

//void Microtubule::wait(double wts)
//{
  //  wait +=wts;
//}
