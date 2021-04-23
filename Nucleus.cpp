/// BARRIER AND NUCLEUS METHODS
//#include "stdafx.h"
#include "defs.h"
#include "pars.h"

extern double mu; // nuclear mobility
extern double pp; // passing probability

Barrier::Barrier(double pos)
{
    position = pos;
    forcesum = 0.0;
}

Nucleus::Nucleus(double pos) : Barrier(pos)
{
    NMT = 0;
}

void Nucleus::move(double frc)
{
    position += mu*frc;
}

void Nucleus::nucleate()
{
    Microtubule* mnew;
    Barrier* trgt;
    enum ORIENTATION ornt = inrand() < 0.5 ? LEFT : RIGHT;   // Choose orientation
    bool pass = inrand() < pp;
    if(!pass)
    {
        trgt = ornt == LEFT ? bPrev : bNext;
    }
    else
    {
        trgt = ornt == LEFT ? leftWall : rightWall;
    }
    mnew = new Microtubule(ornt,trgt);
    MT.push_back(mnew);
    NMT++;
}

unsigned int Nucleus::get_pushing()
{
    unsigned int npush = 0;
    vector<Microtubule*>::iterator mt;
    for(mt = MT.begin(); mt != MT.end(); mt++)
    {
        if((*mt)->force > 0.0) npush++;
    }
    return npush;
}


//unsigned int Nucleus::get_waiting ()

//{
   //unsigned int wait = 0;
   // vector<Microtubule*>::iterator mt;
   //for(mt = MT.begin(); mt != MT.end(); mt++)
   // {
   //     if((*mt)->force > 0.0) wait++;
   // }
    //return wait;
 // }


void Nucleus::cleanup()
{
    vector<Microtubule*>::iterator mt;
    for(mt = MT.begin(); mt != MT.end(); )
    {
        if((*mt)->state == DORMANT)
        {
            delete *mt;
            mt = MT.erase(mt);
        }
        else
        {
            mt++;
        }
    }

}
