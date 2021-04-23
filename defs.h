/// Declarations of globals and common objects

#ifndef DEFS_H_INCLUDED
#define DEFS_H_INCLUDED

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <string>
#include <cmath>
#include <ctime>
#include <fstream>
#include <cstdlib>
#include <cstdio>

using namespace std;

class Barrier;

enum STATE {DORMANT, GROWING, SHRINKING};
enum ORIENTATION {LEFT = -1,RIGHT = 1};

extern double dlp;  // growth per time step
extern double dlm;  // shrinking per time step
extern double pcat; // probability of catastrophe per time step
extern double pres; // probability of rescue per time step
extern double mu; // nuclear mobility
extern Barrier* leftWall;
extern Barrier* rightWall;

inline double inrand(){return (double) rand()/RAND_MAX;}

class Barrier;

class Microtubule
{
   public:
       double length;
       double force;
       double waitingtime;
       enum STATE state;
       enum ORIENTATION orient;
       Barrier* target;

   public:
        Microtubule(void);
        Microtubule(enum ORIENTATION ornt, Barrier* trgt);
        void grow(double dlg);
        void shrink(double dls);
      //  void wait(double wts);
};

class Barrier
{
    public:
        double position;
        double forcesum;

    public:
        Barrier(double pos);
};

class Nucleus : public Barrier
{
    public:
        vector<Microtubule*> MT;
        Barrier* bPrev;
        Barrier* bNext;
        unsigned int NMT;

    public:
        Nucleus(double initpos);
        void move(double frc);
        void nucleate();
        unsigned int get_pushing();
        //unsigned int get_waiting();
        void cleanup();
};




#endif // DEFS_H_INCLUDED


