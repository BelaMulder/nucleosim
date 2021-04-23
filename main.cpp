//#include "stdafx.h"
#include "defs.h"
#include "pars.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <string>
#include <sstream>


//double ppval[KRUN] = {0.0}; // array with passing probabilities
//double rnval[KRUN] = {0.005, 0.01, 0.0275, 0.05, 0.1, 0.275, 0.5, 1.0}; // array with with nucleation rates
double rcval[KRUN] = {0.00005,0.0001,0.0005,0.001,0.005,0.01,0.02,0.03,0.04,0.05}; // array with catastrophe rates

// Similarly for catastrophe rates

unsigned int wSteps = 0;      // Number of waiting steps
unsigned int nSteps = 0;      // Number of simulation steps
unsigned int ksteps = 0;      // Number of simulation steps
string timestamp = "";        // Time-stamp
double dlg = 0.0;             // Unloaded length increase when growing per time step [micron]
double dls = 0.0;             // Length decrease when shrinking per time step [micron]
double rc = 0.0;              // Unloaded catastrophe rate
double pp = 0.0;              // passing probability
double rn = 0.0;              // Nucleation rate [1/second]
double pc = 0.0;              // Unloaded catastrophe probability per time step
double pcc = 0.0;             // Loaded catastrophe probability per time step;
double pn = 0.0;              // Nucleation probability per time step;
double mu = 0.0;              // Nuclear mobility [micron / second / pN]
double km = 0.0;              // Microtubule compressibility [pN / micron]
double rightWallPos = 0.0;
double leftWallPos = 0.0;
double initPos[NUM_NUCLEI] = {0.0};
Barrier* rightWall = NULL;
Barrier* leftWall = NULL;
Nucleus* NC[NUM_NUCLEI] = {NULL};

//unsigned int waiting;
//unsigned int maxWaitingMax=500;
//unsigned int wallWaitingHist=[500];
//int nucWaitingHist =[500];
//int waitingMax=0;

//unsigned int WWaiting[] = {{0}};

double cbinw = CELL_LENGTH/CBINS;
unsigned int NPositions[NUM_NUCLEI][CBINS] = {{0}}; // Histogram with nuclear positions
double fbinw = FMAX/FBINS;
unsigned int MForces[FBINS] = {0};                  // Histogram with forces
double lbinw = LMAX/LBINS;
unsigned int MLengths[LBINS] = {0};
unsigned int MWait[WBINS]= {0};                     // Histogram with length of growing (!) microtubules
double NTrajectories[NUM_NUCLEI][MAX_STEPS] ={{0.0}};        // Trajectories of nuclei
unsigned int MNumbers[NUM_NUCLEI][MAX_STEPS] = {{0}};       // Numbers of microtubules per nucleus
unsigned int MPushing[NUM_NUCLEI][MAX_STEPS] = {{0}};       // Numbers of force producing microtubules per nucleus
double wbinw = WMAX/WBINS;
const double fsens = 1/FS;

string trajfile = "";
ofstream traj_out;
string histfile = "";
ofstream hist_out;
string waitfile = "";
ofstream wait_out;
string lengthfile = "";
ofstream length_out;

void define_parameters();       // Calculate remaining simulation parameters
string create_timestamp();      // Create time-stamp
void write_parameters();        // Output parameter file
void write_trajectory();
void write_histogram();
void write_wait ();
void create_cell();             // Initialize simulation
void write_data(int k,int s);   // Output data
void nbin(int nc,double pos);   // Bin nuclear position
void wbin(double wait);
void lbin(double len);          // Bin microtubule length
void fbin(double frc);          // Bin microtubule force

inline double real_position(int nc,double pos){return pos+(nc+0.5)*NUC_DIAMETER;}     // Calculate real position in cell from reduced position

/// MAIN PROGRAM

int main(void)
{

for (int k=0; k < KRUN; k++)

{
    rc = rcval[k];

    for( int samp = 0; samp < SAMPLES; samp++)
    {


    srand(time(NULL));

    define_parameters();
    timestamp = create_timestamp();
    write_parameters();

    cout << "Run # "<< k << " rc: " << rc << " sample # " << samp << endl;


    trajfile = TIMESTAMPED ? (string) DATADIR_NAME + timestamp + (string) TRAJECTORYFILE_NAME : (string) DATADIR_NAME + (string) TRAJECTORYFILE_NAME;
    traj_out.open(trajfile.c_str());

    histfile = TIMESTAMPED ? (string) DATADIR_NAME + timestamp + (string) HISTOGRAMFILE_NAME : (string) DATADIR_NAME + (string) HISTOGRAMFILE_NAME;
    hist_out.open(histfile.c_str());


    waitfile = TIMESTAMPED ? (string) DATADIR_NAME + timestamp + (string) WAITFILE_NAME : (string) DATADIR_NAME + (string) WAITFILE_NAME;
    wait_out.open(waitfile.c_str());

    lengthfile = TIMESTAMPED ? (string) DATADIR_NAME + timestamp + (string) LENGTHFILE_NAME : (string) DATADIR_NAME + (string) LENGTHFILE_NAME;
    length_out.open(lengthfile.c_str());


    create_cell();

    /// MAIN PROGRAM LOOP

    for(unsigned int step = 0; step < nSteps; step++)
    {

        /// RECORD POSITIONS NUCLEI AND MICROTUBULE NUMBERS

//      cout << "step: " << step << " ";

      for(int nc = 0; nc < NUM_NUCLEI; nc++)
      {
           *(NTrajectories[nc]+step) = real_position(nc,NC[nc]->position);
           *(MNumbers[nc]+step) = NC[nc]->NMT;
           *(MPushing[nc]+step) = NC[nc]->get_pushing();

      }


        /// EVOLVE MICROTUBULES

        for(int nc = 0; nc < NUM_NUCLEI; nc++)
        {
            vector<Microtubule*>:: iterator mt;
            Nucleus* ncur = NC[nc];
            for(mt = ncur->MT.begin(); mt != ncur->MT.end(); mt++)
            {
                Microtubule* mcur = *mt;
                double strain = mcur->length-fabs(mcur->target->position-ncur->position);  // Determine degree of compression
                mcur->force = strain > 0.0 ? km*strain : 0.0;
                fbin(mcur->force);
                lbin(mcur->length);


                if(mcur->state == GROWING)
                {
                    double pick = ((double)rand()/(RAND_MAX));
                    double force = mcur->force;
                    if(force > 0.0)
                    {
                       double risk = exp(fsens*force);   // Multiplicative increase catastrophe rate due to force
                       if(pick < pcc*risk)
                        {
                            mcur->state = SHRINKING;
                            wbin(mcur->waitingtime);
                            mcur->waitingtime=0.0;
                        }
                        else
                        {
                            mcur->grow(dlg/risk);   // Decreased growth speed due to force
                            mcur->waitingtime +=DT;
                        }
                    }
                    else  // MT is unloaded
                    {
                        if(pick < pc)
                        {
                             mcur->state = SHRINKING;
                        }
                        else
                        {
                            mcur->grow(dlg);
                        }
                    }

                }
                else  // MT is shrinking
                {
                    mcur->shrink(dls);
                    if(mcur->length < 0.0)
                    {
                       mcur->state = DORMANT;
                       (ncur->NMT)--;
                    }
                }

            }

        }


        /// NUCLEATE NEW MICROTUBULES

        for(int nc = 0; nc < NUM_NUCLEI; nc++)


        {
            double fract = ((double)rand()/ (RAND_MAX));

            if(fract < pn)
            {
                NC[nc]->nucleate();
            }
        }

        /// CALCULATE FORCES ON NUCLEI

        for(int nc = 0; nc < NUM_NUCLEI; nc++)
        {
            vector<Microtubule*>:: iterator mt;
            Nucleus* ncur = NC[nc];
            for(mt = ncur->MT.begin(); mt != ncur->MT.end(); mt++)
            {
                Microtubule* mcur = *mt;
                ncur->forcesum -= mcur->orient*(mcur->force);              // Reaction force on parent nucleus
                mcur->target->forcesum += mcur->orient*(mcur->force);      // Force on target
            }
        }

        /// MOVE NUCLEI

        for(int nc = 0; nc < NUM_NUCLEI; nc++)
        {
            double newpos = NC[nc]->position+mu*(NC[nc]->forcesum*DT);
            if(newpos > NC[nc]->bPrev->position && newpos < NC[nc]->bNext->position) NC[nc]->position = newpos; // Only move when there is room
            nbin(nc,real_position(nc,NC[nc]->position));
            NC[nc]->forcesum = 0.0;      // Forces are "spent" after the move
        }


        /// REMOVE DORMANT MICROTUBULES

     for(int nc = 0; nc < NUM_NUCLEI; nc++) NC[nc]->cleanup();


    }

    /// END MAIN PROGRAM LOOP

    write_data(k,samp);

  // close open file handles
    traj_out.close();
    hist_out.close();
    wait_out.close();
    length_out.close();

    } // end of sample
} // end of rc

    return 0;
}


/// FUNCTION DEFINITIONS


void define_parameters()
{
    nSteps = (unsigned int) floor((double) TMIN*60/DT);
    dlg = VG*DT;
    dls = VS*DT;
    pc = rc*DT;
    pcc = RCC*DT;
    pn = RN*DT;
    pp = PP;
    leftWallPos = 0.0;
    rightWallPos = CELL_LENGTH-NUM_NUCLEI*NUC_DIAMETER;
    if(INIT_POS == 1)
    {
        for(int n = 0; n < NUM_NUCLEI; n++){initPos[n]=(double) (rightWallPos-leftWallPos)/(NUM_NUCLEI+1)*(n+1);}
    }
    else
    {
        for(int n = 0; n < NUM_NUCLEI; n++){initPos[n]=0.01;}
    }
    km = FS*CF;
    double eps = (CELL_DIAMETER-NUC_DIAMETER)/NUC_DIAMETER;
    mu = 4*pow(eps,2.5)/(9*PI*PI*sqrt(2.0)*ETA*0.5*NUC_DIAMETER);

}

void create_cell()
{
    rightWall = new Barrier(rightWallPos);
    leftWall = new Barrier(leftWallPos);
    for(int n = 0; n < NUM_NUCLEI; n++)
    {
        NC[n] = new Nucleus(initPos[n]);
    }
    NC[0]-> bPrev = leftWall;
    if(NUM_NUCLEI > 1){
    NC[0]->bNext = NC[1];
    for(int nc = 1; nc < NUM_NUCLEI-1; nc++)
    {
        NC[nc]->bPrev = NC[nc-1];
        NC[nc]->bNext = NC[nc+1];
    }
    NC[NUM_NUCLEI-1]->bPrev = NC[NUM_NUCLEI-2];};
    NC[NUM_NUCLEI-1]->bNext = rightWall;
}


string create_timestamp()
{
    time_t curtime;
    struct timeval timev;
    char buffer[80];
    gettimeofday(&timev,NULL);
    curtime = timev.tv_sec;
    strftime(buffer,80,"%d-%m-%Y_%H-%M-%S",localtime(&curtime));
    string stamp = string(buffer);
    ostringstream ss;
    ss << timev.tv_usec;
    string usecs = ss.str();
    stamp += "-"+usecs+"_";
    return stamp;
}

void write_parameters()
{
    string parfile;
    parfile = TIMESTAMPED ? (string) DATADIR_NAME+timestamp+(string) PARFILE_NAME : (string) DATADIR_NAME+ (string) PARFILE_NAME;
    ofstream p_out (parfile.c_str());
    p_out << fixed;
    p_out << "tstamp "  << timestamp << endl;
    p_out << "Nc     "  << NUM_NUCLEI << endl;
    p_out << "vg     "  << VG << endl;
    p_out << "vs     "  << VS << endl;
    p_out << "rcc    "  << RCC <<endl;
    p_out << "rc     "  << rc << endl;
    p_out << "rn     "  << RN << endl;
    p_out << "fs     "  << FS << endl;
    p_out << "cf     "  << CF << endl;
    p_out << "km     "  << km << endl;
    p_out << "pp     "  << pp << endl;
    p_out << "Dn     "  << NUC_DIAMETER << endl;
    p_out << "L      "  << CELL_LENGTH << endl;
    p_out << "Dc     "  << CELL_DIAMETER << endl;
    p_out << "eta    "  << ETA << endl;
    p_out << "mu     "  << mu << endl;
    p_out << "tmin   "  << TMIN << endl;
    p_out << "dt     "  << DT << endl;
    p_out << "nsteps "  << nSteps << endl;
    p_out << "init   "  << INIT_POS << endl;

    p_out.close();
}

void nbin(int nc, double pos)
{
    double realpos = pos;
    if(realpos > 0 && realpos < CELL_LENGTH)
    {
       NPositions[nc][(int) floor(realpos/cbinw)]++;
       //cout << realpos << endl;
    }
   else
   {
       cout << "Position out of range" << endl;
   }
}

void wbin (double waitingtime)
{

if(waitingtime < WMAX)

{
   MWait[(int) floor(waitingtime/wbinw)]++;
}

else
    {
      cout << "too long waiting" << endl;
    }

}

void lbin(double length)
{

if(length < LMAX)

{
    MLengths[(int) floor(length/lbinw)]++;
}
else
    {
      cout << "too long microtubule" << endl;
    }

}


void fbin(double frc)
{
if(frc < FMAX)
    {
       MForces[(int) floor(frc/fbinw)]++;
    }
    else
    {
      cout << "too large force" << endl;
    }
}

void write_data (int k,int s)
{

    // cout << k << "   " << s << endl;
    for(unsigned int step =0; step < nSteps;  step += OUTPUT_STEPSIZE)
    {
        traj_out << step*DT;

        for(int nc = 0; nc < NUM_NUCLEI; nc++)
        {
            traj_out << " " << *(NTrajectories[nc]+step) << " " << *(MNumbers[nc]+step) << " " << *(MPushing[nc]+step);
        }
        traj_out << endl;
    }



    for(int bin = 0; bin < CBINS; bin++)
    {
        hist_out << (bin)*cbinw;
        for(int nc = 0; nc < NUM_NUCLEI; nc++)
        {
            hist_out << " " << NPositions[nc][bin];
        }
         hist_out << endl;
    }


      for(int bin = 0; bin < WBINS; bin++)
      {
        wait_out << (bin)*wbinw << " " << MWait[bin] << endl;
      }

    for(int bin = 0; bin < LBINS; bin++)
    {
       length_out << (bin)*lbinw << " " << MLengths[bin] << endl;
    }


}




