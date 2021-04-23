#ifndef PARS_H_INCLUDED
#define PARS_H_INCLUDED

/// CELL PARAMETERS

#define NUM_NUCLEI 4        // Number of nuclei
#define CELL_LENGTH 32      // Cell length [micron]
#define CELL_DIAMETER 3.5   // Cell diameter [micron] = 2*1.75
#define ETA  0.9            // Cytoplasmic viscosity [pN second / micron^2]


/// NUCLEUS PARAMETERS

#define NUC_DIAMETER  2.6   // Diameter of nucleus [micron] = 2*1.3

/// MICROTUBULE PARAMETERS


#define VG 0.033            // Unloaded growth speed [micron/second]*vg*dt
#define VS 0.15             // Shrinking speed [micron/second]
#define RCC 0.005           // Base catastrophe rate when pushing [1/second]
#define RC 0.005            // Spontaneous catastrophe rate;
#define RN 0.0275           // Nucleation rate (per nucleus?) [1/second]
#define CF 0.56             // Thermal compressibility factor [1/micron]
#define PP 0.725            // Passing probability [0-1]
#define FS 1.67             // Force scale of response [pN]

/// SIMULATION PARAMETERS

#define DT 1.0              // time step [seconds
#define KRUN 10             // Number of independent runs
#define SAMPLES 40          // Number of runs per setting
#define TMIN  100000        // Simulation time [minutes] !!
#define MAX_STEPS 6000000   // Maximum number of time steps (has to be > 60*TMIN!)
#define OUTPUT_STEPSIZE 100 // Report values spaced so many time steps apart
#define PI 3.141592654      // Pi to 10 decimals
#define INIT_POS 0          // Initial positions: 1 equally spaced; 0: in left corner (centrifugation)
#define CBINS 12            // Number of bins for nucleus positions
#define LBINS 15000         // Number of bins for microtubule lengths
#define LMAX  150.0         // Maximum length of microtubules in histogram [micron]
#define WBINS 500
#define WMAX  750.0
#define FBINS 200           // Number of bins for forces
#define FMAX  20.0          // Maximum force in histogram [pN]


#define TIMESTAMPED true                                    // create timestamped files if true
#define DATADIR_NAME "data/tetranucleate/2021/length_runs/"             // Data directory + filename prefix
#define PARFILE_NAME "parameters.txt"                       // Default parameter file name
#define TRAJECTORYFILE_NAME "trajectory.txt"                // Default trajectory file name
#define WAITFILE_NAME "wait.txt"
#define LENGTHFILE_NAME "length.txt"
#define HISTOGRAMFILE_NAME  "histogram.txt"             // Default histogram file name
#define AVERAGESFILE_NAME  "averages.txt"               // Default histogram file name


#endif // PARS_H_INCLUDED
