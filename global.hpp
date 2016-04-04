//
//  global.hpp
//  PLATO: PLAneTesimal locatOr
//
//  Created by Rixin Li on 3/10/16.
//  Copyright © 2016 Rixin Li. All rights reserved.
//

#ifndef global_hpp
#define global_hpp
// C library first
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cassert>
#include <algorithm>
// C++ library
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <bitset>
#include <iomanip>
// other library
#include <unistd.h>
#include <getopt.h>

#define MPI_ON
#ifdef MPI_ON
#include "mpi.h"
#endif // MPI_ON

#ifndef OLDCPP
// Check c++11 support
#if __cplusplus <= 199711L
#error This program needs at least a C++11 compliant compiler (or define OLDCPP as an expedient).
#endif // __cplusplus
#endif // OLDCPP

/*! \class Timer
 *  \brief served as timer */
class Timer {
private:
    /*! \var double begin_time
     *  \brief the time when this class is initiated */
    double begin_time;
    
    /*! \var bool timer_on_flag
     *  \biref a flag indicate if this timer is on or off */
    bool timer_on_flag;
    
    /*! \var double skip_time
     *  \brief if timer is resumed, [stop->resume] is skipped */
    double skip_time;
    
    /*! \fn double GetCurrentTime()
     *  \brief get current time */
    double GetCurrentTime();
    
    /*! \var double end_time
     *  \brief the time when we stop this timer */
    double stop_time;
    
public:
    /*! \var std::vector<double> lap_time;
     *  \brief lap time, just like stopwatch */
    std::vector<double> lap_time;
    
    /*! \fn Timer()
     *  \brief constructor */
    Timer();
    
    /*! \fn void StartTimer()
     *  \brief start the timer */
    void StartTimer();
    
    /*! \fn int void()
     *  \brief record a lap time and return index */
    int Lap();
    
    /*! \fn void StopTimer()
     *  \brief stop the timer */
    void StopTimer();
    
    /*! \fn void ClearTimer()
     *  \brief reset the timer */
    void ClearTimer();
    
    /*! \fn void ResumeTimer()
     *  \brief resume the timer */
    void ResumeTimer();
    
    /*! \fn double GiveTime()
     *  \brief give the time (on) or [start->stop] (off) */
    double GiveTime();
    
    /*! \fn ~Timer()
     *  \brief destructor */
    ~Timer();
};

/*! \var extern std::vector<Timer> timer
 *  \brief declaration of global timer */
extern std::vector<Timer> timer;

/*! \enum TimerTypeIndex
 *  \brief served as index of the timer */
enum TimeTypeIndex {
    // put two underscore at first to avoid possbile name conflict
    __total_elapse_time = 0,      /*!< total elapsed time */
    __waiting_time,               /*!< time for waiting, used when MPI is on */
    __tmp_used_timer,             /*!< calculate time for temporary purpose */
    __time_type_count             /*!< the number of time type */
};

#ifdef MPI_ON
/*! \class MPI_Wrapper
 *  \brief served as wrappers of MPI routines plus related variables */
class MPI_Wrapper {
private:
    
public:
    /*! \var int num_proc
     *  \brief number of processors */
    int num_proc;
    
    /*! \var int rank, master
     *  \brief rank of this cpu / master cpu */
    int myrank, master;
    
    /*! \var int loop_begin, loop_end, loop_step
     *  \brief begin/end/step of the file loop handled by this cpu */
    int loop_begin, loop_end, loop_step;
    
    /*! \var std::vector<Timer> timer
     *  \brief used for recording time usage */
    std::vector<Timer> timer;
    
    /*! \var MPI::Intracomm world
     *  \brief a wrapper of MPI::COMM_WORLD */
    MPI::Intracomm world;
    
    /*! \var MPI::Status status
     *  \brief status used in MPI::COMM_WORLD.Recv function */
    MPI::Status status;
    
    /*! \fn MPI_Wrapper()
     *  \brief constructor */
    MPI_Wrapper();
    
    /*! \fn void Initialization(int argc, const char * argv[])
     *  \brief MPI initializaion */
    void Initialization(int argc, const char * argv[]);
    
    /*! \fn void Determine_Loop(int num_file)
     *  \brief determine the begin/end/step for file loop */
    void DetermineLoop(int num_file);
    
    /*! \fn int Barrier()
     *  \brief a wrapper of MPI::COMM_WORLD.Barrier */
    void Barrier();
    
    /*! \fn double Wtime()
     *  \brief a wrapper of MPI::Wtime */
    double Wtime();
    
    /*! \fn int Finalize()
     *  \brief a wrapper of MPI::Finalize */
    void Finalize();
    
    /*! \fn std::string RankInfo()
     *  \brief return a string contains "Processor myrank: " */
    std::string RankInfo();
    
    /*! \fn ~MPI_Wrapper()
     *  \brief destructor */
    ~MPI_Wrapper();
    
};

/*! \var extern MPI_Wrapper *mpi
 *  \brief declaration of gloal MPI wrapper */
extern MPI_Wrapper *mpi;

#endif // MPI_ON


#endif /* global_hpp */
