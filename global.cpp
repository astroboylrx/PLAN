//
//  global.cpp
//  PLATO: PLAneTesimal locatOr
//
//  Created by Rixin Li on 3/10/16.
//  Copyright © 2016 Rixin Li. All rights reserved.
//

#include "global.hpp"

/***************************/
/********** Timer **********/
/***************************/
/*! \fn Timer()
 *  \brief constructor */
Timer::Timer()
{
    timer_on_flag = 0;
    skip_time = 0;
}

double Timer::GetCurrentTime()
{
#ifdef MPI_ON
    return mpi->Wtime();
#else // MPI_ON
    return double(clock())/CLOCKS_PER_SEC;
#endif // MPI_ON
}

/*! \fn void StartTimer()
 *  \brief start the timer */
void Timer::StartTimer()
{
    ClearTimer();
    begin_time = GetCurrentTime();
    timer_on_flag = 1;
}

/*! \fn int void()
 *  \brief record a lap time and return index */
int Timer::Lap()
{
    if (timer_on_flag) {
        lap_time.push_back(GetCurrentTime() - begin_time - skip_time);
        return int(lap_time.size());
    } else {
        std::cerr << "Error: Timer is off. Cannot lap a time." << std::endl;
        return -1;
    }
}

/*! \fn void StopTimer()
 *  \brief stop the timer and return time */
void Timer::StopTimer()
{
    stop_time = GetCurrentTime() - begin_time - skip_time;
    lap_time.push_back(stop_time);
    timer_on_flag = 0;
}

/*! \fn void ClearTimer()
 *  \brief reset the timer */
void Timer::ClearTimer()
{
    begin_time = 0;
    timer_on_flag = 0;
    skip_time = 0;
    stop_time = 0;
    std::vector<double> TempVector;
    lap_time.swap(TempVector);
}

/*! \fn void ResumeTimer()
 *  \brief resume the timer */
void Timer::ResumeTimer()
{
    skip_time += GetCurrentTime() - stop_time;
    timer_on_flag = 1;
    stop_time = 0;
}

/*! \fn double GiveTime()
 *  \brief show the timer */
double Timer::GiveTime()
{
    if (timer_on_flag) {
        return GetCurrentTime() - begin_time - skip_time;
    } else {
        return stop_time;
    }
}

/*! \fn ~Timer()
 *  \brief destructor */
Timer::~Timer()
{
    std::vector<double> TempVector;
    lap_time.swap(TempVector);
}

/*********************************/
/********** MPI_Wrapper **********/
/*********************************/

#ifdef MPI_ON

/********** Constructor **********/
/*! \fn MPI_Wrapper()
 *  \brief constructor */
MPI_Wrapper::MPI_Wrapper()
{
    ;
}

/********** Initialization **********/
/*! \fn void Initialize(int argc, const char * argv[])
 *  \brief MPI initializaion */
void MPI_Wrapper::Initialization(int argc, const char * argv[])
{
    // common initialization of MPI
    MPI::Init(argc, (char **&)argv);
    world = MPI::COMM_WORLD;
    num_proc = world.Get_size();
    myrank = world.Get_rank();
    master = 0;
    
    // initialize file loop's parameters
    loop_begin = myrank;
    loop_end = myrank;
    loop_step = num_proc;
}

/********** Determine loop parameter **********/
/*! \fn void Determine_Loop(int num_file)
 *  \brief determine the begin/end/offset for file loop */
void MPI_Wrapper::DetermineLoop(int num_file)
{
    if (num_file < num_proc) {
        if (myrank > num_file - 1) {
            loop_end = -1;
        }
        loop_step = 1;
    } else {
        // in order to let master processor become available
        // in fact, no special effect, just for future dev
        loop_begin = num_proc - 1 - myrank;
        loop_end = num_file - 1;
    }
}

/********** Barrier **********/
/*! \fn void Barrier()
 *  \brief a wrapper of MPI Barrier */
void MPI_Wrapper::Barrier()
{
    world.Barrier();
}

/********** Wtime **********/
/*! \fn double Wtime()
 *  \brief a wrapper of MPI Wtime */
double MPI_Wrapper::Wtime()
{
    return MPI::Wtime();
}

/********** Msg Prefix **********/
/*! \fn std::string RankInfo()
 *  \brief return a string contains "Processor myrank: " */
std::string MPI_Wrapper::RankInfo()
{
    std::ostringstream oss;
    oss << "Processor " << myrank << ": ";
    return oss.str();
}

/********** Finalization **********/
/*! \fn void Finalize()
 *  \brief a wrapper of MPI Finalize() */
void MPI_Wrapper::Finalize()
{
    MPI::Finalize();
}

/********** Destructor **********/
/*! \fn ~MPI_info()
 *  \brief a destructor */
MPI_Wrapper::~MPI_Wrapper()
{
    ;
}

#endif // MPI_ON