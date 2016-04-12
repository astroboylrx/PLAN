//
//  global.hpp
//  PLATO: PLAneTesimal locatOr
//
//  Created by Rixin Li on 3/10/16.
//  Copyright © 2016 Rixin Li. All rights reserved.
//

#ifndef global_hpp
#define global_hpp

// Include C libraries first
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cassert>
#include <algorithm>
// Include C++ libraries
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <bitset>
#include <iomanip>
#include <numeric>
// Include other libraries
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
#else // OLDCPP
#define nullptr NULL
#endif // OLDCPP


/******************************/
/********** I/O Part **********/
/******************************/

/*! \class IO_FileName
 *  \brief contains all I/O-related file path & names
 */
class IO_FileName {
private:
    
public:
    /*! \var std::string data_file_dir
     *  \brief path of directory that contains data files */
    std::string data_file_dir;
    
    /*! \var std::string data_file_basename
     *  \brief base name for data files */
    std::string data_file_basename;
    
    /*! \var std::string data_file_postname
     *  \brief post name for data files */
    std::string data_file_postname;
    
    /*! \var std::string output_file_path
     *  \brief path for ostream output */
    std::string output_file_path;
    
    /*! \var std::vector<std::string> lis_data_file_name
     *  \brief file names construced (for particle's lis files) */
    std::vector<std::string> lis_data_file_name;
    
#ifdef SMR_ON
    /*! \var std::string data_level
     *  \brief set level for Static-Mesh-Refinement Run */
    std::string data_level;
    
    /*! \var std::string data_domain
     *  \brief set level for Static-Mesh-Refinement Run */
    std::string data_domain;
#endif // SMR_ON
    
};

/*! \class IO_Operations
 *  \brief Handle all the [file] I/O operations */
class IO_Operations {
private:
    
public:
    /*! \var IO_FileName file_name
     *  \brief all I/O-related file path & names */
    IO_FileName file_name;
    
    /*! \var std::ostringstream out_content
     *  \brief this stores the content to output */
    std::ostringstream out_content;
    
    /*! \var std::ostringostream log_info
     *  \brief this stores the lof information */
    std::ostringstream log_info;
    
    /*! \var std::ostringstream error_message
     *  \brief this stores any temporary error message */
    std::ostringstream error_message;
    
    /*! \var int ostream_level
     *  \brief set level for different amount of output */
    int ostream_level;
    
    /*! \var int debug_flag
     *  \brief set this flag to make log_level = 2 */
    int debug_flag;
    
    /*! \var int start_num, end_num, interval
     *  \brief start/end number and interval fo the entir file loop */
    int start_num, end_num, interval;
    
    /*! \var int num_file
     *  \brief the number of files in total */
    int num_file;
    
    /*! \var int num_cpu
     *  \brief the number of processors used in simulation */
    int num_cpu;
    
    /*! \var int high_log_level_flag;
     *  \brief set this flag to obtain detailed running log */
    int high_log_level_flag;
    
    /*! \enum OutputLevel
     *  \brief served as output_level in Output() */
    enum OutputLevel {
        // put two underscore at first to avoid possbile name conflict
        __normal_output = 0,        /*!< only basic output */
        __more_output,              /*!< print more info for debugging */
        __even_more_output          /*!< give all possible info */
    };
    
    /*! \enum MPI_Level
     *  \brief served as mpi_level in Output() */
    enum MPI_Level {
        // put two underscore at first to avoid possbile name conflict
        __master_only = 0,          /*!< only master processor prints */
        __all_processors            /*!< all processors speak */
    };
    
    /*! \fn FileOperation()
     *  \brief constructor */
    IO_Operations();
    
    /*! \fn int Initialization(int argc, const char * argv[])
     *  \brief initialization */
    int Initialize(int argc, const char * argv[]);
    
    /*! \fn void Output(std::ostream &stream, std::ostringstream &content, OutputLevel &output_level, MPI_Level &mpi_level)
     *  \brief handle output by log level & MPI status */
    void Output(std::ostream &stream, std::ostringstream &content, const OutputLevel &output_level, const MPI_Level &mpi_level);
    
    /*! \fn int PrintUsage(const char *program_name)
     *  \brief print usage */
    void PrintUsage(const char *program_name);
    
    /*! \fn void PrintStars(std::ostream &stream, const OutputLevel &output_level)
     *  \brief print 80 * symbols as a divider line */
    void PrintStars(std::ostream &stream, const OutputLevel &output_level);
    
    /*! \fn void GenerateFilenames()
     *  \brief generate the name of data files for processing */
    void GenerateFilenames();
    
    
    
    /*! \fn void WriteResults()
     *  \brief write calculation results into files */
    void WriteResults();
    
    /*! \fn ~~IO_Operations()
     *  \brief destructor */
    ~IO_Operations();
};

/*! \var FileOperation *io_ops
 *  \brief global class handle all I/O stuff */
extern IO_Operations *io_ops;

/**********************************/
/********** Utility Part **********/
/**********************************/

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
