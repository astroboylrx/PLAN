//
//  global.hpp
//  PLAN: PLantesimal ANalyzer (a better version of project PLATO)
//
//  Created by Rixin Li on 4/26/16.
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
#ifndef OLDCPP
#include <array>
#include <type_traits>
#endif
// Include other libraries
#include <unistd.h>
#include <getopt.h>

#define MPI_ON // Comment out this line before committing!!
#ifdef MPI_ON // Only enable these options during compilation
#include "mpi.h"
#endif // MPI_ON

#ifndef OLDCPP
// Check c++11 support
#if __cplusplus <= 199711L
#error This program needs at least a C++11 compliant compiler (or define OLDCPP as an expedient, which is not guaranteed to work so far).
#endif // __cplusplus
#else // OLDCPP
// In C++, use const instead of defining macros as constants
const void *nullptr = NULL;
#endif // OLDCPP


/*********************************************/
/**********Basic_IO_Operations Part **********/
/*********************************************/

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

/*! \class PhysicalQuantities
 *  \brief all physical quantities that might needed to be calculated */
class PhysicalQuantities {
public:
    /*! \var float time;
     *  \brief simulation time */
    float time;
    
    /*! \var float dt;
     *  \brief simulation time */
    float dt;
    
    /*! \var std::vector<double> particle_scale_height;
     *  \brief particle scale height for all particle sizes */
    std::vector<double> particle_scale_height;
    
    /*! \var double max_particle_density
     *  \brief maximum particle density: $\rho_p$ */
    float max_particle_density {0.0};
};

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

/*! \class IO_Flags
 *  \brief all possible flags used in executation */
class IO_Flags {
public:
    /*! \var int debug_flag
     *  \brief set this flag to make log_level = 2 */
    int debug_flag {0};
    
    /*! \var int verbose_flag
     *  \brief set this flag to make log_level = 1 */
    int verbose_flag {0};
    
    /*! \var int combined_flag
     *  \brief set this flag to read combined data */
    int combined_flag {0};
    
};

/*! \class Basic_IO_Operations
 *  \brief Handle command line and most of I/O functions together with MPI_Wrapper
 *  PS: reading method for lis files is in ParticleSet class */
class Basic_IO_Operations {
private:
    /*! \var int ostream_level
     *  \brief set level for different amount of output */
    int ostream_level {__normal_output};
    
public:
    /*! \var IO_FileName file_name
     *  \brief all I/O-related file path & names */
    IO_FileName file_name;
    
    /*! \var std::vector<PhysicalQuantities> physical_quantities;
     *  \brief all physical quantities that might needed to be calculated */
    std::vector<PhysicalQuantities> physical_quantities;
    
    /*! \var IO_Flags flags
     *  \brief all possible flags used in executation */
    IO_Flags flags;
    
    /*! \var std::ostringstream out_content
     *  \brief this stores the content to output */
    std::ostringstream out_content;
    
    /*! \var std::ostringostream log_info
     *  \brief this stores the lof information */
    std::ostringstream log_info;
    
    /*! \var std::ostringstream error_message
     *  \brief this stores any temporary error message */
    std::ostringstream error_message;
    
    /*! \var int start_num, end_num, interval
     *  \brief start/end number and interval fo the entir file loop */
    int start_num, end_num, interval;
    
    /*! \var int num_file
     *  \brief the number of files in total */
    int num_file;
    
    /*! \var int num_cpu
     *  \brief the number of processors used in simulation */
    int num_cpu;
    
    /*! \var int width {15}
     *  \brief set default width of one data unit */
    int width {15};
    
    /*! \var int column {4}
     *  \brief data column in result file */
    int column {2};
    
    /*! \fn Basic_IO_Operations()
     *  \brief constructor */
    Basic_IO_Operations();
    
    /*! \fn int Initialization(int argc, const char * argv[])
     *  \brief initialization */
    int Initialize(int argc, const char * argv[]);
    
    /*! \fn void Output(std::ostream &stream, std::ostringstream &content, const OutputLevel &output_level, const MPI_Level &mpi_level)
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
    
    /*! \fn inline void Reset(std::ostringstream &content)
     *  \brief reset the content of ostringstream */
    inline void Reset(std::ostringstream &content) {
        content.str(std::string());
        content.clear();
        //std::ostringstream().swap(content); // sometimes need c++14
    }
    
    /*! \fn ~Basic_IO_Operations()
     *  \brief destructor */
    ~Basic_IO_Operations();
};

/*! \var FileOperation *progIO
 *  \brief Handle command line and most of I/O functions together with MPI_Wrapper */
extern Basic_IO_Operations *progIO;

/**********************************/
/********** Utility Part **********/
/**********************************/

/*! \class MPI_Wrapper
 *  \brief served as wrappers of MPI routines plus related variables
 *  Note that even if running as a serial program, this still works since MPI_Wrapper takes care of it. The goal of this encapsulation is to reduce the use of (annoying) "#ifdef MPI_ON" in main function. */
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
    
#ifdef MPI_ON
    /*! \var MPI_Comm world
     *  \brief a wrapper of MPI_COMM_WORLD */
    MPI_Comm world;
    
    /*! \var MPI_Status status
     *  \brief MPI status */
    MPI_Status status;
    
    /*! \var MPI_Offset offset
     *  \brief stream offset used during parallel file I/O */
    MPI_Offset offset {0};
    
    /*! \var MPI_Offset header_offset
     *  \brief stream offset used during parallel file I/O */
    MPI_Offset header_offset {0};
    
#ifdef OLDCPP
    typedef MPI_File file_obj;
#else // OLDCPP
    using file_obj = MPI_File;
#endif // OLDCPP
    
#else // MPI_ON
    
#ifdef OLDCPP
    typedef std::ofstream file_obj
#else // OLDCPP
    using file_obj = std::ofstream;
#endif // OLDCPP
    
#endif // MPI_ON
    
    /*! \var file_obj result_file
     *  \brief default file to output result */
    file_obj result_file;
    
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
     *  \brief a wrapper of MPI_Barrier */
    void Barrier();
    
    /*! \fn int Finalize()
     *  \brief a wrapper of MPI_Finalize */
    void Finalize();
    
    /*! \fn std::string RankInfo()
     *  \brief return a string contains "Processor myrank: " */
    std::string RankInfo();
    
    /*! \fn void OpenFile(file_obj &__file, std::string filename)
     *  \brief open file */
    void OpenFile(file_obj &__file, std::string file_name);
    
    /*! \fn void WriteSingleFile(file_obj &__file, std::ostringstream &content)
     *  \brief all processor write into a file, use with cautions -> read the assumptions in descriptions
     *  When MPI_ON is on, this function assumes that you only write file header if specifying __master_only, and it assumes that every processor are writing the same amout of chunk into the file every time */
    void WriteSingleFile(file_obj &__file, std::ostringstream &content, const MPI_Level &mpi_level);
    
    /*! \fn void WriteItsOwnFile(file_obj &__file, std::ostringstream &content)
     *  \brief all processor write to its own file */
    void WriteItsOwnFile(file_obj &__file, std::ostringstream &content);
    
    /*! \fn void CloseFile(file_obj &__file)
     *  \brief close the file */
    void CloseFile(file_obj &__file);
    
    /*! \fn ~MPI_Wrapper()
     *  \brief destructor */
    ~MPI_Wrapper();
    
};

/*! \var extern MPI_Wrapper *mpi
 *  \brief declaration of gloal MPI wrapper */
extern MPI_Wrapper *mpi;

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

/*! \fn template <typename T> T MaxOf(const T &a, const T &b)
 *  \brief return the larger one of a and b */
template <typename T>
T MaxOf(const T &a, const T &b) {
    return std::max(a, b);
}

/*! \fn template <typename T, typename... Args> T MaxOf(const T &a, const T &b, Args... args)
 *  \brief return the maximum one in a list, downgrade to MaxOf(a, b) */
template <typename T, typename... Args>
T MaxOf(const T &a, const T &b, Args... args) {
    return MaxOf(std::max(a, b), args...);
}

/*! \fn template <typename T> T MinOf(const T &a, const T &b)
 *  \brief return the smaller one of a and b */
template <typename T>
T MinOf(const T &a, const T &b) {
    return std::min(a, b);
}

/*! \fn template <typename T, typename... Args> T MinOf(const T &a, const T &b, Args... args)
 *  \brief return the minimum one in a list, downgrade to MaxOf(a, b) */
template <typename T, typename... Args>
T MinOf(const T &a, const T &b, Args... args) {
    return MaxOf(std::min(a, b), args...);
}

#endif /* global_hpp */
