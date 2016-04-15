//
//  global.cpp
//  PLATO: PLAneTesimal locatOr
//
//  Created by Rixin Li on 3/10/16.
//  Copyright © 2016 Rixin Li. All rights reserved.
//

#include "global.hpp"


/***********************************/
/********** IO_Operations **********/
/***********************************/

/*! \fn IO_Operations()
 *  \brief constructor */
IO_Operations::IO_Operations()
{
    // default log level is 0
    ostream_level = __normal_output;
    // can be changed by options during compilation
#ifdef MORE_OUTPUT
    ostream_level = __more_output;
#endif
#ifdef EVEN_MORE_OUTPUT
    // alternative, set debug_flag to achieve this
    ostream_level = __even_more_output;
#endif
}

/*! \fn int Initialize(int argc, const char * argv[])
 *  \brief initialization */
int IO_Operations::Initialize(int argc, const char * argv[])
{
    // the following manipulation of command line options is from getopb.h
    // an alternative choice is boost library: Boost.Program_options
    
    // set all flags to 0 at first
    debug_flag = 0;
    //Specifying the expected options
    static struct option long_options[] = {
        // These options set a flag
        {"Debug", no_argument, &debug_flag, 1},
        // These options don't set a flag
        {"num_cpu", required_argument, 0, 'c'},
        {"data_dir", required_argument, 0, 'i'},
        {"basename", required_argument, 0, 'b'},
        {"postname", required_argument, 0, 'p'},
        {"file_num", required_argument, 0, 'f'},
        {"output", required_argument, 0, 'o'},
#ifdef SMR_ON
        {"level", required_argument, 0, 'l'},
        {"domain", required_argument, 0, 'd'},
#endif // SMR_ON
        // End
        {0,0,0,0}
    };
    
    log_info << "Verifying command-line-input arguments: \n";
    if (argc < 13) {
        PrintUsage(argv[0]);
    } else {
        int tmp_option;
        while (1) {
            // getopt_long stores the option
            int option_index = 0;
#ifndef SMR_ON
            // remember add ":" after the letter means this option has argument after it
            tmp_option = getopt_long(argc, (char *const *)argv, "c:i:b:p:f:o:", long_options, &option_index);
#else // SMR_ON
            tmp_option = getopt_long(argc, (char *const *)argv, "c:i:b:p:f:o:l:d:", long_options, &option_index);
#endif // SMR_ON
            if (tmp_option == -1) {
                break;
            }
            
            switch (tmp_option) {
                case 0: {
                    // if this option set a flag, do nothing else now
                    if (long_options[option_index].flag != 0) {
                        break;
                    }
                    log_info << "option " << long_options[option_index].name;
                    if (optarg) {
                        out_content << " with arg " << optarg;
                    }
                    log_info << "\n";
                    break;
                }
                case 'c': {
                    std::istringstream tmp_iss;
                    tmp_iss.str(optarg);
                    tmp_iss >> num_cpu;
                    log_info << "num_cpu is " << num_cpu << "\n";
                    break;
                }
                case 'i': {
                    file_name.data_file_dir.assign(optarg);
                    log_info << "data_file_dir is " << file_name.data_file_dir << "\n";
                    break;
                }
                case 'b': {
                    file_name.data_file_basename.assign(optarg);
                    log_info << "data_file_basename is " << file_name.data_file_basename << "\n";
                    break;
                }
                case 'p': {
                    file_name.data_file_postname.assign(optarg);
                    log_info << "data_file_postname is " << file_name.data_file_postname << "\n";
                    break;
                }
                case 'f': {
                    std::string tmp_str;
                    tmp_str.assign(optarg);
                    size_t pos1 = tmp_str.find_first_of(':');
                    size_t pos2 = tmp_str.find_last_of(':');
                    std::istringstream tmp_iss;
                    char tmp_char;
                    if (pos1 == pos2) {
                        tmp_iss.str(optarg);
                        tmp_iss >> start_num >> tmp_char >> end_num;
                        interval = 1;
                    } else {
                        tmp_iss.str(tmp_str.substr(0, pos2));
                        tmp_iss >> start_num >> tmp_char >> end_num;
                        tmp_iss.str(tmp_str.substr(pos2+1));
                        tmp_iss >> interval;
                    }
                    if (start_num < 0) {
                        error_message << "The start number should be positive (Auto fix to 0)" << std::endl;
                        Output(std::cerr, error_message, __normal_output, __master_only);
                        start_num = 0;
                    }
                    if (end_num < start_num) {
                        error_message << "The end number should be larger or equal than the start number. (Auto fix to start number)." << std::endl;
                        Output(std::cerr, error_message, __normal_output, __master_only);
                        end_num = start_num;
                    }
                    if (interval == 0) {
                        error_message << "The interval should be non-zero. (Auto fix to 1)" << std::endl;
                        Output(std::cerr, error_message, __normal_output, __master_only);
                        interval = 1;
                    }
                    num_file = (end_num - start_num) / interval + 1;
                    log_info << "start_num = " << start_num << ", end_num = " << end_num << ", interval = " << interval << ", num_file = " << num_file << "\n";
                    break;
                } // case 'f'
                case 'o': {
                    file_name.output_file_path.assign(optarg);
                    log_info << "output_file_path is " << file_name.output_file_path << "\n";
                    break;
                }
#ifdef SMR_ON
                case 'l': {
                    file_name.data_level.assign(optarg);
                    log_info << "data_level is " << file_name.data_level << "\n";
                    break;
                }
                case 'd': {
                    file_name.data_domain.assign(optarg);
                    log_info << "data_domain is " << file_name.data_domain << "\n";
                    break;
                }
#endif // SMR_ON
                case '?': {
                    
                    if (optopt == 'c' || optopt == 'i' || optopt == 'b' || optopt == 'p' || optopt == 'f' || optopt == 'o'
#ifdef SMR_ON
                        || optopt == 'l' || optopt == 'd'
#endif // SMR_ON
                        ) {
                        error_message << "Error: Option -" << optopt << " requires an arugment.\n";
                        Output(std::cerr, error_message, __normal_output, __master_only);
                    } else if (isprint (optopt)) {
                        error_message << "Error: Unknown option -" << optopt << "\n";
                        Output(std::cerr, error_message, __normal_output, __master_only);
                    } else {
                        error_message << "Error: Unknown option character " << optopt << "\n";
                        Output(std::cerr, error_message, __normal_output, __master_only);
                    }
                    exit(2); // cannot execute
                }
                default: {
                    error_message << tmp_option << "\nError: Argument wrong.\n";
                    Output(std::cerr, error_message, __normal_output, __master_only);
                    exit(2); // cannot execute
                }
            }
        }
        
        if (debug_flag) {
            ostream_level = __even_more_output;
        }
        Output(std::clog, log_info, __more_output, __master_only);
        PrintStars(std::clog, __more_output);
        if (optind < argc) {
            log_info << "Non-option ARGV-elements: ";
            while (optind < argc) {
                log_info << argv[optind++];
            }
            log_info << "\n";
            Output(std::clog, log_info, __normal_output, __master_only);
        }
    } // if (argc < 13)
    
    GenerateFilenames();
    
    return 0;
}

/*! \fn void PrintUsage(const char *program_name)
 *  \brief print usage if required argument is missing */
void IO_Operations::PrintUsage(const char *program_name)
{
    out_content << "USAGE: " << program_name << " -c <num_cpu> -i <data_dir> -b <basename> -p <postname>  -f <range(f1:f2)|range_step(f1:f2:step)> -o <output> [flags]\n" << "Example: ./plato -c 64 -i ./ -b Par_Strat3d -p ds -f 170:227 -o result.txt\n" << "Use --Debug to obtain more information during executation.";
    Output(std::cout, out_content, __normal_output, __master_only);
    exit(2); // cannot execute
}

/*! \fn void Output(std::ostream &stream, std::ostringstream &content, const OutputLevel &output_level, const MPI_Level &mpi_level)
 *  \brief handle output by log level & MPI status */
void IO_Operations::Output(std::ostream &stream, std::ostringstream &content, const OutputLevel &output_level, const MPI_Level &mpi_level)
{
    if (ostream_level >= output_level) {
#ifdef MPI_ON
        if (mpi_level == __master_only) {
            if (mpi->myrank == mpi->master) {
                stream << content.str() << std::flush;
            }
        } else if (mpi_level == __all_processors) {
            stream << mpi->RankInfo() << content.str() << std::flush;
        }
#else // MPI_ON
        stream << content.str() << std::flush;
#endif // MPI_ON
    }
    
    // The most elegant way to clear state and empty content
    // However, this requires those compilers that support c++11
    //std::ostringstream().swap(content);
    
    // alternatively, use
    content.str(std::string()); // empty the content
    content.clear(); // clear error state if any
}

/*! \fn void PrintStars(std::ostream &stream, const OutputLevel &output_level)
 *  \brief print 80 * symbols as a divider line */
void IO_Operations::PrintStars(std::ostream &stream, const OutputLevel &output_level)
{
    out_content << std::setw(80) << std::setfill('*') << "*\n";
    Output(stream, out_content, output_level, __master_only);
}

/*! \fn void GenerateFilenames()
 *  \brief generate the name of data files for processing */
void IO_Operations::GenerateFilenames()
{
    if (*file_name.data_file_dir.rbegin() != '/') {
        file_name.data_file_dir.push_back('/');
    }
    
    file_name.lis_data_file_name.reserve(num_file * num_cpu);
    log_info << "Verifying generated data file names (only id0 and id[max]):\n";
    
    for (int num = start_num; num != end_num+1; num += interval) {
        std::stringstream formatted_num;
        formatted_num << std::setw(4) << std::setfill('0') << num;
        
        file_name.lis_data_file_name.push_back(file_name.data_file_dir+"id0/"+file_name.data_file_basename+"."+formatted_num.str()+"."+file_name.data_file_postname+".lis");
        log_info << file_name.lis_data_file_name.back() << "\n";
        for (int id = 0; id != num_cpu; id++) {
            file_name.lis_data_file_name.push_back(file_name.data_file_dir+"id0/"+file_name.data_file_basename+"-id"+std::to_string(id+1)+"."+formatted_num.str()+"."+file_name.data_file_postname+".lis");
        }
        log_info << file_name.lis_data_file_name.back() << "\n";
    }
    Output(std::clog, log_info, __even_more_output, __master_only);
    PrintStars(std::clog, __even_more_output);
    
}

/*! \fn ~~IO_Operations()
 *  \brief destructor */
IO_Operations::~IO_Operations()
{
    ;
}

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

/*! \fn double GetCurrentTime()
 *  \brief get current time */
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

/*! \fn MPI_Wrapper()
 *  \brief constructor */
MPI_Wrapper::MPI_Wrapper()
{
    ;
}

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

/*! \fn void Barrier()
 *  \brief a wrapper of MPI Barrier */
void MPI_Wrapper::Barrier()
{
    world.Barrier();
}

/*! \fn double Wtime()
 *  \brief a wrapper of MPI Wtime */
double MPI_Wrapper::Wtime()
{
    return MPI::Wtime();
}

/*! \fn std::string RankInfo()
 *  \brief return a string contains "Processor myrank: " */
std::string MPI_Wrapper::RankInfo()
{
    std::ostringstream oss;
    oss << "Processor " << myrank << ": ";
    return oss.str();
}

/*! \fn void Finalize()
 *  \brief a wrapper of MPI Finalize() */
void MPI_Wrapper::Finalize()
{
    MPI::Finalize();
}

/*! \fn ~MPI_info()
 *  \brief a destructor */
MPI_Wrapper::~MPI_Wrapper()
{
    ;
}

#endif // MPI_ON

