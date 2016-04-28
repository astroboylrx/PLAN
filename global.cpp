//
//  global.cpp
//  PLAN: PLantesimal ANalyzer (a better version of project PLATO)
//
//  Created by Rixin Li on 4/26/16.
//  Copyright © 2016 Rixin Li. All rights reserved.
//

#include "global.hpp"


/*****************************************/
/********** Basic_IO_Operations **********/
/*****************************************/

/*! \fn Basic_IO_Operations()
 *  \brief constructor */
Basic_IO_Operations::Basic_IO_Operations()
{
    ;
}

/*! \fn int Initialize(int argc, const char * argv[])
 *  \brief initialization */
int Basic_IO_Operations::Initialize(int argc, const char * argv[])
{
    // the following manipulation of command line options is from getopb.h
    // an alternative choice is boost library: Boost.Program_options
    
    //Specifying the expected options
    static struct option long_options[] = {
        // These options set a flag
        {"Debug", no_argument, &flags.debug_flag, 1},
        {"Verbose", no_argument, &flags.verbose_flag, 1},
        {"Combined", no_argument, &flags.combined_flag, 1},
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
        
        // determine ostream level before first Output()
        if (flags.verbose_flag) {
            ostream_level = __more_output;
        }
        if (flags.debug_flag) {
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
    
    physical_quantities.resize(num_file);
    GenerateFilenames();
    
    return 0;
}

/*! \fn void PrintUsage(const char *program_name)
 *  \brief print usage if required argument is missing */
void Basic_IO_Operations::PrintUsage(const char *program_name)
{
    out_content << "USAGE: " << program_name
    << " -c <num_cpu> -i <data_dir> -b <basename> -p <postname>  -f <range(f1:f2)|range_step(f1:f2:step)> -o <output> [flags]\n"
    << "Example: ./plato -c 64 -i ./ -b Par_Strat3d -p ds -f 170:227 -o result.txt\n"
    << "Use --Verbose to obtain more output during executation"
    << "Use --Debug to obtain all possbile output during executation"
    << "Use --Combined to deal with combined lis files (from all processors";
    out_content << std::endl;
    Output(std::cout, out_content, __normal_output, __master_only);
    exit(2); // cannot execute
}

/*! \fn void Output(std::ostream &stream, std::ostringstream &content, const OutputLevel &output_level, const MPI_Level &mpi_level)
 *  \brief handle output by log level & MPI status */
void Basic_IO_Operations::Output(std::ostream &stream, std::ostringstream &content, const OutputLevel &output_level, const MPI_Level &mpi_level)
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
void Basic_IO_Operations::PrintStars(std::ostream &stream, const OutputLevel &output_level)
{
    out_content << std::setw(80) << std::setfill('*') << "*\n";
    Output(stream, out_content, output_level, __master_only);
}

/*! \fn void GenerateFilenames()
 *  \brief generate the name of data files for processing */
void Basic_IO_Operations::GenerateFilenames()
{
    if (*file_name.data_file_dir.rbegin() != '/') {
        file_name.data_file_dir.push_back('/');
    }
    
    if (flags.combined_flag) {
        file_name.lis_data_file_name.reserve(num_file);
        log_info << "Verifying generated data file names (only the first one and last one):\n";
        
        for (int num = start_num; num != end_num+1; num += interval) {
            std::stringstream formatted_num;
            formatted_num << std::setw(4) << std::setfill('0') << num;
            
            file_name.lis_data_file_name.push_back(file_name.data_file_dir+file_name.data_file_basename+"."+formatted_num.str()+"."+file_name.data_file_postname+".lis");
            if (num == start_num || num == end_num) {
                log_info << file_name.lis_data_file_name.back() << "\n";
            }
        }
        
    } else {
        file_name.lis_data_file_name.reserve(num_file * num_cpu);
        log_info << "Verifying generated data file names (only id0 and id[max]):\n";
        
        for (int num = start_num; num != end_num+1; num += interval) {
            std::stringstream formatted_num;
            formatted_num << std::setw(4) << std::setfill('0') << num;
            
            file_name.lis_data_file_name.push_back(file_name.data_file_dir+"id0/"+file_name.data_file_basename+"."+formatted_num.str()+"."+file_name.data_file_postname+".lis");
            log_info << file_name.lis_data_file_name.back() << "\n";
            for (int id = 1; id != num_cpu; id++) {
                file_name.lis_data_file_name.push_back(file_name.data_file_dir+"id"+std::to_string(id)+"/"+file_name.data_file_basename+"-id"+std::to_string(id)+"."+formatted_num.str()+"."+file_name.data_file_postname+".lis");
            }
            log_info << file_name.lis_data_file_name.back() << "\n";
        }
    }
    
    Output(std::clog, log_info, __even_more_output, __master_only);
    PrintStars(std::clog, __even_more_output);
    
}

/*! \fn ~Basic_IO_Operations()
 *  \brief destructor */
Basic_IO_Operations::~Basic_IO_Operations()
{
    ;
}

/*********************************/
/********** MPI_Wrapper **********/
/*********************************/

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
#ifdef MPI_ON
    // common initialization of MPI
    MPI_Init(&argc, (char***)&argv);
    world = MPI_COMM_WORLD;
    MPI_Comm_size(world, &num_proc);
    MPI_Comm_rank(world, &myrank);
#else // MPI_ON
    num_proc = 1;
    myrank = 0;
#endif // MPI_ON
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
    // use reverse order
    loop_begin = num_proc - 1 - myrank;
    loop_end = num_file - 1;
    if (loop_begin >= num_proc) {
        loop_end = -1;
    }
}

/*! \fn void Barrier()
 *  \brief a wrapper of MPI Barrier */
void MPI_Wrapper::Barrier()
{
#ifdef MPI_ON
    MPI_Barrier(world);
#endif // MPI_ON
}

/*! \fn std::string RankInfo()
 *  \brief return a string contains "Processor myrank: " */
std::string MPI_Wrapper::RankInfo()
{
#ifdef MPI_ON
    std::ostringstream oss;
    oss << "Processor " << myrank << ": ";
    return oss.str();
#else // MPI_ON
    return std::string();
#endif // MPI_ON
}

/*! \fn void Finalize()
 *  \brief a wrapper of MPI Finalize() */
void MPI_Wrapper::Finalize()
{
#ifdef MPI_ON
    MPI_Finalize();
#endif // MPI_ON
}

/*! \fn void OpenFile(file_obj &__file, std::string filename)
 *  \brief open file for data writing */
void MPI_Wrapper::OpenFile(file_obj &__file, std::string file_name)
{
#ifdef MPI_ON
    if (!MPI_File_open(world, file_name.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &__file)) {
        if (__file == MPI::FILE_NULL) {
            progIO->error_message << "Error: Failed to open file " << file_name << std::endl;
            progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
            exit(3); // cannot open file
        }
    }
#else // MPI_ON
    __file.open(file_name.c_str(), std::ofstream::out);
    if (!__file.is_open()) {
        progIO->error_message << "Error: Failed to open file " << file_name << std::endl;
        progIO->Output(std::cerr, progIO->error_message, __normal_output, __all_processors);
        exit(3); // cannot open file
    }
#endif // MPI_ON
}

/*! \fn void WriteFile(file_obj &__file, std::ostringstream &content, const int loop_count)
 *  \brief all processor write into a file, use with cautions -> read the assumptions in descriptions
 *  When MPI_ON is on, this function assumes that you only write file header if specifying __master_only, and it assumes that every processor are writing the same amout of chunk into the file every time */
void MPI_Wrapper::WriteSingleFile(file_obj &__file, std::ostringstream &content, const MPI_Level &mpi_level)
{
#ifdef MPI_ON
    std::string tmp_str = content.str();
    if (mpi_level == __master_only) {
        header_offset = tmp_str.size();
        MPI_Bcast(&header_offset, 1, MPI_OFFSET, master, world);
        if (myrank == master) {
            const char *tmp_char = tmp_str.data();
            MPI_File_write(__file, tmp_char, static_cast<int>(tmp_str.size()), MPI::CHAR, &status);
        } else {
            MPI_File_seek(__file, header_offset, MPI_SEEK_SET);
        }
        MPI_File_sync(__file);
    } else {
        if (loop_end == -1) {
            return; // this should not happen
        }
        MPI_File_get_position(__file, &offset);
        if (offset == header_offset) {
            MPI_File_seek(__file, tmp_str.size()*(loop_begin), MPI_SEEK_CUR);
        } else {
            MPI_File_seek(__file, tmp_str.size()*(loop_step-1), MPI_SEEK_CUR);
        }
        const char *tmp_char = tmp_str.data();
        MPI_File_write(__file, tmp_char, static_cast<int>(tmp_str.size()), MPI::CHAR, &status);
        MPI_File_sync(__file);
    }
#else // MPI_ON
    __file << content.str();
#endif // MPI_ON
    progIO->Reset(content);
}

/*! \fn void CloseFile(file_obj &__file, )
 *  \brief close the file for data writing */
void MPI_Wrapper::CloseFile(file_obj &__file)
{
#ifdef MPI_ON
    MPI_File_close(&__file);
#else // MPI_ON
    __file.close();
#endif // MPI_ON
}

/*! \fn ~MPI_info()
 *  \brief a destructor */
MPI_Wrapper::~MPI_Wrapper()
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
    return MPI::Wtime();
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
