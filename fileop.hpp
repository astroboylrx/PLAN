//
//  fileop.hpp
//  PLATO: PLAneTesimal locatOr
//
//  Created by Rixin Li on 3/10/16.
//  Copyright © 2016 Rixin Li. All rights reserved.
//

#ifndef fileop_hpp
#define fileop_hpp

#include "global.hpp"

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

#endif /* fileop_hpp */
