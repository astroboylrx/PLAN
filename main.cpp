//
//  main.cpp
//  PLATO: PLAneTesimal locatOr
//
//  Created by Rixin Li on 3/8/16.
//  Copyright © 2016 Rixin Li. All rights reserved.
//

#include "global.hpp"
#include "tree.hpp"

// Definitions of Global Varialbes

/*! \var MPI_Wrapper *mpi
 *  \brief global wrapper for MPI routines and related variables
 *  Note that even if running as a serial program, this still works since MPI_Wrapper takes care of it. The goal of this encapsulation is to reduce the use of (annoying) "#ifdef MPI_ON" in main function. */
MPI_Wrapper *mpi = new MPI_Wrapper;

/*! \var std::vector<Timer> timer
 *  \brief global timer */
std::vector<Timer> timer(__time_type_count);

/*! \var FileOperation *io_ops
 *  \brief global class handle all I/O stuff */
IO_Operations *io_ops = new IO_Operations;

/*! \var const int dim = 3
 *  \brief dimension of simulation */
const int dim = 3;

/*! \fn int main(int argc, const char * argv[])
 *  \brief main function. */
int main(int argc, const char * argv[])
{
    /********** Step I: Initialization **********/
    mpi->Initialization(argc, argv);
    timer[__total_elapse_time].StartTimer();
    
    io_ops->log_info << "Program begins now.\n";
    io_ops->Output(std::clog, io_ops->log_info, io_ops->__normal_output, io_ops->__master_only);
    io_ops->PrintStars(std::clog, io_ops->__normal_output);
    
    io_ops->Initialize(argc, argv);
    mpi->DetermineLoop(io_ops->num_file);
    
    ParticleSet<dim> particle_set;
    BHtree<dim> tree;
    
    /********** Step II: Pre-loop Work **********/
    std::vector<std::string>::iterator file_head = io_ops->file_name.lis_data_file_name.begin();
    
    /********** Step III: Loop files, read data and process it **********/
    for (int loop_count = mpi->loop_begin; loop_count <= mpi->loop_end; loop_count += mpi->loop_step) {
        /***** Step III-A, read data from lis files *****/
        
        particle_set.ReadLisFile(file_head + loop_count * io_ops->num_cpu,
                                 file_head + loop_count * io_ops->num_cpu + io_ops->num_cpu);
        
    }

    
    
    
    
    io_ops->PrintStars(std::clog, io_ops->__normal_output);
    /********** Step IV: Post-loop Work **********/
    timer[__total_elapse_time].StopTimer();
    io_ops->log_info << "Program ends now. Elapsed time: " << timer[__total_elapse_time].GiveTime() << "\n";
    io_ops->PrintStars(std::clog, io_ops->__normal_output);
    io_ops->Output(std::clog, io_ops->log_info, io_ops->__normal_output, io_ops->__master_only);
    
    mpi->Finalize();

    return 0;
}
