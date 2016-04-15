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

#ifdef MPI_ON
/*! \var MPI_Wrapper *mpi
 *  \brief global wrapper for MPI routines and related variables */
MPI_Wrapper *mpi = new MPI_Wrapper;
#endif // MPI_ON

/*! \var std::vector<Timer> timer
 *  \brief global timer */
std::vector<Timer> timer(__time_type_count);

/*! \var FileOperation *io_ops
 *  \brief global class handle all I/O stuff */
IO_Operations *io_ops = new IO_Operations;

/*! \fn int main(int argc, const char * argv[])
 *  \brief main function. */
int main(int argc, const char * argv[])
{
    /********** First Step: Initialization **********/
#ifdef MPI_ON
    mpi->Initialization(argc, argv);
#endif // MPI_ON
    timer[__total_elapse_time].StartTimer();
    io_ops->log_info << "Program begins now.\n";
    io_ops->Output(std::clog, io_ops->log_info, io_ops->__normal_output, io_ops->__master_only);
    io_ops->PrintStars(std::clog, io_ops->__normal_output);
    io_ops->Initialize(argc, argv);
    
    SmallVec<double, 3> test(2.5, 3.8, 4.1);
    std::clog << test << std::endl;
    
    
    
    
    
    
    
    
    timer[__total_elapse_time].StopTimer();
    io_ops->log_info << "Program ends now. Elapsed time: " << timer[__total_elapse_time].GiveTime() << "\n";
    io_ops->PrintStars(std::clog, io_ops->__normal_output);
    io_ops->Output(std::clog, io_ops->log_info, io_ops->__normal_output, io_ops->__master_only);
    
#ifdef MPI_ON
    mpi->Finalize();
#endif // MPI_ON

    return 0;
}
