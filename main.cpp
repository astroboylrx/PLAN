//
//  main.cpp
//  PLAN: PLantesimal ANalyzer
//
//  Created by Rixin Li on 3/8/16.
//  Copyright © 2016 Rixin Li. All rights reserved.
//

/*! \file main.cpp
 *  \brief contains definitions of global/const variables and main workflow */

#include "global.hpp"
#include "tree.hpp"
#include "analyses.hpp"

/***********************************************************/
/********** Definitions of Global/Const Varialbes **********/
/***********************************************************/

/*! \var MPI_Wrapper *mpi
 *  \brief global wrapper for MPI routines and related variables
 *  Note that even if running as a serial program, this still works since MPI_Wrapper takes care of it. The goal of this encapsulation is to reduce the use of (annoying) "#ifdef MPI_ON" in main function. */
MPI_Wrapper *mpi = new MPI_Wrapper;

/*! \var std::vector<Timer> timer
 *  \brief global timer */
std::vector<Timer> timer(__time_type_count);

/*! \var FileOperation *progIO
 *  \brief Handle command line and most of I/O functions together with MPI_Wrapper */
Basic_IO_Operations *progIO = new Basic_IO_Operations;

/*! \var const int dim = 3
 *  \brief dimension of simulation */
const int dim = 3;

/***********************************/
/********** Main function **********/
/***********************************/

/*! \fn int main(int argc, const char * argv[])
 *  \brief main function. */
int main(int argc, const char * argv[])
{
    /********** Step I: Initialization **********/
    mpi->Initialization(argc, argv);
    timer[__total_elapse_time].StartTimer();
    
    progIO->log_info << "Program begins now.\n";
    progIO->Output(std::clog, progIO->log_info, __normal_output, __master_only);
    progIO->PrintStars(std::clog, __normal_output);
    
    progIO->Initialize(argc, argv);
    mpi->DetermineLoop(progIO->num_file);
    
    ParticleSet<dim> particle_set;
    BHtree<dim> tree;
    std::vector<Plantesimal<dim>> plantesimals;
    const SmallVec<double, dim> __center(0.0); // default center
    
    /********** Step II: Pre-loop Work **********/
    BasicAnalysesPreWork();
    
    /********** Step III: Loop files, read data and process it **********/
    for (int loop_count = mpi->loop_begin; loop_count <= mpi->loop_end; loop_count += mpi->loop_step) {
        /***** Step III-A, read original data and perform basic analyses *****/
        
        particle_set.ReadLisFile(loop_count);
        BasicAnalyses(particle_set, loop_count);
        
        particle_set.MakeGhostParticles(0.2);
        //tree.BuildTree(__center, 0.1, particle_set, 8);
        
        /***** Step III-B, identity high density region and find planetesimals *****/
        tree.FindPlanetesimals();
        
        
    }

    
    
    /********** Step IV: Post-loop Work **********/
    BasicAnalysesPostWork();
    
    mpi->Barrier();
    timer[__total_elapse_time].StopTimer();
    progIO->log_info << "Program ends now. Elapsed time: " << timer[__total_elapse_time].GiveTime() << "\n";
    progIO->PrintStars(std::clog, __normal_output);
    progIO->Output(std::clog, progIO->log_info, __normal_output, __master_only);
    
    mpi->Finalize();

    return 0;
}
