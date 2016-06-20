//
//  main.cpp
//  PLAN: PLanetesimal ANalyzer
//
//  Created by Rixin Li on 3/8/16.
//  Copyright © 2016 Rixin Li. All rights reserved.
//

/*! \file main.cpp
 *  \brief contains definitions of global/const variables and main workflow */

#include "global.hpp"
#include "tree.hpp"
#include "analyses.hpp"

/*****************************************************/
/********** Definitions of Global Varialbes **********/
/*****************************************************/

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
    progIO->Initialize(argc, argv);
    progIO->log_info << "Program begins now.\n";
    progIO->Output(std::clog, progIO->log_info, __normal_output, __master_only);
    progIO->PrintStars(std::clog, __normal_output);
    
    mpi->DetermineLoop(progIO->num_file);
    ParticleSet<dim> particle_set;
    BHtree<dim> tree;
    std::vector<Planetesimal<dim>> planetesimals;
    
    /********** Step II: Pre-loop Work **********/
    BasicAnalysesPreWork();
    
    /********** Step III: Loop files, read data and process it **********/
    for (int loop_count = mpi->loop_begin; loop_count <= mpi->loop_end; loop_count += mpi->loop_step) {
        
        /***** Step III-A, read original data and perform basic analyses *****/
        particle_set.ReadLisFile(loop_count);
        BasicAnalyses(particle_set, tree, loop_count);
        
        if (progIO->flags.find_clumps_flag || progIO->flags.density_vs_scale_flag) {
            particle_set.MakeGhostParticles(progIO->numerical_parameters);
            tree.BuildTree(progIO->numerical_parameters, particle_set);
            //tree.max_leaf_size = static_cast<unsigned int>(tree.max_leaf_size*(tree.num_nodes/(particle_set.num_total_particle/15.0)));
            //progIO->log_info << "Now change tree.max_leaf_size to " << tree.max_leaf_size << std::endl;
            //progIO->Output(std::clog, progIO->log_info, __even_more_output, __all_processors);
            //tree.BuildTree(progIO->numerical_parameters, particle_set);
            
            tree.CheckTree(tree.root, tree.root_level, tree.root_center, tree.half_width);
            BasicAnalysesWithTree(particle_set, tree, loop_count);

            /***** Step III-B, identity high density region and find planetesimals *****/
            tree.FindPlanetesimals();
            
        } // if (progIO->flags.find_clumps_flag)
        
    }

    /********** Step IV: Post-loop Work **********/
    BasicAnalysesPostWork();
    mpi->Barrier();
    timer[__waiting_time].StopTimer();
    
    progIO->log_info << "Program ends now. Elapsed time: " << timer[__total_elapse_time].GiveTime() << "\n";
    progIO->PrintStars(std::clog, __normal_output);
    progIO->Output(std::clog, progIO->log_info, __normal_output, __master_only);
    
    mpi->Finalize();

    return 0;
}
