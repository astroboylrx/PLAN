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
    progIO->out_content << "Program begins now (local time: " << progIO->LocalTime() << ")." << std::endl;
    progIO->Output(std::cout, progIO->out_content, __normal_output, __master_only);
    progIO->PrintStars(std::cout, __normal_output);
    progIO->Initialize(argc, argv);
    mpi->DetermineLoop(progIO->num_files);
    DataSet<float, dim> ds;
    //sleep(30);
        
    /********** Step II: Pre-loop Work **********/
    BasicAnalysesPreWork();
    
    /********** Step III: Loop files, read data and process it **********/
    for (int loop_count = mpi->loop_begin; loop_count <= mpi->loop_end; loop_count += mpi->loop_step) {
        
        /***** Step III-A, read original data and perform basic analyses *****/
        timer[__tmp_used_timer].StartTimer();
        ds.particle_set.ReadLisFile(loop_count);
        //ds.vtk_data.ReadVtkFile(loop_count);
        timer[__tmp_used_timer].StopTimer();
        progIO->out_content << "Reading loop_count (" << loop_count << ") cost " << timer[__tmp_used_timer].GiveTime() << " seconds\n";
        progIO->Output(std::cout, progIO->out_content, __more_output, __all_processors);

        BasicAnalyses(ds, loop_count);
        
        if (progIO->flags.find_clumps_flag || progIO->flags.density_vs_scale_flag) {
            ds.particle_set.MakeGhostParticles(progIO->numerical_parameters);
            ds.tree.BuildTree(progIO->numerical_parameters, ds.particle_set);

            BasicAnalysesWithTree(ds, loop_count);

            /***** Step III-B, identity high density region and find planetesimals *****/
            ds.tree.FindPlanetesimals(ds, BHtree<dim>::QseudoQuadraticSplinesKernel<float>, loop_count);
            ds.planetesimal_list.WriteBasicResults(loop_count);

        } // if (progIO->flags.find_clumps_flag)
        
    }

    /********** Step IV: Post-loop Work **********/
    BasicAnalysesPostWork();
    mpi->Barrier();
    timer[__waiting_time].StopTimer();
    mpi->Finalize();

    std::flush(std::clog);
    progIO->PrintStars(std::cout, __normal_output);
    progIO->out_content << "Program ends now (local time: " << progIO->LocalTime() << "). Elapsed time: " << timer[__total_elapse_time].GiveTime() << " seconds.\n";
    progIO->Output(std::cout, progIO->out_content, __normal_output, __master_only);

    return 0;
}
