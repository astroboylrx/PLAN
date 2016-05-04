//
//  analyses.hpp
//  PLAN: PLantesimal ANalyzer
//
//  Created by Rixin Li on 5/4/16.
//  Copyright © 2016 Rixin Li. All rights reserved.
//

/*! \file analyses.hpp
 *  \brief provide methods for intriguing analyses */

#ifndef analyses_hpp
#define analyses_hpp

#include "global.hpp"
#include "tree.hpp"

/*! \fn void BasicAnalysesPreWork()
 *  \brief open default result file and write file header */
void BasicAnalysesPreWork();

/*! \fn void template <int dim> BasicAnalyses(ParticleSet<dim> &particle_set, int loop_count)
 *  \brief calculate max($\rho_p$) and $H_p$ */
template <int dim>
void BasicAnalyses(ParticleSet<dim> &particle_set, int loop_count) {
    // N.B. use reference to particle_set to avoid unwanted changes to members in it after calling this function
    
    if (!progIO->flags.basic_analyses_flag) {
        return;
    }
    
    if (loop_count == mpi->loop_begin) {
        progIO->column += particle_set.num_type;
        for (int i = 0; i != particle_set.num_type; i++) {
            std::ostringstream tmp_oss;
            tmp_oss << "H_p[" << i << "]";
            progIO->out_content << std::setw(progIO->width) << std::setfill(' ') << tmp_oss.str();
        }
        progIO->out_content << std::endl;
        mpi->WriteSingleFile(mpi->result_file, progIO->out_content, __master_only);
    }
    
    if (progIO->flags.basic_analyses_flag) {
        // initialize particle scale heights
        progIO->physical_quantities[loop_count].particle_scale_height.resize(particle_set.num_type);
        for (auto &item : progIO->physical_quantities[loop_count].particle_scale_height) {
            item = 0;
        }
        // calculate particle scale heights and find out maximum particle density
        for (__uint32_t i = 0; i != particle_set.num_particle; i++) {
            progIO->physical_quantities[loop_count].particle_scale_height[particle_set.particles[i].property_index] += particle_set.particles[i].pos[2] * particle_set.particles[i].pos[2];
            progIO->physical_quantities[loop_count].max_particle_density = std::max(progIO->physical_quantities[loop_count].max_particle_density, particle_set.particles[i].density);
        }
        progIO->log_info << "loop_count = " << loop_count << ", max(rho_p) = " << progIO->physical_quantities[loop_count].max_particle_density;
        int tmp_index = 0;
        for (auto &item : progIO->physical_quantities[loop_count].particle_scale_height) {
            item = std::sqrt(item/particle_set.num_particle);
            progIO->log_info << ", H_p[" << tmp_index++ << "] = " << item;
        }
        progIO->log_info << std::endl;
        progIO->Output(std::clog, progIO->log_info, __more_output, __all_processors);
    } // if (progIO->flags.basic_analyses_flag)
    
    progIO->out_content << std::setw(progIO->width) << std::setfill(' ') << std::scientific << progIO->physical_quantities[loop_count].time << std::setw(progIO->width) << progIO->physical_quantities[loop_count].max_particle_density;
    for (auto &item : progIO->physical_quantities[loop_count].particle_scale_height) {
        progIO->out_content << std::setw(progIO->width) << std::setfill(' ') << std::scientific << item;
    }
    progIO->out_content << std::endl;
    mpi->WriteSingleFile(mpi->result_file, progIO->out_content, __all_processors);
}

/*! \fn void BasicAnalysesPostWork()
 *  \brief close default result file  */
void BasicAnalysesPostWork();

#endif /* analyses_hpp */
