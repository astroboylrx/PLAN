//
//  analyses.hpp
//  PLAN: PLanetesimal ANalyzer
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

/*! \fn template <int D> void BasicAnalyses(ParticleSet<dim> &particle_set, BHtree<D> &tree, int loop_count)
 *  \brief calculate (but not limited to) max($\rho_p$) and $H_p$ */
template <int D>
void BasicAnalyses(ParticleSet<D> &particle_set, BHtree<D> &tree, int loop_count)
{
    // N.B. use reference to particle_set to avoid unwanted changes to members in it after calling this function
    
    if (progIO->flags.density_vs_scale_flag) {
        progIO->numerical_parameters.ghost_zone_width = progIO->numerical_parameters.max_half_width;
        tree.max_leaf_size = pow(1<<D, 2);
        progIO->flags.no_ghost_particle_flag = 1;
    }
    
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
        mpi->WriteSingleFile(*mpi->files[progIO->file_name.output_file_path], progIO->out_content, __master_only);
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
    mpi->WriteSingleFile(*mpi->files[progIO->file_name.output_file_path], progIO->out_content, __all_processors);
}

/*! \struct template <int D, int d, class Compare, typename std::enable_if<bool(D>d), int>::type = 0> struct CompareParticle 
 *  \brief special object used to compare two Particle objects */
template <int D, int d, class Compare, typename std::enable_if<bool(D>d), int>::type = 0>
struct CompareParticle {
    Compare comp;
    bool operator() (const Particle<D> &lhs, const Particle<D> &rhs) {
        return comp(lhs.pos[d], rhs.pos[d]);
    }
    // usage: define an object then use it
    // CompareParticle<dim, xyz, std::less<double>> partilce_less;
};

/*! \fn template <int D, int d, class Compare> void RecursiveFindExtreme(BHtree<D> &tree, __uint32_t node, __uint32_t &extreme, Compare comp)
 *  \brief recursive find the extreme value of particles' x/y/z coordinates
 *  \nparam D dimension of data
 *  \nparam d which direction to compare, x/y/z
 *  \nparam Compare which extreme to obtain, give std::less for min value, give std::greater for max value */
template <int D, int d, class Compare>
void RecursiveFindExtreme(BHtree<D> &tree, __uint32_t node, __uint32_t &extreme, Compare comp)
{
    if (tree.IsLeaf(node)) { // if leaf, just compare all particles
        for (__uint32_t p = tree.tree[node].begin; p != tree.tree[node].end; p++) {
            if (comp(tree.particle_list[p].pos[d], tree.particle_list[extreme].pos[d])) {
                extreme = p;
            }
        }
    } else { // else, enter daughter nodes
        for (__uint32_t daughter = tree.tree[node].first_daughter; daughter != tree.tree[node].first_daughter + tree.tree[node].num_daughter; daughter++) {
            if (comp(tree.tree[daughter].center[d], tree.tree[node].center[d])) { // if daughter is "comp" than parent, enter it
                RecursiveFindExtreme<D, d, Compare>(tree, daughter, extreme, comp);
            } else { // else, only if parent is "comp" than current extreme, do we enter it
                if (comp(tree.tree[node].center[d], tree.particle_list[extreme].pos[d])) {
                    RecursiveFindExtreme<D, d, Compare>(tree, daughter, extreme, comp);
                }
            }
        }
    }
}

/*! \fn template <int D> void BasicAnalysesWithTree(ParticleSet<dim> &particle_set, BHtree<D> &tree, int loop_count)
 *  \brief calculate (but not limited to) max($\rho_p$) as a function of length scale */
template <int D>
void BasicAnalysesWithTree(ParticleSet<D> &particle_set, BHtree<D> &tree, int loop_count)
{
    if (progIO->flags.density_vs_scale_flag) {
        double z_min, z_max;
        if (D == 3) { // if in 3D, find the z_max and z_min among all particles
            timer[__tmp_used_timer].StartTimer();
            
            /* first method,
            Particle<dim> *p_min, *p_max;
            CompareParticle<dim, 2, std::less<double>> particle_less;
            p_min = std::min_element(particle_set.particles, particle_set.particles+particle_set.num_particle, particle_less);
            p_max = std::max_element(particle_set.particles, particle_set.particles+particle_set.num_particle, particle_less);
            unsigned int t_1st_index = timer[__tmp_used_timer].Lap();
            // second method
            double pos2_min = 0.0, pos2_max = 0.0;
            for (__uint32_t i = 0; i != particle_set.num_particle; i++) {
                if (particle_set[i].pos[2] < pos2_min) {
                    pos2_min = particle_set[i].pos[2];
                }
                if (particle_set[i].pos[2] > pos2_max) {
                    pos2_max = particle_set[i].pos[2];
                }
            }
            unsigned int t_2nd_index = timer[__tmp_used_timer].Lap(); // */
            // third method
            __uint32_t extreme_min = 0, extreme_max = 0;
            RecursiveFindExtreme<dim, 2, std::less<double>>(tree, tree.root, extreme_min, std::less<double>());
            RecursiveFindExtreme<dim, 2, std::greater<double>>(tree, tree.root, extreme_max, std::greater<double>());
            unsigned int t_3rd_index = timer[__tmp_used_timer].Lap();
            timer[__tmp_used_timer].StopTimer();
            
            //progIO->log_info << "First method: z_min = " << p_min->pos[2] << ", z_max = " << p_max->pos[2] << ", time = " << timer[__tmp_used_timer].GiveTime(t_1st_index) << std::endl;
            //progIO->log_info << "Second method: z_min = " << pos2_min << ", z_max = " << pos2_max << ", time = " << timer[__tmp_used_timer].GiveTime(t_2nd_index) << std::endl;
            progIO->log_info << "Third method: z_min = " << tree.particle_list[extreme_min].pos[2] << ", z_max = " << tree.particle_list[extreme_max].pos[2] << ", time = " << timer[__tmp_used_timer].GiveTime(t_3rd_index) << std::endl;
            progIO->Output(std::clog, progIO->log_info, __even_more_output, __master_only);
            
            // force to use cell center
            z_min = progIO->numerical_parameters.box_center[2] - progIO->numerical_parameters.cell_length[2] / 2.0;
            z_max = tree.particle_list[extreme_max].pos[2]; // no need to restrict large end
            while (z_min - progIO->numerical_parameters.cell_length[2] > tree.particle_list[extreme_min].pos[2]) {
                z_min -= progIO->numerical_parameters.cell_length[2];
            }
        } else {
            exit(5); // haven't implement yet
        }
        
        // Calculate how many length scales that we need to calculate, specifically assume that max radius is determined by horizontal box size
        int num_level = round(log10(progIO->numerical_parameters.box_half_width[0]*2/progIO->numerical_parameters.cell_length[0])/log10(2.0)) + 1; // RL: assume the cell has equal size in each dimension
        progIO->physical_quantities[loop_count].max_rhop_vs_scale = std::vector<double>(num_level, 0);
        
        std::vector<double> search_radii(num_level, progIO->numerical_parameters.cell_length[0]/2.0);
        std::vector<SmallVec<double, D>> max_distance(num_level, SmallVec<double, D>(0)); // maximum distance for not applying shear
        double shear_distance = progIO->numerical_parameters.shear_speed * particle_set.time;
        for (unsigned int i = 0; i != search_radii.size(); i++) {
            search_radii[i] *= pow(2.0, i);
            max_distance[i] = progIO->numerical_parameters.box_length - SmallVec<double, D>(search_radii[i]);
            if (D == 3) {
                max_distance[i][2] = 100 * progIO->numerical_parameters.max_half_width; // force to disable periodic search in vertical direction
            }
            for (int j = 0; j != D; j++) {
                assert(max_distance[i][j] > progIO->numerical_parameters.box_half_width[j]);
                // For now, our method only works if max_distance is larger than box_half_width in every dimension
            }
        }
        
        if (loop_count == mpi->loop_begin) {
            progIO->out_content << std::setw(progIO->width) << std::setfill(' ') << "#time";
            for (unsigned int i = 0; i != search_radii.size(); i++) {
                progIO->out_content << std::setw(progIO->width) << std::setfill(' ') << std::scientific << search_radii[i];
                progIO->log_info << "search_radii[" << i << "] = " << search_radii[i] << ", max_distance[" << i << "] = " << max_distance[i] << std::endl;
            }
            progIO->out_content << std::endl;
            mpi->WriteSingleFile(*mpi->files[progIO->file_name.max_rhop_vs_scale_file], progIO->out_content, __master_only);
            progIO->log_info << std::endl;
            progIO->Output(std::clog, progIO->log_info, __even_more_output, __master_only);
            assert(particle_set.num_type == 1);
        }
        progIO->physical_quantities[loop_count].mass_per_particle[0] = progIO->physical_quantities[loop_count].solid_to_gas_ratio * progIO->numerical_parameters.box_length[0] * progIO->numerical_parameters.box_length[1] * std::sqrt(2*3.1415926535897932) / particle_set.num_particle;
        
        // Traverse all the effective cell centers
        SmallVec<double, D> cell_center_min = progIO->numerical_parameters.box_min + 0.5 * progIO->numerical_parameters.cell_length;
        SmallVec<double, D> cell_center_max = progIO->numerical_parameters.box_max;
        SmallVec<double, D> center_step = progIO->numerical_parameters.cell_length, center_substep = center_step / 2.0, center_subsubstep = center_step / 4.0;
        SmallVec<double, D> current_center(0), subcenter(0), subsubcenter(0);
        __uint32_t tmp_count = 0, point_count = 0;
        timer[__tmp_used_timer].StartTimer();
        for (double tmp_z = z_min; tmp_z < z_max; tmp_z += center_step[2]) {
            for (double tmp_y = cell_center_min[1]; tmp_y < cell_center_max[1]; tmp_y += center_step[1]) {
                for (double tmp_x = cell_center_min[0]; tmp_x < cell_center_max[0]; tmp_x += center_step[0]) {
                    point_count++;
                    current_center = SmallVec<double, D>(tmp_x, tmp_y, tmp_z);
                    // loop for all the radii
                    for (unsigned int i = 0; i != search_radii.size(); i++) {
                        tmp_count = 0;
                        //tree.RecursiveBallSearchCount(current_center, tree.root, search_radii[i], tmp_count);
                        tree.RecursiveBallSearchCountWithShear(current_center, tree.root, search_radii[i], tmp_count, max_distance[i], shear_distance);
                        if (tmp_count > progIO->physical_quantities[loop_count].max_rhop_vs_scale[i]) {
                            progIO->physical_quantities[loop_count].max_rhop_vs_scale[i] = tmp_count;
                        }
                        // small radii need more sampling
                        if (i < 3) {
                            for (int orthant = 0; orthant != 1<<D; orthant++) {
                                tmp_count = 0;
                                subcenter = current_center + center_substep.ParaMultiply(Orthant<D>::orthants[orthant]);
                                //tree.RecursiveBallSearchCount(subcenter, tree.root, search_radii[i], tmp_count);
                                tree.RecursiveBallSearchCountWithShear(subcenter, tree.root, search_radii[i], tmp_count, max_distance[i], shear_distance);
                                if (tmp_count > progIO->physical_quantities[loop_count].max_rhop_vs_scale[i]) {
                                    progIO->physical_quantities[loop_count].max_rhop_vs_scale[i] = tmp_count;
                                }
                                if (i == 0) { // the smallest radii needs most sampling
                                    for (int suborthant = 0; suborthant != 1<<D; suborthant++) {
                                        tmp_count = 0;
                                        subsubcenter = subcenter + center_subsubstep.ParaMultiply(Orthant<D>::orthants[orthant]);
                                        //tree.RecursiveBallSearchCount(subsubcenter, tree.root, search_radii[i], tmp_count);
                                        tree.RecursiveBallSearchCountWithShear(subsubcenter, tree.root, search_radii[i], tmp_count, max_distance[i], shear_distance);
                                        if (tmp_count > progIO->physical_quantities[loop_count].max_rhop_vs_scale[i]) {
                                            progIO->physical_quantities[loop_count].max_rhop_vs_scale[i] = tmp_count;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } // RL: this is really not elegent
        timer[__tmp_used_timer].StopTimer();
        progIO->log_info << "loop_count = " << loop_count << ", BallSearch cost " << timer[__tmp_used_timer].GiveTime() << ", points searched = " << point_count << std::endl;
        progIO->Output(std::clog, progIO->log_info, __even_more_output, __all_processors);
        
        progIO->out_content << std::setw(progIO->width) << std::setfill(' ') << std::scientific << progIO->physical_quantities[loop_count].time;
        for (unsigned int i = 0; i != search_radii.size(); i++) {
            progIO->physical_quantities[loop_count].max_rhop_vs_scale[i] *= (progIO->physical_quantities[loop_count].mass_per_particle[0] / (4.1887902047863910 * pow(search_radii[i], 3)));
            progIO->out_content << std::setw(progIO->width) << std::setfill(' ') << std::scientific << progIO->physical_quantities[loop_count].max_rhop_vs_scale[i];
        }
        progIO->out_content << std::endl;
        mpi->WriteSingleFile(*mpi->files[progIO->file_name.max_rhop_vs_scale_file], progIO->out_content, __all_processors);
    }
}


/*! \fn void BasicAnalysesPostWork()
 *  \brief close default result file  */
void BasicAnalysesPostWork();

/*! \fn template <int D> void MinDistanceBetweenParticles(ParticleSet<dim> &particle_set, 
 *  \brief find the minimum distance between particles in data */
template <int D>
void MinDistanceBetweenParticles(ParticleSet<D> &particle_set, BHtree<D> &tree, int loop_count)
{
    if (!progIO->flags.basic_analyses_flag) {
        return;
    }
    
    double min_distance = progIO->numerical_parameters.max_half_width, tmp_min_distance = 0;
    int indices[3];
    for (__uint32_t i = 0; i != particle_set.num_particle; i++) {
        tree.KNN_Search(particle_set[i].pos, 2, tmp_min_distance, indices); // the closest one is itself
        min_distance = std::min(min_distance, tmp_min_distance);
    }
    progIO->log_info << "loop_count = " << loop_count << ", min_distance_between_particles in data = " << min_distance << std::endl;
    progIO->Output(std::clog, progIO->log_info, __more_output, __all_processors);
}

#endif /* analyses_hpp */
