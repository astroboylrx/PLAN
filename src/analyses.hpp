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

/*! \fn template <int D> void BasicAnalyses(DataSet<T, D> &data, int loop_count)
 *  \brief calculate (but not limited to) max($\rho_p$) and $H_p$ */
template <class T, int D>
void BasicAnalyses(DataSet<T, D> &ds, int loop_count)
{
    // N.B. use reference to particle_set to avoid unwanted changes to members in it after calling this function
    
    if (progIO->flags.density_vs_scale_flag) {
        progIO->numerical_parameters.ghost_zone_width = progIO->numerical_parameters.max_half_width;
        ds.tree.max_leaf_size = std::pow(1<<D, 2);
        progIO->flags.no_ghost_particle_flag = 1;
    }

    if (progIO->flags.basic_analyses_flag || progIO->flags.find_clumps_flag) {
        // initialize particle scale heights
        progIO->physical_quantities[loop_count].particle_scale_height.resize(ds.particle_set.num_types);
        for (auto &item : progIO->physical_quantities[loop_count].particle_scale_height) {
            item = 0;
        }

        // calculate particle scale heights and find out maximum particle density
        for (uint32_t i = 0; i != ds.particle_set.num_particles; i++) {
            progIO->physical_quantities[loop_count].particle_scale_height[ds.particle_set.particles[i].property_index] += ds.particle_set.particles[i].pos[2] * ds.particle_set.particles[i].pos[2];
            progIO->physical_quantities[loop_count].max_particle_density = std::max(progIO->physical_quantities[loop_count].max_particle_density, ds.particle_set.particles[i].density);
        }
    }

    if (progIO->flags.basic_analyses_flag) {
        if (loop_count == mpi->loop_begin) {
            // write header for result.txt
            progIO->out_content << std::setw(progIO->width) << std::setfill(' ') << "#time" << std::setw(progIO->width) << "max(rho_p)";
            progIO->column += ds.particle_set.num_types;
            for (unsigned int i = 0; i != ds.particle_set.num_types; i++) {
                std::ostringstream tmp_oss;
                tmp_oss << "H_p[" << i << "]";
                progIO->out_content << std::setw(progIO->width) << std::setfill(' ') << tmp_oss.str();
            }
            if (D == 3) {
                progIO->column += 1; // for vertical flux
                progIO->out_content << std::setw(progIO->width) << std::setfill(' ') << "dSigma" << std::endl;
            }
            mpi->WriteSingleFile(mpi->result_files[mpi->file_pos[progIO->file_name.output_file_path]], progIO->out_content, __master_only);
            
            // write header for MeanSigma.txt
            if (D == 3) {
                progIO->out_content << std::setw(progIO->width) << std::setfill(' ') << 0.0;
                typename VtkDataVector<T, D>::template view_r2d_type<> x_center(ds.vtk_data.cell_center[boost::indices[0][0][sn::b_range()][sn::b_range()]]);
                progIO->out_content << std::scientific;
                int count = 0;
                for (auto item : x_center) {
                    progIO->out_content << std::setw(progIO->width) << std::setfill(' ') << item[0];
                    count++;

                }
                std::cout << count << std::endl;
                for (auto item : x_center) {
                    progIO->out_content << std::setw(progIO->width) << std::setfill(' ') << item[0];
                }
                progIO->out_content << std::endl;
                mpi->WriteSingleFile(mpi->result_files[mpi->file_pos[progIO->file_name.mean_sigma_file]], progIO->out_content, __master_only);
            }
        } // if (loop_count == mpi->loop_begin)

        progIO->log_info << "loop_count = " << loop_count << ", max(rho_p) = " << progIO->physical_quantities[loop_count].max_particle_density;
        int tmp_index = 0;
        for (auto &item : progIO->physical_quantities[loop_count].particle_scale_height) {
            item = std::sqrt(item / ds.particle_set.num_particles);
            progIO->log_info << ", H_p[" << tmp_index++ << "] = " << item;
        }

        // calculate the vertical gas mass flux per unit area (assume 3D cases), <| rho_g c_s |>
        if (D == 3) {
            typename VtkDataVector<T, D>::template view_r1d_type<> slice_z_momentum(ds.vtk_data.vector_data["momentum"].data[boost::indices[sn::b_range().stride(ds.vtk_data.vector_data["momentum"].num_cells[2]-1)][sn::b_range()][sn::b_range()][2]]); // sn::b_range() take all elements in the dimension
            ds.vtk_data.IterateBoostMultiArrayConcept(slice_z_momentum, [](float &element, double &vertical_flux)->void {
                vertical_flux += std::abs(element);
            }, progIO->physical_quantities[loop_count].vertical_flux);
            progIO->physical_quantities[loop_count].vertical_flux /= slice_z_momentum.num_elements();
            progIO->log_info << ", dSigma = " << progIO->physical_quantities[loop_count].vertical_flux;
        }

        progIO->log_info << std::endl;
        progIO->Output(std::clog, progIO->log_info, __more_output, __all_processors);

        // write data to result.txt
        progIO->out_content << std::setw(progIO->width) << std::setfill(' ') << std::scientific << progIO->physical_quantities[loop_count].time << std::setw(progIO->width) << progIO->physical_quantities[loop_count].max_particle_density;
        for (auto &item : progIO->physical_quantities[loop_count].particle_scale_height) {
            progIO->out_content << std::setw(progIO->width) << std::setfill(' ') << item;
        }
        if (D == 3) {
            progIO->out_content << std::setw(progIO->width) << std::setfill(' ')
                                << progIO->physical_quantities[loop_count].vertical_flux;
        }
        progIO->out_content << std::endl;
        mpi->WriteSingleFile(mpi->result_files[mpi->file_pos[progIO->file_name.output_file_path]], progIO->out_content, __all_processors);

        // now calculate <Sigma_g>_{yz}(t) and <Sigma_p>_{yz}(t)
        if (D == 3) {
            double mean_sigma_0 = 0.0;
            unsigned int num_cell_x = ds.vtk_data.scalar_data["density"].num_cells[0];
            progIO->physical_quantities[loop_count].mean_sigma.resize(num_cell_x*2, 0.0);
            for (unsigned int ix = 0; ix != num_cell_x; ix++) {
                typename VtkDataScalar<T, D>::template view_r1d_type<> slice_density(ds.vtk_data.scalar_data["density"].data[boost::indices[sn::b_range()][sn::b_range()][ix]]);
                ds.vtk_data.IterateBoostMultiArrayConcept(slice_density, [](float &element, double &mean_sigma_ix, double &mean_sigma_0)->void {
                    mean_sigma_ix += element;
                    mean_sigma_0 += element;
                }, progIO->physical_quantities[loop_count].mean_sigma[ix], mean_sigma_0);


                typename VtkDataScalar<T, D>::template view_r1d_type<> slice_particle_density(ds.vtk_data.scalar_data["particle_density"].data[boost::indices[sn::b_range()][sn::b_range()][ix]]);
                ds.vtk_data.IterateBoostMultiArrayConcept(slice_particle_density, [](float &element, double &mean_sigma_ix) ->void {
                    mean_sigma_ix += element;
                }, progIO->physical_quantities[loop_count].mean_sigma[ix+num_cell_x]);
            }
            mean_sigma_0 /= num_cell_x;

            for (unsigned int ix = 0; ix != num_cell_x; ix++) {
                progIO->physical_quantities[loop_count].mean_sigma[ix] /= mean_sigma_0;
                // note the following line is different with the previous readvtklis code
                progIO->physical_quantities[loop_count].mean_sigma[ix+num_cell_x] *= (ds.vtk_data.spacing[2]/mean_sigma_0);
            }

            progIO->out_content <<  std::setw(progIO->width) << std::setfill(' ') << std::scientific << progIO->physical_quantities[loop_count].time;
            for (auto &item : progIO->physical_quantities[loop_count].mean_sigma) {
                progIO->out_content << std::setw(progIO->width) << std::setfill(' ') << item;
            }
            progIO->out_content << std::endl;
            mpi->WriteSingleFile(mpi->result_files[mpi->file_pos[progIO->file_name.mean_sigma_file]], progIO->out_content, __all_processors);
        }
    } // if (progIO->flags.basic_analyses_flag)

    if (progIO->flags.tmp_calculation_flag) {

    }

}

/*! \functor template <int D, int d, class Compare, typename std::enable_if<bool(D>d), int>::type = 0> struct CompareParticle 
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

/*! \fn template <int D, int d, class Compare> void RecursiveFindExtreme(BHtree<D> &tree, uint32_t node, uint32_t &extreme, Compare comp)
 *  \brief recursive find the extreme value of particles' x/y/z coordinates
 *  \nparam D dimension of data
 *  \nparam d which direction to compare, x/y/z
 *  \nparam Compare which extreme to obtain, give std::less for min value, give std::greater for max value */
template <int D, int d, class Compare>
void RecursiveFindExtreme(BHtree<D> &tree, uint32_t node, uint32_t &extreme, Compare comp)
{
    if (tree.IsLeaf(node)) { // if leaf, just compare all particles
        for (uint32_t p = tree.tree[node].begin; p != tree.tree[node].end; p++) {
            if (comp(tree.particle_list[p].pos[d], tree.particle_list[extreme].pos[d])) {
                extreme = p;
            }
        }
    } else { // else, enter daughter nodes
        for (uint32_t daughter = tree.tree[node].first_daughter; daughter != tree.tree[node].first_daughter + tree.tree[node].num_daughters; daughter++) {
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

/*! \fn template <int D> void BasicAnalysesWithTree(DataSet<T, D> &ds, int loop_count)
 *  \brief calculate (but not limited to) max($\rho_p$) as a function of length scale */
template <class T, int D>
void BasicAnalysesWithTree(DataSet<T, D> &ds, int loop_count)
{
    if (progIO->flags.density_vs_scale_flag) {
        double z_min, z_max;
        if (D == 3) { // if in 3D, find the z_max and z_min among all particles
            timer[__tmp_used_timer].StartTimer();
            
            /* first method,
            Particle<dim> *p_min, *p_max;
            CompareParticle<dim, 2, std::less<double>> particle_less;
            p_min = std::min_element(ds.particle_set.particles, ds.particle_set.particles+ds.particle_set.num_particles, particle_less);
            p_max = std::max_element(ds.particle_set.particles, ds.particle_set.particles+ds.particle_set.num_particles, particle_less);
            unsigned int t_1st_index = timer[__tmp_used_timer].Lap();
            // second method
            double pos2_min = 0.0, pos2_max = 0.0;
            for (uint32_t i = 0; i != ds.particle_set.num_particles; i++) {
                if (ds.particle_set[i].pos[2] < pos2_min) {
                    pos2_min = ds.particle_set[i].pos[2];
                }
                if (ds.particle_set[i].pos[2] > pos2_max) {
                    pos2_max = ds.particle_set[i].pos[2];
                }
            }
            unsigned int t_2nd_index = timer[__tmp_used_timer].Lap(); // */
            // third method
            uint32_t extreme_min = 0, extreme_max = 0;
            RecursiveFindExtreme<dim, 2, std::less<double>>(ds.tree, ds.tree.root, extreme_min, std::less<double>());
            RecursiveFindExtreme<dim, 2, std::greater<double>>(ds.tree, ds.tree.root, extreme_max, std::greater<double>());
            unsigned int t_3rd_index = timer[__tmp_used_timer].Lap();
            timer[__tmp_used_timer].StopTimer();
            
            //progIO->log_info << "First method: z_min = " << p_min->pos[2] << ", z_max = " << p_max->pos[2] << ", time = " << timer[__tmp_used_timer].GiveTime(t_1st_index) << std::endl;
            //progIO->log_info << "Second method: z_min = " << pos2_min << ", z_max = " << pos2_max << ", time = " << timer[__tmp_used_timer].GiveTime(t_2nd_index) << std::endl;
            progIO->log_info << "Third method: z_min = " << ds.tree.particle_list[extreme_min].pos[2] << ", z_max = " << ds.tree.particle_list[extreme_max].pos[2] << ", time = " << timer[__tmp_used_timer].GiveTime(t_3rd_index) << std::endl;
            progIO->Output(std::clog, progIO->log_info, __even_more_output, __master_only);
            
            // force to use cell center
            z_min = progIO->numerical_parameters.box_center[2] - progIO->numerical_parameters.cell_length[2] / 2.0;
            z_max = ds.tree.particle_list[extreme_max].pos[2]; // no need to restrict large end
            while (z_min - progIO->numerical_parameters.cell_length[2] > ds.tree.particle_list[extreme_min].pos[2]) {
                z_min -= progIO->numerical_parameters.cell_length[2];
            }
        } else {
            exit(5); // haven't implement yet
        }
        
        // Calculate how many length scales that we need to calculate, specifically assume that max radius is determined by horizontal box size
        int num_levels = round(log10(progIO->numerical_parameters.box_half_width[0]*2/progIO->numerical_parameters.cell_length[0])/log10(2.0)) + 1; // RL: assume the cell has equal size in each dimension
        progIO->physical_quantities[loop_count].max_rhop_vs_scale = std::vector<double>(num_levels, 0);
        
        std::vector<double> search_radii(num_levels, progIO->numerical_parameters.cell_length[0]/2.0);
        std::vector<SmallVec<double, D>> max_distance(num_levels, SmallVec<double, D>(0)); // maximum distance for not applying shear
        double shear_distance = progIO->numerical_parameters.shear_speed * ds.particle_set.time;
        for (unsigned int i = 0; i != search_radii.size(); i++) {
            search_radii[i] *= std::pow(2.0, i);
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
            mpi->WriteSingleFile(mpi->result_files[mpi->file_pos[progIO->file_name.max_rhop_vs_scale_file]], progIO->out_content, __master_only);
            progIO->log_info << std::endl;
            progIO->Output(std::clog, progIO->log_info, __even_more_output, __master_only);
            assert(ds.particle_set.num_types == 1);
        }

        // Traverse all the effective cell centers
        SmallVec<double, D> cell_center_min = progIO->numerical_parameters.box_min + 0.5 * progIO->numerical_parameters.cell_length;
        SmallVec<double, D> cell_center_max = progIO->numerical_parameters.box_max;
        SmallVec<double, D> center_step = progIO->numerical_parameters.cell_length, center_substep = center_step / 2.0, center_subsubstep = center_step / 4.0;
        SmallVec<double, D> current_center(0), subcenter(0), subsubcenter(0);
        uint32_t tmp_count = 0, point_count = 0;
        timer[__tmp_used_timer].StartTimer();
        for (double tmp_z = z_min; tmp_z < z_max; tmp_z += center_step[2]) {
            for (double tmp_y = cell_center_min[1]; tmp_y < cell_center_max[1]; tmp_y += center_step[1]) {
                for (double tmp_x = cell_center_min[0]; tmp_x < cell_center_max[0]; tmp_x += center_step[0]) {
                    point_count++;
                    current_center = SmallVec<double, D>(tmp_x, tmp_y, tmp_z);
                    // loop for all the radii
                    for (unsigned int i = 0; i != search_radii.size(); i++) {
                        tmp_count = 0;
                        //ds.tree.RecursiveBallSearchCount(current_center, ds.tree.root, search_radii[i], tmp_count);
                        ds.tree.RecursiveBallSearchCountWithShear(current_center, ds.tree.root, search_radii[i], tmp_count, max_distance[i], shear_distance);
                        if (tmp_count > progIO->physical_quantities[loop_count].max_rhop_vs_scale[i]) {
                            progIO->physical_quantities[loop_count].max_rhop_vs_scale[i] = tmp_count;
                        }
                        // small radii need more sampling
                        if (i < 3) {
                            for (int orthant = 0; orthant != 1<<D; orthant++) {
                                tmp_count = 0;
                                subcenter = current_center + center_substep.ParaMultiply(Orthant<D>::orthants[orthant]);
                                //ds.tree.RecursiveBallSearchCount(subcenter, ds.tree.root, search_radii[i], tmp_count);
                                ds.tree.RecursiveBallSearchCountWithShear(subcenter, ds.tree.root, search_radii[i], tmp_count, max_distance[i], shear_distance);
                                if (tmp_count > progIO->physical_quantities[loop_count].max_rhop_vs_scale[i]) {
                                    progIO->physical_quantities[loop_count].max_rhop_vs_scale[i] = tmp_count;
                                }
                                if (i == 0) { // the smallest radii needs most sampling
                                    for (int suborthant = 0; suborthant != 1<<D; suborthant++) {
                                        tmp_count = 0;
                                        subsubcenter = subcenter + center_subsubstep.ParaMultiply(Orthant<D>::orthants[orthant]);
                                        //ds.tree.RecursiveBallSearchCount(subsubcenter, ds.tree.root, search_radii[i], tmp_count);
                                        ds.tree.RecursiveBallSearchCountWithShear(subsubcenter, ds.tree.root, search_radii[i], tmp_count, max_distance[i], shear_distance);
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
        } // RL: this is really NOT elegant
        timer[__tmp_used_timer].StopTimer();
        progIO->log_info << "loop_count = " << loop_count << ", BallSearch cost " << timer[__tmp_used_timer].GiveTime() << ", points searched = " << point_count << std::endl;
        progIO->Output(std::clog, progIO->log_info, __even_more_output, __all_processors);
        
        progIO->out_content << std::setw(progIO->width) << std::setfill(' ') << std::scientific << progIO->physical_quantities[loop_count].time;
        for (unsigned int i = 0; i != search_radii.size(); i++) {
            progIO->physical_quantities[loop_count].max_rhop_vs_scale[i] *= (progIO->numerical_parameters.mass_per_particle[0] / (4.1887902047863910 * std::pow(search_radii[i], 3))); // 4*PI/3 = 4.1887902047863910
            progIO->out_content << std::setw(progIO->width) << std::setfill(' ') << std::scientific << progIO->physical_quantities[loop_count].max_rhop_vs_scale[i];
        }
        progIO->out_content << std::endl;
        mpi->WriteSingleFile(mpi->result_files[mpi->file_pos[progIO->file_name.max_rhop_vs_scale_file]], progIO->out_content, __all_processors);
    }
}


/*! \fn void BasicAnalysesPostWork()
 *  \brief close default result file  */
void BasicAnalysesPostWork();

/*! \fn template <class T, int D> void MinDistanceBetweenParticles(DataSet<T, D> &ds, int loop_count)
 *  \brief find the minimum distance between particles in data */
template <class T, int D>
void MinDistanceBetweenParticles(DataSet<T, D> &ds, int loop_count)
{
    if (!progIO->flags.basic_analyses_flag) {
        return;
    }
    
    double min_distance = progIO->numerical_parameters.max_half_width, tmp_min_distance = 0;
    uint32_t indices[3];
    for (uint32_t i = 0; i != ds.particle_set.num_particles; i++) {
        ds.tree.KNN_Search(ds.particle_set[i].pos, 2, tmp_min_distance, indices); // the closest one is itself
        min_distance = std::min(min_distance, tmp_min_distance);
    }
    progIO->log_info << "loop_count = " << loop_count << ", min_distance_between_particles in data = " << min_distance << std::endl;
    progIO->Output(std::clog, progIO->log_info, __more_output, __all_processors);
}

#endif /* analyses_hpp */
