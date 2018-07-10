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
        progIO->numerical_parameters.max_ghost_zone_width = progIO->numerical_parameters.max_half_width;
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
    //if (!progIO->flags.basic_analyses_flag) {
    //    return;
    //}
    
    double min_distance = progIO->numerical_parameters.max_half_width, tmp_min_distance = 0;
    uint32_t indices[3];
    for (uint32_t i = 0; i != ds.particle_set.num_particles; i++) {
        ds.tree.KNN_Search(ds.particle_set[i].pos, 2, tmp_min_distance, indices); // the closest one is itself
        min_distance = std::min(min_distance, tmp_min_distance);
    }
    progIO->log_info << "loop_count = " << loop_count << ", min_distance_between_particles in data = " << min_distance << std::endl;
    progIO->Output(std::clog, progIO->log_info, __more_output, __all_processors);
}

template <class T, int D>
void TempCalculation(DataSet<T, D> &ds, int loop_count)
{
    std::string ori_vtk_file_name, ori_lis_file_name;
    if (!progIO->flags.combined_flag) {
        ori_vtk_file_name = *(progIO->file_name.vtk_data_file_name.begin()+loop_count*progIO->num_cpus);
        ori_lis_file_name = *(progIO->file_name.lis_data_file_name.begin()+loop_count*progIO->num_cpus);
    } else {
        ori_vtk_file_name = *(progIO->file_name.vtk_data_file_name.begin()+loop_count);
        ori_lis_file_name = *(progIO->file_name.lis_data_file_name.begin()+loop_count);
    }

    // RL: debug use, output particles by reading particle id from an external file
    /*
    std::ifstream close_b("./red.txt", std::ifstream::in);
    if (!close_b.is_open()) {
        std::cout << "Failed to open red.txt" << std::endl;
        return;
    }
    std::vector<uint32_t> indices;
    uint32_t tmp_idx;
    double tmp_useless;
    while (!close_b.eof()) {
        close_b >> tmp_idx;
        indices.push_back(tmp_idx);
        //for (auto i = 0; i != 7; i++) {
        //    close_b >> tmp_useless;
        //}
    }
    close_b.close();
    std::cout << "Collected " << indices.size() << " particle indices." << std::endl;

    std::sort(ds.particle_set.particles, ds.particle_set.particles+ds.particle_set.num_particles, [](const Particle<D> &a, const Particle<D> &b) {
        return a.id < b.id;
    });
    ds.planetesimal_list.OutputParticles("red_"+std::to_string(loop_count+progIO->start_num)+".txt", indices, ds.particle_set);
    //*/

    // RL: debug use, output simple POINT3D file for ParaView to visualize
    /*
    std::ofstream file_lis2vtk;
    file_lis2vtk.open(ori_lis_file_name.substr(0, ori_lis_file_name.find("lis")) + "vtk", std::ofstream::out);
    if (!file_lis2vtk.is_open()) {
        std::cout << "Failed to open vtk file" << std::endl;
        return;
    } else {
        std::cout << "Writing to " << ori_lis_file_name.substr(0, ori_lis_file_name.find("lis")) + "vtk" << std::endl;
    }

    file_lis2vtk << "# POINT3D file from t="+std::to_string(ds.particle_set.time) << "\n";
    unsigned long sampling_interval = 1;
    file_lis2vtk << std::scientific;
    for (unsigned long i = 0; i < ds.particle_set.num_particles; i += sampling_interval) {
        file_lis2vtk << std::setw(15) << ds.particle_set[i].pos[0] << std::setw(15) << ds.particle_set[i].pos[1] << std::setw(15) << ds.particle_set[i].pos[2] << std::setw(8) << ds.particle_set[i].id << "\n"; // << " " << ds.particle_set[i].id
    }

    //file_lis2vtk << "POINTS " << indices.size() << " float\n";
    //std::sort(ds.particle_set.particles, ds.particle_set.particles+ds.particle_set.num_particles, [](const Particle<D> &a, const Particle<D> &b) {
    //    return a.id < b.id;
    //});
    //for (auto i : indices) {
    //    file_lis2vtk << ds.particle_set[i].pos[0] << " " << ds.particle_set[i].pos[1] << " " << ds.particle_set[i].pos[2] << "\n";
    //}

    file_lis2vtk.close();
    //*/

    /* RL: visualization use, output Sigma_p with custom resolution
    std::ostringstream time; time.precision(3);
    time << std::setfill('0') << std::setw(6) << std::fixed << ds.particle_set.time;
    std::string Finer_Sigma_p_file_name = ori_lis_file_name.substr(0, ori_lis_file_name.find("lis"))+time.str()+".FSigma_p.txt";
    std::ofstream file_Finer_Sigma_p(Finer_Sigma_p_file_name, std::ofstream::out);
    if (!file_Finer_Sigma_p.is_open()) {
        std::cout << "Failed to open "+Finer_Sigma_p_file_name << std::endl;
    } else {
        std::cout << "Writing to " << Finer_Sigma_p_file_name << std::endl;
    }

    double **Sigma_p = ds.particle_set.MakeFinerSurfaceDensityMap(progIO->numerical_parameters.FineSp_Nx[0], progIO->numerical_parameters.FineSp_Nx[1]);
    file_Finer_Sigma_p << std::scientific;
    for (int iy = 0; iy != progIO->numerical_parameters.FineSp_Nx[1]; iy++) {
        for (int ix = 0; ix != progIO->numerical_parameters.FineSp_Nx[0]; ix++) {
            file_Finer_Sigma_p << std::setw(16) << std::setprecision(8) << Sigma_p[iy][ix];
        }
        file_Finer_Sigma_p << std::endl;
    }
    file_Finer_Sigma_p.close();
    if (Sigma_p != nullptr) {
        delete [] Sigma_p[0];
        Sigma_p[0] = nullptr;
    }
    delete [] Sigma_p;
    Sigma_p = nullptr;

    // RL: debug use, output Sigma_p
    /*
    std::string Sigma_p_file_name = ori_lis_file_name.substr(0, ori_lis_file_name.find("lis")) + "Sigma_p.txt";
    std::ofstream file_Sigma_p(Sigma_p_file_name, std::ofstream::out);
    if (!file_Sigma_p.is_open()) {
        std::cout << "Failed to open "+Sigma_p_file_name << std::endl;
        return;
    } else {
        std::cout << "Writing to " << Sigma_p_file_name << std::endl;
    }

    VtkDataScalar<T, 3> tmp_rhop = VtkDataScalar<T, 3>();
    for (int i = 2; i >= 0; i--) {
        tmp_rhop.num_cells[i] = progIO->numerical_parameters.box_resolution[i];
        tmp_rhop.shape[i] = progIO->numerical_parameters.box_resolution[dim-1-i];
    }
    tmp_rhop.data.resize(tmp_rhop.shape);

    Particle<D> *p;
    SmallVec<int, 3> tmp_index;
    for (uint32_t i = 0; i < ds.particle_set.num_particles; i++) {
        p = &ds.particle_set.particles[i];
        for (int j = 0; j != D; j++) {
            tmp_index[j] = static_cast<int>(std::floor((p->pos[j] - progIO->numerical_parameters.box_min[j]) / progIO->numerical_parameters.cell_length[j]));
        }
        if (tmp_index[2] >= progIO->numerical_parameters.box_resolution[2]) {
            continue;
        }
        tmp_rhop.data[tmp_index[2]][tmp_index[1]][tmp_index[0]] = p->density;
    }

    file_Sigma_p << std::scientific;
    for (auto iy = 0; iy != tmp_rhop.num_cells[1]; iy++) {
        for (auto ix = 0; ix != tmp_rhop.num_cells[0]; ix++) {
            double Sigma_p = 0;
            typename VtkDataScalar<T, D>::template view_r2d_type<> column_density(tmp_rhop.data[boost::indices[sn::b_range()][iy][ix]]);
            ds.vtk_data.IterateBoostMultiArrayConcept(column_density, [&Sigma_p](float &element)->void {
                Sigma_p += element;
            });
            file_Sigma_p << std::setw(16) << std::setprecision(8) << Sigma_p;
        }
        file_Sigma_p << std::endl;
    }
    file_Sigma_p.close();

    //*/

    // RL: debug use, output sub-sampled lis file for fast analyses
    /*
    std::ofstream file_subsample;
    unsigned long sampling_in_each_cpu = 4096;
    std::string new_lis_file_name = ori_lis_file_name.substr(0, ori_lis_file_name.find("all")) + "sub.lis";
    std::cout << "Sub-sampling to " << new_lis_file_name << std::endl;
    file_subsample.open(new_lis_file_name, std::ios::binary);
    if (file_subsample.is_open()) {
        int tmp_int;
        long tmp_num_particles, tmp_long;
        float tmp_float_value, tmp_float_vector[D];
        for (int i = 0; i != 12; i++) {
            tmp_float_value = static_cast<float>(ds.particle_set.coor_lim[i]);
            file_subsample.write(reinterpret_cast<char*>(&tmp_float_value), sizeof(float));
        }
        file_subsample.write(reinterpret_cast<char*>(&ds.particle_set.num_types), sizeof(int));
        for (unsigned int i = 0; i != ds.particle_set.num_types; i++) {
            tmp_float_value = static_cast<float>(ds.particle_set.type_info[i]);
            file_subsample.write(reinterpret_cast<char*>(&tmp_float_value), sizeof(float));
        }
        tmp_float_value = static_cast<float>(ds.particle_set.time);
        file_subsample.write(reinterpret_cast<char*>(&tmp_float_value), sizeof(float));
        tmp_float_value = static_cast<float>(ds.particle_set.dt);
        file_subsample.write(reinterpret_cast<char*>(&tmp_float_value), sizeof(float));

        Particle<D> *p;
        uint32_t num_particles_to_output = 0;
        for (uint32_t i = 0; i < ds.particle_set.num_particles; i += 1) {
            p = &ds.particle_set.particles[i];
            if (p->id_in_run < sampling_in_each_cpu) {
                num_particles_to_output++;
            }
        }
        std::cout << "We will subsample " << num_particles_to_output << " particles. " << std::endl;

        tmp_num_particles = static_cast<long>(num_particles_to_output);
        file_subsample.write(reinterpret_cast<char*>(&tmp_num_particles), sizeof(long));

        size_t D_float = D * sizeof(float);
        size_t one_float = sizeof(float);
        size_t one_int = sizeof(int);
        size_t one_long = sizeof(long);

        for (uint32_t i = 0; i < ds.particle_set.num_particles; i += 1) {
            p = &ds.particle_set.particles[i];
            if (p->id_in_run < sampling_in_each_cpu) {
                for (int i = 0; i != D; i++) {
                    tmp_float_vector[i] = static_cast<float>(p->pos[i]);
                }
                file_subsample.write(reinterpret_cast<char *>(&tmp_float_vector), D_float);
                for (int i = 0; i != D; i++) {
                    tmp_float_vector[i] = static_cast<float>(p->vel[i]);
                }
                file_subsample.write(reinterpret_cast<char *>(&tmp_float_vector), D_float);
                tmp_float_value = static_cast<float>(p->density);
                file_subsample.write(reinterpret_cast<char *>(&tmp_float_value), one_float);

                file_subsample.write(reinterpret_cast<char *>(&p->property_index), one_int);
                tmp_long = static_cast<long>(p->id_in_run);
                file_subsample.write(reinterpret_cast<char *>(&tmp_long), one_long);
                tmp_int = static_cast<int>(p->cpu_id);
                file_subsample.write(reinterpret_cast<char *>(&tmp_int), one_int);
            }
        }
    }
    file_subsample.close();
    //*/

    // RL: test use, output particle densities from vtk
    /*
    std::ofstream file_dpar;
    file_dpar.open(progIO->file_name.output_file_path.substr(0, progIO->file_name.output_file_path.find_last_of('/'))+"/dpar.dat", std::ios::binary|std::ios::app);
    // app (append) Set the stream's position indicator to the end of the stream before each output operation
    // so cannot go back to a position in ofstream to change content

    size_t one_float = sizeof(float);
    uint32_t npar_count = 0;

    if (file_dpar.is_open()) {
        file_dpar.write(reinterpret_cast<char*>(&ds.vtk_data.time), sizeof(double));

        ds.vtk_data.IterateBoostMultiArrayConcept(ds.vtk_data.scalar_data["particle_density"].data, [&npar_count](float &element) ->void {
            if (element > 0) {
                npar_count++;
            }
        });
        std::cout << "npar_count = " << npar_count << std::endl;
        file_dpar.write(reinterpret_cast<char*>(&npar_count), sizeof(uint32_t));

        ds.vtk_data.IterateBoostMultiArrayConcept(ds.vtk_data.scalar_data["particle_density"].data, [&file_dpar, &one_float](float &element) ->void {
            if (element > 0) {
                file_dpar.write(reinterpret_cast<char*>(&element), one_float);
            }
        });
        file_dpar << std::endl;
        file_dpar.close();
    } else {
        std::cout << "Cannot open dpar.dat" << std::endl;
    }
    //*/


};

#endif /* analyses_hpp */
