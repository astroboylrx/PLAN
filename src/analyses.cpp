//
//  analysis.cpp
//  PLAN: PLanetesimal ANalyzer
//
//  Created by Rixin Li on 5/4/16.
//  Copyright © 2016 Rixin Li. All rights reserved.
//

/*! \file analyses.cpp
 *  \brief definitions for declarations in analysis.hpp */

#include "analyses.hpp"

/*! \fn void BasicAnalysesPreWork()
 *  \brief open default result file and write file header */
void BasicAnalysesPreWork() {
    // RL: consider creating a map to store flag-file pairs, then we can loop through them to open them
    std::vector<std::string> files_to_open;
    if (progIO->flags.basic_analyses_flag) {
        files_to_open.push_back(progIO->file_name.output_file_path);
        if (dim == 3) {
            files_to_open.push_back(progIO->file_name.mean_sigma_file);
        }
    }
    if (progIO->flags.density_vs_scale_flag) {
        files_to_open.push_back(progIO->file_name.max_rhop_vs_scale_file);
    }
    if (progIO->flags.find_clumps_flag) {
        files_to_open.push_back(progIO->file_name.planetesimals_file);
    }
    for (auto &it : files_to_open) {
        mpi->OpenFile(it);
    }
}

/*! \fn void BasicAnalysesPostWork()
 *  \brief close default result file  */
void BasicAnalysesPostWork() {
    for (auto &it: mpi->result_files) {
        mpi->CloseFile(it);
    }
}
