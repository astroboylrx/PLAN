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
    if (progIO->flags.basic_analyses_flag) {
        mpi->result_files.resize(mpi->result_files.size()+1);
        mpi->files[progIO->file_name.output_file_path] = (--mpi->result_files.end());
        mpi->OpenFile(*mpi->files[progIO->file_name.output_file_path], progIO->file_name.output_file_path);
        progIO->out_content << std::setw(progIO->width) << std::setfill(' ') << "#time" << std::setw(progIO->width) << "max(rho_p)";
        mpi->WriteSingleFile(*mpi->files[progIO->file_name.output_file_path], progIO->out_content, __master_only);
    }
    if (progIO->flags.density_vs_scale_flag) {
        mpi->result_files.resize(mpi->result_files.size()+1);
        mpi->files[progIO->file_name.max_rhop_vs_scale_file] = (--mpi->result_files.end());
        mpi->OpenFile(*mpi->files[progIO->file_name.max_rhop_vs_scale_file], progIO->file_name.max_rhop_vs_scale_file);
    }
}

/*! \fn void BasicAnalysesPostWork()
 *  \brief close default result file  */
void BasicAnalysesPostWork() {
    if (progIO->flags.basic_analyses_flag) {
        mpi->CloseFile(*mpi->files[progIO->file_name.output_file_path]);
    }
    if (progIO->flags.density_vs_scale_flag) {
        mpi->CloseFile(*mpi->files[progIO->file_name.max_rhop_vs_scale_file]);
    }
}