//
//  analysis.cpp
//  PLAN: PLantesimal ANalyzer
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
        mpi->OpenFile(mpi->result_file, progIO->file_name.output_file_path);
        progIO->out_content << std::setw(progIO->width) << std::setfill(' ') << "#time" << std::setw(progIO->width) << "max(rho_p)";
        mpi->WriteSingleFile(mpi->result_file, progIO->out_content, __master_only);
    }
}

/*! \fn void BasicAnalysesPostWork()
 *  \brief close default result file  */
void BasicAnalysesPostWork() {
    if (progIO->flags.basic_analyses_flag) {
        mpi->CloseFile(mpi->result_file);
    }
}