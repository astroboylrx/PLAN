//
//  main.cpp
//  PLATO: PLAnetesimal locaTOr
//
//  Created by Rixin Li on 3/8/16.
//  Copyright © 2016 Rixin Li. All rights reserved.
//

#include "global.hpp"
#include "mpi.h"

int main(int argc, const char * argv[])
{
    int rank = 0, np = 0;
    MPI::Init(argc, (char **&)argv);
    
    np = MPI::COMM_WORLD.Get_size();
    rank = MPI::COMM_WORLD.Get_rank();
    
    for (int i = 0; i != np; i++) {
        if (i == rank) {
            cout << "Hello World! from cpu " << rank << endl;
        }
    }
    
    MPI::Finalize();
    return 0;
}
