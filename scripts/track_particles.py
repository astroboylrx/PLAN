"""
    Run this script to read in a clump ParList output by PLAN and then track all the particles at different snapshots.
    
    usage: track_particles.py [-h] [-c CLUMP] [-i DATA_DIR] [-b BASENAME]
                          [-p POSTNAME] [-f DUMP_ID [DUMP_ID ...]]
                          [-o OUT_DIR] [--Combined]

    optional arguments:
        -h, --help
            show this help message and exit
        -c CLUMP
            path to the clump's particle list file output by PLAN
        -i DATA_DIR
            directory path to the particle data output by Athena, default: ./
        -b BASENAME
            the problem_id used in Athena's output, default: si_trap3d
        -p POSTNAME
            the id name used in Athena's LIS output, default: ds
        -f DUMP_ID [DUMP_ID ...]
            the dump ID of snapshot(s) to track particles, 
            e.g., '-f 270', '-f 270 277', or '-f 250 270 10'
        -o OUT_DIR
            directory path to output new particle lists,
            default: the same directory as the input clump file
        --Combined
            set this flag to deal with combined lis files, default=False

"""

import pyridoxine.utility as rxu
import numpy as np
import os

class ParList:
    """ Read a particle list file """

    def __init__(self, file):
        """  """
        self.id, self.x, self.y, self.z, self.u, self.v, self.w, self.m = rxu.loadtxt(file).T
        self.id = self.id.astype(int)
        self.n = self.id.size

# main function, if you call this as main function
if __name__ == '__main__':

    # dealing with arguments
    import argparse
    parser = argparse.ArgumentParser(description="Run this script to read in a clump ParList output by PLAN and then track all the particles at different snapshots.")
    parser.add_argument('-c', '--clump', dest="clump", help="path to the clump's particle list file output by PLAN", default=None)
    parser.add_argument('-i', '--data_dir', dest="data_dir", help="directory path to the particle data output by Athena, default: ./", default="./")
    parser.add_argument('-b', '--basename', dest="basename", help="the problem_id used in Athena's output, default: si_trap3d", default="si_trap3d")
    parser.add_argument('-p', '--postname', dest="postname", help="the id name used in Athena's LIS output, default: ds", default="ds")
    parser.add_argument('-f', '--dump_id', dest="dump_id", nargs='+', help="the dump ID of snapshot(s) to track particles, e.g., '-f 270', '-f 270 277', or '-f 250 270 10'", type=int, default=[])
    parser.add_argument('-o', '--out_dir', dest="out_dir", help="directory path to output new particle lists, default: the same directory as the input clump file", default=None)
    parser.add_argument("--Combined", dest="combined", help="set this flag to deal with combined lis files, default=False", action='store_const', const=True, default=False)
    args = parser.parse_args()

    if args.clump is None:
        raise ValueError("Please provide a particle list file as the clump of interest to track. Use '-h' to print the usage info.")
    if args.data_dir[-1] is not '/':
        args.data_dir += '/'
    if len(args.dump_id) == 0:
        raise ValueError("Please provide the dump ID of snapshot(s) to track particles. Use '-h' to print the usage info.")
    if len(args.dump_id) == 2:
        if args.dump_id[1] < args.dump_id[0]:
            args.dump_id[0], args.dump_id[1] = args.dump_id[1], args.dump_id[0]
        args.dump_id = np.arange(args.dump_id[0], args.dump_id[1]+1)
    elif len(args.dump_id) == 3:
        if args.dump_id[1] < args.dump_id[0]:
            args.dump_id[0], args.dump_id[1] = args.dump_id[1], args.dump_id[0]
        args.dump_id = np.arange(args.dump_id[0], args.dump_id[1]+1, args.dump_id[2])
    dump_id = ["{:04d}".format(x) for x in args.dump_id]
    if args.out_dir is None:
        args.out_dir = args.clump[:-3] # remove "txt"
    else:
        if args.out_dir[-1] is not '/':
            args.out_dir += '/'
        os.makedirs(args.out_dir, mode=0o755, exist_ok=True)

    print("[info]: reading the particle list file...", end='')
    par_list = ParList(args.clump)
    print("...finished")
    print("[info]: the number of snapshots to track particles: "+str(len(dump_id)))
    
    for i, dump in enumerate(dump_id):
        print("[info]: reading the LIS file with a dump ID: "+dump+"...", end='')
        if args.combined:
            ath_data = rxu.AthenaLIS(args.data_dir + args.basename + '.' + dump + '.' + args.postname + ".lis", sort=True)
        else:
            ath_data = rxu.AthenaMultiLIS(args.data_dir, args.basename, dump+'.'+args.postname+".lis", sort=True)
        print("...finished")
        print("[info]: output selected particles...", end='')
        selected_par = ath_data[ath_data.uni_ids][par_list.id]
        tmp_pos = selected_par['pos']
        tmp_vel = selected_par['vel']
        new_par_list = np.hstack([np.atleast_2d(par_list.id).T, tmp_pos, tmp_vel, np.atleast_2d(par_list.m).T])
        
        np.savetxt(args.out_dir+"at_"+dump+".txt", new_par_list, fmt = "%16d"+"%24.16f"*7)
        print("...finished")
        
        
