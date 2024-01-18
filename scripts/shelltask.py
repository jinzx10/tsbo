import subprocess

def xabacus(abacus_path, jobdir, nthreads, nprocs, stdout, stderr):
    '''
    Executes ABACUS in a directory.

    Parameters
    ----------
    abacus_path : str
        Path to the ABACUS executable.
    jobdir : str
        Path to the directory where ABACUS is executed.
    nthreads : int
        Number of openmp threads per MPI process.
    nprocs : int
        Number of MPI processes.
    stdout, stderr :
        See the documentation of subprocess.

    '''
    subprocess.run("cd {jobdir}; " \
                   "OMP_NUM_THREADS={nthreads} mpirun -np {nprocs} {abacus_path}" \
                   .format(jobdir=jobdir, nthreads=nthreads, nprocs=nprocs, abacus_path=abacus_path), \
                   shell=True, stdout=stdout, stderr=stderr)


def grep_energy(jobdir, suffix='ABACUS'):
    '''
    Extracts the total energy from the ABACUS output.

    Parameters
    ----------
    jobdir : str
        Path to the directory where ABACUS was executed.
    suffix : str
        Suffix of the ABACUS output directory.

    '''
    result = subprocess.run("grep '!FINAL' {jobdir}/OUT.{suffix}/running_scf.log | awk '{{print $2}}'" \
                            .format(jobdir=jobdir, suffix=suffix),
                            shell=True, capture_output=True, text=True)
    return float(result.stdout)



############################################################
#                       Test
############################################################
import unittest
import shutil
import os
from pathlib import Path
from inputio import write_input
from struio import write_stru

class TestShellTask(unittest.TestCase):

    def test_xabacus_and_grep(self):
        abacus_path = '/home/zuxin/abacus-develop/bin/abacus'
        jobdir = './testfiles/tmp/'
        suffix = 'In2'
        Path(jobdir).mkdir(parents=True, exist_ok=True)

        write_input(jobdir,
                    calculation='scf',
                    suffix=suffix,
                    ecutwfc=50.0,
                    scf_thr=1.0e-8,
                    pseudo_dir='../', 
                    orbital_dir='../',
                    basis_type='lcao',
                    gamma_only=1,
                    )

        species = [
                {'symbol': 'In', 'mass': 1.0, 'pp_file':'In_ONCV_PBE-1.0.upf'},
                ]

        orbitals = [
                'In_gga_10au_100Ry_3s3p3d2f.orb',
                ]

        lattice = {
                'latconst': 20.0,
                'latvec': [
                            [1.0, 0.0, 0.0],
                            [0.0, 1.0, 0.0],
                            [0.0, 0.0, 1.0],
                            ],
                }

        atoms = ['Cartesian_angstrom',
                 {
                    'In': {'mag_each': 0.0,
                            'num': 2,
                            'coord': [
                                [0.0, 0.0, 0.0],
                                [0.0, 0.0, 3.5],
                                ],
                           },
                    }
                 ]

        write_stru(jobdir, species, lattice, atoms, orbitals)

        nthreads = 2
        nprocs = 4
        stdout = subprocess.DEVNULL
        #stdout = None
        stderr = subprocess.DEVNULL

        xabacus(abacus_path, jobdir, nthreads, nprocs, stdout, stderr)
        self.assertTrue(os.path.exists(jobdir + '/OUT.%s/running_scf.log'%(suffix)))

        energy = grep_energy(jobdir, suffix)
        self.assertTrue(abs(energy - (-2912.92846754)) < 1e-7)
        shutil.rmtree(jobdir, ignore_errors=True)


if __name__ == '__main__':
    unittest.main()

