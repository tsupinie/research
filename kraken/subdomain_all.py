
import subprocess
import os

from util import editNamelistFile

def generateBatchFile(base_path, job_name, job_n_cpus, commands, wall_time):

    file_text = """#!/bin/csh

#PBS -A TG-MCA95C006
#PBS -N %s 
#PBS -o %s/$PBS_JOBNAME.output
#PBS -e %s/$PBS_JOBNAME.error
#PBS -j oe
#PBS -l walltime=%s:00,size=%d

set base=%s
cd $base

%s

""" % (job_name, base_path, base_path, wall_time, job_n_cpus, base_path, commands)

    return file_text

def main():
    n_ens_members = 40
    t_ens_start = 15000
    t_ens_end = 15000
    dt_ens_step = 300

    t_ens = 13500

    debug = False
    join = True

    base_path = "/lustre/scratch/tsupinie/05June2009/"
    job_name = "1kmf-zs25-offtime-05XP"

    input_path = "%s/input" % base_path
    debug_path = "%s/debug" % base_path
    work_path = "%s/%s" % (base_path, job_name)

    if not join:
        commands = []

    for n_ens in xrange(n_ens_members + 1):
        ena_name = "ena%03d" % (n_ens + 1)
        ena_directory = "EN%03d" % (n_ens + 1)

        if n_ens == n_ens_members:
            ena_name = "enmean"
            ena_directory = "."

        if join:
            commands = []

            start_time = t_ens_start
            for t_ens in xrange(t_ens_start, t_ens_end + dt_ens_step, dt_ens_step):
                if os.path.exists("%s/%s.hdf%06d" % (work_path, ena_name, t_ens)):
                    start_time = t_ens + dt_ens_step

            if start_time > t_ens_end:
                print "Ensemble member %s is done, moving on ..." % ena_name
                continue

        if not join:
           kwargs = { 'hisfile(1)':"%s/%s.hdf%06d" % (work_path, ena_name, t_ens) }

        if join:
            editNamelistFile("%s/hdfsubdomain.input" % base_path, "%s/%s.hdfsubdomain.input" % (input_path, ena_name),
                hdmpinopt=1,
                runname=ena_name,
                dirname="%s" % work_path,
                tbgn_dmpin=start_time,
                tend_dmpin=t_ens_end,
                hdmpfheader="%s/%s/%s" % (work_path, ena_directory, ena_name),
            )
        else:
            editNamelistFile("%s/hdfsubdomain.input" % base_path, "%s/%s.hdfsubdomain.input" % (input_path, ena_name),
                hdmpinopt=2,
                nproc_x=1, nproc_y=1,
                runname=ena_name,
                dirname="%s/%s" % (work_path, ena_directory),
                hdmpfheader="%s/%s/%s" % (work_path, ena_directory, ena_name),
                grdbasfn="%s/%s.hdfgrdbas" % (work_path, ena_name),
                nprocx_in=1, nprocy_in=1,
                nprocx_out=2, nprocy_out=12,
                **kwargs
            )

        if join:
            commands.append("aprun -n 24 $base/hdfsubdomain < %s/%s.hdfsubdomain.input > %s/%s.hdfsubdomain.debug" % (input_path, ena_name, debug_path, ena_name))
        else:
            commands.append("( aprun -n 1 $base/hdfsubdomain < %s/%s.hdfsubdomain.input > %s/%s.hdfsubdomain.debug )&" % (input_path, ena_name, debug_path, ena_name))

        if join:
            text = generateBatchFile(base_path, job_name, 24, "\n".join(commands), "00:10")

            batch_file_name = "/lustre/scratch/tsupinie/%s-hdfsubdomain.csh" % job_name

            if n_ens == 0:
                file = open(batch_file_name, 'w')
                file.write(text)
                file.close()

            if not debug:
                echo = subprocess.Popen([ "echo", text ], stdout=subprocess.PIPE)
                pid = subprocess.check_output([ "qsub" ], stdin=echo.stdout).strip()
                print pid

    if not join:
        commands.append("wait")

        text = generateBatchFile(base_path, "%s" % job_name, 48, "\n".join(commands), "00:30")

        batch_file_name = "/lustre/scratch/tsupinie/%s-hdfsubdomain.csh" % job_name

        file = open(batch_file_name, 'w')
        file.write(text)
        file.close()

        if not debug:
            echo = subprocess.Popen([ "echo", text ], stdout=subprocess.PIPE)
            pid = subprocess.check_output([ "qsub" ], stdin=echo.stdout).strip()
            print pid

    return

if __name__ == "__main__":
    main()
