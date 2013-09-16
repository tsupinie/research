
import os
import subprocess 
import time
import argparse
from math import ceil

from util import editNamelistFile

def parseQStat(text):
    field_widths = []
    qstat = []
    for line in text:
        if line[:3] == '---':
            field_widths = [ len(f) for f in line.split() ]
            continue

        if field_widths != []:
            fields = []
            offset = 0
            for width in field_widths:
                fields.append(line[(offset):(offset + width)])
                offset += width + 1

            print fields
            qstat.append(fields)
    return qstat

def isDivisible(dividend, divisor):
    return float(dividend) / int(divisor) == int(dividend) / int(divisor)

def doForEnsembleDict(commands, member_list, current_commands=None):

    if type(commands) not in [ list, tuple ]:
        commands = [ commands ]

    if current_commands is None:
        current_commands = {}

    for n_ens in member_list:
        key = "ena%03d" % (n_ens + 1)
        for cmd in commands:
            try:
                current_commands[key].append(cmd % {'ens':(n_ens + 1)})
            except KeyError:
                current_commands[key] = [ cmd % {'ens':(n_ens + 1)} ]

    return current_commands

def doForEnsemble(commands, member_list, footprint, n_cores_available):
    command_lines = []

    if type(commands) not in [ list, tuple ]:
        commands = [ commands ]

    even_division = False

    for mem_number, n_ens in enumerate(member_list):
        even_division = False

        for cmd in commands:
            full_command = cmd % {'ens':(n_ens + 1)}
            command_lines.append("( %s )&" % full_command)

        if not ((mem_number + 1) % (n_cores_available / (footprint * len(commands)))):
            even_division = True

            command_lines.append("")
            command_lines.append("wait")
            command_lines.append("")

    if not even_division:
        command_lines.append("")
        command_lines.append("wait")
        command_lines.append("")

    return command_lines

def doForMPIConfig(command, mpi_config):
    nproc_x, nproc_y = mpi_config

    commands = []

    for idx in range(nproc_x):
        for jdy in range(nproc_y):
            commands.append(command % {'x_proc':idx + 1, 'y_proc':jdy + 1, 'ens':'%(ens)03d' })

    return commands

def generateBatchFile(base_path, job_name, job_n_cpus, commands, wall_time):

    file_text = """#!/bin/csh

#PBS -A TG-MCA95C006
#PBS -N %s 
#PBS -o %s/debug/$PBS_JOBNAME.output
#PBS -e %s/debug/$PBS_JOBNAME.error
#PBS -j oe
#PBS -l walltime=%s:00,size=%d

set base=%s
cd $base

%s

""" % (job_name, base_path, base_path, wall_time, job_n_cpus, base_path, commands)

    return file_text

def generateEnsembleIntegration(base_path, job_name, member_list, mpi_config, n_cores, start_time, end_time, dump_time, split_files='neither', move_for_assim=True, current_commands=None, input_file_name="arps.input"):
    work_path = "%s/%s/" % (base_path, job_name)
    input_path = "%s/input/" % base_path
    debug_path = "%s/debug/" % base_path
    batch_path = "%s/batch/" % base_path
    bc_path = "%s/boundary/" % base_path

    nproc_x, nproc_y = mpi_config

    extraneous_files = [
        "%s/%s.hdf%06d.01" % (work_path, 'ena%(ens)03d', start_time),
        "%s/%s.hdfgrdbas.01" % (work_path, 'ena%(ens)03d'),
        "%s/%s.log" % (work_path, 'ena%(ens)03d'),
        "%s/%s.maxmin" % (work_path, 'ena%(ens)03d'),
    ]

    read_split = split_files in ['read', 'both']
    dump_split = split_files in ['dump', 'both']

    if move_for_assim:
        if dump_split:
            epilogue = doForMPIConfig("mv %s/%s/%s.hdf%06d_%s %s/%s/%s.hdf%06d_%s" % (work_path, 'EN%(ens)s', 'ena%(ens)s', end_time, '%(x_proc)03d%(y_proc)03d', work_path, 'ENF%(ens)s', 'enf%(ens)s', end_time, '%(x_proc)03d%(y_proc)03d'), mpi_config)
        else:
            epilogue = [
                "mv %s/%s.hdf%06d %s/%s.hdf%06d" % (work_path, 'ena%(ens)03d', end_time, work_path, 'enf%(ens)03d', end_time)
            ]
    else:
        epilogue = []

    for n_ens in member_list:
        ens_member_name = "ena%03d" % (n_ens + 1)
        ens_member_directory = "EN%03d" % (n_ens + 1)

        arps_input_file_name = "%s/%s.%d-%d.arps.input" % (input_path, ens_member_name, start_time, end_time)

        if read_split:
            init_file = "%s/%s/%s.hdf%06d" % (work_path, ens_member_directory, ens_member_name, start_time),
            init_grdbas = "%s/%s/%s.hdfgrdbas" % (work_path, ens_member_directory, ens_member_name),
        else:
            init_file = "%s/%s.hdf%06d" % (work_path, ens_member_name, start_time),
            init_grdbas = "%s/%s.hdfgrdbas" % (work_path, ens_member_name),

        if dump_split:
            out_dir="%s/%s/" % (work_path, ens_member_directory)
        else:
            out_dir="%s/" % work_path

        editNamelistFile("%s/%s" % (base_path, input_file_name), arps_input_file_name,
            nproc_x=nproc_x, nproc_y=nproc_y,
            runname=ens_member_name,
            initopt=3,
            inifile=init_file,
            inigbf=init_grdbas,
            exbcname="%s/%s" % (bc_path, ens_member_name),
            tstart=start_time,
            tstop=end_time,
            tstrtdmp=start_time + dump_time,
            thisdmp=dump_time,
            dmp_out_joined=int(not dump_split),
            inisplited=int(read_split),
            dirname=out_dir)

    arps_input_file_name = "%s/%s.%d-%d.arps.input" % (input_path, 'ena%(ens)03d', start_time, end_time)
    arps_debug_file_name = "%s/%s.%d-%d.arps.debug" % (debug_path, 'ena%(ens)03d', start_time, end_time)

    command = "aprun -n %d $base/arps %s > %s ; rm %s ; %s" % \
        (nproc_x * nproc_y, arps_input_file_name, arps_debug_file_name, " ".join(extraneous_files), " ; ".join(epilogue))

    command_lines = doForEnsembleDict(command, member_list, current_commands)

    return command_lines

def generateEnsemblePerturbations(base_path, job_name, member_list, n_cores, start_time, current_commands=None, input_file_name="arpsenkfic.input"):
    work_path = "%s/%s/" % (base_path, job_name)
    input_path = "%s/input/" % base_path
    debug_path = "%s/debug/" % base_path

    for n_ens in member_list:
        ens_member_name = "ena%03d" % (n_ens + 1)

        arpsenkfic_input_file_name = "%s/%s.arpsenkfic.input" % (input_path, ens_member_name)

        editNamelistFile("%s/%s" % (base_path, input_file_name), arpsenkfic_input_file_name,
            seeds=-n_ens,
            dirnamp="%s/" % work_path,
            outdumpdir="%s/" % work_path,
            outname=ens_member_name,
            tfgs=start_time)

    arps_input_file_name = "%s/arps.input" % input_path
    arpsenkfic_input_file_name = "%s/%s.arpsenkfic.input" % (input_path, 'ena%(ens)03d')
    arpsenkfic_debug_file_name = "%s/%s.arpsenkfic.debug" % (debug_path, 'ena%(ens)03d')

    command = "aprun -n 1 $base/arpsenkfic %s < %s > %s" % \
        (arps_input_file_name, arpsenkfic_input_file_name, arpsenkfic_debug_file_name)

    command_lines = doForEnsembleDict(command, member_list, current_commands)

    return command_lines

def generateEnKFAssimilation(base_path, job_name, n_ens_members, mpi_config, n_cores, assim_time, radar_data_flag, split_files, covariance_inflation, assim_all=False, input_file_name="arpsenkf.input"):
    work_path = "%s/%s/" % (base_path, job_name)
    input_path = "%s/input/" % base_path
    debug_path = "%s/debug/" % base_path
    batch_path = "%s/batch/" % base_path

    nproc_x, nproc_y = mpi_config

    arps_input_file_name = "%s/arps.input" % input_path
    enkf_input_file_name = "%s/%d.arpsenkf.input" % (input_path, assim_time)
    enkf_debug_file_name = "%s/%d.arpsenkf.debug" % (debug_path, assim_time)
    batch_file_name = "%s/%d.csh" % (batch_path, assim_time)
    submit_job_name = "%s.%d" % (job_name, assim_time)

    if assim_all:
        kwargs = {'sndassim':1, 'proassim':1}
    else:
        kwargs = {}

    if split_files:
        kwargs['nproc_x_in'] = nproc_x
        kwargs['nproc_y_in'] = nproc_y
        kwargs['inidirname'] = "%s/ENF%s" % (work_path, "%3N")
    else:
        kwargs['inidirname'] = work_path

    if assim_time == 5400:
        kwargs['vcpmode(2)'] = 32
        kwargs['vcpmode(3)'] = 32

    kwargs['mult_inflat'] = 0
    kwargs['adapt_inflat'] = 0
    kwargs['relax_inflat'] = 0

    try:
        covariance_inflation = covariance_inflation.split(',')
    except ValueError:
        covariance_inflation = [ covariance_inflation ]

    for cov_infl in covariance_inflation:
        inflation_method, inflation_factor = cov_infl.split('=')
        inflation_factor = float(inflation_factor)

        if inflation_method == "mult":
#           if kwargs['adapt_inflat'] > 0:
#               kwargs['mult_inflat'] = 1
#           else:
#               kwargs['mult_inflat'] = 2

            if inflation_factor > 1.05:
                kwargs['mult_inflat'] = 1
            else:
                kwargs['mult_inflat'] = 2

            kwargs['cinf'] = inflation_factor
        elif inflation_method == "adapt":
            kwargs['adapt_inflat'] = 1
            kwargs['rlxf'] = inflation_factor

#           if kwargs['mult_inflat'] == 2:
#               kwargs['mult_inflat'] = 1

        elif inflation_method == "relax":
            kwargs['relax_inflat'] = 1
            kwargs['rlxf'] = inflation_factor

    n_radars = len([ f for f in radar_data_flag.values() if f ])
    radardaopt = 0
    if n_radars > 0:
        radardaopt = 1

    editNamelistFile("%s/%s" % (base_path, input_file_name), "%s/%d.arpsenkf.input" % (input_path, assim_time), 
        nen=n_ens_members,
        enkfdtadir="%s/" % work_path,
        assim_time=assim_time,
        radardaopt=radardaopt,
        nrdrused=n_radars,
        rmsfcst=2,
        **kwargs
    )

    command_lines = []

    command_lines.append("cd %s" % work_path)
    command_lines.append("aprun -n %d -d 12 $base/arpsenkf %s < %s > %s" % (nproc_x * nproc_y, arps_input_file_name, enkf_input_file_name, enkf_debug_file_name))
    command_lines.append("cd -")
    command_lines.append("")

    return command_lines

def generateDomainSubset(base_path, job_name, ic_path, member_list, n_cores, start_time, end_time, step_time, current_commands=None, perturb_ic=True):
    work_path = "%s/%s/" % (base_path, job_name)
    input_path = "%s/input/" % base_path
    debug_path = "%s/debug/" % base_path
    bc_path = "%s/boundary/" % base_path

    for n_ens in member_list:
        ens_member_name = "ena%03d" % (n_ens + 1)

        interp_input_file_name = "%s/%s.arpsintrp.input" % (input_path, ens_member_name)
        arpsenkfic_input_file_name = "%s/%s.arpsenkfic.input" % (input_path, ens_member_name)
        arps_input_file_name = "%s/%s.arps.input" % (input_path, ens_member_name)

        editNamelistFile("%s/arpsintrp.input" % base_path, interp_input_file_name,
            runname="%s" % ens_member_name,
            hdmpfheader="%s/%s" % (ic_path, ens_member_name),
            dirname=bc_path,
            tbgn_dmpin=start_time,
            tend_dmpin=end_time,
            tintv_dmpin=step_time,
        )

        editNamelistFile("%s/arps.input" % base_path, arps_input_file_name,
            runname="%s" % ens_member_name,
            initopt=3,
            inifile="%s/%s.hdf%06d" % (bc_path, ens_member_name, start_time),
            inigbf="%s/%s.hdfgrdbas" % (bc_path, ens_member_name),
            tstart=start_time,
            tstop=end_time,
            dirname="%s/" % work_path
        )

        if perturb_ic:
            editNamelistFile("%s/arpsenkfic.input" % base_path, arpsenkfic_input_file_name,
                seeds=-n_ens,
                dirnamp="%s/" % work_path,
                outdumpdir="%s/" % work_path,
                outname=ens_member_name,
                tfgs=start_time
            )

    interp_input_file_name = "%s/ena%s.arpsintrp.input" % (input_path, '%(ens)03d')
    interp_debug_file_name = "%s/ena%s.arpsintrp.debug" % (debug_path, '%(ens)03d')
    arps_input_file_name = "%s/ena%s.arps.input" % (input_path, '%(ens)03d')
    arpsenkfic_input_file_name = "%s/ena%s.arpsenkfic.input" % (input_path, '%(ens)03d')
    arpsenkfic_debug_file_name = "%s/ena%s.arpsenkfic.debug" % (debug_path, '%(ens)03d')

    commands = [ "rm %s/%sicbc.*" % (bc_path, "ena%(ens)03d"), "" ]
    if perturb_ic:
        commands.append("aprun -n 1 $base/arpsintrp %s/input/arps.input < %s > %s ; aprun -n 1 -d 12 $base/arpsenkfic %s < %s > %s" % (base_path, interp_input_file_name, 
            interp_debug_file_name, arps_input_file_name, arpsenkfic_input_file_name, arpsenkfic_debug_file_name))
    else:
        commands.append("aprun -n 1 $base/arpsintrp %s/input/arps.input < %s > %s" % (base_path, interp_input_file_name, interp_debug_file_name))
        commands.extend([
           "cp %s/%s.hdf%06d %s" % (bc_path, 'ena%(ens)03d', start_time, work_path),
           "cp %s/%s.hdfgrdbas %s" % (bc_path, 'ena%(ens)03d', work_path),
           "cp %s/%s.hdfgrdbas %s/%s.hdfgrdbas" % (bc_path, 'ena%(ens)03d', work_path, 'enf%(ens)03d')
        ])

    command_lines = doForEnsembleDict(commands, member_list, current_commands)

    return command_lines

def submit(base_path, job_name, n_cores_available, command_lines, wall_time, submit, sleep):
    job_suffixes = []
    job_batch = []
    if type(command_lines) in [ dict ]:
        written_example = False
        for key, commands in command_lines.iteritems():
            file_text = generateBatchFile(base_path, "%s-%s" % (key, job_name), n_cores_available / len(command_lines.keys()), "\n".join(commands), wall_time)
            job_batch.append(file_text)
            job_suffixes.append(key)

            if not written_example:
                file = open("%s.csh" % job_name, 'w')
                file.write(file_text)
                file.close()
                written_example = True

    elif type(command_lines) in [ list, tuple ]:
        file_text = generateBatchFile(base_path, "enkf-%s" % job_name, n_cores_available, "\n".join(command_lines), wall_time)
        job_batch.append(file_text)
        job_suffixes.append("enkf")

        file = open("%s.csh" % job_name, 'w')
        file.write(file_text)
        file.close()

    if sleep:
        job_completed = []
        if submit:
            for batch in job_batch:
                echo = subprocess.Popen([ "echo", batch ], stdout=subprocess.PIPE)
                pid = subprocess.check_output([ "qsub" ], stdin=echo.stdout).strip()
                print pid
                job_completed.append(False)
        else:
            for batch in job_batch:
                job_completed.append(False)

            print "Done submitting ..."

        while True:
            stat = subprocess.check_output(['qstat', '-u', 'tsupinie']).split("\n")
            stat = parseQStat(stat)
            if len(stat) > 0:
                job_name_len = len(stat[0][3])
                jobs_queued = [ r[3].strip() for r in stat ]

                for idx, suffix in enumerate(job_suffixes):
                    key = "%s-%s" % (suffix, job_name)
                    if key[:job_name_len] in jobs_queued:
                        jdy = jobs_queued.index(key[:job_name_len])
                        if stat[jdy][9] == 'C':
                            job_completed[idx] = True
            else:
                for idx in range(len(job_completed)):
                    job_completed[idx] = True

            print job_completed
            if all(job_completed):
                print "Jobs are completed, returning for the next cycle ..."
#               time.sleep(6 * 60)
                return        

            time.sleep(4 * 60)
    elif submit:   
        job_id = subprocess.check_output([ 'qsub', "%s.csh" % job_name ])
        print job_id

    return

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n-ens', dest='n_ens_members', default=4, type=int)
    ap.add_argument('--members', dest='members', nargs='+', default=[], type=int)
    ap.add_argument('--base-path', dest='base_path', default=os.getcwd())
    ap.add_argument('--job-name', dest='job_name', default="run_osse_test")
    ap.add_argument('--n-cores', dest='n_cores_available', default=1, type=int)
    ap.add_argument('--mpi-config', dest='mpi_config', nargs=2, default=(3, 4), type=int)

    ap.add_argument('--ens-start', dest='t_ens_start', default=1200, type=int)
    ap.add_argument('--ens-end', dest='t_ens_end', default=1500, type=int)
    ap.add_argument('--assim-step', dest='dt_assim_step', default=300, type=int)
    ap.add_argument('--ens-step', dest='dt_ens_step', default=300, type=int)

    ap.add_argument('--assimilate', dest='osse_assimilate', action='store_true')
    ap.add_argument('--initial-conditions', dest='init_cond', default='')
    ap.add_argument('--subset-ic', dest='subset', action='store_true')
    ap.add_argument('--free-forecast', dest='free_forecast', action='store_true')
    ap.add_argument('--covariance-inflation', dest='cov_infl', nargs='+', default=["mult=1.1"])
    ap.add_argument('--boundary-conditions', dest='bound_cond', default='')

    ap.add_argument('--init-fcst-req', dest='init_fcst_req', default='0:40')
    ap.add_argument('--fcst-req', dest='fcst_req', default='0:20')
    ap.add_argument('--init-free-fcst-req', dest='init_free_fcst_req', default='1:45')
    ap.add_argument('--free-fcst-req', dest='free_fcst_req', default='1:30')
    ap.add_argument('--assim-on-req', dest='assim_on_req', default='1:00')
    ap.add_argument('--assim-off-req', dest='assim_off_req', default='0:45')

    ap.add_argument('--piecewise', dest='piecewise', action='store_true')
    ap.add_argument('--chunk-size', dest='chunk_size', default=1200, type=int)
    ap.add_argument('--split-files', dest='split_files', action='store_true')
    ap.add_argument('--split-init', dest='split_init', choices=['auto', 'yes', 'no'], default='auto')

    ap.add_argument('--no-submit', dest='no_submit', action='store_true')
    ap.add_argument('--restart', dest='restart', action='store_true')
    ap.add_argument('--debug', dest='debug', action='store_true')

    args = ap.parse_args()

    work_path = "%s/%s" % (args.base_path, args.job_name)
    input_path = "%s/input" % args.base_path

    member_list = [ m - 1 for m in args.members ]
    if member_list == []: member_list = range(args.n_ens_members)

    boundary_path = args.bound_cond
    if args.bound_cond == "":
        boundary_path = "%s/boundary/" % args.base_path

    nproc_x, nproc_y = args.mpi_config

    command_lines = {}

    prologue = []

    prologue.append("setenv OMP_NUM_THREADS 12") # % (args.n_cores_available / (nproc_x * nproc_y)))
    prologue.append("")

    prologue.append("echo \"\" > %s/debug/ena%s-%s.output" % (args.base_path, '%(ens)03d', args.job_name))
    prologue.append("")

    command_lines = doForEnsembleDict(prologue, member_list, current_commands=command_lines)

    joined = 1
    if args.split_files: joined = 0

    editNamelistFile("%s/arps.input" % args.base_path, "%s/arps.input" % input_path,
        nproc_x=nproc_x, nproc_y=nproc_y,
        dmp_out_joined=joined,
        inisplited=3 * (1 - joined),
        sfcdat=3,
        dirname=work_path)

    do_first_integration = True
    ens_chunk_start = args.t_ens_start

    if not args.free_forecast and args.chunk_size > args.dt_assim_step:
        args.chunk_size = args.dt_assim_step

    radar_data_flag_3km = {
        1800: {'KCYS':False, 'KFTG':False, 'KRIW':True, '05XP':False }, #1830
        3600: {'KCYS':False, 'KFTG':False, 'KRIW':True, '05XP':False }, #1900
        5400: {'KCYS':True,  'KFTG':True,  'KRIW':True, '05XP':False }, #1930
        7200: {'KCYS':True,  'KFTG':True,  'KRIW':True, '05XP':False }, #2000
        9000: {'KCYS':True,  'KFTG':True,  'KRIW':True, '05XP':False }, #2030
        9600: {'KCYS':True,  'KFTG':True,  'KRIW':True, '05XP':False }, #2040
        10200:{'KCYS':True,  'KFTG':True,  'KRIW':True, '05XP':False }, #2050
        10800:{'KCYS':True,  'KFTG':True,  'KRIW':True, '05XP':False }, #2100

        11100:{'KCYS':True,  'KFTG':True,  'KRIW':True, '05XP':False }, #2105
        11400:{'KCYS':True,  'KFTG':True,  'KRIW':True, '05XP':False }, #2110
        11700:{'KCYS':True,  'KFTG':True,  'KRIW':True, '05XP':False }, #2115
        12000:{'KCYS':True,  'KFTG':True,  'KRIW':True, '05XP':False }, #2120
        12300:{'KCYS':True,  'KFTG':True,  'KRIW':True, '05XP':False }, #2125
        12600:{'KCYS':True,  'KFTG':True,  'KRIW':True, '05XP':False }, #2130
        12900:{'KCYS':True,  'KFTG':True,  'KRIW':True, '05XP':False }, #2135
        13200:{'KCYS':True,  'KFTG':True,  'KRIW':True, '05XP':False }, #2140
        13500:{'KCYS':True,  'KFTG':True,  'KRIW':True, '05XP':False }, #2145
        13800:{'KCYS':True,  'KFTG':True,  'KRIW':True, '05XP':False }, #2150
        14100:{'KCYS':True,  'KFTG':True,  'KRIW':True, '05XP':False }, #2155
        14400:{'KCYS':True,  'KFTG':True,  'KRIW':True, '05XP':False }, #2200
    }

    radar_data_flag_1km = {
        9600: {'KCYS':True,  'KFTG':True,  'KRIW':False, '05XP':False}, #2040
        10200:{'KCYS':True,  'KFTG':True,  'KRIW':False, '05XP':False}, #2050
        10800:{'KCYS':True,  'KFTG':True,  'KRIW':False, '05XP':False}, #2100

        11100:{'KCYS':True,  'KFTG':False, 'KRIW':False, '05XP':False}, #2105
        11400:{'KCYS':True,  'KFTG':False, 'KRIW':False, '05XP':False}, #2110
        11700:{'KCYS':True,  'KFTG':True,  'KRIW':False, '05XP':False}, #2115
        12000:{'KCYS':True,  'KFTG':True,  'KRIW':False, '05XP':False}, #2120
        12300:{'KCYS':True,  'KFTG':True,  'KRIW':False, '05XP':False}, #2125
        12600:{'KCYS':True,  'KFTG':True,  'KRIW':False, '05XP':False}, #2130
        12900:{'KCYS':True,  'KFTG':True,  'KRIW':False, '05XP':False}, #2135
        13200:{'KCYS':True,  'KFTG':True,  'KRIW':False, '05XP':False}, #2140
        13500:{'KCYS':True,  'KFTG':True,  'KRIW':False, '05XP':True }, #2145
        13800:{'KCYS':True,  'KFTG':True,  'KRIW':False, '05XP':True }, #2150
        14100:{'KCYS':True,  'KFTG':True,  'KRIW':False, '05XP':True }, #2155
        14400:{'KCYS':True,  'KFTG':True,  'KRIW':False, '05XP':True }, #2200
    }

    radar_data_flag = radar_data_flag_1km

    exp_start = args.t_ens_start

    # Copy the configuration information to the working directory, so we'll **ALWAYS HAVE IT IF WE NEED TO GO BACK AND LOOK AT IT**
    if radar_data_flag == radar_data_flag_1km:
        input_files = ['arps.1km.input', 'arpsenkf.1km.input', 'arpsenkf.1km.05XP.input', 'arpsenkfic.1km.input']
    elif radar_data_flag == radar_data_flag_3km:
        input_files = ['arps.3km.input', 'arpsenkf.3km.input', 'arpsenkfic.3km.input']

    config_files = [ "%s/%s" % (args.base_path, f) for f in input_files ]
    config_files.extend(['run_real_data_case.py', 'run_real_data_case.csh'])

    for file in config_files:
        subprocess.Popen(['cp', file, work_path])

    if args.restart:
        for t_ens in xrange(args.t_ens_start, args.t_ens_end + args.dt_assim_step, args.dt_assim_step):

            if args.split_files:
                ena_exist = [ os.path.exists("%s/EN%03d/ena%03d.hdf%06d_001001" % (work_path, n_ens + 1, n_ens + 1, t_ens)) for n_ens in member_list ]
                enf_exist = [ os.path.exists("%s/ENF%03d/enf%03d.hdf%06d_001001" % (work_path, n_ens + 1, n_ens + 1, t_ens)) for n_ens in member_list ]
            else:
                ena_exist = [ os.path.exists("%s/ena%03d.hdf%06d" % (work_path, n_ens + 1, t_ens)) for n_ens in member_list ]
                enf_exist = [ os.path.exists("%s/enf%03d.hdf%06d" % (work_path, n_ens + 1, t_ens)) for n_ens in member_list ]

            all_ena_exist = all(ena_exist)
            all_enf_exist = all(enf_exist)

            if all_ena_exist:
                args.t_ens_start = t_ens

                if args.piecewise:
                    ens_chunk_start = t_ens
            elif all_enf_exist and not all_ena_exist:
                args.t_ens_start = t_ens - args.dt_assim_step
                do_first_integration = False

            if args.piecewise:
                for t_chunk in xrange(t_ens, t_ens + args.dt_assim_step, args.chunk_size):
                    if args.split_files:
                        ena_exist = [ os.path.exists("%s/EN%03d/ena%03d.hdf%06d_001001" % (work_path, n_ens + 1, n_ens + 1, t_chunk)) for n_ens in member_list ]
                    else:
                        ena_exist = [ os.path.exists("%s/ena%03d.hdf%06d" % (work_path, n_ens + 1, t_chunk)) for n_ens in member_list ]

                    if all(ena_exist):
                        ens_chunk_start = t_chunk

        if do_first_integration:
            print "Restarting from time %d (with integration) ..." % (args.t_ens_start)
        else:
            print "Restarting from time %d (no integration) ..." % (args.t_ens_start + args.dt_assim_step)

    else:
        print "No restart ..."
#       command_lines.append("rm %s/*" % work_path)
#       command_lines.append("")

        if args.split_files:
            command = "mkdir %s/%s ; mkdir %s/%s" % (work_path, 'EN%(ens)03d', work_path, 'ENF%(ens)03d')
            command_lines = doForEnsembleDict(command, member_list, current_commands=command_lines)

        if args.init_cond == "":
            print "Generate random initial conditions ..."
            command_lines = generateEnsemblePerturbations(args.base_path, args.job_name, member_list, args.n_cores_available, args.t_ens_start, current_commands=command_lines)

            command = "cp %s/%s.hdfgrdbas %s/%s.hdfgrdbas" % (work_path, 'ena%(ens)03d', work_path, 'enf%(ens)03d')
            command_lines = doForEnsembleDict(command, member_list, current_commands=command_lines)
        else:
            print "Use supplied initial conditions ..."
            if args.subset:
                print "Subset and perturb the domain ..."
                command_lines = generateDomainSubset(args.base_path, args.job_name, args.init_cond, member_list, args.n_cores_available, args.t_ens_start, args.t_ens_end, args.dt_ens_step, current_commands=command_lines, perturb_ic=True)

                command = "cp %s/%s.hdfgrdbas %s/%s.hdfgrdbas" % (boundary_path, 'ena%(ens)03d', work_path, 'enf%(ens)03d')
                command_lines = doForEnsembleDict(command, member_list, current_commands=command_lines)

            else:
                print "No domain subset ..."
                command = [
                    "cp %s/%s.hdf%06d %s" % (boundary_path, 'ena%(ens)03d', args.t_ens_start, work_path),
                    "cp %s/%s.hdfgrdbas %s" % (boundary_path, 'ena%(ens)03d', work_path),
                    "cp %s/%s.hdfgrdbas %s/%s.hdfgrdbas" % (boundary_path, 'ena%(ens)03d', work_path, 'enf%(ens)03d')
                ]
                command_lines = doForEnsembleDict(command, member_list, current_commands=command_lines)

    if args.free_forecast:
        if args.subset and args.t_ens_start == exp_start:
            print "Subset the boundary conditions ..."
            ic_path = args.init_cond
            if ic_path == "":
                ic_path = boundary_path
            command_lines = generateDomainSubset(args.base_path, args.job_name, ic_path, member_list, args.n_cores_available, args.t_ens_start, args.t_ens_end, args.dt_ens_step, current_commands=command_lines, perturb_ic=False)

            command = "cp %s/%s.hdfgrdbas %s/%s.hdfgrdbas" % (boundary_path, 'ena%(ens)03d', work_path, 'enf%(ens)03d')
            command_lines = doForEnsembleDict(command, member_list, current_commands=command_lines)

        if args.piecewise:
            n_chunks = int(ceil(float(args.t_ens_end - args.t_ens_start) / args.chunk_size))
            n_chunk_start = 0

            n_chunk_start = (ens_chunk_start - args.t_ens_start) / args.chunk_size

            for n_chunk, t_chunk in enumerate(range(ens_chunk_start, args.t_ens_end, args.chunk_size)):
                print "Generating free forecast from %d to %d (chunk %d of %d) ..." % (args.t_ens_start, args.t_ens_end, n_chunk + n_chunk_start + 1, n_chunks)
                chunk_start = t_chunk
                chunk_end = t_chunk + args.chunk_size

                which_split = 'neither'
                if args.split_files and (args.split_init == 'auto' or args.split_init == 'yes'):
                    which_split = 'both'
                elif args.split_files and args.split_init == 'no':
                    which_split = 'dump'

                if chunk_end > args.t_ens_end:
                    chunk_end = args.t_ens_end

                command_lines = generateEnsembleIntegration(args.base_path, args.job_name, member_list, args.mpi_config, args.n_cores_available, chunk_start, chunk_end, args.dt_ens_step, which_split, False, current_commands=command_lines)

                req_time = args.free_fcst_req
                if args.subset and t_chunk == exp_start:
                    req_time = args.init_free_fcst_req

                submit(args.base_path, args.job_name, args.n_cores_available, command_lines, req_time, True, True)
                command_lines = {}
                
        else:
            print "Generating free forecast from %d to %d ..." % (args.t_ens_start, args.t_ens_end)
            command_lines = generateEnsembleIntegration(args.base_path, args.job_name, member_list, args.mpi_config, args.n_cores_available, args.t_ens_start, args.t_ens_end, args.dt_assim_step, args.dt_ens_step, current_commands=command_lines)
    else:
        for t_ens in xrange(args.t_ens_start, args.t_ens_end, args.dt_assim_step):
            print "Generating timestep %d ..." % t_ens

            start_time = t_ens
            end_time = t_ens + args.dt_assim_step

            if do_first_integration or t_ens > args.t_ens_start:
                if args.piecewise:
                    n_chunks = int(ceil(float(end_time - start_time) / args.chunk_size))
                    n_chunk_start = 0

                    if start_time == args.t_ens_start:
                        n_chunk_start = (ens_chunk_start - start_time) / args.chunk_size
                        start_time = ens_chunk_start

                    for n_chunk, t_chunk in enumerate(range(start_time, end_time, args.chunk_size)):
                        print "Submitting ensemble integration for timestep %d (chunk %d of %d) ..." % (t_ens, n_chunk + n_chunk_start + 1, n_chunks)

                        chunk_start = t_chunk
                        chunk_end = t_chunk + args.chunk_size
                        if chunk_end > end_time:
                            chunk_end = end_time

                        which_split = 'neither'
                        if args.split_files:
                            if chunk_start == exp_start:
                                # For the first chunk
                                if args.restart:
                                    # We're restarting
                                    if args.split_init == 'auto' or args.split_init == 'yes':
                                        which_split = 'both'
                                    elif args.split_init == 'no':
                                        which_split = 'dump'
                                else:
                                    # No restart
                                    if args.split_init == 'auto' or args.split_init == 'no':
                                        which_split = 'dump'
                                    elif args.split_init == 'yes':
                                        which_split = 'both'
                            else:
                                # Everything after the first chunk
                                which_split = 'both'
              
                        command_lines = generateEnsembleIntegration(args.base_path, args.job_name, member_list, args.mpi_config, args.n_cores_available, chunk_start, chunk_end, args.dt_ens_step, which_split, chunk_end == end_time, current_commands=command_lines)

                        if args.split_files and t_chunk == exp_start:
                            command = doForMPIConfig("cp %s/%s/%s/%s.hdfgrdbas_%s %s/%s/%s/%s.hdfgrdbas_%s" % (args.base_path, args.job_name, 'EN%(ens)s', 'ena%(ens)s', '%(x_proc)03d%(y_proc)03d', 
                                args.base_path, args.job_name, 'ENF%(ens)s', 'enf%(ens)s', '%(x_proc)03d%(y_proc)03d'), args.mpi_config)
                            command_lines = doForEnsembleDict(command, member_list, current_commands=command_lines)

                        req_time = args.fcst_req # 0:45 for 30-minute forecasts
                        if (args.subset or not args.restart) and t_chunk == exp_start:
                            req_time = args.init_fcst_req # 1:15 for 30-minute forecasts

                        submit(args.base_path, args.job_name, args.n_cores_available, command_lines, req_time, True, True)
                        command_lines = {}
                else:
                    command_lines = generateEnsembleIntegration(args.base_path, args.job_name, member_list, args.mpi_config, args.n_cores_available, start_time, end_time, args.dt_ens_step, args.split_files, True, current_commands=command_lines)

            assim_all = False
            req_time = args.assim_off_req
            if isDivisible(end_time, 3600):
                print "Assimilate all data ..."
                assim_all = True
                req_time = args.assim_on_req

            for cov_infl in args.cov_infl:
                if ':' in cov_infl:
                    time, factors = cov_infl.split(':')
                    if end_time >= int(time):
                        covariance_inflation = factors
                else:
                    covariance_inflation = cov_infl

            print "Covariance inflation for this timestep is", covariance_inflation

            if radar_data_flag[end_time]['05XP']:
                input_file = "arpsenkf.1km.05XP.input"
            else:
                input_file = "arpsenkf.input"

            assimilation_lines = generateEnKFAssimilation(args.base_path, args.job_name, args.n_ens_members, args.mpi_config, args.n_cores_available, end_time, radar_data_flag[end_time], args.split_files, covariance_inflation, assim_all, input_file_name=input_file)

            if args.piecewise:
                print "Submitting assimilation for timestep %d ..." % t_ens
                submit(args.base_path, args.job_name, nproc_x * nproc_y * 12, assimilation_lines, req_time, True, True)

    if not args.piecewise:
        submit(args.base_path, args.job_name, args.n_cores_available, command_lines, args.free_fcst_req, not args.no_submit, False)

    return

if __name__ == "__main__":
    main()
