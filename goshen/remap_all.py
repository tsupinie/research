
import subprocess
import os
import glob
import shutil

from editNamelist import editNamelistFile

def toSeconds(time_string):
    return int(time_string[:2]) * 3600 + int(time_string[2:-2]) * 60 + int(time_string[-2:])

def findTags(base_path, radar_id):
    if radar_id == "MWR05XP":
        tags = set([ f.split(".")[2] for f in glob.glob("%s/%s/*/swp.*" % (base_path, radar_id)) ])
    elif radar_id == "KRIW":
        tags = set([ f[-10:-4] for f in glob.glob("%s/%s/*" % (base_path, radar_id)) ])
    else:
        tags = set([ f[-3:] for f in glob.glob("%s/%s/swp.*.%s.*v413" % (base_path, radar_id, radar_id)) ])

    return sorted(list(tags))

def makeSweepFile(base_path, radar_id, tag):
    sweep_path = "sweeps"
    if radar_id == "MWR05XP":
        files = glob.glob("%s/%s/*/swp.%s.*" % (base_path, radar_id, tag))
        def vol(file_name):
            index = file_name.rfind("v")
            return int(file_name[(index + 1):])
        files.sort(key=vol)
    elif radar_id == "KRIW":
        files = sorted(glob.glob("%s/%s/*%s*" % (base_path, radar_id, tag)))
    else:
        files = sorted(glob.glob("%s/%s/swp.*_v%s" % (base_path, radar_id, tag)))

    sweep_file = open("%s/sweeps.%s.txt" % (sweep_path, tag), 'w')
    sweep_file.write("SOLO\n")
    for file in files:
        sweep_file.write("%s\n" % file)
    sweep_file.close()
    return files

def main():
    radar_location_KCYS = (41.15194, -104.80611)
    radar_location_KFTG = (39.78667, -104.54583)
    radar_location_KRIW = (43.06611, -108.47722)
    radar_location_05XP = (41.56150, -104.298996)

    domain_center_1km = (41.61795, -104.34843)
    domain_center_3km = (40.61795, -107.344)
    domain_std_lon = -107.344

    # Set these!  ALL OF THESE!
    radar_id = "KCYS"
    path_to_data = "raw"
    radar_location = radar_location_KCYS
    domain_center = domain_center_3km
    domain_grid_spacing = 3000
    domain_grid_size = 411
    manual_qc_flag = 'manual'

    tags = findTags(path_to_data, radar_id)

    initial_time = toSeconds("180000")

    max_pids = 8
    pids = []

    if radar_id == "KRIW":
        apparent_lat = domain_center[0]
        apparent_lon = domain_center[1]
        apparent_std_lon = domain_std_lon

        if os.path.exists("88d2arps"):
            os.unlink("88d2arps")
        os.link("/data6/tsupinie/88d2arps_KRIW", "88d2arps")
    else:
#       apparent_lat = domain_center[0]
#       apparent_lon = -domain_center[1] + 2 * (domain_center[1] - radar_location[1])
#       apparent_std_lon = -domain_std_lon + 2 * (domain_std_lon - radar_location[1])

        apparent_lat = domain_center[0]
        apparent_lon = domain_center[1]
        apparent_std_lon = domain_std_lon

        if manual_qc_flag == "manual":
            good_flag = "good"
        else:
            good_flag = "bad"

        if os.path.exists("solo2arps"):
            os.unlink("solo2arps")
        os.link("/data6/tsupinie/solo2arps_%s" % (radar_id[-4:],), "solo2arps")

    for tag in tags:
        print "Extracting tag %s ..." % tag

        files = makeSweepFile(path_to_data, radar_id, tag)
        if radar_id != "KRIW":
            time_lbound = files[0].find("swp") + 11
            time_ubound = time_lbound + 6
            time = files[0][time_lbound:time_ubound]

            if initial_time == -1:
                initial_time = toSeconds(time)
        
            offset_seconds = toSeconds(time) - initial_time
        else:
            offset_seconds = toSeconds(tag) - initial_time

        kwargs = {}
        if radar_id == "KRIW":
            kwargs['radfname'] = files[0]
            kwargs['refvarname'] = "DBZ"
            kwargs['velvarname'] = "VEL"
            src_remap_file = "radremap_88D.input"
        elif radar_id == "KCYS":
            if manual_qc_flag == "manual":
                kwargs['refvarname'] = "DZ"
                kwargs['velvarname'] = "DV"
            else:
                kwargs['refvarname'] = "DBZ"
                kwargs['velvarname'] = "VEL"
            src_remap_file = "radremap_88D.input"
        elif radar_id == "KFTG":
            if manual_qc_flag == "manual":
                kwargs['refvarname'] = "DZ"
                kwargs['velvarname'] = "DV"
            else:
                kwargs['refvarname'] = "DBZ"
                kwargs['velvarname'] = "VEL"
            src_remap_file = "radremap_88D.input"
        elif radar_id == "MWR05XP":
            src_remap_file = "radremap_MWR05XP.input"
            kwargs['refvarname'] = "DZ"
            kwargs['velvarname'] = "DV"

        if manual_qc_flag == "automated":
            kwargs['rngmin'] = 100000.0

        editNamelistFile(src_remap_file, "input/radremap.%s.input" % tag,
            nx=domain_grid_size, ny=domain_grid_size,
            dx=domain_grid_spacing, dy=domain_grid_spacing,
            radname=radar_id[-4:],
            ctrlat=apparent_lat,
            ctrlon=apparent_lon,
            trulon=apparent_std_lon,
            outtime=offset_seconds,
            **kwargs
        )

        if radar_id == "KRIW":
            pid = subprocess.Popen([ "./88d2arps" ], stdin=open("input/radremap.%s.input" % tag, 'r'), stdout=open("debug/radremap.%s.debug" % tag, 'w')).pid
        else:
            pid = subprocess.Popen([ "./solo2arps", "-fdisk", "sweeps/sweeps.%s.txt" % tag ], stdin=open("input/radremap.%s.input" % tag, 'r'), stdout=open("debug/radremap.%s.debug" % tag, 'w')).pid

        pids.append(pid)

        if len(pids) == max_pids or tag == tags[-1]:
            for pid in pids:
                os.waitpid(pid, 0)
            pids = []

    for file in glob.glob("%s.20090605.*" % radar_id):
        os.rename(file, "qc/%s/%dkm/." % (manual_qc_flag, domain_grid_spacing / 1000))

    for file in glob.glob("goshen.hdfr*"):
        os.rename(file, "hdf/%s/%dkm/." % (radar_id, domain_grid_spacing / 1000))

    return

if __name__ == "__main__":
    main()
