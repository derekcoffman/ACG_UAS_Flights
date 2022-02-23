import json
from datetime import datetime, timedelta
import pytz
import time
import sys
from netCDF4 import Dataset, num2date
import os


class Data:
    def __init__(self):

        self._data = dict()
        self._metadata = None

        self.path = "./"
        self.base_filename = None

    def load_jsonl_file(self, fname, path=None):

        filename = fname
        fname_parts = fname.split(".")
        self.base_filename = fname_parts[0]
        # print(self.base_filename)

        if path:
            # prepend path to fname
            self.path = path
            filename = path + fname
        # print(f"filename = {filename}, {self.path}")
        try:
            with open(filename, mode="r") as dfile:
                # print(f"{dfile}")
                for line in list(dfile):
                    # line = dfile.readline()
                    # print(f'line = {line}\n\n')

                    try:
                        data_line = json.loads(line)
                        self.append(data_line)
                    except json.JSONDecodeError:
                        # skip bad line
                        print("skipping badly formed json entry")
                        pass
        # try:
        #     dfile = open(filename, mode='r')
        except FileNotFoundError as e:
            print(f"File not found: {e}")
            return None

        # print("file loaded")

    def load_csl_lidar_nc_file(self, fname, path=None):
        filename = fname
        fname_parts = fname.split(".")
        self.base_filename = fname_parts[0]
        # print(self.base_filename)

        if path:
            # prepend path to fname
            self.path = path
            filename = path + fname
        # print(f"filename = {filename}, {self.path}")

        try:
            nc = Dataset(filename, "r")
        except Exception as e:
            print(f"Error opening nc file: {e}")
            return None

        axes = tuple(["TIME"])
        data_map = {axes: {}}

        for name in nc.variables:
            if name in ["version", "year", "borderHeights"]:
                pass
            else:
                data_map[axes][name] = dict()
                self._data[name] = nc[name][:]

        # make datetime
        self._data["datetime"] = []
        year = int(nc["year"][0])
        yday = nc["yDay"][:]
        dt = []
        for doy in yday:
            dt.append(doy_to_dt(year, doy))

        # dt = num2date(nc[time_var][:], units=nc[time_var].units)
        for x in dt:
            self._data["datetime"].append(dt_to_string(x))

        # print("file loaded")
        return data_map

    def load_nc_file(self, fname, path=None, time_var="time"):
        filename = fname
        fname_parts = fname.split(".")
        self.base_filename = fname_parts[0]
        # print(self.base_filename)

        if path:
            # prepend path to fname
            self.path = path
            filename = path + fname
        # print(f"filename = {filename}, {self.path}")
        # print(os.getcwd())
        try:
            nc = Dataset(filename, "r")
        except Exception as e:
            print(f"Error opening nc file: {e}")
            return None

        axes = tuple(["TIME"])
        data_map = {axes: {}}

        for name in nc.variables:
            data_map[axes][name] = dict()
            self._data[name] = nc[name][:]

        # make datetime
        self._data["datetime"] = []
        if time_var not in nc.variables:
            print(f"time_var not found in nc file")
            return None

        dt = num2date(nc[time_var][:], units=nc[time_var].units)
        for x in dt:
            self._data["datetime"].append(x.strftime(format=isofmt))

        print("file loaded")
        return data_map

    def load_psd_file(self, fname, path=None):

        # load psd header
        header_name = "./psd_met/header_ATOMIC_2020_1min_met_sea_data.txt"
        try:
            with open(header_name, mode="r") as dfile:
                header = dfile.readline()
                col_names = header.split()
        except FileNotFoundError:
            print("psd header not found!")
            return None

        axes = tuple(["TIME"])
        data_map = {axes: {}}
        for name in col_names:
            data_map[axes][name] = dict()

        filename = fname
        fname_parts = fname.split(".")
        self.base_filename = fname_parts[0]
        # print(self.base_filename)

        if path:
            # prepend path to fname
            self.path = path
            filename = path + fname
        # print(f"filename = {filename}, {self.path}")
        try:
            with open(filename, mode="r") as dfile:
                # print(f'{dfile}')
                for line in list(dfile):
                    # line = dfile.readline()
                    # print(f'line = {line}\n\n')
                    data_line = line.split()

                    for name, data in zip(col_names, data_line):
                        if name not in self._data:
                            self._data[name] = []

                        self._data[name].append(float(data))

        except FileNotFoundError as e:
            print(f"File not found: {e}")
            return None

        # make datetime
        self._data["datetime"] = []
        for i, year in enumerate(self._data["year"]):
            # year = self._data['year'][i]
            month = self._data["month"][i]
            day = self._data["day"][i]
            hour = self._data["hour"][i]
            if not hour:
                hour = 0
            minute = self._data["minute"][i]
            if not minute:
                minute = 0
            second = 0

            dt_str = f"{int(year)}-{int(month):02}-{int(day):02}T{int(hour):02}:{int(minute):02}:00Z"
            self._data["datetime"].append(dt_str)

        print("file loaded")
        return data_map

    def append(self, json_data):
        # print(f'json_data: {json_data}')
        if "DATA" in json_data:
            dt = json_data["DATA"]["DATETIME"]
            if "datetime" not in self._data:
                self._data["datetime"] = []
            self._data["datetime"].append(dt)

            measurements = json_data["DATA"]["MEASUREMENTS"]
            for meas, dat in measurements.items():
                if meas not in self._data:
                    self._data[meas] = []
                self._data[meas].append(dat["VALUE"])

        if not self._metadata and "METADATA" in json_data:
            self._metadata = json_data["METADATA"]
            # print(f'meta = {self._metadata}')

        # print(f'_data: {self._data}')

    def convert(self, output_format, output_path=None, data_map=None):
        # print("here")
        if output_path:
            self.output_path = output_path
        else:
            self.output_path = self.path

        # print(f"{self.output_path} {self.path} {output_format}")
        if output_format == "AeroSizing":
            self.write_AS_file()
        elif output_format == "msems":
            self.write_msems_file()
        elif output_format == "json":
            self.write_json_file()
        elif output_format == "igor":
            self.write_igor_file(data_map=data_map)

    def write_AS_file(self):

        with open(
            "/home/derek/Software/python/envDataSystem_analysis/utilities/AS_header.txt",
            mode="r",
        ) as h:
            header = h.read()
            # print(header)

        outfile = self.output_path + "AS" + self.base_filename + ".dat"

        print(f"output file: {outfile}")
        try:
            asfile = open(outfile, mode="w")
        except Exception as e:
            print(f"write error AS: {e}")
            return None

        asfile.write(header)

        # for name, dat in self._data.items():
        #     print(f"name: {name}")
        dt_string = self._data["datetime"]
        for i, val in enumerate(dt_string):
            # print(f'dt[{i}] = {string_to_dt(dt[i])}')
            # print(string_to_dt(val).timetuple())
            dt = string_to_dt(val).timetuple()
            year = dt.tm_year
            month = dt.tm_mon
            mday = dt.tm_mday
            HH = dt.tm_hour
            MM = dt.tm_min
            SS = dt.tm_sec

            doy = dt.tm_yday

            duration = 300

            sec_per_day = 24 * 60 * 60
            start_sec = HH * 60 * 60 + MM * 60 + SS
            start_doy = round(doy + (start_sec / sec_per_day), 5)
            stop_sec = start_sec + duration

            stop_doy = round(doy + (stop_sec / sec_per_day), 5)

            line = f"{year}\t{month}\t{mday}\t{HH}\t{MM}\t{SS}\t{duration}\t{start_doy}\t{stop_doy}"

            dmps_dp = self._data["dmps_diameter_um"][i]
            dmps_dndlogdp = self._data["dmps_dndlogdp"][i]
            line += f"\t{len(dmps_dp)}\t{len(dmps_dndlogdp)}"

            aps_dp = self._data["aps_diameter_um"][i]
            aps_dndlogdp = self._data["aps_dndlogdp"][i]
            line += f"\t{len(aps_dp)}\t{len(aps_dndlogdp)}"

            for dp in dmps_dp:
                line += f"\t{dp}"
            for dn in dmps_dndlogdp:
                line += f"\t{dn}"

            for dp in aps_dp:
                line += f"\t{dp}"
            for dn in aps_dndlogdp:
                line += f"\t{dn}"

            line += "\n"

            asfile.write(line)

        if asfile:
            asfile.close()

    def write_msems_file(self):

        num_bins = 30

        testing = False
        if testing:
            pass
        else:

            with open(
                "/home/derek/Software/python/envDataSystem_analysis/utilities/msems_header.txt",
                mode="r",
            ) as h:
                header = h.read()
                # print(header)

            col_names1 = [
                "#YY/MM/DD",
                "HR:MN:SC",
                "scan_direction",
                "actual_max_dia",
                "scan_max_volts",
                "scan_min_volts",
                "sheath_flw_avg",
                "sheath_flw_stdev",
                "mcpc_smpf_avg",
                "mcpc_smpf_stdev",
                "press_avg",
                "press_stdev",
                "temp_avg",
                "temp_stdev",
            ]
            # bin1...binN
            col_names2 = [
                "msems_errs",
                "mcpc_smpf",
                "mcpc_satf",
                "mcpc_cndt",
                "mcpc_satt",
                # "mcpc_errs"
            ]
            col_last = "mcpc_errs"
            # col_last = "\n"

            # InstrumentType_SN_YYMMDD_HHMMSS_SCAN.dat
            # eg mSEMS_1_190913_163127_SCAN.dat
            # outfile = self.output_path + 'msems_raw_' + self.base_filename + '.dat'
            outfile = (
                self.output_path
                + "mSEMS_1_"
                + self.base_filename
                + "_000000_SCAN"
                + ".dat"
            )

            print(f"output file: {outfile}")
            try:
                ofile = open(outfile, mode="w")
            except Exception as e:
                print(f"write error AS: {e}")
                return None

            ofile.write(header)

            # write column headers
            for col in col_names1:
                ofile.write(col + "\t")

            # bin counts
            bins = self._data["number_bins"][0]
            for i in range(1, bins + 1):
                ofile.write("bin" + str(i) + "\t")

            for col in col_names2:
                ofile.write(col + "\t")

            ofile.write(col_last + "\n")
            # ofile.write("\n")

            # # bin diameters
            # for i in range(1, bins + 1):
            #     ofile.write("bin_dia" + str(i) + "\t")

            # # bin concentrations
            # for i in range(1, bins + 1):
            #     ofile.write("bin_conc" + str(i) + "\t")

            # ofile.write("\n")

            # start writing data
            # for name, dat in self._data.items():
            #     print(f"name: {name}")

            dt_string = self._data["datetime"]
            for i, val in enumerate(dt_string):
                # print(f'dt[{i}] = {string_to_dt(dt[i])}')
                # print(string_to_dt(val).timetuple())
                dt = string_to_dt(val).timetuple()
                year = dt.tm_year - 2000  # msems only uses last 2 digits
                month = dt.tm_mon
                mday = dt.tm_mday
                HH = dt.tm_hour
                MM = dt.tm_min
                SS = dt.tm_sec

                # doy = dt.tm_yday

                # duration = 300

                # sec_per_day = 24*60*60
                # start_sec = HH*60*60 + MM*60 + SS
                # start_doy = round(doy + (start_sec/sec_per_day), 5)
                # stop_sec = start_sec + duration

                # stop_doy = round(doy + (stop_sec/sec_per_day), 5)

                line = f"{year:02}/{month:02}/{mday:02}\t{HH:02}:{MM:02}:{SS:02}\t"

                line += f'{self._data["scan_direction"][i]}\t'

                line += f'{self._data["actual_max_dp"][i]}\t{self._data["max_volts"][i]}\t{self._data["min_volts"][i]}\t'
                line += f'{self._data["sheath_flow_avg"][i]}\t{self._data["sheath_flow_sd"][i]}\t'
                line += f'{self._data["sample_flow_avg"][i]}\t{self._data["sample_flow_sd"][i]}\t'
                line += (
                    f'{self._data["pressure_avg"][i]}\t{self._data["pressure_sd"][i]}\t'
                )
                line += f'{self._data["temperature_avg"][i]}\t{self._data["temperature_sd"][i]}\t'
                # print(f'{self._data['bin_concentration'][i]}')

                # changed data file. this allows both
                dn = self._data["bin_concentration"][i]
                if 'bin_counts' in self._data:
                    dn = self._data["bin_counts"][i]
                    # print(f'bin_counts: {i}, {val}')

                if not dn:
                    continue

                # eff_factor = [.4, .415, .43, .445, .46, .475, .49, .505, .52, .535, .55, .7, .85]
                # for i, n in enumerate(dn):
                #     if i<= 12:
                #         line += f"{n/(eff_factor[i]*.8)}\t"
                #     else:
                #         line += f"{n}\t"
                # print(f"line: {line}")
                for n in dn:
                    line += f"{n}\t"
                # for n in self._data['bin_concentration'][i]:
                #     f'{n}\t'

                line += f'{self._data["msems_error"][i]}\t'

                line += f'{self._data["mcpc_sample_flow"][i]}\t{self._data["mcpc_saturator_flow"][i]}\t'
                line += f'{self._data["mcpc_condenser_temp"][i]}\t{45}\t'
                line += f'{self._data["mcpc_error"][i]}\n'

                # if 'diameter_nm' in self._data:
                #     dps = self._data["diameter_nm"][i]
                #     # print(f'bin_counts: {i}, {val}')
                # if not dps:
                #     continue
                # for dp in dps:
                #     line += f"{dp}\t"

                # if 'bin_concentration' in self._data:
                #     dn = self._data["bin_concentration"][i]
                #     # print(f'bin_counts: {i}, {val}')
                # if not dn:
                #     continue
                # for bini, n in enumerate(dn):
                #     if bini == (bins-1):
                #         line += f"{n}\n"
                #     else:
                #         line += f"{n}\t"

                # line += "\n"

                if "None" in line:
                    # print(line)
                    continue

                ofile.write(line)

            if ofile:
                ofile.close()

    def write_json_file(self):

        outfile = self.output_path + self.base_filename + ".json"
        with open(outfile, "w") as out:
            json.dump(self._data, out)

    def write_igor_file(self, data_map=None):

        # with open('/home/horton/derek/Software/python/envDataSystem_analysis/utilities/msems_header.txt', mode='r') as h:
        #     header = h.read()
        #     # print(header)

        # col_names1 = [
        #     '#YY/MM/DD', 'HR:MN:SC', 'scan_direction', 'actual_max_dia', 'scan_max_volts', 'scan_min_volts',
        #     'sheath_flw_avg', 'sheath_flw_stdev', 'mcpc_smpf_avg', 'mcpc_smpf_stdev	press_avg', 'press_stdev', 'temp_avg',
        #     'temp_stdev'
        # ]
        # # bin1...binN
        # col_names2 = ['sems_errs', 'mcpc_smpf', 'mcpc_satf', 'mcpc_cndt', 'mcpc_satt']
        # col_last = 'mcpc_errs'

        # InstrumentType_SN_YYMMDD_HHMMSS_SCAN.dat
        # eg mSEMS_1_190913_163127_SCAN.dat
        # outfile = self.output_path + 'msems_raw_' + self.base_filename + '.dat'

        if not data_map:
            data_map = self.parse_meta()

        # need to have a way to parse if no meta - later
        if not data_map:
            return None

        outfile = self.output_path + self.base_filename + ".itx"

        print(f"output file: {outfile}")
        try:
            ofile = open(outfile, mode="w")
        except Exception as e:
            print(f"write error igor: {e}")
            return None

        ofile.write("IGOR\n")

        # dt_mapped = False
        # for axes, meas in data_map['1d']
        # ofile.write(header)

        include_dt = False
        dt_done = False
        dt_string = self._data["datetime"]
        igor_epoch = pytz.utc.localize(datetime(1904, 1, 1, 0, 0, 0))

        for axes, meas_map in data_map.items():

            if len(axes) == 1 and axes[0] == "TIME":

                # do timeseries waves
                #
                # if axes[0] == 'TIME' and not dt_done:
                #     include_dt = True

                line = "WAVES/D\tStartDateTime"

                for name, meas in meas_map.items():
                    if name == "pops_datetime" or name == "magic_datetime" or name == "error_string":
                        pass
                    elif name == "time":
                        line += "\t" + "time_native"
                    elif name == "5v_monitor":
                        line += "\t" + "hk_5v_monitor"
                    else:
                        line += "\t" + name
                line += "\n"
                ofile.write(line)

                ofile.write("BEGIN\n")

                for i, val in enumerate(dt_string):
                    # print(f'{string_to_dt(val)} {igor_epoch} {(string_to_dt(val) - igor_epoch).total_seconds()}')
                    igor_dt = int((string_to_dt(val) - igor_epoch).total_seconds())
                    line = f"\t{igor_dt}"

                    for name, meas in meas_map.items():
                        if name == "pops_datetime" or name == "magic_datetime" or name == "error_string":
                            pass
                        else:
                            if self._data[name][i]:
                                line += f"\t{self._data[name][i]}"
                            else:
                                line += "\tNAN"
                    line += "\n"
                    ofile.write(line)

                ofile.write("END\n")
                line = f'X SetScale/P x 0,1,"", StartDateTime; SetScale y 0,0,"", StartDateTime\n'
                ofile.write(line)

                for name, meas in meas_map.items():
                    if name == "pops_datetime" or name == "magic_datetime" or name == "error_string":
                        pass
                    elif name == "time":
                        time_native = "time_native"
                        line = f'X SetScale/P x 0,1,"", {time_native}; SetScale y 0,0,"", {time_native}\n'
                        ofile.write(line)
                    elif name == "5v_monitor":
                        label = "hk_5v_monitor"
                        line = f'X SetScale/P x 0,1,"", {label}; SetScale y 0,0,"", {label}\n'
                    else:
                        line = f'X SetScale/P x 0,1,"", {name}; SetScale y 0,0,"", {name}\n'
                        ofile.write(line)

                dt_done = True

            else:
                if len(axes) == 1:
                    for name, meas in meas_map.items():
                        ofile.write(f"WAVES/D {name}\n")
                        ofile.write("BEGIN\n")

                        for i in range(0, len(self._data[name])):
                            ofile.write(f"{self._data[name]}\n")

                        ofile.write("END\n")
                        line = f'X SetScale/P x 0,1,"", {name}; SetScale y 0,0,"", {name}\n\n'
                        ofile.write(line)

                elif len(axes) == 2:

                    for name, meas in meas_map.items():
                        rows = len(self._data[name])
                        cols = len(self._data[name][0])
                        if name == "pops_datetime":
                            pass
                        else:
                            line = f"WAVES/D/N=({rows}, {cols})"
                            line += "\t" + name
                        line += "\n"
                        ofile.write(line)
                        ofile.write("BEGIN\n")

                        for r in range(0, rows):
                            line = ""
                            row = self._data[name][r]
                            for col in row:
                                line += f"\t{col}"
                            line += "\n"
                            ofile.write(line)
                        ofile.write("END\n")

                        line = f'X SetScale/P x 0 ,1,"dat", {name}; SetScale/P y 0,1,"", {name}; SetScale d 0,0,"", {name}'
                        line += "\n\n"
                        ofile.write(line)

        if ofile:
            ofile.close()

    def parse_meta(self):

        if not self._metadata:
            return None

        if "measurement_meta" in self._metadata:
            data_map = dict()
            # data_map = {
            #     '1d': dict(),
            #     '2d': dict()
            # }
            for t, group in self._metadata["measurement_meta"].items():
                for name, meas in group.items():
                    # print(f"name: {name}")
                    axes = meas["dimensions"]["axes"]
                    axes_id = tuple(axes)
                    if axes_id not in data_map:
                        data_map[axes_id] = dict()
                    data_map[axes_id][name] = meas
                    # if len(axes) == 1:
                    #     data_map['1d'][axes] = {
                    #         name: meas
                    #     }
                    # elif len(axes) == 2:
                    #     data_map['2d'][axes] = {
                    #         name: meas
                    #     }

            # print(f"data_map: {data_map}")
            return data_map
        return None

    # def pops2igor(self, fname, path='./'):

    def run(self, fname, **kw):

        conv_type = "json2itx"
        if "convert_type" in kw:
            conv_type = kw["convert_type"]

        fpath = "./"
        if "path" in kw:
            fpath = kw["path"]

        output_path = None
        if "output_path" in kw:
            output_path = kw["output_path"]

        # print(f'convert_path = {conv_type}')
        if conv_type == "jsonl2itx":
            out_fmt = "igor"
            self.load_jsonl_file(fname, path=fpath)
            self.convert(out_fmt)

        elif conv_type == "jsonl2as":
            out_fmt = "AeroSizing"
            self.load_jsonl_file(fname, path=fpath)
            self.convert(out_fmt)

        elif conv_type == "jsonl2msems":
            out_fmt = "msems"
            self.load_jsonl_file(fname, path=fpath)
            self.convert(out_fmt, output_path=output_path)

        elif conv_type == "jsonl2json":
            out_fmt = "json"
            self.load_jsonl_file(fname, path=fpath)
            self.convert(out_fmt)

        elif conv_type == "psd2itx":
            out_fmt = "igor"
            data_map = self.load_psd_file(fname, path=fpath)
            if data_map:
                self.convert(out_fmt, data_map=data_map)

        elif conv_type == "nc2itx":
            out_fmt = "igor"

            if "time_var" in kw:
                time_var = kw["time_var"]
                data_map = self.load_nc_file(fname, path=fpath, time_var=time_var,)
            else:
                data_map = self.load_nc_file(fname, path=fpath,)

            if data_map:
                self.convert(out_fmt, data_map=data_map)

        elif conv_type == "csl_lidar2itx":
            out_fmt = 'igor'
            data_map = self.load_csl_lidar_nc_file(fname, path=fpath)
            # print(f'{data_map}')
            if data_map:
                self.convert(out_fmt, data_map=data_map)

        # # d.load_jsonl_file('2020-01-08.jsonl', path='./')
        # fpath = None
        # if 'path' in kw:
        #     fpath = kw['path']
        # # fpath = './utilities/'
        # d.load_jsonl_file(fname, path=fpath)

        # # out_fmt = 'igor_text'
        # # fpath = './utilities/'
        # out_fmt = 'igor'
        # if 'output_format' in kw:
        #     out_fmt = kw['output_format']
        # d.convert(out_fmt)


isofmt = "%Y-%m-%dT%H:%M:%SZ"


def string_to_dt(dt_string):
    dt = datetime.strptime(dt_string, isofmt)
    # print(type(dt))
    # print(type(pytz.utc.localize(dt)))
    return pytz.utc.localize(dt)


def dt_to_string(dt):
    utc = pytz.utc.localize(dt)
    return utc.strftime(isofmt)


def doy_to_dt(year, doy):
    dt = datetime(year, 1, 1) + timedelta(doy - 1)
    return dt


if __name__ == "__main__":

    # ---- uncomment from here
    # if len(sys.argv) < 2:
    #     print("specify file name to convert")
    #     sys.exit()

    # print(f"args : {sys.argv}")

    # fname = sys.argv[1]
    # kw = {}
    # # if len(sys.argv) > 2:
    # # kw = {}
    # for arg in sys.argv[2:]:
    #     print(arg)
    #     parts = arg.split("=")
    #     kw[parts[0]] = parts[1]
    # ----- to here for normal use

    path = "/home/derek/Data/envDataSystem/from_cloudbase/UIServer/cloudy1.acg.pmel.noaa.gov/cloudysky"

    # for testing
    # fname = 'ATOMIC_2020_1min_met_sea_data_01_20_jd020.txt'
    # kw = {'convert_type': 'psd2itx', 'path': './utilities/psd_met/'}
    # fpath = './utilities/psd_met/'
    #
    # fname = 'ATOMIC_CN_v1.nc'
    # kw = {'convert_type': 'nc2itx', 'path': './utilities/', 'time_var': 'time'}
    # print(f"fname={fname}, kw={kw}")

    d = Data()

    # kw options:
    #   path=<path to input/output files> default='./'
    #   convert_type=<conversion routine to run> default='json2itx'
    #       options: json2as, psd2itx
    d.run(fname, **kw)
    d.load_jsonl_file
    # d.load_jsonl_file('2020-01-08.jsonl', path='./')
    # fpath = None
    # if 'path' in kw:
    #     fpath = kw['path']
    # # fpath = './utilities/'
    # d.load_jsonl_file(fname, path=fpath)

    # # out_fmt = 'igor_text'
    # # fpath = './utilities/'
    # out_fmt = 'igor'
    # if 'output_format' in kw:
    #     out_fmt = kw['output_format']
    # d.convert(out_fmt)
