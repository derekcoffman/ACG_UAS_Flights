import glob
from operator import inv
import os
import json
import struct

# from tracemalloc import start
from click import clear
import numpy as np
import xarray as xr
import pandas as pd

# import hvplot
# import hvplot.pandas
# import hvplot.xarray
from datetime import datetime, timedelta

from convert import Data, dt_to_string, string_to_dt


class DataEvent:

    isofmt = "%Y-%m-%dT%H:%M:%SZ"
    isofmt_day = "%Y-%m-%d"

    def __init__(self, event_id, config=None, auto_create=True) -> None:

        self._data = dict()
        self.event_id = event_id
        self.handle_config(config)

        if config and auto_create:
            self.create_event()

    def handle_config(self, config):

        if config is None:
            return

        self._data["kind"] = config["kind"]
        self._data["metadata"] = dict()
        for name, meta in config["metadata"].items():
            self._data["metadata"][name] = meta

        try:
            event = config["events"][self.event_id]
        except KeyError:
            print("bad event_id")
            return

        start_dt = datetime.strptime(event["start_time"], DataEvent.isofmt)
        end_dt = datetime.strptime(event["end_time"], DataEvent.isofmt)
        self._data["times"] = dict()
        self._data["times"]["start"] = {"time": event["start_time"], "dt": start_dt}
        self._data["times"]["end"] = {"time": event["end_time"], "dt": end_dt}

        # self._data["times"]["duration"] = (end_dt-start_dt).seconds
        self._data["datasystem"] = dict()
        for ds in event["datasystems"]:
            self._data["datasystem"][ds] = {"spec": config["datasystems"][ds]}

        # for ds, spec in config["datasystem"].items():
        #     self._data["datasystem"][ds] = {"spec": spec}

    def get_data(self):
        return self._data

    def merge(self, tb=1, ds_list=None):

        # return self.get_data()["datasystem"][datasystem][controller][instrument]
        print(f"ds_list = {ds_list}" )
        if ds_list is None:
            ds_list = []
            for dsystem_name, dsystem in self._data["datasystem"].items():
                for cont_name, controller in dsystem.items():
                    if cont_name == "spec":
                        continue
                    for dsname, ds in controller.items():
                        # print(dsystem_name, cont_name, dsname)
                        if ds.attrs["timebase"] == tb:
                            ds_list.append(dsname)

        print(f"ds_list2 = {ds_list}" )
        start_dt = self._data["times"]["start"]["dt"]
        end_dt = self._data["times"]["end"]["dt"]
        freq = f"{tb}s"
        tol = None
        # tol = "500ms"
        if tb > 1:
            tol = f"{tb/2}s"
        dr = pd.date_range(start_dt, end_dt, freq=freq, closed="left")
        ds_out = xr.Dataset()
        ds_out_list = []
        for dsystem_name, dsystem in self._data["datasystem"].items():
            print(dsystem_name)
            for cont_name, controller in dsystem.items():
                if cont_name == "spec":
                    continue
                print(cont_name)
                for dsname, ds in controller.items():
                    print(f"dsname = {dsname}")
                    if dsname in ds_list:
                        print(f"dsname in list = {dsname}")
                        # print(dsname)
                        ds_1 = ds.sel(time=~ds.get_index("time").duplicated())
                        if dsname == "navigation":
                            print(f"ds_1 = {ds_1}")
                        ds_2 = (
                            ds_1.reindex(time=sorted(ds_1.time.values))
                            .resample(time=freq)
                            .mean(keep_attrs=True)
                        )
                        ds_out_list.append(
                            # ds.reindex(time=sorted(ds.time.values)).sel(time=~ds.get_index("time").duplicated()).sel(time=slice(start_dt, end_dt)).reindex({"time": dr}, method="nearest", tolerance=tol)
                            ds_2.sel(time=slice(start_dt, end_dt)).reindex(
                                {"time": dr}, method="nearest"  # , tolerance=tol
                            )
                        )
                        if dsname == "navigation":
                            print(ds_2)
        ds_out = xr.combine_by_coords(ds_out_list, combine_attrs="drop")

        # # add extra time arrays: mid, end
        tb_ms = int(tb * 1000)  # use ms
        tb_mid = int(tb_ms / 2)
        # print(tb)
        ds_out["time_mid"] = ds_out.time + np.timedelta64(tb_mid, "ms")
        ds_out["time_end"] = ds_out.time + np.timedelta64(tb_ms, "ms")
        ds_out["duration"] = xr.DataArray(
            np.full(ds_out.dims["time"], tb),
            coords={"time": ds_out.time.values},
            dims=["time"],
        )
        meta = self._data["metadata"]
        for key, val in meta.items():
            ds_out.attrs[key] = val
        ds_out.attrs["event_id"] = self.event_id
        ds_out.attrs["timebase"] = tb
        return ds_out

    def process_msems(self, instrument, datasystem=None, controller=None):
        ds = self.get_dataset(instrument, datasystem=datasystem, controller=controller)
        process_msems(ds)

    def process_cdp(self, instrument, datasystem=None, controller=None, **kw):
        ds = self.get_dataset(instrument, datasystem=datasystem, controller=controller)
        process_cdp(ds, **kw)

    def process_pops(self, instrument, datasystem=None, controller=None, **kw):
        ds = self.get_dataset(instrument, datasystem=datasystem, controller=controller)
        ds = process_pops(ds, **kw)
        self.set_dataset(ds, instrument, datasystem=datasystem, controller=controller)

    def process_payload(self, instrument, datasystem=None, controller=None, **kw):
        ds = self.get_dataset(instrument, datasystem=datasystem, controller=controller)
        ds = process_payload(ds, **kw)
        # self.set_dataset(ds, instrument, datasystem=datasystem, controller=controller)

    def process_piccolo(self, instrument, datasystem=None, controller=None):
        ds = self.get_dataset(instrument, datasystem=datasystem, controller=controller)
        ds = process_piccolo(ds)
        self.set_dataset(ds, instrument, datasystem=datasystem, controller=controller)

    def get_dataset(self, instrument, datasystem=None, controller=None):
        # eventually, this should look for the first instrument in the data structure, but for now, require all
        if datasystem is None or controller is None:
            return None

        try:
            return self.get_data()["datasystem"][datasystem][controller][instrument]
        except KeyError:
            return None

    def set_dataset(self, ds, instrument, datasystem=None, controller=None):
        # eventually, this should look for the first instrument in the data structure, but for now, require all
        if datasystem is None or controller is None:
            return None

        try:
            self.get_data()["datasystem"][datasystem][controller][instrument] = ds
        except KeyError:
            return None

    def get_start_dt(self):
        return self._data["times"]["start"]["dt"]

    def get_end_dt(self):
        return self._data["times"]["end"]["dt"]

    def create_event(self) -> dict:

        # get list of day filenames that should be loaded
        # file_list = []
        # for dt in pd.date_range(start=self.get_start_dt(),end=self.get_end_dt(), freq="D").to_pydatetime().tolist():
        #     file_list.append(f"{dt.strftime('%Y-%m-%d')}.jsonl")

        # file_list = get_file_list(start_time=self.get_start_dt(), end_time=self.get_end_dt())

        # load and create xr.Dataset for each instrument in each datasystem
        for ds_name, datasystem in self._data["datasystem"].items():
            spec = datasystem["spec"]
            for cont_name, controller in spec.items():
                path = controller["base_path"]
                for name, cfg in controller["instruments"].items():
                    data_format = "envds"
                    try:
                        data_format = cfg["format"]
                    except KeyError:
                        pass

                    data = dict()
                    if data_format in ["envdsys", "envds"]:
                        file_list = get_file_list(
                            start_time=self.get_start_dt(), end_time=self.get_end_dt()
                        )

                        for dayfile in file_list:
                            fname = os.path.join(path, cont_name, name, dayfile)
                            data = self.load_envds_datafile(fname, data=data)

                    elif data_format in ["piccolo-log"]:
                        base_file_list = get_piccolo_file_list(
                            start_time=self.get_start_dt(), end_time=self.get_end_dt()
                        )
                        file_list = []
                        for bf in base_file_list:
                            # print(f"{bf}, {path}")
                            pf = glob.glob(os.path.join(path, bf))
                            # pf = glob.glob(bf)
                            for f in pf:
                                file_list.append(f)
                            # print(f"filelist: {file_list}")
                        for dayfile in file_list:
                            fname = os.path.join(path, dayfile)
                            fname = dayfile
                            print(fname)
                            data = self.load_piccolo_datafile(
                                fname,
                                # start_dt=self.get_start_dt(),
                                # end_dt=self.get_end_dt(),
                                data=data,
                            )

                    elif data_format in ["clear-payload-dat"]:
                        # base_file_list = get_clear_dat_file_list(
                        #     start_time=self.get_start_dt(), end_time=self.get_end_dt()
                        # )
                        file_list = []
                        # for bf in base_file_list:
                        # bf = "payload/*.DAT"
                        bf = "*.DAT"
                        pf = glob.glob(os.path.join(path, bf))
                        for f in pf:
                            file_list.append(f)
                        for dayfile in file_list:
                            # fname = os.path.join(path, dayfile)
                            fname = dayfile
                            data = self.load_clear_payload_dat_datafile(
                                fname,
                                # start_dt=self.get_start_dt(),
                                # end_dt=self.get_end_dt(),
                                data=data,
                            )

                    elif data_format in ["pops-bin"]:
                        # base_file_list = get_clear_dat_file_list(
                        #     start_time=self.get_start_dt(), end_time=self.get_end_dt()
                        # )
                        try:
                            cal_file = cfg["calibration_file"]
                            caldf = pd.read_csv(cal_file)
                            calds = xr.Dataset(caldf)
                            # calds = calds.rename_dims({"dim_0": "cal_bins"}).drop(
                            #     "dim_0"
                            # )
                            calds = (
                                calds.rename({"Amplitude": "amplitude"})
                                .set_coords("amplitude")
                                .swap_dims({"dim_0": "amplitude"})
                                .drop("dim_0")
                            )

                        except (KeyError, FileNotFoundError):
                            print(
                                "Missing POPS calibration file, not able to load pops-bin data"
                            )
                            continue

                        try:
                            bin_count = cfg["bin_count"]
                        except KeyError:
                            bin_count = 26

                        file_list = []
                        # for bf in base_file_list:
                        bf = "*.b"
                        pf = glob.glob(os.path.join(path, bf))
                        for f in pf:
                            file_list.append(f)
                        for dayfile in file_list:
                            # fname = os.path.join(path, dayfile)
                            fname = dayfile
                            data = self.load_pops_bin_datafile(
                                fname,
                                calds,
                                bin_count,
                                # start_dt=self.get_start_dt(),
                                # end_dt=self.get_end_dt(),
                                data=data,
                            )

                    else:
                        return None

                    vars = None
                    if "variables" in cfg:
                        vars = cfg["variables"]

                    ds = self.create_dataset(data, dims=cfg["dims"], vars=vars)
                    # ds = self.prepare_dataset(ds)
                    if ds:
                        ds.attrs["start_dt"] = self.get_start_dt()
                        ds.attrs["end_dt"] = self.get_end_dt()
                        # add timebase to ds - default is 1s
                        try:
                            ds.attrs["timebase"] = cfg["timebase"]
                        except KeyError:
                            ds.attrs["timebase"] = 1

                        # # add extra time arrays: mid, end
                        # tb = int(ds.attrs["timebase"] * 1000) # use ms
                        # tb_mid = int(tb/2)
                        # # print(tb)
                        # ds.coords["time_mid"] = ds.time + np.timedelta64(tb_mid, "ms")
                        # ds.coords["time_end"] = ds.time + np.timedelta64(tb, "ms")
                        # ds["duration"] = xr.DataArray(np.full(ds.dims["time"], tb), coords={"time": ds.time.values}, dims=["time"], name="duration")

                        if cont_name not in datasystem:
                            datasystem[cont_name] = dict()
                        if "variables" in cfg:
                            # datasystem[cont_name][name] = ds[cfg["variables"]]
                            var_list = []
                            for v in cfg["variables"]:
                                var_list.append(v[1])
                            datasystem[cont_name][name] = ds[var_list]
                        else:
                            datasystem[cont_name][name] = ds
                        pass

    def load_envds_datafile(self, file_name, data=None):
        if data is None:
            data = dict()
        # data["extra_coords"] = dict()
        need_meta = True
        metdata = None
        # with open("/home/derek/Data/envDataSystem/from_cloudbase/UIServer/cloudy1.acg.pmel.noaa.gov/cloudysky/uas_cloudy/fast_trh/2021-06-17.jsonl") as f:

        start_dt = self.get_start_dt()
        end_dt = self.get_end_dt()
        try:
            with open(file_name) as f:
                # line = f.readline()
                for line in f:
                    # print(line)
                    try:
                        entry = json.loads(line)
                    except json.JSONDecodeError:
                        continue

                    if "time" not in data:
                        data["time"] = []

                    dup_index = None
                    try:
                        dt = datetime.strptime(
                            entry["DATA"]["DATETIME"], "%Y-%m-%dT%H:%M:%SZ"
                        )
                        if dt < start_dt or dt >= end_dt:
                            continue

                        if dt in data["time"]:
                            dup_index = data["time"].index(dt)
                            # print(f"len of data[time] = {len(data['time'])}")
                            data["time"].pop(dup_index)
                            # print(f"\tlen of data[time] = {len(data['time'])}")

                        data["time"].append(dt)
                    except KeyError:
                        continue

                    for key, dat in entry["DATA"]["MEASUREMENTS"].items():
                        if key not in data:
                            data[key] = []

                        try:
                            if dup_index:
                                data[key].pop(dup_index)
                            record = entry["DATA"]["MEASUREMENTS"][key]["VALUE"]
                            data[key].append(record)
                        except KeyError:
                            data[key].append(None)
                        except TypeError:
                            print(f"{key}, {dat}")
                            data[key].append(None)

                    if "METADATA" in entry:
                        data["metadata"] = entry["METADATA"]
        except FileNotFoundError:
            print(f"load evnds: file not found: {file_name}")

        return data

    def load_piccolo_datafile(self, file_name, data=None):
        if data is None:
            data = dict()
        # data["extra_coords"] = dict()
        need_meta = True
        metdata = None
        # with open("/home/derek/Data/envDataSystem/from_cloudbase/UIServer/cloudy1.acg.pmel.noaa.gov/cloudysky/uas_cloudy/fast_trh/2021-06-17.jsonl") as f:

        isofmt = "%Y-%m-%dT%H:%M:%S.%fZ"
        ds_start_dt = self.get_start_dt()
        ds_end_dt = self.get_end_dt()
        varlist = []
        print(file_name)
        try:
            with open(file_name) as f:
                # get header and variable list
                line = f.readline().rstrip()
                params = line.split(" ")
                for par in params:
                    # print(f"par: {par}, {par.split('>[')}, {par.split('>[')[0].replace('<','')}")
                    parts = par.split(">")
                    # print(f"{parts}, {len(parts)}")
                    name = parts[0].replace("<", "")
                    varlist.append(name)
                    data[name] = []
                    if len(parts) > 1 and parts[1]:
                        units = parts[1].replace("[", "").replace("]", "")
                        if "metadata" not in data:
                            data["metadata"] = {"measurement_meta": {"primary": {}}}
                        data["metadata"]["measurement_meta"]["primary"][name] = {
                            "units": units
                        }
                        #     "measurement_meta": {"primary": {name: {"units": units}}}
                        # }

                for line in f:
                    parts = line.rstrip().split(" ")

                    for name, val in zip(varlist, parts):
                        # sanity check for date
                        # print(parts[1], start_dt)
                        if (
                            int(parts[1]) < ds_start_dt.year
                            or int(parts[1]) > ds_end_dt.year
                        ):
                            continue
                        try:
                            data[name].append(float(val))
                        except ValueError:
                            data[name].append(None)

                base_ms = data["Clock"][0]
                # fulldate = fulldate + datetime.timedelta(milliseconds=500)
                yr = data["Year"][0]
                mo = data["Month"][0]
                da = data["Day"][0]
                hr = data["Hours"][0]
                mi = data["Minutes"][0]
                se = data["Seconds"][0]
                millis = (se - int(se)) * 1000.0
                isofmt = "%Y-%m-%dT%H:%M:%S.%fZ"
                dtstr = f"{int(yr):02}-{int(mo):02}-{int(da):02}T{int(hr):02}:{int(mi):02}:{int(se):02}.000000Z"
                # start_dt = datetime.strptime(dtstr, isofmt).replace(microsecond=0)
                data["time"] = []
                # for (clock, yr,mo,da,hr,mi,se) in zip(data["Clock"]["data"],data["Year"]["data"],data["Month"]["data"],data["Day"]["data"],data["Hours"]["data"],data["Minutes"]["data"],data["Seconds"]["data"]):
                for i, clck in enumerate(data["Clock"]):
                    # dtstr = f"{int(yr)}-{int(mo)}-{int(da)}T{int(hr)}:{int(mi)}:{int(se)}Z"
                    # ms = int(clock)
                    # dt = (start_dt + timedelta(milliseconds=int(ms)))#.replace(microsecond=0)
                    yr = data["Year"][i]
                    mo = data["Month"][i]
                    da = data["Day"][i]
                    hr = data["Hours"][i]
                    mi = data["Minutes"][i]
                    se = data["Seconds"][i]
                    if se < 0:
                        se = 0
                    elif se > 59.99:
                        se = 59.99
                    millis = (se - int(se)) * 1000.0
                    dtstr = f"{int(yr):02}-{int(mo):02}-{int(da):02}T{int(hr):02}:{int(mi):02}:{int(se):02}.{int(millis):<06}Z"
                    dt = datetime.strptime(dtstr, isofmt)
                    # dup_index = 99999
                    # while dup_index > 0:
                    #     # dup_index = -1
                    #     if dt in data["time"]:
                    #         dup_index = data["time"].index(dt)
                    #         for key, dat in data.items():
                    #             if key == "metadata":
                    #                 continue
                    #             data[key].pop(dup_index)
                    #     else:
                    #         dup_index = -1

                    # data["time"].append((start_dt + timedelta(milliseconds=int(ms))).replace(microsecond=0))
                    data["time"].append(dt)
            # res = [idx for idx, val in enumerate(data["time"]) if val in data["time"][:idx]]
        except FileNotFoundError:
            print(f"load piccolo: file not found: {file_name}")

        return data

    def load_clear_payload_dat_datafile(self, file_name, data=None):
        need_meta = True
        metdata = None

        raw_data = dict()
        isofmt = "%Y-%m-%dT%H:%M:%S.%fZ"
        ds_start_dt = self.get_start_dt()
        ds_end_dt = self.get_end_dt()
        varlist = []

        curr_date = ds_start_dt.strftime("%Y-%m-%d")
        # print(file_name)
        try:
            with open(file_name) as f:

                current_time = None
                for line in f:

                    # split on whitespace
                    parts = line.split()

                    for p in parts:
                        if ":" in p:
                            # if "time" not in raw_data:
                            #     raw_data["time"] = []
                            # raw_data["time"].append(p)
                            current_time = f"{curr_date}T{p}Z"
                            if current_time not in raw_data:
                                raw_data[current_time] = dict()
                        elif "=" in p:
                            try:
                                name, val = p.split("=")
                            except ValueError:
                                print(f"split payload error: {p}")
                                raw_data[current_time][name] = None
                                continue
                            # if name not in raw_data[current_time]:
                            #     raw_data[current_time][name] = []
                            if name == "POP":
                                pops_data = val.split(",")
                                raw_data[current_time][name] = [
                                    x for x in pops_data[1:]
                                ]
                            elif name == "Date":
                                curr_date = val
                                raw_data[current_time][name] = val
                            elif name == "AT":
                                try:
                                    raw_data[current_time][name] = float(val)
                                except ValueError:
                                    print(f"AT(string) = {val}")
                                    raw_data[current_time][name] = None
                            else:
                                raw_data[current_time][name] = val

            par_list = []
            if data is None or len(data) == 0:
                data = dict()
                # build key list
                parameters = ["time"]
                for name, record in raw_data.items():
                    # if name == "time":
                    #     continue
                    for parname in record.keys():
                        if parname not in data:
                            data[parname] = []
                            par_list.append(parname)
                        # if parname not in parameters:
                        #     parameters.append(parname)
            else:
                par_list = list(data.keys())
                if "time" in par_list:
                    time_index = par_list.index("time")
                    par_list.pop(time_index)

            # print(par_list)
            for dtstr, record in raw_data.items():

                # dt = datetime.strptime(f"{curr_date}T{dtstr}Z", "%Y-%m-%dT%H:%M:%SZ")
                dt = datetime.strptime(dtstr, "%Y-%m-%dT%H:%M:%SZ")
                if "time" not in data:
                    data["time"] = []

                dup_index = None
                try:
                    if dt in data["time"]:
                        dup_index = data["time"].index(dt)
                        # print(f"time dup_index = {dup_index}")
                        # print(f"len of data[time] = {len(data['time'])}")
                        data["time"].pop(dup_index)
                        # print(f"\tlen of data[time] = {len(data['time'])}")
                    if dt >= ds_start_dt and dt <= ds_end_dt:
                        data["time"].append(dt)
                    else:
                        continue
                except KeyError:
                    continue

                for key in par_list:
                    if key not in data:
                        print(f"make array: {key}")
                        data[key] = []

                    try:
                        if dup_index:
                            data[key].pop(dup_index)
                            # print(f"var[{key}] dup_index = {dup_index}")
                        rec = record[key]
                        data[key].append(rec)
                    except KeyError:
                        data[key].append(None)
                    except TypeError:
                        # print(f"{key}, {dat}")
                        data[key].append(None)
        except FileNotFoundError:
            print(f"load payload: file not found: {file_name}")

        return data

    def load_pops_bin_datafile(self, file_name, cal, bin_count=26, data=None):
        if data is None:
            data = dict()
        # data["extra_coords"] = dict()
        need_meta = True
        metdata = None
        # with open("/home/derek/Data/envDataSystem/from_cloudbase/UIServer/cloudy1.acg.pmel.noaa.gov/cloudysky/uas_cloudy/fast_trh/2021-06-17.jsonl") as f:

        bin_bounds = []
        bin_bounds.append(140)
        dlogDp = np.log10(3000 / 140) / bin_count
        # dlogDp = 0.05
        pops_dp_um = []
        for i in range(1, 27):
            bin_bounds.append(round(bin_bounds[i - 1] * (10 ** dlogDp), 3))
            pops_dp_um.append(
                round(np.sqrt(bin_bounds[i - 1] * bin_bounds[i]) / 1000.0, 3)
            )

        ds_start_dt = self.get_start_dt()
        ds_end_dt = self.get_end_dt()

        try:
            with open(file_name, "rb") as f:
                while True:
                    buf = f.read(12)
                    if not buf or len(buf) != 12:
                        break
                    num_part, tstamp = struct.unpack("<Ld", buf)
                    try:
                        dt = datetime.utcfromtimestamp(np.floor(tstamp))
                    except OverflowError:
                        dt = datetime.utcfromtimestamp(0)

                    # print(num_part, tstamp)
                    if num_part > 0:  # and dt >= ds_start_dt and dt < ds_end_dt:
                        data_peaks = []
                        for i in range(0, num_part):
                            buf = f.read(12)
                            if len(buf) != 12:
                                break
                            pk, _, _ = struct.unpack("<3L", buf)
                            data_peaks.append(pk)

                        if dt >= ds_start_dt and dt < ds_end_dt:
                            amplitude = xr.DataArray(
                                np.log10(data_peaks), dims=["count"]
                            )
                            data_dp = cal.interp(
                                # cal_bins=amplitude,
                                amplitude=amplitude,
                                method="cubic",
                                kwargs={"fill_value": "extrapolate"},
                            ).Size.to_dataset()

                            try:
                                if "time" not in data:
                                    data["time"] = []
                                if "bin_counts" not in data:
                                    data["bin_counts"] = []
                                if "diameter_um" not in data:
                                    data["diameter_um"] = []

                                data["bin_counts"].append(
                                    list(
                                        data_dp.groupby_bins(
                                            group=data_dp.Size, bins=bin_bounds
                                        )
                                        .count()
                                        .Size.values
                                    )
                                )
                                data["diameter_um"].append(pops_dp_um)
                                data["time"].append(dt)
                                # hist_data.append(list(data_dp.groupby_bins(group=data_dp.Size, bins=bin_bounds).count().Size.values))
                                # dt_data.append(datetime.fromtimestamp(np.floor(tstamp)))
                            except ValueError:
                                pass
                            # if len(hist_data) > 10:
                            #     break

            # ('bin_counts', 'pops_bin_counts', 'clear_pops_bin_counts'),
            # ('diameter_um', 'pops_dp_um_2d', 'cloudy_msems_diameter_um'),

            # raw_data = dict()
            # isofmt = "%Y-%m-%dT%H:%M:%S.%fZ"
            # ds_start_dt = self.get_start_dt()
            # ds_end_dt = self.get_end_dt()
            # varlist = []

            # curr_date = ds_start_dt.strftime("%Y-%m-%d")
            # print(file_name)
            # with open(file_name) as f:

            #     current_time = None
            #     for line in f:

            #         # split on whitespace
            #         parts = line.split()

            #         for p in parts:
            #             if ":" in p:
            #                 # if "time" not in raw_data:
            #                 #     raw_data["time"] = []
            #                 # raw_data["time"].append(p)
            #                 current_time = f"{curr_date}T{p}Z"
            #                 if current_time not in raw_data:
            #                     raw_data[current_time] = dict()
            #             elif "=" in p:
            #                 name,val = p.split("=")
            #                 # if name not in raw_data[current_time]:
            #                 #     raw_data[current_time][name] = []
            #                 if name == "POPS":
            #                     pops_data = val.split(",")
            #                     raw_data[current_time][name] = [x for x in pops_data[1:]]
            #                 elif name == "Date":
            #                     curr_date = val
            #                     raw_data[current_time][name] = val
            #                 else:
            #                     raw_data[current_time][name] = val

            # par_list = []
            # if data is None or len(data) == 0:
            #     data = dict()
            #     # build key list
            #     parameters = ["time"]
            #     for name, record in raw_data.items():
            #         # if name == "time":
            #         #     continue
            #         for parname in record.keys():
            #             if parname not in data:
            #                 data[parname] = []
            #                 par_list.append(parname)
            #             # if parname not in parameters:
            #             #     parameters.append(parname)

            # for dtstr, record in raw_data.items():

            #     # dt = datetime.strptime(f"{curr_date}T{dtstr}Z", "%Y-%m-%dT%H:%M:%SZ")
            #     dt = datetime.strptime(dtstr, "%Y-%m-%dT%H:%M:%SZ")
            #     if "time" not in data:
            #         data["time"] = []

            #     dup_index =  None
            #     try:
            #         if dt in data["time"]:
            #             dup_index = data["time"].index(dt)
            #             # print(f"len of data[time] = {len(data['time'])}")
            #             data["time"].pop(dup_index)
            #             # print(f"\tlen of data[time] = {len(data['time'])}")

            #         data["time"].append(dt)
            #     except KeyError:
            #         continue

            #     for key in par_list:
            #         if key not in data:
            #             data[key] = []

            #         try:
            #             if dup_index:
            #                 data[key].pop(dup_index)
            #             rec = record[key]
            #             data[key].append(rec)
            #         except KeyError:
            #             data[key].append(None)
            #         except TypeError:
            #             # print(f"{key}, {dat}")
            #             data[key].append(None)
        except FileNotFoundError:
            print(f"load pops: file not found: {file_name}")

        return data

    def get_dimensions(self, data, parameter):
        dims = []
        if "metadata" not in data:
            dims = ["time"]

        for key, val in data["metadata"].items():
            if key == "measurement_meta":
                for mtype, mdata in val.items():
                    if parameter in mdata and "dimensions" in mdata[parameter]:
                        dims = [
                            d.lower() for d in mdata[parameter]["dimensions"]["axes"]
                        ]
                        return dims

        return dims

    def create_dataset(self, data, dims=None, vars=None):
        # print(data)
        # coords = self.parse_coordinates(data)
        # print(len(data["time"]), len(data["CONCN"]))
        if data is None:
            return None

        ds = xr.Dataset()
        # add coordinates
        coords = ["time"]
        try:
            ds.coords["time"] = data["time"]
        except KeyError:
            return None

        if dims is None:
            dims = ["time"]
        # for name, coord in data["extra_coords"].items():
        #     par = coord["parameter"]
        #     if isinstance(data[par][0], list):
        #         ds.coords[name] = data[par][0]
        #     else:
        #         ds.coords[name] = data[par]
        #     coords.append(name)
        # if len(ds)
        # for name, coord in coords.items():
        #     if name == "time":
        #         ds.coords[name] = data[coord["parameter"]]

        if vars:
            var_map = dict()
            for v in vars:
                var_map[v[0]] = {"name": v[1], "long_name": v[2]}
            # print(var_map)
            for key, entry in data.items():
                if key not in coords and key != "metadata":
                    # if key == "diameter_um":
                    #     print('here')
                    # if key != "time" and key != "metadata":
                    # get dimensions
                    # dims = self.get_dimensions(data, key)
                    long_name = ""
                    if key in var_map:
                        long_name = var_map[key]["long_name"]
                        key = var_map[key]["name"]

                        if isinstance(entry[0], list):
                            ds[key] = (dims, entry)
                            try:
                                if any(isinstance(x, str) for x in entry[0]):
                                    ds[key] = (
                                        dims,
                                        np.asarray(entry, dtype=np.float64),
                                    )
                                else:
                                    ds[key] = (dims, np.asarray(entry))
                            except (TypeError, ValueError):
                                ds[key] = (dims, np.asarray(entry))
                        else:
                            if any(isinstance(x, str) for x in entry):
                                try:
                                    ds[key] = (
                                        ["time"],
                                        np.asarray(entry, dtype=np.float64),
                                    )
                                except (TypeError, ValueError):
                                    ds[key] = (["time"], np.asarray(entry))
                            else:
                                ds[key] = (["time"], np.asarray(entry))

                    # don't use long_name at the moment
                    # if long_name:
                    #     ds[key].attrs["long_name"] = long_name

            if "metadata" in data:
                # print(data["metadata"])
                for key, val in data["metadata"].items():
                    # print(key)
                    # if key in var_map:
                    #     key = var_map[key]

                    if key == "measurement_meta":  # add parameter metadata
                        for mtype, mdata in val.items():
                            for param, meta in mdata.items():
                                if param in var_map:
                                    param = var_map[param]["name"]
                                for att, attval in meta.items():
                                    try:
                                        ds[param].attrs[att] = attval
                                    except KeyError:
                                        pass
                    else:
                        try:
                            if key == "plot_meta":
                                continue
                            ds.attrs[key] = val
                        except KeyError:
                            pass
        # print(len(ds.time))
        return ds

    def to_float(self, val) -> float:
        if val is None:
            return None  # is not None: # and type(val) is str:
        return float(val)
        #     try:
        #         return np.float64(val)
        #     except TypeError:
        #         raise TypeError
        # return None

    def prepare_dataset(self, ds):
        # print(ds)
        # sort on time in case things got wonky
        # tsorted = xr.DataArray(sorted(ds.time.values))
        ds_sort = ds.reindex(time=sorted(ds.time.values))

        # remove any duplicates
        ds_final = ds_sort.sel(time=~ds_sort.get_index("time").duplicated())

        # select time slice of interest and put on common time base
        ds_all = (
            ds_final.sel(time=slice(self.get_start_dt(), self.get_end_dt()))
            .resample(time="S")
            .fillna(value=None)
        )
        return ds_all


def get_file_list(start_time, end_time):
    # get list of day filenames that should be loaded
    file_list = []
    for dt in (
        pd.date_range(start=start_time, end=end_time, freq="D").to_pydatetime().tolist()
    ):
        file_list.append(f"{dt.strftime('%Y-%m-%d')}.jsonl")
    return file_list


def get_piccolo_file_list(start_time, end_time):
    # get list of day filenames that should be loaded
    file_list = []
    for dt in (
        pd.date_range(start=start_time, end=end_time, freq="D").to_pydatetime().tolist()
    ):
        file_list.append(f"Piccolo_*_{dt.strftime('%Y-%m-%d')}_*-*-*.log")
    # print(file_list)
    return file_list


def get_inverted_file_list(start_time, end_time):
    # get list of day filenames that should be loaded
    file_list = []
    for dt in (
        pd.date_range(start=start_time, end=end_time, freq="D").to_pydatetime().tolist()
    ):
        file_list.append(f"mSEMS_1_{dt.strftime('%Y-%m-%d')}_000000_INVERTED.dat")
    return file_list


def save_for_invert_msems(instrument, event_id, datasystem, controller, config):

    try:
        event = config["events"][event_id]
    except KeyError:
        print("bad event_id")
        return

    # start_dt = datetime.strptime(event["start_time"], DataEvent.isofmt)
    # end_dt = datetime.strptime(event["end_time"], DataEvent.isofmt)

    file_list = get_file_list(event["start_time"], event["end_time"])
    try:
        controller_config = config["datasystems"][datasystem][controller]
        inv_config = controller_config["instruments"][instrument]
        inv_path = os.path.join(
            inv_config["invert_path"], ""
        )  # add trailing slash if necessary
        base_path = os.path.join(
            controller_config["base_path"], controller, instrument, ""
        )

        # create inversion path if necessary
        try:
            os.makedirs(inv_path)
        except FileExistsError:
            pass
        except PermissionError:
            print(f"Do not have permission to create inversion directory: {inv_path}")
            return

        d = Data()
        for dayfile in file_list:
            # path = "/home/derek/Data/envDataSystem/from_cloudbase/UIServer/cloudy1.acg.pmel.noaa.gov/cloudysky/uas_cloudy/cloudy_msems/"
            # fname = os.path.join(,dayfile)
            # print(f"inv path: {inv_path}")
            d.run(
                fname=dayfile,
                path=base_path,
                convert_type="jsonl2msems",
                output_path=inv_path,
            )
    except KeyError:
        print(f"could not find msems")
        return None
    # for ds_name, datasystem in config["datasystem"].items():
    #     for cont_name, controller in datasystem.items():
    #         if instrument in controller_config["instruments"]:
    #             inv_config = controller["instruments"][instrument]
    #             inv_path = inv_config["invert_path"]
    #             base_path = os.path.join(controller["base_path"], cont_name, instrument, "")

    #             d = Data()
    #             for dayfile in file_list:
    #                 path = "/home/derek/Data/envDataSystem/from_cloudbase/UIServer/cloudy1.acg.pmel.noaa.gov/cloudysky/uas_cloudy/cloudy_msems/"
    #                 # fname = os.path.join(,dayfile)
    #                 d.run(fname=dayfile, path=base_path, convert_type="jsonl2msems", output_path=inv_path)


# def load_inverted_msems(name, config, data=None):
def load_inverted_msems(
    instrument, event_id, datasystem, controller, config, data=None
):

    try:
        event = config["events"][event_id]
    except KeyError:
        print("bad event_id")
        return

    # start_dt = datetime.strptime(event["start_time"], DataEvent.isofmt)
    # end_dt = datetime.strptime(event["end_time"], DataEvent.isofmt)

    file_list = get_inverted_file_list(event["start_time"], event["end_time"])
    try:
        controller_config = config["datasystems"][datasystem][controller]
        inv_config = controller_config["instruments"][instrument]
        inv_path = os.path.join(
            inv_config["invert_path"], ""
        )  # add trailing slash if necessary
        base_path = os.path.join(
            controller_config["base_path"], controller, instrument, ""
        )

        dt = []
        conc = []
        for invfile in file_list:
            invname = os.path.join(inv_path, invfile)
            # print(invname)
            with open(invname) as f:
                next(f)  # skip header
                for line in f:
                    parts = line.split()
                    dt.append(
                        datetime.strptime(
                            " ".join([parts[0], parts[1]]), "%y/%m/%d %H:%M:%S"
                        )
                    )
                    conc.append(
                        [float(x) for x in parts[5:]]
                    )  # serial numbers are missing from data so header is wrong - conc start at col 5
        da = xr.DataArray(conc, coords={"time": dt}, dims=inv_config["dims"])
        if data:
            try:
                msems = data["datasystem"][datasystem][controller][instrument]
                msems["msems_dNdlogDp"] = da
            except KeyError:
                pass
        else:
            return da

    except KeyError:
        print("error in load inverted file")
        return None

    # file_list = get_inverted_file_list(config["start_time"], config["end_time"])
    # for ds_name, datasystem in config["datasystem"].items():
    #     for cont_name, controller in datasystem.items():
    #         if name in controller["instruments"]:
    #             inv_config = controller["instruments"][name]
    #             inv_path = inv_config["invert_path"]
    #             # base_path = os.path.join(controller["base_path"], cont_name, name, "")
    #             dt = []
    #             conc = []
    #             for invfile in file_list:
    #                 invname = os.path.join(inv_path, invfile)
    #                 with open(invname) as f:
    #                     next(f) # skip header
    #                     for line in f:
    #                         parts = line.split()
    #                         dt.append(datetime.strptime(" ".join([parts[0], parts[1]]), "%y/%m/%d %H:%M:%S"))
    #                         conc.append([float(x) for x in parts[5:]]) #serial numbers are missing from data so header is wrong - conc start at col 5
    #             da = xr.DataArray(conc, coords={"time": dt}, dims=inv_config["dims"])
    #             if data:
    #                 try:
    #                     msems = data["datasystem"][ds_name][cont_name][name]
    #                     msems["dNdlogDp"] = da
    #                 except KeyError:
    #                     pass
    #             else:
    #                 return da


def calc_dp_bounds(dp_um):
    dlogDp = []
    # dp_um = msems.msems_dp_um.values
    dlogDp = np.ones(len(dp_um))
    dp_bounds_um = np.ones(len(dp_um) + 1)
    # print(len(dp_bounds_um))
    for i in range(1, len(dp_um)):
        dp_bounds_um[i] = np.sqrt(dp_um[i] * dp_um[i - 1])
        # print(i, dp_bounds_um[i])
    dp_bounds_um[0] = dp_um[0] - (dp_bounds_um[1] - dp_um[0])
    dp_bounds_um[len(dp_um)] = dp_um[-1] + (dp_um[-1] - dp_bounds_um[len(dp_um) - 1])
    for i in range(0, len(dp_um)):
        dlogDp[i] = np.log10(dp_bounds_um[i + 1] / dp_bounds_um[i])

    return dlogDp, dp_bounds_um


def process_msems(msems):
    # calc dlogDp, dp_bounds, dp_um, intNSV,
    msems.coords["msems_dp_um"] = (["msems_bins"], msems.msems_dp_um_2d.values[0])
    msems.msems_dp_um.attrs["units"] = "um"
    (dlogDp, dp_bounds_um) = calc_dp_bounds(msems.msems_dp_um.values)
    # dlogDp = []
    # dp_um = msems.msems_dp_um.values
    # dlogDp = np.ones(len(dp_um))
    # dp_bounds_um = np.ones(len(dp_um)+1)
    # # print(len(dp_bounds_um))
    # for i in range(1,len(dp_um)):
    #     dp_bounds_um[i] = np.sqrt(dp_um[i]*dp_um[i-1])
    #     # print(i, dp_bounds_um[i])
    # dp_bounds_um[0] = dp_um[0] - (dp_bounds_um[1]-dp_um[0])
    # dp_bounds_um[len(dp_um)] = dp_um[-1] + (dp_um[-1] - dp_bounds_um[len(dp_um)-1])
    # for i in range(0,len(dp_um)):
    #     dlogDp[i] = np.log(dp_bounds_um[i+1]/dp_bounds_um[i])

    msems["msems_dlogDp"] = (["msems_bins"], dlogDp)
    msems.coords["msems_dp_bounds_um"] = (["msems_bin_bounds"], dp_bounds_um)
    msems.msems_dp_bounds_um.attrs["units"] = "um"

    msems["msems_dN"] = (
        ["time", "msems_bins"],
        msems.msems_dNdlogDp.data * msems.msems_dlogDp.data,
    )
    msems["msems_dSdlogDp"] = (
        ["time", "msems_bins"],
        msems.msems_dNdlogDp.data
        * 4.0
        * np.pi
        * np.power(msems.msems_dp_um.data / 2, 2),
    )
    msems["msems_dS"] = (
        ["time", "msems_bins"],
        msems.msems_dSdlogDp.data * msems.msems_dlogDp.data,
    )
    msems["msems_dVdlogDp"] = (
        ["time", "msems_bins"],
        msems.msems_dNdlogDp.data
        * 4.0
        / 3.0
        * np.pi
        * np.power(msems.msems_dp_um.data / 2, 3),
    )
    msems["msems_dV"] = (
        ["time", "msems_bins"],
        msems.msems_dVdlogDp.data * msems.msems_dlogDp.data,
    )

    msems["msems_intN"] = msems.msems_dN.sum(dim="msems_bins")
    msems.msems_intN.attrs["units"] = "cm-3"
    msems["msems_intS"] = msems.msems_dS.sum(dim="msems_bins")
    msems.msems_intS.attrs["units"] = "um2/cm3"
    msems["msems_intV"] = msems.msems_dV.sum(dim="msems_bins")
    msems.msems_intV.attrs["units"] = "um3/cm3"


def calc_cdp_dN(cdp, speed=None):
    if speed is None:
        speed = xr.DataArray(np.full(cdp.dims["time"], 20.0), dims=["time"])

    # print(cdp.cdp_bin_counts)
    # cdp["cdp_dN"] = cdp.cdp_bin_counts * (0.24 * (1 / 100)) * speed * 100
    # fixed 12 Aug 2022
    cdp["cdp_dN"] = cdp.cdp_bin_counts / ((0.24 * (1 / 100)) * speed * 100)
    # print(cdp.cdp_dN)
    # cdp.cdp_dN.attrs["units"] = "m/s"


def process_cdp(cdp, speed=None, recode=False):

    if recode:  # fix bad binary decoding of earlier data
        # skip if already recoded
        if "recode_status" in cdp.attrs and cdp.attrs["recode_status"]:
            pass

        else:
            import json
            import struct
            import csv

            counts = cdp["cdp_bin_counts"]
            new_counts = []

            for scan in counts:
                new_scan = []
                for bin in scan:
                    if bin is not None:
                        rebin = struct.pack("<I", int(bin))
                        newbin = struct.pack(">2h", *struct.unpack("<2h", rebin))
                        new_scan.append(struct.unpack(">I", newbin)[0])
                    else:
                        new_scan.append(None)
                new_counts.append(new_scan)

            cdp["cdp_bin_counts"] = xr.DataArray(new_counts, dims=["time", "cdp_bins"])
            cdp.attrs["recode_status"] = True

    # if speed is None:
    #     speed = 20 # m/s

    # print(f"speed = {speed}")
    cdp.coords["cdp_dp_um"] = (["cdp_bins"], cdp.cdp_dp_um_2d.values[0])
    cdp.cdp_dp_um.attrs["units"] = "um"

    (dlogDp, dp_bounds_um) = calc_dp_bounds(cdp.cdp_dp_um.values)

    cdp["cdp_dlogDp"] = (["cdp_bins"], dlogDp)
    cdp.coords["cdp_dp_bounds_um"] = (["cdp_bin_bounds"], dp_bounds_um)
    cdp.cdp_dp_bounds_um.attrs["units"] = "um"

    # print(cdp.cdp_bin_counts)
    calc_cdp_dN(cdp, speed)

    cdp["cdp_dNdlogDp"] = (
        ["time", "cdp_bins"],
        cdp.cdp_dN.data / cdp.cdp_dlogDp.data,
    )

    cdp["cdp_dSdlogDp"] = (
        ["time", "cdp_bins"],
        cdp.cdp_dNdlogDp.data * 4.0 * np.pi * np.power(cdp.cdp_dp_um.data / 2, 2),
    )
    cdp["cdp_dS"] = (
        ["time", "cdp_bins"],
        cdp.cdp_dSdlogDp.data * cdp.cdp_dlogDp.data,
    )
    cdp["cdp_dVdlogDp"] = (
        ["time", "cdp_bins"],
        cdp.cdp_dNdlogDp.data * 4.0 / 3.0 * np.pi * np.power(cdp.cdp_dp_um.data / 2, 3),
    )
    cdp["cdp_dV"] = (
        ["time", "cdp_bins"],
        cdp.cdp_dVdlogDp.data * cdp.cdp_dlogDp.data,
    )

    cdp["cdp_intN"] = cdp.cdp_dN.sum(dim="cdp_bins")
    cdp.cdp_intN.attrs["units"] = "cm-3"
    cdp["cdp_intS"] = cdp.cdp_dS.sum(dim="cdp_bins")
    cdp.cdp_intS.attrs["units"] = "um2/cm3"
    cdp["cdp_intV"] = cdp.cdp_dV.sum(dim="cdp_bins")
    cdp.cdp_intV.attrs["units"] = "um3/cm3"

    cdp["cdp_lwc"] = cdp.cdp_intV / 100 ** 3
    cdp.cdp_lwc.attrs["units"] = "g/m3"


def resample_pops(pops, resample_tb=30):

    time_raw = list(pops.time.values)
    bin_counts_raw = list(pops.pops_bin_counts)
    dp_um_2d_raw = list(pops.pops_dp_um_2d)

    freq = f"{resample_tb}s"
    pops = pops.resample(time=freq).sum(keep_attrs=True)
    pops.coords["time_raw"] = (["time_raw"], time_raw)
    pops["pops_bin_counts_raw"] = (["time_raw", "pops_bins"], bin_counts_raw)
    pops["pops_dp_um_2d"] = (
        ["time", "pops_bins"],
        np.full((pops.dims["time"], pops.dims["pops_bins"]), dp_um_2d_raw[0]),
    )
    pops.attrs["resample_tb"] = resample_tb
    pops.attrs["timebase"] = resample_tb
    return pops


def process_pops(pops, flow_rate=None, firstDp=0.140, lastDp=3.0, avg_tb=30):

    # if recode: # fix bad binary decoding of earlier data
    #     import json
    #     import struct
    #     import csv

    #     counts = cdp["cdp_bin_counts"]
    #     new_counts = []

    #     for scan in counts:
    #         new_scan = []
    #         for bin in scan:
    #             if bin is not None:
    #                 rebin = struct.pack("<I", int(bin))
    #                 newbin = struct.pack(">2h", *struct.unpack("<2h", rebin))
    #                 new_scan.append(struct.unpack(">I", newbin)[0])
    #             else:
    #                 new_scan.append(None)
    #         new_counts.append(new_scan)

    # cdp["cdp_bin_counts"] = xr.DataArray(new_counts, dims=["time", "cdp_bins"])

    if "time_raw" not in pops.dims:
        pops = resample_pops(pops, avg_tb)
        # time_raw = list(pops.time.values)
        # bin_counts_raw = list(pops.pops_bin_counts)
        # freq = f"{resample_tb}s"
        # pops.resample(time=freq).sum()
        # pops.coords["time_raw"] = (["time_raw"], time_raw)
        # pops.coords["bin_counts_raw"] = (["time_raw", bin_counts_raw], bin_counts_raw)
        # pops.attrs["resample_tb"] = resample_tb

    elif "resample_tb" in pops.attrs and pops.attrs["resample_tb"] != avg_tb:
        time_raw = list(pops.time_raw.values)
        bin_counts_raw = list(pops.pops_bin_counts_raw)
        # freq = f"{resample_tb}s"
        pops = xr.Dataset()
        pops.coords["time"] = (["time"], time_raw)
        pops["pops_bin_counts"] = (["time", "pops_bins"], bin_counts_raw)
        pops = resample_pops(pops, avg_tb)

    pops.coords["pops_dp_um"] = (["pops_bins"], pops.pops_dp_um_2d.values[0])
    pops.pops_dp_um.attrs["units"] = "um"

    dlogDp = np.log10(lastDp / firstDp) / pops.dims["pops_bins"]
    dp_bounds_um = [firstDp]
    for i in range(1, (pops.dims["pops_bins"] + 1)):
        dp_bounds_um.append(dp_bounds_um[i - 1] * 10 ** dlogDp)

    pops["pops_dlogDp"] = (["pops_bins"], np.full(pops.dims["pops_bins"], dlogDp))
    pops.coords["pops_dp_bounds_um"] = (["pops_bin_bounds"], dp_bounds_um)
    pops.pops_dp_bounds_um.attrs["units"] = "um"

    if flow_rate is None:
        flow_rate = xr.DataArray(np.full([pops.dims["time"]], 3.0), dims=["time"])
    # flow = xr.DataArray(np.full([pops.dims["time"], pops.dims["pops_bins"]], flow_rate), dims=["time", "pops_bins"])
    flow_rate = flow_rate.where(flow_rate > 0)
    # print(flow)
    # pops["pops_dN"] = (["time", "pops_bins"], pops.pops_bin_counts / (avg_tb*flow_rate))
    pops["pops_dN"] = pops.pops_bin_counts / (avg_tb * flow_rate)
    # print(pops.pops_dN)
    pops["pops_dNdlogDp"] = (
        ["time", "pops_bins"],
        pops.pops_dN.data / pops.pops_dlogDp.data,
    )
    pops["pops_dSdlogDp"] = (
        ["time", "pops_bins"],
        pops.pops_dNdlogDp.data * 4.0 * np.pi * np.power(pops.pops_dp_um.data / 2, 2),
    )
    pops["pops_dS"] = (
        ["time", "pops_bins"],
        pops.pops_dSdlogDp.data * pops.pops_dlogDp.data,
    )
    pops["pops_dVdlogDp"] = (
        ["time", "pops_bins"],
        pops.pops_dNdlogDp.data
        * 4.0
        / 3.0
        * np.pi
        * np.power(pops.pops_dp_um.data / 2, 3),
    )
    pops["pops_dV"] = (
        ["time", "pops_bins"],
        pops.pops_dVdlogDp.data * pops.pops_dlogDp.data,
    )

    pops["pops_intN"] = pops.pops_dN.sum(dim="pops_bins")
    pops.pops_intN.attrs["units"] = "cm-3"
    pops["pops_intS"] = pops.pops_dS.sum(dim="pops_bins")
    pops.pops_intS.attrs["units"] = "um2/cm-3"
    pops["pops_intV"] = pops.pops_dV.sum(dim="pops_bins")
    pops.pops_intV.attrs["units"] = "um3/cm-3"

    return pops


def process_payload(payload, **kwargs):
    process_payload_abs(payload, **kwargs)

    return payload


def process_payload_abs(payload, ref_init=None, avg_time=30):
    if payload is None:
        return None

    # process abs/psap
    if ref_init is None:
        ref_init = {
            payload.time.values[0]: {
                "init_450": (payload.ABSBA.values[0] / payload.ABSBB.values[0]),
                "init_525": (payload.ABSGA.values[0] / payload.ABSGB.values[0]),
                "init_624": (payload.ABSRA.values[0] / payload.ABSRB.values[0]),
            },
        }

    # print(ref_init)
    tr_init_450 = xr.DataArray(
        np.full(payload.dims["time"], None, dtype=np.float64), dims=["time"]
    )
    tr_init_525 = xr.DataArray(
        np.full(payload.dims["time"], None, dtype=np.float64), dims=["time"]
    )
    tr_init_624 = xr.DataArray(
        np.full(payload.dims["time"], None, dtype=np.float64), dims=["time"]
    )

    # spot_size = 19.3  # mm
    spot_size = 18.96
    # interpolate ABSFL to remove nans
    # payload["ABSFL"] = payload.ABSFL.interpolate_na(dim="time")
    flow_rate = (payload.ABSFL.interpolate_na(dim="time") / 60.0).rename("flow_rate")

    wavelength = [450, 525, 624]
    # wavelength = [525]
    for wl in wavelength:
        tr_init = xr.DataArray(
            np.full(payload.dims["time"], None, dtype=np.float64), dims=["time"]
        )
        if wl == 450:
            color = "B"
        elif wl == 525:
            color = "G"
        elif wl == 624:
            color = "R"
        else:
            continue

        tr_init_name = f"tr_init_{wl}"
        tr_name = f"tr_{wl}"
        bap_name = f"bap_{wl}"

        for dt, entry in ref_init.items():
            # init_name = f"tr_init_{wl}"
            if isinstance(dt, str):
                dt = np.datetime64(dt)

            tr_init = xr.where(
                payload.time >= dt,
                entry[f"init_{wl}"],
                tr_init.where(payload.time < dt),
            ).rename("tr_init")

        np.timedelta64()
        payload[tr_init_name] = tr_init
        tr = xr.DataArray(
            (payload[f"ABS{color}A"] / payload[f"ABS{color}B"] / tr_init),
            coords=[payload.time.values],
            dims=["time"],
            name=tr_name,
        )
        ba = xr.DataArray(
            np.full(payload.dims["time"], np.nan, dtype=np.float64),
            coords=[payload.time.values],
            dims=["time"],
            name=bap_name,
        )

        for i in range(0, payload.dims["time"]):
            dti = payload.time[i].values
            delta2x = np.timedelta64((2 * avg_time), "s")
            delta1x = np.timedelta64(avg_time, "s")
            dr2 = pd.date_range(
                start=(dti - delta2x), end=(dti - delta1x), freq="S", closed="left"
            )
            dr1 = pd.date_range(
                start=(dti - delta1x), end=(dti), freq="S", closed="left"
            )
            try:
                tr2 = tr.sel(time=dr2).mean(skipna=True)
                tr1 = tr.sel(time=dr1).mean(skipna=True)
                trdiff = np.log(tr2 / tr1)
                flow_avg = (
                    flow_rate.sel(time=dr1)
                    .where(tr.sel(time=dr1).notnull())
                    .mean(skipna=True)
                )
                td = (
                    payload.time.sel(time=dr1).where(tr.sel(time=dr1).notnull()).mean()
                    - payload.time.sel(time=dr2)
                    .where(tr.sel(time=dr2).notnull())
                    .mean()
                )
                new_avg_time = (td / (np.timedelta64(1, "s"))).values
                ba[i] = trdiff * spot_size / flow_avg / new_avg_time
                ba[i] *= 0.873 / (1.317 * tr1 + 0.866)
                ba[i] *= 1e6
                # print(tr.sel(time=dr1))
                # print(trdiff[i], tr.sel(time=dr1).mean(), tr.sel(time=dr2).mean())
            except (ValueError, KeyError):
                ba[i] = np.nan

        # # tr60 = np.log(tr.rolling(time=60).mean())
        # tr60 = tr.rolling(time=60).mean()
        # # tr60 = 10**lntr60
        # # trdiff = tr60.diff("time")
        # lntr60 = np.log(tr60)
        # trdiff = lntr60.diff("time")*(-1.0)
        # # trdiff = 10**lntrdiff
        # # ba = trdiff * 19.3 / (1600.0 / 60.0) / 60 * 1e6 * (1.317 * tr60[1:] + 0.866)
        # ba = (
        #     trdiff * spot_size / flow_rate[1:] / 60 * 1e6 * (1.317 * tr60[1:] + 0.866)
        # ).rename(bap_name)
        # ba.reindex(time=tr60.time.values)
        payload[tr_name] = tr
        payload[bap_name] = ba
        payload[bap_name].attrs["units"] = "Mm-1"

    return payload


def process_piccolo(piccolo):
    if piccolo is None:
        return None

    # clean up duplicates and sort
    ds_1 = piccolo.sel(time=~piccolo.get_index("time").duplicated())
    piccolo = (
        ds_1.reindex(time=sorted(ds_1.time.values))
        # .resample(time=freq)
        # .mean(keep_attrs=True)
    )

    if "start_dt" in piccolo.attrs and "end_dt" in piccolo.attrs:
        piccolo = piccolo.sel(
            time=slice(piccolo.attrs["start_dt"], piccolo.attrs["end_dt"])
        )

    # convert radians to degrees
    for name, var in piccolo.items():
        if var.attrs["units"] == "rad":
            piccolo[name].values = var * 180 / np.pi
            if name == "latitude":
                piccolo[name].attrs["units"] = "degrees_north"
            elif name == "longitude":
                piccolo[name].attrs["units"] = "degree_east"
            else:
                piccolo[name].attrs["units"] = "degree"

    # add height variables in feet for convenience
    heights = ["altitude", "pressure_altitude", "height_agl"]
    for h in heights:
        try:
            piccolo[f"{h}_ft"] = xr.DataArray(
                piccolo[h].values * 3.28084, dims=["time"]
            )
            piccolo[f"{h}_ft"].attrs["units"] = "feet"
        except KeyError:
            pass

    return piccolo


def save_flight_itx(ds):
    try:
        payload_id = ds.attrs["payload_id"]
        project = ds.attrs["project"]
        flight_id = ds.attrs["flight_id"]
        event_id = ds.attrs["event_id"]
        tb = ds.attrs["timebase"]
    except KeyError:
        return

    # create output path if necessary
    try:
        os.makedirs(os.path.join("data", "igor"))
    except FileExistsError:
        pass
    except PermissionError:
        print(
            f"Do not have permission to create inversion directory: {os.path.join('data', 'igor')}"
        )
        return

    fname = os.path.join(
        "data", "igor", f"{project}_{flight_id}_{payload_id}_{event_id}_{tb}s.itx"
    )
    # print(fname)
    igor_epoch = (
        np.datetime64("1904-01-01T00:00:00").astype("datetime64[s]").astype("int")
    )

    # print(ds.coords)

    out = open(fname, "w")
    # print(out)

    out.write("IGOR\n")

    # name, coord = ds.items()
    # print(f"name={name}")
    for name, coord in ds.reset_coords().coords.items():
        name = name.replace("-", "_")
        if name in ["time", "time_mid", "time_end", "time_raw"]:
            vname = f"{event_id}_{name}"
            line = f"WAVES/O/D\t{vname}\n"
            out.write(line)
            out.write("BEGIN\n")
            line = ""
            for ts in coord.values.astype("datetime64[s]").astype("int"):
                line = "\t".join([line, str(ts - igor_epoch)])
            out.write(f"{line}\n")
            out.write("END\n")

        else:
            line = f"WAVES/O/D\t{name}\n"
            out.write(line)
            out.write("BEGIN\n")
            line = ""
            for val in coord.values:
                line = "\t".join([line, str(val)])
            out.write(f"{line}\n")
            out.write("END\n")

    for name, param in ds.reset_coords().items():
        name = name.replace("-", "_")
        if name in ["time", "time_mid", "time_end"]:
            vname = f"{event_id}_{name}"
            line = f"WAVES/O/D\t{vname}\n"
            out.write(line)
            out.write("BEGIN\n")
            line = ""
            for ts in param.values.astype("datetime64[s]").astype("int"):
                line = "\t".join([line, str(ts - igor_epoch)])
            out.write(f"{line}\n")
            out.write("END\n")
        elif len(param.shape) > 1:
            line = f"WAVES/O/D/N={param.shape}\t{name}\n"
            out.write(line)
            out.write("BEGIN\n")
            line = ""
            for row in param.values:
                line = ""
                for val in row:
                    line = "\t".join([line, str(val)])
                out.write(f"{line}\n")
            out.write("END\n")
            line = f'X SetScale/P x {ds.time[0].values.astype("datetime64[s]").astype("int")-igor_epoch} ,{tb},"dat", {name}\n'
            out.write(line)
        else:
            line = f"WAVES/O/D\t{name}\n"
            out.write(line)
            out.write("BEGIN\n")
            line = ""
            for val in param.values:
                line = "\t".join([line, str(val)])
            out.write(f"{line}\n")
            out.write("END\n")
            line = f'X SetScale/P x {ds.time[0].values.astype("datetime64[s]").astype("int")-igor_epoch} ,{tb},"dat", {name}\n'
            out.write(line)

    line = f"X make/o/D/n=(numpnts({event_id}_time)+1) {event_id}_time_im; {event_id}_time_im[0]={event_id}_time[0]; {event_id}_time_im[1,]={event_id}_time_im[p-1]+duration[p-1]\n"
    out.write(line)
    out.close()


if __name__ == "__main__":

    config = {
        # "version": "?",
        "kind": "Flight",  # Cruise, Cruise-Leg
        "metadata": {
            "project": "PMEL VP2022",
            "platform": "FVR-55-001",
            "experiment_id": "VP2022_SimFlight_001",
        },
        "start_time": "2022-02-01T21:00:00Z",
        "end_time": "2022-02-01T22:00:00Z",
        "datasystem": {
            "CloudySky": {
                "uas_cloudy": {
                    "base_path": "/home/derek/Data/envDataSystem/from_cloudbase/UIServer/cloudy1.acg.pmel.noaa.gov/cloudysky",
                    "instruments": [
                        "ambient_trh",
                        "cloudy_cdp",
                        "cloudy_msems",
                        "fast_trh",
                        "msems_sheath_trh",
                    ],
                }
            },
            "GroundStation": {
                "uasground": {
                    "base_path": "/home/derek/Data/envDataSystem/from_cloudbase",
                    "instruments": ["amcpc_cn"],
                }
            },
        },
    }
    config2 = {
        "kind": "CloudySkyFlight",
        "metadata": {
            "project": "Tucson2021",
            "platform": "FVR-55",
            "flight_id": "Tucson2021_05_CloudySky",
        },
        "events": {
            "preflight": {},
            "flight": {
                "start_time": "2021-03-26T20:45:00Z",
                "end_time": "2021-03-26T21:55:00Z",
                "datasystems": ["CloudySky", "AutoPilot"],
            },
            "postflight": {},
        },
        "datasystems": {
            "CloudySky": {
                "uas_cloudy": {
                    "base_path": "/home/derek/Data/envDataSystem/from_cloudbase/UIServer/cloudy1-cloudysky",
                    "instruments": {
                        "ambient_trh": {
                            "format": "envdsys",
                            "timebase": 1,
                            "dims": ["time"],
                            "variables": [
                                ("temperature", "ambient_T", "ambient_temperature"),
                                (
                                    "relative_humidity",
                                    "ambient_rh",
                                    "ambient_relative_humidity",
                                ),
                            ],
                        },
                        "msems_sheath_trh": {
                            "format": "envdsys",
                            "timebase": 1,
                            "dims": ["time"],
                            "variables": [
                                (
                                    "temperature",
                                    "msems_sh_T",
                                    "msems_sheath_temperature",
                                ),
                                (
                                    "relative_humidity",
                                    "msems_sh_rh",
                                    "msems_sheath_relative_humidity",
                                ),
                            ],
                        },
                        "cloudy_cdp": {
                            "format": "envdsys",
                            "timebase": 1,
                            "dims": ["time", "cdp_bins"],
                            "process_options": {"type": "cdp"},
                            "variables": [
                                (
                                    "bin_counts",
                                    "cdp_bin_counts",
                                    "cloudy_cdp_bin_counts",
                                ),
                                (
                                    "diameter_um",
                                    "cdp_dp_um_2d",
                                    "cloudy_cdp_diameter_um_2d",
                                ),
                                (
                                    "integral_counts",
                                    "cdp_intN",
                                    "cloudy_cdp_integral_counts",
                                ),
                                (
                                    "laser_current",
                                    "cdp_laser_current",
                                    "cloudy_cdp_laser_current",
                                ),
                                (
                                    "wingboard_temperature",
                                    "cdp_wingboard_T",
                                    "cloudy_cdp_wingboard_temperature",
                                ),
                                (
                                    "laser_temperature",
                                    "cdp_laser_T",
                                    "cloudy_cdp_laser_temperature",
                                ),
                                (
                                    "control_board_temperature",
                                    "cdp_control_board_T",
                                    "cloudy_cdp_control_board_temperature",
                                ),
                            ],
                        },
                        "cloudy_msems": {
                            "format": "envdsys",
                            "timebase": 30,
                            "resample": False,
                            "dims": ["time", "msems_bins"],
                            "invert_path": "/home/derek/Data/UAS/msems_inversion/",
                            "process_options": {"type": "msems"},
                            "variables": [
                                (
                                    "bin_counts",
                                    "msems_bin_counts",
                                    "cloudy_msems_bin_counts",
                                ),
                                (
                                    "diameter_um",
                                    "msems_dp_um_2d",
                                    "cloudy_msems_diameter_um",
                                ),
                                (
                                    "scan_direction",
                                    "msems_scan_direction",
                                    "cloudy_msems_scan_direction",
                                ),
                                (
                                    "sheath_flow_avg",
                                    "msems_sheath_flow_avg",
                                    "cloudy_msems_sheath_flow_avg",
                                ),
                                (
                                    "sheath_flow_sd",
                                    "msems_sheath_flow_sd",
                                    "cloudy_msems_sheath_flow_sd",
                                ),
                                (
                                    "sample_flow_avg",
                                    "msems_sample_flow_avg",
                                    "cloudy_msems_sample_flow_avg",
                                ),
                                (
                                    "sample_flow_sd",
                                    "msems_sample_flow_sd",
                                    "cloudy_msems_sample_flow_sd",
                                ),
                                (
                                    "pressure_avg",
                                    "msems_pressure_avg",
                                    "cloudy_msems_pressure_avg",
                                ),
                                (
                                    "pressure_sd",
                                    "msems_pressure_sd",
                                    "cloudy_msems_pressure_sd",
                                ),
                                (
                                    "temperature_avg",
                                    "msems_T_avg",
                                    "cloudy_msems_temperature_avg",
                                ),
                                (
                                    "temperature_sd",
                                    "msems_T_sd",
                                    "cloudy_msems_temprature_sd",
                                ),
                                (
                                    "mcpc_sample_flow",
                                    "mcpc_sample_flow",
                                    "cloudy_msems_mcpc_sample_flow",
                                ),
                                (
                                    "mcpc_saturator_flow",
                                    "mcpc_saturator_flow",
                                    "cloudy_msems_mcpc_saturator_flow",
                                ),
                                (
                                    "mcpc_condenser_temp",
                                    "mcpc_condenser_T",
                                    "cloudy_msems_mcpc_condenser_temp",
                                ),
                                ("bin_time", "msems_bin_time", "cloudy_msems_bin_time"),
                                (
                                    "scan_type",
                                    "msems_scan_type",
                                    "cloudy_msems_scan_type",
                                ),
                                (
                                    "plumbing_time",
                                    "msems_plumbing_time",
                                    "cloudy_msems_plumbing_time",
                                ),
                            ],
                        },
                    },
                }
            },
            "GroundStation": {
                "uasground": {
                    "base_path": "/home/derek/Data/envDataSystem/from_cloudbase",
                    "instruments": {
                        "amcpc_cn": {
                            "format": "envdsys",
                            "timebase": 1,
                            "dims": ["time"],
                            "variables": [
                                ("concentration", "ground_cn", "amcpc_cn_concentration")
                            ],
                        }
                    },
                }
            },
            "AutoPilot": {
                "navigation": {
                    "base_path": "/home/derek/Data/UAS/autopilot/piccolo",
                    "instruments": {
                        "piccolo": {
                            "format": "piccolo-log",
                            "timebase": 1,
                            "dims": ["time"],
                            "variables": [
                                ("Lat", "latitude", ""),
                                ("Lon", "longitude", ""),
                                ("Height", "altitude", ""),
                                ("GroundSpeed", "ground_speed", ""),
                                ("Direction", "heading", ""),
                                ("BaroAlt", "pressure_altitude", ""),
                                ("TAS", "true_air_speed", ""),
                                ("Roll", "roll", ""),
                                ("Pitch", "pitch", ""),
                                ("Yaw", "yaw", ""),
                                ("MagHdg", "heading_mag", ""),
                                ("AGL", "height_agl", ""),
                            ],
                        }
                    },
                }
            },
        },
    }

    clear_config = {
        "kind": "ClearSkyFlight",
        "metadata": {
            "project": "VP2022",
            "platform": "FVR-55",
            "flight_id": "Flight_05",
            "payload_id": "ClearSky",
        },
        "events": {
            "preflight": {},
            "flight": {
                "start_time": "2022-02-16T20:00:00Z",
                "end_time": "2022-02-16T22:15:00Z",
                "datasystems": ["ClearSky", "AutoPilot", "GroundStation"],
            },
            "postflight": {},
        },
        "datasystems": {
            "ClearSky": {
                "payload": {
                    "base_path": "./data/payload",
                    "instruments": {
                        "payload": {
                            "format": "clear-payload-dat",
                            "timebase": 1,
                            "dims": ["time"],
                            "variables": [
                                ("CONCN", "CONCN", "CN concentration"),
                                ("ABSRA", "ABSRA", "PSAP red sample"),
                                ("ABSRB", "ABSRB", "PSAP red reference"),
                                ("ABSGA", "ABSGA", "PSAP green sample"),
                                ("ABSGB", "ABSGB", "PSAP green reference"),
                                ("ABSBA", "ABSBA", "PSAP blue sample"),
                                ("ABSBB", "ABSBB", "PSAP blue reference"),
                                ("CHMPS", "CHMPS", "Chem filter number"),
                                ("AT", "AT", "Air temp from slow probe"),
                                ("RH", "RH", "RH from slow probe"),
                                ("SMPFL", "SMPFL", "MCPC sample flow"),
                                ("SMPFP", "SMPFP", "MCPC sample pump power settings"),
                                ("SATFL", "SATFL", "MCPC Saturator flow"),
                                (
                                    "SATFP",
                                    "SATFP",
                                    "MCPC saturator pump power settings",
                                ),
                                ("ABSFL", "ABSFL", "PSAP flow"),
                                (
                                    "ABSFP",
                                    "ABSFP",
                                    "MCPC saturator pump power settings",
                                ),
                                ("CHMFL", "CHMFL", "Chem filter flow"),
                                ("CHMFP", "CHMFP", "Chem filter pump power settings"),
                                ("OPCFL", "POPS_FL", "POPS flow"),
                                ("OPCFP", "OPCFP", "POPS pump power settings"),
                                ("OPTCT", "OPTCT", "MCPC optics block temp"),
                                ("OPTCP", "OPTCP", "MCPC optics block power"),
                                ("CONDT", "CONDT", "MCPC condenser temperature"),
                                ("CONDP", "CONDP", "MCPC condenser power settings"),
                                ("SATTT", "SATTT", "MCPC Saturator top temp"),
                                ("SATTP", "SATTP", "MCPC Saturator power setting"),
                                ("SATBT", "SATBT", "MCPC Saturator bottom temp"),
                                (
                                    "SATBP",
                                    "SATBP",
                                    "MCPC Saurator bottom power setting",
                                ),
                                ("INLTT", "INLTT", "Temp at the LEF manifold"),
                                (
                                    "FILLC",
                                    "FILLC",
                                    "BuOH indicator, starts incrementing when BuOH is below full",
                                ),
                                ("CABNT", "CABNT", "Temp on the MCPC side of payload"),
                                (
                                    "PRESS",
                                    "PRESS",
                                    "Absolute pressure at in the flow manifold",
                                ),
                                ("PSAP-T", "PSAP_T", "PSAP Temp"),
                                ("PSAP-RH", "PSAP-RH", "PSAP RH"),
                                ("POPS-T", "POPS_T", "POPS Temp"),
                                ("POPS-RH", "POPS_RH", "POPS RH"),
                                ("FastT", "FastT", "Fast temp sensor"),
                                ("FastRH", "FastRH", "Fast RH sensor"),
                            ],
                        }
                    },
                },
                "pops": {
                    "base_path": "/home/derek/Data/UAS/UAS_Flights/ClearSky/VP2022/Flight_05/data/pops",
                    "instruments": {
                        "pops": {
                            "format": "pops-bin",
                            "timebase": 30,
                            "resample": False,
                            "dims": ["time", "pops_bins"],
                            "bin_count": 26,
                            "calibration_file": "/home/derek/Data/UAS/pops_calibration/Cal_curve.csv",
                            "process_options": {"type": "pops"},
                            "variables": [
                                (
                                    "bin_counts",
                                    "pops_bin_counts",
                                    "clear_pops_bin_counts",
                                ),
                                (
                                    "diameter_um",
                                    "pops_dp_um_2d",
                                    "cloudy_msems_diameter_um",
                                ),
                            ],
                        }
                    },
                },
            },
            "GroundStation": {
                "uasground": {
                    "base_path": "/home/derek/Data/envDataSystem/from_cloudbase",
                    "instruments": {
                        "amcpc_cn": {
                            "format": "envdsys",
                            "timebase": 1,
                            "dims": ["time"],
                            "variables": [
                                ("concentration", "ground_cn", "amcpc_cn_concentration")
                            ],
                        }
                    },
                }
            },
            "AutoPilot": {
                "navigation": {
                    "base_path": "./data/piccolo",
                    "instruments": {
                        "piccolo": {
                            "format": "piccolo-log",
                            "timebase": 1,
                            "dims": ["time"],
                            "variables": [
                                ("Lat", "latitude", ""),
                                ("Lon", "longitude", ""),
                                ("Height", "altitude", ""),
                                ("GroundSpeed", "ground_speed", ""),
                                ("Direction", "heading", ""),
                                ("BaroAlt", "pressure_altitude", ""),
                                ("TAS", "true_air_speed", ""),
                                ("Roll", "roll", ""),
                                ("Pitch", "pitch", ""),
                                ("Yaw", "yaw", ""),
                                ("MagHdg", "heading_mag", ""),
                                ("AGL", "height_agl", ""),
                            ],
                        }
                    },
                }
            },
        },
    }
    # de = DataEvent(preflight_config)
    de = DataEvent("flight", config=clear_config)
    print(de._data)

    # save_for_invert_msems("cloudy_msems", event_id="preflight", datasystem="CloudySky", controller="uas_cloudy", config=config2)

