import os
import shutil
import sys
import json

if __name__ == "__main__":

    base_path = os.path.abspath("")
    if not os.path.isdir("templates"):
        base_path = os.path.dirname(os.path.abspath(""))
        os.chdir("..")
    # print(os.path.isdir("templates"))
    # print(os.getcwd())
    if base_path not in sys.path:
        sys.path.insert(0, base_path)
    # base_path = sys.path.insert(0, os.path.dirname(os.path.abspath("")))

    # print(sys.path)
    from config.settings import data_paths

    # sys.argv.append("create_flight")
    # sys.argv.append("--project=TESTTESTTEST")
    # print(sys.argv)

    if len(sys.argv) < 2:
        print("nothing to do")
        exit
    action = sys.argv[1]
    # print(f"action = {action}")

    # print(data_paths)

    if action == "create_flight" or action == "update_flight":
        if len(sys.argv) < 3:
            print("don't know what to create")
            exit

        args = sys.argv[2:]

        payload = "CloudySky"
        project = "Default"
        platform = "FVR-55"
        flight_id = "Flight_XX"

        for arg in args:
            option, value = arg.split("=")
            if option == "--payload":
                payload = value
            elif option == "--project":
                project = value
            elif option == "--platform":
                platform = value
            elif option == "--flight_id":
                flight_id = value

        tmpl_name = "CloudySky_data_tmpl.ipynb"
        data_folders = ["piccolo", "msems_inversion", "igor"]
        if payload == "ClearSky":
            tmpl_name = "ClearSky_data_tmpl.ipynb"
            data_folders = ["payload", "pops", "sasp", "piccolo", "igor"]

        nb_name = tmpl_name.replace("tmpl", f"{project}_{flight_id}")
        # print(nb_name)

        # create directories
        try:
            for fld in data_folders:
                pth = os.path.join(payload, project, flight_id, "data", fld)
                os.makedirs(pth)
        except FileExistsError:
            pass

        try:
            pth = os.path.join(payload, project, flight_id, "config")
            os.makedirs(pth)
        except FileExistsError:
            pass

        success = "updated"
        if action == "create_flight":
            # copy template to flight folder
            shutil.copy(
                os.path.join("templates", tmpl_name),
                os.path.join(payload, project, flight_id, nb_name),
            )
            success = "created"

        # create flight config
        config = {
            "payload_id": payload,
            "project": project,
            "platform": platform,
            "flight_id": flight_id,
            "data_paths": data_paths
        }

        with open(os.path.join(payload, project, flight_id, "config", "settings.json"), "w") as f:
            json.dump(config, f)

        print(f"Flight {success}: {payload}/{project}/{flight_id}")
        # os.path.join("templates", "CloudySky_data_tmpl.ipynb")
