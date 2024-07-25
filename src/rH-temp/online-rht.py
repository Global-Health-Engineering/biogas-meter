import os
import git
import time
import datetime
import AHT20

"""
source of AHT20 library: https://github.com/Chouffy/python_sensor_aht20?tab=readme-ov-file
"""

def get_git_root(path):
    git_repo = git.Repo(path, search_parent_directories=True)
    return git_repo.working_dir


def main():
    aht20 = AHT20.AHT20(BusNum=1)

    now = datetime.datetime.now()
    now_string = now.strftime("%Y-%m-%d %H:%M:%S")

    with open(os.path.join(get_git_root(os.getcwd()), "data", "raw_data", "rH-temp", f"{now_string}.csv"), "w") as fh:
        fh.write("Time,Relative humidity [%RH],Temperature [degC]\n")
        while True:
            # get time
            now = datetime.datetime.now()
            nowstr = now.strftime("%Y-%m-%d %H:%M:%S")
            # get humidity
            rH = aht20.get_humidity()
            # get temperature
            temp = aht20.get_temperature()

            # Fill a string with date, humidity and temperature
            print(f"{nowstr}: {rH:10.2f} %RH, {temp:10.2f} °C")
            
            # crc8 check (cyclic redundancy check)
            # rH_crc8 = aht20.get_humidity_crc8()
            # T_crc8 = aht20.get_temperature_crc8()
            
            # Append in a file
            fh.write(f"{nowstr},{rH},{temp}\n")

            # Wait 2 sec to avoid sensor drift
            time.sleep(2)


if __name__ == "__main__":
    main()
