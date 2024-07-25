import os
import git
import time
import board
import datetime
import adafruit_ahtx0


def get_git_root(path):
    git_repo = git.Repo(path, search_parent_directories=True)
    return git_repo.working_dir


def main():
    # Create sensor object, communicating over the board's default I2C bus
    i2c = board.I2C() # uses board.SCL and board.SDA
    aht20 = adafruit_ahtx0.AHTx0(i2c)

    now = datetime.datetime.now()
    nowStr = now.strftime("%Y-%m-%d %H:%M:%S")

    with open(os.path.join(get_git_root(os.getcwd()),
                           "data",
                           "raw_data",
                           "rH-temp",
                           f"{nowStr}.csv"),
              "w") as fh:
        fh.write("Time,Relative humidity [%RH],Temperature [degC]\n")
        while True:
            # get time
            now = datetime.datetime.now()
            nowStr = now.strftime("%Y-%m-%d %H:%M:%S")

            # Fill a string with date, humidity and temperature
            print(f"{nowStr}: {sensor.relative_humidity:10.2f} %RH, {sensor.temperature:10.2f} °C")
            
            # Append to the end of the file
            fh.write(f"{nowStr},{sensor.relative_humidity},{sensor.temperature}\n")

            # Wait 2 sec to avoid sensor drift
            time.sleep(2)


if __name__ == "__main__":
    main()
