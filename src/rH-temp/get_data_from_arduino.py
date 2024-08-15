import os
import serial
from datetime import datetime


def readserial(comport, baudrate):
    ser = serial.Serial(comport, baudrate, timeout=0.1)
    now = datetime.now()
    fname = os.path.join(os.getcwd(), "rH_data", f"{now.strftime('%Y-%m-%d_%H%M%S')}.log")
    with open(fname, "a") as fh:
        fh.write("Datetime, Temp [deg C], Humidity [rH]\n")
        while True:
            now = datetime.now()
            nowstr = now.strftime('%Y-%m-%d %H:%M:%S')
            data = ser.readline().decode().strip()
            if data:
                line = f"{nowstr}, {data}"
                print(line)
                fh.write(f"{line}\n")


def main():
    readserial(comport='COM5', baudrate=9600)


if __name__ == '__main__':
    main()
