# Temperature and humidity sensor AHT20

CircuitPython driver for the Adafruit AHT10 or AHT20 Humidity and Temperature Sensor is available here:

```
pip3 install adafruit-circuitpython-ahtx0
```

## Usage example

```
import time
import board
import adafruit_ahtx0

# Create sensor object, communicating over the board's default I2C bus
i2c = board.I2C()  # uses board.SCL and board.SDA
sensor = adafruit_ahtx0.AHTx0(i2c)

while True:
    print("\nTemperature: %0.1f C" % sensor.temperature)
    print("Humidity: %0.1f %%" % sensor.relative_humidity)
    time.sleep(2)
```