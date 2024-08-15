#include <Wire.h>
#include <Adafruit_AHTX0.h>

// declare objects
Adafruit_AHTX0 aht;

void setup() {
  // Initialize serial communication at 9600 baud
  Serial.begin(9600);

  // initialize AHT20 sensor
  if (!aht.begin()) {
      Serial.println("Could not find AHT20");
      Serial.flush();
      while (1);
    }
}

void loop() {
  sensors_event_t humidity, temp;
  aht.getEvent(&humidity, &temp);
  Serial.print(temp.temperature);
  Serial.print(", ");
  Serial.println(humidity.relative_humidity);

  // delay 2 seconds
  delay(2000);
}
  