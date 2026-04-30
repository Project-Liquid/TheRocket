# GSE Code Full Integration

## TO USE:
1. Download the entire gse_control folder (or clone into repository)
2. Open only gse_control.ino in the Arduino IDE. The other files will open automatically.
3. For field testing, run gui.py on a separate computer connected to a radio.

Do not add any other files to this folder. To add a feature, create a new .cpp and .h file to be implemented in the same main .ino file. 

## LIBRARIES TO DOWNLOAD:
Arduino:
- [HX711 Arduino Library](https://github.com/bogde/HX711) by Bogdan Necula
- [Encoder](https://www.pjrc.com/teensy/td_libs_Encoder.html) by Paul Stoffregen
- [Adafruit_MCP9600](https://github.com/adafruit/Adafruit_MCP9600) by Adafruit
- [SparkFun Qwiic Scale NAU7802 Library](https://github.com/sparkfun/SparkFun_Qwiic_Scale_NAU7802_Arduino_Library) by sparkfun
- [StandardCplusplus](https://github.com/maniacbug/StandardCplusplus/blob/master/README.md) by maniacbug. NOTE: this library needs to be manually imported. Download the StandardCplusplus.zip folder in this repository (TheRocket/GSE) and import it with Sketch>Include Library>Add .ZIP Library in the Arduino IDE.  
- [ADS1118]() by Alvaro Salazar

Python:
- [pySerial](https://github.com/pyserial/pyserial)
- [PyQt5](https://doc.qt.io/archives/qtforpython-5/)
- [PyQtGraph](https://github.com/pyqtgraph/pyqtgraph)
