#include <Wire.h>
#include <Adafruit_MLX90614.h>
#include <i2cmaster.h>

// Pins: Standard: SDA:A4  SCL:A5


byte MLXAddr  = 0x00;
float emissivity = .8;
//Universal address: byte MLXAddr = 0;
Adafruit_MLX90614 mlx = Adafruit_MLX90614();

void setup() {
  Serial.begin(9600);

  Serial.println("Setup...");
 
  i2c_init();                              //Initialise the i2c bus
  //Using exteneral pullups //PORTC = (1 << PORTC4) | (1 << PORTC5);   //enable pullups

  unsigned int emiss = round(65535 * emissivity );
  Serial.print("Emissivity: ");
  Serial.println(emiss);
  Serial.print("Lower Byte: ");
  Serial.println(emiss & 0xFF);
  Serial.print("Higher Byte: ");
  Serial.println(emiss >> 8);
  
  delay(5000);                    // Wait to allow serial connection
  ReadEmissivity(0);                    // Read current Emissivity bytes
  ReadTemp(0);                   // Read temperature using current emissivity
  ChangeEmissivity(emiss & 0xFF, emiss >> 8);         // Change address to new value
  ReadEmissivity(0);                    // Read address bytes
  delay(10000);                    // Cycle power to MLX during this pause
  ReadTemp(0);                    // Read temperature using new emissivity
}

void loop(){
   delay(1000); // wait a second
}

word ChangeEmissivity(byte NewEmis1, byte NewEmis2) {

 Serial.println("> Change emissivity");

 i2c_start_wait(0 + I2C_WRITE);    //send start condition and write bit
 i2c_write(0x24);                  //send command for device to return emissivity
 i2c_write(0x00);                  // send low byte zero to erase
 i2c_write(0x00);                  //send high byte zero to erase
 if (i2c_write(0xE8) == 0) {
   i2c_stop();                     //Release bus, end transaction
   Serial.println("  Data erased.");
 }
 else {
   i2c_stop();                     //Release bus, end transaction
   Serial.println("  Failed to erase data");
   return -1;
 }

 Serial.print("  Writing data: ");
 Serial.print(NewEmis1, HEX);
 Serial.print(", ");
 Serial.println(NewEmis2, HEX);

 for (int a = 0; a != 256; a++) {
   i2c_start_wait(0 + I2C_WRITE);  //send start condition and write bit
   i2c_write(0x24);                //send command for device to return emissivity
   i2c_write(NewEmis1);            // send low byte zero to erase
   i2c_write(NewEmis2);            //send high byte zero to erase
   if (i2c_write(a) == 0) {
     i2c_stop();                   //Release bus, end transaction
     delay(100);                   // then wait 10ms
     Serial.print("Found correct CRC: 0x");
     Serial.println(a, HEX);
     return a;
   }
 }
 i2c_stop();                       //Release bus, end transaction
 Serial.println("Correct CRC not found");
 return -1;
}

void ReadEmissivity(byte Address) {

 Serial.println("> Read Emissivity");

 Serial.print("  MLX address: ");
 Serial.print(Address, HEX);
 Serial.print(", Data: ");

 i2c_start_wait(Address + I2C_WRITE);  //send start condition and write bit
 i2c_write(0x24);                  //send command for device to return emissivity
 i2c_rep_start(Address + I2C_READ);
 
 Serial.print(i2c_readAck(), HEX); //Read 1 byte and then send ack
 Serial.print(", ");
 Serial.print(i2c_readAck(), HEX); //Read 1 byte and then send ack
 Serial.print(", ");
 Serial.println(i2c_readNak(), HEX);
 i2c_stop();
}

float ReadTemp(byte Address) {
 int data_low = 0;
 int data_high = 0;
 int pec = 0;

 Serial.println("> Read temperature");

 Serial.print("  MLX address: ");
 Serial.print(Address, HEX);
 Serial.print(", ");

 i2c_start_wait(Address + I2C_WRITE);
 i2c_write(0x07);                  // Address of temp bytes
 
 // read
 i2c_rep_start(Address + I2C_READ);
 data_low = i2c_readAck();         //Read 1 byte and then send ack
 data_high = i2c_readAck();        //Read 1 byte and then send ack
 pec = i2c_readNak();
 i2c_stop();
 
 //This converts high and low bytes together and processes temperature, MSB is a error bit and is ignored for temps
 float Temperature = 0x0000;       // zero out the data
 
 // This masks off the error bit of the high byte, then moves it left 8 bits and adds the low byte.
 Temperature = (float)(((data_high & 0x007F) << 8) + data_low);
 Temperature = (Temperature * 0.02) - 273.16;
 
 Serial.print(Temperature);
 Serial.println(" C");
 return Temperature;
}
