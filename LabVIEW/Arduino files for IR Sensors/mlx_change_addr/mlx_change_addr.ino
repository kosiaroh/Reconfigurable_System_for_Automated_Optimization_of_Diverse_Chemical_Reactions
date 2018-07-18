#include <Wire.h>
#include <Adafruit_MLX90614.h>
#include <i2cmaster.h>

// Pins: Standard: SDA:A4  SCL:A5


byte MLXAddr  = 0x04; //Current Address
//Universal address: byte MLXAddr = 0;
Adafruit_MLX90614 mlx = Adafruit_MLX90614();
const int button1Pin = 2;  // pushbutton 1 pin

void setup() {
  pinMode(button1Pin, INPUT);
  Serial.begin(9600);

  Serial.println("Setup...");
 
  i2c_init();                              //Initialise the i2c bus
  //Using exteneral pullups //PORTC = (1 << PORTC4) | (1 << PORTC5);   //enable pullups
 
  delay(5000);                    // Wait to allow serial connection
  ReadAddr(0);                    // Read current address bytes
  ChangeAddr(MLXAddr, 0x00);         // Change address to new value
 //ChangeAddr(0x5A, 0xBE);       // Change address to default value
  ReadAddr(0);                    // Read address bytes
  delay(10000);                    // Cycle power to MLX during this pause
  ReadTemp(0);                    // Read temperature using default address
  ReadTemp(MLXAddr);              // Read temperature using new address
}

void loop(){
   delay(1000); // wait a second
}

word ChangeAddr(byte NewAddr1, byte NewAddr2) {

 Serial.println("> Change address");

 i2c_start_wait(0 + I2C_WRITE);    //send start condition and write bit
 i2c_write(0x2E);                  //send command for device to return address
 i2c_write(0x00);                  // send low byte zero to erase
 i2c_write(0x00);                  //send high byte zero to erase
 if (i2c_write(0x6F) == 0) {
   i2c_stop();                     //Release bus, end transaction
   Serial.println("  Data erased.");
 }
 else {
   i2c_stop();                     //Release bus, end transaction
   Serial.println("  Failed to erase data");
   return -1;
 }

 Serial.print("  Writing data: ");
 Serial.print(NewAddr1, HEX);
 Serial.print(", ");
 Serial.println(NewAddr2, HEX);

 for (int a = 0; a != 256; a++) {
   i2c_start_wait(0 + I2C_WRITE);  //send start condition and write bit
   i2c_write(0x2E);                //send command for device to return address
   i2c_write(NewAddr1);            // send low byte zero to erase
   i2c_write(NewAddr2);            //send high byte zero to erase
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

void ReadAddr(byte Address) {

 Serial.println("> Read address");

 Serial.print("  MLX address: ");
 Serial.print(Address, HEX);
 Serial.print(", Data: ");

 i2c_start_wait(Address + I2C_WRITE);  //send start condition and write bit
 i2c_write(0x2E);                  //send command for device to return address
 i2c_rep_start(Address + I2C_READ);
 
 Serial.print(i2c_readAck(), HEX); //Read 1 byte and then send ack
 Serial.print(", ");
 Serial.print(i2c_readAck(), HEX); //Read 1 byte and then send ack
 Serial.print(", ");
 Serial.println(i2c_readNak(), HEX);
 i2c_stop();
}

byte calcCrc8(byte X, byte crc) {
 X = X ^ crc;
 for (byte i=0; i<8; i++) {
   if (X > 0x7F) {
     X <<= 1;
     X = X ^ 0x07;
   }
   else {
     X <<= 1;
   }
 }
 return X;
}

byte calcPec(byte slaveAddr, byte memAddr, byte loByte, byte hiByte) {
 byte crc = 0x00;
 byte slaveAddr2 = slaveAddr << 1;
 crc = calcCrc8(slaveAddr2, crc);
 crc = calcCrc8(memAddr, crc);
 crc = calcCrc8(slaveAddr2+1, crc);
 crc = calcCrc8(loByte, crc);
 crc = calcCrc8(hiByte, crc);
 crc = calcCrc8(0x19, crc);
 return crc;
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


/*  Serial.println("Adafruit MLX90614 test");  

  mlx.begin();  
}

void loop() {
  int button1State;

  button1State = digitalRead(button1Pin);

  if (button1State == LOW){
    Serial.print("Ambient = "); Serial.print(mlx.readAmbientTempC()); 
    Serial.print("*C\tObject = "); Serial.print(mlx.readObjectTempC()); Serial.println("*C");
    Serial.print("Ambient = "); Serial.print(mlx.readAmbientTempF()); 
    Serial.print("*F\tObject = "); Serial.print(mlx.readObjectTempF()); Serial.println("*F");
  
    Serial.println();
    delay(500);
  }else{
    
  }
}*/
