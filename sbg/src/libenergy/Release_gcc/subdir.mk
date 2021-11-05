################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../icosa.cpp \
../korpe.cpp \
../pd2.cpp \
../rama.cpp \
../tobi.cpp 

OBJS += \
./icosa.o \
./korpe.o \
./pd2.o \
./rama.o \
./tobi.o 

CPP_DEPS += \
./icosa.d \
./korpe.d \
./pd2.d \
./rama.d \
./tobi.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I../include -I../../../src -O3 -Wall -c -fmessage-length=0 -fpermissive -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


