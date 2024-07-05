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
	@echo 'Invoking: Intel C++ Compiler'
	icpx -g -O3 -inline-level=2 -I../include -I../../../src -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


