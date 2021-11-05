################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Io.cpp \
../Memstuff.cpp \
../Surface.cpp 

OBJS += \
./Io.o \
./Memstuff.o \
./Surface.o 

CPP_DEPS += \
./Io.d \
./Memstuff.d \
./Surface.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel C++ Compiler'
	icpc -O3 -inline-level=2 -I../include -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


