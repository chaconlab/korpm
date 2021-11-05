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
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I../include -O2 -U_FORTIFY_SOURCE -Wall -c -fmessage-length=0 -fno-omit-frame-pointer -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


