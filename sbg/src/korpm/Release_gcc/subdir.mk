################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../korpm.cpp 

OBJS += \
./korpm.o 

CPP_DEPS += \
./korpm.d 


# Each subdirectory must supply rules for building sources it contributes
korpm.o: ../korpm.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I../../../src -O3 -g -Wall -c -fmessage-length=0 -fpermissive -Wmaybe-uninitialized -MMD -MP -MF"$(@:%.o=%.d)" -MT"korpm.d" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


