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
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel C++ Compiler'
	icpc -g -O3 -mkl=parallel -I../../../src -no-multibyte-chars -fargument-noalias -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


