//Generated from template file: utopia_TemplateReduce.cl

#include "utopia_Operations.cl"

#pragma OPENCL EXTENSION cl_khr_fp64 : enable

#define STRIDE 1

__kernel void <@KERNEL_NAME@>(
	const SizeType n_entries, 
	__global Scalar * result, 
	__local Scalar  * local_work,
	//auto generated variables here below:
	<@[list]ARG_IN@>)
{
	Scalar private_reduction[STRIDE];

	for(SizeType d = 0; d < STRIDE; ++d) {
		private_reduction[0] = 0; //TODO !@INITIAL_VALUE@>;
	}
	
	for(SizeType index =  get_global_id(0); index < n_entries; index += get_global_size(0)) {
		//the composite expression here below:
		for(SizeType d = 0; d < STRIDE; ++d) {
			private_reduction[0] = <@REDUCE_OP@>(private_reduction[0], <@OP@>);
		}
	}

	const SizeType l = get_local_id(0);                                                                
	const SizeType localSize = get_local_size(0);                                                      
	SizeType offset = localSize >> 1;                                                                  

	generic_copy(STRIDE, private_reduction, &local_work[l*STRIDE]);                                   
	for(SizeType k = 0; offset > 0; ++k) {                                                             

		SizeType index = (k + offset) * STRIDE; 

		barrier(CLK_LOCAL_MEM_FENCE);                                                                

		if(index < localSize) { 
			generic_copy(STRIDE, (&local_work[index]), private_reduction);  
		}

		barrier(CLK_LOCAL_MEM_FENCE);                                                                

		for(SizeType i = 0; i < STRIDE; ++i) {
			local_work[index + i] = <@REDUCE_OP@>(private_reduction[i], local_work[index+i]);
		}

		offset >>= 1;                                                                                
	}                                                                                                  

	if(l == 0) {                                                                                       
		const SizeType global_offset = get_group_id(0) * STRIDE;										            
		generic_copy(STRIDE, &local_work[0], &result[global_offset]);			         				   
	}	                              
}

