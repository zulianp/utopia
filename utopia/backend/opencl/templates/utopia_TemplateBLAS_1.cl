//Generated from template file: utopia_TemplateBLAS_1.cl

#include "utopia_Operations.cl"

<@[list]FUNCTIONS@>

//auto generated kernel
__kernel void <@KERNEL_NAME@>(
	const SizeType n_entries, 
	//auto generated variables here below:
	<@[list]ARG_IN@>)
{
	for(SizeType index =  get_global_id(0); index < n_entries; index += get_global_size(0)) {
		//beginning of the composite expression
		<@OP@>;
		//end of the composite expression
	}

#ifdef VERBOSE
	if(get_global_id(0) == 0) {
		printf("[KERNEL_MESSAGE]: <@KERNEL_NAME@> terminated\n");
	}
#endif //VERBOSE	
}
