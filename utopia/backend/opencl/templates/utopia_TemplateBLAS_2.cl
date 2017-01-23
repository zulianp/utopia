//Generated from template file: utopia_TemplateBLAS_2.cl
#include "utopia_Operations.cl"

////////////////////////////////////////////////////////////////////////////////////////////////////

__kernel void <@KERNEL_NAME@>( 
	const SizeType n_entries,
	<@[list]ARG_IN@>)
{
	const SizeType rows    = <@ROWS@>;
	const SizeType columns = <@COLUMNS@>;

	for(SizeType index0 = get_global_id(0); index0 < rows; index0 += get_global_size(0)) {
		const SizeType offset = index0 * columns;
		for(SizeType index1 = get_global_id(1); index1 < columns; index1 += get_global_size(1)) {
			const SizeType index = offset + index1;

			//see if it is necessary to enable it through macros
			const SizeType index_transposed = index1 * rows + index0;
			<@OP@>;
		}
	}	
}

///////////////////////////////////////////////////////////////////////////////////////////////////

