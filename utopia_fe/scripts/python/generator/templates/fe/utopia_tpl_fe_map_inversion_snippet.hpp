{{
	static constexpr int max_iter = 100;
	bool converged = false;
	GeoT residual_norm;
	for(int i = 0; i < max_iter; ++i) {{
		{update}

		// FIXME not SIMD friendly
		if(residual_norm < tol) {{
			converged = true;
			break;
		}}
	}}

	assert(converged);
}}