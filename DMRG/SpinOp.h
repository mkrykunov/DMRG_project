#if __SPIN == 2

double Sn_up[2][2] = {{ 0, 1 },
					  { 0, 0 }};

double Sn_z[2][2] = {{ 0.5,   0 },
					 {   0,-0.5 }};

#elif __SPIN == 3

double Sn_up[3][3] = {{ 0, 1, 0 },
					  { 0, 0, 1 },
					  { 0, 0, 0 }};

double Sn_z[3][3] = {{ 1, 0, 0 },
					 { 0, 0, 0 },
					 { 0, 0,-1 }};

#endif
