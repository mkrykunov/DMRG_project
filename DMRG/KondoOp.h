
static double Lna[4][4] = {{ 0, 1, 0, 0 },
						   { 0, 0, 0, 0 },
						   { 0, 0, 0, 1 },
						   { 0, 0, 0, 0 }};

static double Lnb[4][4] = {{ 0, 0, 1, 0 },
						   { 0, 0, 0,-1 },
						   { 0, 0, 0, 0 },
						   { 0, 0, 0, 0 }};

static int nLna[4] = { 0, 1, 0, 1 };

static int nLnb[4] = { 0, 0, 1, 1 };


static double SLn_up[4][4] = {{ 0, 0, 0, 0 },
							  { 0, 0, 1, 0 },
							  { 0, 0, 0, 0 },
							  { 0, 0, 0, 0 }};

static double SLn_z[4][4] = {{ 0,   0,   0, 0 },
							 { 0, 0.5,   0, 0 },
							 { 0,   0,-0.5, 0 },
							 { 0,   0,   0, 0 }};


static double Ana[4][4] = {{ 0, 1, 0, 0 },
						   { 0, 0, 0, 0 },
	    				   { 0, 0, 0,-1 },
						   { 0, 0, 0, 0 }};

static double Anb[4][4] = {{ 0, 0, 1, 0 },
						   { 0, 0, 0, 1 },
						   { 0, 0, 0, 0 },
						   { 0, 0, 0, 0 }};

static int nAna[4] = { 0, 1, 0, 1 };

static int nAnb[4] = { 0, 0, 1, 1 };


static double dna[4][4] = {{ 0, 1, 0, 0 },
						   { 0, 0, 0, 0 },
						   { 0, 0, 0, 1 },
						   { 0, 0, 0, 0 }};

static double dnb[4][4] = {{ 0, 0, 1, 0 },
						   { 0, 0, 0,-1 },
						   { 0, 0, 0, 0 },
						   { 0, 0, 0, 0 }};

static int ndna[2] = { 1, 0 };

static int ndnb[2] = { 0, 1 };

static double Sdn_up[2][2] = {{ 0, 1 },
							  { 0, 0 }};

static double Sdn_z[2][2] = {{ 0.5,   0 },
							 {   0,-0.5 }};
