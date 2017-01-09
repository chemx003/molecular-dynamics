void bCond(double x[], double y[], double z[], double l, int n);

double dRand(double dMin, double dMax);

double errorFile(double x[], double y[], double z[], double ex[], double ey[], double ez[], double gx1, double gy1, double gz1, double gx2, double gy2, double gz2, double dpot_dci, double dpot_dcij, double deps_dci, double deps_dcij, double eps1, double e1Mag, double e2Mag, double rij_mag, double rho, double drhoterm, double epsilon, int i, int j);

void halfstep(double x[], double y[], double z[], double vx[], double vy[], double vz[], double fx[], double fy[], double fz[], double ex[], double ey[], double ez[], double ux[], double uy[], double uz[],double gx[], double gy[], double gz[], double mass, double dt, int n, double I);

void gb(double x[], double y[], double z[], double fx[], double fy[], double fz[], double ex[], double ey[], double ez[], double gx[], double gy[], double gz[], double& V, double l,double& P, double kB, double T, int n, double sigE, int loop);

void gbTest(double x[], double y[], double z[], double fx[], double fy[], double fz[], double ex[], double ey[], double ez[], double gx[], double gy[], double gz[], double& V, double l, double& P, double kB, double T, int n, double sigE, int loop);

void init(double x[], double y[], double z[], double vx[], double vy[], double vz[], double ux[], double uy[], double uz[], double ex[], double ey[], double ez[], double m[], double mass, double I, double l, double dt, double temp, double kB, int n);

void initVerlet(double x[], double y[], double z[], double vx[], double vy[], double vz[], double ux[], double uy[], double uz[], double ex[], double ey[], double ez[], double m[], double mass, double I, double l, double dt, double temp, double kB, int n);

void orientationInfo(double ex[], double ey[], double ez[], int n);

double pairCor(double x[], double y[], double z[], int n, double l);

void leapfrog(double x[], double y[], double z[], double vx[], double vy[], double vz[], double fx[], double fy[], double fz[], double ex[], double ey[], double ez[], double ux[], double uy[], double uz[], double gx[], double gy[], double gz[], double mass,double I, double& K, double dt, int n, double& sumvx, double& sumvy, double& sumvz, double l, int loop);

void orientMag(double ex[], double ey[], double ez[], int n);

void verlet(double x[], double y[], double z[], double vx[], double vy[], double vz[], double fx[], double fy[], double fz[], double ex[], double ey[], double ez[], double ux[], double uy[], double uz[], double gx[], double gy[], double gz[], double mass, double I, double dT, double& K,int n);

void writeEnergy(double V, double K, double E, int time);

void writeVectors(double x[], double y[], double z[], double ex[], double ey[], double ez[], double i, int n);
