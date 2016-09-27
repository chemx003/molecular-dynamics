/*I'm going to transcribe the example code from the Allen and Tildesley on 
github.
	
	I've proof read it once, but I just need to create the function specify
the arguments and declare all of the shit*/

double mu=2.0, nu=1.0;
double kappa=3.0, xappa=5.0;

double chi=(pow(kappa,2)-1.0)/(pow(kappa,2)-1.0);
double xhi=(pow(xappa,1.0/mu)-1.0)/(pow(xappa,1/mu)+1.0);

double rc=4.0;

double dx,dy,dz;
double rij_sq,rij_mag;
double dx_hat,dy_hat,dz_hat;
double ci,cj,cij;
double cp, cm;
double cpchi, cmchi, sigma;
double eps1, cpxhi, cmxhi,eps2,epsilon;
double rho, rho6, rho12, rhoterm,drhoterm,pot;
double prefac, dsig_dci, dsig_dcj, dsig_dcij;
double deps_dci, deps_dcj, deps_cij;
double dpot_drij, dpot_dci, dpot_dcj, dpot_dcij;
double fxi, fyi, fzi;
double g1x,g1y,g1z,g2x,g2y,g2z;

for(int i=0; i<n-1; i++){
	for(int j=i+1; j<n; j++){
		dx=x[i]-x[j];
		dy=y[i]-y[j];
		dz=z[i]-z[j];

		rij_sq=dx*dx+dy*dy+dz*dz;
		rij_mag=pow(rij_sq,0.5);

		dx_hat=dx/rij_mag;
		dy_hat=dy/rij_mag;
		dz_hat=dz/rij_mag;

		ci=dx_hat*ex[i]+dy_hat*ey[i]+dz_hat*ez[i];
		cj=dx_hat*ex[j]+dy_hat*ey[j]+dz_hat*ez[j];
		cij=ex[i]*ex[j]+ey[i]*ey[j]+ez[i]*ez[j];

		cp=ci+cj;
		cm=ci-cj;

		cpchi=cp/(1.0+chi*cij);
		cmchi=cm/(1.0-chi*cij);
		sigma=1.0/pow(1.0-0.5*chi*(cp*cpchi+cm*cmchi),0.5);

		eps1=1.0/pow(1.0-(chi*chi*cij*cij),0.5);
		cpxhi=cp/(1.0+xhi*cij);
		cmxhi=cm/(1.0-xhi*cij);
		eps2=1.0-0.5*xhi*(cp*cpxhi+cm*cmxhi);
		epsilon=(pow(eps1,nu))*(pow(eps2,mu));

		rho=rij_mag-sigma+1.0;
		rho6=1.0/pow(rho,6);
		rho12=rho6*rho6;
		rhoterm=4.0*(rho12-rho6);
		drhoterm=-24.0*(2.0*rho12-rho6)/rho;
		pot=epsilon*rhoterm;
		
		prefac=0.5*chi*pow(sigma,3);
		dsig_dci=prefac*(cpchi+cmchi);
		dsig_dcj=prefac*(cpchi-cmchi);
		prefac=prefac*(0.5*chi);
		dsig_dcij=-prefac*(cpchi*cpchi-cmchi*cmchi);

		prefac=-mu*xhi*pow(eps1,nu)*pow(eps2,mu-1);
		deps_dci=prefac*(cpxhi+cmxhi);
		deps_dcj=prefac*(cpxhi-cmxhi);
		prefac=preac*(0.5*xhi);
		deps_dcij=-prefac*(cpxhi*cpxhi-cmxhi*cmxhi);
		deps_dcij=deps_dcij+nu*chi*chi*pow(eps1,nu+2)*pow(eps2,mu)*cij;
		
		dpot_drij=epsilon*drhoterm;
		dpot_dci=rhoterm*deps_dci-epsilon*drhoterm*dsig_dci;
		dpot_dcj=rhoterm*deps_dcj-epsilon*drhoterm*dsig+dcj;
		dpot_dcij=rhoterm*deps_dcij-epsilon*drhoterm*dsig_dcij;

		fxi=-dpot_drij*dx_hat+dpot_dci*(ex[i]-ci*dx_hat)/rij_mag
			-dpot_dcj*(ex[j]-cj*dx_hat)/rij_mag;
		fyi=-dpot_drij*dy_hat+dpot_dci*(ey[i]-ci*dy_hat)/rij_mag
			-dpot_dcj*(ey[j]-cj*dy_hat)/rij_mag;
		fzi=-dpot_drij*dz_hat+dpot_dci*(ez[i]-ci*dz_hat)/rij_mag
			-dpot_dcj*(ez[j]-cj*dz_hat)/rij_mag;

		g1x=dpot_dci*dx_hat+dpot_dcij*ex[j];
		g1y=dpot_dci*dy_hat+dpot_dcij*ey[j];
		g1z=dpot_dci*dz_hat+dpot_dcij*ez[j];

		g2x=dpot_dci*dx_hat+dpot_dcij*ex[i];
		g2y=dpot_dci*dy_hat+dpot_dcij*ey[i];
		g2z=dpot_dci*dz_hat+dpot_dcij*ez[i];

		fx[i]=fx[i]+fxi;
		fy[i]=fy[i]+fyi;
		fz[i]=fz[i]+fzi;

		gx[i]=gx[i]-g1x;
		gy[i]=gy[i]-g1y;
		gz[i]=gz[i]-g1z;

		gx[j]=gx[j]-g2x;
		gy[j]=gy[j]-g2y;
		gz[j]=gz[j]-g2z;
	}
}
