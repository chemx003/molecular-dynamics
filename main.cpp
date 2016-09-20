/*
 * File:   main.cpp,
 * Author: Bryan
 *
 */

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <time.h>
#include <sstream>

using namespace std;

/*Implementation of periodic boundary conditions. Updates positions
  after particle positions are incremented each timestep*/
void bCond(double x[], double y[], double z[], double l, int n){

       for(int i=0; i<n; i++){
        if(x[i]<0){
            if(fmod(abs(x[i]),l)<l/2){
                x[i]=x[i]-l*round(x[i]/l)+l;
            }
            else{
                x[i]=x[i]-l*round(x[i]/l);
            }
        }

        else if(x[i]>l){
            if(fmod(abs(x[i]),l)<l/2){
                x[i]=x[i]-l*round(x[i]/l);
            }
            else{
                x[i]=x[i]-l*round(x[i]/l)+l;
            }
        }
        if(y[i]<0){

           if(fmod(abs(y[i]),l)<l/2){
                y[i]=y[i]-l*round(y[i]/l)+l;
            }
            else{
                y[i]=y[i]-l*round(y[i]/l);
            }
        }
        else if(y[i]>l){
            if(fmod(abs(y[i]),l)<l/2){
                y[i]=y[i]-l*round(y[i]/l);
            }
            else{
                y[i]=y[i]-l*round(y[i]/l)+l;
            }
        }

        if(z[i]<0){
           if(fmod(abs(z[i]),l)<l/2){
                z[i]=z[i]-l*round(z[i]/l)+l;
            }
            else{
                z[i]=z[i]-l*round(z[i]/l);
            }
        }
        else if(z[i]>l){
            if(fmod(abs(z[i]),l)<l/2){
                z[i]=z[i]-l*round(z[i]/l);
            }
            else{
                z[i]=z[i]-l*round(z[i]/l)+l;
            }
        }
    }
}//from ubbcluj

/*Generates a random double between dMin and dMax*/
double dRand(double dMin, double dMax){

    double d = (double)rand()/RAND_MAX;
    return dMin + d*(dMax-dMin);
}

/*Increments the velocity and angular velocity half a timestep and the 
position and orientation a full timestep*/
void halfstep(double x[], double y[], double z[], double vx[], double vy[],
        double vz[], double fx[], double fy[], double fz[], double ex[],
		double ey[], double ez[], double ux[], double uy[], double uz[],
		double gx[], double gy[], double gz[], double mass, double dt, 
		int n, double I){
	
	double lm;
	double dot;

    for(int i=0; i<n; i++){
		
		//Advance velocities 1/2 timestep
        vx[i]=vx[i]+0.5*dt*fx[i]/mass;
        vy[i]=vy[i]+0.5*dt*fy[i]/mass;
        vz[i]=vz[i]+0.5*dt*fz[i]/mass;
    }

    for(int i=0; i<n; i++){

		//Advance positions 1 timestep
        x[i]=x[i]+dt*vx[i];
        y[i]=y[i]+dt*vy[i];
        z[i]=z[i]+dt*vz[i];
    }
	
	for(int i=0; i<n; i++){

		//Calculate lagrange multiplier to constrain length
		lm=-2.0*(ux[i]*ex[i]+uy[i]*ey[i]+uz[i]*ez[i]);
		
		//dot product
		dot=gx[i]*ex[i]+gy[i]*ey[i]+gz[i]*ez[i];
			
		//Take perpendicular component of the gorque
		gx[i]=gx[i]-(dot)*ex[i];
		gy[i]=gy[i]-(dot)*ey[i];
		gz[i]=gz[i]-(dot)*ez[i];

		/*Advance angular velocities 1/2 timestep.*/
        ux[i]=ux[i]+0.5*dt*(gx[i]/I)+0.5*lm*ex[i];
        uy[i]=uy[i]+0.5*dt*(gy[i]/I)+0.5*lm*ey[i];
        uz[i]=uz[i]+0.5*dt*(gz[i]/I)+0.5*lm*ez[i];	
	}

	for(int i=0; i<n; i++){

		//Advance orientations 1 timestep
		ex[i]=ex[i]+dt*ux[i];
        ey[i]=ey[i]+dt*uy[i];
        ez[i]=ez[i]+dt*uz[i];

		/*double mag=pow(ex[i]*ex[i]+ey[i]*ey[i]+ez[i]*ez[i], 0.5);

        ex[i]=ex[i]/mag;
        ey[i]=ey[i]/mag;
        ez[i]=ez[i]/mag;*/
	
	}
}

/*Calculation of the intermolecular forces and torques*/
void gb(double x[], double y[], double z[], double fx[],
        double fy[], double fz[], double ex[], double ey[], double ez[],
        double gx[], double gy[], double gz[], double& V, double l,
        double& P, double kB, double T, int n, double sigE, int loop){

	//Simulation parameters
    double mu=2.0, nu=1.0; 
    double sigmaE=sigE, sigmaS=1.0, epsilonE=5.0, epsilonS=1.0;
    double kappa=sigmaE/sigmaS, kappaPrime=epsilonE/epsilonS;
    double chi=(pow(kappa,2.0)-1.0)/(pow(kappa,2.0)+1.0);
    double chiPrime=(pow(kappaPrime,1.0/mu)-1.0)/(pow(kappaPrime,1.0/mu)+1.0);
    double rc=3.80*sigmaS, rc2=rc*rc; //cuttoff

	//Calculated quantities
	double dx, dy, dz;
    double dot1, dot2, dot12, dot122, dotSum, dotSum2, dotDif, dotDif2;
    double g, gPrime, gHalf, dgx, dgy, dgz, dgxPrime, dgyPrime, dgzPrime;
    double R, R_1, R_2, R_6, R_7, R_12, R_13, distF;
    double ePrime, ePn, ePn1, gPm, gPm1;
    double fxi, fyi, fzi, gx1, gy1, gz1, gx2, gy2, gz2;
    double dotSByChi, dotDByChi, dotSByChip, dotDByChip;
	double dotSByChi2, dotDByChi2, dotSByChip2, dotDByChip2;
    double dex1, dey1, dez1, dex2, dey2, dez2;
    double depx1, depy1, depz1, depx2, depy2, depz2;
    double drx1,dry1,drz1,drx2,dry2,drz2;
    double dgx1,dgy1,dgz1,dgx2,dgy2,dgz2;
    double dgpx1,dgpy1,dgpz1,dgpx2,dgpy2,dgpz2;
	double r, r2, e1Mag, e2Mag;

	//Resetting quantities
    for(int i=0; i<n; i++){
        fx[i]=0; gx[i]=0;
        fy[i]=0; gy[i]=0;
        fz[i]=0; gz[i]=0;
    }
	V=0;
    P=0;

	//Main loop
    for(int i=0; i<n-1; i++){
        for(int j=i+1; j<n; j++){

			//components of separation vector
            dx=x[i]-x[j]; 
            dy=y[i]-y[j];
            dz=z[i]-z[j];
			
			//Minimum image convention
            dx=dx-l*round(dx/l);
            dy=dy-l*round(dy/l); 
            dz=dz-l*round(dz/l);
			
			//Magnitudes of separation and orientation vectors
            r2=dx*dx+dy*dy+dz*dz;
            r=pow(r2,0.5);
			e1Mag=sqrt(ex[i]*ex[i]+ey[i]*ey[i]+ez[i]*ez[i]);
			e2Mag=sqrt(ex[j]*ex[j]+ey[j]*ey[j]+ez[j]*ez[j]);

            if(r2<rc2){

                //dot products
                dot1=dx*ex[i]+dy*ey[i]+dz*ez[i];
                dot2=dx*ex[j]+dy*ey[j]+dz*ez[j];
                dot12=ex[i]*ex[j]+ey[i]*ey[j]+ez[i]*ez[j]; dot122=dot12*dot12;
                dotSum=dot1+dot2; dotSum2=dotSum*dotSum;
                dotDif=dot1-dot2; dotDif2=dotDif*dotDif;

                dotSByChi=dotSum/(1+chi*dot12);
                dotDByChi=dotDif/(1-chi*dot12);
                dotSByChip=dotSum/(1+chiPrime*dot12);
                dotDByChip=dotDif/(1-chiPrime*dot12);

				dotSByChi2=dotSum2/(1+chi*dot122);
                dotDByChi2=dotDif2/(1-chi*dot122);
                dotSByChip2=dotSum2/(1+chiPrime*dot122);
                dotDByChip2=dotDif2/(1-chiPrime*dot122);

                //calculation of g, gPrime, and others
                g=1.0-(chi/(2.0*r2))*(dotSum*dotSByChi+dotDif*dotDByChi);
                gPrime=1.0-(chiPrime/(2.0*r2))*(dotSum*dotSByChip+dotDif*dotDByChip); //epsilon
                ePrime=1.0/pow(1.0-chi*chi*dot122,0.5);//epsilon single prime

                gPm=pow(gPrime,mu);
                gPm1=pow(gPrime,mu-1);
                ePn=pow(ePrime,nu);
                ePn1=pow(ePrime,nu-1);

                gHalf=pow(g,0.5);

                distF=sigmaS/pow(g,0.5);

                R=(r-distF+sigmaS)/sigmaS;
                R_1=1.0/R;
                R_6=R_1*R_1*R_1*R_1*R_1*R_1;
                R_7=R_6*R_1;                
				R_12=R_6*R_6;
				R_13=R_12*R_1;

                //derivatives of g and gPrime wrt separation
                dgx=-(chi/r2)*(dotSByChi*(ex[i]+ex[j])+dotDByChi*(ex[i]-ex[j]))
                        +(dx*chi/(r2*r2))*(dotSum*dotSByChi+dotDif*dotDByChi);
                dgy=-(chi/r2)*(dotSByChi*(ey[i]+ey[j])+dotDByChi*(ey[i]-ey[j]))
                        +(dy*chi/(r2*r2))*(dotSum*dotSByChi+dotDif*dotDByChi);
                dgz=-(chi/r2)*(dotSByChi*(ez[i]+ez[j])+dotDByChi*(ez[i]-ez[j]))
                        +(dz*chi/(r2*r2))*(dotSum*dotSByChi+dotDif*dotDByChi);

                dgxPrime=-(chiPrime/r2)*(dotSByChip*(ex[i]+ex[j])+dotDByChip*(ex[i]-ex[j]))
                        +(dx*chiPrime/(r2*r2))*(dotSum*dotSByChip+dotDif*dotDByChip);
                dgyPrime=-(chiPrime/r2)*(dotSByChip*(ey[i]+ey[j])+dotDByChip*(ey[i]-ey[j]))
                        +(dy*chiPrime/(r2*r2))*(dotSum*dotSByChip+dotDif*dotDByChip);
                dgzPrime=-(chiPrime/r2)*(dotSByChip*(ez[i]+ez[j])+dotDByChip*(ez[i]-ez[j]))
                        +(dz*chiPrime/(r2*r2))*(dotSum*dotSByChip+dotDif*dotDByChip);

                //derivatives of g and gPrime with respect to orientation
                dgx1=-(chi/2)*((dx/r)*(2*dotSByChi+2*dotDByChi)+chi*ex[j]
                        *(dotDByChi2+dotSByChi2));
                dgy1=-(chi/2)*((dy/r)*(2*dotSByChi+2*dotDByChi)+chi*ey[j]
                        *(dotDByChi2+dotSByChi2));
                dgz1=-(chi/2)*((dz/r)*(2*dotSByChi+2*dotDByChi)+chi*ez[j]
                        *(dotDByChi2+dotSByChi2));

                dgx2=-(chi/2)*((dx/r)*(2*dotSByChi+2*dotDByChi)+chi*ex[i]
                        *(dotDByChi2+dotSByChi2));
                dgy2=-(chi/2)*((dy/r)*(2*dotSByChi+2*dotDByChi)+chi*ey[i]
                        *(dotDByChi2+dotSByChi2));
                dgz2=-(chi/2)*((dz/r)*(2*dotSByChi+2*dotDByChi)+chi*ez[i]
                        *(dotDByChi2+dotSByChi2));

                dgpx1=-(chiPrime/2)*((dx/r)*(2*dotSByChip+2*dotDByChip)+chiPrime*ex[j]
                        *(dotDByChip2+dotSByChip2));
                dgpy1=-(chiPrime/2)*((dy/r)*(2*dotSByChip+2*dotDByChip)+chiPrime*ey[j]
                        *(dotDByChip2+dotSByChip2));
                dgpz1=-(chiPrime/2)*((dz/r)*(2*dotSByChip+2*dotDByChip)+chiPrime*ez[j]
                        *(dotDByChip2+dotSByChip2));

                dgpx2=-(chiPrime/2)*((dx/r)*(2*dotSByChip+2*dotDByChip)+chiPrime*ex[i]
                        *(dotDByChip2+dotSByChip2));
                dgpy2=-(chiPrime/2)*((dy/r)*(2*dotSByChip+2*dotDByChip)+chiPrime*ey[i]
                        *(dotDByChip2+dotSByChip2));
                dgpz2=-(chiPrime/2)*((dz/r)*(2*dotSByChip+2*dotDByChip)+chiPrime*ez[i]
                        *(dotDByChip2+dotSByChip2));

                //derivatives of R with respect to orientations
                drx1=0.5*pow(distF/sigmaS,3)*dgx1;
                dry1=0.5*pow(distF/sigmaS,3)*dgy1;
                drz1=0.5*pow(distF/sigmaS,3)*dgz1;

                drx2=0.5*pow(distF/sigmaS,3)*dgx2;
                dry2=0.5*pow(distF/sigmaS,3)*dgy2;
                drz2=0.5*pow(distF/sigmaS,3)*dgz2;

                //derivatives of epsilon single prime wrt orientation
                depx1=chi*chi*pow(ePrime,3)*ex[j];
                depy1=chi*chi*pow(ePrime,3)*ey[j];
                depz1=chi*chi*pow(ePrime,3)*ez[j];

                depx2=chi*chi*pow(ePrime,3)*ex[i];
                depy2=chi*chi*pow(ePrime,3)*ey[i];
                depz2=chi*chi*pow(ePrime,3)*ez[i];

                //derivatives of epsilon double prime wrt orientation
                dex1=epsilonS*(ePn*mu*gPm1*dgpx1+gPm*nu*ePn1*depx1);
                dey1=epsilonS*(ePn*mu*gPm1*dgpy1+gPm*nu*ePn1*depy1);
                dez1=epsilonS*(ePn*mu*gPm1*dgpz1+gPm*nu*ePn1*depz1);

                dex2=epsilonS*(ePn*mu*gPm1*dgpx2+gPm*nu*ePn1*depx2);
                dey2=epsilonS*(ePn*mu*gPm1*dgpy2+gPm*nu*ePn1*depy2);
                dez2=epsilonS*(ePn*mu*gPm1*dgpz2+gPm*nu*ePn1*depz2);

                //force components
                fxi=-epsilonS*(ePn*gPm*(6*R_7-12*R_13)*(dx/r+(sigmaS/2)
                    /(pow(g,1.5))*dgx)+mu*gPm1*(R_12-R_6)*dgxPrime); 
                fyi=-epsilonS*(ePn*gPm*(6*R_7-12*R_13)*(dy/r+(sigmaS/2)
                    /(pow(g,1.5))*dgy)+mu*gPm1*(R_12-R_6)*dgyPrime);
                fzi=-epsilonS*(ePn*gPm*(6*R_7-12*R_13)*(dz/r+(sigmaS/2)
                    /(pow(g,1.5))*dgz)+mu*gPm1*(R_12-R_6)*dgzPrime);

                //torque components
                gx1=(R_12-R_6)*dex1+g*(6*R_7-12*R_13)*drx1;
                gy1=(R_12-R_6)*dey1+g*(6*R_7-12*R_13)*dry1;
                gz1=(R_12-R_6)*dez1+g*(6*R_7-12*R_13)*drz1;

                gx2=(R_12-R_6)*dex2+g*(6*R_7-12*R_13)*drx2;
                gy2=(R_12-R_6)*dey2+g*(6*R_7-12*R_13)*dry2;
                gz2=(R_12-R_6)*dez2+g*(6*R_7-12*R_13)*drz2;
				
				//Summing forces
                fx[i]=fx[i]+fxi; fx[j]=fx[j]-fxi; 
                fy[i]=fy[i]+fyi; fy[j]=fy[j]-fyi;
                fz[i]=fz[i]+fzi; fz[j]=fz[j]-fzi;

				//Summing torques
                gx[i]=gx[i]-gx1; gx[j]=gx[j]-gx2;
                gy[i]=gy[i]-gy1; gy[j]=gy[j]-gy2;
                gz[i]=gz[i]-gz1; gz[j]=gz[j]-gz2;
				
				//check if particles are over lapping or gorques are NaN
           	    if(r<=0.8 || isnan(gx1)==1 || isnan(gx2)==1){
						cout<<"r="<<r<<"    i="<<i<<"    j="<<j<<endl;
                        cout<<"r1(" << x[i] <<"," << y[i] << ", " << z[i] << ")"<<endl;
                        cout<<"r2(" << x[j] <<"," << y[j] << ", " << z[j] << ")"<<endl;
                        cout<<"e1(" << ex[i] <<","<< ey[i] <<", " << ez[i] << ") MAG: "<<e1Mag<<endl;
                        cout<<"e2(" << ex[j] <<","<< ey[j] << ", " << ez[j]<< ") MAG: "<<e2Mag<<endl;
                        cout<<"g1("<< gx1 <<"," << gy1 << ", " << gz1 << ")"<<endl;
                        cout<<"g2("<< gx2 <<"," << gy2 << ", " << gz2 << ")"<<endl;
						cout<<"g: "<<g<<" dex1: "<<dex1<<" depx1: "<<depx1<<" eprime: " << ePrime <<endl<<endl;
           	    }

				//Calculate potential
                V=V+4.0*epsilonS*ePn*gPm*R_6*(R_6-1.0);

				//Calculate Pressure
                P=P+fxi*dx+fyi*dy+fzi*dz;
            }
        }
    }
	//Pressure converted to Pa -should move this
    P=n*kB*pow(10,-21)*T+P*pow(10,-21)/3;
	P=P/(l*l*l*pow(10,-18));
}

/*Places particles on a cubic lattice with random deviations from 
lattice sites, random orientations, velocities, etc. according to 
specified temperature*/
void init(double x[], double y[], double z[], double vx[],
        double vy[], double vz[], double ux[], double uy[], double uz[],
        double ex[], double ey[], double ez[],
        double m[], double mass, double I, double l, double dt, double temp,
        double kB, int n){

    double sumvx=0.0, sumvy=0.0, sumvz=0.0; //used to set lin mtm = 0
    double sumux=0.0, sumuy=0.0, sumuz=0.0;
    double sumx=0.0, sumy=0.0, sumz=0.0; // for debugging get rid of later
    double sumv2x=0.0, sumv2y=0.0, sumv2z=0.0; //set kinetic energy
    double sumu2x=0.0, sumu2y=0.0, sumu2z=0.0;
    double mag;//unit vector magnitude
	int N, p; double a;
	double fsvx, fsvy, fsvz, fsux, fsuy, fsuz;

    for(int i=0; i<n; i++){

        m[i]=mass;
		
		//Assign random orientation to each molecule
        ex[i]=1.0;//dRand(0,1);
        ey[i]=1.0;//dRand(0,1);
        ez[i]=1.0;//dRand(0,1);
		
		//Make orientation vector unit length
        mag=pow(ex[i]*ex[i]+ey[i]*ey[i]+ez[i]*ez[i], 0.5);
        ex[i]=ex[i]/mag;
        ey[i]=ey[i]/mag;
        ez[i]=ez[i]/mag;
    }

    N=ceil(pow(n,1.0/3.0)); //Third root of n to find # of particles in a direction
    a=l/N; //spacing

    p=0; //number of particles placed

    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            for(int k=0; k<N; k++){
                if(p<n){
					
					//place particles on lattice sites with random deviations
                    x[p]=(i+0.5+dRand(-0.1,0.1))*a;
                    y[p]=(j+0.5+dRand(-0.1,0.1))*a;
                    z[p]=(k+0.5+dRand(-0.1,0.1))*a;
					
					//assign random velocities and ang. velocities
                    vx[p]=dRand(-0.5,0.5);
                    vy[p]=dRand(-0.5,0.5);
                    vz[p]=dRand(-0.5,0.5);

                    ux[p]=dRand(-0.5,0.5);
                    uy[p]=dRand(-0.5,0.5);
                    uz[p]=dRand(-0.5,0.5);
					
					//Sum velocities and squares for energy and mtm
                    sumvx=sumvx+vx[p];
                    sumvy=sumvy+vy[p];
                    sumvz=sumvz+vz[p];

                    sumv2x=sumv2x+pow(vx[p],2);
                    sumv2y=sumv2y+pow(vy[p],2);
                    sumv2z=sumv2z+pow(vz[p],2);

                    sumu2x=sumu2x+pow(ux[p],2);
                    sumu2y=sumu2y+pow(uy[p],2);
                    sumu2z=sumu2z+pow(uz[p],2);
                }
                p++;
            }
        }
    }
	
	//cm velocity of system
    sumvx=sumvx/n; sumvy=sumvy/n; sumvz=sumvz/n;

	//mean squared velocities
    sumv2x=sumv2x/n; sumv2y=sumv2y/n; sumv2z=sumv2z/n;
    sumu2x=sumu2x/n; sumu2y=sumu2y/n; sumu2z=sumu2z/n;

	//calculate scaling factors to set correct temperature
    fsvx=sqrt(kB*temp/sumv2x);
    fsvy=sqrt(kB*temp/sumv2y);
    fsvz=sqrt(kB*temp/sumv2z);

    fsux=sqrt(kB*temp/sumu2x);
    fsuy=sqrt(kB*temp/sumu2y);
    fsuz=sqrt(kB*temp/sumu2z);

    for(int i=0; i<n; i++){

		//scale velocites and ang. velocities
        vx[i]=(vx[i]-sumvx)*fsvx;
        vy[i]=(vy[i]-sumvy)*fsvy;
        vz[i]=(vz[i]-sumvz)*fsvz;

        ux[i]=(ux[i])*fsux;
        uy[i]=(uy[i])*fsuy;
        uz[i]=(uz[i])*fsuz;
    }
}

/*Places particles on a cubic lattice with random deviations from 
lattice sites, random orientations, velocities, etc. according to 
specified temperature THIS INITIALIZES ACCORDING TO VERLET INTEGRATION
ALGORITHM*/
void initVerlet(double x[], double y[], double z[], double vx[],
        double vy[], double vz[], double ux[], double uy[], double uz[],
        double ex[], double ey[], double ez[],
        double m[], double mass, double I, double l, double dt, double temp,
        double kB, int n){

    double sumvx=0.0, sumvy=0.0, sumvz=0.0; //used to set lin mtm = 0
    double sumux=0.0, sumuy=0.0, sumuz=0.0;
    double sumx=0.0, sumy=0.0, sumz=0.0; // for debugging get rid of later
    double sumv2x=0.0, sumv2y=0.0, sumv2z=0.0; //set kinetic energy
    double sumu2x=0.0, sumu2y=0.0, sumu2z=0.0;
    double mag;//unit vector magnitude
	int N, p; double a;
	double fsvx, fsvy, fsvz, fsux, fsuy, fsuz;
	double tvx[n], tvy[n], tvz[n], tux[n], tuy[n], tuz[n];

    for(int i=0; i<n; i++){

        m[i]=mass;
		
		//Assign random orientation to each molecule
        ex[i]=1.0;//dRand(0,1);
        ey[i]=1.0;//dRand(0,1);
        ez[i]=1.0;//dRand(0,1);
		
		//Make orientation vector unit length
        mag=pow(ex[i]*ex[i]+ey[i]*ey[i]+ez[i]*ez[i], 0.5);
        ex[i]=ex[i]/mag;
        ey[i]=ey[i]/mag;
        ez[i]=ez[i]/mag;
    }

    N=ceil(pow(n,1.0/3.0)); //Third root of n to find # of particles in a direction
    a=l/N; //spacing

    p=0; //number of particles placed

    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            for(int k=0; k<N; k++){
                if(p<n){
					
					//place particles on lattice sites with random deviations
                    x[p]=(i+0.5+dRand(-0.1,0.1))*a;
                    y[p]=(j+0.5+dRand(-0.1,0.1))*a;
                    z[p]=(k+0.5+dRand(-0.1,0.1))*a;
					
					//assign random velocities and ang. velocities
                    tvx[p]=dRand(-0.5,0.5);
                    tvy[p]=dRand(-0.5,0.5);
                    tvz[p]=dRand(-0.5,0.5);

                    tux[p]=dRand(-0.5,0.5);
                    tuy[p]=dRand(-0.5,0.5);
                    tuz[p]=dRand(-0.5,0.5);
					
					//Sum velocities and squares for energy and mtm
                    sumvx=sumvx+vx[p];
                    sumvy=sumvy+vy[p];
                    sumvz=sumvz+vz[p];

                    sumv2x=sumv2x+pow(vx[p],2);
                    sumv2y=sumv2y+pow(vy[p],2);
                    sumv2z=sumv2z+pow(vz[p],2);

                    sumu2x=sumu2x+pow(ux[p],2);
                    sumu2y=sumu2y+pow(uy[p],2);
                    sumu2z=sumu2z+pow(uz[p],2);
                }
                p++;
            }
        }
    }
	
	//cm velocity of system
    sumvx=sumvx/n; sumvy=sumvy/n; sumvz=sumvz/n;

	//mean squared velocities
    sumv2x=sumv2x/n; sumv2y=sumv2y/n; sumv2z=sumv2z/n;
    sumu2x=sumu2x/n; sumu2y=sumu2y/n; sumu2z=sumu2z/n;

	//calculate scaling factors to set correct temperature
    fsvx=sqrt(kB*temp/sumv2x);
    fsvy=sqrt(kB*temp/sumv2y);
    fsvz=sqrt(kB*temp/sumv2z);

    fsux=sqrt(kB*temp/sumu2x);
    fsuy=sqrt(kB*temp/sumu2y);
    fsuz=sqrt(kB*temp/sumu2z);

    for(int i=0; i<n; i++){

		//scale velocites and ang. velocities
        tvx[i]=(vx[i]-sumvx)*fsvx;
        tvy[i]=(vy[i]-sumvy)*fsvy;
        tvz[i]=(vz[i]-sumvz)*fsvz;

        tux[i]=(ux[i])*fsux;
        tuy[i]=(uy[i])*fsuy;
        tuz[i]=(uz[i])*fsuz;

		vx[i]=x[i]-dt*tvx[i];
		vy[i]=y[i]-dt*tvy[i];
		vz[i]=z[i]-dt*tvz[i];
		
		//I feel like there might be some problems with
		//this not being unit length but we'll see
		ux[i]=ex[i]-dt*tux[i];
		uy[i]=ey[i]-dt*tuy[i];
		uz[i]=ez[i]-dt*tuz[i];
    }
}

/*Calculates average orientation and scalar order parameter once
I implement that*/
void orientationInfo(double ex[], double ey[], double ez[], int n){
	
	double exAvg, eyAvg, ezAvg;

    for(int i=0; i<n; i++){
        exAvg=exAvg+ex[i];
        eyAvg=eyAvg+ey[i];
        ezAvg=ezAvg+ez[i];
    }

    exAvg=exAvg/n; eyAvg=eyAvg/n; ezAvg=ezAvg/n;
	cout<<"(exAvg,eyAvg,ezAvg)=("<<exAvg<<","<<eyAvg<<","<<ezAvg<<")"<<endl;
}

/*Calculates the pair correlation function of the system -need to 
implement to take the average over the last couple timesteps. 
Example in Frenkel*/
double pairCor(double x[], double y[], double z[], int n, double l){
    int bins=180;
    double histo[bins][2]; //forty bins
    double R=0;
    double dR=1.0/30;//need to change numerator depending on particle!!!
    double dx, dy, dz, r;

    for(int b=0; b<bins; b++){
        histo[b][0]=R;

        for(int i=0; i<n-1; i++){
            for(int j=i+1; j<n; j++){

				//components of separation vector
                dx=x[i]-x[j];
                dy=y[i]-y[j];
                dz=z[i]-z[j];

				//correct for min image convention
                dx=dx-l*round(dx/l); 
                dy=dy-l*round(dy/l);
                dz=dz-l*round(dz/l);
				
				//magnitude of separation
                r=sqrt(dx*dx+dy*dy+dz*dz);
				
				//check if particle in shell
                if(r>R && r<R+dR){
                    histo[b][1]++;
                }
            }
        }
        R=R+dR;
    }

    ofstream o;
    o.open("rpc.data");
    R=0;
    for(int i=1; i<bins; i++){
        histo[i][1]=histo[i][1]*2/(4*3.1415*pow(R,2)*dR*256*256/(l*l*l)); //added factor of two.. need to count each particle
        if(histo[i][1]<1000 && histo[i][1]>0)
            o << histo[i][0] << "\t" << histo[i][1] << "\n";
        R=R+dR;
    }
    o.close();

}

/*Integrates motion of the particles*/
void leapfrog(double x[], double y[], double z[], double vx[],
        double vy[], double vz[], double fx[], double fy[], double fz[],
		double ex[], double ey[], double ez[], double ux[], double uy[],
		double uz[], double gx[], double gy[], double gz[], double mass,
        double I, double& K, double dt, int n, double& sumvx,
        double& sumvy, double& sumvz, double l, int loop){

    double dtSqr=dt*dt;
    double dt2=2*dt;
    double xNEW, yNEW, zNEW, vxi, vyi, vzi;
	double lm, uxi, uyi, uzi;
	
	//reset quantities
    K=0; sumvx=0; sumvy=0; sumvz=0;

	//dot product
	double dot;
	
    for(int i=0; i<n; i++){

		//save old velocities for energy calculation
        vxi=vx[i]; 
        vyi=vy[i];
        vzi=vz[i];
		
		//advance velocites and positions
        vx[i]=vx[i]+dt*fx[i]/mass;
        vy[i]=vy[i]+dt*fy[i]/mass;
        vz[i]=vz[i]+dt*fz[i]/mass;

        x[i]=x[i]+dt*vx[i];
        y[i]=y[i]+dt*vy[i];
        z[i]=z[i]+dt*vz[i];

		//calculate velocity at t
	    vxi=0.5*(vxi+vx[i]);
        vyi=0.5*(vyi+vy[i]);
        vzi=0.5*(vzi+vz[i]);

		//Calculate lagrange multiplier to constrain length
        lm=-2.0*(ux[i]*ex[i]+uy[i]*ey[i]+uz[i]*ez[i]);

		//dot product
		dot=gx[i]*ex[i]+gy[i]*ey[i]+gz[i]*ez[i];
		
		/*Save old velocities for energy calculation
		actually need to calculate the average velocity*/
        uxi=ux[i];
        uyi=uy[i];
        uzi=uz[i];
	
		//Take perpendicular component of the gorque
		gx[i]=gx[i]-(dot)*ex[i];
		gy[i]=gy[i]-(dot)*ey[i];
		gz[i]=gz[i]-(dot)*ez[i];
		
		//Advance bond vector derivatives
        ux[i]=ux[i]+dt*(gx[i]/I)+lm*ex[i];
        uy[i]=uy[i]+dt*(gy[i]/I)+lm*ey[i];
        uz[i]=uz[i]+dt*(gz[i]/I)+lm*ez[i];

		//Advance orentaions
        ex[i]=ex[i]+dt*ux[i];
        ey[i]=ey[i]+dt*uy[i];
        ez[i]=ez[i]+dt*uz[i];

		/*double mag=pow(ex[i]*ex[i]+ey[i]*ey[i]+ez[i]*ez[i], 0.5);

        ex[i]=ex[i]/mag;
        ey[i]=ey[i]/mag;
        ez[i]=ez[i]/mag;*/

		//Calculate total kinetic energy
		K=K+0.5*mass*(vxi*vxi+vyi*vyi+vzi*vzi)+0.5*I*(uxi*uxi+uyi*uyi+uzi*uzi);
    }
}

void orientMag(double ex[], double ey[], double ez[], int n){

	double mag;
	
	for(int i=0;i<n;i++){

		//calculate magnitude of orientation vector
		mag=ex[i]*ex[i]+ey[i]*ey[i]+ez[i]*ez[i];
		
		//print to terminal if mag>1
		if(mag>1.0000001){
			cout<<"warning i: "<<i<<" eMag: "<< mag <<endl;
		}
	}
}

/*Integrates the equations of motion using the verlet algorithm
instead in this case let vx, vy, .. and ux, uy, .. be the old velocities*/
void verlet(double x[], double y[], double z[], double vx[], double vy[],
			double vz[], double fx[], double fy[], double fz[], 
			double ex[], double ey[], double ez[], double ux[],
			double uy[], double uz[], double gx[], double gy[],
			double gz[], double mass, double I, double dT, double& K,
			 int n){

	double xNEW, yNEW, zNEW, exNEW, eyNEW, ezNEW;
	double vxi, vyi, vzi, uxi, uyi, uzi;
	double lm, dot1, dot2;

	for(int i=0; i<n; i++){
		
		xNEW=2*x[i]-vx[i]+dT*dT*fx[i]/mass;
		yNEW=2*y[i]-vy[i]+dT*dT*fy[i]/mass;
		zNEW=2*z[i]-vz[i]+dT*dT*fz[i]/mass;

		vx[i]=x[i];
		vy[i]=y[i];
		vz[i]=z[i];

		x[i]=xNEW;
		y[i]=yNEW;
		z[i]=zNEW;

		vxi=(x[i]-vx[i])/(2*dT)
		vyi=(y[i]-vy[i])/(2*dT)
		vzi=(z[i]-vz[i])/(2*dT)

		exNEW=2*ex[i]-ux[i]+dT*dT*gx[i]/I;
		eyNEW=2*ey[i]-uy[i]+dT*dT*gy[i]/I;
		ezNEW=2*ez[i]-uz[i]+dT*dT*gz[i]/I;

		dot1=ex[i]*exNEW+ey[i]*eyNEW+ez[i]*ezNEW;
		dot2=exNEW*exNEW+eyNEW*eyNEW+ezNEW*ezNEW;

		lm=-dot1+pow(dot1*dot1-dot2+1,0.5);

		exNEW=exNEW+ex[i]*lm;
		eyNEW=eyNEW+ey[i]*lm;
		ezNEW=ezNEW+ez[i]*lm;

		ux[i]=ex[i];
		uy[i]=ey[i];
		uz[i]=ez[i];

		ex[i]=exNEW;
		ey[i]=eyNEW;
		ez[i]=ezNEW;
		
		uxi=(ex[i]-ux[i])/(2*dT);
		uyi=(ey[i]-uy[i])/(2*dT);
		uzi=(ez[i]-uz[i])/(2*dT);

		K=K+0.5*mass*(vxi*vxi+vyi*vyi*+vzi*vzi)+0.5*I*(uxi*uxi+uyi*uyi+uzi*uzi);
	}
}

/*Writes out potential, kinetic, and total energy to a
3-column tab separated data file*/
void writeEnergy(double V, double K, double E, int time){
    ofstream o;
    o.open("energy.data", ios::app);
    o << time << "\t" << V << "\t" << K << "\t" << E << "\n";
    o.close();
}

/*Writes out positions and orientations to a 6-column
tab seperated data file*/
void writeVectors(double x[], double y[], double z[],
    double ex[], double ey[], double ez[], double i, int n){

    ofstream o;
    std::stringstream ss;
    ss<<"vector"<<i<<".data";
    o.open(ss.str().c_str());
    for(int i=0; i<n; i++){
        o<<x[i]<<"\t"<<y[i]<<"\t"<<z[i]<<"\t"<<ex[i]<<
            "\t"<<ey[i]<<"\t"<<ez[i]<<"\n";
    }
}

int main(int argc, char** argv) {
    //Number of particles
    int n=256;
    //Time information
    int tau=2000;//10*pow(10,3); //Number of time steps
    double dT=0.0015;//pow(10,-4); //Length of time step ** used a smaller step
    double T=tau*dT; //Total time
    //Particle info
    double mass=1.0;
    //Storage
    double x[n],y[n],z[n],vx[n],vy[n],vz[n],ex[n],ey[n],ez[n], ux[n], uy[n],
            uz[n],m[n],fx[n],fy[n],fz[n],gx[n],gy[n],gz[n];
    //Simulation box length
    double l=14.92472;
    //Kinetic/Potential/Total Energy;
    double K,V; double E;
    //Temperature
    double temp=1.7;
    //Boltzmann Cons
    double kB=0.0025;
    //momentum
    double sumvx, sumvy, sumvz;
    //pressure
    double P;
    //moment of inertia
   `double I=1.0;

    double sigE=3.0;
	int rand=1;
	cout.precision(17);
	
	/*Initialize and reinitialize until the temperature is acceptable*/
    do {
        srand(time(NULL));
        temp=1.7;
        init(x, y, z, vx, vy, vz, ux, uy, uz, ex, ey, ez, m, mass, I, l, dT, temp, kB,n);
        gb(x, y, z, fx, fy, fz, ex, ey, ez, gx, gy, gz, V, l, P, kB, T, n, sigE, 0);
        halfstep(x, y, z, vx, vy, vz, fx, fy, fz, ex, ey, ez, ux, uy, uz, gx, gy, gz, mass, dT, n, I);
		orientMag(ex,ey,ez,n);
        bCond(x, y, z, l, n);
        gb(x, y, z, fx, fy, fz, ex, ey, ez, gx, gy, gz, V, l, P, kB, T, n, sigE, 1);
        leapfrog(x, y, z, vx, vy, vz, fx, fy, fz, ex, ey, ez, ux, uy, uz, gx, gy, gz, mass, I,
				 K, dT, n, sumvx, sumvy, sumvz, l, 1);
        bCond(x, y, z, l, n);
        temp=2*K/(5*n*kB);
        rand++;
    } while(temp>300);

    for(int i=2; i<tau; i++){

		//Calculate forces and torques, translate and rotate
        gb(x, y, z, fx, fy, fz, ex, ey, ez, gx, gy, gz, V, l, P, kB, T, n, sigE, i);
        leapfrog(x, y, z, vx, vy, vz, fx, fy, fz, ex, ey, ez, ux, uy, uz, gx, gy, gz, mass, I,
				 K, dT, n, sumvx, sumvy, sumvz, l, i); 
        bCond(x, y, z, l, n);
		
		//Calculate total energy
        E=K+V;
	
		//Calculate temperature
        temp=2*K/(5*n*kB); 

        if(i%10==0 || i==2){
            cout << "Loop# " << i << endl;
            cout << "V: " << V << endl;
            cout << "K: " << K << endl;
            cout << "E: " << E << endl;
            cout << "T: " << temp << endl;
            cout << "P: " << P << endl << endl;
			
			orientMag(ex,ey,ez,n);
            writeEnergy(V,K,E,i);
         }
         writeVectors(x,y,z,ex,ey,ez,i,n);
     }
	orientationInfo(ex,ey,ez,n);
    pairCor(x,y,z,n,l);
}
