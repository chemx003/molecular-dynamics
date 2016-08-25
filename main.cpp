/*
 * File:   main.cpp,
 * Author: Bryan
 *
 * Created on January 6, 2016, 3:17 PM
 */

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <time.h>
#include <sstream>

using namespace std;

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

void largeForce(double x[], double y[], double z[], double fx[], double fy[],
        double fz[], double E, int n){
    ofstream o;
    o.open("large_force.dat");
    o << "#" << "\t" << "x" << "\t" << "fx" << "\t" << E << endl;

    for(int i=0; i<n; i++){
        if(abs(fx[i])>0){
            o << i << "\t" << x[i] << "\t" << fx[i] << endl;
        }
    }
    o.close();
}

void lfOrient(double ex[], double ey[], double ez[], double ux[],
        double uy[], double uz[], double gx[], double gy[], double gz[],
        double x[], double y[], double z[], double I, double& K, double dt,
        int n, double l, int loop){
    double lm, uxi, uyi, uzi;


    for(int i=0; i<n; i++){

        lm=-2*(ux[i]*ex[i]+uy[i]*ey[i]+uz[i]*ez[i]);

        uxi=ux[i]; //saving old velocities for energy calc.
        uyi=uy[i];
        uzi=uz[i];

		gx[i]=gx[i]-(gx[i]*ex[i]+gy[i]*ey[i]+gz[i]*ez[i])*ex[i];
		gy[i]=gy[i]-(gx[i]*ex[i]+gy[i]*ey[i]+gz[i]*ez[i])*ey[i];
		gz[i]=gz[i]-(gx[i]*ex[i]+gy[i]*ey[i]+gz[i]*ez[i])*ez[i];

        ux[i]=ux[i]+dt*(gx[i]/I)+lm*ex[i];
        uy[i]=uy[i]+dt*(gy[i]/I)+lm*ey[i];
        uz[i]=uz[i]+dt*(gz[i]/I)+lm*ez[i];

        ex[i]=ex[i]+dt*ux[i];
        ey[i]=ey[i]+dt*uy[i];
        ez[i]=ez[i]+dt*uz[i];

        double mag=pow(ex[i]*ex[i]+ey[i]*ey[i]+ez[i]*ez[i], 0.5);

        ex[i]=ex[i]/mag;
        ey[i]=ey[i]/mag;
        ez[i]=ez[i]/mag;//just here to test remove and debug latera
		
		K=K+0.5*I*(uxi*uxi+uyi*uyi+uzi*uzi);

        //if(i==1)
        //cout<<"("<<ex[i]<<","<<ey[i]<<","<<ez[i]<<")"<<endl<<endl;

        //if(mag>=1){
        //   cout<< i << " (" << mag << ")" << endl;
        //}
        //don't forget to calculate energy from rotation
    }
}

double dRand(double dMin, double dMax){
    double d = (double)rand()/RAND_MAX;
    return dMin + d*(dMax-dMin);
}

void halfstep(double x[], double y[], double z[], double vx[], double vy[],
        double vz[], double fx[], double fy[], double fz[], double mass,
        double dt, int n){

    for(int i=0; i<n; i++){
        vx[i]=vx[i]+0.5*dt*fx[i]/mass;
        vy[i]=vy[i]+0.5*dt*fy[i]/mass;
        vz[i]=vz[i]+0.5*dt*fz[i]/mass;
    }

    for(int i=0; i<n; i++){
        x[i]=x[i]+dt*vx[i];
        y[i]=y[i]+dt*vy[i];
        z[i]=z[i]+dt*vz[i];
    }
}

void gb(double x[], double y[], double z[], double fx[],
        double fy[], double fz[], double ex[], double ey[], double ez[],
        double gx[], double gy[], double gz[], double& V, double l,
        double& P, double kB, double T, int n, double sigE, int loop){
    double mu=2.0, nu=1.0;
    double dx, dy, dz;
    double sigmaE=sigE, sigmaS=1.0, epsilonE=5.0, epsilonS=1.0;
    double kappa=sigmaE/sigmaS, kappaPrime=epsilonE/epsilonS;
    double chi=(pow(kappa,2.0)-1.0)/(pow(kappa,2.0)+1.0);
    double chiPrime=(pow(kappaPrime,1.0/mu)-1.0)/(pow(kappaPrime,1.0/mu)+1.0);
    double rc=3.25*sigmaS, rc2=rc*rc; //cuttoff
    double dot1, dot2, dot12, dot122, dotSum, dotSum2, dotDif, dotDif2;
    double g, gPrime, gHalf, dgx, dgy, dgz, dgxPrime, dgyPrime, dgzPrime;
    double R, R_1, R_2, R_6, R_7, R_12, R_13, distF;
    double ePrime, ePn, ePn1, gPm, gPm1;
    double fxi, fyi, fzi, gx1, gy1, gz1, gx2, gy2, gz2;
    double dotSByChi, dotDByChi, dotSByChip, dotDByChip;
	double dotSByChi2, dotDByChi2, dotSByChip2, dotDByChip2;
    double dex1, dey1, dez1, dex2, dey2, dez2;
    double depx1, depy1, depz1, depx2, depy2, depz2;//derivatives of epsilon single prime wrt orientation
    double drx1,dry1,drz1,drx2,dry2,drz2;//derivatives of r wrt orientation
    double dgx1,dgy1,dgz1,dgx2,dgy2,dgz2;
    double dgpx1,dgpy1,dgpz1,dgpx2,dgpy2,dgpz2;

    V=0;
    P=0;

    for(int i=0; i<n; i++){
        fx[i]=0; gx[i]=0;
        fy[i]=0; gy[i]=0;
        fz[i]=0; gz[i]=0;
    }

    int lrg=0;

    for(int i=0; i<n-1; i++){
        for(int j=i+1; j<n; j++){

            dx=x[i]-x[j]; //components of distance vector
            dy=y[i]-y[j];
            dz=z[i]-z[j];

            dx=dx-l*round(dx/l); //correct for min image convention
            dy=dy-l*round(dy/l); //from frenkel... we'll see how this goes
            dz=dz-l*round(dz/l);

            double r2=dx*dx+dy*dy+dz*dz;
            double r=pow(r2,0.5);

			double e1Mag=sqrt(ex[i]*ex[i]+ey[i]*ey[i]+ez[i]*ez[i]);
			double e2Mag=sqrt(ex[j]*ex[j]+ey[j]*ey[j]+ez[j]*ez[j]);

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
                R_6=R_1*R_1*R_1*R_1*R_1*R_1;    R_12=R_6*R_6;
                R_7=R_6*R_1;                    R_13=R_12*R_1;

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

                dgpx2=-(chiPrime/2)*((dx/r)*(2*dotSByChip+2*dotDByChip)+chiPrime*ex[j]
                        *(dotDByChip2+dotSByChip2));
                dgpy2=-(chiPrime/2)*((dy/r)*(2*dotSByChip+2*dotDByChip)+chiPrime*ey[j]
                        *(dotDByChip2+dotSByChip2));
                dgpz2=-(chiPrime/2)*((dz/r)*(2*dotSByChip+2*dotDByChip)+chiPrime*ez[j]
                        *(dotDByChip2+dotSByChip2));// I need to make this look cleaner!!!

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
                    /(pow(g,1.5))*dgx)+mu*gPm1*(R_12-R_6)*dgxPrime); //forces between the pairs
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

                fx[i]=fx[i]+fxi; fx[j]=fx[j]-fxi; //total force on particle
                fy[i]=fy[i]+fyi; fy[j]=fy[j]-fyi;
                fz[i]=fz[i]+fzi; fz[j]=fz[j]-fzi;

                gx[i]=gx[i]-gx1; gx[j]=gx[j]-gx2;
                gy[i]=gy[i]-gy1; gy[j]=gy[j]-gy2;
                gz[i]=gz[i]-gz1; gz[j]=gz[j]-gz2;

                //if(i==1){
                //    cout<<gx[i]<<endl;
                //}

           	    if(r<=1 || isnan(gx1)==1 || isnan(gx2)==1){
						cout<<"r="<<r<<"    i="<<i<<"    j="<<j<<endl;
                        cout<<"r1(" << x[i] <<"," << y[i] << ", " << z[i] << ")"<<endl;
                        cout<<"r2(" << x[j] <<"," << y[j] << ", " << z[j] << ")"<<endl;
                        cout<<"e1(" << ex[i] <<","<< ey[i] <<", " << ez[i] << ") MAG: "<<e1Mag<<endl;
                        cout<<"e2(" << ex[j] <<","<< ey[j] << ", " << ez[j]<< ") MAG: "<<e2Mag<<endl;
                        cout<<"g1("<< gx1 <<"," << gy1 << ", " << gz1 << ")"<<endl;
                        cout<<"g2("<< gx2 <<"," << gy2 << ", " << gz2 << ")"<<endl;
						cout<<"g: "<<g<<" dex1: "<<dex1<<" depx1: "<<depx1<<" eprime: " << ePrime <<endl<<endl;
           	    }

                V=V+4.0*epsilonS*ePn*gPm*R_6*(R_6-1.0);
                //cout<<V<<endl;
                P=P+fxi*dx+fyi*dy+fzi*dz;//pressure
            }
        }
    }
//    cout<<"LARGE INIDENTS: "<<lrg<<endl<< endl;
//    double sumx=0;
//    for(int i=0; i<n ; i++){
//        sumx=sumx+fx[i];
//    }
//    cout<<sumx<<endl; //Checking if sum of forces = 0;
    P=n*kB*pow(10,-21)*T+P*pow(10,-21)/3; P=P/(l*l*l*pow(10,-18));
    //converted to Pa

}

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

    for(int i=0; i<n; i++){
        m[i]=mass;
        ex[i]=dRand(0,1);
        ey[i]=dRand(0,1);
        ez[i]=dRand(0,1);
        mag=pow(ex[i]*ex[i]+ey[i]*ey[i]+ez[i]*ez[i], 0.5);
        ex[i]=ex[i]/mag;
        ey[i]=ey[i]/mag;
        ez[i]=ez[i]/mag;
    }

    int N=ceil(pow(n,1.0/3.0)); //Third root of n to find # of particles in a direction
    double a=l/N; //spacing

    int p=0; //number of particles placed

    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            for(int k=0; k<N; k++){
                if(p<n){
                    x[p]=(i+0.5+dRand(-0.1,0.1))*a;
                    y[p]=(j+0.5+dRand(-0.1,0.1))*a;
                    z[p]=(k+0.5+dRand(-0.1,0.1))*a;

                    vx[p]=dRand(-0.5,0.5);
                    vy[p]=dRand(-0.5,0.5);
                    vz[p]=dRand(-0.5,0.5);

                    ux[p]=dRand(-0.5,0.5);
                    uy[p]=dRand(-0.5,0.5);
                    uz[p]=dRand(-0.5,0.5);

                    sumvx=sumvx+vx[p];
                    sumvy=sumvy+vy[p];
                    sumvz=sumvz+vz[p];

                    sumux=sumux+ux[p];
                    sumuy=sumuy+uy[p];
                    sumuz=sumuz+uz[p];

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
    }//Places particle in lattice & assigns random velocities

    sumvx=sumvx/n; sumvy=sumvy/n; sumvz=sumvz/n; //cm velocities
    sumv2x=sumv2x/n; sumv2y=sumv2y/n; sumv2z=sumv2z/n; //mean-squared velocities
    sumu2x=sumu2x/n; sumu2y=sumu2y/n; sumu2z=sumu2z/n;

    double fsvx=sqrt(kB*temp/sumv2x);
    double fsvy=sqrt(kB*temp/sumv2y);
    double fsvz=sqrt(kB*temp/sumv2z);

    double fsux=sqrt(kB*temp/sumu2x);
    double fsuy=sqrt(kB*temp/sumu2y);
    double fsuz=sqrt(kB*temp/sumu2z);

    for(int i=0; i<n; i++){
        vx[i]=(vx[i]-sumvx)*fsvx;
        vy[i]=(vy[i]-sumvy)*fsvy;
        vz[i]=(vz[i]-sumvz)*fsvz;

        ux[i]=(ux[i])*fsux;
        uy[i]=(uy[i])*fsuy;
        uz[i]=(uz[i])*fsuz;
    }
}

void rescale(double vx[],double vy[], double vz[], double ux[],
        double uy[], double uz[], double temp, double kB, int n){

    double sumvx=0.0, sumvy=0.0, sumvz=0.0; //used to set lin mtm = 0
    double sumux=0.0, sumuy=0.0, sumuz=0.0;
    double sumx=0.0, sumy=0.0, sumz=0.0; // for debugging get rid of later
    double sumv2x=0.0, sumv2y=0.0, sumv2z=0.0; //set kinetic energy
    double sumu2x=0.0, sumu2y=0.0, sumu2z=0.0;
    double mag;//unit vector magnitude

    for(int p=0; p<n; p++){

                    sumvx=sumvx+vx[p];
                    sumvy=sumvy+vy[p];
                    sumvz=sumvz+vz[p];

                    sumux=sumux+ux[p];
                    sumuy=sumuy+uy[p];
                    sumuz=sumuz+uz[p];

                    sumv2x=sumv2x+pow(vx[p],2);
                    sumv2y=sumv2y+pow(vy[p],2);
                    sumv2z=sumv2z+pow(vz[p],2);

                    sumu2x=sumu2x+pow(ux[p],2);
                    sumu2y=sumu2y+pow(uy[p],2);
                    sumu2z=sumu2z+pow(uz[p],2);
                }

    sumvx=sumvx/n; sumvy=sumvy/n; sumvz=sumvz/n; //cm velocities
    sumv2x=sumv2x/n; sumv2y=sumv2y/n; sumv2z=sumv2z/n; //mean-squared velocities
    sumu2x=sumu2x/n; sumu2y=sumu2y/n; sumu2z=sumu2z/n;

    double fsvx=sqrt(kB*temp/sumv2x);
    double fsvy=sqrt(kB*temp/sumv2y);
    double fsvz=sqrt(kB*temp/sumv2z);

    double fsux=sqrt(kB*temp/sumu2x);
    double fsuy=sqrt(kB*temp/sumu2y);
    double fsuz=sqrt(kB*temp/sumu2z);

    for(int i=0; i<n; i++){
        vx[i]=(vx[i]-sumvx)*fsvx;
        vy[i]=(vy[i]-sumvy)*fsvy;
        vz[i]=(vz[i]-sumvz)*fsvz;

        ux[i]=(ux[i])*fsux;
        uy[i]=(uy[i])*fsuy;
        uz[i]=(uz[i])*fsuz;
    }
}

double pairCor(double x[], double y[], double z[], int n, double l){
    int bins=180;
    double histo[bins][2]; //forty bins
    double R=0;
    double dR=1.0/30;//need to change numerator depending on particle!!!
    double dx, dy, dz;

    for(int b=0; b<bins; b++){
        histo[b][0]=R;

        for(int i=0; i<n-1; i++){
            for(int j=i+1; j<n; j++){


                dx=x[i]-x[j]; //components of distance vector
                dy=y[i]-y[j];
                dz=z[i]-z[j];

                dx=dx-l*round(dx/l); //correct for min image convention
                dy=dy-l*round(dy/l);
                dz=dz-l*round(dz/l);

                double r=sqrt(dx*dx+dy*dy+dz*dz);

                if(r>R and r<R+dR){
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

void leapfrog(double x[], double y[], double z[], double vx[],
        double vy[], double vz[], double fx[], double fy[], double fz[],
        double mass, double& K, double dt, int n, double& sumvx,
        double& sumvy, double& sumvz, double l, int loop){

    double dtSqr=dt*dt;
    double dt2=2*dt;
    double xNEW, yNEW, zNEW, vxi, vyi, vzi;
    K=0; sumvx=0; sumvy=0; sumvz=0;

    for(int i=0; i<n; i++){

        vxi=vx[i]; //save old velocities for energy calculation
        vyi=vy[i];
        vzi=vz[i];

        vx[i]=vx[i]+dt*fx[i]/mass;
        vy[i]=vy[i]+dt*fy[i]/mass;
        vz[i]=vz[i]+dt*fz[i]/mass;

        vxi=0.5*(vxi+vx[i]);
        vyi=0.5*(vyi+vy[i]);
        vzi=0.5*(vzi+vz[i]);

        x[i]=x[i]+dt*vx[i];
        y[i]=y[i]+dt*vy[i];
        z[i]=z[i]+dt*vz[i];

        K=K+vxi*vxi+vyi*vyi+vzi*vzi;

        sumvx=sumvx+vxi;//velocity of the center of mass
        sumvy=sumvy+vyi;
        sumvz=sumvz+vzi;
    }
    K=0.5*mass*K; //Just assuming one mass for now
}

void track(int p, double x[], double y[], double z[],double loop){
    cout<<loop<<" ("<<x[p]<<","<< y[p]<< "," << z[p]<<")"<<endl;
}

void writeXYZ(double x[], double y[], double z[], int n){
    ofstream o;
    o.open("random-testing.xyz",ios::app); //I should make this a setting in main())
    double t,j,k;

    o << 255 << endl;

    for(int i=0; i<n; i++){

        t=x[i]*100;
        j=y[i]*100;
    }
    o.close();
}

void writeEnergy(double V, double K, double E, int time){
    ofstream o;
    o.open("energy.data", ios::app);
    o << time << "\t" << V << "\t" << K << "\t" << E << "\n";
    o.close();
    }
    o.close();
}

void writeEnergy(double V, double K, double E, int time){
    ofstream o;
    o.open("energy.data", ios::app);
    o << time << "\t" << V << "\t" << K << "\t" << E << "\n";
    o.close();
}

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
    int tau=10000;//10*pow(10,3); //Number of time steps
    double dT=0.0015;//pow(10,-4); //Length of time step ** used a smaller step
    double T=tau*dT; //Total time
    //Particle info
    double mass=1;
    //Storage
    double x[n],y[n],z[n],vx[n],vy[n],vz[n],ex[n],ey[n],ez[n], ux[n], uy[n],
            uz[n],m[n],fx[n],fy[n],fz[n],gx[n],gy[n],gz[n];
    //Simulation box length
    double l=14.92472; //scaled density of 0.2
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
    double I=1.0;
    double sigE=3.0;

    int rand=1.0;//
    do {//
        //Random seed;
        srand(rand*time(NULL));//time(NULL)
        temp=1.7;//
        init(x, y, z, vx, vy, vz, ux, uy, uz, ex, ey, ez, m, mass, I, l, dT, temp, kB,n);
        //writeXYZ(x, y, z, n);
        gb(x, y, z, fx, fy, fz, ex, ey, ez, gx, gy, gz, V, l, P, kB, T, n, sigE, 0);
        halfstep(x, y, z, vx, vy, vz,fx, fy, fz, mass, dT, n);
        bCond(x, y, z, l, n);
        //writeXYZ(x, y, z, n);
        gb(x, y, z, fx, fy, fz, ex, ey, ez, gx, gy, gz, V, l, P, kB, T, n, sigE, 1);//
        leapfrog(x, y, z, vx, vy, vz, fx, fy, fz, mass, K, dT, n, sumvx, sumvy, sumvz, l, 1);//
        lfOrient(ex,ey,ez,ux,uy,uz,gx,gy,gz,x,y,z,I,K,dT,n,l,1);
        bCond(x, y, z, l, n);//
        temp=2*K/(5*n*kB);//
        cout<<temp<<endl;//
        rand++;//
    } while(temp>300);//


    for(int i=2; i<tau; i++){
        gb(x, y, z, fx, fy, fz, ex, ey, ez, gx, gy, gz, V, l, P, kB, T, n, sigE, i);
		if(i<5000){
        	leapfrog(x, y, z, vx, vy, vz, fx, fy, fz, mass, K, dT, n, sumvx, sumvy, sumvz, l, i);
		}
        lfOrient(ex,ey,ez,ux,uy,uz,gx,gy,gz,x,y,z,I,K,dT,n,l,i);
        bCond(x, y, z, l, n);
        //rescale(vx,vy,vz,ux,uy,uz,1.7,kB,n);
        //writeXYZ(x,y,z,n);
        E=K+V; //in scaled units
        temp=2*K/(5*n*kB); //in kelvin kg*m^2/s^2 -15*-6^2/-3^2  /-21

        if(i%100==0 || i==2){
            cout << "Loop# " << i << endl;
            cout << "V: " << V << endl;
            cout << "K: " << K << endl;
            cout << "E: " << E << endl;
            cout << "T: " << temp << endl;
            cout << "P: " << P << endl << endl;

            writeEnergy(V,K,E,i);
            // if(E>0){
            //      largeForce(x,y,z,fx,fy,fz,E,n);
            // }
         }
         writeVectors(x,y,z,ex,ey,ez,i,n);
     }

    double exAvg, eyAvg, ezAvg;
    for(int i=0; i<n; i++){
        exAvg=exAvg+ex[i];
        eyAvg=eyAvg+ey[i];
        ezAvg=ezAvg+ez[i];
    }
    exAvg=exAvg/n; eyAvg=eyAvg/n; ezAvg=ezAvg/n;
    cout<<"(exAvg,eyAvg,ezAvg)=("<<exAvg<<","<<eyAvg<<","<<ezAvg<<")";
    pairCor(x,y,z,n,l);
}
