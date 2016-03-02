/* 
 * File:   main.cpp,
 * Author: bryan
 *
 * Created on January 6, 2016, 3:17 PM
 */

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>

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

void gbForces(double x[], double y[], double z[], double fx[],
        double fy[], double fz[], double ex[], double ey[], double ez[],
        double& V, double l, double& P, double kB, double T, int n){
    
    double mu=2, nu=1;
    double dx, dy, dz;
    double sigmaE=0.0081, sigmaS=0.00027, epsilonE=1.38064851, epsilonS=6.9032426;
    double kappa=sigmaE/sigmaS, kappaPrime=epsilonS/epsilonE;
    double chi=(pow(kappa,2)-1)/(pow(kappa,2)+1);
    double chiPrime=(pow(kappaPrime,1/mu)-1)/(pow(kappaPrime,1/mu)+1);
    double rc=3.25*sigmaS, rc2=rc*rc; //cuttoff
    double dot1, dot2, dot12, dot122, dotSum, dotSum2, dotDif, dotDif2;
    double g, gPrime, gHalf, dgx, dgy, dgz, dgxPrime, dgyPrime, dgzPrime;
    double R, R_1, R_2, R_6, distF;
    double ePrime=0;
    double fxi, fyi, fzi;
    
    V=0;
    P=0;
    
    for(int i=0; i<n; i++){
        fx[i]=0;
        fy[i]=0;
        fz[i]=0;
    }
    
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
            
            if(r2<rc2){
                dot1=dx*ex[i]+dy*ey[i]+dz*ez[i]; 
                dot2=dx*ex[j]+dy*ey[j]+dz*ez[j];
                dot12=ex[i]*ex[j]+ey[i]*ey[j]+ez[i]*ez[j]; dot122=pow(dot12,2);
                dotSum=dot1+dot2; dotSum2=pow(dotSum,2);
                dotDif=dot1-dot2; dotDif2=pow(dotDif,2);

                g=1-(chi/(2*r2))*((dotSum2/(1+chi*dot12))+(dotDif2/(1-chi*dot12)));
                gPrime=1-(chiPrime/(2*r2))*((dotSum2/(1+chiPrime*dot12))+(dotDif2/(1-chiPrime*dot12))); //epsilon
                gHalf=pow(g,0.5);

                distF=sigmaS/gHalf;

                R=(r-distF+sigmaS)/sigmaS;
                R_1=1/R;
                R_2=R_1*R_1;
                R_6=R_2*R_2*R_2;

                ePrime=1/pow(1-chi*chi*dot122,0.5);

                dgx=-(chi/r2)*((dotSum/(1+chi*dot12))*(ex[i]+ex[j])+(dotDif/(1-chi*dot12))
                    *(ex[i]-ex[j]))+(dx*chi/(r2*r2))*(dotSum2/(1+chi*dot12)+dotDif2/(1-chi*dot12));
                dgy=-(chi/r2)*((dotSum/(1+chi*dot12))*(ey[i]+ey[j])+(dotDif/(1-chi*dot12))
                    *(ey[i]-ey[j]))+(dy*chi/(r2*r2))*(dotSum2/(1+chi*dot12)+dotDif2/(1-chi*dot12));
                dgz=-(chi/r2)*((dotSum/(1+chi*dot12))*(ez[i]+ez[j])+(dotDif/(1-chi*dot12))
                    *(ez[i]-ez[j]))+(dz*chi/(r2*r2))*(dotSum2/(1+chi*dot12)+dotDif2/(1-chi*dot12));

                dgxPrime=-(chiPrime/r2)*((dotSum/(1+chiPrime*dot12))*(ex[i]+ex[j])+(dotDif/(1-chiPrime*dot12))
                    *(ex[i]-ex[j]))+(dx*chiPrime/(r2*r2))*(dotSum2/(1+chiPrime*dot12)+dotDif2/(1-chiPrime*dot12));
                dgyPrime=-(chiPrime/r2)*((dotSum/(1+chiPrime*dot12))*(ey[i]+ey[j])+(dotDif/(1-chiPrime*dot12))
                    *(ey[i]-ey[j]))+(dy*chiPrime/(r2*r2))*(dotSum2/(1+chiPrime*dot12)+dotDif2/(1-chiPrime*dot12));
                dgzPrime=-(chiPrime/r2)*((dotSum/(1+chiPrime*dot12))*(ez[i]+ez[j])+(dotDif/(1-chiPrime*dot12))
                    *(ez[i]-ez[j]))+(dz*chiPrime/(r2*r2))*(dotSum2/(1+chiPrime*dot12)+dotDif2/(1-chiPrime*dot12));
                
                fxi=-epsilonS*(pow(ePrime,nu)*pow(gPrime,mu)*R_6*R_1*(6-12*R_6)*(dx/r+(sigmaS/2)/(gHalf*gHalf*gHalf)*dgx)+mu*pow(gPrime,mu-1)*R_6*(R_6-1)*dgxPrime); //forces between the pairs
                fyi=-epsilonS*(pow(ePrime,nu)*pow(gPrime,mu)*R_6*R_1*(6-12*R_6)*(dy/r+(sigmaS/2)/(gHalf*gHalf*gHalf)*dgy)+mu*pow(gPrime,mu-1)*R_6*(R_6-1)*dgyPrime);
                fzi=-epsilonS*(pow(ePrime,nu)*pow(gPrime,mu)*R_6*R_1*(6-12*R_6)*(dz/r+(sigmaS/2)/(gHalf*gHalf*gHalf)*dgz)+mu*pow(gPrime,mu-1)*R_6*(R_6-1)*dgzPrime);
                
                fx[i]=fx[i]+fxi; fx[j]=fx[j]-fxi; //total force on particle
                fy[i]=fy[i]+fyi; fy[j]=fy[j]-fyi;
                fz[i]=fz[i]+fzi; fz[j]=fz[j]-fzi;
                V=V+4.0*epsilonS*pow(ePrime,nu)*pow(gPrime,mu)*R_6*(R_6-1.0);
                P=P+fxi*dx+fyi*dy+fzi*dz;//pressure
            }
        }
    }
//    cout<<fx[135]<<endl;
    P=n*kB*pow(10,-21)*T+P*pow(10,-21)/3; P=P/(l*l*l*pow(10,-18));
    //converted to Pa
    
}

void init(double x[], double y[], double z[], double vx[],
        double vy[], double vz[], double ex[], double ey[], double ez[],
        double m[], double mass, double l, double dt, double temp, int n){
    
    double sumvx=0.0, sumvy=0.0, sumvz=0.0; //used to set lin mtm = 0
    double sumx=0.0, sumy=0.0, sumz=0.0; // for debugging get rid of later
    double sumv2x=0.0, sumv2y=0.0, sumv2z=0.0; //set kinetic energy
    
    for(int i=0; i<n; i++){
        m[i]=mass;
        ex[i]=0;
        ey[i]=0;
        ez[i]=1;
    }
    
    int N=ceil(pow(n,1.0/3.0)); //Third root of n to find # of particles in a direction
    double a=l/N; //spacing
        
    int p=0; //number of particles placed
    
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            for(int k=0; k<N; k++){
                if(p<n){
                    x[p]=(i+0.5)*a;
                    y[p]=(j+0.5)*a;
                    z[p]=(k+0.5)*a;
                
                    vx[p]=dRand(-0.5,0.5); 
                    vy[p]=dRand(-0.5,0.5); 
                    vz[p]=dRand(-0.5,0.5); 
                    
                    sumvx=sumvx+vx[p];
                    sumvy=sumvy+vy[p];
                    sumvz=sumvz+vz[p];
                    
                    sumv2x=sumv2x+pow(vx[p],2);
                    sumv2y=sumv2y+pow(vy[p],2);
                    sumv2z=sumv2z+pow(vz[p],2);
                }
                p++;
            }
        }
    }//Places particle in lattice & assigns random velocities
    
    sumvx=sumvx/n; sumvy=sumvy/n; sumvz=sumvz/n; //cm velocities
    sumv2x=sumv2x/n; sumv2y=sumv2y/n; sumv2z=sumv2z/n; //mean-squared velocities
    
    
    double fsx=sqrt(temp/sumv2x);
    double fsy=sqrt(temp/sumv2y);
    double fsz=sqrt(temp/sumv2z);
    

    for(int i=0; i<n; i++){
        vx[i]=(vx[i]-sumvx)*fsx;
        vy[i]=(vy[i]-sumvy)*fsy;
        vz[i]=(vz[i]-sumvz)*fsz;
//        cout<<vx[i]<<endl;
    }
}

double pairCor(double x[], double y[], double z[], int n, double l){
    int bins=160;
    double histo[bins][2]; //forty bins
    double R,dR=0.00034/40;
    double dx, dy, dz;
    
    for(int b=0; b<bins; b++){
        histo[b][0]=R/0.00034;
        
        for(int i=0; i<n-1; i++){
            for(int j=i+1; j<n; j++){


                dx=x[i]-x[j]; //components of distance vector
                dy=y[i]-y[j];
                dz=z[i]-z[j];

                dx=dx-l*round(dx/l); //correct for min image convention
                dy=dy-l*round(dy/l); //from frenkel... we'll see how this goes
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
    o.open("pair-correlation.data");
    R=dR;
    for(int i=1; i<bins; i++){
        histo[i][1]=histo[i][1]*2/(4*3.1415*pow(R,2)*dR*864*864/(pow(10.229*0.00034,3))); //added factor of two to scale... but i'll need to look at this again
        if(histo[i][1]<20 && histo[i][1]>0)
            o << histo[i][0] << "\t" << histo[i][1] << "\n";
        R=R+dR;
    }
    o.close();
    
}

void verletLeapfrog(double x[], double y[], double z[], double vx[],
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
//    cout<< "x: " << x[135]<< " vx: " << vx[135] << " fx: " << fx[135] <<endl<<endl;
    //cout << "(" << sumvx << "," << sumvy << "," << sumvz << ")" << endl;
    K=0.5*mass*K; //Just assuming one mass for now
}

void writeXYZ(double x[], double y[], double z[], int n){
    ofstream o;
    o.open("test.xyz",ios::app);
    double t,j,k;
    
    o << 863 << endl;
    
    for(int i=0; i<n; i++){

        t=x[i]*10000;
        j=y[i]*10000;
        k=z[i]*10000;
        o << "Ar" << "\t" <<  t << "\t" << j << "\t" << k << "\n";
        
    }
    o.close();
}

int main(int argc, char** argv) {
    //Number of particles
    int n=864;
    //Time information
    int tau=pow(10,4); //Number of time steps
    double dT=pow(10,-11); //Length of time step ** used a smaller step
    double T=tau*dT; //Total time
    //Particle info
    double mass=4.27*pow(10,-10);
    //Storage
    double x[n],y[n],z[n],vx[n],vy[n],vz[n],ex[n],ey[n],ez[n],m[n],
            fx[n], fy[n], fz[n];
    //Simulation box length
    double l=13.92477*0.00027;
    //Kinetic/Potential/Total Energy;
    double K,V; double E;
    //Temperature
    double temp=390;
    //Boltzmann Cons     
    double kB=0.0138064852;
    //momentum
    double sumvx, sumvy, sumvz;
    //pressure
    double P;
    
    init(x, y, z, vx, vy, vz, ex, ey, ez, m, mass, l, dT, temp, n); 
    writeXYZ(x, y, z, n);
    gbForces(x, y, z, fx, fy, fz, ex, ey, ez, V, l, P, kB, T, n);
//    cout<< "x: " << x[135]<< " vx: " << vx[135] << " fx: " << fx[135] <<endl<<endl;
    halfstep(x, y, z, vx, vy, vz,fx, fy, fz, mass, dT, n);
//    cout<< "x: " << x[135]<< " vx: " << vx[135] << " fx: " << fx[135] <<endl<<endl;
    bCond(x, y, z, l, n);
    writeXYZ(x, y, z, n);

    for(int i=2; i<tau; i++){
        gbForces(x, y, z, fx, fy, fz, ex, ey, ez, V, l, P, kB, T, n);
        
        verletLeapfrog(x, y, z, vx, vy, vz, fx, fy, fz, mass, K, dT, n, 
                    sumvx, sumvy, sumvz, l, i);
        
        bCond(x, y, z, l, n);
        
        writeXYZ(x,y,z,n);
        
        E=K+V; //in scaled units
        temp=2*K/(3*n*kB); //in kelvin
        
        if(i%100==0){            
            cout << "Loop# " << i << endl;
            cout << "V: " << V << endl;
            cout << "K: " << K << endl;
            cout << "E: " << E << endl;
            cout << "T: " << temp << endl;
            cout << "P: " << P << endl << endl;
            
            // if(E>0){
            //      largeForce(x,y,z,fx,fy,fz,E,n);
            // }
         }
     }
     pairCor(x,y,z,n,l);
}

