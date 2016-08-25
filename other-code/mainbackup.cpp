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

void bCond(double x[], double y[], double z[], double xOLD[],
        double yOLD[], double zOLD[], double l, int n){
    
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

double dRand(double dMin, double dMax){
    double d = (double)rand()/RAND_MAX;
    return dMin + d*(dMax-dMin);
}

void forces(double x[], double y[], double z[], double fx[],
        double fy[], double fz[], double& V, double l, int n){
    
    double sigma=2.341, epsilon=0.1, rc=3*sigma;
    double sigma2=sigma*sigma, rc2=rc*rc;
    double ecut;//discrepancy in sources
    double dx, dy, dz;
    double fr2, fr6, frp;
    double fxi, fyi, fzi;
    
    V=0;
    
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
            
            if(r2<rc2){
                fr2=sigma2/r2;
                fr6=pow(fr2,3);
                frp=48*epsilon*fr6*(fr6-0.5)/r2;
                
                fxi=frp*dx;
                fyi=frp*dy;
                fzi=frp*dz;
                
                fx[i]=fx[i]+fxi; fx[j]=fx[j]-fxi;
                fy[i]=fy[i]+fyi; fy[j]=fy[j]-fyi;
                fz[i]=fz[i]+fzi; fz[j]=fz[j]-fzi;
                
                V=V+4.0*epsilon*fr6*(fr6-1.0);
            }
        }
    }
    
}

void init(double x[], double y[], double z[], double xOLD[],
        double yOLD[], double zOLD[], double m[], double mass, 
        double l, double dt, double temp, int n){
    
    double sumvx=0.0, sumvy=0.0, sumvz=0.0; //used to set lin mtm = 0
    double sumv2x=0.0, sumv2y=0.0, sumv2z=0.0; //set kinetic energy
    
    for(int i=0; i<n; i++){
        m[i]=mass;
    }//initialize mass for a system of 1 gas
    
    int N=ceil(pow(n,1.0/3.0)); //Third root of n to find # of particles in a direction
    double a=l/N; //spacing
        
    int p=0; //number of particles placed
    double vx[n]; //temporary holders for velocity
    double vy[n];
    double vz[n];
    
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
        vz[i]=(vy[i]-sumvz)*fsz;

        xOLD[i]=x[i]-vx[i]*dt;
        yOLD[i]=y[i]-vy[i]*dt;
        zOLD[i]=z[i]-vz[i]*dt;
    }           
}

void verlet(double x[], double y[], double z[], double xOLD[],
        double yOLD[], double zOLD[], double fx[], double fy[], double fz[],
        double m[], double& K, double dt, int n, double& sumvx, 
        double& sumvy, double& sumvz){
    
    double dtSqr=dt*dt;
    double dt2=2*dt;
    double xNEW, yNEW, zNEW, vxi, vyi, vzi;
    K=0; sumvx=0; sumvy=0; sumvz=0;
    
    for(int i=0; i<n; i++){
        xNEW=2.0*x[i]-xOLD[i]+dtSqr*fx[i]/m[i];
        yNEW=2.0*y[i]-yOLD[i]+dtSqr*fy[i]/m[i];
        zNEW=2.0*z[i]-zOLD[i]+dtSqr*fz[i]/m[i];
        
        vxi=(xNEW-xOLD[i])/dt2;
        vyi=(yNEW-yOLD[i])/dt2;
        vzi=(zNEW-zOLD[i])/dt2;
                       
        K=K+vxi*vxi+vyi*vyi+vzi*vzi;
        
        sumvx=sumvx+vxi;//velocity of the center of mass
        sumvy=sumvy+vyi;
        sumvz=sumvz+vzi;
        
        xOLD[i]=x[i]; yOLD[i]=y[i]; zOLD[i]=z[i];
        x[i]=xNEW; y[i]=yNEW; z[i]=zNEW;
    }

    K=0.5*m[0]*K; //Just assuming one mass for now
}

int main(int argc, char** argv) {
    //file io
    ofstream data;
    data.open("data.txt");
    //Number of particles
    int n=100;
    //Time information
    int tau=pow(10,6); //Number of time steps
    double dT=pow(10,-8); //Length of time step
    double T=tau*dT; //Total time
    //Particle info
    double mass=22.99;
    //Storage
    double x[n],y[n],z[n],xOLD[n],yOLD[n],zOLD[n],m[n],
            fx[n], fy[n], fz[n];
    //Simulation box length
    double l=10;
    //Kinetic/Potential/Total Energy;
    double K,V; double E;
    //Temperature
    double temp=300;
    //Boltzmann Const;
    double kB=1/temp;
    //momentum
    double sumvx, sumvy, sumvz;
    
    init(x, y, z, xOLD, yOLD, zOLD, m, mass, l, dT, temp, n); 
    
    for(int i=0; i<tau; i++){
        forces(x, y, z, fx, fy, fz, V, l, n);
        
//        if(i%1==0 & i>3659 & i<3665){
//            for(int f=0; f<n/2; f++){
//                if(f==41){
//                    data << "Loop#: " << i << endl;
//                    data << "fx" << f << "(verlet): " << fx[f] << endl;
//                }
//            }
//        }
        
        verlet(x, y, z, xOLD, yOLD, zOLD, fx, fy, fz, m, K, dT, n, 
               sumvx, sumvy, sumvz);
        
//        if(i%1==0 & i>3659 & i<3665){
//            for(int f=0; f<n/2; f++){
//                if(f==41){
//                    data << "x" << f << "(verlet): " << x[f] << endl;
//                }
//            }
//        }
        
        bCond(x, y, z, xOLD, yOLD, zOLD, l, n);
        
//        if(i%1==0 & i>3659 & i<3665){
//            for(int f=0; f<n/2; f++){
//                if(f==41){
//                    data << "x" << f << "(BC): " << x[f] << endl;
//                    data << "xOLD" << f << "(BC): " << xOLD[f] << endl << endl;
//                }
//            }
//        }
// 
        E=K+V;
        
        if(i%1000==0){            
            cout << "Loop# " << i << endl;
            cout << "V: " << V << endl;
            cout << "K: " << K << endl;
            cout << "E: " << E<< endl << endl;
        }
    }
    data.close();
}

