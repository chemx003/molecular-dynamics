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

void largeForce(double x[], double y[], double z[], double fx[], double fy[],
        double fz[], double E, int n){
    ofstream o;
    o.open("large_pos_debug.dat");
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

void forces(double x[], double y[], double z[], double fx[],
        double fy[], double fz[], double& V, double l, double& P, 
        double kB, double T, int n){
    
    double sigma=.00034, epsilon=1.656768, rc=2.25*sigma;
    double sigma2=sigma*sigma, rc2=rc*rc;
    double ecut;//discrepancy in sources
    double dx, dy, dz;
    double fr2, fr6, frp;
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
            
            //cout << sqrt(r2) << endl;
            
            if(r2<rc2){
                fr2=sigma2/r2;
                fr6=pow(fr2,3);
                frp=48*epsilon*fr6*(fr6-0.5)/r2;
                
 
                fxi=frp*dx; //forces between the pairs
                fyi=frp*dy;
                fzi=frp*dz;
                
                fx[i]=fx[i]+fxi; fx[j]=fx[j]-fxi; //total force on particle
                fy[i]=fy[i]+fyi; fy[j]=fy[j]-fyi;
                fz[i]=fz[i]+fzi; fz[j]=fz[j]-fzi;
                
                V=V+4.0*epsilon*fr6*(fr6-1.0);
                P=P+fxi*dx+fyi*dy+fzi*dz;
            }
        }
    }
    
    P=n*kB*pow(10,-21)*T+P*pow(10,-21)/3; P=P/(l*l*l*pow(10,-18));
    //converted to Pa
    
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

void verlet(double x[], double y[], double z[], double xOLD[],
        double yOLD[], double zOLD[], double fx[], double fy[], double fz[],
        double m[], double& K, double dt, int n, double& sumvx, 
        double& sumvy, double& sumvz, double l, int loop){
    
    double dtSqr=dt*dt;
    double dt2=2*dt;
    double xNEW, yNEW, zNEW, vxi, vyi, vzi;
    K=0; sumvx=0; sumvy=0; sumvz=0;
    
    for(int i=0; i<n; i++){
        if(x[i]-xOLD[i]>l-0.00001){
            x[i]=x[i]-l;
        }//kind of works...
        
        xNEW=2.0*x[i]-xOLD[i]+dtSqr*fx[i]/m[i];
        yNEW=2.0*y[i]-yOLD[i]+dtSqr*fy[i]/m[i];
        zNEW=2.0*z[i]-zOLD[i]+dtSqr*fz[i]/m[i];
        
        vxi=(xNEW-xOLD[i])/dt2;
        vyi=(yNEW-yOLD[i])/dt2;
        vzi=(zNEW-zOLD[i])/dt2;
        
        //cout << vxi << endl;
                       
        K=K+vxi*vxi+vyi*vyi+vzi*vzi;
        
        if(vxi>2000000){
            cout << i<< ": " << vxi <<"\t"<< loop << endl;
            cout << xNEW << "\t" << xOLD[i] <<"\t"<<x[i]<< endl << endl;
        }
        
        sumvx=sumvx+vxi;//velocity of the center of mass
        sumvy=sumvy+vyi;
        sumvz=sumvz+vzi;
        
        xOLD[i]=x[i]; yOLD[i]=y[i]; zOLD[i]=z[i];
        x[i]=xNEW; y[i]=yNEW; z[i]=zNEW;
	if(i==0){
	//cout<< "x: " << x[0]<< " vx: " << vxi << " fx: " << fx[0] <<endl<<endl;
    }}

    K=0.5*m[0]*K; //Just assuming one mass for now
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
    int tau=2000; //Number of time steps
    double dT=pow(10,-11); //Length of time step ** used a smaller step
    double T=tau*dT; //Total time
    //Particle info
    double mass=6.6335209*pow(10,-11);
    //Storage
    double x[n],y[n],z[n],xOLD[n],yOLD[n],zOLD[n],m[n],
            fx[n], fy[n], fz[n];
    //Simulation box length
    double l=10.229*0.00034;
    //Kinetic/Potential/Total Energy;
    double K,V; double E;
    //Temperature
    double temp=90;
    //Boltzmann Const;
    double kB=0.0138064852;
    //momentum
    double sumvx, sumvy, sumvz;
    //pressure
    double P;
    
    init(x, y, z, xOLD, yOLD, zOLD, m, mass, l, dT, temp, n); 
    writeXYZ(x,y,z,n);
//    for(int f=0; f<n; f++){
//        cout << "x" << f<< ": " << x[f] << " y: " << y[f] << " z: " << z[f] <<endl; 
//    }
    
    for(int i=0; i<tau; i++){
        forces(x, y, z, fx, fy, fz, V, l, P, kB, T, n);
        
        verlet(x, y, z, xOLD, yOLD, zOLD, fx, fy, fz, m, K, dT, n, 
               sumvx, sumvy, sumvz, l, i);
        
        bCond(x, y, z, xOLD, yOLD, zOLD, l, n);
        
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
            
            if(E>0){
                largeForce(x,y,z,fx,fy,fz,E,n);
            }
        }
    }
    pairCor(x,y,z,n,l);
}
