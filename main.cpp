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
#include <time.h>

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
        
        ux[i]=ux[i]+dt*(gx[i]/I+lm*ex[i]);
        uy[i]=uy[i]+dt*(gy[i]/I+lm*ey[i]);
        uz[i]=uz[i]+dt*(gz[i]/I+lm*ez[i]);
        
        ex[i]=ex[i]+dt*ux[i];
        ey[i]=ey[i]+dt*uy[i];
        ez[i]=ez[i]+dt*uz[i];
        
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
        double& V, double l, double& P, double kB, double T, int n,
        double sigE, int loop){
    //1.1045188 0.00081
    double mu=2.0, nu=1.0;
    double dx, dy, dz;
    double sigmaE=sigE, sigmaS=1.0, epsilonE=0.2, epsilonS=1.0;
    double kappa=sigmaE/sigmaS, kappaPrime=epsilonS/epsilonE;
    double chi=(pow(kappa,2.0)-1.0)/(pow(kappa,2.0)+1.0);
    double chiPrime=(pow(kappaPrime,1.0/mu)-1.0)/(pow(kappaPrime,1.0/mu)+1.0);
    double rc=3.25*sigmaS, rc2=rc*rc; //cuttoff
    double dot1, dot2, dot12, dot122, dotSum, dotSum2, dotDif, dotDif2;
    double g, gPrime, gHalf, dgx, dgy, dgz, dgxPrime, dgyPrime, dgzPrime;
    double R, R_1, R_2, R_6, R_7, R_12, R_13, distF;
    double ePrime, ePn, gPm, gPm1;
    double fxi, fyi, fzi;
    double dotSByChi,dotDByChi,dotSByChip, dotDByChip;
    
    V=0;
    P=0;
    
    for(int i=0; i<n; i++){
        fx[i]=0;
        fy[i]=0;
        fz[i]=0;
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
            
            if(r2<rc2){
                dot1=dx*ex[i]+dy*ey[i]+dz*ez[i]; 
                dot2=dx*ex[j]+dy*ey[j]+dz*ez[j];
                dot12=ex[i]*ex[j]+ey[i]*ey[j]+ez[i]*ez[j]; dot122=dot12*dot12;
                dotSum=dot1+dot2; dotSum2=dotSum*dotSum;
                dotDif=dot1-dot2; dotDif2=dotDif*dotDif;
                
                dotSByChi=dotSum/(1+chi*dot12);
                dotDByChi=dotDif/(1-chi*dot12);
                dotSByChip=dotSum/(1+chiPrime*dot12);
                dotDByChip=dotDif/(1-chiPrime*dot12);
                
                g=1.0-(chi/(2*r2))*(dotSum*dotSByChi+dotDif*dotDByChi);
                gPrime=1.0-(chiPrime/(2*r2))*(dotSum*dotSByChip+dotDif*dotDByChip); //epsilon
                ePrime=1/pow(1-chi*chi*dot122,0.5);
                
                gPm=pow(gPrime,mu);
                gPm1=pow(gPrime,mu-1);
                ePn=pow(ePrime,nu);
                
                gHalf=pow(g,0.5);
                
                distF=sigmaS/pow(g,0.5);

                R=(r-distF+sigmaS)/sigmaS;
                R_1=1/R;
                R_6=R_1*R_1*R_1*R_1*R_1*R_1;
                R_7=R_6*R_1;
                R_12=R_6*R_6;
                R_13=R_12*R_1;

                ePrime=1/pow(1-chi*chi*dot122,0.5);

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
                
                fxi=-epsilonS*(ePn*gPm*(6*R_7-12*R_13)*(dx/r+(sigmaS/2)
                    /(pow(g,1.5))*dgx)+mu*gPm1*(R_12-R_6)*dgxPrime); //forces between the pairs
                fyi=-epsilonS*(ePn*gPm*(6*R_7-12*R_13)*(dy/r+(sigmaS/2)
                    /(pow(g,1.5))*dgy)+mu*gPm1*(R_12-R_6)*dgyPrime);
                fzi=-epsilonS*(ePn*gPm*(6*R_7-12*R_13)*(dz/r+(sigmaS/2)
                    /(pow(g,1.5))*dgz)+mu*gPm1*(R_12-R_6)*dgzPrime);
                
//                if(fxi>=10 || fyi>=10 || fzi>=10){
//                        cout<<"fx("<< loop <<"," <<i << ", " << j<< "): " <<fxi<<endl;
//                        cout<<"fy("<< loop <<"," <<i << ", " << j<< "): " <<fyi<<endl;
//                        cout<<"fz("<< loop <<"," <<i << ", " << j<< "): " <<fzi<<endl;
//                        cout<<"r: "<< r << " (" << dx << "," << dy << "," << dz << ")" << endl;
//                        cout<<"dot 1: "<< dot1 << " dot 2: " << dot2 << " dot12: " << dot12 << endl;
//                        cout<<"distF: "<< distF << " chi: " << chi << " chiPrime " << chiPrime << endl;
//                        cout<<"g: " << g << " gPrime: " << gPrime << " ePrime: " << ePn << endl;
//                        cout<<"gHalf: " << gHalf << endl;
//                        cout<<"dg: ("<< dgx <<"," << dgy << ", " << dgz << ")" << endl; 
//                        cout<<"dgPrime: (" << dgxPrime << "," << dgyPrime << "," << dgzPrime << ")" << endl;
//                        cout<<"R: " << R << " R_6: " << R_6 << endl<<endl;
//                        lrg++;
//                }

                fx[i]=fx[i]+fxi; fx[j]=fx[j]-fxi; //total force on particle
                fy[i]=fy[i]+fyi; fy[j]=fy[j]-fyi;
                fz[i]=fz[i]+fzi; fz[j]=fz[j]-fzi;
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
        double vy[], double vz[], double ex[], double ey[], double ez[],
        double m[], double mass, double l, double dt, double temp, int n){
    
    double sumvx=0.0, sumvy=0.0, sumvz=0.0; //used to set lin mtm = 0
    double sumx=0.0, sumy=0.0, sumz=0.0; // for debugging get rid of later
    double sumv2x=0.0, sumv2y=0.0, sumv2z=0.0; //set kinetic energy
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
                    x[p]=(i+0.5+dRand(-0.25,0.25))*a;
                    y[p]=(j+0.5+dRand(-0.25,0.25))*a;
                    z[p]=(k+0.5+dRand(-0.25,0.25))*a;
                
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
        if(histo[i][1]<100 && histo[i][1]>0)
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
        
        //cout << vxi << endl;
        
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
    o.open("random-testing.xyz",ios::app); //I should make this a setting in main())
    double t,j,k;
    
    o << 255 << endl;
    
    for(int i=0; i<n; i++){

        t=x[i]*100;
        j=y[i]*100;
        k=z[i]*100;
        o << "Ar" << "\t" <<  t << "\t" << j << "\t" << k << "\n";
        
    }
    o.close();
}



int main(int argc, char** argv) {
    //Number of particles
    int n=256;
    //Time information
    int tau=500;//10*pow(10,3); //Number of time steps
    double dT=0.0015;//pow(10,-4); //Length of time step ** used a smaller step
    double T=tau*dT; //Total time
    //Particle info
    double mass=1;
    //Storage
    double x[n],y[n],z[n],vx[n],vy[n],vz[n],ex[n],ey[n],ez[n], ux[n], uy[n], 
            uz[n],m[n],fx[n],fy[n],fz[n],gx[n],gy[n],gz[n];
    //Simulation box length
    double l=10.857670466; //scaled density of 0.2
    //Kinetic/Potential/Total Energy;
    double K,V; double E;
    //Temperature
    double temp=3.0;
    //Boltzmann Cons     
    double kB=0.0025;
    //momentum
    double sumvx, sumvy, sumvz;
    //pressure
    double P;
    //moment of inertia
    double I;
    //SigmaE
    double sigE=1.60;

    int rand;//
    do {//
        //Random seed;
        srand(rand*time(NULL));//time(NULL)
        temp=3.0;//
        init(x, y, z, vx, vy, vz, ex, ey, ez, m, mass, l, dT, temp, n); 
        writeXYZ(x, y, z, n);
        gb(x, y, z, fx, fy, fz, ex, ey, ez, V, l, P, kB, T, n, sigE, 0);
        halfstep(x, y, z, vx, vy, vz,fx, fy, fz, mass, dT, n);
        bCond(x, y, z, l, n);
        writeXYZ(x, y, z, n);
        gb(x, y, z, fx, fy, fz, ex, ey, ez, V, l, P, kB, T, n, sigE, 1);//
        leapfrog(x, y, z, vx, vy, vz, fx, fy, fz, mass, K, dT, n, //
                        sumvx, sumvy, sumvz, l, 1);//
        bCond(x, y, z, l, n);//
        temp=2*K/(3*n*kB);//
        cout<<temp<<endl;//
        rand++;//
    } while(temp>100000);//


    for(int i=2; i<tau; i++){
        gb(x, y, z, fx, fy, fz, ex, ey, ez, V, l, P, kB, T, n, sigE, i);
        
        leapfrog(x, y, z, vx, vy, vz, fx, fy, fz, mass, K, dT, n, 
                    sumvx, sumvy, sumvz, l, i);
        
        bCond(x, y, z, l, n);
        
        writeXYZ(x,y,z,n);
        
        E=K+V; //in scaled units
        temp=2*K/(3*n*kB); //in kelvin kg*m^2/s^2 -15*-6^2/-3^2  /-21
        
        if(i%100==0 || i==2){            
            cout << "Loop# " << i << endl;
            cout << "V: " << V << endl;
            cout << "K: " << K << endl;
            cout << "E: " << E << endl;
            cout << "T: " << temp << endl;
            cout << "P: " << P << endl;
            double fxavg;
            double fxmax=0;
            for(int f=0; f<n; f++){
                fxavg=fxavg+fx[f];
                if(abs(fx[f])>abs(fxmax)){
                    fxmax=fx[f];
                }
            }
            
            fxavg=fxavg/n;
            cout<< "Maximum force: " << fxmax << endl;
            cout<< "Average x-force: " << fxavg << endl <<endl;
            // if(E>0){
            //      largeForce(x,y,z,fx,fy,fz,E,n);
            // }
         }
     }
     pairCor(x,y,z,n,l);
}

