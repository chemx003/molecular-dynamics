#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <time.h>

using namespace std;

//while performing the simulation the pcf does not have a wall at 1... which
//is expected since sigma s=1.0... i'm going to plot the forces.

int main(int argc, char** argv){
    double mu=2.0, nu=1.0;
    double dx, dy, dz;
    double sigmaE=3.0, sigmaS=1.0, epsilonE=5.0, epsilonS=1.0;
    double kappa=sigmaE/sigmaS, kappaPrime=epsilonS/epsilonE;
    double chi=(pow(kappa,2.0)-1.0)/(pow(kappa,2.0)+1.0);
    double chiPrime=(1.0-pow(kappaPrime,1.0/mu))/(1.0+pow(kappaPrime,1.0/mu));
    double rc=20.0*sigmaS, rc2=rc*rc; //cuttoff
    double dot1, dot2, dot12, dot122, dotSum, dotSum2, dotDif, dotDif2;
    double g, gPrime, gHalf, dgx, dgy, dgz, dgxPrime, dgyPrime, dgzPrime;
    double R, R_1, R_2, R_6, R_7, R_12, R_13, distF;
    double ePrime, ePn, gPm, gPm1;
    double fxi, fyi, fzi;
    double dotSByChi,dotDByChi,dotSByChip, dotDByChip;
    double x[2],y[2],z[2];
    double ex[2],ey[2],ez[2];
    double V, distance;
    double sumx=0;

    ofstream o,f,rp;
    o.open("p2.data",ios::app);
    f.open("f.data",ios::app);
    rp.open("rp.data",ios::app);

    for(int k=1; k<2000; k++){
        distance=k/100.0;
        x[0]=distance; ex[0]=1;
        y[0]=0; ey[0]=0;
        z[0]=0; ez[0]=0;
        x[1]=0; ex[1]=1;
        y[1]=0; ey[1]=0;
        z[1]=0; ez[1]=0;

        dx=x[0]-x[1]; //components of distance vector
        dy=y[0]-y[1];
        dz=z[0]-z[1];

        double r2=dx*dx+dy*dy+dz*dz;
        double r=pow(r2,0.5);

        dot1=dx*ex[0]+dy*ey[0]+dz*ez[0];
        dot2=dx*ex[1]+dy*ey[1]+dz*ez[1];
        dot12=ex[0]*ex[1]+ey[0]*ey[1]+ez[0]*ez[1]; dot122=dot12*dot12;
        dotSum=dot1+dot2; dotSum2=dotSum*dotSum;
        dotDif=dot1-dot2; dotDif2=dotDif*dotDif;

        dotSByChi=dotSum/(1.0+chi*dot12);
        dotDByChi=dotDif/(1.0-chi*dot12);
        dotSByChip=dotSum/(1.0+chiPrime*dot12);
        dotDByChip=dotDif/(1.0-chiPrime*dot12);

        g=1.0-(chi/(2.0*r2))*(dotSum*dotSByChi+dotDif*dotDByChi);
        gPrime=1.0-(chiPrime/(2.0*r2))*(dotSum*dotSByChip+dotDif*dotDByChip); //epsilon
        ePrime=1.0/pow(1.0-chi*chi*dot122,0.5);

        gPm=pow(gPrime,mu);
        gPm1=pow(gPrime,mu-1);
        ePn=pow(ePrime,nu);

        distF=sigmaS/pow(g,0.5);

        if(r>2.5){
        R=(r-distF+sigmaS)/sigmaS;
        R_1=1.0/R;
        R_6=R_1*R_1*R_1*R_1*R_1*R_1;
        R_7=R_6*R_1;
        R_12=R_6*R_6;
        R_13=R_12*R_1;

        dgx=-(chi/r2)*(dotSByChi*(ex[0]+ex[1])+dotDByChi*(ex[0]-ex[1]))
                +(dx*chi/(r2*r2))*(dotSum*dotSByChi+dotDif*dotDByChi);
        dgy=-(chi/r2)*(dotSByChi*(ey[0]+ey[1])+dotDByChi*(ey[0]-ey[1]))
                +(dy*chi/(r2*r2))*(dotSum*dotSByChi+dotDif*dotDByChi);
        dgz=-(chi/r2)*(dotSByChi*(ez[0]+ez[1])+dotDByChi*(ez[0]-ez[1]))
                +(dz*chi/(r2*r2))*(dotSum*dotSByChi+dotDif*dotDByChi);

        dgxPrime=-(chiPrime/r2)*(dotSByChip*(ex[0]+ex[1])+dotDByChip*(ex[0]-ex[1]))
                +(dx*chiPrime/(r2*r2))*(dotSum*dotSByChip+dotDif*dotDByChip);
        dgyPrime=-(chiPrime/r2)*(dotSByChip*(ey[0]+ey[1])+dotDByChip*(ey[0]-ey[1]))
                +(dy*chiPrime/(r2*r2))*(dotSum*dotSByChip+dotDif*dotDByChip);
        dgzPrime=-(chiPrime/r2)*(dotSByChip*(ez[0]+ez[1])+dotDByChip*(ez[0]-ez[1]))
                +(dz*chiPrime/(r2*r2))*(dotSum*dotSByChip+dotDif*dotDByChip);

        fxi=-epsilonS*(ePn*gPm*(6*R_7-12*R_13)*(dx/r+(sigmaS/2)
            /(pow(g,1.5))*dgx)+mu*gPm1*(R_12-R_6)*dgxPrime); //forces between the pairs
        fyi=-epsilonS*(ePn*gPm*(6*R_7-12*R_13)*(dy/r+(sigmaS/2)
            /(pow(g,1.5))*dgy)+mu*gPm1*(R_12-R_6)*dgyPrime);
        fzi=-epsilonS*(ePn*gPm*(6*R_7-12*R_13)*(dz/r+(sigmaS/2)
            /(pow(g,1.5))*dgz)+mu*gPm1*(R_12-R_6)*dgzPrime);}
        else{
            R=r-3;
            R_1=1.0/R;
            R_6=R_1*R_1*R_1*R_1*R_1*R_1;
            R_7=R_6*R_1;
            R_12=R_6*R_6;
            R_13=R_12*R_1;

            dgx=-(chi/r2)*(dotSByChi*(ex[0]+ex[1])+dotDByChi*(ex[0]-ex[1]))
                    +(dx*chi/(r2*r2))*(dotSum*dotSByChi+dotDif*dotDByChi);
            dgy=-(chi/r2)*(dotSByChi*(ey[0]+ey[1])+dotDByChi*(ey[0]-ey[1]))
                    +(dy*chi/(r2*r2))*(dotSum*dotSByChi+dotDif*dotDByChi);
            dgz=-(chi/r2)*(dotSByChi*(ez[0]+ez[1])+dotDByChi*(ez[0]-ez[1]))
                    +(dz*chi/(r2*r2))*(dotSum*dotSByChi+dotDif*dotDByChi);

            dgxPrime=-(chiPrime/r2)*(dotSByChip*(ex[0]+ex[1])+dotDByChip*(ex[0]-ex[1]))
                    +(dx*chiPrime/(r2*r2))*(dotSum*dotSByChip+dotDif*dotDByChip);
            dgyPrime=-(chiPrime/r2)*(dotSByChip*(ey[0]+ey[1])+dotDByChip*(ey[0]-ey[1]))
                    +(dy*chiPrime/(r2*r2))*(dotSum*dotSByChip+dotDif*dotDByChip);
            dgzPrime=-(chiPrime/r2)*(dotSByChip*(ez[0]+ez[1])+dotDByChip*(ez[0]-ez[1]))
                    +(dz*chiPrime/(r2*r2))*(dotSum*dotSByChip+dotDif*dotDByChip);

            fxi=-epsilonS*(ePn*gPm*(6*R_7-12*R_13)*(dx/r+(sigmaS/2)
                /(pow(g,1.5))*dgx)+mu*gPm1*(R_12-R_6)*dgxPrime); //forces between the pairs
            fyi=-epsilonS*(ePn*gPm*(6*R_7-12*R_13)*(dy/r+(sigmaS/2)
                /(pow(g,1.5))*dgy)+mu*gPm1*(R_12-R_6)*dgyPrime);
            fzi=-epsilonS*(ePn*gPm*(6*R_7-12*R_13)*(dz/r+(sigmaS/2)
                /(pow(g,1.5))*dgz)+mu*gPm1*(R_12-R_6)*dgzPrime);
        }

        V=4.0*epsilonS*ePn*gPm*(R_12-R_6);

        double twelveSix=R_12-R_6;

        cout<<twelveSix<<endl;

        o << distance << "\t" << V << "\n";
        f << distance << "\t" << fxi << "\n";
        rp << distance << "\t" << twelveSix << "\n";
        sumx=sumx+fxi;
    }

    cout<<sumx<<endl;
    sumx=0;

        for(int k=1; k<2000; k++){
        distance=k/100.0;
        x[0]=0; ex[0]=1;
        y[0]=distance; ey[0]=0;
        z[0]=0; ez[0]=0;
        x[1]=0; ex[1]=1;
        y[1]=0; ey[1]=0;
        z[1]=0; ez[1]=0;

        dx=x[0]-x[1]; //components of distance vector
        dy=y[0]-y[1];
        dz=z[0]-z[1];

        double r2=dx*dx+dy*dy+dz*dz;
        double r=pow(r2,0.5);

        dot1=dx*ex[0]+dy*ey[0]+dz*ez[0];
        dot2=dx*ex[1]+dy*ey[1]+dz*ez[1];
        dot12=ex[0]*ex[1]+ey[0]*ey[1]+ez[0]*ez[1]; dot122=dot12*dot12;
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

        distF=sigmaS/pow(g,0.5);

        if(r2>1.0){
        R=(r-distF+sigmaS)/sigmaS;
        R_1=1/R;
        R_6=R_1*R_1*R_1*R_1*R_1*R_1;
        R_7=R_6*R_1;
        R_12=R_6*R_6;
        R_13=R_12*R_1;

        dgx=-(chi/r2)*(dotSByChi*(ex[0]+ex[1])+dotDByChi*(ex[0]-ex[1]))
                +(dx*chi/(r2*r2))*(dotSum*dotSByChi+dotDif*dotDByChi);
        dgy=-(chi/r2)*(dotSByChi*(ey[0]+ey[1])+dotDByChi*(ey[0]-ey[1]))
                +(dy*chi/(r2*r2))*(dotSum*dotSByChi+dotDif*dotDByChi);
        dgz=-(chi/r2)*(dotSByChi*(ez[0]+ez[1])+dotDByChi*(ez[0]-ez[1]))
                +(dz*chi/(r2*r2))*(dotSum*dotSByChi+dotDif*dotDByChi);

        dgxPrime=-(chiPrime/r2)*(dotSByChip*(ex[0]+ex[1])+dotDByChip*(ex[0]-ex[1]))
                +(dx*chiPrime/(r2*r2))*(dotSum*dotSByChip+dotDif*dotDByChip);
        dgyPrime=-(chiPrime/r2)*(dotSByChip*(ey[0]+ey[1])+dotDByChip*(ey[0]-ey[1]))
                +(dy*chiPrime/(r2*r2))*(dotSum*dotSByChip+dotDif*dotDByChip);
        dgzPrime=-(chiPrime/r2)*(dotSByChip*(ez[0]+ez[1])+dotDByChip*(ez[0]-ez[1]))
                +(dz*chiPrime/(r2*r2))*(dotSum*dotSByChip+dotDif*dotDByChip);

        fxi=-epsilonS*(ePn*gPm*(6*R_7-12*R_13)*(dx/r+(sigmaS/2)
            /(pow(g,1.5))*dgx)+mu*gPm1*(R_12-R_6)*dgxPrime); //forces between the pairs
        fyi=-epsilonS*(ePn*gPm*(6*R_7-12*R_13)*(dy/r+(sigmaS/2)
            /(pow(g,1.5))*dgy)+mu*gPm1*(R_12-R_6)*dgyPrime);
        fzi=-epsilonS*(ePn*gPm*(6*R_7-12*R_13)*(dz/r+(sigmaS/2)
            /(pow(g,1.5))*dgz)+mu*gPm1*(R_12-R_6)*dgzPrime);}

        else{
            R=r;
            R_1=1.0/R;
            R_6=R_1*R_1*R_1*R_1*R_1*R_1;
            R_7=R_6*R_1;
            R_12=R_6*R_6;
            R_13=R_12*R_1;

            dgx=-(chi/r2)*(dotSByChi*(ex[0]+ex[1])+dotDByChi*(ex[0]-ex[1]))
                    +(dx*chi/(r2*r2))*(dotSum*dotSByChi+dotDif*dotDByChi);
            dgy=-(chi/r2)*(dotSByChi*(ey[0]+ey[1])+dotDByChi*(ey[0]-ey[1]))
                    +(dy*chi/(r2*r2))*(dotSum*dotSByChi+dotDif*dotDByChi);
            dgz=-(chi/r2)*(dotSByChi*(ez[0]+ez[1])+dotDByChi*(ez[0]-ez[1]))
                    +(dz*chi/(r2*r2))*(dotSum*dotSByChi+dotDif*dotDByChi);

            dgxPrime=-(chiPrime/r2)*(dotSByChip*(ex[0]+ex[1])+dotDByChip*(ex[0]-ex[1]))
                    +(dx*chiPrime/(r2*r2))*(dotSum*dotSByChip+dotDif*dotDByChip);
            dgyPrime=-(chiPrime/r2)*(dotSByChip*(ey[0]+ey[1])+dotDByChip*(ey[0]-ey[1]))
                    +(dy*chiPrime/(r2*r2))*(dotSum*dotSByChip+dotDif*dotDByChip);
            dgzPrime=-(chiPrime/r2)*(dotSByChip*(ez[0]+ez[1])+dotDByChip*(ez[0]-ez[1]))
                    +(dz*chiPrime/(r2*r2))*(dotSum*dotSByChip+dotDif*dotDByChip);

            fxi=-epsilonS*(ePn*gPm*(6*R_7-12*R_13)*(dx/r+(sigmaS/2)
                /(pow(g,1.5))*dgx)+mu*gPm1*(R_12-R_6)*dgxPrime); //forces between the pairs
            fyi=-epsilonS*(ePn*gPm*(6*R_7-12*R_13)*(dy/r+(sigmaS/2)
                /(pow(g,1.5))*dgy)+mu*gPm1*(R_12-R_6)*dgyPrime);
            fzi=-epsilonS*(ePn*gPm*(6*R_7-12*R_13)*(dz/r+(sigmaS/2)
                /(pow(g,1.5))*dgz)+mu*gPm1*(R_12-R_6)*dgzPrime);
        }

        V=4.0*epsilonS*ePn*gPm*(R_12-R_6);

        double twelveSix=R_12-R_6;

        o << distance << "\t" << V << "\n";
        f << distance << "\t" << fyi << "\n";
        rp << distance << "\t" << twelveSix << "\n";
        sumx=sumx+fxi;
    }

    cout<<sumx<<endl;
    sumx=0;

        for(int k=1; k<2000; k++){
        distance=k/100.0;
        x[0]=distance; ex[0]=0;
        y[0]=0; ey[0]=0;
        z[0]=0; ez[0]=0;
        x[1]=0; ex[1]=1;
        y[1]=0; ey[1]=0;
        z[1]=0; ez[1]=0;

        dx=x[0]-x[1]; //components of distance vector
        dy=y[0]-y[1];
        dz=z[0]-z[1];

        double r2=dx*dx+dy*dy+dz*dz;
        double r=pow(r2,0.5);

        dot1=dx*ex[0]+dy*ey[0]+dz*ez[0];
        dot2=dx*ex[1]+dy*ey[1]+dz*ez[1];
        dot12=ex[0]*ex[1]+ey[0]*ey[1]+ez[0]*ez[1]; dot122=dot12*dot12;
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

        distF=sigmaS/pow(g,0.5);

        if(r2>1.0){
        R=(r-distF+sigmaS)/sigmaS;
        R_1=1/R;
        R_6=R_1*R_1*R_1*R_1*R_1*R_1;
        R_7=R_6*R_1;
        R_12=R_6*R_6;
        R_13=R_12*R_1;

        dgx=-(chi/r2)*(dotSByChi*(ex[0]+ex[1])+dotDByChi*(ex[0]-ex[1]))
                +(dx*chi/(r2*r2))*(dotSum*dotSByChi+dotDif*dotDByChi);
        dgy=-(chi/r2)*(dotSByChi*(ey[0]+ey[1])+dotDByChi*(ey[0]-ey[1]))
                +(dy*chi/(r2*r2))*(dotSum*dotSByChi+dotDif*dotDByChi);
        dgz=-(chi/r2)*(dotSByChi*(ez[0]+ez[1])+dotDByChi*(ez[0]-ez[1]))
                +(dz*chi/(r2*r2))*(dotSum*dotSByChi+dotDif*dotDByChi);

        dgxPrime=-(chiPrime/r2)*(dotSByChip*(ex[0]+ex[1])+dotDByChip*(ex[0]-ex[1]))
                +(dx*chiPrime/(r2*r2))*(dotSum*dotSByChip+dotDif*dotDByChip);
        dgyPrime=-(chiPrime/r2)*(dotSByChip*(ey[0]+ey[1])+dotDByChip*(ey[0]-ey[1]))
                +(dy*chiPrime/(r2*r2))*(dotSum*dotSByChip+dotDif*dotDByChip);
        dgzPrime=-(chiPrime/r2)*(dotSByChip*(ez[0]+ez[1])+dotDByChip*(ez[0]-ez[1]))
                +(dz*chiPrime/(r2*r2))*(dotSum*dotSByChip+dotDif*dotDByChip);

        fxi=-epsilonS*(ePn*gPm*(6*R_7-12*R_13)*(dx/r+(sigmaS/2)
            /(pow(g,1.5))*dgx)+mu*gPm1*(R_12-R_6)*dgxPrime); //forces between the pairs
        fyi=-epsilonS*(ePn*gPm*(6*R_7-12*R_13)*(dy/r+(sigmaS/2)
            /(pow(g,1.5))*dgy)+mu*gPm1*(R_12-R_6)*dgyPrime);
        fzi=-epsilonS*(ePn*gPm*(6*R_7-12*R_13)*(dz/r+(sigmaS/2)
            /(pow(g,1.5))*dgz)+mu*gPm1*(R_12-R_6)*dgzPrime);}

        else{
            R=r;
            R_1=1.0/R;
            R_6=R_1*R_1*R_1*R_1*R_1*R_1;
            R_7=R_6*R_1;
            R_12=R_6*R_6;
            R_13=R_12*R_1;

            dgx=-(chi/r2)*(dotSByChi*(ex[0]+ex[1])+dotDByChi*(ex[0]-ex[1]))
                    +(dx*chi/(r2*r2))*(dotSum*dotSByChi+dotDif*dotDByChi);
            dgy=-(chi/r2)*(dotSByChi*(ey[0]+ey[1])+dotDByChi*(ey[0]-ey[1]))
                    +(dy*chi/(r2*r2))*(dotSum*dotSByChi+dotDif*dotDByChi);
            dgz=-(chi/r2)*(dotSByChi*(ez[0]+ez[1])+dotDByChi*(ez[0]-ez[1]))
                    +(dz*chi/(r2*r2))*(dotSum*dotSByChi+dotDif*dotDByChi);

            dgxPrime=-(chiPrime/r2)*(dotSByChip*(ex[0]+ex[1])+dotDByChip*(ex[0]-ex[1]))
                    +(dx*chiPrime/(r2*r2))*(dotSum*dotSByChip+dotDif*dotDByChip);
            dgyPrime=-(chiPrime/r2)*(dotSByChip*(ey[0]+ey[1])+dotDByChip*(ey[0]-ey[1]))
                    +(dy*chiPrime/(r2*r2))*(dotSum*dotSByChip+dotDif*dotDByChip);
            dgzPrime=-(chiPrime/r2)*(dotSByChip*(ez[0]+ez[1])+dotDByChip*(ez[0]-ez[1]))
                    +(dz*chiPrime/(r2*r2))*(dotSum*dotSByChip+dotDif*dotDByChip);

            fxi=-epsilonS*(ePn*gPm*(6*R_7-12*R_13)*(dx/r+(sigmaS/2)
                /(pow(g,1.5))*dgx)+mu*gPm1*(R_12-R_6)*dgxPrime); //forces between the pairs
            fyi=-epsilonS*(ePn*gPm*(6*R_7-12*R_13)*(dy/r+(sigmaS/2)
                /(pow(g,1.5))*dgy)+mu*gPm1*(R_12-R_6)*dgyPrime);
            fzi=-epsilonS*(ePn*gPm*(6*R_7-12*R_13)*(dz/r+(sigmaS/2)
                /(pow(g,1.5))*dgz)+mu*gPm1*(R_12-R_6)*dgzPrime);
        }

        V=4.0*epsilonS*ePn*gPm*(R_12-R_6);

        double twelveSix=R_12-R_6;

        o << distance << "\t" << V << "\n";
        f << distance << "\t" << fxi << "\n";
        rp << distance << "\t" << twelveSix << "\n";
        sumx=sumx+fxi;

    }

    cout<<sumx<<endl;
    sumx=0;

        for(int k=1; k<2000; k++){
        distance=k/100.0;
        x[0]=0; ex[0]=0;
        y[0]=distance; ey[0]=0;
        z[0]=0; ez[0]=0;
        x[1]=0; ex[1]=1;
        y[1]=0; ey[1]=0;
        z[1]=0; ez[1]=0;

        dx=x[0]-x[1]; //components of distance vector
        dy=y[0]-y[1];
        dz=z[0]-z[1];

        double r2=dx*dx+dy*dy+dz*dz;
        double r=pow(r2,0.5);

        dot1=dx*ex[0]+dy*ey[0]+dz*ez[0];
        dot2=dx*ex[1]+dy*ey[1]+dz*ez[1];
        dot12=ex[0]*ex[1]+ey[0]*ey[1]+ez[0]*ez[1]; dot122=dot12*dot12;
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

        distF=sigmaS/pow(g,0.5);

        if(r2>1.0){
        R=(r-distF+sigmaS)/sigmaS;
        R_1=1/R;
        R_6=R_1*R_1*R_1*R_1*R_1*R_1;
        R_7=R_6*R_1;
        R_12=R_6*R_6;
        R_13=R_12*R_1;

        dgx=-(chi/r2)*(dotSByChi*(ex[0]+ex[1])+dotDByChi*(ex[0]-ex[1]))
                +(dx*chi/(r2*r2))*(dotSum*dotSByChi+dotDif*dotDByChi);
        dgy=-(chi/r2)*(dotSByChi*(ey[0]+ey[1])+dotDByChi*(ey[0]-ey[1]))
                +(dy*chi/(r2*r2))*(dotSum*dotSByChi+dotDif*dotDByChi);
        dgz=-(chi/r2)*(dotSByChi*(ez[0]+ez[1])+dotDByChi*(ez[0]-ez[1]))
                +(dz*chi/(r2*r2))*(dotSum*dotSByChi+dotDif*dotDByChi);

        dgxPrime=-(chiPrime/r2)*(dotSByChip*(ex[0]+ex[1])+dotDByChip*(ex[0]-ex[1]))
                +(dx*chiPrime/(r2*r2))*(dotSum*dotSByChip+dotDif*dotDByChip);
        dgyPrime=-(chiPrime/r2)*(dotSByChip*(ey[0]+ey[1])+dotDByChip*(ey[0]-ey[1]))
                +(dy*chiPrime/(r2*r2))*(dotSum*dotSByChip+dotDif*dotDByChip);
        dgzPrime=-(chiPrime/r2)*(dotSByChip*(ez[0]+ez[1])+dotDByChip*(ez[0]-ez[1]))
                +(dz*chiPrime/(r2*r2))*(dotSum*dotSByChip+dotDif*dotDByChip);

        fxi=-epsilonS*(ePn*gPm*(6*R_7-12*R_13)*(dx/r+(sigmaS/2)
            /(pow(g,1.5))*dgx)+mu*gPm1*(R_12-R_6)*dgxPrime); //forces between the pairs
        fyi=-epsilonS*(ePn*gPm*(6*R_7-12*R_13)*(dy/r+(sigmaS/2)
            /(pow(g,1.5))*dgy)+mu*gPm1*(R_12-R_6)*dgyPrime);
        fzi=-epsilonS*(ePn*gPm*(6*R_7-12*R_13)*(dz/r+(sigmaS/2)
            /(pow(g,1.5))*dgz)+mu*gPm1*(R_12-R_6)*dgzPrime);}

        else{
            R=r;
            R_1=1.0/R;
            R_6=R_1*R_1*R_1*R_1*R_1*R_1;
            R_7=R_6*R_1;
            R_12=R_6*R_6;
            R_13=R_12*R_1;

            dgx=-(chi/r2)*(dotSByChi*(ex[0]+ex[1])+dotDByChi*(ex[0]-ex[1]))
                    +(dx*chi/(r2*r2))*(dotSum*dotSByChi+dotDif*dotDByChi);
            dgy=-(chi/r2)*(dotSByChi*(ey[0]+ey[1])+dotDByChi*(ey[0]-ey[1]))
                    +(dy*chi/(r2*r2))*(dotSum*dotSByChi+dotDif*dotDByChi);
            dgz=-(chi/r2)*(dotSByChi*(ez[0]+ez[1])+dotDByChi*(ez[0]-ez[1]))
                    +(dz*chi/(r2*r2))*(dotSum*dotSByChi+dotDif*dotDByChi);

            dgxPrime=-(chiPrime/r2)*(dotSByChip*(ex[0]+ex[1])+dotDByChip*(ex[0]-ex[1]))
                    +(dx*chiPrime/(r2*r2))*(dotSum*dotSByChip+dotDif*dotDByChip);
            dgyPrime=-(chiPrime/r2)*(dotSByChip*(ey[0]+ey[1])+dotDByChip*(ey[0]-ey[1]))
                    +(dy*chiPrime/(r2*r2))*(dotSum*dotSByChip+dotDif*dotDByChip);
            dgzPrime=-(chiPrime/r2)*(dotSByChip*(ez[0]+ez[1])+dotDByChip*(ez[0]-ez[1]))
                    +(dz*chiPrime/(r2*r2))*(dotSum*dotSByChip+dotDif*dotDByChip);

            fxi=-epsilonS*(ePn*gPm*(6*R_7-12*R_13)*(dx/r+(sigmaS/2)
                /(pow(g,1.5))*dgx)+mu*gPm1*(R_12-R_6)*dgxPrime); //forces between the pairs
            fyi=-epsilonS*(ePn*gPm*(6*R_7-12*R_13)*(dy/r+(sigmaS/2)
                /(pow(g,1.5))*dgy)+mu*gPm1*(R_12-R_6)*dgyPrime);
            fzi=-epsilonS*(ePn*gPm*(6*R_7-12*R_13)*(dz/r+(sigmaS/2)
                /(pow(g,1.5))*dgz)+mu*gPm1*(R_12-R_6)*dgzPrime);
        }

        V=4.0*epsilonS*ePn*gPm*(R_12-R_6);

        double twelveSix=R_12-R_6;

        o << distance << "\t" << V << "\n";
        f << distance << "\t" << fyi << "\n";
        rp << distance << "\t" << twelveSix << "\n";
        sumx=sumx+fxi;
    }

    cout<<sumx<<endl;
    sumx=0;
}
