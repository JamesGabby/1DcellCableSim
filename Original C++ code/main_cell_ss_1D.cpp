/*
 *  O'Hara-Rudy Human model with RK method
 *  Created by Mao-Tsuen Jeng on 6/13/11.
 *
 *	Copyright 2011 __Colleen Clancy Lab__. All rights reserved.
 *
 *  One-D simulations
 */




#include <fstream>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <ctime>
#pragma warning(disable:4996)
//#include <omp.h>
#include <cstdlib>
using namespace std;




const double BCL = 1000;

const int stimuli = 2;

const int tl = 165;  // length of the tissue
const int tw = 1;  // width of the tissue

const double base_dt = 0.005;
const int fold = 1.0/base_dt;

const double dx= 0.01; //cm


const double gj_use = 2.5;// uS
const double rad = 0.0011;
const double Rmyo = 150.0;       //ohm*cm
const double Rg = 3.14159*pow(rad,2)/(gj_use*1.0e-6);
const double Df = (1e3*rad)/(4.0*(Rmyo+Rg/dx));


typedef struct cell {
	double y0[41];
	double highestV, lowestV, Apd90;
    double highestdvdt, V90;
    double DI, t1, t2;
	double  highestV2, lowestV2, t3, t4, s2APD;
	double dvdt, v_new, v, dvdt2;
	double currents[20];
    int celltype;
    double I_inj;
    double coef_PKA;
    
} Cell;

#include "ord_model_betaStim.h"
#include "integrate_rk2rk1.h"

typedef struct simState {
	double t, tt;
	int tstep, counter, beat;
	Cell cellData[tw][tl];
} SimState;

SimState theState;
SimState *S = &theState;

double Calcu_I_Total(Cell *theCell, double dt );


int main () {
	
	double p = 0;//par;  // Parameter array for passing nondefault conditions
	
	//// Initial conditions
	double v=-87.;
	double nai=7.;
	double nass=nai;
	double ki=145.;
	double kss=ki;
	double cai=1.0e-4;
	double cass=cai;
	double cansr=1.2;
	double cajsr=cansr;
	double m=0.;
	double hf=1.;
	double hs=1.;
	double j=1.;
	double hsp=1.;
	double jp=1.;
	double mL=0.;
	double hL=1.;
	double hLp=1.;
	double a=0.;
	double iF=1.;
	double iS=1.;
	double ap=0.;
	double iFp=1.;
	double iSp=1.;
	double d=0.;
	double ff=1.;
	double fs=1.;
	double fcaf=1.;
	double fcas=1.;
	double jca=1.;
	double nca=0.;
	double ffp=1.;
	double fcafp=1.;
	double xrf=0.;
	double xrs=0.;
	double xs1=0.;
	double xs2=0.;
	double xk1=1.;
	double Jrelnp=0.;
	double Jrelp=0.;
	double CaMKt=0.;
	// %X0 is the vector for initial sconditions for state variables
	//	X0=[v nai nass ki kss cai cass cansr cajsr m hf hs j hsp jp mL hL hLp a iF iS ap iFp iSp d ff fs fcaf fcas jca nca ffp fcafp xrf xrs xs1 xs2 xk1 Jrelnp Jrelp CaMKt]';
    
    int sy0 = 41;
    
    double y0[] = { v, nai, nass, ki, kss,
        cai, cass, cansr, cajsr, m,
        hf, hs, j, hsp, jp,
        mL, hL, hLp, a, iF,
        iS, ap, iFp, iSp, d,
        ff, fs, fcaf, fcas, jca,
        nca, ffp, fcafp, xrf, xrs,
        xs1, xs2, xk1, Jrelnp, Jrelp,
        CaMKt };
	
	
	char name[30];
	FILE *output;
	FILE *output_m;
	FILE *output1;
	FILE *output2;
	FILE *output3;
	int ll, ww;
	
	double V, dvdt, dvdt2;
	Cell *theCell;
	double dt;
	double I_Total;
	double dv;
	double Vm = -86;
	double cycle_length;
	
	double I_inj;
	int done = 0;
	const char *SSstateFileName = "model_INPUTSteadyStateControl.1DInitials";
    const char *stateFileName = "Save1D.currentState"; // For current simulated state output name

	double waitTime;
	
	time_t startTime;
	time_t previousTime;
	
	const time_t timeSave = 0.5*60*60;
	startTime = time(NULL);
	previousTime = startTime;
    
    //OUTPUT results files
	output_m = fopen("model_OUTPUT_ecg.txt", "w");
	output=fopen("model_OUTPUT_ap.txt", "w");
	
    
    //READ in initial conditions from 1D steady-state
	FILE *fp = fopen(SSstateFileName, "r");
    
	if (fp == NULL) {
		cout << "New start" << endl;
		
		S->counter=0;
        
        for ( ll = 0; ll < tl; ll++) {
            
            theCell = &(S->cellData[0][ll]);
            
            for ( int i = 0; i < sy0 ; i++ ) {
                theCell->y0[i] = y0[i];
            }
            
            theCell->Apd90 = 0;
            theCell->DI = 0;
            theCell->s2APD = 0;
            theCell->highestV = -86;
            theCell->lowestV = -86;
            theCell->highestV2 = -86;
            theCell->lowestV2 = -86;
            theCell->t1 = 1E10;
            theCell->t2 = 1E10;
            theCell->t3 = 1E10;
            theCell->t4 = 1E10;
            theCell->dvdt = 0;
            theCell->dvdt2 = 0;
            theCell->v_new = theCell->y0[0];
            theCell->I_inj = 0;
            
            //Set transmural cable
            if ( ll < 60 ) {
                theCell->celltype = 0; // %endo = 0, epi = 1, M = 2
                
            } else if ( ll < 105 ) {
                theCell->celltype = 2; // %endo = 0, epi = 1, M = 2
                
            } else {
                theCell->celltype = 1; // %endo = 0, epi = 1, M = 2
            }
            
            theCell->coef_PKA = 0; // beta-stimulation input range from 0 to 1 FOR linear change
            //For Control case, theCell->coef_PKA = 0 (no beta-stimulation)
        }
        
		S->beat = 1;
		S->t = 0.0;
		S->tt = 0.0;
		S->tstep = 0;
	} else {
		fread(S, sizeof(SimState), 1, fp);
		cout << "restarting from 1D steady-state " <<  endl;
		
		fclose(fp);
	}
	
	
	cycle_length = BCL;
	
    dt = base_dt;
    
	while (!done) {
		
		
        
        ww = 0;
        
#pragma omp parallel default(none) private( ll, theCell,  I_Total ) shared( S, dt, cout  )
        {
#pragma omp for schedule(static, 1)
            for ( ll = 0; ll < tl; ll++) {
                
                theCell = &(S->cellData[0][ll]);
            
                theCell->v = theCell->v_new;
                
                if ( S->t < 0.5 && ( ll < 2 ) ) {
                    theCell->I_inj = - 300;
                } else {
                    theCell->I_inj=0.0;
                }
                
                
                I_Total = Calcu_I_Total( theCell, dt );
                
                
                theCell->v = theCell->y0[0];
                
            }
        }
        
#pragma omp parallel default(none) private( ll, theCell, dv ) shared( cout , S, dt )
        {
#pragma omp for // schedule(static, 1)
            
            for ( ll = 0; ll < tl; ll++ ) {
                
                theCell = &(S->cellData[0][ll]);
                
                if( ll > 0 && ll < (tl-1) ) {
                    dv = dt * ( Df * ( S->cellData[0][ll-1].v - 2 * theCell->v + S->cellData[0][ll+1].v ) / ( dx * dx ) );
                }
                else if ( ll == 0 ) {
                    dv = dt * ( Df * ( -theCell->v + S->cellData[0][1].v ) / ( dx * dx ) );
                }
                else if ( ll == (tl-1) ) {
                    dv = dt * ( Df * ( -theCell->v + S->cellData[0][ll-1].v ) / ( dx * dx ) );
                }
                
                theCell->y0[0] += dv;
            }
            
        } // end openmp
        
        
        
        S->t += dt;
		S->tt += dt;
		S->counter += 1;
        
        
        
        
		
		if(  S->beat > 0 && ( S->counter % (fold*1) == 0 ) ) {
            
            
			fprintf ( output, "%10.3f\t",
					 S->tt );
            
            for ( ll = 0; ll < tl; ll+=1) {

              fprintf ( output, "%8.5f\t",S->cellData[0][ll].y0[0] );
            }
            
            fprintf ( output, "\n" );
            
            double Ev = 0;
            
            fprintf ( output_m, "%10.3f\t",
                     S->tt );
            
            for (ll=15; ll<=(tl-16); ll++){
                
                Ev = Ev + ( S->cellData[ww][ll].v - S->cellData[ww][ll+1].v )*( 1 / ((tl - ll - 1)*dx+2)-1/((tl - ll)*dx+2));
                
            }
            
            fprintf ( output_m, "%9.7g\t", Ev );
            
            
            fprintf ( output_m, "\n" );
            
		}
		
		
        if (S->counter % (fold*100) == 0) {
            cout<< "tt: " << S->tt << ", runtime: " << (time(NULL) - startTime)/60 << " min " << ", tstep: " << S->tstep << ", counter; " << S->counter << endl;
        }
		
		
		
		if(S->t >= cycle_length) {
			S->t = 0;
			S->beat++;
			cout << "now the beat is: " << S->beat << endl;
		}
		
		if (S->beat > stimuli) {
			S->beat = 1;
			done = 1;
		}
		
        //%%%%%%%%% For long simulations: save current state at a certain period of time in case of simulations stop
        
        
		if (time(NULL) - previousTime > timeSave) {
			previousTime = time(NULL);
			cout << "Saving current states." << endl;
			fp = fopen(stateFileName, "w");
			fwrite(S, sizeof(SimState), 1, fp);
			fclose(fp);
			cout << "Saved." << endl;
		}
		//%%%%%%%%%
        
        
	} // end while
	
    // Save cells for 2D reading
    cout << "Saving 1D cells for 2D.";
    const char *from1D = "model_INPUT.from1D";
    fp = fopen( from1D, "w");
    for ( ll = 0; ll < tl; ll ++ ) {
        theCell = &(S->cellData[0][ll]);
        fwrite( theCell, sizeof(Cell), 1, fp);
    }
    fclose(fp);
    cout << "Saved." << endl;

    
	S->beat=1;
    S->t = 0.0;
	S->tt = 0.0;
	S->counter = 0;
	S->tstep =0;
	
    fp = fopen(stateFileName, "w");
    fwrite(S, sizeof(SimState), 1, fp);
    fclose(fp);
    
    
	fclose ( output );
	
	fclose (output_m);
	return 0;
}


double Calcu_I_Total( Cell *theCell, double dt ) {
	
	double tspan[2];
	double tspan_total[] = { 0, dt }; // ms
	double partition = dt;  // ms
	int parts_total = floor( 0.5 + (tspan_total[1] - tspan_total[0]) / partition );
	int part;
	const int sy = 1;
	const int sy0 = 41;
	
	double y[sy][sy0], t1[sy];
	double * y1[sy];
	for( int i = 0; i < sy; i++ ) {
		y1[i] = y[i];
	}
	
	time_t t_start, t_stop;
	t_start = time(NULL);
	
	
	double V = theCell->y0[0];
	
	for( part = 0; part < parts_total; part++ ) {
		
		
		tspan[0] = partition * part ;
		tspan[1] = tspan[0] + partition ;
        
		integrate_rk2rk1( ord_model_betaStim , tspan, sy0, theCell->y0, dt, sy, y1 ,t1 , theCell );
		
		for( int i = 0; i < sy0; i++ ) {
			theCell->y0[i] = y1[sy-1][i];
		}
		
		
        
	}
	
	return ( - (theCell->y0[0] - V) / dt  );
	
    
}
