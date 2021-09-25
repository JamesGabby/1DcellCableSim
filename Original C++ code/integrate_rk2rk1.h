/*
 *  integrate_rk2rk1.h
 *
 *
 *  Created by Mao-Tsuen Jeng on 4/12/2016.
 *
 */

#ifndef integrate_rk2rk1_H
#define integrate_rk2rk1_H

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <vector>
// #include <omp.h>

void integrate_rk2rk1( void (*f)( double *, double * , Cell * , double ), double * tspan, int sy0, double * y0, double dt,  int sy, double * y[] , double * t1 , Cell * theCell );
// f : function that evaluate ydots ( 1st derivatives ).
// tspan : time interval for integration
// dt : time step size
// y0 : initial value
// t1 : array of saved time.
// sy0 : size of y0
// sy : size of y

using namespace std;

void integrate_rk2rk1( void (*f)( double *, double *, Cell * , double ), double * tspan, int sy0, double * y0, double dt,  int sy, double * y[] , double * t1 , Cell * theCell ){
    
	cout.precision(16);
    
	double tol = 0.001; 
    
	int counter, counter_sy ;
	double t, tn, err_y;
	std::vector<double> y1(sy0), y2(sy0), z(sy0), z2(sy0);
	double s, sm; 
	std::vector<double> dy1(sy0), dy2(sy0);
	std::vector<double> ydot(sy0); 
	char runType[] = "ydot";
	int i, j, k;
	int st = floor( 0.5 + ( tspan[1] - tspan[0] ) / ( dt * sy ) );
	
    
	for( i = 0; i < sy0; i++ ){
		y1[i] = y0[i];
		y2[i] = y0[i];
	}
	t = tspan[0];
    
    
	for( counter_sy = 0; counter_sy < sy; counter_sy ++ ) {
        
		for( counter = 0; counter < st; counter ++ ) {
            
			
			f(y1.data(), ydot.data(), theCell, dt);

            
			for( i = 0; i < sy0; i++ ) {
				dy1[i] = ydot[i] * dt; // k1
				y2[i] = y1[i] + dy1[i];
                z2[i] = y2[i];
			}
            
			f(y2.data(), ydot.data(), theCell, dt);

			for( i = 0; i < sy0; i++ ) {
				dy2[i] = ydot[i] * dt; // k2
                                       // Result using RK2
				z[i] = y1[i] + 0.5 * ( dy1[i] + dy2[i] );
			}
            
            
            
			sm = 1;
            
			for( i = 0; i < sy0; i++ ) {
                
				if ( z[i] != z2[i] && z[i] == z[i] ) {
					s = pow( tol * dt * ( fabs(z[i]) + fabs(z2[i]) ) / ( 2 * fabs( z[i] - z2[i] ) ) , 0.5 );
					
                    
					if( sm > s && s == s ) {
                        
						sm = s;
						if( 1E-1 > sm ) {
							sm = 1E-1;
							
						}
						
                        
					} else if ( s != s ) { // if s is NAN
                        
						sm = 1E-1;
					}
					//cout << endl;
				} else if ( z[i] != z[i] ){
					cout << z[i] << endl;  // z[i] is NAN
                    exit(1);
					sm = 1E-1;
				}
			}
			if( sm >= 1 ){
				// Keep RK2
				t += dt;
				for( i = 0; i < sy0; i++ ) {
					
					y1[i] = z[i];
					
				}
			} else {
                
				// Use RK2
                
                
                int invsm = 8;//1 / sm;
				double dt2 = dt / invsm; //sm;
				for ( int id1 = 1; id1 <= invsm ; id1 ++ ) {
                    
					f(y1.data(), ydot.data(), theCell, dt2);

                    
                    for( i = 0; i < sy0; i++ ) {
                        dy1[i] = ydot[i] * dt2; // k1
                        y2[i] = y1[i] + dy1[i];
                        
                    }
                    
					f(y2.data(), ydot.data(), theCell, dt2);

                    for( i = 0; i < sy0; i++ ) {
                        dy2[i] = ydot[i] * dt2; // k2
                                               // Result using RK2
                        y1[i] = y1[i] + 0.5 * ( dy1[i] + dy2[i] );
                    }
                    
					t += dt2;
					
				}
			}
            
            
		}
        
		
		t1[counter_sy] = t;
        
		for( i = 0; i < sy0; i++ ) {
			y[counter_sy][i] = y1[i];
					
		}
        
        
	}
	
    
}

#endif
