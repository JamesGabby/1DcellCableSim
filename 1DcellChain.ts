/*
 *  O'Hara-Rudy Human model with RK method
 *  TypeScript version based on the C++ implementation created by Mao-Tsuen Jeng on 6/13/11.
 *
 *	Copyright 2011 __Colleen Clancy Lab__. All rights reserved.
 *
 *  One-D simulations
 */

const BCL = 1000;
const stimuli = 2;
const tl = 3; //length of the tissue, 165 for full length
const tw = 1; //width of the tissue
const base_dt = 0.005;
const fold = 1.0/base_dt;
const dx= 0.01; //cm
const gj_use = 2.5; //uS
const rad = 0.0011;
const Rmyo = 150.0; //ohm*cm
const Rg = 3.14159*Math.pow(rad,2)/(gj_use*1.0e-6);
const Df = (1e3*rad)/(4.0*(Rmyo+Rg/dx));

// Initial conditions
let v=-87.;
let nai=7.;
let nass=nai;
let ki=145.;
let kss=ki;
let cai=1.0e-4;
let cass=cai;
let cansr=1.2;
let cajsr=cansr;
let m=0.;
let hf=1.;
let hs=1.;
let j=1.;
let hsp=1.;
let jp=1.;
let mL=0.;
let hL=1.;
let hLp=1.;
let a=0.;
let iF=1.;
let iS=1.;
let ap=0.;
let iFp=1.;
let iSp=1.;
let d=0.;
let ff=1.;
let fs=1.;
let fcaf=1.;
let fcas=1.;
let jca=1.;
let nca=0.;
let ffp=1.;
let fcafp=1.;
let xrf=0.;
let xrs=0.;
let xs1=0.;
let xs2=0.;
let xk1=1.;
let Jrelnp=0.;
let Jrelp=0.;
let CaMKt=0.;
let sy0 = 41;
let y0 = [ v, nai, nass, ki, kss,
	cai, cass, cansr, cajsr, m,
	hf, hs, j, hsp, jp,
	mL, hL, hLp, a, iF,
	iS, ap, iFp, iSp, d,
	ff, fs, fcaf, fcas, jca,
	nca, ffp, fcafp, xrf, xrs,
	xs1, xs2, xk1, Jrelnp, Jrelp,
	CaMKt ];
let p = 0;
let ll, ww;
let V, dvdt, dvdt2;
let dt;
let I_Total;
let dv;
let Vm = -86;
let cycle_length;
let I_inj;
let done = 0;
let waitTime;

class Cell {
	public y0 = [];
	public highestV; public lowestV; public Apd90;
	public highestdvdt; public V90;
	public DI; public t1; public t2;
	public highestV2; public lowestV2; public t3; public t4; public s2APD;
	public dvdt; public v_new; public v; public dvdt2;
	public currents = [];
	public celltype;
	public I_inj;
	public coef_PKA;

	constructor(y0: number[]) {
		this.y0 = y0;
	}
} 

class SimState {
	public t; public tt;
	public tstep; public counter; public beat;
	public cellData = new Array();

	constructor(cellData: Cell[]) {
		this.cellData = cellData;
		for (var i = 0; i <= tl; i++) {
			this.cellData[i] = [new Cell([])];
		}
	}
}

let cell = new Cell(y0);
let theState = new SimState([cell]);
let theCell: Cell = new Cell(y0);

function main() {
    console.log("New start");
    theState.counter = 0;
    for (ll = 0; ll < tl; ll++) {
        theCell = theState.cellData[ll][0];
        for (let i = 0; i < sy0; i++) {
            theCell.y0[i] = y0[i];
        }
        theCell.v = v;
        theCell.Apd90 = 0;
        theCell.DI = 0;
        theCell.s2APD = 0;
        theCell.highestV = -86;
        theCell.lowestV = -86;
        theCell.highestV2 = -86;
        theCell.lowestV2 = -86;
        theCell.t1 = 1E10;
        theCell.t2 = 1E10;
        theCell.t3 = 1E10;
        theCell.t4 = 1E10;
        theCell.dvdt = 0;
        theCell.dvdt2 = 0;
        theCell.v_new = theCell.y0[0];
        theCell.I_inj = 0;
		
        //Set transmural cable
        if (ll < 60) {
            theCell.celltype = 0; // %endo = 0, epi = 1, M = 2
        } else if (ll < 105) {
            theCell.celltype = 2; // %endo = 0, epi = 1, M = 2
        } else {
            theCell.celltype = 1; // %endo = 0, epi = 1, M = 2
        }
        theCell.coef_PKA = 0; // beta-stimulation input range from 0 to 1 FOR linear change
        //For Control case, theCell.coef_PKA = 0 (no beta-stimulation)
    }
	theState.beat = 1;
	theState.t = 0.0;
	theState.tt = 0.0;
	theState.tstep = 0;
	cycle_length = BCL;
    dt = base_dt;
	while (!done) {
		ww = 0;
		for (ll = 0; ll < tl; ll++) {
			theCell = theState.cellData[ll][0];
			theCell.v = theCell.v_new;
			if (theState.t < 0.5 && ( ll < 2 )) {
				theCell.I_inj = -300;
			} else {
				theCell.I_inj = 0.0;
			}
			I_Total = Calcu_I_Total(theCell, dt);
			theCell.v = theCell.y0[0];
		}

		for (ll = 0; ll < tl; ll++) {
			//theCell = theState.cellData[ll][0];
			if (ll > 0 && ll < (tl-1)) {
				dv = dt * ( Df * ( theState.cellData[ll-1][0].v - 2 * theCell.v + theState.cellData[ll+1][0].v ) / ( dx * dx ) );
			} else if ( ll == 0 ) {
				dv = dt * ( Df * ( -theCell.v + theState.cellData[1][0].v ) / ( dx * dx ) );
			} else if ( ll == (tl-1) ) {
				dv = dt * ( Df * ( -theCell.v + theState.cellData[ll-1][0].v ) / ( dx * dx ) );
			}	
			theCell.y0[0] += dv;
		}

		theState.t += dt;
		theState.tt += dt;
		theState.counter += 1;

		if (theState.beat > 0 && (theState.counter % (fold*1) == 0 )) {
			// fprintf ( output, "%10.3f\t",
			// 		 theState.tt );
		//	console.log("theState.tt: "+theState.tt.toFixed());
            
            for (ll = 0; ll < tl; ll+=1) {
              //fprintf ( output, "%8.5f\t",theState.cellData[0][ll].y0[0] );
			  console.log("theState.cellData[ll][0].y0[0]: "+theState.cellData[ll][0].y0[0]);
            }
            
            //fprintf ( output, "\n" );
            let Ev = 0;
            // fprintf ( output_m, "%10.3f\t",
            //          theState.tt );
            
            for (ll = 15; ll <= (tl-16); ll++) {
                Ev = Ev + (theState.cellData[ll][ww].v - theState.cellData[ll+1][ww].v ) * ( 1 / ((tl - ll - 1) * dx+2)-1 / ((tl - ll) * dx+2));
            }
            
            // fprintf ( output_m, "%9.7g\t", Ev );
            // fprintf ( output_m, "\n" );
		}
		
        if (theState.counter % (fold*100) == 0) {
            console.log("tt: " + theState.tt /*+ ", runtime: " + (time(NULL) - startTime)/60 + " min "*/ + ", tstep: " + theState.tstep + ", counter: " + theState.counter);
        }
		if (theState.t >= cycle_length) {
			theState.t = 0;
			theState.beat++;
			console.log("now the beat is: " + theState.beat);
		}
		if (theState.beat > stimuli) {
			theState.beat = 1;
			done = 1;
		}
	} //end of while
	
	// Save cells for 2D reading
    console.log("Saving 1D cells for 2D.");
    //const char *from1D = "model_INPUT.from1D";
    //fp = fopen( from1D, "w");
    for ( ll = 0; ll < tl; ll ++ ) {
        //theCell = theState.cellData[ll][0];
        //fwrite( theCell, sizeof(Cell), 1, fp);
		console.log("theCell: "+theCell);
    }
    //fclose(fp);
    console.log("Saved.");
	theState.beat=1;
    theState.t = 0.0;
	theState.tt = 0.0;
	theState.counter = 0;
	theState.tstep =0;
	console.log("theState: "+theState);
}

function Calcu_I_Total(theCell: Cell, dt: number): number {
	theCell = new Cell(y0);
	let tspan = [];
	let tspan_total = [ 0, dt ]; // ms
	let partition = dt;  // ms
	let parts_total = Math.floor( 0.5 + (tspan_total[1] - tspan_total[0]) / partition );
	let part;
	const sy = 1;
	const sy0 = 41;
	let y = new Array(sy0);
	let t1 = new Array(sy);
	for (let j=0; j < y.length; j++) {
		y[j] = new Array(sy);
	}
	let y1 = new Array(sy);
	for( let i = 0; i < sy; i++ ) {
		y1[i] = y[i];
	}
	V = theCell.y0[0];
	for (part = 0; part < parts_total; part++) {
		tspan[0] = partition * part ;
		tspan[1] = tspan[0] + partition ;
		integrate_rk2rk1(ord_model_betaStim, tspan, sy0, theCell.y0, dt, sy, y1, t1, theCell);
		for (let i = 0; i < sy0; i++) {
			 theCell.y0[i] = y1[sy-1][i];
			// console.log(theCell.y0[i]);
		}
	}
	//console.log(-(theCell.y0[0] - V) / dt);	
	return (-(theCell.y0[0] - V) / dt);	
}

//=====================================================================================================

function integrate_rk2rk1(f: (a, b, c: Cell , d: number) => void, tspan: any[], sy0: number, y0: any[], dt: number, sy: number, y: number[], t1: any[], theCell: Cell) {
	let tol = 0.001; 
	let counter, counter_sy;
	let t, tn, err_y;
	let y1=[], y2=[], z=[], z2=[];
	let s, sm; 
	let dy1=[], dy2=[];
	let ydot=[];
	let runType = ['y','d','o','t'];
	let i, j, k;
	let st = Math.floor(0.5 + ( tspan[1] - tspan[0] ) / ( dt * sy ));
	for (i = 0; i < sy0; i++) {
		y1[i] = y0[i];
		y2[i] = y0[i];
	}
	t = tspan[0];
	for (counter_sy = 0; counter_sy < sy; counter_sy++) {
		for (counter = 0; counter < st; counter++) {
			f(y1, ydot, theCell, dt); // [?]
			for (i = 0; i < sy0; i++) {
				dy1[i] = ydot[i] * dt; // k1
				y2[i] = y1[i] + dy1[i];
                z2[i] = y2[i];
			}
			f(y2, ydot, theCell, dt);
			for( i = 0; i < sy0; i++ ) {
				dy2[i] = ydot[i] * dt; // k2
				z[i] = y1[i] + 0.5 * ( dy1[i] + dy2[i] ); // Result using RK2
			}
			sm = 1;
			for (i = 0; i < sy0; i++) {
				if (z[i] != z2[i] && z[i] == z[i]) {
					s = Math.pow( tol * dt * ( Math.abs(z[i]) + Math.abs(z2[i]) ) / ( 2 * Math.abs( z[i] - z2[i] ) ) , 0.5 );
					if (sm > s && s == s) {     
						sm = s;
						if (1E-1 > sm) {
							sm = 1E-1;
						}
					} else if (s != s) { // if s is NAN   
						sm = 1E-1;
					}
					//console.log( << endl;
				} else if ( z[i] != z[i] ){
					//console.log(z[i]);  // z[i] is NAN
					//break; // [?]
					sm = 1E-1;
				}
			}
			if (sm >= 1) {
				// Keep RK2
				t += dt;
				for (i = 0; i < sy0; i++ ) {
					y1[i] = z[i];
				}
			} else {
				// Use RK2
                let invsm = 8;//1 / sm;
				let dt2 = dt / invsm; //sm;
				for ( let id1 = 1; id1 <= invsm ; id1 ++ ) {
					f(y1, ydot, theCell, dt2);
                    
                    for ( i = 0; i < sy0; i++ ) {
                        dy1[i] = ydot[i] * dt2; // k1
                        y2[i] = y1[i] + dy1[i];
                    }
                    
					f(y2, ydot, theCell, dt2);

                    for ( i = 0; i < sy0; i++ ) {
                        dy2[i] = ydot[i] * dt2; // k2
                                               // Result using RK2
                        y1[i] = y1[i] + 0.5 * ( dy1[i] + dy2[i] );
                    }  
					t += dt2;
				}
			}
		}
		t1[counter_sy] = t;
		for (i = 0; i < sy0; i++ ) {
			y[counter_sy][i] = y1[i];			
		}   
	}	
}

//========================================================================================================

function ord_model_betaStim(X, Xdot, theCell: Cell, dt: number) {
	let celltype = theCell.celltype; // %endo = 0, epi = 1, M = 2
    let Istim = theCell.I_inj;
    let coef_PKA = theCell.coef_PKA;
    let currents = theCell.currents;
	X = [];
	Xdot = [];
    
	// %extracellular ionic concentrations
	const nao = 140.0;
	const cao = 1.8;
	const ko = 5.4;
	
	// %physical constants
	const R = 8314.0;
	const T = 310.0;
	const F = 96485.0;
	
	// %cell geometry
	const L = 0.01;
	const rad = 0.0011;
	const vcell = 1000*3.14*rad*rad*L;
	const Ageo = 2*3.14*rad*rad+2*3.14*rad*L;
	const Acap = 2*Ageo;
	const vmyo = 0.68*vcell;
	const vnsr = 0.0552*vcell;
	const vjsr = 0.0048*vcell;
	const vss = 0.02*vcell;
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
    if (X[6] > 0.03) {  // cass <= 0.03 [mM]
        X[6] = 0.03;
    }
    
	// %give names to the state vector values
	const v=X[0];
	const nai=X[1];
	const nass=X[2];
	const ki=X[3];
	const kss=X[4];
	const cai=X[5];
	const cass=X[6];
	const cansr=X[7];
	const cajsr=X[8];
	const m=X[9];
	const hf=X[10];
	const hs=X[11];
	const j=X[12];
	const hsp=X[13];
	const jp=X[14];
	const mL=X[15];
	const hL=X[16];
	const hLp=X[17];
	const a=X[18];
	const iF=X[19];
	const iS=X[20];
	const ap=X[21];
	const iFp=X[22];
	const iSp=X[23];
	const d=X[24];
	const ff=X[25];
	const fs=X[26];
	const fcaf=X[27];
	const fcas=X[28];
	const jca=X[29];
	const nca=X[30];
	const ffp=X[31];
	const fcafp=X[32];
	const xrf=X[33];
	const xrs=X[34];
	const xs1=X[35];
	const xs2=X[36];
	const xk1=X[37];
	const Jrelnp=X[38];
	const Jrelp=X[39];
	const CaMKt=X[40];
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	// %CaMK constants
	const KmCaMK=0.15;
	
	const aCaMK=0.05;
	const bCaMK=0.00068;
	const CaMKo=0.05;
	const KmCaM=0.0015;
	// %update CaMK
	const CaMKb=CaMKo*(1.0-CaMKt)/(1.0+KmCaM/cass);
	const CaMKa=CaMKb+CaMKt;
	const dCaMKt=aCaMK*CaMKb*(CaMKb+CaMKt)-bCaMK*CaMKt;
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	// %reversal potentials
	let ENa=(R*T/F)*Math.log(nao/nai);
	let EK=(R*T/F)*Math.log(ko/ki);
	let PKNa=0.01833;
	let EKs=(R*T/F)*Math.log((ko+PKNa*nao)/(ki+PKNa*nai));
	
	// %convenient shorthand calculations
	let vffrt=v*F*F/(R*T);
	let vfrt=v*F/(R*T);
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	// %calculate INa
	let mss=1.0/(1.0+Math.exp((-(v+39.57))/9.871));
	let tm=1.0/(6.765*Math.exp((v+11.64)/34.77)+8.552*Math.exp(-(v+77.42)/5.955));
	let dm=(mss-m)/tm;
	let hss = 1.0 / ( 1. + Math.exp( ( v + 82.90 + coef_PKA * 5. ) / 6.086 ) ); //!PKA effect (+0~5)
	let thf=1.0/(1.432e-5*Math.exp(-(v+1.196)/6.285)+6.149*Math.exp((v+0.5096)/20.27));
	let ths=1.0/(0.009794*Math.exp(-(v+17.95)/28.05)+0.3343*Math.exp((v+5.730)/56.66));
	let Ahf=0.99;
	let Ahs=1.0-Ahf;
	let dhf=(hss-hf)/thf;
	let dhs=(hss-hs)/ths;
	let h=Ahf*hf+Ahs*hs;
	let jss=hss;
	let tj=2.038+1.0/(0.02136*Math.exp(-(v+100.6)/8.281)+0.3052*Math.exp((v+0.9941)/38.45));
	let dj=(jss-j)/tj;
	let hssp = 1.0 / ( 1. + Math.exp( ( v + 89.1 + coef_PKA * 5. ) / 6.086 ) ); //!PKA effect (+0~5)
	let thsp=3.0*ths;
	let dhsp=(hssp-hsp)/thsp;
	let hp=Ahf*hf+Ahs*hsp;
	let tjp=1.46*tj;
	let djp=(jss-jp)/tjp;
	let GNa = 75. * ( 1. + coef_PKA * 1.7 ); //!PKA effect (*1~2.7)
	let fINap=(1.0/(1.0+KmCaMK/CaMKa));
	let INa=GNa*(v-ENa)*m*m*m*((1.0-fINap)*h*j+fINap*hp*jp);
	
	// %calculate INaL
	let mLss=1.0/(1.0+Math.exp((-(v+42.85))/5.264));
	let tmL=tm;
	let dmL=(mLss-mL)/tmL;
	let hLss=1.0/(1.0+Math.exp((v+87.61)/7.488));
	let thL=200.0;
	let dhL=(hLss-hL)/thL;
	let hLssp=1.0/(1.0+Math.exp((v+93.81)/7.488));
	let thLp=3.0*thL;
	let dhLp=(hLssp-hLp)/thLp;
	let GNaL=0.0075;
	if ( celltype==1 ) {
		GNaL=GNaL*0.6;
	}
	let fINaLp=(1.0/(1.0+KmCaMK/CaMKa));
	let INaL=GNaL*(v-ENa)*mL*((1.0-fINaLp)*hL+fINaLp*hLp);
	
	// %calculate Ito
	let ass=1.0/(1.0+Math.exp((-(v-14.34))/14.82));
	let ta=1.0515/(1.0/(1.2089*(1.0+Math.exp(-(v-18.4099)/29.3814)))+3.5/(1.0+Math.exp((v+100.0)/29.3814)));
	let da=(ass-a)/ta;
	let iss=1.0/(1.0+Math.exp((v+43.94)/5.711));
	let delta_epi;
	if ( celltype==1 ) {
		delta_epi=1.0-(0.95/(1.0+Math.exp((v+70.0)/5.0)));
	} else {
		delta_epi=1.0;
	}
	let tiF=4.562+1/(0.3933*Math.exp((-(v+100.0))/100.0)+0.08004*Math.exp((v+50.0)/16.59));
	let tiS=23.62+1/(0.001416*Math.exp((-(v+96.52))/59.05)+1.780e-8*Math.exp((v+114.1)/8.079));
	tiF=tiF*delta_epi;
	tiS=tiS*delta_epi;
	let AiF=1.0/(1.0+Math.exp((v-213.6)/151.2));
	let AiS=1.0-AiF;
	let diF=(iss-iF)/tiF;
	let diS=(iss-iS)/tiS;
	let i=AiF*iF+AiS*iS;
	let assp=1.0/(1.0+Math.exp((-(v-24.34))/14.82));
	let dap=(assp-ap)/ta;
	let dti_develop=1.354+1.0e-4/(Math.exp((v-167.4)/15.89)+Math.exp(-(v-12.23)/0.2154));
	let dti_recover=1.0-0.5/(1.0+Math.exp((v+70.0)/20.0));
	let tiFp=dti_develop*dti_recover*tiF;
	let tiSp=dti_develop*dti_recover*tiS;
	let diFp=(iss-iFp)/tiFp;
	let diSp=(iss-iSp)/tiSp;
	let ip=AiF*iFp+AiS*iSp;
	let Gto=0.02;
	if( celltype==1 ){
		Gto=Gto*4.0;
	} else if( celltype==2 ){
		Gto=Gto*4.0;
	}
	let fItop=(1.0/(1.0+KmCaMK/CaMKa));
	let Ito=Gto*(v-EK)*((1.0-fItop)*a*i+fItop*ap*ip);
	
	// %calculate ICaL, ICaNa, ICaK
	let dss = 1.0 / ( 1.0 + Math.exp( ( -( v + 3.940 + coef_PKA * 16. ) ) / 4.230 ) ); //!PKA effect (+0~16)
	let td=0.6+1.0/(Math.exp(-0.05*(v+6.0))+Math.exp(0.09*(v+14.0)));
	let dd=(dss-d)/td;
	let fss = 1.0 / ( 1.0 + Math.exp( ( v + 19.58 + coef_PKA * 8. ) / 3.696 ) ); //!PKA effect (+0~8)
	let tff=7.0+1.0/(0.0045*Math.exp(-(v+20.0)/10.0)+0.0045*Math.exp((v+20.0)/10.0));
	let tfs=1000.0+1.0/(0.000035*Math.exp(-(v+5.0)/4.0)+0.000035*Math.exp((v+5.0)/6.0));
	let Aff=0.6;
	let Afs=1.0-Aff;
	let dff=(fss-ff)/tff;
	let dfs=(fss-fs)/tfs;
	let f=Aff*ff+Afs*fs;
	let fcass=fss;
	let tfcaf=7.0+1.0/(0.04*Math.exp(-(v-4.0)/7.0)+0.04*Math.exp((v-4.0)/7.0));
	let tfcas=100.0+1.0/(0.00012*Math.exp(-v/3.0)+0.00012*Math.exp(v/7.0));
	let Afcaf=0.3+0.6/(1.0+Math.exp((v-10.0)/10.0));
	let Afcas=1.0-Afcaf;
	let dfcaf=(fcass-fcaf)/tfcaf;
	let dfcas=(fcass-fcas)/tfcas;
	let fca=Afcaf*fcaf+Afcas*fcas;
	let tjca=75.0;
	let djca=(fcass-jca)/tjca;
	let tffp=2.5*tff;
	let dffp=(fss-ffp)/tffp;
	let fp=Aff*ffp+Afs*fs;
	let tfcafp=2.5*tfcaf;
	let dfcafp=(fcass-fcafp)/tfcafp;
	let fcap=Afcaf*fcafp+Afcas*fcas;
	let Kmn=0.002;
	let k2n=1000.0;
	let km2n=jca*1.0;
	let anca=1.0 / ( k2n/km2n+ Math.pow((1.0+Kmn/cass),4) );
	let dnca=anca*k2n-nca*km2n;
	
    //!PKA effect cass upto 0.03 mM
    let PhiCaL=4.0*vffrt*(cass*Math.exp(2.0*vfrt)-0.341*cao)/(Math.exp(2.0*vfrt)-1.0);
    
	let PhiCaNa=1.0*vffrt*(0.75*nass*Math.exp(1.0*vfrt)-0.75*nao)/(Math.exp(1.0*vfrt)-1.0);
	let PhiCaK=1.0*vffrt*(0.75*kss*Math.exp(1.0*vfrt)-0.75*ko)/(Math.exp(1.0*vfrt)-1.0);
	let zca=2.0;
    
	let PCa = 0.0001 * ( 1. + coef_PKA * 1.5 ) ; //!PKA effect (*1~2.5)
    
	if ( celltype==1 ) {
		PCa=PCa*1.2;
	} else if ( celltype==2 ) {
		PCa=PCa*2.5;
	}
	let PCap=1.1*PCa;
	let PCaNa=0.00125*PCa;
	let PCaK=3.574e-4*PCa;
	let PCaNap=0.00125*PCap;
	let PCaKp=3.574e-4*PCap;
	let fICaLp=(1.0/(1.0+KmCaMK/CaMKa));
	let ICaL=(1.0-fICaLp)*PCa*PhiCaL*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCap*PhiCaL*d*(fp*(1.0-nca)+jca*fcap*nca);
	let ICaNa=(1.0-fICaLp)*PCaNa*PhiCaNa*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCaNap*PhiCaNa*d*(fp*(1.0-nca)+jca*fcap*nca);
	let ICaK=(1.0-fICaLp)*PCaK*PhiCaK*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCaKp*PhiCaK*d*(fp*(1.0-nca)+jca*fcap*nca);
	
	// %calculate IKr
	let xrss=1.0/(1.0+Math.exp((-(v+8.337))/6.789));
	let txrf=12.98+1.0/(0.3652*Math.exp((v-31.66)/3.869)+4.123e-5*Math.exp((-(v-47.78))/20.38));
	let txrs=1.865+1.0/(0.06629*Math.exp((v-34.70)/7.355)+1.128e-5*Math.exp((-(v-29.74))/25.94));
	let Axrf=1.0/(1.0+Math.exp((v+54.81)/38.21));
	let Axrs=1.0-Axrf;
	let dxrf=(xrss-xrf)/txrf;
	let dxrs=(xrss-xrs)/txrs;
	let xr=Axrf*xrf+Axrs*xrs;
	let rkr=1.0/(1.0+Math.exp((v+55.0)/75.0))*1.0/(1.0+Math.exp((v-10.0)/30.0));
	let GKr=0.046;
	if( celltype==1 ) {
		GKr=GKr*1.3;
	} else if( celltype==2 ) {
		GKr=GKr*0.8;
	}
	let IKr=GKr*Math.sqrt(ko/5.4)*xr*rkr*(v-EK);
	
	// %calculate IKs
	let xs1ss = 1.0 / ( 1.0 + Math.exp( ( -( v + 11.60 ) ) / 8.932 ) ) * ( 1. + coef_PKA * 0.4 ); //!PKA effect (*1~1.2)
	let txs1=817.3+1.0/(2.326e-4*Math.exp((v+48.28)/17.80)+0.001292*Math.exp((-(v+210.0))/230.0));
	let dxs1=(xs1ss-xs1)/txs1;
	let xs2ss = xs1ss;
	let txs2=1.0/(0.01*Math.exp((v-50.0)/20.0)+0.0193*Math.exp((-(v+66.54))/31.0));
	let dxs2=(xs2ss-xs2)/txs2;
	let KsCa=1.0+0.6/(1.0+ Math.pow((3.8e-5/cai),1.4) );
    
	let GKs = 0.0034 * ( 1. + coef_PKA * 2.2 ); //!PKA effect (*1~3.2)
    
	if ( celltype==1 ) {
		GKs=GKs*1.4;
	}
	let IKs=GKs*KsCa*xs1*xs2*(v-EKs);
	
	let xk1ss=1.0/(1.0+Math.exp(-(v+2.5538*ko+144.59)/(1.5692*ko+3.8115)));
	let txk1=122.2/(Math.exp((-(v+127.2))/20.36)+Math.exp((v+236.8)/69.33));
	let dxk1=(xk1ss-xk1)/txk1;
	let rk1=1.0/(1.0+Math.exp((v+105.8-2.6*ko)/9.493));
	let GK1=0.1908;
	if ( celltype==1 ) {
		GK1=GK1*1.2;
	} else if( celltype==2 ) {
		GK1=GK1*1.3;
	}
	let IK1=GK1*Math.sqrt(ko)*rk1*xk1*(v-EK);
	
	// %calculate INaCa_i
	let kna1=15.0;
	let kna2=5.0;
	let kna3=88.12;
	let kasymm=12.5;
	let wna=6.0e4;
	let wca=6.0e4;
	let wnaca=5.0e3;
	let kcaon=1.5e6;
	let kcaoff=5.0e3;
	let qna=0.5224;
	let qca=0.1670;
	let hca=Math.exp((qca*v*F)/(R*T));
	let hna=Math.exp((qna*v*F)/(R*T));
	let h1=1+nai/kna3*(1+hna);
	let h2=(nai*hna)/(kna3*h1);
	let h3=1.0/h1;
	let h4=1.0+nai/kna1*(1+nai/kna2);
	let h5=nai*nai/(h4*kna1*kna2);
	let h6=1.0/h4;
	let h7=1.0+nao/kna3*(1.0+1.0/hna);
	let h8=nao/(kna3*hna*h7);
	let h9=1.0/h7;
	let h10=kasymm+1.0+nao/kna1*(1.0+nao/kna2);
	let h11=nao*nao/(h10*kna1*kna2);
	let h12=1.0/h10;
	let k1=h12*cao*kcaon;
	let k2=kcaoff;
	let k3p=h9*wca;
	let k3pp=h8*wnaca;
	let k3=k3p+k3pp;
	let k4p=h3*wca/hca;
	let k4pp=h2*wnaca;
	let k4=k4p+k4pp;
	let k5=kcaoff;
	let k6=h6*cai*kcaon;
	let k7=h5*h2*wna;
	let k8=h8*h11*wna;
	let x1=k2*k4*(k7+k6)+k5*k7*(k2+k3);
	let x2=k1*k7*(k4+k5)+k4*k6*(k1+k8);
	let x3=k1*k3*(k7+k6)+k8*k6*(k2+k3);
	let x4=k2*k8*(k4+k5)+k3*k5*(k1+k8);
	let E1=x1/(x1+x2+x3+x4);
	let E2=x2/(x1+x2+x3+x4);
	let E3=x3/(x1+x2+x3+x4);
	let E4=x4/(x1+x2+x3+x4);
	let KmCaAct=150.0e-6;
	let allo=1.0/(1.0+ Math.pow((KmCaAct/cai),2));
	let zna=1.0;
	let JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
	let JncxCa=E2*k2-E1*k1;
	let Gncx=0.0008;
	if ( celltype==1 ) {
		Gncx=Gncx*1.1;
	} else if( celltype==2 ) {
		Gncx=Gncx*1.4;
	}
	let INaCa_i=0.8*Gncx*allo*(zna*JncxNa+zca*JncxCa);
	
	// %calculate INaCa_ss
	h1=1+nass/kna3*(1+hna);
	h2=(nass*hna)/(kna3*h1);
	h3=1.0/h1;
	h4=1.0+nass/kna1*(1+nass/kna2);
	h5=nass*nass/(h4*kna1*kna2);
	h6=1.0/h4;
	h7=1.0+nao/kna3*(1.0+1.0/hna);
	h8=nao/(kna3*hna*h7);
	h9=1.0/h7;
	h10=kasymm+1.0+nao/kna1*(1+nao/kna2);
	h11=nao*nao/(h10*kna1*kna2);
	h12=1.0/h10;
	k1=h12*cao*kcaon;
	k2=kcaoff;
	k3p=h9*wca;
	k3pp=h8*wnaca;
	k3=k3p+k3pp;
	k4p=h3*wca/hca;
	k4pp=h2*wnaca;
	k4=k4p+k4pp;
	k5=kcaoff;
	k6=h6*cass*kcaon;
	k7=h5*h2*wna;
	k8=h8*h11*wna;
	x1=k2*k4*(k7+k6)+k5*k7*(k2+k3);
	x2=k1*k7*(k4+k5)+k4*k6*(k1+k8);
	x3=k1*k3*(k7+k6)+k8*k6*(k2+k3);
	x4=k2*k8*(k4+k5)+k3*k5*(k1+k8);
	E1=x1/(x1+x2+x3+x4);
	E2=x2/(x1+x2+x3+x4);
	E3=x3/(x1+x2+x3+x4);
	E4=x4/(x1+x2+x3+x4);
	KmCaAct=150.0e-6;
	allo=1.0/(1.0+ Math.pow((KmCaAct/cass),2));
	JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
	JncxCa=E2*k2-E1*k1;
	let INaCa_ss=0.2*Gncx*allo*(zna*JncxNa+zca*JncxCa);
	
	// %calculate INaK
	let k1p=949.5;
	let k1m=182.4;
	let k2p=687.2;
	let k2m=39.4;
	k3p=1899.0;
	let k3m=79300.0;
	k4p=639.0;
	let k4m=40.0;
	let Knai0=9.073 * ( 1. - coef_PKA * 0.3 ); //!PKA effect (*0.7~1);
	let Knao0=27.78;
	let delta=-0.1550;
    let Knai = Knai0 * Math.exp( ( delta * v * F ) / ( 3.0 * R * T ) );
	let Knao=Knao0*Math.exp(((1.0-delta)*v*F)/(3.0*R*T));
	let Kki=0.5;
	let Kko=0.3582;
	let MgADP=0.05;
	let MgATP=9.8;
	let Kmgatp=1.698e-7;
	let H=1.0e-7;
	let eP=4.2;
	let Khp=1.698e-7;
	let Knap=224.0;
	let Kxkur=292.0;
	let P=eP/(1.0+H/Khp+nai/Knap+ki/Kxkur);
	let a1=(k1p * Math.pow((nai/Knai),3) )/( Math.pow((1.0+nai/Knai),3) + Math.pow((1.0+ki/Kki),2) - 1.0);
	let b1=k1m*MgADP;
	let a2=k2p;
	let b2=(k2m * Math.pow((nao/Knao),3) )/( Math.pow((1.0+nao/Knao),3) + Math.pow((1.0+ko/Kko),2) - 1.0);
	let a3=(k3p * Math.pow((ko/Kko),2) )/( Math.pow((1.0+nao/Knao),3) + Math.pow((1.0+ko/Kko),2) - 1.0);
	let b3=(k3m*P*H)/(1.0+MgATP/Kmgatp);
	let a4=(k4p*MgATP/Kmgatp)/(1.0+MgATP/Kmgatp);
	let b4=(k4m * Math.pow((ki/Kki),2) )/( Math.pow((1.0+nai/Knai),3) + Math.pow((1.0+ki/Kki),2) - 1.0);
	x1=a4*a1*a2+b2*b4*b3+a2*b4*b3+b3*a1*a2;
	x2=b2*b1*b4+a1*a2*a3+a3*b1*b4+a2*a3*b4;
	x3=a2*a3*a4+b3*b2*b1+b2*b1*a4+a3*a4*b1;
	x4=b4*b3*b2+a3*a4*a1+b2*a4*a1+b3*b2*a1;
	E1=x1/(x1+x2+x3+x4);
	E2=x2/(x1+x2+x3+x4);
	E3=x3/(x1+x2+x3+x4);
	E4=x4/(x1+x2+x3+x4);
	let zk=1.0;
	let JnakNa=3.0*(E1*a3-E2*b3);
	let JnakK=2.0*(E4*b1-E3*a1);
	let Pnak=30;
	if ( celltype==1 ) {
		Pnak=Pnak*0.9;
	} else if( celltype==2 ) {
		Pnak=Pnak*0.7;
	}
	let INaK=Pnak*(zna*JnakNa+zk*JnakK);
	
	// %calculate IKb
	let xkb=1.0/(1.0+Math.exp(-(v-14.48)/18.34));
	let GKb = 0.003 * ( 1. + coef_PKA * 1.5 ); //!PKA effect cass  (*1~2.5)
	if ( celltype==1 ) {
		GKb=GKb*0.6;
	}
	let IKb=GKb*xkb*(v-EK);
	
	// %calculate INab
	let PNab=3.75e-10;
	let INab=PNab*vffrt*(nai*Math.exp(vfrt)-nao)/(Math.exp(vfrt)-1.0);
	
	// %calculate ICab
	let PCab=2.5e-8;
	let ICab=PCab*4.0*vffrt*(cai*Math.exp(2.0*vfrt)-0.341*cao)/(Math.exp(2.0*vfrt)-1.0);
	
	// %calculate IpCa
	let GpCa=0.0005;
	let IpCa=GpCa*cai/(0.0005+cai);
	
		
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	// %udtate the membrane voltage
	let dv=-(INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa_i+INaCa_ss+INaK+INab+IKb+IpCa+ICab+Istim);
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	// %calculate diffusion fluxes
	let JdiffNa=(nass-nai)/2.0;
	let JdiffK=(kss-ki)/2.0;
	let Jdiff=(cass-cai)/0.2;
	
	// %calculate ryanodione receptor calcium induced calcium release from the jsr
	let bt=4.75;
	let a_rel = 0.5 * bt * ( 1. + coef_PKA * 0.75 ); //!PKA effect (*1~1.75)
    
	let Jrel_inf=a_rel*(-ICaL)/(1.0 + Math.pow((1.5/cajsr),8) );
	if(  celltype==2 ){
		Jrel_inf=Jrel_inf*1.7;
	}
	let tau_rel = bt / ( 1.0 + 0.0123 / cajsr ) * ( 1. - coef_PKA * 0.5 ); //!PKA effect (*0.5~1)
	
	if ( tau_rel<0.001 ) {
		tau_rel=0.001; 
	}
	
	let dJrelnp=(Jrel_inf-Jrelnp)/tau_rel;
	let btp=1.25*bt;
	let a_relp = 0.5 * btp * ( 1. + coef_PKA * 0.75 ) ; //!PKA effect (*1~1.75)
	let Jrel_infp=a_relp*(-ICaL)/(1.0 + Math.pow((1.5/cajsr),8) );
	if( celltype==2 ) {
		Jrel_infp=Jrel_infp*1.7;
	}
	let tau_relp = btp / ( 1.0 + 0.0123 / cajsr ) * ( 1. - coef_PKA * 0.5 ); //!PKA effect (*0.5~1)
	
	if ( tau_relp<0.001 ) {
		tau_relp=0.001; 
	}
	
	let dJrelp=(Jrel_infp-Jrelp)/tau_relp;
	let fJrelp=(1.0/(1.0+KmCaMK/CaMKa));
	let Jrel=(1.0-fJrelp)*Jrelnp+fJrelp*Jrelp;
	
	// %calculate serca pump, ca uptake flux
	let Jupnp = 0.004375 * cai / ( cai
                                     + 0.00092 * ( 1. - coef_PKA * 0.46 ) ) ; //!PKA effect (*0.54~1)
	let Jupp = 2.75 * 0.004375 * cai / ( cai
                                           + ( ( 0.00092 - 0.00017 )
                                              * ( 1. - coef_PKA * 0.46 ) ) ); //!PKA effect (*0.54~1)
	if ( celltype==1 ) {
		Jupnp=Jupnp*1.3;
		Jupp=Jupp*1.3;
	}
	let fJupp=(1.0/(1.0+KmCaMK/CaMKa));
	let Jleak=0.0039375*cansr/15.0;
	let Jup=(1.0-fJupp)*Jupnp+fJupp*Jupp-Jleak;
	
	// %calculate tranlocation flux
	let Jtr=(cansr-cajsr)/100.0;
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	// %calcium buffer constants
	let cmdnmax=0.05;
	if ( celltype==1 ) {
		cmdnmax=cmdnmax*1.3;
	}
	let kmcmdn=0.00238;
	let trpnmax=0.07;
	let kmtrpn = 0.0005 * ( 1. + coef_PKA * 0.6 ); //!PKA effect (*1~1.6)
	let BSRmax=0.047;
	let KmBSR=0.00087;
	let BSLmax=1.124;
	let KmBSL=0.0087;
	let csqnmax=10.0;
	let kmcsqn=0.8;
	
	// %update intracellular concentrations, using buffers for cai, cass, cajsr
	let dnai=-(INa+INaL+3.0*INaCa_i+3.0*INaK+INab)*Acap/(F*vmyo)+JdiffNa*vss/vmyo;
	let dnass=-(ICaNa+3.0*INaCa_ss)*Acap/(F*vss)-JdiffNa;
	
	let dki=-(Ito+IKr+IKs+IK1+IKb+Istim-2.0*INaK)*Acap/(F*vmyo)+JdiffK*vss/vmyo;
	let dkss=-(ICaK)*Acap/(F*vss)-JdiffK;
	
	let Bcai=1.0/(1.0 + cmdnmax * kmcmdn / Math.pow((kmcmdn+cai),2) + trpnmax * kmtrpn / Math.pow((kmtrpn+cai),2) );
	let dcai=Bcai*(-(IpCa+ICab-2.0*INaCa_i)*Acap/(2.0*F*vmyo)-Jup*vnsr/vmyo+Jdiff*vss/vmyo);
	
	let Bcass=1.0/(1.0 + BSRmax * KmBSR / Math.pow((KmBSR+cass),2) + BSLmax * KmBSL / Math.pow((KmBSL+cass),2) );
	let dcass=Bcass*(-(ICaL-2.0*INaCa_ss)*Acap/(2.0*F*vss)+Jrel*vjsr/vss-Jdiff);
	
	let dcansr=Jup-Jtr*vjsr/vnsr;
	
	let Bcajsr=1.0/(1.0 + csqnmax * kmcsqn / Math.pow((kmcsqn+cajsr),2) );
	let dcajsr=Bcajsr*(Jtr-Jrel);
	
    
    if ( cass >= 0.03 ) { // cass <= 0.03 [mM]
        if ( dcass > 0 ) {
            dcass = 0;
        }
    }
    
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// %output the state venctor when ode_flag==1, and the calculated currents and fluxes otherwise
	//if flag_ode==1
	//    output=[dv dnai dnass dki dkss dcai dcass dcansr dcajsr dm dhf dhs dj dhsp djp dmL dhL dhLp da diF diS dap diFp diSp dd dff dfs dfcaf dfcas djca dnca dffp dfcafp dxrf dxrs dxs1 dxs2 dxk1 dJrelnp dJrelp dCaMKt]';
	//else
	//    output=[INa INaL Ito ICaL IKr IKs IK1 INaCa_i INaCa_ss INaK  IKb INab ICab IpCa Jdiff JdiffNa JdiffK Jup Jleak Jtr Jrel CaMKa Istim];
	//end
	//
	Xdot[0] = dv;
	Xdot[1] = dnai;
	Xdot[2] = dnass;
	Xdot[3] = dki;
	Xdot[4] = dkss;
	Xdot[5] = dcai;
	Xdot[6] = dcass;
	Xdot[7] = dcansr;
	Xdot[8] = dcajsr;
	Xdot[9] = dm;
	Xdot[10] = dhf;
	Xdot[11] = dhs;
	Xdot[12] = dj;
	Xdot[13] = dhsp;
	Xdot[14] = djp;
	Xdot[15] = dmL;
	Xdot[16] = dhL;
	Xdot[17] = dhLp;
	Xdot[18] = da;
	Xdot[19] = diF;
	Xdot[20] = diS;
	Xdot[21] = dap;
	Xdot[22] = diFp;
	Xdot[23] = diSp;
	Xdot[24] = dd;
	Xdot[25] = dff;
	Xdot[26] = dfs;
	Xdot[27] = dfcaf;
	Xdot[28] = dfcas;
	Xdot[29] = djca;
	Xdot[30] = dnca;
	Xdot[31] = dffp;
	Xdot[32] = dfcafp;
	Xdot[33] = dxrf;
	Xdot[34] = dxrs;
	Xdot[35] = dxs1;
	Xdot[36] = dxs2;
	Xdot[37] = dxk1;
	Xdot[38] = dJrelnp;
	Xdot[39] = dJrelp;
	Xdot[40] = dCaMKt;
	currents[0] = INa;
    currents[1] = ICaL;
    currents[2] = IKs;
    currents[3] = IKr;
    currents[4] = INaK;
    currents[5] = IKb;
    currents[6] = Ito;
    currents[7] = INaCa_i;
    currents[8] = INaCa_ss;
    currents[9] = IK1;
    currents[10] = INaK;
    // currents[11] = I;
    /// currents[12] = I;
	// INaL Ito ICaL IKr IKs IK1 INaCa_i INaCa_ss INaK  IKb INab ICab IpCa Jdiff JdiffNa JdiffK Jup Jleak Jtr Jrel CaMKa Istim
}
