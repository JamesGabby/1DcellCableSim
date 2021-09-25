/*
 *  O'Hara-Rudy Human model with RK method
 *  TypeScript version based on the C++ implementation created by Mao-Tsuen Jeng on 6/13/11.
 *
 *	Copyright 2011 __Colleen Clancy Lab__. All rights reserved.
 *
 *  One-D simulations
 */
var BCL = 1000;
var stimuli = 2;
var tl = 3; //length of the tissue, 165 for full length
var tw = 1; //width of the tissue
var base_dt = 0.005;
var fold = 1.0 / base_dt;
var dx = 0.01; //cm
var gj_use = 2.5; //uS
var rad = 0.0011;
var Rmyo = 150.0; //ohm*cm
var Rg = 3.14159 * Math.pow(rad, 2) / (gj_use * 1.0e-6);
var Df = (1e3 * rad) / (4.0 * (Rmyo + Rg / dx));
// Initial conditions
var v = -87.;
var nai = 7.;
var nass = nai;
var ki = 145.;
var kss = ki;
var cai = 1.0e-4;
var cass = cai;
var cansr = 1.2;
var cajsr = cansr;
var m = 0.;
var hf = 1.;
var hs = 1.;
var j = 1.;
var hsp = 1.;
var jp = 1.;
var mL = 0.;
var hL = 1.;
var hLp = 1.;
var a = 0.;
var iF = 1.;
var iS = 1.;
var ap = 0.;
var iFp = 1.;
var iSp = 1.;
var d = 0.;
var ff = 1.;
var fs = 1.;
var fcaf = 1.;
var fcas = 1.;
var jca = 1.;
var nca = 0.;
var ffp = 1.;
var fcafp = 1.;
var xrf = 0.;
var xrs = 0.;
var xs1 = 0.;
var xs2 = 0.;
var xk1 = 1.;
var Jrelnp = 0.;
var Jrelp = 0.;
var CaMKt = 0.;
var sy0 = 41;
var y0 = [v, nai, nass, ki, kss,
    cai, cass, cansr, cajsr, m,
    hf, hs, j, hsp, jp,
    mL, hL, hLp, a, iF,
    iS, ap, iFp, iSp, d,
    ff, fs, fcaf, fcas, jca,
    nca, ffp, fcafp, xrf, xrs,
    xs1, xs2, xk1, Jrelnp, Jrelp,
    CaMKt];
var p = 0;
var ll, ww;
var V, dvdt, dvdt2;
var dt;
var I_Total;
var dv;
var Vm = -86;
var cycle_length;
var I_inj;
var done = 0;
var waitTime;
var Cell = /** @class */ (function () {
    function Cell(y0) {
        this.y0 = [];
        this.currents = [];
        this.y0 = y0;
    }
    return Cell;
}());
var SimState = /** @class */ (function () {
    function SimState(cellData) {
        this.cellData = new Array();
        this.cellData = cellData;
        for (var i = 0; i <= tl; i++) {
            this.cellData[i] = [new Cell([])];
        }
    }
    return SimState;
}());
var cell = new Cell(y0);
var theState = new SimState([cell]);
var theCell = new Cell(y0);
function main() {
    console.log("New start");
    theState.counter = 0;
    for (ll = 0; ll < tl; ll++) {
        theCell = theState.cellData[ll][0];
        for (var i = 0; i < sy0; i++) {
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
        }
        else if (ll < 105) {
            theCell.celltype = 2; // %endo = 0, epi = 1, M = 2
        }
        else {
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
            if (theState.t < 0.5 && (ll < 2)) {
                theCell.I_inj = -300;
            }
            else {
                theCell.I_inj = 0.0;
            }
            I_Total = Calcu_I_Total(theCell, dt);
            theCell.v = theCell.y0[0];
        }
        for (ll = 0; ll < tl; ll++) {
            //theCell = theState.cellData[ll][0];
            if (ll > 0 && ll < (tl - 1)) {
                dv = dt * (Df * (theState.cellData[ll - 1][0].v - 2 * theCell.v + theState.cellData[ll + 1][0].v) / (dx * dx));
            }
            else if (ll == 0) {
                dv = dt * (Df * (-theCell.v + theState.cellData[1][0].v) / (dx * dx));
            }
            else if (ll == (tl - 1)) {
                dv = dt * (Df * (-theCell.v + theState.cellData[ll - 1][0].v) / (dx * dx));
            }
            theCell.y0[0] += dv;
        }
        theState.t += dt;
        theState.tt += dt;
        theState.counter += 1;
        if (theState.beat > 0 && (theState.counter % (fold * 1) == 0)) {
            // fprintf ( output, "%10.3f\t",
            // 		 theState.tt );
            //	console.log("theState.tt: "+theState.tt.toFixed());
            for (ll = 0; ll < tl; ll += 1) {
                //fprintf ( output, "%8.5f\t",theState.cellData[0][ll].y0[0] );
                console.log("theState.cellData[ll][0].y0[0]: " + theState.cellData[ll][0].y0[0]);
            }
            //fprintf ( output, "\n" );
            var Ev = 0;
            // fprintf ( output_m, "%10.3f\t",
            //          theState.tt );
            for (ll = 15; ll <= (tl - 16); ll++) {
                Ev = Ev + (theState.cellData[ll][ww].v - theState.cellData[ll + 1][ww].v) * (1 / ((tl - ll - 1) * dx + 2) - 1 / ((tl - ll) * dx + 2));
            }
            // fprintf ( output_m, "%9.7g\t", Ev );
            // fprintf ( output_m, "\n" );
        }
        if (theState.counter % (fold * 100) == 0) {
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
    for (ll = 0; ll < tl; ll++) {
        //theCell = theState.cellData[ll][0];
        //fwrite( theCell, sizeof(Cell), 1, fp);
        console.log("theCell: " + theCell);
    }
    //fclose(fp);
    console.log("Saved.");
    theState.beat = 1;
    theState.t = 0.0;
    theState.tt = 0.0;
    theState.counter = 0;
    theState.tstep = 0;
    console.log("theState: " + theState);
}
function Calcu_I_Total(theCell, dt) {
    theCell = new Cell(y0);
    var tspan = [];
    var tspan_total = [0, dt]; // ms
    var partition = dt; // ms
    var parts_total = Math.floor(0.5 + (tspan_total[1] - tspan_total[0]) / partition);
    var part;
    var sy = 1;
    var sy0 = 41;
    var y = new Array(sy0);
    var t1 = new Array(sy);
    for (var j_1 = 0; j_1 < y.length; j_1++) {
        y[j_1] = new Array(sy);
    }
    var y1 = new Array(sy);
    for (var i = 0; i < sy; i++) {
        y1[i] = y[i];
    }
    V = theCell.y0[0];
    for (part = 0; part < parts_total; part++) {
        tspan[0] = partition * part;
        tspan[1] = tspan[0] + partition;
        integrate_rk2rk1(ord_model_betaStim, tspan, sy0, theCell.y0, dt, sy, y1, t1, theCell);
        for (var i = 0; i < sy0; i++) {
            theCell.y0[i] = y1[sy - 1][i];
            // console.log(theCell.y0[i]);
        }
    }
    //console.log(-(theCell.y0[0] - V) / dt);	
    return (-(theCell.y0[0] - V) / dt);
}
//=====================================================================================================
function integrate_rk2rk1(f, tspan, sy0, y0, dt, sy, y, t1, theCell) {
    var tol = 0.001;
    var counter, counter_sy;
    var t, tn, err_y;
    var y1 = [], y2 = [], z = [], z2 = [];
    var s, sm;
    var dy1 = [], dy2 = [];
    var ydot = [];
    var runType = ['y', 'd', 'o', 't'];
    var i, j, k;
    var st = Math.floor(0.5 + (tspan[1] - tspan[0]) / (dt * sy));
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
            for (i = 0; i < sy0; i++) {
                dy2[i] = ydot[i] * dt; // k2
                z[i] = y1[i] + 0.5 * (dy1[i] + dy2[i]); // Result using RK2
            }
            sm = 1;
            for (i = 0; i < sy0; i++) {
                if (z[i] != z2[i] && z[i] == z[i]) {
                    s = Math.pow(tol * dt * (Math.abs(z[i]) + Math.abs(z2[i])) / (2 * Math.abs(z[i] - z2[i])), 0.5);
                    if (sm > s && s == s) {
                        sm = s;
                        if (1E-1 > sm) {
                            sm = 1E-1;
                        }
                    }
                    else if (s != s) { // if s is NAN   
                        sm = 1E-1;
                    }
                    //console.log( << endl;
                }
                else if (z[i] != z[i]) {
                    //console.log(z[i]);  // z[i] is NAN
                    //break; // [?]
                    sm = 1E-1;
                }
            }
            if (sm >= 1) {
                // Keep RK2
                t += dt;
                for (i = 0; i < sy0; i++) {
                    y1[i] = z[i];
                }
            }
            else {
                // Use RK2
                var invsm = 8; //1 / sm;
                var dt2 = dt / invsm; //sm;
                for (var id1 = 1; id1 <= invsm; id1++) {
                    f(y1, ydot, theCell, dt2);
                    for (i = 0; i < sy0; i++) {
                        dy1[i] = ydot[i] * dt2; // k1
                        y2[i] = y1[i] + dy1[i];
                    }
                    f(y2, ydot, theCell, dt2);
                    for (i = 0; i < sy0; i++) {
                        dy2[i] = ydot[i] * dt2; // k2
                        // Result using RK2
                        y1[i] = y1[i] + 0.5 * (dy1[i] + dy2[i]);
                    }
                    t += dt2;
                }
            }
        }
        t1[counter_sy] = t;
        for (i = 0; i < sy0; i++) {
            y[counter_sy][i] = y1[i];
        }
    }
}
//========================================================================================================
function ord_model_betaStim(X, Xdot, theCell, dt) {
    var celltype = theCell.celltype; // %endo = 0, epi = 1, M = 2
    var Istim = theCell.I_inj;
    var coef_PKA = theCell.coef_PKA;
    var currents = theCell.currents;
    X = [];
    Xdot = [];
    // %extracellular ionic concentrations
    var nao = 140.0;
    var cao = 1.8;
    var ko = 5.4;
    // %physical constants
    var R = 8314.0;
    var T = 310.0;
    var F = 96485.0;
    // %cell geometry
    var L = 0.01;
    var rad = 0.0011;
    var vcell = 1000 * 3.14 * rad * rad * L;
    var Ageo = 2 * 3.14 * rad * rad + 2 * 3.14 * rad * L;
    var Acap = 2 * Ageo;
    var vmyo = 0.68 * vcell;
    var vnsr = 0.0552 * vcell;
    var vjsr = 0.0048 * vcell;
    var vss = 0.02 * vcell;
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (X[6] > 0.03) { // cass <= 0.03 [mM]
        X[6] = 0.03;
    }
    // %give names to the state vector values
    var v = X[0];
    var nai = X[1];
    var nass = X[2];
    var ki = X[3];
    var kss = X[4];
    var cai = X[5];
    var cass = X[6];
    var cansr = X[7];
    var cajsr = X[8];
    var m = X[9];
    var hf = X[10];
    var hs = X[11];
    var j = X[12];
    var hsp = X[13];
    var jp = X[14];
    var mL = X[15];
    var hL = X[16];
    var hLp = X[17];
    var a = X[18];
    var iF = X[19];
    var iS = X[20];
    var ap = X[21];
    var iFp = X[22];
    var iSp = X[23];
    var d = X[24];
    var ff = X[25];
    var fs = X[26];
    var fcaf = X[27];
    var fcas = X[28];
    var jca = X[29];
    var nca = X[30];
    var ffp = X[31];
    var fcafp = X[32];
    var xrf = X[33];
    var xrs = X[34];
    var xs1 = X[35];
    var xs2 = X[36];
    var xk1 = X[37];
    var Jrelnp = X[38];
    var Jrelp = X[39];
    var CaMKt = X[40];
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // %CaMK constants
    var KmCaMK = 0.15;
    var aCaMK = 0.05;
    var bCaMK = 0.00068;
    var CaMKo = 0.05;
    var KmCaM = 0.0015;
    // %update CaMK
    var CaMKb = CaMKo * (1.0 - CaMKt) / (1.0 + KmCaM / cass);
    var CaMKa = CaMKb + CaMKt;
    var dCaMKt = aCaMK * CaMKb * (CaMKb + CaMKt) - bCaMK * CaMKt;
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // %reversal potentials
    var ENa = (R * T / F) * Math.log(nao / nai);
    var EK = (R * T / F) * Math.log(ko / ki);
    var PKNa = 0.01833;
    var EKs = (R * T / F) * Math.log((ko + PKNa * nao) / (ki + PKNa * nai));
    // %convenient shorthand calculations
    var vffrt = v * F * F / (R * T);
    var vfrt = v * F / (R * T);
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // %calculate INa
    var mss = 1.0 / (1.0 + Math.exp((-(v + 39.57)) / 9.871));
    var tm = 1.0 / (6.765 * Math.exp((v + 11.64) / 34.77) + 8.552 * Math.exp(-(v + 77.42) / 5.955));
    var dm = (mss - m) / tm;
    var hss = 1.0 / (1. + Math.exp((v + 82.90 + coef_PKA * 5.) / 6.086)); //!PKA effect (+0~5)
    var thf = 1.0 / (1.432e-5 * Math.exp(-(v + 1.196) / 6.285) + 6.149 * Math.exp((v + 0.5096) / 20.27));
    var ths = 1.0 / (0.009794 * Math.exp(-(v + 17.95) / 28.05) + 0.3343 * Math.exp((v + 5.730) / 56.66));
    var Ahf = 0.99;
    var Ahs = 1.0 - Ahf;
    var dhf = (hss - hf) / thf;
    var dhs = (hss - hs) / ths;
    var h = Ahf * hf + Ahs * hs;
    var jss = hss;
    var tj = 2.038 + 1.0 / (0.02136 * Math.exp(-(v + 100.6) / 8.281) + 0.3052 * Math.exp((v + 0.9941) / 38.45));
    var dj = (jss - j) / tj;
    var hssp = 1.0 / (1. + Math.exp((v + 89.1 + coef_PKA * 5.) / 6.086)); //!PKA effect (+0~5)
    var thsp = 3.0 * ths;
    var dhsp = (hssp - hsp) / thsp;
    var hp = Ahf * hf + Ahs * hsp;
    var tjp = 1.46 * tj;
    var djp = (jss - jp) / tjp;
    var GNa = 75. * (1. + coef_PKA * 1.7); //!PKA effect (*1~2.7)
    var fINap = (1.0 / (1.0 + KmCaMK / CaMKa));
    var INa = GNa * (v - ENa) * m * m * m * ((1.0 - fINap) * h * j + fINap * hp * jp);
    // %calculate INaL
    var mLss = 1.0 / (1.0 + Math.exp((-(v + 42.85)) / 5.264));
    var tmL = tm;
    var dmL = (mLss - mL) / tmL;
    var hLss = 1.0 / (1.0 + Math.exp((v + 87.61) / 7.488));
    var thL = 200.0;
    var dhL = (hLss - hL) / thL;
    var hLssp = 1.0 / (1.0 + Math.exp((v + 93.81) / 7.488));
    var thLp = 3.0 * thL;
    var dhLp = (hLssp - hLp) / thLp;
    var GNaL = 0.0075;
    if (celltype == 1) {
        GNaL = GNaL * 0.6;
    }
    var fINaLp = (1.0 / (1.0 + KmCaMK / CaMKa));
    var INaL = GNaL * (v - ENa) * mL * ((1.0 - fINaLp) * hL + fINaLp * hLp);
    // %calculate Ito
    var ass = 1.0 / (1.0 + Math.exp((-(v - 14.34)) / 14.82));
    var ta = 1.0515 / (1.0 / (1.2089 * (1.0 + Math.exp(-(v - 18.4099) / 29.3814))) + 3.5 / (1.0 + Math.exp((v + 100.0) / 29.3814)));
    var da = (ass - a) / ta;
    var iss = 1.0 / (1.0 + Math.exp((v + 43.94) / 5.711));
    var delta_epi;
    if (celltype == 1) {
        delta_epi = 1.0 - (0.95 / (1.0 + Math.exp((v + 70.0) / 5.0)));
    }
    else {
        delta_epi = 1.0;
    }
    var tiF = 4.562 + 1 / (0.3933 * Math.exp((-(v + 100.0)) / 100.0) + 0.08004 * Math.exp((v + 50.0) / 16.59));
    var tiS = 23.62 + 1 / (0.001416 * Math.exp((-(v + 96.52)) / 59.05) + 1.780e-8 * Math.exp((v + 114.1) / 8.079));
    tiF = tiF * delta_epi;
    tiS = tiS * delta_epi;
    var AiF = 1.0 / (1.0 + Math.exp((v - 213.6) / 151.2));
    var AiS = 1.0 - AiF;
    var diF = (iss - iF) / tiF;
    var diS = (iss - iS) / tiS;
    var i = AiF * iF + AiS * iS;
    var assp = 1.0 / (1.0 + Math.exp((-(v - 24.34)) / 14.82));
    var dap = (assp - ap) / ta;
    var dti_develop = 1.354 + 1.0e-4 / (Math.exp((v - 167.4) / 15.89) + Math.exp(-(v - 12.23) / 0.2154));
    var dti_recover = 1.0 - 0.5 / (1.0 + Math.exp((v + 70.0) / 20.0));
    var tiFp = dti_develop * dti_recover * tiF;
    var tiSp = dti_develop * dti_recover * tiS;
    var diFp = (iss - iFp) / tiFp;
    var diSp = (iss - iSp) / tiSp;
    var ip = AiF * iFp + AiS * iSp;
    var Gto = 0.02;
    if (celltype == 1) {
        Gto = Gto * 4.0;
    }
    else if (celltype == 2) {
        Gto = Gto * 4.0;
    }
    var fItop = (1.0 / (1.0 + KmCaMK / CaMKa));
    var Ito = Gto * (v - EK) * ((1.0 - fItop) * a * i + fItop * ap * ip);
    // %calculate ICaL, ICaNa, ICaK
    var dss = 1.0 / (1.0 + Math.exp((-(v + 3.940 + coef_PKA * 16.)) / 4.230)); //!PKA effect (+0~16)
    var td = 0.6 + 1.0 / (Math.exp(-0.05 * (v + 6.0)) + Math.exp(0.09 * (v + 14.0)));
    var dd = (dss - d) / td;
    var fss = 1.0 / (1.0 + Math.exp((v + 19.58 + coef_PKA * 8.) / 3.696)); //!PKA effect (+0~8)
    var tff = 7.0 + 1.0 / (0.0045 * Math.exp(-(v + 20.0) / 10.0) + 0.0045 * Math.exp((v + 20.0) / 10.0));
    var tfs = 1000.0 + 1.0 / (0.000035 * Math.exp(-(v + 5.0) / 4.0) + 0.000035 * Math.exp((v + 5.0) / 6.0));
    var Aff = 0.6;
    var Afs = 1.0 - Aff;
    var dff = (fss - ff) / tff;
    var dfs = (fss - fs) / tfs;
    var f = Aff * ff + Afs * fs;
    var fcass = fss;
    var tfcaf = 7.0 + 1.0 / (0.04 * Math.exp(-(v - 4.0) / 7.0) + 0.04 * Math.exp((v - 4.0) / 7.0));
    var tfcas = 100.0 + 1.0 / (0.00012 * Math.exp(-v / 3.0) + 0.00012 * Math.exp(v / 7.0));
    var Afcaf = 0.3 + 0.6 / (1.0 + Math.exp((v - 10.0) / 10.0));
    var Afcas = 1.0 - Afcaf;
    var dfcaf = (fcass - fcaf) / tfcaf;
    var dfcas = (fcass - fcas) / tfcas;
    var fca = Afcaf * fcaf + Afcas * fcas;
    var tjca = 75.0;
    var djca = (fcass - jca) / tjca;
    var tffp = 2.5 * tff;
    var dffp = (fss - ffp) / tffp;
    var fp = Aff * ffp + Afs * fs;
    var tfcafp = 2.5 * tfcaf;
    var dfcafp = (fcass - fcafp) / tfcafp;
    var fcap = Afcaf * fcafp + Afcas * fcas;
    var Kmn = 0.002;
    var k2n = 1000.0;
    var km2n = jca * 1.0;
    var anca = 1.0 / (k2n / km2n + Math.pow((1.0 + Kmn / cass), 4));
    var dnca = anca * k2n - nca * km2n;
    //!PKA effect cass upto 0.03 mM
    var PhiCaL = 4.0 * vffrt * (cass * Math.exp(2.0 * vfrt) - 0.341 * cao) / (Math.exp(2.0 * vfrt) - 1.0);
    var PhiCaNa = 1.0 * vffrt * (0.75 * nass * Math.exp(1.0 * vfrt) - 0.75 * nao) / (Math.exp(1.0 * vfrt) - 1.0);
    var PhiCaK = 1.0 * vffrt * (0.75 * kss * Math.exp(1.0 * vfrt) - 0.75 * ko) / (Math.exp(1.0 * vfrt) - 1.0);
    var zca = 2.0;
    var PCa = 0.0001 * (1. + coef_PKA * 1.5); //!PKA effect (*1~2.5)
    if (celltype == 1) {
        PCa = PCa * 1.2;
    }
    else if (celltype == 2) {
        PCa = PCa * 2.5;
    }
    var PCap = 1.1 * PCa;
    var PCaNa = 0.00125 * PCa;
    var PCaK = 3.574e-4 * PCa;
    var PCaNap = 0.00125 * PCap;
    var PCaKp = 3.574e-4 * PCap;
    var fICaLp = (1.0 / (1.0 + KmCaMK / CaMKa));
    var ICaL = (1.0 - fICaLp) * PCa * PhiCaL * d * (f * (1.0 - nca) + jca * fca * nca) + fICaLp * PCap * PhiCaL * d * (fp * (1.0 - nca) + jca * fcap * nca);
    var ICaNa = (1.0 - fICaLp) * PCaNa * PhiCaNa * d * (f * (1.0 - nca) + jca * fca * nca) + fICaLp * PCaNap * PhiCaNa * d * (fp * (1.0 - nca) + jca * fcap * nca);
    var ICaK = (1.0 - fICaLp) * PCaK * PhiCaK * d * (f * (1.0 - nca) + jca * fca * nca) + fICaLp * PCaKp * PhiCaK * d * (fp * (1.0 - nca) + jca * fcap * nca);
    // %calculate IKr
    var xrss = 1.0 / (1.0 + Math.exp((-(v + 8.337)) / 6.789));
    var txrf = 12.98 + 1.0 / (0.3652 * Math.exp((v - 31.66) / 3.869) + 4.123e-5 * Math.exp((-(v - 47.78)) / 20.38));
    var txrs = 1.865 + 1.0 / (0.06629 * Math.exp((v - 34.70) / 7.355) + 1.128e-5 * Math.exp((-(v - 29.74)) / 25.94));
    var Axrf = 1.0 / (1.0 + Math.exp((v + 54.81) / 38.21));
    var Axrs = 1.0 - Axrf;
    var dxrf = (xrss - xrf) / txrf;
    var dxrs = (xrss - xrs) / txrs;
    var xr = Axrf * xrf + Axrs * xrs;
    var rkr = 1.0 / (1.0 + Math.exp((v + 55.0) / 75.0)) * 1.0 / (1.0 + Math.exp((v - 10.0) / 30.0));
    var GKr = 0.046;
    if (celltype == 1) {
        GKr = GKr * 1.3;
    }
    else if (celltype == 2) {
        GKr = GKr * 0.8;
    }
    var IKr = GKr * Math.sqrt(ko / 5.4) * xr * rkr * (v - EK);
    // %calculate IKs
    var xs1ss = 1.0 / (1.0 + Math.exp((-(v + 11.60)) / 8.932)) * (1. + coef_PKA * 0.4); //!PKA effect (*1~1.2)
    var txs1 = 817.3 + 1.0 / (2.326e-4 * Math.exp((v + 48.28) / 17.80) + 0.001292 * Math.exp((-(v + 210.0)) / 230.0));
    var dxs1 = (xs1ss - xs1) / txs1;
    var xs2ss = xs1ss;
    var txs2 = 1.0 / (0.01 * Math.exp((v - 50.0) / 20.0) + 0.0193 * Math.exp((-(v + 66.54)) / 31.0));
    var dxs2 = (xs2ss - xs2) / txs2;
    var KsCa = 1.0 + 0.6 / (1.0 + Math.pow((3.8e-5 / cai), 1.4));
    var GKs = 0.0034 * (1. + coef_PKA * 2.2); //!PKA effect (*1~3.2)
    if (celltype == 1) {
        GKs = GKs * 1.4;
    }
    var IKs = GKs * KsCa * xs1 * xs2 * (v - EKs);
    var xk1ss = 1.0 / (1.0 + Math.exp(-(v + 2.5538 * ko + 144.59) / (1.5692 * ko + 3.8115)));
    var txk1 = 122.2 / (Math.exp((-(v + 127.2)) / 20.36) + Math.exp((v + 236.8) / 69.33));
    var dxk1 = (xk1ss - xk1) / txk1;
    var rk1 = 1.0 / (1.0 + Math.exp((v + 105.8 - 2.6 * ko) / 9.493));
    var GK1 = 0.1908;
    if (celltype == 1) {
        GK1 = GK1 * 1.2;
    }
    else if (celltype == 2) {
        GK1 = GK1 * 1.3;
    }
    var IK1 = GK1 * Math.sqrt(ko) * rk1 * xk1 * (v - EK);
    // %calculate INaCa_i
    var kna1 = 15.0;
    var kna2 = 5.0;
    var kna3 = 88.12;
    var kasymm = 12.5;
    var wna = 6.0e4;
    var wca = 6.0e4;
    var wnaca = 5.0e3;
    var kcaon = 1.5e6;
    var kcaoff = 5.0e3;
    var qna = 0.5224;
    var qca = 0.1670;
    var hca = Math.exp((qca * v * F) / (R * T));
    var hna = Math.exp((qna * v * F) / (R * T));
    var h1 = 1 + nai / kna3 * (1 + hna);
    var h2 = (nai * hna) / (kna3 * h1);
    var h3 = 1.0 / h1;
    var h4 = 1.0 + nai / kna1 * (1 + nai / kna2);
    var h5 = nai * nai / (h4 * kna1 * kna2);
    var h6 = 1.0 / h4;
    var h7 = 1.0 + nao / kna3 * (1.0 + 1.0 / hna);
    var h8 = nao / (kna3 * hna * h7);
    var h9 = 1.0 / h7;
    var h10 = kasymm + 1.0 + nao / kna1 * (1.0 + nao / kna2);
    var h11 = nao * nao / (h10 * kna1 * kna2);
    var h12 = 1.0 / h10;
    var k1 = h12 * cao * kcaon;
    var k2 = kcaoff;
    var k3p = h9 * wca;
    var k3pp = h8 * wnaca;
    var k3 = k3p + k3pp;
    var k4p = h3 * wca / hca;
    var k4pp = h2 * wnaca;
    var k4 = k4p + k4pp;
    var k5 = kcaoff;
    var k6 = h6 * cai * kcaon;
    var k7 = h5 * h2 * wna;
    var k8 = h8 * h11 * wna;
    var x1 = k2 * k4 * (k7 + k6) + k5 * k7 * (k2 + k3);
    var x2 = k1 * k7 * (k4 + k5) + k4 * k6 * (k1 + k8);
    var x3 = k1 * k3 * (k7 + k6) + k8 * k6 * (k2 + k3);
    var x4 = k2 * k8 * (k4 + k5) + k3 * k5 * (k1 + k8);
    var E1 = x1 / (x1 + x2 + x3 + x4);
    var E2 = x2 / (x1 + x2 + x3 + x4);
    var E3 = x3 / (x1 + x2 + x3 + x4);
    var E4 = x4 / (x1 + x2 + x3 + x4);
    var KmCaAct = 150.0e-6;
    var allo = 1.0 / (1.0 + Math.pow((KmCaAct / cai), 2));
    var zna = 1.0;
    var JncxNa = 3.0 * (E4 * k7 - E1 * k8) + E3 * k4pp - E2 * k3pp;
    var JncxCa = E2 * k2 - E1 * k1;
    var Gncx = 0.0008;
    if (celltype == 1) {
        Gncx = Gncx * 1.1;
    }
    else if (celltype == 2) {
        Gncx = Gncx * 1.4;
    }
    var INaCa_i = 0.8 * Gncx * allo * (zna * JncxNa + zca * JncxCa);
    // %calculate INaCa_ss
    h1 = 1 + nass / kna3 * (1 + hna);
    h2 = (nass * hna) / (kna3 * h1);
    h3 = 1.0 / h1;
    h4 = 1.0 + nass / kna1 * (1 + nass / kna2);
    h5 = nass * nass / (h4 * kna1 * kna2);
    h6 = 1.0 / h4;
    h7 = 1.0 + nao / kna3 * (1.0 + 1.0 / hna);
    h8 = nao / (kna3 * hna * h7);
    h9 = 1.0 / h7;
    h10 = kasymm + 1.0 + nao / kna1 * (1 + nao / kna2);
    h11 = nao * nao / (h10 * kna1 * kna2);
    h12 = 1.0 / h10;
    k1 = h12 * cao * kcaon;
    k2 = kcaoff;
    k3p = h9 * wca;
    k3pp = h8 * wnaca;
    k3 = k3p + k3pp;
    k4p = h3 * wca / hca;
    k4pp = h2 * wnaca;
    k4 = k4p + k4pp;
    k5 = kcaoff;
    k6 = h6 * cass * kcaon;
    k7 = h5 * h2 * wna;
    k8 = h8 * h11 * wna;
    x1 = k2 * k4 * (k7 + k6) + k5 * k7 * (k2 + k3);
    x2 = k1 * k7 * (k4 + k5) + k4 * k6 * (k1 + k8);
    x3 = k1 * k3 * (k7 + k6) + k8 * k6 * (k2 + k3);
    x4 = k2 * k8 * (k4 + k5) + k3 * k5 * (k1 + k8);
    E1 = x1 / (x1 + x2 + x3 + x4);
    E2 = x2 / (x1 + x2 + x3 + x4);
    E3 = x3 / (x1 + x2 + x3 + x4);
    E4 = x4 / (x1 + x2 + x3 + x4);
    KmCaAct = 150.0e-6;
    allo = 1.0 / (1.0 + Math.pow((KmCaAct / cass), 2));
    JncxNa = 3.0 * (E4 * k7 - E1 * k8) + E3 * k4pp - E2 * k3pp;
    JncxCa = E2 * k2 - E1 * k1;
    var INaCa_ss = 0.2 * Gncx * allo * (zna * JncxNa + zca * JncxCa);
    // %calculate INaK
    var k1p = 949.5;
    var k1m = 182.4;
    var k2p = 687.2;
    var k2m = 39.4;
    k3p = 1899.0;
    var k3m = 79300.0;
    k4p = 639.0;
    var k4m = 40.0;
    var Knai0 = 9.073 * (1. - coef_PKA * 0.3); //!PKA effect (*0.7~1);
    var Knao0 = 27.78;
    var delta = -0.1550;
    var Knai = Knai0 * Math.exp((delta * v * F) / (3.0 * R * T));
    var Knao = Knao0 * Math.exp(((1.0 - delta) * v * F) / (3.0 * R * T));
    var Kki = 0.5;
    var Kko = 0.3582;
    var MgADP = 0.05;
    var MgATP = 9.8;
    var Kmgatp = 1.698e-7;
    var H = 1.0e-7;
    var eP = 4.2;
    var Khp = 1.698e-7;
    var Knap = 224.0;
    var Kxkur = 292.0;
    var P = eP / (1.0 + H / Khp + nai / Knap + ki / Kxkur);
    var a1 = (k1p * Math.pow((nai / Knai), 3)) / (Math.pow((1.0 + nai / Knai), 3) + Math.pow((1.0 + ki / Kki), 2) - 1.0);
    var b1 = k1m * MgADP;
    var a2 = k2p;
    var b2 = (k2m * Math.pow((nao / Knao), 3)) / (Math.pow((1.0 + nao / Knao), 3) + Math.pow((1.0 + ko / Kko), 2) - 1.0);
    var a3 = (k3p * Math.pow((ko / Kko), 2)) / (Math.pow((1.0 + nao / Knao), 3) + Math.pow((1.0 + ko / Kko), 2) - 1.0);
    var b3 = (k3m * P * H) / (1.0 + MgATP / Kmgatp);
    var a4 = (k4p * MgATP / Kmgatp) / (1.0 + MgATP / Kmgatp);
    var b4 = (k4m * Math.pow((ki / Kki), 2)) / (Math.pow((1.0 + nai / Knai), 3) + Math.pow((1.0 + ki / Kki), 2) - 1.0);
    x1 = a4 * a1 * a2 + b2 * b4 * b3 + a2 * b4 * b3 + b3 * a1 * a2;
    x2 = b2 * b1 * b4 + a1 * a2 * a3 + a3 * b1 * b4 + a2 * a3 * b4;
    x3 = a2 * a3 * a4 + b3 * b2 * b1 + b2 * b1 * a4 + a3 * a4 * b1;
    x4 = b4 * b3 * b2 + a3 * a4 * a1 + b2 * a4 * a1 + b3 * b2 * a1;
    E1 = x1 / (x1 + x2 + x3 + x4);
    E2 = x2 / (x1 + x2 + x3 + x4);
    E3 = x3 / (x1 + x2 + x3 + x4);
    E4 = x4 / (x1 + x2 + x3 + x4);
    var zk = 1.0;
    var JnakNa = 3.0 * (E1 * a3 - E2 * b3);
    var JnakK = 2.0 * (E4 * b1 - E3 * a1);
    var Pnak = 30;
    if (celltype == 1) {
        Pnak = Pnak * 0.9;
    }
    else if (celltype == 2) {
        Pnak = Pnak * 0.7;
    }
    var INaK = Pnak * (zna * JnakNa + zk * JnakK);
    // %calculate IKb
    var xkb = 1.0 / (1.0 + Math.exp(-(v - 14.48) / 18.34));
    var GKb = 0.003 * (1. + coef_PKA * 1.5); //!PKA effect cass  (*1~2.5)
    if (celltype == 1) {
        GKb = GKb * 0.6;
    }
    var IKb = GKb * xkb * (v - EK);
    // %calculate INab
    var PNab = 3.75e-10;
    var INab = PNab * vffrt * (nai * Math.exp(vfrt) - nao) / (Math.exp(vfrt) - 1.0);
    // %calculate ICab
    var PCab = 2.5e-8;
    var ICab = PCab * 4.0 * vffrt * (cai * Math.exp(2.0 * vfrt) - 0.341 * cao) / (Math.exp(2.0 * vfrt) - 1.0);
    // %calculate IpCa
    var GpCa = 0.0005;
    var IpCa = GpCa * cai / (0.0005 + cai);
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // %udtate the membrane voltage
    var dv = -(INa + INaL + Ito + ICaL + ICaNa + ICaK + IKr + IKs + IK1 + INaCa_i + INaCa_ss + INaK + INab + IKb + IpCa + ICab + Istim);
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // %calculate diffusion fluxes
    var JdiffNa = (nass - nai) / 2.0;
    var JdiffK = (kss - ki) / 2.0;
    var Jdiff = (cass - cai) / 0.2;
    // %calculate ryanodione receptor calcium induced calcium release from the jsr
    var bt = 4.75;
    var a_rel = 0.5 * bt * (1. + coef_PKA * 0.75); //!PKA effect (*1~1.75)
    var Jrel_inf = a_rel * (-ICaL) / (1.0 + Math.pow((1.5 / cajsr), 8));
    if (celltype == 2) {
        Jrel_inf = Jrel_inf * 1.7;
    }
    var tau_rel = bt / (1.0 + 0.0123 / cajsr) * (1. - coef_PKA * 0.5); //!PKA effect (*0.5~1)
    if (tau_rel < 0.001) {
        tau_rel = 0.001;
    }
    var dJrelnp = (Jrel_inf - Jrelnp) / tau_rel;
    var btp = 1.25 * bt;
    var a_relp = 0.5 * btp * (1. + coef_PKA * 0.75); //!PKA effect (*1~1.75)
    var Jrel_infp = a_relp * (-ICaL) / (1.0 + Math.pow((1.5 / cajsr), 8));
    if (celltype == 2) {
        Jrel_infp = Jrel_infp * 1.7;
    }
    var tau_relp = btp / (1.0 + 0.0123 / cajsr) * (1. - coef_PKA * 0.5); //!PKA effect (*0.5~1)
    if (tau_relp < 0.001) {
        tau_relp = 0.001;
    }
    var dJrelp = (Jrel_infp - Jrelp) / tau_relp;
    var fJrelp = (1.0 / (1.0 + KmCaMK / CaMKa));
    var Jrel = (1.0 - fJrelp) * Jrelnp + fJrelp * Jrelp;
    // %calculate serca pump, ca uptake flux
    var Jupnp = 0.004375 * cai / (cai
        + 0.00092 * (1. - coef_PKA * 0.46)); //!PKA effect (*0.54~1)
    var Jupp = 2.75 * 0.004375 * cai / (cai
        + ((0.00092 - 0.00017)
            * (1. - coef_PKA * 0.46))); //!PKA effect (*0.54~1)
    if (celltype == 1) {
        Jupnp = Jupnp * 1.3;
        Jupp = Jupp * 1.3;
    }
    var fJupp = (1.0 / (1.0 + KmCaMK / CaMKa));
    var Jleak = 0.0039375 * cansr / 15.0;
    var Jup = (1.0 - fJupp) * Jupnp + fJupp * Jupp - Jleak;
    // %calculate tranlocation flux
    var Jtr = (cansr - cajsr) / 100.0;
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // %calcium buffer constants
    var cmdnmax = 0.05;
    if (celltype == 1) {
        cmdnmax = cmdnmax * 1.3;
    }
    var kmcmdn = 0.00238;
    var trpnmax = 0.07;
    var kmtrpn = 0.0005 * (1. + coef_PKA * 0.6); //!PKA effect (*1~1.6)
    var BSRmax = 0.047;
    var KmBSR = 0.00087;
    var BSLmax = 1.124;
    var KmBSL = 0.0087;
    var csqnmax = 10.0;
    var kmcsqn = 0.8;
    // %update intracellular concentrations, using buffers for cai, cass, cajsr
    var dnai = -(INa + INaL + 3.0 * INaCa_i + 3.0 * INaK + INab) * Acap / (F * vmyo) + JdiffNa * vss / vmyo;
    var dnass = -(ICaNa + 3.0 * INaCa_ss) * Acap / (F * vss) - JdiffNa;
    var dki = -(Ito + IKr + IKs + IK1 + IKb + Istim - 2.0 * INaK) * Acap / (F * vmyo) + JdiffK * vss / vmyo;
    var dkss = -(ICaK) * Acap / (F * vss) - JdiffK;
    var Bcai = 1.0 / (1.0 + cmdnmax * kmcmdn / Math.pow((kmcmdn + cai), 2) + trpnmax * kmtrpn / Math.pow((kmtrpn + cai), 2));
    var dcai = Bcai * (-(IpCa + ICab - 2.0 * INaCa_i) * Acap / (2.0 * F * vmyo) - Jup * vnsr / vmyo + Jdiff * vss / vmyo);
    var Bcass = 1.0 / (1.0 + BSRmax * KmBSR / Math.pow((KmBSR + cass), 2) + BSLmax * KmBSL / Math.pow((KmBSL + cass), 2));
    var dcass = Bcass * (-(ICaL - 2.0 * INaCa_ss) * Acap / (2.0 * F * vss) + Jrel * vjsr / vss - Jdiff);
    var dcansr = Jup - Jtr * vjsr / vnsr;
    var Bcajsr = 1.0 / (1.0 + csqnmax * kmcsqn / Math.pow((kmcsqn + cajsr), 2));
    var dcajsr = Bcajsr * (Jtr - Jrel);
    if (cass >= 0.03) { // cass <= 0.03 [mM]
        if (dcass > 0) {
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
