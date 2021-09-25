Source code contains: main_cell_ss_1D.cpp (main program), ord_model_betaStim.h (O’Hara-Rudy Human model) and integrate_rk2rk1.h (Numerical method)

######################################################################################

model_INPUTSteadyStateControl.1DInitials were arrived at by starting initial conditions in the main program (line 77 to 118 and line 182 to 212) and allowing the transmural cells to run for 1000 beats (dt = 0.005, I_inj = -300 at the indicated frequency without beta-stimulation.

Figures of individual AP and ECG were generated using model_INPUTSteadyStateControl.Initials, pacing for 2 beats without beta-stimulation.

######################################################################################
Transmural cable: Cells 1–60 were endo, 61–105 were M, and 106– 165 were epi.

           
######################################################################################


Key Parameters in main_cell_ss_1D.cpp:

1. INPUT initial steady-state condition: Change initials input name at line 153

2. OUTPUT file names: Change OUTPUT file names at line 166 to line 167

3. Set basic cycle length (double BCL) at line 24
 
4. Set number of beat (double stimuli) at line 26

5. Set length of the cable (double tl) at line 28

6. Set cell size (double dx) at lind 34

7. Set gap junction (double gj_use) at line 37

8. Set cell type (theCell->celltype) at line 202 to 212 

9. Set beta-stimulation (theCell->coef_PKA) at line 213 

NOTE: Beta-stimulation inputs for linear change: range from 0 to 1 (0: basal level, 1: maximum level)

For example, set theCell->coef_PKA = 0 (no beta-stimulation) for control case.

10. Set inject current ( theCell->I_inj) at line 249



