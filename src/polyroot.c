#include "polyroot.h"

/*
C     ALGORITHM 419 COLLECTED ALGORITHMS FROM ACM.
C     ALGORITHM APPEARED IN COMM. ACM, VOL. 15, NO. 02,
C     P. 097.
 Used under ACM Software License Agreement for Non-commercial Use
 Read it in Full  - http://www.acm.org/publications/policies/softwarecrnotice
*/

static void mcon(double *ETA, double *INFINY, double *SMALNO, double *BASE) {
	*ETA = macheps();
	*INFINY = DBL_MAX;
	*SMALNO = DBL_MIN;
	*BASE = FLT_RADIX;
}

static double CMOD(double *R, double *I) {
	double cmod;
	double AR, AI,temp;
	AR = fabs(*R);
	AI = fabs(*I);

	if (AR < AI) {
		temp = AR / AI;
		cmod = AI * sqrt(1.0 + temp*temp);
		return cmod;
	}

	if (AR > AI) {
		temp = AI / AR;
		cmod = AR * sqrt(1.0 + temp*temp);
		return cmod;
	}

	cmod = AR * sqrt(2.0);

	return cmod;
}

static void CDIVID(double *AR, double *AI, double *BR, double *BI, double INFINY, double *CR, double *CI) {
	// C COMPLEX DIVISION C = A/B, AVOIDING OVERFLOW.
	double INFIN,R,D;
	if ( !(*BR != 0.0 || *BI != 0.0) ) {
		INFIN = INFINY;
		*CR = *CI = INFIN;
	}
	else {
		if (fabs(*BR) < fabs(*BI)) {
			R = *BR / *BI;
			D = *BI + R * *BR;
			*CR = (*AR*R + *AI) / D;
			*CI = (*AI*R - *AR) / D;
		}
		else {
			R = *BI / *BR;
			D = *BR + R * *BI;
			*CR = (*AR + *AI*R) / D;
			*CI = (*AI - *AR*R) / D;
		}
	}

}

static double RESCALE(int NN, double *PT, double ETA, double INFIN, double SMALNO, double BASE) {
	/*
C RETURNS A SCALE FACTOR TO MULTIPLY THE COEFFICIENTS OF THE
C POLYNOMIAL. THE SCALING IS DONE TO AVOID OVERFLOW AND TO AVOID
C UNDETECTED UNDERFLOW INTERFERING WITH THE CONVERGENCE
C CRITERION.  THE FACTOR IS A POWER OF THE BASE.
C PT - MODULUS OF COEFFICIENTS OF P
C ETA,INFIN,SMALNO,BASE - CONSTANTS DESCRIBING THE
C FLOATING POINT ARITHMETIC.*/
	double HI, LO, MAX, MIN, X, rescale,SC,L;
	int i;

	HI = sqrt(INFIN);
	LO = SMALNO / ETA;
	MAX = 0.0;
	MIN = INFIN;

	for (i = 0; i < NN; ++i) {
		X = PT[i];
		if (X > MAX) {
			MAX = X;
		}

		if (X != 0.0 && X < MIN) {
			MIN = X;
		}
	}
	// C SCALE ONLY IF THERE ARE VERY LARGE OR VERY SMALL COMPONENTS.

	rescale = 1.0;

	if (MIN >= LO && MAX <= HI) {
		return rescale;
	}
	X = LO / MIN;
	if (X <= 1.0) {
		SC = 1.0 / (sqrt(MAX)*sqrt(MIN));
	}
	else {
		SC = X;

		if (INFIN / SC > MAX) {
			SC = 1.0;
		}
	}

	L = log(SC) / log(BASE) + .500;
	rescale = pow(BASE, L);
	return rescale;
}

static double CAUCHY(int NN, double *PT, double *Q) {
	/*
	C CAUCHY COMPUTES A LOWER BOUND ON THE MODULI OF THE ZEROS OF A
	C POLYNOMIAL - PT IS THE MODULUS OF THE COEFFICIENTS.
	*/
	double cauchy, X, XM, F, DX, DF;
	int i,N;
	PT[NN - 1] = -PT[NN - 1];
	N = NN - 1;
	X = exp((log(-PT[NN - 1]) - log(PT[0])) / (double) N);

	if (PT[N - 1] != 0.0) {
		XM = -PT[NN - 1] / PT[N-1];
		if (XM < X) {
			X = XM;
		}
	}

	XM = 0.1 * X;
    F = PT[0];
    int max_iterations = 1000; // Maximum number of iterations to prevent infinite loop
    int iterations = 0; // Initialize the iteration counter
    while (F > 0.0 && iterations < max_iterations) {
        F = PT[0];
        for (i = 1; i < NN; ++i) {
            F = F*XM + PT[i];
        }
        X = XM;
        XM *= 0.1; // Reduce XM for the next iteration
        iterations++; // Increment the iteration counter
    }

	DX = X;
	iterations = 0; // Initialize the iteration counter
	while (fabs(DX / X) > 0.005 && iterations < max_iterations) {
		Q[0] = PT[0];
		for (i = 1; i < NN; ++i) {
			Q[i] = Q[i - 1]*X + PT[i];
		}
		F = Q[NN - 1];
		DF = Q[0];
		for (i = 1; i < N; ++i) {
			DF = DF*X + Q[i];
		}
		DX = F / DF;
		X = X - DX;
        iterations++; // Increment the iteration counter
	}

	cauchy = X;

	return cauchy;
}

static double ERREV(int NN, double *QR, double *QI, double MS, double MP, double ARE, double MRE) {
	/*
	C BOUNDS THE ERROR IN EVALUATING THE POLYNOMIAL BY THE HORNER
	C RECURRENCE.
	C QR,QI - THE PARTIAL SUMS
	C MS    -MODULUS OF THE POINT
	C MP    -MODULUS OF POLYNOMIAL VALUE
	C ARE, MRE -ERROR BOUNDS ON COMPLEX ADDITION AND MULTIPLICATION
	*/
	double errev, E;
	int i;

	E = CMOD(&QR[0], &QI[0]) * MRE / (ARE + MRE);

	for (i = 0; i < NN; ++i) {
		E = E*MS + CMOD(&QR[i], &QI[i]);
	}
	errev = E*(ARE + MRE) - MP*MRE;

	return errev;
}

static void POLYEV(int NN, double *SR, double *SI,double *PR, double *PI, double *QR, double *QI, double *PVR, double *PVI) {
	double T;
	int i;
	QR[0] = PR[0];
	QI[0] = PI[0];
	*PVR = QR[0];
	*PVI = QI[0];

	for (i = 1; i < NN; ++i) {
		T = *PVR * *SR - *PVI * *SI + PR[i];
		*PVI = *PVR * *SI + *PVI * *SR + PI[i];
		*PVR = T;
		QR[i] = *PVR;
		QI[i] = *PVI;
	}
}

/*      COMMON/GLOBAL/PR,PI,HR,HI,QPR,QPI,QHR,QHI,SHR,SHI,
     *    SR,SI,TR,TI,PVR,PVI,ARE,MRE,ETA,INFIN,NN
	 
*/

static void NEXTH(int *BOOL, int NN, double *QHR, double *QHI, double *QPR, double *QPI, double *TR, double *TI,double *HR, double *HI) {
	/*
	C CALCULATES THE NEXT SHIFTED H POLYNOMIAL.
	C BOOL   -  INT , IF 1 H(S) IS ESSENTIALLY ZERO
	*/

	double T1, T2;
	int N, NM1,j;
	N = NN - 1;
	NM1 = N - 1;

	if (*BOOL == 0) {
		for (j = 1; j < N; ++j) {
			T1 = QHR[j - 1];
			T2 = QHI[j - 1];
			HR[j] = *TR*T1 - *TI*T2 + QPR[j];
			HI[j] = *TR*T2 + *TI*T1 + QPI[j];
		}
		HR[0] = QPR[0];
		HI[0] = QPI[0];
	}
	else {
		for (j = 1; j < N; ++j) {
			HR[j] = QHR[j - 1];
			HI[j] = QHI[j - 1];
		}
		HR[0] = 0.0;
		HI[0] = 0.0;
	}
}

static void CALCT(int *BOOL, int NN, double *SR, double *SI, double *HR, double *HI,double *QHR, double *QHI,double ARE, double INFINY,double *PVR, double *PVI,double *TR, double *TI ) {
	/*
	C COMPUTES  T = -P(S)/H(S).
	C BOOL   - LOGICAL, SET TRUE IF H(S) IS ESSENTIALLY ZERO.
	*/
	double HVR, HVI, pvrt,pvit;
	int N;
	N = NN - 1;
	POLYEV(N, SR, SI, HR, HI, QHR, QHI, &HVR, &HVI);

	if (CMOD(&HVR, &HVI) <= ARE*10.0*CMOD(&HR[N - 1], &HI[N - 1])) {
		*BOOL = 1;
	}
	else {
		*BOOL = 0;
	}
	pvrt = -1.0 * *PVR;
	pvit = -1.0 * *PVI;
	if (*BOOL == 0) {
		CDIVID(&pvrt, &pvit, &HVR, &HVI, INFINY,TR, TI);
	}
	else {
		*TR = 0.0;
		*TI = 0.0;
	}

}

static void VRSHFT(int NN, int L3, double *ZR, double *ZI, int *CONV,double *SR, double *SI,double *PR, double *PI, double *QPR, double *QPI, double *PVR, double *PVI,
	double ARE, double MRE, double ETA, double INFIN, double *HR, double *HI, double *QHR, double *QHI, double *SHR, double *SHI, double *TR, double *TI) {
	/*
	C CARRIES OUT THE THIRD STAGE ITERATION.
	C L3 - LIMIT OF STEPS IN STAGE 3.
	C ZR,ZI   - ON ENTRY CONTAINS THE INITIAL ITERATE, IF THE
	C ITERATION CONVERGES IT CONTAINS THE FINAL ITERATE
	C ON EXIT.
	C CONV    -  .TRUE. IF ITERATION CONVERGES
	*/
	double MP, MS, OMP, RELSTP, R1, R2, TP;
	int B,BOOL,i,j;
	*SR = *ZR;
	*SI = *ZI;
	B = 0;
	*CONV = 0;


	for (i = 0; i < L3; ++i) {
		POLYEV(NN, SR, SI, PR, PI, QPR, QPI, PVR, PVI);
		MP = CMOD(PVR,PVI);
		MS = CMOD(SR, SI);
		if (MP <= 20.0*ERREV(NN, QPR, QPI, MS, MP, ARE, MRE)) {
			*CONV = 1;
			*ZR = *SR;
			*ZI = *SI;
			break;
		}
		else {//10
			if (i == 0) {
				OMP = MP;
			}
			else {
				if (!(B == 1 || MP < OMP || RELSTP >= .05)) {
					TP = RELSTP;
					B = 1;
					if (RELSTP < ETA) {
						TP = ETA;
					}
					R1 = sqrt(TP);
					R2 = *SR * (1.0 + R1) - *SI*R1;
					*SI = *SR*R1 + *SI*(1.0 + R1);
					*SR = R2;
					POLYEV(NN, SR, SI, PR, PI, QPR, QPI, PVR, PVI);
					for (j = 0; j < 5;++j) {
						CALCT(&BOOL,NN,SR,SI,HR,HI,QHR, QHI,ARE,INFIN,PVR, PVI,TR,TI);
						NEXTH(&BOOL,NN,QHR,QHI,QPR,QPI,TR, TI, HR,HI);
					}
					OMP = INFIN;
				}
				else {
					if (MP*0.1 > OMP) {
						break;
					}
					OMP = MP;
				}

			}

			CALCT(&BOOL, NN, SR, SI, HR, HI, QHR, QHI, ARE, INFIN, PVR, PVI, TR, TI);
			NEXTH(&BOOL, NN, QHR, QHI, QPR, QPI, TR, TI, HR, HI);
			CALCT(&BOOL, NN, SR, SI, HR, HI, QHR, QHI, ARE, INFIN, PVR, PVI, TR, TI);

			if (BOOL == 0) {
				RELSTP = CMOD(TR, TI) / CMOD(SR, SI);
				*SR = *SR + *TR;
				*SI = *SI + *TI;
			}

		}
	}//60
}

static void FXSHFT(int NN, int L2, double *ZR, double *ZI, int *CONV, double *SR, double *SI, double *PR, double *PI, double *QPR, double *QPI, double *PVR, double *PVI,
	double ARE, double MRE, double ETA, double INFIN, double *HR, double *HI, double *QHR, double *QHI, double *SHR, double *SHI, double *TR, double *TI) {
	/* 
	C COMPUTES L2 FIXED-SHIFT H POLYNOMIALS AND TESTS FOR
	C CONVERGENCE.
	C INITIATES A VARIABLE-SHIFT ITERATION AND RETURNS WITH THE
	C APPROXIMATE ZERO IF SUCCESSFUL.
	C L2 - LIMIT OF FIXED SHIFT STEPS
	C ZR,ZI - APPROXIMATE ZERO IF CONV IS .TRUE.
	C CONV  - LOGICAL INDICATING CONVERGENCE OF STAGE 3 ITERATION
	*/
	double OTR, OTI, SVSR, SVSI;
	int N, TEST, PASD, BOOL,j,i;
	double tp1, tp2;
	N = NN - 1;
	
	POLYEV(NN, SR, SI, PR, PI, QPR, QPI, PVR, PVI);
	TEST = 1;
	PASD = 0;
	CALCT(&BOOL, NN, SR, SI, HR, HI, QHR, QHI, ARE, INFIN, PVR, PVI, TR, TI);
	for (j = 0; j < L2; ++j) {
		OTR = *TR;
		OTI = *TI;
		NEXTH(&BOOL, NN, QHR, QHI, QPR, QPI, TR, TI, HR, HI);
		CALCT(&BOOL, NN, SR, SI, HR, HI, QHR, QHI, ARE, INFIN, PVR, PVI, TR, TI);
		//printf("\n S %g %g", *SR, *SI);
		*ZR = *SR + *TR;
		*ZI = *SI + *TI;
		if (!(BOOL == 1 || TEST == 0 || j == L2 - 1)) {
			tp1 = *TR - OTR;
			tp2 = *TI - OTI;
			if (CMOD(&tp1, &tp2) < 0.5*CMOD(ZR, ZI)) {
				if (PASD == 1) {
					for (i = 0; i < N; ++i) {
						SHR[i] = HR[i];
						SHI[i] = HI[i];
					}
					SVSR = *SR;
					SVSI = *SI;
					VRSHFT(NN, 10, ZR, ZI, CONV, SR, SI, PR, PI, QPR, QPI, PVR, PVI, ARE, MRE, ETA, INFIN, HR, HI, QHR, QHI, SHR, SHI, TR, TI);
					if (*CONV == 1) {
						break; //C THE ITERATION FAILED TO CONVERGE. TURN OFF TESTING AND RESTORE
					}
					TEST = 0;
					for (i = 0; i < N; ++i) {
						HR[i] = SHR[i];
						HI[i] = SHI[i];
					}
					*SR = SVSR;
					*SI = SVSI;
					POLYEV(NN, SR, SI, PR, PI, QPR, QPI, PVR, PVI);
					CALCT(&BOOL, NN, SR, SI, HR, HI, QHR, QHI, ARE, INFIN, PVR, PVI, TR, TI);
				}
				else {
					PASD = 1;
				}

			} else {
				PASD = 0;
			}
		}
	}

	VRSHFT(NN, 10, ZR, ZI, CONV, SR, SI, PR, PI, QPR, QPI, PVR, PVI, ARE, MRE, ETA, INFIN, HR, HI, QHR, QHI, SHR, SHI, TR, TI);
}

static void NOSHFT(int NN,int L1, double *SR, double *SI, double *PR, double *PI, double *QPR, double *QPI, double *PVR, double *PVI,
	double ARE, double MRE, double ETA, double INFIN, double *HR, double *HI, double *QHR, double *QHI, double *SHR, double *SHI, double *TR, double *TI) {
	/*
	C COMPUTES  THE DERIVATIVE  POLYNOMIAL AS THE INITIAL H
	C POLYNOMIAL AND COMPUTES L1 NO-SHIFT H POLYNOMIALS.
	*/
	double XNI, T1, T2;
	int N, NM1, i, i1,j, jj,j1;
	double prt, pit;

	N = NN - 1;
	NM1 = N - 1;

	for (i = 1; i <= N; ++i) {
		i1 = i - 1;
		XNI = (double)(NN - i);
		HR[i1] = XNI*PR[i1] / (double)N;
		HI[i1] = XNI*PI[i1] / (double)N;
	}

	for (jj = 0; jj < L1; ++jj) {
		if (CMOD(&HR[N-1], &HI[N-1]) > ETA*10.0*CMOD(&PR[N-1], &PI[N-1])) {
			prt = -1.0 * *PR;
			pit = -1.0 * *PI;
			CDIVID(&prt, &pit, &HR[N - 1], &HI[N - 1], INFIN, TR, TI);
			for (i = 1; i <= NM1; ++i) {
				j = NN - i - 1;
				j1 = j - 1;
				T1 = HR[j1];
				T2 = HI[j1];
				HR[j] = *TR*T1 - *TI*T2 + PR[j];
				HI[j] = *TR*T2 + *TI*T1 + PI[j];
			}
			HR[0] = PR[0];
			HI[0] = PI[0];
		}
		else {
			for (i = 1; i <= NM1; ++i) {
				j = NN - i - 1;
				HR[j] = HR[j - 1];
				HI[j] = HI[j - 1];
			}
			HR[0] = 0.0;
			HI[0] = 0.0;
		}
	}

}

int cpoly(double *OPR, double *OPI, int DEGREE, double *ZEROR, double *ZEROI) {
	/* FINDS THE ZEROS OF A COMPLEX POLYNOMIAL.
	C OPR, OPI  -  DOUBLE PRECISION VECTORS OF REAL AND
	C IMAGINARY PARTS OF THE COEFFICIENTS IN
	C ORDER OF DECREASING POWERS.
	C DEGREE    -  INTEGER DEGREE OF POLYNOMIAL.
	C ZEROR, ZEROI  -  OUTPUT DOUBLE PRECISION VECTORS OF
	C REAL AND IMAGINARY PARTS OF THE ZEROS.
	C FAIL      -  OUTPUT LOGICAL PARAMETER,  TRUE  ONLY IF
	C LEADING COEFFICIENT IS ZERO OR IF CPOLY
	C HAS FOUND FEWER THAN DEGREE ZEROS.
	C THE PROGRAM HAS BEEN WRITTEN TO REDUCE THE CHANCE OF OVERFLOW
	C OCCURRING. IF IT DOES OCCUR, THERE IS STILL A POSSIBILITY THAT
	C THE ZEROFINDER WILL WORK PROVIDED THE OVERFLOWED QUANTITY IS
	C REPLACED BY A LARGE NUMBER.*/
	int fail,CONV;
	double *PR, *PI, *HR, *HI, *QPR, *QPI, *QHR, *QHI, *SHR, *SHI;
	double SR, SI, TR, TI, PVR, PVI, ARE, MRE, ETA, INFIN;
	double XX, YY, COSR, SINR, SMALNO, BASE, XXX, ZR, ZI, BND;
	int NN,CNT1,CNT2,IDNN2,i;
	double prt, pit;

	PR = (double*)malloc(sizeof(double)* 100);
	PI = (double*)malloc(sizeof(double)* 100);
	HR = (double*)malloc(sizeof(double)* 100);
	HI = (double*)malloc(sizeof(double)* 100);
	QPR = (double*)malloc(sizeof(double)* 100);
	QPI = (double*)malloc(sizeof(double)* 100);
	QHR = (double*)malloc(sizeof(double)* 100);
	QHI = (double*)malloc(sizeof(double)* 100);
	SHR = (double*)malloc(sizeof(double)* 100);
	SHI = (double*)malloc(sizeof(double)* 100);

	mcon(&ETA, &INFIN, &SMALNO, &BASE);

	ARE = ETA;
	MRE = 2.0*sqrt(2.0)*ETA;
	XX = .70710678;
	YY = -XX;
	COSR = -.060756474;
	SINR = .99756405;
	fail = 0;
	NN = DEGREE + 1;

	if (!(OPR[0] != 0.0 || OPI[0] != 0.0)) {
		fail = 1;
		free(PR);
		free(PI);
		free(HR);
		free(HI);
		free(QPR);
		free(QPI);
		free(QHR);
		free(QHI);
		free(SHR);
		free(SHI);
		return fail;
	}

	while (!(OPR[NN - 1] != 0.0 || OPI[NN - 1] != 0.0)) {
		IDNN2 = DEGREE - NN + 2;
		ZEROR[IDNN2 - 1] = 0.0;
		ZEROI[IDNN2 - 1] = 0.0;
		NN = NN - 1;
	}

	for (i = 0; i < NN; ++i) {
		PR[i] = OPR[i];
		PI[i] = OPI[i];
		SHR[i] = CMOD(&PR[i], &PI[i]);
	}

	BND = RESCALE(NN, SHR, ETA, INFIN, SMALNO, BASE);

	if (BND != 1.0) {
		for (i = 0; i < NN; ++i) {
			PR[i] = BND*PR[i];
			PI[i] = BND*PI[i];
		}
	}

	NNITER: if (NN <= 2) {
				prt = -1.0 * PR[1];
				pit = -1.0 * PI[1];
				CDIVID(&prt, &pit, &PR[0], &PI[0], INFIN, &ZEROR[DEGREE - 1], &ZEROI[DEGREE - 1]);
				free(PR);
				free(PI);
				free(HR);
				free(HI);
				free(QPR);
				free(QPI);
				free(QHR);
				free(QHI);
				free(SHR);
				free(SHI);
				return fail;
			}

	for (i = 0; i < NN; ++i) {
		SHR[i] = CMOD(&PR[i], &PI[i]);
	}

	BND = CAUCHY(NN, SHR, SHI);

	/*
	C OUTER LOOP TO CONTROL 2 MAJOR PASSES WITH DIFFERENT SEQUENCES
	C OF SHIFTS.
	*/

	for (CNT1 = 1; CNT1 <= 2; ++CNT1) {
		//C FIRST STAGE CALCULATION, NO SHIFT.
		NOSHFT(NN, 5, &SR, &SI, PR, PI, QPR, QPI, &PVR, &PVI, ARE, MRE, ETA, INFIN, HR, HI, QHR, QHI, SHR, SHI, &TR, &TI);
		//C INNER LOOP TO SELECT A SHIFT.
		for (CNT2 = 1; CNT2 <= 9; ++CNT2) {
			//C SHIFT IS CHOSEN WITH MODULUS BND AND AMPLITUDE ROTATED BY
			//C 94 DEGREES FROM THE PREVIOUS SHIFT
			XXX = COSR*XX - SINR*YY;
			YY = SINR*XX + COSR*YY;
			XX = XXX;
			SR = BND*XX;
			SI = BND*YY;
			//C SECOND STAGE CALCULATION, FIXED SHIFT.
			FXSHFT(NN, 10*CNT2, &ZR, &ZI, &CONV, &SR, &SI, PR, PI, QPR, QPI, &PVR, &PVI, ARE, MRE, ETA, INFIN, HR, HI, QHR, QHI, SHR, SHI, &TR, &TI);
			if (CONV == 1) {
				//C THE SECOND STAGE JUMPS DIRECTLY TO THE THIRD STAGE ITERATION.
				//C IF SUCCESSFUL THE ZERO IS STORED AND THE POLYNOMIAL DEFLATED.
				IDNN2 = DEGREE - NN + 2;
				ZEROR[IDNN2 - 1] = ZR;
				ZEROI[IDNN2 - 1] = ZI;
				NN = NN - 1;
				for (i = 0; i < NN; ++i) {
					PR[i] = QPR[i];
					PI[i] = QPI[i];
				}
				goto NNITER;
			}//80
		}
	}

	fail = 1;

	free(PR);
	free(PI);
	free(HR);
	free(HI);
	free(QPR);
	free(QPI);
	free(QHR);
	free(QHI);
	free(SHR);
	free(SHI);
	return fail;
}

int polyroot(double *coeff, int DEGREE, double *ZEROR, double *ZEROI) {
	/*
	Input - coeff[0] + coeff[1] * x^1 + coeff[2] * x^2 + ---- + coeff[DEGREE+1] * x^DEGREE 

	The program fails if there are leading zeros. Check return value.

	Return Value = 0 -> Probable Success
	Return Value = 1 -> Likely Failure
	*/
	int fail,i,N;
	double *OPR, *OPI;

	N = DEGREE + 1;
	OPR = (double*)malloc(sizeof(double)* N);
	OPI = (double*)malloc(sizeof(double)* N);

	for (i = 0; i < N; ++i) {
		OPR[i] = coeff[N - 1 - i];
		OPI[i] = 0.0;
	}
	fail = cpoly(OPR, OPI, DEGREE, ZEROR, ZEROI);

	free(OPR);
	free(OPI);
	return fail;
}

int cpolyroot(double *rcoeff,double *icoeff, int DEGREE, double *ZEROR, double *ZEROI) {
	/*
	Input - rcoeff[0]+i*icoeff[0] + (rcoeff[1] + i*icoeff[1]) * x^1 + (rcoeff[2] + i*icoeff[2]) * x^2 + ---- + (rcoeff[DEGREE+1] + i*icoeff[DEGREE+1]) * x^DEGREE

	The program fails if there are leading zeros [At least one coefficient of the highest degree should be non-zero. Reduce degree by one if both are zero]. Check return value.

	Return Value = 0 -> Probable Success
	Return Value = 1 -> Likely Failure
	*/
	int fail, i, N;
	double *OPR, *OPI;

	N = DEGREE + 1;
	OPR = (double*)malloc(sizeof(double)* N);
	OPI = (double*)malloc(sizeof(double)* N);

	for (i = 0; i < N; ++i) {
		OPR[i] = rcoeff[N - 1 - i];
		OPI[i] = icoeff[N - 1 - i];
	}
	fail = cpoly(OPR, OPI, DEGREE, ZEROR, ZEROI);

	free(OPR);
	free(OPI);
	return fail;
}