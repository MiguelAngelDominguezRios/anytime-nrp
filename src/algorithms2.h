
using namespace std;

void Abrir_problema_y_fijar_parametros_CPLEX(BOILP *P1 ,const char *titulo) {
	int status;
	string a;
	
	CPXENVptr env = *P1->env;
	CPXLPptr lp = *P1->lp;

	env = CPXopenCPLEX(&status);
	lp = CPXcreateprob(env, &status, titulo);
	
	FILE *fp;
	char cad[80];
	string strcad;
	int cplexparameter;
	double value;
	
	fp = fopen("parameters_NRP", "r");
	if (fp == NULL) {
		printf("Unable to open parameters_NRP file");
		exit (1);
	}
	while (!feof(fp)) {
		fscanf(fp, "%s", cad);	
		strcad = cad;
		fscanf(fp, "%lf", &value);

		if (strcad == "CPX_PARAM_PARALLELMODE") {	cplexparameter = CPX_PARAM_PARALLELMODE;		status = CPXsetintparam(env, cplexparameter, int(value));	}
		else if (strcad == "CPX_PARAM_THREADS") {	cplexparameter = CPX_PARAM_THREADS;				status = CPXsetintparam(env, cplexparameter, int(value));	}
		else if (strcad == "CPX_PARAM_EPINT") {		cplexparameter = CPX_PARAM_EPINT;				status = CPXsetdblparam(env, cplexparameter, value);		}
		else if (strcad == "CPX_PARAM_EPGAP") {		cplexparameter = CPX_PARAM_EPGAP;				status = CPXsetdblparam(env, cplexparameter, value);		}
		else if (strcad == "CPX_PARAM_EPAGAP") {	cplexparameter = CPX_PARAM_EPAGAP;				status = CPXsetdblparam(env, cplexparameter, value);		}
		else if (strcad == "CPX_PARAM_EPRHS") {		cplexparameter = CPX_PARAM_EPRHS;				status = CPXsetdblparam(env, cplexparameter, value);		}
		else if (strcad == "CPX_PARAM_EPOPT") {		cplexparameter = CPX_PARAM_EPOPT;				status = CPXsetdblparam(env, cplexparameter, value);		}
		else if (strcad == "CPX_PARAM_ADVIND") {	cplexparameter = CPX_PARAM_ADVIND;				status = CPXsetintparam(env, cplexparameter, int(value));	}
		else if (strcad == "Tune") {
			int x = int(value);
			P1->Tune = (x == 0) ? false : true;
		}
	}
	fclose(fp);

	*P1->env = env;
	*P1->lp = lp;
}

int comparar_vectores(double *a, double *b, int tamano) {
	//Compare if two vectors have the same components
	int i;
	for (i = 0; i < tamano; i++) {
		if (a[i] != b[i]) return 0;
	}
	return (1);
}

double mult(double *x, double *y, int n) {
	//Given the vectors   x = (x1 , x2 , ... , xn)     y = (y1 , y2 , ... , yn), return the value :   sum(xi*yi)
	double valor = 0.0;

	for (int i = 0; i < n; i++) {
		valor += x[i] * y[i];
	}
	return valor;
}

void Solve_Econstraint_1(int row_to_modify, double new_eps, CPXENVptr *env, CPXLPptr *lp, int *stat, double *objval, int *numiteraciones) {
	int m = row_to_modify;
	CPXchgrhs(*env, *lp, 1, &m, &new_eps);	//Update constraint   f1 <= eps  (f2 <= eps)
	CPXmipopt(*env, *lp);
	CPXgetobjval(*env, *lp, objval);
	*stat = CPXgetstat(*env, *lp);
	(*numiteraciones)++;
}

int Solve_Econstraint_2(int row_to_modify, double new_eps, CPXENVptr env, CPXLPptr lp, CPXLPptr lp2 ,  int *indice_col, int n, int f_obj, double **F, double *objval,
	int *numiteraciones, double **nuevo_y) {

	int m = row_to_modify;
	int stat, f_rest = 1 - f_obj;

	//First subproblem
	CPXchgrhs(env, lp, 1, &m, &new_eps);	//s.a. f2 <= eps
	CPXmipopt(env, lp);
	CPXgetobjval(env, lp, objval);	//The objective will be the new eps
	stat = CPXgetstat(env, lp);
	(*numiteraciones)++;
	if (stat == CPXMIP_INFEASIBLE) return 0;
	(*nuevo_y)[f_obj] = round(*objval); 
	
	//Second subproblem
	new_eps = round(*objval);
	CPXchgrhs(env, lp2, 1, &m, &new_eps);
	
	CPXmipopt(env, lp2);
	CPXgetobjval(env, lp2, objval);
	stat = CPXgetstat(env, lp2);
	(*numiteraciones)++;
	(*nuevo_y)[f_rest] = round(*objval); 
	return 1;
}


int Solve_Econstraint_22(int row_to_modify, double new_eps, CPXENVptr env, CPXLPptr lp, int *indice_col, int n, int f_obj, double **F, double *objval,
	int *numiteraciones, double **nuevo_y) {

	int m = row_to_modify, j;
	int stat, f_rest = 1 - f_obj;

	
	//First subproblem
	CPXchgobj(env, lp, n, indice_col, F[f_obj]); //min f1
	for (j = 0; j < n; j++) CPXchgcoef(env, lp, m, j, F[f_rest][j]);  //s.a f2 <= eps
	CPXchgrhs(env, lp, 1, &m, &new_eps);	//s.a. f2 <= eps
	CPXmipopt(env, lp);
	CPXgetobjval(env, lp, objval);	//The objective will be the new eps
	stat = CPXgetstat(env, lp);
	(*numiteraciones)++;
	if (stat == CPXMIP_INFEASIBLE) return 0;
	
	(*nuevo_y)[f_obj] = round(*objval);

	//Second subproblem
	CPXchgobj(env, lp, n, indice_col, F[f_rest]); //min f2
	for (j = 0; j < n; j++)	CPXchgcoef(env, lp, m, j, F[f_obj][j]);  //We modify the mth-constraint in order to create the constraint f1 <= f2(x*)
	new_eps = round(*objval);
	CPXchgrhs(env, lp, 1, &m, &new_eps);
	CPXmipopt(env, lp);
	CPXgetobjval(env, lp, objval);
	stat = CPXgetstat(env, lp);
	(*numiteraciones)++;

	
	(*nuevo_y)[f_rest] = round(*objval);
	return 1;
}


int Solve_Augmecon(double *F_obj, double *F_rest, double lambda, int row_to_modify, double new_eps, CPXENVptr env, CPXLPptr lp, int *indice_col, int n, 
	double *objval,	int *numiteraciones) {

	//Execute biobjective AUGMEGON algorithm using CPLEX
	//Input: f1-costs, f2-costs, lambda, index of constraint f2+s2 <= eps2, value eps2
	//Output: objval , numiteraciones + 1
	
	//AUGMECON algorithm
	//In every iteration we do one CPLEX call, and obtain one efficient solution
	//
	// min (f1(x) - lambda * s2)   
	// s.t.  f2(x) + s2 <= eps2
	//			s2 >= 0
	//			x en X
	// (it can be f1(x) <--> f2(x))

	int stat;
	double oo = -lambda;
	int m = row_to_modify;
	CPXchgobj(env, lp, 1, &n, &oo); //min f2 - lambda*s
	CPXchgrhs(env, lp, 1, &m, &new_eps);	// s.a f1 + s <= new_eps
	CPXmipopt(env, lp);
	CPXgetobjval(env, lp, objval);
	stat = CPXgetstat(env, lp);
	(*numiteraciones)++;
	
	if ( (stat == CPXMIP_OPTIMAL) || (stat == 102) ) return 1;
	return 0;
}

bool Solve_SPF_approximate(double *z1, double *z2, int n, double **F, CPXENVptr env, CPXLPptr lp, int *indice_col,
	int *numiteraciones, double **nuevo_x, double **z3) {
	//Search a solution which is paremeterized by z1 and z2, and must be z1, z2 or a new z3
	//The output is never infeasible 
	//Return true if a new solution is found
	//Return false if the solution is previously found (z1 or z2)

	int i;
	double *nuevo_c = (double *)malloc(n * sizeof(double));
	double objval;

	//Normal vector of the line determined by z1 and z2
	double lambda1 = z1[1] - z2[1];
	double lambda2 = z2[0] - z1[0];

	for (i = 0; i < n; i++) nuevo_c[i] = lambda1 * F[0][i] + lambda2 * F[1][i];

	CPXchgobj(env, lp, n, indice_col, nuevo_c); //min k1*f1 + k2*f2
	CPXmipopt(env, lp);
	CPXgetobjval(env, lp, &objval);
	(*numiteraciones)++;

	free(nuevo_c);

	CPXgetx(env, lp, *nuevo_x, 0, n - 1);
	(*z3)[0] = round(mult(*nuevo_x, F[0], n));
	(*z3)[1] = round(mult(*nuevo_x, F[1], n));

	if ((comparar_vectores(z1, *z3, 2) || (comparar_vectores(z2, *z3, 2)))) {
		return false;
	}
	else {
		return true;
	}
}

bool Solve_SPF(double *z1, double *z2, int n, int row_reference, double **F, CPXENVptr env, CPXLPptr lp, int *indice_col,
	int *numiteraciones, double **nuevo_x, double **z3) {
	// Look for a new solution parameterized by z1 and z2, which might be z1, z2, or a new z3.
	
	int i;
	double *nuevo_c = (double *)malloc(n * sizeof(double));
	double objval;

	//Normal vector of the line determined by z1 and z2
	double lambda1 = z1[1] - z2[1];
	double lambda2 = z2[0] - z1[0];

	double objcaja = lambda1 * z1[0] + lambda2 * z1[1]; //Objective value in the extreme of the box

	for (i = 0; i < n; i++)	nuevo_c[i] = lambda1 * F[0][i] + lambda2 * F[1][i];

	int indices[2] = { row_reference , row_reference + 1 };
	double values[2] = { z2[0] - 1 , z1[1] - 1 };

	CPXchgobj(env, lp, n, indice_col, nuevo_c); //min k1*f1 + k2*f2
	CPXchgrhs(env, lp, 2, indices, values); //s.a.  f1 <= k1,   f2 <= k2

	CPXmipopt(env, lp);
	CPXgetobjval(env, lp, &objval);
	int stat = CPXgetstat(env, lp);
	(*numiteraciones)++;

	free(nuevo_c);
	
	if ( ((stat == CPXMIP_OPTIMAL) || (stat == 102)) && (objval <= objcaja)) {
		CPXgetx(env, lp, *nuevo_x, 0, n - 1);
		(*z3)[0] = round(mult(*nuevo_x, F[0], n));
		(*z3)[1] = round(mult(*nuevo_x, F[1], n));
		return true;
	}
	return false;
}

int calcula_parametros_augmented_Tchebycheff(double z11, double z12, double z21, double z22, double *w1, double *w2, double *ro) {
	double x = z21 - z11;
	double y = z12 - z22;

	double u = 0.1;
	
	if ( (x > y) && (y >= 2) ){
		*w1 = (x * y - x - y + u *(2 - u)) / (x * y - y - 3 * x + x * x + 2 * u * (2 - u));
		*w2 = 1 - *w1;
		*ro = (x - u)*(1 - u) / (x * y - y - 3 * x + x * x + 2 * u * (2 - u));
	}
	else if ((x == y) && (y >= 2)) {
		*w1 = *w2 = 0.5;
		*ro = (1 - u) / (2 * (x + u - 2));
	}
	else if ((x < y) && (x >= 2)) {
		*w1 = ((y - u) * (y + u - 2)) / (x * y - x - 3 * y + y * y + 2 * u * (2 - u));
		*w2 = 1 - *w1;
		*ro = (y - u)*(1 - u) / (x * y - x - 3 * y + y * y + 2 * u * (2 - u));
	}
	else {
		return 0;
	}
	return 1;
}

bool Solve_Tchebycheff(double *z1, double *z2, double *ideal ,int n, int row_reference, double **F , CPXENVptr env, CPXLPptr lp, int *indice_col,
	int *numiteraciones, double **nuevo_x, double **z3) {
	//Execute biobjective AUGMENTED TCHEBYCHEFF algorithm usin CPLEX
	//Input: extreme values of the box (z1 and z2)
	//Output: an efficient solution, the non-dominated point, numiterations + 1
	
	int i;
	double w1, w2, ro;
	calcula_parametros_augmented_Tchebycheff(z1[0], z1[1], z2[0], z2[1], &w1, &w2, &ro);

	double objval;

	ideal[0] = z1[0];	ideal[1] = z2[1];
	
	//min (lambda + ro * (z1-s1 + z2-s2))
	double *values = (double *)malloc(n * sizeof(double));
	for (i = 0; i < n; i++) values[i] = ro * (F[0][i] + F[1][i]);
	CPXchgobj(env, lp, n, indice_col, values);
	CPXchgobjoffset(env, lp, ro * (-ideal[0] - ideal[1]));
	
	//s.t. wi* zi -lambda <= wi * si    i = 1,2
	double rhs;	int indice;

	indice = row_reference; rhs = w1 * ideal[0];
	for (i = 0; i < n; i++) CPXchgcoef(env, lp, row_reference, i, w1 * F[0][i]);
	CPXchgrhs(env, lp, 1, &indice, &rhs); 
	
	indice = row_reference + 1; rhs = w2 * ideal[1];
	for (i = 0; i < n; i++) CPXchgcoef(env, lp, row_reference + 1, i, w2 * F[1][i]);
	CPXchgrhs(env, lp, 1, &indice, &rhs);
	
	indice = row_reference + 2; rhs = z2[0];
	CPXchgrhs(env, lp, 1, &indice, &rhs);

	indice = row_reference + 3; rhs = z1[1];
	CPXchgrhs(env, lp, 1, &indice, &rhs);
	
	CPXmipopt(env, lp);
	CPXgetobjval(env, lp, &objval);
	int stat = CPXgetstat(env, lp);
	(*numiteraciones)++;

	if ((stat == CPXMIP_OPTIMAL) || (stat == 102) ){
		CPXgetx(env, lp, *nuevo_x, 0, n);

		(*z3)[0] = round(mult(*nuevo_x, F[0], n));
		(*z3)[1] = round(mult(*nuevo_x, F[1], n));
		if ((comparar_vectores(z1, *z3, 2)) || (comparar_vectores(z2, *z3 , 2))) {
			return false;
		}
		else {
			return true;
		}
	}
	else { 
		printf("ERROR FOUND");
		return false;
	}
}

bool Solve_Hibrido_cuadrantes(double zz0 , double zz1, int n, int row_reference, double **F, CPXENVptr env, CPXLPptr lp, 
	int *numiteraciones, double **nuevo_x, double **z3) {
	//Look for a new solution parameterized by z1 and z2, which might be z1, z2 or a new z3
	
	int indices[2] = { row_reference , row_reference + 1 };
	double objval, values[2] = { zz0 - 1 , zz1 - 1 };
	
	CPXchgrhs(env, lp, 2, indices, values); //s.a.  f1 <= k1,   f2 <= k2
	CPXmipopt(env, lp);
	CPXgetobjval(env, lp, &objval);
	int stat = CPXgetstat(env, lp);
	(*numiteraciones)++;

	if ((stat == CPXMIP_OPTIMAL) || (stat == 102)) {
		CPXgetx(env, lp, *nuevo_x, 0, n - 1);
		(*z3)[0] = mult(*nuevo_x, F[0], n);
		(*z3)[1] = mult(*nuevo_x, F[1], n);
		return true;
	}
	else {
		return false;
	}
}

bool Solve_Hibrido(double *z1, double *z2, int n, int row_reference, double **F, CPXENVptr env, CPXLPptr lp, int *indice_col,
	int *numiteraciones, double **nuevo_x, double **z3) {
	//Look for a new solution parameterized by z1 and z2, which might be z1, z2 or a new z3
	
	int i;
	double *nuevo_c = (double *)malloc(n * sizeof(double));
	double objval;

	//Normal vector of the line determined by z1 and z3
	double lambda1 = z1[1] - z2[1];
	double lambda2 = z2[0] - z1[0];

	for (i = 0; i < n; i++) nuevo_c[i] = lambda1 * F[0][i] + lambda2 * F[1][i];

	int indices[2] = { row_reference , row_reference + 1 };
	double values[2] = { z2[0] - 1 , z1[1] - 1 };

	CPXchgobj(env, lp, n, indice_col, nuevo_c); //min k1*f1 + k2*f2
	CPXchgrhs(env, lp, 2, indices, values); //s.a.  f1 <= k1,   f2 <= k2
	CPXmipopt(env, lp);
	CPXgetobjval(env, lp, &objval);
	int stat = CPXgetstat(env, lp);
	(*numiteraciones)++;

	free(nuevo_c);

	if ((stat == CPXMIP_OPTIMAL) || (stat == 102)) {
		CPXgetx(env, lp, *nuevo_x, 0, n - 1);

		(*z3)[0] = round(mult(*nuevo_x, F[0], n));
		(*z3)[1] = round(mult(*nuevo_x, F[1], n));
		return true;
	}
	else {
		return false;
	}
}

