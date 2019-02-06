// Bi-objective Next Release Problem (NRP) Version 3
//Date of last modification: may 2017 (december 2018-translation)
//Author: Miguel Angel Dominguez Rios
//Institution: University of Malaga


#include <stdio.h>
#include <stdlib.h>
#include <ilcplex/cplex.h>
#include <iostream>
#include <string>
#include <list>
#include <time.h>
#include <math.h>

using namespace std;


struct solution { //A "solution" contains the efficient solution, the nondominated point, their dimensions, the time when the solution is found and its accumulated hypervolume
	int dimension_x, dimension_y;
	double *vector_x, *vector_y;
	double time , acum_hyper;
};

struct BOILP {
	CPXENVptr *env;	//CPLEX environtment
	CPXLPptr *lp;	//CPLEX problem
	list<solution> *Efic; //List of solutions
	double **F;		//Coefficients of objective functions
	int n;			//Number of variables
	int m;			//Number of constraints
	int n_iterations;	//Number of calls to CPLEX
	double hypervolume;	//Hypervolume
	double max_time;	//Max time of execution (0 if no time limit)
	double total_time;	//Total time of execution
	char sense;			//Sense of the model: MAX or MIN
	bool Tune;			//Option of tune the CPLEX model

};

struct caja_bidimensional { //Bidimensional box
	double *z1, *z2;
	double value;
	int algorithm_to_use;
};


#include "algorithms2.h"
#include "biobjective_algorithms.h"


void mensaje_de_error() {
	printf("\n\nNRP (Input)(Max time execution)(Algorithm)(Option_1)[Option_2][Option_3]");
	printf("\n\n Input: Name of the archive to execute");
	printf("\n Max time execution: Measured in seconds. When no limit time, type '0'");
	printf("\n Algorithm: Choose between the following algorithms");
	printf("\n  'econst1': Epsilon constraint with 1 ILP per iteration");
	printf("\n  'econst2': Epsilon constraint with 2 ILPs per iteration");
	printf("\n  'augmecon': Improved Epsilon constraint, AUGMEGON");
	printf("\n  'spf': Supported Pareto Front");
	printf("\n  'hybrid': Hybrid (parametrization + epsilon constraint)");
	printf("\n  'tchebycheff': Augmented Thebycheff algorithm");
	printf("\n  'mixed': Mixed two types of algorithms in {hybrid, tchebycheff}");

	printf("\n\nALGORITHMS AND OPTIONS");
	printf("\n ALGORITHM  OPTION_1  OPTION_2  OPTION_3");
	printf("\n econst1 {f1,f2}");
	printf("\n econst2 {f1,f2}");
	printf("\n augmecon exact {f1,f2} {'fix','variable'}");
	printf("\n augmecon anytime {f1 , f2} {'fix','variable'}");
	printf("\n spf {exact,approximate}	");
	printf("\n hybrid exact");
	printf("\n hybrid anytime {normal,quadrants} {area,unexplored}");
	printf("\n tchebycheff {exact,anytime}	");
	printf("\n mixed {hybrid,tchebycheff} {tchebycheff,hybrid} [spf]");

	printf("\n\n  All options are strings except 'fix' (double) or 'variable' (0)\n");
}

int leer_problema_desde_archivo_NRP(CPXENVptr env, CPXLPptr lp, int *nn, int *mm, double ***F, char *nombre_fichero) {
	int status, i, j, n, m;
	double *b, *lb, *ub;
	char *tipo_rest, *tipo_var;
	FILE *fp;

	
	//Open file 
	fp = fopen(nombre_fichero, "r");
	if (fp == NULL) {
		printf("Impossible to open file %s \n", nombre_fichero);
		return 0;
	}

	//Reading file
	int n_levels, n_dep, n_cli, x;
	int req = 0; //Number of requirements
	std::list<int> listaF1, listaF2, listaRest; //Lists of F1 costs, F2 costs and the coefficient of constraints

	fscanf(fp, "%d", &n_levels); //Numer of requirements' levels

	for (i = 0; i < n_levels; i++) { //Number of requirements in every level
		fscanf(fp, "%d", &req);
		for (j = 0; j < req; j++) { //Costs of requirements in level i
			fscanf(fp, "%d", &x);
			listaF1.push_back(x);
		}
	}

	fscanf(fp, "%d", &n_dep);//Number of dependencies

	for (i = 0; i < n_dep; i++) { //Saving dependencies (constraint xi - xj >= 0)
		fscanf(fp, "%d", &x);
		listaRest.push_back(x);
		fscanf(fp, "%d", &x);
		listaRest.push_back(x);
	}

	fscanf(fp, "%d", &n_cli);//Number of stakeholders

	n = int(listaF1.size());	//Current number of variables

	for (i = 0; i < n_cli; i++) {
		n++;
		fscanf(fp, "%d", &x);
		listaF2.push_back(x); //Weight of the stakeholder
		fscanf(fp, "%d", &req); //Number of requirements for the stakeholder
		for (j = 0; j < req; j++) {
			fscanf(fp, "%d", &x);
			listaRest.push_back(x); //Add constraint xi - xk >= 0
			listaRest.push_back(n);
		}
	}

	//Close file
	fclose(fp);
	
	//Objective functions
	*F = (double **)malloc(2 * sizeof(double *));
	for (int i = 0; i < 2; i++) {
		(*F)[i] = (double *)malloc(n * sizeof(double));
	}

	i = 0;
	j = int(listaF1.size());
	for (i = 0; i < j; i++) {
		(*F)[0][i] = double(listaF1.front());		listaF1.pop_front();
		(*F)[1][i] = 0.0;
	}
	for (i = j; i < n; i++) {
		(*F)[0][i] = 0.0;
		(*F)[1][i] = double(-listaF2.front());		listaF2.pop_front();
	}

	if ((!listaF1.empty() || (!listaF2.empty()))) {
		printf("Reading data error. ");
		return 0;
	}

	m = int(listaRest.size()) / 2;

	int *indice_col = (int *)malloc(n * sizeof(int));
	lb = (double *)malloc(n * sizeof(double));
	ub = (double *)malloc(n * sizeof(double));
	b = (double *)malloc(m * sizeof(double));
	tipo_rest = (char *)malloc(m * sizeof(char));
	tipo_var = (char *)malloc(n * sizeof(char));
	int *matbeg;
	matbeg = (int *)malloc(m * sizeof(int));

	for (i = 0; i < m; i++) {
		tipo_rest[i] = 'G';
		b[i] = 0.0;
		matbeg[i] = 0;
	}


	for (j = 0; j < n; j++) {
		lb[j] = 0.0;
		ub[j] = 1.0; 
		tipo_var[j] = CPX_BINARY;
		indice_col[j] = j;
	}

	status = CPXaddrows(env, lp, n, m, 0, NULL, tipo_rest, matbeg, NULL, NULL, NULL, NULL);

	for (i = 0; i < m; i++) {
		x = int(listaRest.front());		listaRest.pop_front();
		status = CPXchgcoef(env, lp, i, x - 1, 1.0);
		x = int(listaRest.front());		listaRest.pop_front();
		status = CPXchgcoef(env, lp, i, x - 1, -1.0);
	}

	if (!listaRest.empty()) {
		printf("Error reading constraints. ");
		return 0;
	}

	char *tipo_lu = (char *)malloc(n * sizeof(char));
	for (i = 0; i < n; i++) {
		tipo_lu[i] = 'U';
	}

	status = CPXchgbds(env, lp, n, indice_col, tipo_lu, ub); //We change upper bounds for variables
	status = CPXchgprobtype(env, lp, CPXPROB_MILP); //Change CPLEX problem into mixed integer linear problem
	status = CPXchgctype(env, lp, n, indice_col, tipo_var); //Change variable types 
	status = CPXchgobj(env, lp, n, indice_col, (*F)[0]); //By default, the objective is F1.

	*nn = n;
	*mm = m;
	return (1);
}

int control_de_entrada(int argc, char **argv, double *tmax, std::string *opcion_elegida) {

	bool ejecutar = true;
	*opcion_elegida = "";

	printf("\nExecuting NRP");
	for (int g = 1; g < argc; g++) printf(" %s", argv[g]);
	printf("...");

	std::string a3, a4;
	if (argc >= 5) {
		a3 = argv[3];
		a4 = argv[4];
	}
	
	if (argc < 5) ejecutar = false;
	else if (argc == 5) {
		if (((a3 == "econst1") && (a4 == "f1")) ||			((a3 == "econst1") && (a4 == "f2")) ||
			((a3 == "econst2") && (a4 == "f1")) ||			((a3 == "econst2") && (a4 == "f2")) ||
			((a3 == "spf") && (a4 == "exact")) ||			((a3 == "spf") && (a4 == "approximate")) ||
			((a3 == "hybrid") && (a4 == "exact")) ||
			((a3 == "tchebycheff") && (a4 == "exact")) ||	((a3 == "tchebycheff") && (a4 == "anytime")) )	 {
			*opcion_elegida += a3; *opcion_elegida += "_";	  //algorithm
			*opcion_elegida += a4;   //option 1
		}
		else ejecutar = false;
	}
	else if (argc == 6) {
		std::string a5 = argv[5];
		
		if (((a3 == "hybrid") && (a4 == "anytime") && (a5 == "quadrants")) ||
			((a3 == "mixed") && (a4 == "hybrid") && (a5 == "tchebycheff")) ||
			((a3 == "mixed") && (a4 == "tchebycheff") && (a5 == "hybrid")) ){
			
			*opcion_elegida += a3; *opcion_elegida += "_";	  //algorithm
			*opcion_elegida += a4; *opcion_elegida += "_";	  //option 1
			*opcion_elegida += a5; *opcion_elegida ;	  //option 2
		}
		else ejecutar = false;
	}
	else if (argc == 7) {
		std::string a5 = argv[5];
		std::string a6 = argv[6];
		if (((a3 == "augmecon") && (a4 == "exact") && (a5 == "f1")) ||
			((a3 == "augmecon") && (a4 == "exact") && (a5 == "f2")) ||
			((a3 == "augmecon") && (a4 == "anytime") && (a5 == "f1")) ||
			((a3 == "augmecon") && (a4 == "anytime") && (a5 == "f2")) ||
			((a3 == "hybrid") && (a4 == "anytime") && (a5 == "normal") && (a6 == "area")) ||
			((a3 == "hybrid") && (a4 == "anytime") && (a5 == "normal") && (a6 == "unexplored")) ||
			((a3 == "mixed") && (a4 == "hybrid") && (a5 == "tchebycheff") && (a6 == "spf")) ||
			((a3 == "mixed") && (a4 == "tchebycheff") && (a5 == "hybrid") && (a6 == "spf")) ) {
			*opcion_elegida += a3; *opcion_elegida += "_";	  //algorithm
			*opcion_elegida += a4; *opcion_elegida += "_";	  //option 1
			*opcion_elegida += a5; *opcion_elegida += "_";	  //option 2
			*opcion_elegida += a6;	  //option 3
		}
		else ejecutar = false;
	}
	else ejecutar = false;

	if (!ejecutar) {
		mensaje_de_error();
		printf("\nUnable to complete the execution");
		return 0;
	}
	
	std::string a2 = argv[2];
	*tmax = atof(argv[2]); //Convert string to double (max time of execution)
	return 1;
}

void liberar_datos_main(BOILP *P1) {
	CPXfreeprob(*P1->env, P1->lp);
	CPXcloseCPLEX(P1->env);
	free(P1->F[0]);
	free(P1->F[1]);
	P1->Efic->clear();
	delete(P1->env);
	delete(P1->lp);
	delete(P1->Efic);
	delete(P1->F);
	delete(P1);
}

void Write_NRP_results_in_file(BOILP *P1, std::string *out, int argc, char **argv) {
	FILE *fp;
	std::string archive = argv[1];
	std::string algoritmo = argv[3];

	size_t eficsize = P1->Efic->size();
	double hipervolumen = P1->hypervolume;
	int numiteraciones = P1->n_iterations;
	double tmax = P1->max_time;
	double total_time = P1->total_time;
		
	const char *name = "NRP_results.txt";

	fp = fopen(name, "a+");
	char c;
	fread(&c, sizeof(c), 1, fp);	//We read one character to check if the file is empty

	if (feof(fp)) { //Si el fichero está vacio, creamos cabecera
		fprintf(fp, "NRP_Results\nDate      \tHour  \tInstance_name \tMaxt\tAlgorithm	 Opt1\tOpt2\tOpt3\t|N|\tIter\tHypervolume\tExecution Time");
		fprintf(fp, "\t\tPARALLELMODE\tTHREADS\tEPINT\tTune\tEPGAP\tEPAGAP\tEPRHS\tEPOPT\tADVIND\n");
	}
		
	rewind(fp); //Go to the beginning of file

	//Reading date and system's hour
	struct tm *newtime;
	time_t long_time;
	time(&long_time);
	newtime = localtime(&long_time);

	std::string fecha;
	if (newtime->tm_mday < 10) 		fecha += ("0");
	fecha += std::to_string(newtime->tm_mday);
	if (newtime->tm_mon < 9)		fecha += ("0");
	fecha += std::to_string(newtime->tm_mon + 1);		//Month
	fecha += std::to_string(newtime->tm_year + 1900);	//Year

	std::string hora;
	if (newtime->tm_hour < 10) hora += ("0");
	hora += std::to_string(newtime->tm_hour);			//Hours
	if (newtime->tm_min < 10) hora += ("0");
	hora += std::to_string(newtime->tm_min);			//Minutes
	if (newtime->tm_sec < 10) hora += ("0");
	hora += std::to_string(newtime->tm_sec);			//Seconds

	//Reading algorithm and options
	std::string Algoritmo;
	std::string opt1, opt2, opt3;

	opt1 = opt2 = opt3 = "----";

	if (algoritmo == "econst1") {
		Algoritmo = "E-const 1ILP";
		opt1 = "exacto";
		opt2 = argv[4];
		opt3 = "----";
	}
	else if (algoritmo == "econst2") {
		Algoritmo = "E-const 2ILP";
		opt1 = "exacto";
		opt2 = argv[4];
		opt3 = "----";
	}
	else if (algoritmo == "augmecon") {
		Algoritmo = "Augmecon";
		opt1 = argv[4];
		opt2 = argv[5];
		opt3 = argv[6];
		double lambda = atof(argv[6]);	//lambda value (0, if lambda is variable)
		if (lambda != 0) {
			opt3 = "l=" + opt3.substr(0, 9);	//Substring, we take the first 8 digits
		}
		else {
			opt3 = "----";
		}
	}
	else if (algoritmo == "spf") {
		Algoritmo = "F.Pareto Sop";
		opt1 = argv[4];
		opt2 = "----";
		opt3 = "----";
	}
	else if (algoritmo == "hybrid") {
		Algoritmo = "Hibrido";
		opt1 = argv[4];
		if (argc >= 6) opt2 = argv[5];
		if (argc == 7) opt3 = argv[6];
	}
	else if (algoritmo == "benson") {
		Algoritmo = "Benson";
		opt1 = argv[4];
		if (argc >= 6) opt2 = argv[5];
		if (argc == 7) opt3 = argv[6];
	}
	else if (algoritmo == "tchebycheff") {
		Algoritmo = "Tchebycheff";
		opt1 = argv[4];
	}
	else if (algoritmo == "mixed") {
		Algoritmo = "mixed";
		opt1 = argv[4];
		opt2 = argv[5];
		if (argc ==7) opt3 = argv[6];
	}

	//Write data in file
	int i_param[4];
	double d_param[5];
	CPXgetintparam(*P1->env, CPX_PARAM_PARALLELMODE, &i_param[0]);	
	CPXgetintparam(*P1->env, CPX_PARAM_THREADS, &i_param[1]);		
	CPXgetdblparam(*P1->env, CPX_PARAM_EPINT, &d_param[0]);
	i_param[2] = (P1->Tune) ? 1 : 0;
	CPXgetdblparam(*P1->env, CPX_PARAM_EPGAP, &d_param[1]);		
	CPXgetdblparam(*P1->env, CPX_PARAM_EPAGAP, &d_param[2]);		
	CPXgetdblparam(*P1->env, CPX_PARAM_EPRHS, &d_param[3]);
	CPXgetdblparam(*P1->env, CPX_PARAM_EPOPT, &d_param[4]);
	CPXgetintparam(*P1->env, CPX_PARAM_ADVIND, &i_param[3]);

	fprintf(fp, "%s\t%s\t", fecha.c_str(), hora.c_str());  //Date and hour
	fprintf(fp, "%s\t", archive.c_str());					//Input data name
	if (tmax == 0)
		fprintf(fp, "----\t");
	else	fprintf(fp, "%0.2f\t", tmax);

	fprintf(fp, "%s\t%s\t%s\t%s\t%d\t%d\t%0.1f\t%0.2f", Algoritmo.c_str(), opt1.c_str(), opt2.c_str(), opt3.c_str(), (int)eficsize, numiteraciones, hipervolumen, total_time); //Results
	fprintf(fp, "\t\t%d\t%d\t%2.9f\t%d\t%2.9f\t%2.9f\t%2.9f\t%2.9f\t%d", i_param[0], i_param[1], d_param[0], i_param[2], d_param[1], d_param[2], d_param[3], d_param[4], i_param[3]);
	fprintf(fp, "\n");

	fclose(fp);
	
	std::string aux = argv[1];
	aux += "_" + *out + "_" + fecha + hora;
	*out = aux;
}

void Print_solution_in_file(BOILP *P1, std::string nombre_fichero, std::list<solution> *L, double tiempo, double hipervolumen) {
	int contador = 0;
	std::list<solution>::iterator it;
	FILE *fp , *fp2;

#ifdef _WIN32
	std::string a = ".\\PARETO_FRONTS\\";
	_mkdir(a.c_str()); //New folder. Return -1 if it exists, 0 if don't
#elif __linux__
	std::string a = "./PARETO_FRONTS/";
	//const int dir_err = mkdir(a.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	mkdir(a.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif

	L->sort(ordenar_lex_2); //Arrange lexicographally the non-dominated points

	nombre_fichero = a + nombre_fichero;
	fp = fopen(nombre_fichero.c_str(), "w");
	nombre_fichero = nombre_fichero + "_time_hypervolume";
	fp2 = fopen(nombre_fichero.c_str(), "w");

	for (it = L->begin(); it != L->end(); it++) {
		fprintf(fp, "(");
		int i;
		for (i = 0; i < L->front().dimension_y - 1; i++)	fprintf(fp, "%1.0f, ", it->vector_y[i]);
		fprintf(fp, "%1.0f)\n", it->vector_y[i]);
		contador++;
		fprintf(fp2, "%0.5f\t%0.1f\n", it->time, it->acum_hyper);
	}
	int param_int;
	double param_double;
	CPXgetintparam(*P1->env, CPX_PARAM_PARALLELMODE, &param_int);	fprintf(fp, "\n\n\nCPX_PARAM_PARALLELMODE = %d        ", param_int);
	CPXgetintparam(*P1->env, CPX_PARAM_THREADS, &param_int);		fprintf(fp, "\nCPX_PARAM_THREADS = %d", param_int);
	CPXgetdblparam(*P1->env, CPX_PARAM_EPINT, &param_double);		fprintf(fp, "\nCPX_PARAM_EPINT = %2.9f", param_double);
	param_int = (P1->Tune) ? 1 : 0;									fprintf(fp, "\nTune = %d", param_int);
	CPXgetdblparam(*P1->env, CPX_PARAM_EPGAP, &param_double);		fprintf(fp, "\nCPX_PARAM_EPGAP = %2.9f;   ", param_double);
	CPXgetdblparam(*P1->env, CPX_PARAM_EPAGAP, &param_double);		fprintf(fp, "CPX_PARAM_EPAGAP = %2.9f;   ", param_double);
	CPXgetdblparam(*P1->env, CPX_PARAM_EPRHS, &param_double);		fprintf(fp, "CPX_PARAM_EPRHS = %2.9f;   ", param_double);
	CPXgetdblparam(*P1->env, CPX_PARAM_EPOPT, &param_double);		fprintf(fp, "CPX_PARAM_EPOPT = %2.9f;   ", param_double);
	CPXgetintparam(*P1->env, CPX_PARAM_ADVIND, &param_int);			fprintf(fp, "CPX_PARAM_ADVIND = %d", param_int);

	fprintf(fp, "\n\nTotal non-dominated points: %d ", contador);
	fprintf(fp, "\nTime : %2.2f seconds", tiempo);
	fprintf(fp, "\n\nHypervolume : %2.1f\n", hipervolumen);

	fclose(fp);
	fclose(fp2);
	
}

void Execute_NRP(BOILP *P1, int argc, char **argv) {

	double tmax = P1->max_time;
	std::string algoritmo = argv[3], a4 = argv[4], a5, a6;

	if (tmax == 0.0) P1->max_time = CPX_INFBOUND;
	if (argc > 5) a5 = argv[5];
	if (argc > 6) a6 = argv[6];

	if (algoritmo == "econst1") {
		int f_obj = (a4 == "f1") ? 0 : 1;
		Ejecutar_algoritmo_econstraint_1ILP(P1, f_obj);
	}
	else if (algoritmo == "econst2") {
		int f_obj = (a4 == "f1") ? 0 : 1;
		Ejecutar_algoritmo_econstraint_2ILP(P1, f_obj);
	}
	else if (algoritmo == "augmecon") {
		int f_obj = (a5 == "f1") ? 0 : 1;
		double lambda = atof(argv[6]);	//lambda value (0, if it is variable)

		if (a4 == "exact")			Ejecutar_algoritmo_augmecon_exact(P1 , f_obj, lambda);
		else if (a4 == "anytime")	Ejecutar_algoritmo_augmecon_anytime(P1, f_obj, lambda);
	}
	else if (algoritmo == "spf") {
		if (a4 == "exact") 				Ejecutar_algoritmo_SPF_exact(P1);
		else if (a4 == "approximate")	Ejecutar_algoritmo_SPF_approximate(P1);
	}
	else if (algoritmo == "hybrid") {
		enum { zero, box_area, unexplored_area };
		int tipo_valor_caja = (a4 == "exact") ? zero : ((a6 == "area") ? box_area : unexplored_area);
		
		if ((a4 == "anytime") && (a5 == "quadrants")) 	Ejecutar_algoritmo_hibrido_por_cuadrantes(P1);
		else											
			Ejecutar_algoritmo_hibrido(P1 , tipo_valor_caja);
	}
	else if (algoritmo == "tchebycheff") {		//Options EXACT O ANYTIME
		if (a4 == "exact") 			Ejecutar_algoritmo_tchebycheff(P1 , 0);
		else if (a4 == "anytime") 	Ejecutar_algoritmo_tchebycheff(P1 , 1);
	}
	else if (algoritmo == "mixed") {
		if (argc == 6) 			Ejecutar_algoritmo_mixed(P1, a4, a5);
		else if (argc == 7) 	Ejecutar_algoritmo_mixed_with_SPF(P1, a4, a5);
	}
	else {
		mensaje_de_error();
		exit(1);
	}
	if (tmax == 0.0) P1->max_time = 0.0;
}

int Cargar_datos_NRP(BOILP *P , double *tmax , int argc , char **argv) {

	CPXENVptr *env = new(CPXENVptr); //CPLEX environtment
	CPXLPptr *NRP = new(CPXLPptr);	//CPLEX problem
	std::list<solution> *Efic = new(std::list<solution>);
	int n = 0; 	//Number of variables (columns)
	int m = 0; //Number of constraints (rows)
	int *numiteraciones = new(int); //Number of iterations (calls to CPLEX)
	double **F = new(double *); //F[0] is the vector costs of f1;   F[1] is the vector costs of f2
	
	P->env = env;
	P->lp = NRP;
	P->Efic = Efic;
	P->F = F;
	P->n = P->m = 0;
	P->hypervolume = 0.0;
	P->max_time = *tmax;
	P->total_time = 0.0;
	P->n_iterations = 0;
	P->sense = 'L';
	P->Tune = false;
	
	const char *problem_name = "Next Release Problem";
	Abrir_problema_y_fijar_parametros_CPLEX(P ,problem_name);

	if (!leer_problema_desde_archivo_NRP(*P->env, *P->lp, &P->n, &P->m, &P->F, argv[1])) {
		printf("\n Impossible to read the file %s.\nEND\n", argv[1]);
		return 0;
	};
	
	return 1;
}

int main(int argc, char **argv) {
	
	BOILP *P1 = new(BOILP); //New bi-objective Integer Linear Program
	double tmax;			//Maximum execution time
	std::string out;		//Output file name (Pareto Front)
	TIEMPO t;				//Object time

	
	//Control of the input parameters
	if (!control_de_entrada(argc, argv , &tmax , &out)) return 0;
	if (!Cargar_datos_NRP(P1, &tmax, argc, argv)) 		return 0;
	

	t.init();
	Execute_NRP(P1 , argc , argv);
	t.interval(); //Total execution time	


	P1->total_time = t.value();

	Write_NRP_results_in_file(P1 ,&out , argc , argv);
	Print_solution_in_file(P1 , out, P1->Efic, P1->total_time, P1->hypervolume);
	printf("\n           Total time : %0.2f. |N| = %d. Hypervolume %0.1f\n", P1->total_time , int(P1->Efic->size()) , P1->hypervolume);
	

	//Free memory
	liberar_datos_main(P1);
	
	return 0;	
}


