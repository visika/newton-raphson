#include<math.h>
#include<stdio.h>

double f(double x);
double f1(double x);
int sign(double x);
double derivata_numerica(double x, double h);
double trova_radice(double x_0, double epsilon);

int main() {
	// trovare le radici dell'equazione f(x)=0; metodo di Newton-Raphson
	// punto di partenza x=0; ci sposteremo in avanti man mano; mi aspetto 4 radici reali
	double x_0 = 0; // punto di partenza della mia iterazione
	double epsilon = 0.00001;
	double radici[] = {0, 0, 0, 0}; // so che 0 non è radice, lo uso come segnaposto
	double step = 0.05; // avanzo di questo passo per trovare radici successive
	double radice;
	// radici successive; finché è la stessa di prima, ignora ed avanza il punto iniziale
	// non passo direttamente alla prima radice, poiché potrei averne totalmente saltata una nel mezzo
	for (int i=0; i<4; i++) {
		do {
			x_0 += step;
			printf("Newton parte da: %lf", x_0);
			radice = trova_radice(x_0, epsilon);
			printf("Radice trovata: %lf", radice);
		} while ((fabs(radice - radici[2]) < epsilon) || (fabs(radice - radici[1]) < epsilon) || (fabs(radice - radici[0]) < epsilon));
		radici[i] = radice;
	}
	for (int i=0; i<4; i++) {
		printf("Radice: x=%10.8lf\n", radici[i]);
	}
	
	
	
	// SECONDA PARTE
	
	
	
	// integrale definito tra xmin = 1 e xmax=4
	double xmin=1, a=1;
	double xmax=4, b=4;
	// metodo dei trapezi
	int n = 4;
	double h = (b-a)/n;
	double integrale = 0;
	double vecchio_integrale = 0;
	double somma = 0;
	somma += f(a);
	for (int i=1; i<=n-1; i++) {
		somma += 2*f(a+i*h);
	}
	somma += f(b);
	integrale = h/2*somma;
	while (fabs(integrale - vecchio_integrale) > epsilon) {
		vecchio_integrale = integrale;
		somma = 0;
		n *= 2;
		h = (b-a)/n;
		for (int i=1; i<=n-1; i+=2) {
			somma += f(a+i*h);
		}
		integrale = vecchio_integrale/2 + h*somma;
	}
	printf("L'integrale definito vale %lf\n", integrale);
	
	
	
	// TERZA PARTE
	// determinare la derivata numericamente
	
	int N = 100;
	double errore = 1E-3;
	double dn[N];
	double da[N];
	h = (b-a)/N;
	int coincidono;
	do {
		coincidono = 1;
		for (int i=0; i<N; i++) {
			dn[i] = derivata_numerica(a+i*h, errore);
			da[i] = f1(a+i*h);
			if (fabs(dn[i] - da[i]) >= errore) {
				coincidono = 0;
				break;
			}
		}
		h /= 2;
	} while (coincidono = 0);
	for (int i=0; i<N; i++) {
		printf("%3d) DN: %.15lf\t:::\tDA: %.15lf\n", i+1, dn[i], da[i]);
	}
	
	return 0;
}

double f(double x) {
	return x*x*x*x - 10*x*x*x + 35*x*x - 50*x + 24;
}

double f1(double x) {
	return 4*x*x*x - 30*x*x + 70*x - 50;
}

int sign(double x) {
	if (x > 0) return 1;
	if (x < 0) return -1;
	return 0;
}

double derivata_numerica(double x, double h) {
	return (f(x-2*h) - 8*f(x-h) + 8*f(x+h) - f(x+2*h))/(12*h);
}


double trova_radice(double x_0, double epsilon) {
	double x_i=x_0, x_i1=x_0;
	do {
		x_i = x_i1;
		x_i1 = x_i - f(x_i) / f1(x_i);
		 printf("x=%lf\n", x_i1); // per debugging
	} while (fabs(x_i - x_i1) > epsilon);
	double radice = x_i1;
	printf("Radice %lf\n", radice);
	return radice;
}


/* Tengo la vecchia versione di trova_radice, prima ripetuta ogni volta come qui sotto
do {
	x_i = x_i1;
	x_i1 = x_i - f(x_i) / f1(x_i);
	printf("x=%lf\n", x_i1);
} while (fabs(x_i - x_i1) > epsilon);
double radice = x_i1;
printf("La radice 1 è %lf\n", radice);
radici[0] = radice;
*/
/*
do {
	x_0 += step;
	radice = trova_radice(x_0, epsilon);
} while (fabs(radice - radici[0]) < 2*epsilon);
printf("La radice 2 è %lf\n", radice);
radici[1] = radice;
do {
	x_0+=step;
	radice = trova_radice(x_0, epsilon);
} while ((fabs(radice - radici[1]) < 2*epsilon) || (fabs(radice - radici[0]) < 2*epsilon));
printf("La radice 3 è %lf\n", radice);
radici[2] = radice;
do {
	x_0+=step;
	radice = trova_radice(x_0, epsilon);
} while ((fabs(radice - radici[2]) < 2*epsilon) || (fabs(radice - radici[1]) < 2*epsilon) || (fabs(radice - radici[0]) < 2*epsilon));
printf("La radice 4 è %lf\n", radice);
radici[3] = radice;
*/
