#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>
#include<time.h>
#include<stdlib.h>
#include "winbgi2.h"
#include "interface.h"
#include<stdbool.h>

#define WIN_WIDTH	800
#define WIN_HEIGHT	600
#define G			50


void init(int* x, int* y, double* vx, double* vy, double* m, int N, int* cl);
void init2(int* x, int* y, double* vx, double* vy, double* m, int N, int* cl);

double get_random(double min, double max);
void status(int* x, int* y, double* vx, double* vy, double* m, int N);
void displacement(int* x, int* y, double* vx, double* vy, int N, double t);
void detect_wall_collision(int* x, int* y, double* vx, double* vy, double* m, int N, int* cl);
void find_max(double* x, double* y, int N);
void kinetic_en(double* vx, double* vy, double* m, double* kin, int N);
void find_closest(int* x, int* y, int N);
void detect_ball_collision(int* x, int* y, double* vx, double* vy, double* m, int N);
void detect_sphere_collision(int* x, int* y, double* vx, double* vy, double* m, int N);
void init3(int* x, int* y, double* vx, double* vy, double* m, int N);
void dmuchawa(int* x, int* y, double* vx, double* vy, double* m, int N);



void sim_init()
{
	//graphics(WIN_WIDTH, WIN_HEIGHT);

	int N = 10;

	int* x, * y;						// wsp pilek
	double* vx, * vy;				// skladowe predkosci pilek
	double* m;						// masy pilek
	int* cl;						// collision logic
	double* kin;					// energia kinetyczna 

	x = (int*)malloc(N * sizeof(int));
	y = (int*)malloc(N * sizeof(int));
	vx = (double*)malloc(N * sizeof(double));
	vy = (double*)malloc(N * sizeof(double));
	m = (double*)malloc(N * sizeof(double));
	cl = (int*)malloc(N * sizeof(int));
	kin = (double*)malloc(N * sizeof(double));



	// nadawanie pi�kom randomowych parametrów
	//init(x,	y,vx,vy,m,N, cl);	//cl
	//init2(x, y, vx, vy, m, N, cl);	//cl

	init3(x, y, vx, vy, m, N);	//cl
	int mouse_x = 0;
	int mouse_y = 0;

	int blowing = 0;
	int start = 0;

	double t = 0.0167;				// czas jaki zabiera skok
	for (;;) {
		clear();
		mouse_x = mouseclickx();
		mouse_y = mouseclicky();



		// start 
		rectangle(40, 40, 240, 100);
		settextjustify(1, 1);
		settextstyle(0, HORIZ_DIR, 3);
		outtextxy(140, 80, "START");

		// PRZYCISKi
		// menu
		rectangle(660, 40, 760, 90);
		settextjustify(1, 1);
		settextstyle(0, HORIZ_DIR, 1);
		outtextxy(710, 65, "menu");

		if ((mouse_x > 660) && (mouse_x < 760) && (mouse_y > 40) && (mouse_y < 90)) {
			mouse_x = 0;
			mouse_y = 0;

			clear();
			free(x);	free(y);	free(vx);	free(vy);	free(m);	free(cl);	 free(kin);

			landing_page();
			break;
		}

		// dmuchawa on and off 
		settextjustify(1, 1);
		settextstyle(0, HORIZ_DIR, 1);
		outtextxy(710, 150, "Dmuchawa");

		rectangle(660, 175, 705, 225);
		settextjustify(1, 1);
		settextstyle(0, HORIZ_DIR, 1);
		outtextxy(680, 200, "on");


		rectangle(710, 175, 760, 225);
		settextjustify(1, 1);
		settextstyle(0, HORIZ_DIR, 1);
		outtextxy(735, 200, "off");

		// on
		if ((mouse_x > 660) && (mouse_x < 705) && (mouse_y > 175) && (mouse_y < 225)) {
			
			//printf("dmucahwa\n");
			mouse_x = 0;
			mouse_y = 0;
			blowing = 1;
		}
		if ((mouse_x > 710) && (mouse_x < 760) && (mouse_y > 175) && (mouse_y < 225)) {

			
			mouse_x = 0;
			mouse_y = 0;
			blowing = 0;
		}

		// załączanie dmuchawy
		if (blowing == 1) {
			circle(657, 147, 5);
			dmuchawa(x, y, vx, vy, m, N);
		}

		detect_ball_collision(x, y, vx, vy, m, N);
		displacement(x, y, vx, vy, N, t);
		detect_sphere_collision(x, y, vx, vy, m, N);
		status(x, y, vx, vy, m, N);


		delay(20);
	}
	free(x);	free(y);	free(vx);	free(vy);	free(m);	free(cl);	 free(kin);

	wait();
}


// randomowe rozmieszczenie z randomowymi prędkościami
void init(int* x, int* y, double* vx, double* vy, double* m, int N, int* cl) {
	srand(time(NULL));
	for (int i = 0; i < N; i++) {

		int middle_x = WIN_WIDTH / 2;
		int middle_y = WIN_HEIGHT / 2;
		int radius = 250;



		x[i] = (int)get_random(250, 250 + 300);
		y[i] = (int)get_random(150, 150 + 300);
		vx[i] = get_random(-10, 10);
		vy[i] = get_random(-10, 10);
		m[i] = 30.0;
		//m[i] = get_random(15, 25);

		cl[i] = 0;
	}

	double overlap;
	double centers_unit_vector[2];

	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {

			double distance_between_centers = sqrt(pow(x[i] - x[j], 2) + pow(y[i] - y[j], 2));

			if (distance_between_centers < (int)(m[i] + m[j])) {

				// przeniesienie kulek tak aby sie tylko stykaly
				// uzyty do tego zostanie wektor jednostkowy, 
				centers_unit_vector[0] = ((double)(x[i] - x[j])) / distance_between_centers;
				centers_unit_vector[1] = ((double)(y[i] - y[j])) / distance_between_centers;

				// przemieszczenie piłek tak aby tylko stykały się ze sobą
				overlap = 0.5 * (distance_between_centers - m[j] - m[i]);

				x[i] -= overlap * centers_unit_vector[0];
				y[i] -= overlap * centers_unit_vector[1];

				x[j] += overlap * centers_unit_vector[0];
				y[j] += overlap * centers_unit_vector[1];
			}
		}
	}
}

// do testowania zderzeń idealnie sprężystych 
void init2(int* x, int* y, double* vx, double* vy, double* m, int N, int* cl) {
	srand(time(NULL));
	for (int i = 0; i < N; i++) {

		// 1.
		//x[0] = 400;
		//y[0] = 300;
		x[0] = (int)get_random(250, 250 + 300);
		y[0] = (int)get_random(150, 150 + 300);
		//vx[0] = -10;
		//vy[0] = 10;
		m[0] = 30.0;
		vx[i] = 20;
		vy[i] = 0;

		// 2.
		/*
		x[1] = 500;
		y[1] = 300;
		vx[1] = -10;
		vy[1] = 0;
		m[1] = 30.0;
		*/
		//m[i] = get_random(15, 25);

		/*
		x[2] = 250;
		y[2] = 100;
		vx[2] = 0;
		vy[2] = 0;
		m[2] = 30.0;
		*/
		cl[i] = 0;
	}
}

// rozmieszcza piłeczki w centrum okręgu, bez prędkości początkowych
void init3(int* x, int* y, double* vx, double* vy, double* m, int N) {
	// rozmieszczenie piłek tak jak w instrukcji
	for (int i = 0; i < N; i++) {

		int middle_x = WIN_WIDTH / 2;
		int middle_y = WIN_HEIGHT / 2;
		int radius = 250;


		x[i] = (int)get_random(250, 250 + 300);
		y[i] = (int)get_random(150, 150 + 300);
		vx[i] = 0.0;
		vy[i] = 0.0;
		m[i] = 30.0;
		//m[i] = get_random(15, 25);

	}
}

// zwraca randomowego doubla z zakresu [min, max] (gownie wykorzystuje ja init())
double get_random(double min, double max) {
	double random = ((double)rand()) / RAND_MAX;
	double range = max - min;
	return random * range + min;
}

// zwraca wartosci obecny stan wszystkich kulek + wypisuje numery kulek
void status(int* x, int* y, double* vx, double* vy, double* m, int N) {
	for (int i = 0; i < N; i++) {

		// printf("%lf %lf %lf %lf %lf\n", x[i], y[i], vx[i], vy[i], m[i]);	// parametry kazdego z 
		circle(x[i], y[i], (int)m[i]);										// rysowanie kolek

		char score_s[4] = "0";
		sprintf(score_s, "%d", i + 1);
		settextjustify(1, 1);
		settextstyle(0, HORIZ_DIR, 2);
		outtextxy(x[i], y[i] + 5, score_s);

		//line(x[i], y[i], x[i] + vx[i] * 10, y[i] + vy[i] * 10);						//rysowanie wektora predkosci
	}
	int middle_x = WIN_WIDTH / 2;
	int middle_y = WIN_HEIGHT / 2;
	int radius = 250;

	circle(middle_x, middle_y, radius);
	circle(middle_x, middle_y, 5);
}

// przemieszczenie 
void displacement(int* x, int* y, double* vx, double* vy, int N, double t) {

	for (int i = 0; i < N; i++) {

		// zmiana pr�dko�ci sk�adowych
		vx[i] = vx[i];
		vy[i] = vy[i] + (G * t);
		// zmiana przemieszczenia na podstawie pr�dko�ci sk�adowych
		x[i] = x[i] + (int)vx[i];
		y[i] = y[i] + (int)vy[i];

		//printf("vx[%d] = %lf \t ", i, vx[i]);
	}
	//printf("\n");



}

// wykrywanie kolizji z podloga i scianami 
// NOTE - nie potrzebne przy finalnej wersji
void detect_wall_collision(int* x, int* y, double* vx, double* vy, double* m, int N, int* cl) {
	for (int i = 0; i < N; i++) {
		//podloga
		if (((y[i] + (int)m[i]) > WIN_HEIGHT)) {
			vy[i] = 0.8 * (-1) * vy[i];
			y[i] = WIN_HEIGHT - (int)m[i];
		}
		//sufit
		if (((y[i]) - (int)m[i] < 0) && cl[i] == 0) {
			vy[i] = 0.8 * (-1) * vy[i];
			y[i] = y[i] + 1;
		}
		// prawa sciana
		if (((x[i]) + (int)m[i] > WIN_WIDTH)) {
			vx[i] = 0.8 * (-1) * vx[i];
			x[i] = x[i] - 1;
		}
		// lewa	sciana
		if (((x[i] - (int)m[i]) < 0)) {
			vx[i] = 0.8 * (-1) * vx[i];
			x[i] = 1 + (int)m[i];

		}
	}
}

// TODO - naprawic zderzenia z okręgiem
void detect_sphere_collision(int* x, int* y, double* vx, double* vy, double* m, int N) {

	int middle_x = WIN_WIDTH / 2;
	int middle_y = WIN_HEIGHT / 2;
	int radius = 250;

	double overlap;
	double centers_unit_vector[2] = { 0.0, 0.0 };
	double vel_unit_vector[2];

	for (int i = 0; i < N; i++) {

		// odleglosc od srodka sfery
		double center_dist = sqrtf(pow(x[i] - middle_x, 2) + pow(y[i] - middle_y, 2));

		centers_unit_vector[0] = ((double)(x[i] - middle_x)) / center_dist;
		centers_unit_vector[1] = ((double)(y[i] - middle_y)) / center_dist;

		if ((int)(center_dist + m[i]) > radius) {

			double temp_x = (double)x[i];
			double temp_y = (double)y[i];

			double offset = radius - center_dist - m[i];
			
			x[i] += (int)(offset * centers_unit_vector[0]);
			y[i] += (int)(offset * centers_unit_vector[1]);


			centers_unit_vector[0] = ((double)(x[i] - middle_x)) / center_dist;
			centers_unit_vector[1] = ((double)(y[i] - middle_y)) / center_dist;

			double vel_value = sqrtf(vx[i] * vx[i] + vy[i] * vy[i]);
			vel_unit_vector[0] = vx[i] / vel_value;
			vel_unit_vector[1] = vy[i] / vel_value;

			// narzedzia do rysowania
			// wektor promienia
			//line(middle_x, middle_y, middle_x + centers_unit_vector[0] * 800, middle_y + centers_unit_vector[1] * 800);
			//line(x[i], y[i], x[i] + vel_unit_vector[0] * 500, y[i] + vel_unit_vector[1] * 500);


			// zmiania wektora predkosci
			//vx[i] = 0.8*(-1) * vx[i];
			//vy[i] = 0.8*(-1) * vy[i];

			double cos_alpha = (vel_unit_vector[0] * centers_unit_vector[0] + vel_unit_vector[1] * centers_unit_vector[1]) / (1 + 1);
			double cos_beta = (centers_unit_vector[0] * (-1) + centers_unit_vector[1] * 0) / (2);

			//printf("arccos alpha: %lf\n", acos(cos_alpha));
			//printf("arccos beta: %lf\n", acos(cos_beta));

			/* zmiana wektora predkosci - interpretacja
				1. traktujemy zderzenie piłki o ściane tak jakby to były dwie piłki
			*/

			// n w tutorialu

			// tangent
			double tx = centers_unit_vector[0];
			double ty = centers_unit_vector[1];

			// Dot Porduct Tangent
			double dpTan1 = vx[i] * tx + vy[i] * ty;
			//double dpTan2 = vx[j] * tx + vy[j] * ty;

			line(x[i], y[i], x[i] + vx[i] * 10, y[i] + vy[i] * 10);

			vx[i] = -0.8 * tx * dpTan1;
			vy[i] = -0.8 * ty * dpTan1;

			//printf("vx %lf vy %lf \n", vx[i], vy[i]);


			//printf("out of sphere\n");
			//status(x, y, vx, vy, m, N);
			//delay(1000);
		}
	}
}

// Działa dobrze
void detect_ball_collision(int* x, int* y, double* vx, double* vy, double* m, int N) {

	double relative_vel[2];
	double distance_between_centers;
	double overlap;
	double centers_unit_vector[2];
	double vn[2];
	double v2[2];
	double v1[2];
	double check_sum1;
	double check_sum2;


	relative_vel[0] = 0.0;
	relative_vel[1] = 0.0;

	centers_unit_vector[0] = 0.0;
	centers_unit_vector[1] = 0.0;

	vn[0] = 0.0; vn[1] = 0.0;
	v2[0] = 0.0; v2[1] = 0.0;
	v1[0] = 0.0; v1[1] = 0.0;

	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {

			distance_between_centers = sqrtf(pow(x[i] - x[j], 2) + pow(y[i] - y[j], 2));

			if (distance_between_centers < (int)(m[i] + m[j])) {

				printf("yeet\n");
				// przeniesienie kulek tak aby sie tylko stykaly
				// uzyty do tego zostanie wektor jednostkowy, 
				centers_unit_vector[0] = ((double)(x[i] - x[j])) / distance_between_centers;
				centers_unit_vector[1] = ((double)(y[i] - y[j])) / distance_between_centers;

				// przemieszczenie piłek tak aby tylko stykały się ze sobą
				overlap = 0.5 * (distance_between_centers - m[j] - m[i]);

				x[i] -= overlap * centers_unit_vector[0];
				y[i] -= overlap * centers_unit_vector[1];

				x[j] += overlap * centers_unit_vector[0];
				y[j] += overlap * centers_unit_vector[1];

				// v1 w tutorialu
				relative_vel[0] = vx[i] - vx[j];
				relative_vel[1] = vy[i] - vy[j];

				// L w tutorialu
				distance_between_centers = sqrtf(pow(x[i] - x[j], 2) + pow(y[i] - y[j], 2));

				// n w tutorialu
				centers_unit_vector[0] = ((double)(x[i] - x[j])) / distance_between_centers;
				centers_unit_vector[1] = ((double)(y[i] - y[j])) / distance_between_centers;

				// tangent
				double tx = -centers_unit_vector[1];
				double ty = centers_unit_vector[0];

				// Dot Porduct Tangent
				double dpTan1 = vx[i] * tx + vy[i] * ty;
				double dpTan2 = vx[j] * tx + vy[j] * ty;

				// Dot Product Normal
				double dpNorm1 = vx[i] * centers_unit_vector[0] + vy[i] * centers_unit_vector[1];
				double dpNorm2 = vx[j] * centers_unit_vector[0] + vy[j] * centers_unit_vector[1];

				double m1 = 0.8 * (dpNorm1 * (m[i] - m[j]) + 2.0 * m[j] * dpNorm2) / (m[i] + m[j]);
				double m2 = 0.8 * (dpNorm2 * (m[j] - m[i]) + 2.0 * m[i] * dpNorm2) / (m[i] + m[j]);

				vx[i] = tx * dpTan1 + centers_unit_vector[0] * m1;
				vy[i] = ty * dpTan1 + centers_unit_vector[1] * m1;

				vx[j] = -(tx * dpTan2 + centers_unit_vector[0] * m2);
				vy[j] = ty * dpTan2 + centers_unit_vector[1] * m2;

				/*
				printf("ball 2 status:\n tx %lf \t dpTan2: %lf \t centers_unit_vector: %lf \t m1: %lf \n",
					tx,
					dpTan2,
					centers_unit_vector[0],
					m1
				);
				*/

			}
		}
	}
}

void dmuchawa(int* x, int* y, double* vx, double* vy, double* m, int N) {
	rectangle(10, 400, 790, 600);
	for (int i = 0; i < N; i++) {
		// definiowanie obaszaru gdzie wieje dmucha wiatr a gdzie nie
		if ((x[i] > 10 && x[i] < 790) && (y[i] > 400 && y[i] < 600)) {
			vy[i] += -5.0;
			vx[i] += -2.0;
		}
	}
}


// (cw. z lab: Wykyrwanie minimalnej i maksymalnej wsp x oraz y wraz indeksami tych pilek)
void find_max(double* x, double* y, int N) {
	double indeks_x = WIN_WIDTH;
	double indeks_y = WIN_HEIGHT;

	for (int i = 0; i < N; i++) {
		if (x[i] < indeks_x) {
			indeks_x = x[i];
		}
		if (y[i] < indeks_y) {
			indeks_y = y[i];
		}
	}
	printf("min x: %lf, y: %lf\n ", indeks_x, indeks_y);
}

// (cw. z lab: Obliczanie minimalnej, sredniej i maksymalnej energii kinetycznej pilek, wypisac te energie na ekranie
void kinetic_en(double* vx, double* vy, double* m, double* kin, int N) {

	long double kin_sum = 0.0;

	// obliczanie energii kinetycznej dla ka�dej pi�eczki
	for (int i = 0; i < N; i++) {
		double kinetic = (0.5) * m[i] * (pow(vx[i], 2) + pow(vy[i], 2));			// energia kinetyczna pojedynczej pi�eczki
		kin[i] = kinetic;
		kin_sum = kin_sum + kinetic;
	}

	// znajdowanie najmniejszej 
	double smallest_kin = 10000;
	double highest_kin = 0;


	for (int i = 0; i < N; i++) {
		if (kin[i] > highest_kin) {
			highest_kin = kin[i];
		}
		if (kin[i] < smallest_kin) {
			smallest_kin = kin[i];
		}
	}
	printf("highest kin: %lf\t lowest kin: %lf\t avg: %lf\n", highest_kin, smallest_kin, kin_sum / (double)N);
	// znajdowanie najwi�kszej

	// wyznaczanie �redniej

}

// (cw. z lab: Znajdz najmniejsza i najwieksza old miedzy pilkami. Wypisz te odl oraz rysuje linie pomiedzy srodkami tych pilek
void find_closest(int* x, int* y, int N) {

	double smallest_dist = 1000.0;
	double longest_dist = 0.0;

	double s_x1, s_x2, s_y1, s_y2;
	double l_x1, l_x2, l_y1, l_y2;


	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {

			long double distance = sqrt(pow(x[i] - x[j], 2) + pow(y[i] - y[j], 2));

			if (distance > longest_dist) {
				longest_dist = distance;
				l_x1 = x[i];
				l_y1 = y[i];
				l_x2 = x[j];
				l_y2 = y[j];
			}
			if (distance < smallest_dist) {
				smallest_dist = distance;
				s_x1 = x[i];
				s_y1 = y[i];
				s_x2 = x[j];
				s_y2 = y[j];
			}
		}
	}

	// rysowanie lini pomi�dzy �rodkami tych k�ek
	line(s_x1, s_y1, s_x2, s_y2);
	line(l_x1, l_y1, l_x2, l_y2);

}

