#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>
#include<time.h>
#include<stdlib.h>
#include "winbgi2.h"
#include "board.h"
#include<stdbool.h>

#define WIN_WIDTH	800
#define WIN_HEIGHT	600
#define G			50

void landing_page();

void graphics_init() {
	graphics(WIN_WIDTH, WIN_HEIGHT);
	landing_page();
}


void landing_page() {

	int mouse_x = 0;
	int mouse_y = 0;

	int loop = true;
	while (loop) {
		rectangle(10, 10, 780, 580);

		// Title
		settextjustify(1, 1);
		settextstyle(9, HORIZ_DIR, 5);
		outtextxy(WIN_WIDTH/2, WIN_HEIGHT/3, "Totolotek");

		// Start
		rectangle(280, 250, 380, 300);
		settextjustify(1, 1);
		settextstyle(0, HORIZ_DIR, 1);
		outtextxy(330, 280, "Start");

		// QUIT
		rectangle(420, 250, 520, 300);
		settextjustify(1, 1);
		settextstyle(0, HORIZ_DIR, 1);
		outtextxy(470, 280, "Exit");

		int loop2 = true;
		while (loop2) {
			mouse_x = mouseclickx();
			mouse_y = mouseclicky();

			if ((mouse_x > 280) && (mouse_x < 380) && (mouse_y > 250) && (mouse_y < 300)) {
				printf("START\n");
				mouse_x = 0;
				mouse_y = 0;
				loop2 = false;
				clear();
				sim_init();
				
			}
			//QUIT
			if ((mouse_x > 420) && (mouse_x < 520) && (mouse_y > 250) && (mouse_y < 300)) {
				printf("QUIT\n");
				mouse_x = 0;
				mouse_y = 0;
				loop2 = false;
				loop = false;
			}

			delay(100);
		}

	}
}