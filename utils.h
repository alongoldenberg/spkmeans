/**
* Created by Alon Goldenberg & Ofir Nissan.
*/

#ifndef SPKMEANS_UTILS_H
#define SPKMEANS_UTILS_H
#include <stdio.h>


double **allocate_data(int n, int d);
double distance(double *arr1, double *arr2, int d);
void print_error(void);
void print_invalid_input(void);
int check_is_num(char *str);
void print_matrix(double **matrix, int n, int d);
void print_arr(double *arr, int n);
int count_d(FILE *pFile);
int count_lines(FILE *file);


#endif
