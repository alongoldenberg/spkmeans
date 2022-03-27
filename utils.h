//
// Created by Alon Goldenberg on 07/03/2022.
//

#ifndef SPKMEANS_UTILS_H
#define SPKMEANS_UTILS_H

double **allocate_data(int n, int d);
double distance(double *arr1, double *arr2, int d);
static void print_error();
void print_invalid_input();
int check_is_num(char *str);
void print_matrix(double **matrix, int n, int d);

#endif //SPKMEANS_UTILS_H
