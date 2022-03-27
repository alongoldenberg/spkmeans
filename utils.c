//
// Created by Alon Goldenberg on 07/03/2022.
//

#include "utils.h"


double distance(double *arr1, double *arr2, int d){
    double res = 0; int i; double delta; double sqr_delta;
    for (i=0; i<d; i++){
        delta = arr1[i] - arr2[i];
        sqr_delta = delta*delta;
        res += sqr_delta;
    }
    return sqrt(res);
}


double **allocate_data(int n, int d) {
    double *p; double **a; int i;
    p = (double*)calloc(n*d, sizeof(double));
    a = (double**) calloc(n,sizeof(double *));
    if(a == NULL || p == NULL){
        print_error();
    }
    for(i=0; i<n; i++) {
        a[i] = p + i*d;
    }
    return a;
}

static void print_error(){
    printf("An Error Has Occurred!\n");
    exit(1);
}

void print_invalid_input(){
    printf("Invalid Input!\n");
    exit(1);
}

int check_is_num(char *str) {
    int i;
    for(i=0;str[i]!=0;i++){
        if(!isdigit(str[i])){
            print_invalid_input();
        }
    }
    return 1;
}

void print_matrix(double **matrix, int n, int d){
    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < d; j++)
        {
            printf("%d\t", matrix[i][j]);
        }
        printf("\n");
    }
}