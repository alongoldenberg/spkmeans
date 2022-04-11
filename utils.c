/**
* Created by Alon Goldenberg & Ofir Nissan.
*/
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

void print_matrix(double **matrix, int n, int d) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < d; j++) {
            if(j==(d-1)){
                printf("%.4f", matrix[i][j]);
            }
            else{
                printf("%.4f,", matrix[i][j]);
            }
        }
        printf("\n");
    }
}

void print_arr(double *arr, int n) {
    int i;
    for (i = 0; i < n; i++) {
            if(i==(n-1)){
                printf("%.4f", arr[i]);
            }
            else{
                printf("%.4f,", arr[i]);
            }
        }
    printf("\n");
}



int count_lines(FILE *file) {
    char c = 0;
    int n = 0;
    while ((c = fgetc(file) ) != EOF){
        if(c=='\n')
            n++;
    }
    return n;
}

int count_d(FILE *pFile) {
    int d, c=0;
    rewind(pFile);
    d = 1;
    while ((c = fgetc(pFile) ) != '\n' ){
        if(c==',')
            d++;
    }
    rewind(pFile);
    return d;
}

