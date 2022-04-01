//
// Created by Alon Goldenberg on 07/03/2022.
//

#include "utils.h"

int n ,d;
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
            printf("%fl\t", matrix[i][j]);
        }
        printf("\n");
    }
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
    int c=0;
    rewind(pFile);
    d = 1;
    while ((c = fgetc(pFile) ) != '\n' ){
        if(c==',')
            d++;
    }
    rewind(pFile);
    return d;
}
double **parse_file(char *filename){
    /* Allocate memory and put data from file to a 2d array */
    FILE *data = NULL; double **arr; double num; int i,j;
    
    data = fopen(filename, "r");
    if(data == NULL){
        print_invalid_input();
    }
    n = count_lines(data);
    d = count_d(data);
    arr = allocate_data(n, d);
    rewind(data);
    for(i = 0; i<n; i++){
        for(j=0; j<d; j++) {
            fscanf(data, "%lf,", &num);
            arr[i][j] = num;
        }
    }
    fclose(data);
    return arr;
}