/*************************************************************************
	> File Name: ctc_func.h
	> Author: Chenggaofeng
	> Mail: Chenggaofeng@hccl.ioa.ac.c 
	> Created Time: Sat 05 Sep 2015 12:37:39 PM CST
 ************************************************************************/
#ifndef CTC_FUNC_H
#define CTC_FUNC_H
#include<iostream>
#include<vector>
using namespace std;
//sub_max() is func used to find the max value of every column,
//and substract the max value from the Matrix
void sub_max(float **(&params),int array_n,int array_m);

//array_exp() is func used to get the exp value of every element of the Matrix
void array_exp(float **(&params),int array_n,int array_m);

//array_sum_column() is func used to sum every column of the Matrix and put the summation of the column to
//the parameter sum_of_column.
void array_sum_column(float (**params),float *(&sum_of_column),int array_n,int array_m);

//array_norm() is func used to normalise the Matrix in the column direction
void array_norm(float **(&params),int array_n,int array_m);

//array_sum_sincolumn() is func used to sum the Matrix in the column direction,
//but you have to declare the num of the column and the beginning and ending of 
//the summation
float array_sum_sincolumn(float **array,int num_of_columns,int t=0,int start=0,int end=0);

//delete_bidimen_array() is the func used to delete the two dimension dynamic ptr_
void delete_bidimen_array(float **array,int n,int m);

//ctc_loss() is func used to calculate the ctcloss and the grad of the ctc NN
float ctc_loss(float **(&grad),float **params,int array_n,int array_m,int *seq,int seqLen,int blank=0,bool is_prob=true);

//best_path_decode() is func used to decode the NN using the best-path method
void best_path_decode(float **(&params),int array_n,int array_m,vector<int> &seq,int blank=0);

#endif

