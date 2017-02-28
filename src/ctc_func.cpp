/*************************************************************************
	> File Name: ctc_func.cpp
	> Author: Chenggaofeng
	> Mail: Chenggaofeng@hccl.ioa.ac.c 
	> Created Time: Sat 05 Sep 2015 12:37:26 PM CST
 ************************************************************************/

#include<cmath>
#include<iostream>
#include<vector>
using namespace std;
void sub_max(float **(&params),int array_n,int array_m)
{
    //search the largest value in the column of the array/matrix
	fprintf(stdout, "rows(array_m) : %d\ncols(array_n) : %d\n", array_m, array_n);
    for(int i=0;i<array_m;i++)
    {
		fprintf(stdout, "Check the %dth frame to find the best prob of current phoneme state\n", i);
        float colum_max=params[0][i];
        for(int j=0;j<array_n;j++)
        {
			fprintf(stdout, "Check the %dth state to find the best prob of phoneme state at this frame\n", j);
            if(params[j][i]>colum_max)
                colum_max=params[j][i];
        }
        for(int j=0;j<array_n;j++)
        {
            params[j][i]=params[j][i]-colum_max;
        }
    }
}
void array_exp(float **(&params),int array_n,int array_m)
{
    for(int i=0;i<array_n;i++)
        for(int j=0;j<array_m;j++)
            params[i][j]=exp(params[i][j]);
}
void array_sum_column(float (**params),float *(&sum_of_column),int array_n,int array_m)
{
    float tempt=0;
    for(int j=0;j<array_m;j++)
    {
        tempt=0;
        for(int i=0;i<array_n;i++)
        {
        tempt+=params[i][j];    
        }
        sum_of_column[j]=tempt;
    }
}
void array_norm(float **(&params),int array_n,int array_m)
{
	fprintf(stdout, "Do softmax computation here\n");
    for(int j=0;j<array_m;j++)
    {
        float col_sum=0;
        for(int i=0;i<array_n;i++)
            col_sum+=params[i][j];
        for(int i=0;i<array_n;i++)
            params[i][j]=params[i][j]/col_sum;
    }
}

//*******************there may exit algorithm problem*******************************
float array_sum_sincolumn(float **array,int num_of_columns,int t=0,int start=0,int end=0)
{
	//here
	fprintf(stdout, "Here, we add all the probs(softmax outputs) of u at time(frame) t=%d\n",t);
	fprintf(stdout, "u from %d to %d(do not include %d) on u axis at time t : %d\n", start, end, end, t);
    float c=0;
    for(int i=start;i<end;i++)
    {
        c+=array[i][t];
    }
    return c;
}
//*********************************************************
void delete_bidimen_array(float **array,int n,int m)
{
    for(int i=0;i<n;i++)
        delete []array[i];
    delete []array;
}
//array_n stands for numPhones,array_m stands for T,2*seqLen+1 for u
float ctc_loss(float **(&grad),float **params,int array_n,int array_m,
	int *seq,int seqLen,int blank=0,bool is_prob=true)
{
    //1.initializing the value 
	fprintf(stdout, "Phone number(Softmax output number of NN, array_n) : %d\nFrame number(the length of training observations, array_m) : %d\nSequence length(seqLen, purely symbols comes from symbol set A) : %d\n",
		array_n, array_m, seqLen);
	fprintf(stdout, "The sequence is the current example for training\n\n");
	
    int L=seqLen*2+1;
    int T=array_m;

	fprintf(stdout, "CTC trellis : frame length : %d, labellings length (symbols from A'): %d\n",
		T, L);
    
    float **alphas;
    float **beta;
    //2.new the space for alphas and beta
    alphas=new float *[L];
    for(int i=0;i<L;i++)
    {
        alphas[i]=new float[array_m];
        for(int j=0;j<array_m;j++)
        {
            alphas[i][j]=0;
        }
    }

    beta=new float *[L];
    for(int i=0;i<L;i++)
    {
        beta[i]=new float[array_m];
        for(int j=0;j<array_m;j++)
        {
            beta[i][j]=0;   
        }
    }

#ifdef _DEBUG
	std::cout << "CTC loss trellis have been allocated and initialized here\n";
#endif
    //3.transfer the params into prob form if necessary 
	std::cout << "如果必要，将参数转换成概率形式" << std::endl;
    if (is_prob==false)
    {  
        //substracting the max
        sub_max(params,array_n,array_m);    
        //exp the vals of the array
        array_exp(params,array_n,array_m);
        //normalise the array in the column direction
        array_norm(params,array_n,array_m);
    }
    //4.initialize alphas and the forward pass
    alphas[0][0]=params[blank][0];//init blank path start
    alphas[1][0]=params[seq[0]][0];//init normal label path start
    float c=0;
    //communicate the column num of alphas to sum function
    c=array_sum_sincolumn(alphas,L,0,0,L-1);

    alphas[0][0]=params[blank][0]/c;
    alphas[1][0]=params[seq[0]][0]/c;

    static float llforward=log(c);
    //used for test
    int count=0;
    int testnum=0;
    //used for test
    //5.calculating the val of alphas
	fprintf(stdout, "\n\n\nHere, we compute all values of alpha and store them to params \n");
    for(int t=1;t<T;t++)
    {
		fprintf(stdout, "computing alpha at time (frame index) : %d from %dth frame\n", t, t-1);
        int start=(0>(L-2*(T-t))) ? 0:(L-2*(T-t));
        int end=((2*t+2)>L) ? L:(2*t+2);
#ifdef _DEBUG
		fprintf(stdout, "current u from %d to %d(do not include %d) \n", start, end, end);
#endif
        for(int s=start;s<L;s++)
        {
            int l=(s-1)/2;
            if(s%2 == 0)
            {
                if(s==0)
                    alphas[s][t]=alphas[s][t-1]*params[blank][t];
                else
                    alphas[s][t]=(alphas[s][t-1]+alphas[s-1][t-1])*params[blank][t];
            }
            else if(s==1 || seq[l]==seq[l-1] )
            {
                alphas[s][t]=(alphas[s][t-1]+alphas[s-1][t-1])*params[seq[l]][t];
            }
            else
            {
                alphas[s][t]=(alphas[s][t-1]+alphas[s-1][t-1]+alphas[s-2][t-1])*params[seq[l]][t];
            }
        }
        //used for debug>
        cout<<"Col index : "<<  count+1  <<"   of CTC computation matrix : \n";
        count++;
        for(int i=0;i<L;i++)
			fprintf(stdout, "u = %d, t = %d, alpha[%d][%d] = %f\n", i, t, i, t, alphas[i][t]);
		std::cout << std::endl;
        //used for debug<
        c=array_sum_sincolumn(alphas,L,t,start,end);
        //used for debug>
        cout<<endl;
		fprintf(stdout, "The c VALUE : %f at time (frame index) : %d \n", c, t);
        
        //used for debug<
        for(int i=start;i<end;i++)
        {
            alphas[i][t]=alphas[i][t]/c;
        }
        cout<<endl;
        
        llforward+=log(c);
        //used for debug
        std::cout<<"llforward : "<<  llforward  <<std::endl;
        //used for debug
    }
        //used for debug>
        for(int t=0;t<T;t++)
        {
            for(int i=0;i<L;i++)
				fprintf(stdout, "u = %d, t = %d, alpha[%d][%d] = %f\n", i, t, i, t, alphas[i][t]);
            cout<<endl;
        }
        //used for debug<
    
    //6.initialize the betas and the backward pass computation
    beta[L-1][T-1]=params[blank][T-1];
    beta[L-2][T-1]=params[seq[seqLen-1]][T-1];
    c=beta[L-1][T-1]+beta[L-2][T-1];
    for(int i=0;i<L;i++)
        beta[i][T-1]=beta[i][T-1]/c;
    static float llbackward=log(c);
    //debug module
    count=T-2;
    //debug module
    //7.calculating the val of the beta
#ifdef _DEBUG
	fprintf(stdout, "\n\n\nHere, we compute all values of beta and store them to params \n");
#endif
    for(int t=T-2;t>=0;t--)
    {
        int start=(0>=(L-(T-t)*2)) ? 0:(L-(T-t)*2);
        int end=(L<=(2*t+2)) ? L:(2*t+2);
		fprintf(stdout, "current u from %d to %d \n", start, end);
        for(int s=(end-1);s>=0;s--)
        {
            int l=(s-1)/2;
            if(s%2==0)
            {
                if (s==(L-1))
                    beta[s][t]=beta[s][t+1]*params[blank][t];
                else
                    beta[s][t]=(beta[s][t+1]+beta[s+1][t+1])*params[blank][t];
            }
            else if((s==(L-2))|| (seq[l]==seq[l+1]))
            beta[s][t]=(beta[s][t+1]+beta[s+1][t+1])*params[seq[l]][t];
            else
            beta[s][t]=(beta[s][t+1]+beta[s+1][t+1]+beta[s+2][t+1])*params[seq[l]][t];
        }
        //*used for debug>
        cout<<"Col index :"<<  count  <<  "   *** : \n";
        count--;
        for(int i=0;i<L;i++)
			fprintf(stdout, "u = %d, t = %d, beta[%d][%d] = %f\n", i, t, i, t, beta[i][t]);
        cout<<endl;
        ///used for debug<
        c=array_sum_sincolumn(beta,L,t,start,end);
        for(int i=start;i<end;i++)
        {
            beta[i][t]=beta[i][t]/c;
        }
        llbackward+=log(c);
		std::cout<<"llbackward : "  <<  llbackward  << std::endl;
    }
//###############################################################################################
//###calculate the grad#########################################################################

    grad=new float *[array_n];
    for(int i=0;i<array_n;i++)
    {
        grad[i]=new float [T];
        for(int j=0;j<T;j++)
            grad[i][j]=0;
    }
	fprintf(stdout, "Grad computation, grad matrix : rows(Softmax output) : %d, cols(frame length) : %d\n", array_n, T);
//claculating the alphas*beta>
	fprintf(stdout, "Computating alphas * beta here\n");
    float **ab;
    ab=new float *[L];
    for(int i=0;i<L;i++)
    {
        ab[i]=new float [T];
            for(int t=0;t<T;t++)
			{
				ab[i][t]=alphas[i][t]*beta[i][t];
				fprintf(stdout, "u = %d, t = %d, alpha*beta[%d][%d] = %f\n", i, t, i, t, ab[i][t]);
			}
    }
//claculating the alphas*beta<
    for(int s=0;s<L;s++)
    {
        if(s%2==0)
        {
        for(int t=0;t<T;t++)
        {
            grad[blank][t]+=ab[s][t];
            ab[s][t]=ab[s][t]/params[blank][t];
        }
        }
        else
        {
        for(int t=0;t<T;t++)
        {
            grad[seq[(s-1)/2]][t]+=ab[s][t];
            ab[s][t]=ab[s][t]/(params[seq[(s-1)/2]][t]);
        }
        }
    }
	fprintf(stdout, "Done alpha*beta normalization\n");
    float *absum;
    absum=new float[T];
    array_sum_column(ab,absum,L,T);
//module used to judge the output>
    float llDiff=abs(llforward-llbackward);
    bool absum_0=false;
    for(int t=0;t<T;t++)
        if(absum[t]==0)
        {
            absum_0=true;
        }
//>
//cout<<"the raw grad##############:"<<endl;
//for(int i=0;i<array_n;i++)
//{
//for(int t=0;t<T;t++)
//  cout<<grad[i][t]<<" ";
//cout<<endl;
//}
//<
    if(llDiff>1e-5  || absum_0)
        {
            cout<<"Difference in forward/backward LogLikelyhood:"<<llDiff<<endl;
            cout<<"zero found in absum"<<endl;
            delete_bidimen_array(alphas,L,T);
            delete_bidimen_array(beta,L,T);

            return -llforward;
        }
    else
    {
//module used to judge the output<
        for(int t=0;t<T;t++)
            for(int i=0;i<array_n;i++)
            {
                grad[i][t]=params[i][t]-grad[i][t]/(params[i][t]*absum[i]);
            }
        delete_bidimen_array(alphas,L,T);
        delete_bidimen_array(beta,L,T);
        return -llforward;
    }
}
//here the input of the func is params
void best_path_decode(float **(&params),int array_n,int array_m,vector<int> &seq,int blank=0)
{
    int *seq_decode;
    seq_decode=new int[array_m];
    for(int t=0;t<array_m;t++)
    {
        int max_phone_of_time_t=0;
        float max_val_of_time_t=params[0][t];
        for(int u=0;u<array_n;u++)
        {
        if(params[u][t]>max_val_of_time_t)
        {
            max_val_of_time_t=params[u][t];
            max_phone_of_time_t=u;
        }
        }
        seq_decode[t]=max_phone_of_time_t;
    }
    seq.push_back(0);//push a blank symbols and delete it later
    for(int t=0;t<array_m;t++)
    {
        int &tempt=seq.back();
        if((t==0)&&(seq_decode[t]!=0))
            seq.push_back(seq_decode[t]);
        else if(((t>=1)&&(seq_decode[t]!=0))&&(seq_decode[t]!=tempt))
            seq.push_back(seq_decode[t]);
    }
}
