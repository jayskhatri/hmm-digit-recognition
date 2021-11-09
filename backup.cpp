#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <vector>

#include <algorithm>
#include <sstream>

#define CONVERGE_ITERATIONS 200
#define M 32 //Number of obsevation symbols per state
#define N 5 //Number of states
#define P 12
#define LIMIT 5000
#define CB_SIZE 32
#define PI 3.142857142857
#define FRAME_SIZE 320
#define FRAMES 60

const int WIN_SIZE  = (FRAME_SIZE * FRAMES);

const int T = 60; //Time sequence length
using namespace std;

//Global variables
int count_samples = 0;
int index_max;
int digit_count = 0;			//variable to store digits number, which is currently being tested.
int Num_Recordings = 20;	
long double dcShift, nFactor, mx;
long double const threshold = 1e-30;   //Min threshold to be assigned to zero values in matrix B.
long int sSize = 0;
long double max_pobs_model = 0;
int test_ans = 0;

//Globally defined arrays
int O[T+1];	//Observation sequence
int Q[T+1];	//state sequence.
long double pstar = 0, prev_p_star = -1;
long double Alpha[T+1][N+1];
long double Beta[T+1][N+1];
long double Gamma[T+1][N+1];
long double Delta[T+1][N+1];
int Psi[T+1][N+1]; 
long double Xi[T+1][N+1][N+1];

long double codeBook[CB_SIZE][P];

long int sample[1000000];
long double energy[100000];
long double Ai[P+1], Ri[P+1], Ci[P+1];
//tokhura weights
double tokhuraWeights[]={1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};

//Model parameters A, B and Pi
long double A[N+1][N+1] = {0};
long double B[N+1][M+1] = {0};
long double Pi[N+1] = {0};

long double Abar[N+1][N+1] = {0};
long double Bbar[N+1][M+1] = {0};
long double Pibar[N+1] = {0};

long double a_average[N+1][N+1] = {0};
long double b_average[N+1][M+1] = {0};
long double pi_average[N+1] = {0};

//files
char* A_file = "a_i_j.txt";
char* B_file = "b_i_j.txt";
char* PI_file = "pi.txt";
char *avg_a_file = "avg_A.txt";
char *avg_b_file = "avg_B.txt";
char *avg_pi_file = "avg_Pi.txt";

string alpha_file = "alpha";
string beta_file = "beta";
string gamma_file = "gamma";
string qstart_file = "state_seq_viterbi.txt";

int cnt = 1;
long double P_Obs_for_Model = 0;

//Calculation of alpha variable to find the solution of problem number 1.
void forward_procedure(int iteration){
	int i , j , t;
	long double sum ;
	int index = O[1];
	P_Obs_for_Model = 0;

	for(i=1;i<=N;i++){
		Alpha[1][i] = Pi[i]*B[i][index];
	}
	
	for (t = 1; t < T; t++){
		for (j = 1; j <= N; j++){
			sum = 0;
			for (i = 1; i <= N; i++){
				sum += Alpha[t][i] * A[i][j];
			}
			Alpha[t + 1][j] = sum * B[j][O[t + 1]];
		}
	}
	for(i=1;i<=N;i++){
		P_Obs_for_Model = P_Obs_for_Model + Alpha[T][i];
	}
	// prob_seq[digit_count]=P_Obs_for_Model;
	//cout << P_Obs_for_Model <<"\n";
	// if(P_Obs_for_Model > max_pobs_model){
	// 	max_pobs_model = P_Obs_for_Model;
	// 	test_ans = iteration;
	// }

	// cout << "Digit:"<<iteration<<"\tP(obs/model) : " << P_Obs_for_Model <<endl;
}


void solution_to_prob1(int iteration){
	forward_procedure(iteration);
}

//Calculation of Beta variable.
void backward_procedure(){
	int i , j , t;
	long double sum;
	int index = 0;
	for(i=1;i<=N;i++){
		Beta[T][i] = 1.0;
	}
	for(t=T-1;t>=1;t--){
		index = O[t+1];
		for(i=1;i<=N;i++){
			sum = 0;
			for(j=1;j<=N;j++){
				sum = sum + B[j][index]*A[i][j]*Beta[t+1][j];
			}
			Beta[t][i]=sum;
		}
	}
}

//calculating gamma values using xita
void calculate_gamma_values(){
	int i, j, t;
	for (t = 1; t <= T - 1; t++){
		for (i = 1; i <= N; i++){
			Gamma[t][i] = 0;
			for (j = 1; j <= N; j++){
				Gamma[t][i] += Xi[t][i][j];
			}
		}
	}
}

//calculating alpha and beta values
void calculate_gamma(){
	for(int t=1;t<=T;t++){
		for(int j=1;j<=N;j++){
			long double summation=0;
			for(int k=1;k<=N;k++){
				summation += Alpha[t][k] * Beta[t][k];
			}
			Gamma[t][j]=(Alpha[t][j] * Beta[t][j])/summation;
		}
	}
}

//loading the model parameters with new calculated values
void load_calculated_model(){
	int i, j;
	for(i=1;i<=N;i++){
		Pi[i]=Pibar[i];
	}
	
	for(i=1;i<=N;i++){
		for(j=1;j<=N;j++){
			A[i][j]= Abar[i][j];
		}
	}
	for(i=1;i<=N;i++){
		for(j=1;j<=M;j++){
			B[i][j] = Bbar[i][j];
		}
	}
}

//TO APPLY THRESHOLD VALUE TO MAKE VALUES NON-ZERO
void apply_threshold_to_Bij(){
	int i, j;
	long double diff;
	long double max;
	int max_i=0;
	int max_j=0;
	for (i = 1; i <= N; i++){
		diff = 0;
		max = 0;
		for (j = 1; j <= M; j++){
			if (Bbar[i][j] > max){
				max = Bbar[i][j];
				max_i = i;
				max_j = j;
			}
			if (Bbar[i][j] < threshold){
				diff += Bbar[i][j] - threshold;
				Bbar[i][j] = threshold;
			}
		}
		Bbar[max_i][max_j] = max;
	}
}

void reevaluate_model_parameters(){
	int i, j, k, t;
	long double sum1=0 , sum2 =0;
	//Re-evaluating Pi
	for(i=1;i<=N;i++){
		Pibar[i] = Gamma[1][i];
	}
	
	for(i = 1; i<=N; i++){
		for(j = 1; j <= N; j++){
			long double t1 = 0, t2 = 0;
			for(t = 1; t <= T-1; t++){
				t1 += Xi[t][i][j];
				t2 += Gamma[t][i];
				//cout<<"Xi["<<t<<"]["<<i<<"]["<<j<<"]: "<<Xi[t][i][j]<<", Gamma["<<t<<"]["<<j<<"]: "<<Gamma[t][j]<<endl;
			}
			//cout<<"t1 "<<t1<<" t2: "<<t2<<endl;
			//cout<<"t1/t2: "<<t1/t2<<endl;
			//system("pause");
			Abar[i][j] = t1/t2;
			// cout<<"Abar: "<<Abar[i][j]<<endl;
		}
	}
	//Re-evaluating B
	for(j=1;j<=N;j++){
		int count=0;
		long double max=0;
		int ind_j=0, ind_k=0;
		
		for(k=1;k<=M;k++){
			sum1 =0 , sum2 =0;
			for(t=1;t<T;t++){
				sum1 = sum1 + Gamma[t][j];
				if(O[t]==k){
					sum2 = sum2 + Gamma[t][j];				
				}
			}
			Bbar[j][k] = sum2/sum1;
			
			//finding max
			if(Bbar[j][k]>max){
				max=Bbar[j][k];
				ind_j = j;
				ind_k = k;
			}
			
			//updating new bij with threshold value if it is zero
			if(Bbar[j][k] == 0){
				Bbar[j][k]=threshold;
				count++;
			}
		}
		Bbar[ind_j][ind_k] = max - count*threshold;
	}
	load_calculated_model();
}

//Adjusting Model Parameters
void calculate_xi(){

	int i , j , t , index;
	long double summation[FRAMES + 1];

	for(t=1;t<=T;t++){
		// index = ;
		summation[t] = 0;
		for(i=1;i<=N;i++){
			for(j=1;j<=N;j++){
				summation[t] += Alpha[t][i]*A[i][j]*B[j][O[t+1]]*Beta[t+1][j];
			}
		}

		for(i=1;i<=N;i++){
			long double x;
			for(j=1;j<=N;j++){
				x = Alpha[t][i]*A[i][j]*B[j][O[t+1]]*Beta[t+1][j];
				Xi[t][i][j]= x/summation[t];
			}
		}
	}
}

//viterbi algorithm
void viterbi(){
    for(int i=1; i<=N; i++){
        Delta[1][i] = Pi[i] * B[i][O[1]];
        Psi[1][i] = 0;
    }

	for(int j=1; j<=N; j++){
		for(int t=2; t<=T; t++){
            long double max = 0, ti = 0;
            int ind = 0;
            
            for(int i=1; i<=N; i++){
                ti = Delta[t-1][i] * A[i][j];
                if(ti > max){
					max = ti;
					ind = i;
				}
            }

            Delta[t][j] = max * B[j][O[t]];
			Psi[t][j] = ind;
        }
    }

    long double max = 0;
    for(int i=1; i<=N; i++){
        if(Delta[T][i] > max) {
			max = Delta[T][i];
			Q[T] = i;
		}

        pstar = max;
    }

    for(int t = T-1; t>0; t--){
        Q[t] = Psi[t+1][Q[t+1]];
    }
}

//writing updated A matrix to file
void write_final_A_matrix(FILE *fp){
	// FILE *fp;
	// fp = fopen(filename,"w");
	int i, j;
	fprintf(fp, "---------------A Matrix----------------\n");
	for (i = 1; i <= N; i++){
		for (j = 1; j <= N; j++){
			fprintf(fp,"%Le   ",A[i][j]);
		}
		fprintf(fp,"\n");
	}
	// fclose(fp);
}

//writing updated B matrix to file
void write_final_B_matrix(FILE *fp){
	// FILE *fp;
	// fp = fopen(filename, "w");
	int i, j;
	fprintf(fp, "---------------B Matrix----------------\n");
	for (i = 1; i <= N; i++){
		for (j = 1; j <= M; j++){
			fprintf(fp, "%Le   ", B[i][j]);
		}
		fprintf(fp, "\n");
	}
	// fclose(fp);
}

//writing updated pi values to file
void write_final_pi_matrix(FILE *fp){
	// FILE *fp;
	// fp = fopen(filename, "w");
	fprintf(fp, "---------------Pi values----------------\n");
	int i, j;
	for (i = 1; i <= N; i++){
		fprintf(fp, "%Le   ", Pi[i]);
	}
	// fclose(fp);
}

//dump the model
void dump_converged_model(FILE *fp){
	write_final_A_matrix(fp);
	write_final_B_matrix(fp);
	write_final_pi_matrix(fp);
}

//read A
bool readA(char *filename){
    fstream fin;
	fin.open(filename);    

    //file does not exist
	if(!fin){
		cout<<"Couldn't open file: "<<filename<<"\n";
		return false;
	}
	long double word;

    int row = 1, col = 1;
    //until input is available
	while(fin >> word){
        col = 1;
        A[row][col++] = word;

        for(int i=2; i<=N; i++){
            fin>>word;
            A[row][col++] = word;
        }
        row++;
    }

	fin.close();
	return true;
}

//read B
bool readB(string filename){
	fstream fin;
	fin.open(filename);    

    //file does not exist
	if(!fin){
		cout<<"Couldn't open file: "<<filename<<"\n";
		return false;
	}
	long double words;

    int row = 1, col = 1;

	while(fin>>words){
		col = 1;
		B[row][col++] = words;

		for(int i=1; i<M; i++){
			fin>>words;
			B[row][col++] = words;
		}
		row++;
	}
	//cout<<"row: "<<row<<endl;
	fin.close();
	return true;
}

//read Pi
bool readPi(string filename){
	fstream fin;
	fin.open(filename);    

    //file does not exist
	if(!fin){
		cout<<"Couldn't open file: "<<filename<<"\n";
		return false;
	}
	long double word;

    int col = 1;
    //until input is available
	while(fin >> word){
		col = 1;
        Pi[col++] = word;

        //save whole row
		for(int i=1;i<N;i++){
			fin>>word;
            Pi[col++] = word;
		}
	}

	fin.close();
	return true;
}

//saving the model
void save(ofstream &op, char digit, FILE *fp){
	// alpha_file += digit +"_" + to_string(cnt) + ".txt";
	// alpha_file.append(digit);
	// beta_file += digit +"_" + to_string(cnt) + ".txt";
	// gamma_file += digit +"_" + to_string(cnt++) + ".txt";
	// char alphafn[100], betafn[100], gammafn[100];
	// sprintf(alphafn, "alpha_%c_%d.txt", digit, cnt);
	// sprintf(betafn, "beta_%c_%d.txt", digit, cnt);
	// sprintf(gammafn, "gamma_%c_%d.txt", digit, cnt++);

	// if(cnt == 6) cnt = 0;

	//solution to prob 1
	// ofstream op(alphafn);
	// op<<"\n\n------------------Solution to problem 1------------------\n";
	// cout<<"Alpha - ";
	// op<<"Probability: "<<P_Obs_for_Model<<endl;
	// if(!op) cout<<"Problem with "<<alpha_file<<endl;
	
	// op<<"----------------------Forward (alpha matrix)----------------------------\n";
	// for(int i=1; i<=T; i++){
	// 	for(int j=1; j<=N; j++){
	// 		op<<setw(15)<<Alpha[i][j]<<"   ";
	// 		// cout<<Alpha[i][j]<<"   ";
	// 	}
	// 	op<<endl;
	// 	// cout<<endl;
	// }
	
	// cout<<"Alpha written to file "<<endl;

	// cout<<"Beta - ";
	
	// op<<"----------------------Backward (beta matrix)----------------------------\n";
	// for(int i=1; i<=T; i++){
	// 	for(int j=1; j<=N; j++){
	// 		// cout<<Beta[i][j]<<"   ";
	// 		op<<setw(15)<<Beta[i][j]<<"   ";
	// 	}
	// 	// cout<<endl;
	// 	op<<endl;
	// }
	
	// cout<<"Beta written to file "<<endl;


	// op<<"\n\n--------------------------Solution to problem 2-------------------------\n";
	// cout<<"Probability : "<<pstar<<endl;
	fprintf(fp, "Probability P*: %Le\n",pstar);

	fprintf(fp, "obsseq: ");
	for(int i=1;i<T;i++){
		// cout<<O[i]<<" ";
		// op<<setw(3)<<O[i]<<"   ";
		fprintf(fp, "%3d   ", O[i]);
	}
	// op<<endl;
	// fprintf(fp, "\n");
	
	// op<<"qstar : ";
	fprintf(fp, "\nq star: ");
	for(int i=1;i<=T;i++){
		// cout<<Q[i]<<" ";
		// op<<setw(3)<<Q[i]<<"   ";
		fprintf(fp, "%3d   ", Q[i]);
	}
	fprintf(fp, "\n\n");
	// op<<"\n\n--------------------------Solution to problem 3-------------------------\n";
	// cout<<"Gamma - ";
	// // op.open(gammafn);
	// for(int i=1; i<=T; i++){
	// 	for(int j=1; j<=N; j++){
	// 		// cout<<Gamma[i][j]<<"   ";
	// 		op<<setw(15)<<Gamma[i][j]<<"   ";
	// 	}
	// 	// cout<<endl;
	// 	op<<endl;
	// }
	
	// op<<endl;
	// cout<<"Gamma written to file "<<endl;
	int i, j;
	// cout<<"---------------------Pi---------------------------\n";
	// op<<"---------------------Pi---------------------------\n";
	// for(i=1;i<=N;i++){
	// 	// Pi[i]=Pibar[i];
	// 	// cout<<Pibar[i]<<"   ";
	// 	op<<Pibar[i]<<"   ";
	// }
	// // cout<<endl;
	// op<<endl;
	// op<<"-------------------A matrix---------------------\n";
	// cout<<"-------------------A matrix---------------------\n";
	// for(i=1;i<=N;i++){
	// 	for(j=1;j<=N;j++){
	// 		// A[i][j]= Abar[i][j];
	// 		// cout<<Abar[i][j]<<"   ";
	// 		op<<Abar[i][j]<<"   ";
	// 	}
	// 	// cout<<endl;
	// 	op<<endl;
	// }
	// // cout<<"-------------------B matrix----------------------\n";
	// op<<"-------------------B matrix----------------------\n";
	// fprintf(fp, "-------------------B matrix----------------------\n");
	// for(i=1;i<=N;i++){
	// 	for(j=1;j<=M;j++){
	// 		// B[i][j] = Bbar[i][j];
	// 		// cout<<Bbar[i][j]<<"   ";
	// 		op<<Bbar[i][j]<<"   ";
	// 		fprintf(fp, "%Le   ", Bbar[i][j]);
	// 	}
	// 	fprintf(fp, "\n");
	// 	// cout<<endl;
	// 	op<<endl;
	// }
	cout<<"New model (PI, A, B) written to file"<<endl;
}

// make the model values, average model values and bar model values -  0
void erase_all_model(){
	for(int i=1; i<=N; i++){
		for(int j=1; j<=N; j++){
			A[i][j] = 0;
			a_average[i][j] = 0;
			Abar[i][j] = 0;
		}
	}

	for(int i=1; i<=N; i++){
		for(int j=1; j<=M; j++){
			B[i][j] = 0;
			b_average[i][j] = 0;
			Bbar[i][j] = 0;
		}
	}

	for(int i=1; i<=N; i++){
		Pi[i] = 0;
		Pibar[i] = 0;
		pi_average[i] = 0;
	}
}

//erasing the current model
void erase_model(){
	for(int i=1; i<=N; i++){
		for(int j=1; j<=N; j++){
			A[i][j] = 0;
		}
	}

	for(int i=1; i<=N; i++){
		for(int j=1; j<=M; j++){
			B[i][j] = 0;
		}
	}

	for(int i=1; i<=N; i++){
		Pi[i] = 0;
	}
}

// erasing average model
void erase_avg_model(){
	for(int i=1; i<=N; i++){
		for(int j=1; j<=N; j++){
			a_average[i][j] = 0;
		}
	}

	for(int i=1; i<=N; i++){
		for(int j=1; j<=M; j++){
			b_average[i][j] = 0;
		}
	}

	for(int i=1; i<=N; i++){
		pi_average[i] = 0;
	}
}

//initialize model according to parameters
void initialize_model(int digit, int seq, char *filename = "--"){
	char a_file[100], b_file[100], pi_file[100], obs_file[100];

	if(filename == "--"){
		readA(A_file);
		readB(B_file);
		readPi(PI_file);
	}else if(filename  == "avg"){
		sprintf(a_file, "models/model_A_%d_%d.txt", digit, seq);
		sprintf(b_file, "models/model_B_%d_%d.txt", digit, seq);
		sprintf(pi_file, "models/model_pi_%d_%d.txt", digit, seq);
		
		readA(a_file);
		readB(b_file);
		readPi(pi_file);
	}else if(filename == "init"){
		sprintf(a_file, "Model_Obs_Seq/Digit %d/A_%d.txt", digit, digit);
		sprintf(b_file, "Model_Obs_Seq/Digit %d/B_%d.txt", digit, digit);
		sprintf(pi_file, "Model_Obs_Seq/Digit %d/pi_%d.txt", digit, digit);

		readA(a_file);
		readB(b_file);
		readPi(pi_file);
	}
}

//adding current model values to avg model
void add_to_avg_model(){
	int i, j;
	for (i = 1; i <= N; i++){
		for (j = 1; j <= N; j++){
			a_average[i][j] += A[i][j];
		}
	}
	for (i = 1; i <= N; i++){
		pi_average[i] += Pi[i];
	}
	for (int i = 1; i <= N; i++){
		for (int j = 1; j <= M; j++){
			b_average[i][j] += B[i][j];
		}
	}
}

//writing updated a values
void write_final_A_values(char filename[]){
	FILE *fp = fopen(filename, "w");
	int i, j;
	for (i = 1; i <= N; i++){
		for (j = 1; j <= N; j++){
			// out << a_i_j[i][j] << " ";
			fprintf(fp, "%Le   ", A[i][j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

//print alpha beta gamma on screen
void print(){
	cout<<"Alpha values - \n";
	for(int i = 1; i<=T; i++){
		for(int j = 1; j<=N; j++){
			cout<<Alpha[i][j]<<"   ";
		}
		cout<<endl;
	}

	cout<<"Beta values -\n";
	for(int i = 1; i<=T; i++){
		for(int j = 1; j<=N; j++){
			cout<<Beta[i][j]<<"   ";
		}
		cout<<endl;
	}

	cout<<"Gamma values -\n";
	for(int i = 1; i<=T; i++){
		for(int j = 1; j<=N; j++){
			cout<<Gamma[i][j]<<"   ";
		}
		cout<<endl;
	}
}

//writing updated b values
void write_final_B_values(char filename[]){
	// out.open(filename);
	ofstream out(filename);
	for(int i=1; i<=N; i++){
		for(int j=1; j<=M; j++){
			out<<B[i][j]<<"   ";
		}
		out<<endl;
	}
	out.close();
}

//writing updated pi values
void write_final_pi_values(char filename[]){
	FILE *fp = fopen(filename, "w");
	int i, j;
	for (i = 1; i <= N; i++){
		// out << pi[i] << " ";
		fprintf(fp, "%Le   ", Pi[i]);
	}
	// out.close();
	fclose(fp);
}

// dump the final model in output/digit d/ folder wise
void dump_final_model(int seq, int digit){
	char index[10];
	char a_file_final[100], b_file_final[100], pi_file_final[100];

	sprintf(a_file_final, "output/digit %d/model_A_%d.txt", digit, seq);
	write_final_A_values(a_file_final);

	sprintf(b_file_final, "output/digit %d/model_B_%d.txt", digit, seq);
	write_final_B_values(b_file_final);

	sprintf(pi_file_final, "output/digit %d/model_PI_%d.txt", digit, seq);
	write_final_pi_values(pi_file_final);

	cout<<"digit "<<digit<<", sequence : "<<seq<<" model dumped successfully\n";
}

// taking average of the avg model
void average_of_avg_model(int total_iterations){
	int i, j;
	for (i = 1; i <= N; i++){
		for (j = 1; j <= N; j++){
			a_average[i][j] /= total_iterations;

		}
	}
	for (i = 1; i <= N; i++){
		for (j = 1; j <= M; j++){
			b_average[i][j] /= total_iterations;
		}
	}
	for (i = 1; i <= N; i++){
		pi_average[i] /= total_iterations;
	}
}

// dumping average model to file
void dump_avg_model(int digit){
	char a_file_avg[100], b_file_avg[100], pi_file_avg[100], ind[3];

	sprintf(a_file_avg, "output/avgmodels/digit_%d_A.txt", digit);
	FILE *fp = fopen(a_file_avg, "w");
	for(int i=1; i<=N; i++){
		for(int j=1; j<=N; j++){
			fprintf(fp, "%Le   ", a_average[i][j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	
	sprintf(b_file_avg, "output/avgmodels/digit_%d_B.txt", digit);
	ofstream fout(b_file_avg);
	for(int i=1; i<=N; i++){
		for(int j=1; j<=M; j++){
			//fprintf(fp, "%Le   ", b_average[i][j]);
			fout<<b_average[i][j]<<"   ";
		}
		fout<<endl;
		//fprintf(fp, "\n");
	}
	fout.close();

	sprintf(pi_file_avg, "output/avgmodels/digit_%d_PI.txt", digit);
	fp = fopen(pi_file_avg, "w");
	for(int i=1; i<=N; i++){
		fprintf(fp, "%Le   ", pi_average[i]);
	}
	fclose(fp);
}

//finding dc shift
void get_DC_shift(){
	long int sample_count = 0;
    FILE *fp;
    char line[80];

    //reading dc_shift.txt file
    fp = fopen("silence.txt", "r");
    
    if(fp == NULL){
        printf("Silence File not found\n");
        exit(1);
    }
    
	dcShift = 0;
    while(!feof(fp)){
        fgets(line, 80, fp);
        dcShift += atof(line);
        sample_count++;
    }
    dcShift /= sample_count;
    
    fclose(fp);

}

//function to setup the global variable like, max and nFactor
//max and nFactor depends on the vowel recording file and are used to do the normalization
void setupGlobal(char *filename){
    FILE *fp;
    long int totalSample = 0;
    char line[100];

    fp = fopen(filename, "r");
    if(fp == NULL){
        printf("Error opening file\n");
    }

    //get max value
    mx = 0;
    while(!feof(fp)){
        fgets(line, 100, fp);
        if(!isalpha(line[0])){
            totalSample++;
            if(mx < abs(atoi(line)))
                mx = abs(atoi(line));
        }
    }
    
    nFactor = (double)LIMIT/mx;
   
    //setup dcShift
    get_DC_shift();
    fclose(fp);
}

//Calculating Tokhura's Distance Using Code Book
void calculate_tokhura_distance(long double c[12], int index, FILE *fp){
	int  min_index = 0;
	long double min = DBL_MAX;
	long double sum[CB_SIZE] = { 0 };
	string temp, temp1;

	for (int j = 0; j < CB_SIZE; j++){
		for (int i = 0; i < P; i++){
			sum[j] += tokhuraWeights[i] * (c[i] - codeBook[j][i])*(c[i] - codeBook[j][i]);
		}
		if (sum[j] < min){
			min = sum[j];
			min_index = j;
		}
	}

	O[index] = min_index + 1;

	//cout << O[index] << " ";
	fprintf(fp, "%4d ", O[index]);
}

//This function calulate the cepstral coefficients Ci's
void calculate_Cis(){
	double sum=0;
	Ci[0]=log(Ri[0]*Ri[0]);

	for(int m=1;m<=P;m++){
		sum=0;
		for(int k=1;k<m;k++){
			sum += (k*Ci[k]*Ai[m-k])/(m*1.0);
		}
		Ci[m]=Ai[m]+sum;
	}
}

// Function to apply Durbin Algorithm And Find The value of ai's 
void durbinAlgo(){
	double alpha[13][13],E[13],K[13];
	double sum = 0;
	E[0] = Ri[0];
	//loop for p from 1 to 12
	for(int i=1;i<=P;i++){
		sum=0;
		for(int j=1;j<=i-1;j++){
			sum += alpha[i-1][j]*Ri[i-j];	
		}
		
		K[i]=(Ri[i]-sum)/E[i-1];
			
		alpha[i][i]=K[i];
	
		for(int j=1;j<=i-1;j++){
			alpha[i][j]=alpha[i-1][j] - K[i]*alpha[i-1][i-j];
		}
	
		E[i]=(1-(K[i]*K[i]))*E[i-1];
	}

	//storing the ai values
	for(int i=1;i<=P;i++){
		Ai[i]= alpha[P][i];
	}	
}

//calculating the Ris values
void calculate_Ris(double *samp){
	for(int m =0; m<=P; m++){
		Ri[m] = 0;
		for(int k=0; k<N-m; k++){
			Ri[m] += samp[k]*samp[k+m];
		}
	}
}

//function to apply the Raised Sine Window in Ci of each frame
void applyingRaisedSinWindow(){
	long double sum=0;
	for(int m=1;m<=P;m++){
		sum = (P/2)*sin((PI*m)/P);
		Ci[m]*=sum;	
	}	
}

//calculating c prime values
void calculate_c_prime(double *samp){
	calculate_Ris(samp);
	//calling durbinAlgo to find ai values
	durbinAlgo();
	//finding cepstral constants
	calculate_Cis();
	//applying raised sin window on cis
	applyingRaisedSinWindow();
}

double* trim_waveform1(){
	int i;
	long double total_energy, current_value;
	char trimmed_file[100];
	int start_marker = 0, end_marker = 0;
	
	total_energy = 0;
	int interval = WIN_SIZE, index_count = 0, arr_index = 0;
	double samp[WIN_SIZE] = {0};
	long double total_max_energy = 0;
	int initial_shift_pointer = 0, shift_pointer_count = 0, writing_count = 0;
	long double arr_energy[50];
	long count1 = 0;
	string temp = "";
	//cout << "\n................Trimming Waveform..................." << endl;
	int count = interval;
	//in.open(normalized_file);
	ofstream out("energy_trim.txt");

	//window shifting logic
	for (i = 0; i < sSize; i++){
		count1++;
		//if (shift_pointer_count >= initial_shift_pointer){

		//current_value = stold(temp);
		current_value = sample[i];
		current_value *= 0.1;
		total_energy += (current_value*current_value);
		count--;
		if (!count){

			out << to_string(total_energy) << endl;

			if (total_energy > total_max_energy){
				total_max_energy = total_energy;
				start_marker = initial_shift_pointer;
				end_marker = initial_shift_pointer + interval;

			}

			total_energy = 0;
			count = interval;
			initial_shift_pointer += 1000;
			i = initial_shift_pointer - 1;
			shift_pointer_count = 0;
		}
	}
	out.close();
	
	//cout << "Start Marker : " << start_marker << endl;
	//cout << "End Marker : " << end_marker << endl;

	if(end_marker > sSize){
		start_marker -= (end_marker - sSize);
		end_marker = sSize;
	}
	
	
	sprintf(trimmed_file, "trim.txt");
	ofstream out1(trimmed_file);
	for(int i=0; i<sSize; i++){
		index_count++;
		current_value = sample[i]; //stold(temp);
		//current_value = current_value - dc_shift_value;
		//current_value = current_value * normalization_ratio;

		if (index_count >= start_marker && index_count < end_marker){
			samp[arr_index++] = current_value;
			out1 << samp[arr_index - 1] << endl;
		}
	}

	out1.close();
	return samp;
}

//generate observation sequence
void generate_obs_sequence(char *filename){
	int obs_ind = 1;
	// char obs_file[100];
	// //framing
	// sprintf(obs_file, "output/live_test/obs_seq.txt", d, u);
	FILE *op = fopen(filename, "a+");
	if(op == NULL) {
		printf("locha hai idhar bhaiiiii.. \n");
		exit(1);
	}
	
	double *frame_samples = trim_waveform1();
	double fsamp[FRAME_SIZE];

	for(int i=0; i<WIN_SIZE; i++){
		for(int j = 0; j<320; j++)
			fsamp[j] = frame_samples[i++];

		calculate_c_prime(fsamp);
		calculate_tokhura_distance(Ci, obs_ind++, op);
	}
	fprintf(op, "\n");
	fclose(op);
	cout<<"wrote observation seq in file: "<<filename<<"\n";
}

//generate observation sequence
void generate_obs_sequence(char *filename, FILE *op){
	int obs_ind = 1;
	// char obs_file[100];
	// //framing
	// sprintf(obs_file, "output/live_test/obs_seq.txt", d, u);
	// FILE *op = fopen(filename, "w");
	// if(op == NULL) {
	// 	printf("locha hai idhar bhaiiiii.. \n");
	// 	exit(1);
	// }
	
	double *frame_samples = trim_waveform1();
	double fsamp[FRAME_SIZE];

	for(int i=0; i<WIN_SIZE; i++){
		for(int j = 0; j<320; j++)
			fsamp[j] = frame_samples[i++];

		calculate_c_prime(fsamp);
		calculate_tokhura_distance(Ci, obs_ind++, op);
	}
	cout<<"wrote observation seq in file: "<<filename<<"\n";
}

void train_file(char *filename, int digit){
	char line[100], obs_file[100];
	// for(int d = 0; d<=1; d++){
		erase_model();
		
		// for(int u = 1; u<=20; u++){
		// sprintf(filename, "input/recordings/Digit %d/rec_%d.txt", d, u);

		FILE *f = fopen(filename, "r");

		if(f == NULL){
			printf("Issue in opening file %s", filename);
			exit(1);
		}
		
		//setting dcshift and nfactor
		setupGlobal(filename);

		sSize = 0;
		//reading the samples and normalizing them
		while(!feof(f)){
			fgets(line, 100, f);
			
			//input file may contain header, so we skip it
			if(!isalpha(line[0])){
				int y = atof(line);
				double normalizedX = floor((y-dcShift)*nFactor);
				//if(abs(normalizedX) > 1)
				sample[sSize++] = normalizedX;
			}
		}
		fclose(f);

		//framing
		sprintf(obs_file, "output/delete/digit_%d_obs_seq.txt", digit);
		generate_obs_sequence(obs_file);
		initialize_model(digit, 1, "--");

		int iteration = 1;
		while(iteration <= CONVERGE_ITERATIONS){
			//cout<<"iteration: "<<iteration++<<endl;
			iteration++;
			forward_procedure(0);
			backward_procedure();
			viterbi();
			calculate_xi();
			calculate_gamma();
			reevaluate_model_parameters();
			//print();
		}

		char a_file_final[100], b_file_final[100], pi_file_final[100];

		sprintf(a_file_final, "output/delete/digit_%d_model_A.txt", digit);
		write_final_A_values(a_file_final);
		sprintf(b_file_final, "output/digit_%d_model_B.txt", digit);
		write_final_B_values(b_file_final);
		sprintf(pi_file_final, "output/digit_%d_model_PI.txt", digit);
		write_final_pi_values(pi_file_final);

		system("pause");
	// }
}

void make_all_obs_seq(){
	char filename[100], line[100], obs_file[100];
	for(int d = 0; d<=9; d++){

		sprintf(obs_file, "output/obs_seq/HMM_OBS_SEQ_%d.txt", d);
		FILE *op = fopen(obs_file, "w");
		
		for(int u = 1; u<=20; u++){
			sprintf(filename, "input/recordings/Digit %d/rec_%d.txt", d, u);

			FILE *f = fopen(filename, "r");

			if(f == NULL){
				printf("Issue in opening file %s", filename);
				exit(1);
			}
			
			//setting dcshift and nfactor
			setupGlobal(filename);

			sSize = 0;
			//reading the samples and normalizing them
			while(!feof(f)){
				fgets(line, 100, f);
				
				//input file may contain header, so we skip it
				if(!isalpha(line[0])){
					int y = atof(line);
					double normalizedX = floor((y-dcShift)*nFactor);
					//if(abs(normalizedX) > 1)
					sample[sSize++] = normalizedX;
				}
			}
			fclose(f);
			//calling helper function to generate observation seq
			generate_obs_sequence(obs_file, op);
			fprintf(op, "\n");
		}
		fclose(op);
	}
}

void training(){
	char filename[100], line[100], obs_file[100];
	erase_all_model();
	for(int d = 0; d<=9; d++){
		erase_model();

		// sprintf(obs_file, "output/obs_seq/HMM_OBS_SEQ_%d.txt", d);
		// FILE *op = fopen(obs_file, "w");
		
		for(int u = 1; u<=20; u++){
			sprintf(filename, "input/recordings/Digit %d/rec_%d.txt", d, u);

			FILE *f = fopen(filename, "r");

			if(f == NULL){
				printf("Issue in opening file %s", filename);
				exit(1);
			}
			
			//setting dcshift and nfactor
			setupGlobal(filename);

			sSize = 0;
			//reading the samples and normalizing them
			while(!feof(f)){
				fgets(line, 100, f);
				
				//input file may contain header, so we skip it
				if(!isalpha(line[0])){
					int y = atof(line);
					double normalizedX = floor((y-dcShift)*nFactor);
					//if(abs(normalizedX) > 1)
					sample[sSize++] = normalizedX;
				}
			}
			fclose(f);

			//framing
			// double frame_sample[320];
			// int obs_ind = 1;
			// if(op == NULL) {
			// 	printf("locha hai idhar bhaii.. \n");
			// 	exit(1);
			// }
			// // for(int i=1,k=0;k<=(240*85) && k<sSize-240;i++,k+=240){ //calculates ci's for 320 samples frame
			// // 	num_frames++;
			// // 	for(int lmk = 0; lmk<320; lmk++)
			// // 		samp[lmk] = sample[lmk+k];
			// // 	//tn=start+k;
			// // 	// calculate_ri(samp,r);
 			// // 	// durbin(r,a);
 			// // 	// cepstral(a,c);
			// // 	calculate_c_prime(samp);

			// // 	calculate_tokhura_distance(Ci, obs_ind++, op);
			// // }

			// double *frame_samples = trim_waveform1();
			// double fsamp[FRAME_SIZE];

			// for(int i=0; i<WIN_SIZE; i++){
			// 	for(int j = 0; j<320; j++)
			// 		fsamp[j] = frame_samples[i++];

			// 	calculate_c_prime(fsamp);
			// 	calculate_tokhura_distance(Ci, obs_ind++, op);
			// }
			// fclose(op);
			// cout<<"wrote output/obs_seq/digit_"<<d<<"_obs_"<<u<<".txt file\n";
			sprintf(obs_file, "output/obs_seq/HMM_OBS_SEQ_%d_%d.txt", d, u);
			generate_obs_sequence(obs_file);
			// fprintf(op, "\n");
			// calculate_c_prime(samp);
			// calculate_tokhura_distance(Ci, obs_ind++, op);
			// fclose(op);


			initialize_model(d, 1, "--");

			int iteration = 1;
			while(iteration <= CONVERGE_ITERATIONS){
				//cout<<"iteration: "<<iteration++<<endl;
				iteration++;
				forward_procedure(0);
				backward_procedure();
				viterbi();
				calculate_xi();
				calculate_gamma();
				reevaluate_model_parameters();
				//print();
			}

			add_to_avg_model();
			dump_final_model(u, d);

			//cout<<"number of frames: "<<obs_ind<<endl;
			// system("pause");
		}
		// fclose(op);
		average_of_avg_model(20);
		dump_avg_model(d);
		erase_avg_model();
		
		system("pause");
	}
}

//TO READ CODEBOOK FROM FILE
void read_codebook(){
	ifstream in("provided_codebook.txt");
	for (int i = 0; i < CB_SIZE; i++){
		for (int j = 0; j < P; j++){
			in >> codeBook[i][j];
		}
	}
	in.close();
}

void read_average_model(int iteration){
	
	char filename[100];
	sprintf(filename, "output/avgmodels/digit_%d_A.txt", iteration);
	readA(filename);

	sprintf(filename, "output/avgmodels/digit_%d_B.txt", iteration);
	readB(filename);

	sprintf(filename, "output/avgmodels/digit_%d_PI.txt", iteration);
	readPi(filename);
}

void test_prediction(){
	test_ans = 0;
	max_pobs_model = 0;
	for(int k = 0; k<=9; k++){
		read_average_model(k);
		solution_to_prob1(k);
	}

	printf("Detected digit %d\n", test_ans);
}

void live_testing(){
	char obs_file[100], line[100];
	printf("\n----------Live testing----------\n");

	system("Recording_Module.exe 3 input.wav input_file.txt");

	initialize_model(0, 0, "--");

	FILE *f = fopen("input_file.txt", "r");
	if(f == NULL){
		printf("Issue in opening file input_file.txt");
		exit(1);
	}

	//setting dcshift and nfactor
	setupGlobal("input_file.txt");

	sSize = 0;
	//reading the samples and normalizing them
	while(!feof(f)){
		fgets(line, 100, f);
		
		//input file may contain header, so we skip it
		if(!isalpha(line[0])){
			int y = atof(line);
			double normalizedX = floor((y-dcShift)*nFactor);
			//if(abs(normalizedX) > 1)
			sample[sSize++] = normalizedX;
		}
	}
	fclose(f);
	generate_obs_sequence("output/live_test/obs_seq.txt");
	
	test_prediction();
}


int _tmain(int argc, _TCHAR* argv[]){

	printf("--------------Recording silence-------------\n");
	//system("Recording_Module.exe 3 silence.wav silence_file.txt");

	read_codebook();
	training();
	// train_file("input/recordings/Digit 0/rec_11.txt", 0);
	//live_testing();
	//initialize_model('1', 1, "avg");
	
	return 0;
}
