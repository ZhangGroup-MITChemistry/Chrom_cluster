#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include <cmath>
using namespace std;


// ~~~~~~~~~~~~~~~~~ //
// Declare Constants //
// ~~~~~~~~~~~~~~~~~ //

const double K_ON=0.002;
const double K_OFF=0.04;
const int K=1;
const double THRES=5E-7;
const int BEAD=600;
const int CTCF=24;
const int COL=2;
const int T_TOT=1;
const double T_STEP=1E-2;
const int STEP=T_TOT/T_STEP;
const int BIND_CUTOFF=1;
const int EQUSTEP=4E6;
const int RECSTEP=4E4;
const int MAXSTEP=EQUSTEP+RECSTEP;
const int OUTGAP=20;

#define RES 1E8
static vector<vector<int> > cohes_lst;

struct vectorCol
{
    vector<double> tern_all;
    vector<vector<int> > zero_seq_iter;
};


// ~~~~~~~~~~~~~~~~~ //
// Declare Functions //
// ~~~~~~~~~~~~~~~~~ //

void main_proc(int** dna_seq_iter,vector<int> ctcf_lst,int CTCF,double **ctcf_perm,double **prob_mat_raw,int **iter_num);
void update_step(int** dna_seq_iter, int spot_t3, vector<vector<int> > zero_seq_iter);
int update_center_side(int** dna_seq_iter, int spot_t1, int tern_index, int center_type, int left_type, int right_type, int left_change, int right_change);
void update_cohes(int tern_index, int spot_t1, int side_type);
int update_side(int tern_index,int side_type);
void gillespie(vector<double> react_lst, int size, double *ti, int *step);
vectorCol ternry_all(int **seq_iter, int nbead, double **ctcf_perm);
void ternry_each(int center_type, int center_id, int left_type, int right_type, int nbead, double *tern, double **ctcf_perm);
void each_step(int center_type, int left_type, int right_type, double p_left, double p_right, double *tern);
int Inarray(int val, vector<int> vec);
int Index2d(int spot_cohes, vector<vector<int> > cohes_lst, int col);
int Index1d(int spot_cohes, vector<int> ctcf_lst);
double max_diff(double **lst_2d_1,double **lst_2d_2,int di);
double **readdblfile(const char* fi_name, int row, int col);
int **readintfile(const char* fi_name, int row, int col);
void print_vec(vector<vector<int> > cohes_lst);
void print_mat(double **mat, int size);
int rowoffile(const char* fi_name);


// ~~~~~~~~~~~~~~ //
// Main Functions //
// ~~~~~~~~~~~~~~ //

void main_proc(int** dna_seq_iter,vector<int> ctcf_lst,int CTCF,double **ctcf_perm,double **prob_mat_raw,int **iter_num)
{
	// int step = T_TOT/T_STEP;
	double time = 0.0;
	double time_store = 0.0;
	double count = iter_num[0][0];
	double tau = 0.0;
	int stepi = 0;
	int index_0 = 0;
	int index_1 = 0;
	int index_t_fst = 0;
	int index_t_lst = 0;

	int kk,ic,ip,iq;
	double diff=1.0;
	int num_cohes;

	vectorCol newVectors;

	double **prob_mat_last;
	prob_mat_last = (double**)malloc((CTCF)*sizeof(double*));
	for (int ii=0;ii<CTCF;ii++)
	{
		prob_mat_last[ii] = (double*)malloc((CTCF)*sizeof(double));
	}
	for (ip=0;ip<CTCF;ip++) {
		for (iq=0;iq<CTCF;iq++) {
			prob_mat_last[ip][iq] = prob_mat_raw[ip][iq];
		}
	}

	double **prob_mat_ta_all;
	prob_mat_ta_all = (double**)malloc((CTCF)*sizeof(double*));
	for (int jj=0;jj<CTCF;jj++)
	{
		prob_mat_ta_all[jj] = (double*)malloc((CTCF)*sizeof(double));
	}
	for (ip=0;ip<CTCF;ip++) {
		for (iq=0;iq<CTCF;iq++) {
			prob_mat_ta_all[ip][iq] = prob_mat_raw[ip][iq];
		}
	}

	double **prob_mat;
	prob_mat = (double**)malloc((CTCF)*sizeof(double*));
	for (kk=0;kk<CTCF;kk++)
	{
		prob_mat[kk] = (double*)malloc((CTCF)*sizeof(double));
	}


	ofstream fo_diff("diff_record.txt");
	ofstream fo_cohes("cohes_record.txt");
	ofstream fo_prob("prob_mat_raw_restart.txt");
	ofstream fo_dna_seq_iter("ctcf_position_restart.txt");
	ofstream fo_cohes_lst("cohes_lst_restart.txt");
	ofstream fo_iter_num("iter_num_restart.txt");


	while (count < MAXSTEP)
	{
		time = 0.0;
		time_store = 0.0;
		count += 1.0;

		fo_prob.clear();
		fo_prob.seekp(0);

		fo_dna_seq_iter.clear();
		fo_dna_seq_iter.seekp(0);

		// fo_cohes_lst.clear();
		// fo_cohes_lst.seekp(0);

		fo_iter_num.clear();
		fo_iter_num.seekp(0);


		for (ip=0;ip<CTCF;ip++) {
			for (iq=0;iq<CTCF;iq++) {
				prob_mat[ip][iq] = 0.0;
			}
		}

		while (time < T_TOT)
		{
			// zero_seq_iter.clear();
			newVectors = ternry_all(dna_seq_iter,BEAD,ctcf_perm);

			// for (int ip=1;ip<newVectors.zero_seq_iter.size();ip++)
			// {
			// 	cout << newVectors.zero_seq_iter[ip][0] << '\t' << newVectors.zero_seq_iter[ip][1] << '\n';
			// }

			gillespie(newVectors.tern_all, newVectors.tern_all.size(), &tau, &stepi);
			time += tau;

			int l = cohes_lst.size();
			for (ic=0;ic<l;ic++)
			{
				if ( Inarray(cohes_lst[ic][0], ctcf_lst) && Inarray(cohes_lst[ic][1], ctcf_lst) )
				{
					index_0 = Index1d(cohes_lst[ic][0],ctcf_lst);
					index_1 = Index1d(cohes_lst[ic][1],ctcf_lst);
					index_t_fst = (int) time_store / T_STEP;
					index_t_lst = (int) time / T_STEP;
					prob_mat[index_0][index_1] += (index_t_lst - index_t_fst);
				}
			}

			// cout << stepi << endl;
			time_store = time;
			update_step(dna_seq_iter,stepi,newVectors.zero_seq_iter);
		}

		for (ip=0;ip<CTCF;ip++){
			for (iq=0;iq<CTCF;iq++) {
				prob_mat[ip][iq] /= STEP;
			}
		}

		for (ip=0;ip<CTCF;ip++){
			for (iq=0;iq<CTCF;iq++){
				prob_mat_ta_all[ip][iq] = (prob_mat_ta_all[ip][iq]*(count-1.0) + prob_mat[ip][iq]) / count;
			}
		}

		// print_mat(prob_mat,CTCF);
		// cout << endl;
		// print_mat(prob_mat_ta_all,CTCF);
		// cout << endl;
		// print_mat(prob_mat_last,CTCF);
		// cout << endl;

		diff = max_diff(prob_mat_ta_all,prob_mat_last,CTCF);
		// cout << "here the diff is: " << diff << endl;

		if (diff != 0.0)
		{
			fo_diff << diff << endl;
			fo_diff.flush();

			num_cohes = cohes_lst.size();
			fo_cohes << num_cohes << endl;
			fo_cohes.flush();
			//cout << "The no. of iteration: " << count << endl;
			cout << count << endl;
			//cout << "The difference is: " << diff << endl;
			cout << diff << endl;
			cout << endl;

			if (diff < THRES){
				fo_diff << diff << endl;
				fo_diff.flush();
				// break;
			}
			else{
				// memcpy(prob_mat_last,prob_mat_ta_all,CTCF*CTCF*2);
				for (ip=0;ip<CTCF;ip++){
					for (iq=0;iq<CTCF;iq++){
						prob_mat_last[ip][iq] = prob_mat_ta_all[ip][iq];
					}
				}
			}

		}

		if ((int)count%OUTGAP == 0 && (int)count>EQUSTEP)
		{
			for (ip=0;ip<CTCF;ip++){
				for (iq=0;iq<CTCF;iq++){
					fo_prob << prob_mat_ta_all[ip][iq] << '\t';
				}
				fo_prob << endl;
			}
			fo_prob.flush();

			for (ip=1;ip<BEAD+1;ip++){
				for (iq=0;iq<2;iq++){
					fo_dna_seq_iter << dna_seq_iter[ip][iq] << '\t';
				}
				fo_dna_seq_iter << endl;
			}
			fo_dna_seq_iter.flush();

			for (ip=0;ip<cohes_lst.size();ip++){
				for (iq=0;iq<2;iq++){
					fo_cohes_lst << cohes_lst[ip][iq] << '\t';
				}
				fo_cohes_lst << endl;
			}
			fo_cohes_lst << endl;
			fo_cohes_lst << endl;
			
			fo_cohes_lst.flush();

			fo_iter_num << count;
			fo_iter_num.flush();

		}

	}

	for (ip=0;ip<CTCF;ip++){
		for (iq=0;iq<CTCF;iq++){
			fo_prob << prob_mat_ta_all[ip][iq] << '\t';
		}
		fo_prob << endl;
	}

	fo_diff.close();
	fo_cohes.close();
	fo_prob.close();
}


void update_step(int** dna_seq_iter, int spot_t3, vector<vector<int> > zero_seq_iter)
{
	if (spot_t3 < BEAD*3)
	{
		int spot_t1 = (int) spot_t3/3;
		spot_t1 += 1;
		int spot_t1_2 = -2;
		int tern_index = (spot_t3+1)%3;

		int center_type = dna_seq_iter[spot_t1][1];
		int left_change;
		int right_change;
		
		if (center_type == 0) {left_change=32; right_change=31;}
		else if (center_type == 4) {left_change=74; right_change=64;}
		else if (center_type == 1) {left_change=81; right_change=51;}
		else if (center_type == 2) {left_change=52; right_change=82;}
		else if (center_type == 3) {left_change=0; right_change=0;}
		else if (center_type == 31) {left_change=0; right_change=0;}
		else if (center_type == 32) {left_change=0; right_change=0;}
		else if (center_type == 51) {left_change=0; right_change=0;}
		else if (center_type == 52) {left_change=0; right_change=0;}
		else if (center_type == 64) {left_change=0; right_change=0;}
		else if (center_type == 74) {left_change=0; right_change=0;}
		else if (center_type == 81) {left_change=0; right_change=0;}
		else if (center_type == 82) {left_change=0; right_change=0;}
		else {cout << "Something is wrong." << endl;}

		int left_type = dna_seq_iter[spot_t1-1][1];
		int right_type = dna_seq_iter[spot_t1+1][1];

		spot_t1_2 = update_center_side(dna_seq_iter,spot_t1,tern_index,center_type,left_type,right_type,left_change,right_change);

		vector<int> unbind2_1 = {32};
		vector<int> unbind2_2 = {52,74,81};

		if (spot_t1_2 != -2)
		{
			int spot_type = dna_seq_iter[spot_t1_2][1];
			// cout << "spot_type of dna_seq_iter[spot_t1_2][1] is: " << spot_type << endl;

			if (Inarray(spot_type,unbind2_1))
			{
				dna_seq_iter[spot_t1_2][1] = 0;
			}
			else if (Inarray(spot_type,unbind2_2))
			{
				dna_seq_iter[spot_t1_2][1] = spot_type%10;
			}
			else
			{
				cout << "The second cohesin head is wrong." << endl;
			}
		}
		else;
	}

	else
	{
		// for (int ip=0;ip<zero_seq_iter.size();ip++) {cout << zero_seq_iter[ip][0] << '\t' << zero_seq_iter[ip][1] << endl;}
		// cout << zero_seq_iter[40][0] << endl;

		vector<int> temp;
		int fst,scd;
		int react_id;
		react_id = spot_t3 - BEAD*3;
		// cout << react_id << endl;

		fst = zero_seq_iter[react_id][0];
		scd = zero_seq_iter[react_id][1];

		temp.push_back(fst);
		temp.push_back(scd);
		cohes_lst.push_back(temp);
		temp.clear();

		dna_seq_iter[fst][1] = 31;
		dna_seq_iter[scd][1] = 32;
	}
	
}


int update_center_side(int** dna_seq_iter, int spot_t1, int tern_index, int center_type, int left_type, int right_type, int left_change, int right_change)
{
	int left_spot_type = left_type;
	int right_spot_type = right_type;

	vector<int> bind = {0};
	vector<int> unbind1 = {3};
	vector<int> unbind2 = {31,51,64,82};
	vector<int> nothing = {1,2,4};

	int spot_t1_2 = -2;
	int index;
	int spot_type;

	if (tern_index == 1)
	{
		spot_type = left_change;
		left_spot_type = update_side(tern_index,left_type);
		update_cohes(tern_index,spot_t1,left_type);
	}
	else if (tern_index == 2)
	{
		spot_type = right_change;
		right_spot_type = update_side(tern_index,right_type);
		update_cohes(tern_index,spot_t1,right_type);
	}
	else
	{
		if (Inarray(center_type,bind))
		{
			spot_type = 3;
			cohes_lst.push_back(vector<int>(2,spot_t1));
		}
		else if (Inarray(center_type,unbind1))
		{
			spot_type = 0;
			index = Index2d(spot_t1,cohes_lst,0);
			cohes_lst.erase(cohes_lst.begin()+index);
		}
		else if (Inarray(center_type,unbind2))
		{
			if (center_type == 31) {spot_type = 0;}
			else {spot_type = center_type%10;}

			index = Index2d(spot_t1,cohes_lst,0);
			spot_t1_2 = cohes_lst[index][1];
			// cout << "spot_t1_2 is: " << spot_t1_2 << endl;
			cohes_lst.erase(cohes_lst.begin()+index);
		}
		else if (Inarray(center_type,nothing))
		{
			cout << "Error here." << endl;
		}

	}

	// cout << "spot_t1 is: " << spot_t1 << endl;

	dna_seq_iter[spot_t1-1][1] = left_spot_type;
	dna_seq_iter[spot_t1+1][1] = right_spot_type;
	dna_seq_iter[spot_t1][1] = spot_type;
	
	// for (int ip=0;ip<26;ip++){
	// 	for (int iq=0;iq<2;iq++){
	// 		cout << dna_seq_iter[ip][iq] << '\t';
	// 	}
	// 	cout << endl;
	// }

	// print_vec(cohes_lst);
	// cout << endl;

	return spot_t1_2;

}


void update_cohes(int tern_index, int spot_t1, int side_type)
{
	vector<int> cohes = {3};
	vector<int> left_cohes = {32,52,74,81};
	vector<int> right_cohes = {31,51,64,82};
	int index;
	int spot_cohes;
	
	if (tern_index == 1)
	{
		spot_cohes = spot_t1-1;
		if (Inarray(side_type,cohes))
		{
			index = Index2d(spot_cohes,cohes_lst,0);
			cohes_lst[index][1] += 1;
		}
		else if (Inarray(side_type,left_cohes))
		{
			index = Index2d(spot_cohes,cohes_lst,1);
			cohes_lst[index][1] += 1;
		}
		else;
	}

	else if (tern_index == 2)
	{
		spot_cohes = spot_t1+1;
		if (Inarray(side_type,cohes))
		{
			index = Index2d(spot_cohes,cohes_lst,1);
			cohes_lst[index][0] -= 1;
		}
		else if (Inarray(side_type,right_cohes))
		{
			index = Index2d(spot_cohes,cohes_lst,0);
			// cout << endl;
			// cout << "old cohes_lst[index][0] is: " << cohes_lst[index][0] << endl;
			cohes_lst[index][0] -= 1;
		}
		else;
	}
			// cout << "tern_index is: " << tern_index << endl;
			// cout << "spot_t1 is: " << spot_t1 << endl;
			// cout << "spot_cohes is: " << spot_cohes << endl;
			// cout << "side_type is: " << side_type << endl;
			// cout << "Index is: " << index << endl;
			// cout << "new cohes_lst[index][0] is: " << cohes_lst[index][0] << endl;
			// int l = cohes_lst.size();
			// cout << "the length of the cohes_lst is: " << l << endl;
			// print_vec(cohes_lst);
			// cout << endl;
}


int update_side(int tern_index,int side_type)
{
	if (tern_index == 1)
	{
		if (side_type == 3) {side_type = 31;}
		else if (side_type == 32) {side_type = 0;}
		else if (side_type == 52 || side_type == 74 || side_type == 81) {side_type = side_type%10;}
		else;
	}
	else if (tern_index == 2)
	{
		if (side_type == 3) {side_type = 32;}
		else if (side_type == 31) {side_type = 0;}
		else if (side_type == 51 || side_type == 64 || side_type == 82) {side_type = side_type%10;}
		else;
	}
	else {cout << "The orientation is wrongly assigned." << endl;}

	return side_type;
}



// ~~~~~~~~~ //
// Gillespie //
// ~~~~~~~~~ //

void gillespie(vector<double> react_lst, int size, double *ti, int *step)
{
	double seed_1 = 0.0;
	double seed_2 = 0.0;
	seed_1 = rand()%(int(RES)+1)/(float)(RES+1);
	seed_2 = rand()%(int(RES)+1)/(float)(RES+1);
	// cout << seed_1 << endl;
	// cout << seed_2 << endl;

	double alpha = 0.0;
	for (int i=0; i<size; i++)
	{
		alpha += react_lst[i];
	}

	double tau = 0.0;
	tau = 1.0/alpha * log(1.0/seed_1);

	double cum = 0.0;
	for (int j=0; j<size; j++)
	{
		if ((seed_2 >= cum/alpha) && (seed_2 < (react_lst[j]+cum)/alpha))
		{
			*ti = tau;
			*step = j;
			break;
		}
		else
		{
			cum += react_lst[j];
		}
	}
}



// ~~~~~~~ //
// Ternary //
// ~~~~~~~ //

vectorCol ternry_all(int **seq_iter, int nbead, double **ctcf_perm)
{
	int i,j;
	double prob;
	// vector<double> tern_all;
	vector<int> temp;
	vectorCol newVectors;

	for (i=1;i<nbead+1;i++)
	{
		double tern[3] = {0.0};
		// ternry_each( *((*seq_iter+i)+1), i, *((*seq_iter+i-1)+1), *((*seq_iter+i+1)+1), nbead, tern, ctcf_perm );
		ternry_each( seq_iter[i][1] , i, seq_iter[i-1][1], seq_iter[i+1][1], nbead, tern, ctcf_perm );
		for (j=0;j<3;j++)
		{
			// cout << tern[j] << endl;
			newVectors.tern_all.push_back(tern[j]);
		}
	}

	// remote binding
	for (i=1;i<nbead+1;i++)
	{
		for (j=i+1;j<nbead+1;j++)
		{
			if ( (seq_iter[i][1] == 0) && (seq_iter[j][1] == 0) && (seq_iter[j][0]-seq_iter[i][0] <= BIND_CUTOFF) )
			{
				temp.push_back(seq_iter[i][0]);
				temp.push_back(seq_iter[j][0]);
				newVectors.zero_seq_iter.push_back(temp);
				temp.clear();
			}
		}
	}

	for (i=0;i<newVectors.zero_seq_iter.size();i++)
	{
		prob = K_ON * pow((newVectors.zero_seq_iter[i][1] - newVectors.zero_seq_iter[i][0]),-0.75);
		newVectors.tern_all.push_back(prob);
	}

	// for (int ip=0;ip<newVectors.tern_all.size();ip++) {cout << newVectors.tern_all[ip] << '\t';}
	// cout << endl;

	return newVectors;
}


void ternry_each(int center_type, int center_id, int left_type, int right_type, int nbead, double *tern, double **ctcf_perm)
{
	double p_left = 0.0;
	double p_right = 0.0;

	if (center_id == 1)
	{
		p_right = ctcf_perm[center_id][1];
	}
	else if (center_id == nbead)
	{
		p_left = ctcf_perm[center_id-2][1];
	}
	else
	{
		p_left = ctcf_perm[center_id-2][1];
		p_right = ctcf_perm[center_id][1];
	}
	each_step(center_type, left_type, right_type, p_left, p_right, tern);
}

 
void each_step(int center_type, int left_type, int right_type, double p_left, double p_right, double *tern)
{
	vector <int> ctcf = {1,2,4};
	vector <int> cohesin_1 = {3,31,51,64,82};

	// first bead
	if (left_type == -1) {
		if (p_left != 0.0) {cout << "Error with the first bead" << endl;}
		*tern = 0;
		*(tern+1) = K*p_right;
		if (center_type == 0) {
			*(tern+2) = K_ON/(pow(0.5,-0.75));
		} else if ( Inarray(center_type, cohesin_1) ) {
			*(tern+2) = K_OFF;
		} else {
			*(tern+2) = 0.0;
		}
	}

	// last bead
	else if (right_type == -1) {
		if (p_right != 0.0) {cout << "Error with the last bead" << endl;}
		*tern = K*p_left;
		*(tern+1) = 0;
		if (center_type == 0) {
			*(tern+2) = K_ON/(pow(0.5,-0.75));
		} else if ( Inarray(center_type, cohesin_1) ) {
			*(tern+2) = K_OFF;
		} else {
			*(tern+2) = 0.0;
		}
	}

	// beads in the middle
	else {
		if (center_type == 0) {
			*(tern+2) = K_ON/(pow(0.5,-0.75));
		} else if ( Inarray(center_type, cohesin_1) ) {
			*(tern+2) = K_OFF;
		} else {
			*(tern+2) = 0.0;
		}

		if (center_type == 0 || Inarray(center_type,ctcf))
		{
			if (left_type == 3 || left_type == 32 || left_type == 81) {
				*tern = K;
			} else if (left_type == 52 || left_type == 74) {
				if (p_left != 0) {*tern = K*p_left;}
				else {cout << "CTCF position is wrong." << endl;}
			} else {
				*tern = 0.0;
			}

			if (right_type == 3 || right_type == 31 || right_type == 82) {
				*(tern+1) = K;
			} else if (right_type == 51 || right_type == 64) {
				if (p_right != 0) {*(tern+1) = K*p_right;}
				else {cout << "CTCF position is wrong." << endl;}
			} else {
				*(tern+1) = 0.0;
			}
		}

	}
}



//~~~~~~~~~~~~~~~~~~//
// Assist Functions //
//~~~~~~~~~~~~~~~~~~//

int Inarray(int val, vector<int> vec)
{
	bool found = ( find(vec.begin(), vec.end(), val) != vec.end() );
	if (found) {return 1;}
	else {return 0;}
}


int Index2d(int spot_cohes, vector<vector<int> > cohes_lst, int col)
{
	int index;
	for (int ir=0;ir<cohes_lst.size();ir++)
	{
		if (cohes_lst[ir][col] == spot_cohes)
		{
            index = ir;
            break;
		}
	}

	return index;
}


int Index1d(int spot_cohes, vector<int> ctcf_lst)
{
	int index;
	for (int ir=0;ir<ctcf_lst.size();ir++)
	{
		if (ctcf_lst[ir] == spot_cohes)
		{
            index = ir;
            break;
		}
	}

	return index;
}


double max_diff(double **lst_2d_1,double **lst_2d_2,int di)
{
	int ip,iq;
	double diff = 0.0;
	double diff_max = 0.0;

	for (ip=0;ip<di;ip++)
	{
		for (iq=0;iq<di;iq++)
		{
			diff = abs(lst_2d_1[ip][iq] - lst_2d_2[ip][iq]);
			if (diff > diff_max){
				diff_max = diff;
			}
		}
	}

	return diff_max;
}


double **readdblfile(const char* fi_name, int row, int col)
{
	ifstream fin(fi_name);
	double **mat;
	mat = (double**)malloc(row*sizeof(double*));
	for (int ii=0;ii<row;ii++) {mat[ii] = (double*)malloc(col*sizeof(double));}

	while (!fin.eof())
	{
		for (int i=0;i<row;i++)
		{
			for (int j=0;j<col;j++)
			{
				fin >> mat[i][j];
			}
		}
	}

	fin.close();
	return mat;
}


int **readintfile(const char* fi_name, int row, int col)
{
	if (row > 0)
	{
		ifstream fin(fi_name);
		int **mat;
		mat = (int**)malloc(row*sizeof(int*));
		for (int ii=0;ii<row;ii++) {mat[ii] = (int*)malloc(col*sizeof(int));}

		while (!fin.eof())
		{
			for (int i=0;i<row;i++)
			{
				for (int j=0;j<col;j++)
				{
					fin >> mat[i][j];
				}
			}
		}
		
		fin.close();
		return mat;
	}
	else {return NULL;}
}


int rowoffile(const char* fi_name)
{
	int num_of_lines = 0;
    string line;
    ifstream fin(fi_name);

    while (getline(fin, line))
        num_of_lines++;

    return num_of_lines;
}



//~~~~~~//
// Temp //
//~~~~~~//

void print_vec(vector<vector<int> > cohes_lst)
{
    for (int ip = 0; ip <cohes_lst.size(); ip++)
    {
        for (int iq = 0; iq <cohes_lst[0].size(); iq++)
        {
            cout << cohes_lst[ip][iq] << '\t';
        }
        cout << endl;
    }
}


void print_mat(double **mat, int size)
{
	for (int ip=0; ip <size; ip++)
	{
		for (int iq=0; iq <size; iq++)
		{
			cout << mat[ip][iq] << '\t';
		}
		cout << endl;
	}
}



//~~~~~~//
// Main //
//~~~~~~//

int main(int argc, char **argv)
{
	// const char* dna_seq_fi = "ctcf_position_chr21_Gm12878_From20MbTo45Mb.txt";
	// read_restart
	const char* dna_seq_fi = "ctcf_position.txt";
	const char* ctcf_seq_fi = "raw_ctcf_position_artificial.txt";
	const char* ctcf_perm_fi = "scale_param.txt";
	const char* cohes_lst_fi = "cohes_lst.txt";
	const char* prob_mat_raw_fi = "prob_mat_raw.txt";
	const char* iter_num_fi = "iter_num.txt";

	// int **dna_seq = readintfile(dna_seq_fi,BEAD,COL);
	// read_restart
	int **dna_seq = readintfile(dna_seq_fi,BEAD,COL);
	int **ctcf_seq = readintfile(ctcf_seq_fi,CTCF,COL);
	double **ctcf_perm = readdblfile(ctcf_perm_fi,BEAD,COL);

	int cohes_lst_lines = rowoffile(cohes_lst_fi);
	int **cohes_lst_old = readintfile(cohes_lst_fi,cohes_lst_lines,COL);

	double **prob_mat_raw = readdblfile(prob_mat_raw_fi,CTCF,CTCF);
	int **iter_num = readintfile(iter_num_fi,1,1);

	// read_restart
	// for (int id=0;id<BEAD;id++)
	// {
	// 	if (dna_seq[id][1] == 3)
	// 	{
	// 		dna_seq[id][1] = 0;
	// 	}
	// }

	vector<int> ctcf_lst;
	// ctcf_lst = (int*)malloc(CTCF*sizeof(int));

	for (int ic=0;ic<CTCF;ic++)
	{
		ctcf_lst.push_back(ctcf_seq[ic][0]);
	}

	int **dna_seq_iter;
	dna_seq_iter = (int**)malloc((BEAD+2)*sizeof(int*));

	for (int ii=0;ii<BEAD+2;ii++)
	{
	 	dna_seq_iter[ii] = (int*)malloc(2*sizeof(int));
	}

	// dna_seq_iter = new int *[BEAD+2];
	// for( int i=0; i<BEAD+2; i++ )
	// {
	// 	dna_seq_iter[i] = new int[2];
 // 	}


	// initializa
	dna_seq_iter[0][0] = 0;
	dna_seq_iter[0][1] = -1;
	dna_seq_iter[BEAD+1][0] = BEAD+1;
	dna_seq_iter[BEAD+1][1] = -1;

	for (int it=1;it<BEAD+1;it++)
	{
		for (int jt=0;jt<2;jt++)
		{
			dna_seq_iter[it][jt] = dna_seq[it-1][jt];
		}
	}


	// static vector<vector<int> > cohes_lst;
	if (cohes_lst_lines > 0)
	{
		for (int ico=0;ico<cohes_lst_lines;ico++)
		{
			vector<int> temp;
			temp.push_back(cohes_lst_old[ico][0]);
			temp.push_back(cohes_lst_old[ico][1]);
			cohes_lst.push_back(temp);
			temp.clear();
		}
	}
	else;


	// debugging print
	// for (int ip=0;ip<BEAD+2;ip++)
	// {
	// 	for (int ip2=0;ip2<2;ip2++)
	// 	{
	// 		cout << dna_seq_iter[ip][ip2] << '\t';
	// 	}
	// 	cout << endl;
	// }
	
	srand((unsigned) time(NULL));
	main_proc(dna_seq_iter,ctcf_lst,CTCF,ctcf_perm,prob_mat_raw,iter_num);

	// debugging print
	// for (int it=0;it<tern_all.size();it++)
	// {
	// 	cout << tern_all[it] << '\t' << it << endl;

	// }
	// 

	return 0;
}
