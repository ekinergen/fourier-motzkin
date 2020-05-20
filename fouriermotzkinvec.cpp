#include <fstream> // for file-access
#include <string>
#include<sstream>
#include <iostream>
#include<vector>
#include<list>
const double small=-10000;
const double big = 10000; //these can be changed, our test instances are small enough anyway. I deliberately didn't choose them too big in order to prevent overflows.

void readfile(size_t num_rows, size_t num_cols, std::list<double> &c, std::vector< std::list<double> > &A, char* filename){
    std::ifstream infile(filename); //open the file
    std::string line;
    if (infile.is_open() && infile.good()) {
        int counter = 0;
        std::string line;
        while (std::getline(infile, line)) {
            std::istringstream row(line);
            double field;
            if(counter == 0){
                row>> num_rows >> num_cols;
            }
            if(counter ==1){
                while(row>>field){
                    c.push_back(field);
                }
            }
            if(counter==2){
                while(row>>field){
                    std::list<double> column;
					column.push_back(field);
					A.push_back(column);
                }
            }
            else{
                while(row>>field){
					A[counter-3].push_back(field);
                }
            }
            counter++;
        }
    } else {
        std::cout << "Failed to open file..";
    }
}

void normalize(std::list<double> &vec){ //normalizes at the last variable
	if (vec.back()==0) return;
	double divide=vec.back();
	for(auto it=vec.begin(); it!=vec.end(); it++){
		*(it)=*(it)/divide;
		
	}
}

void printmatrix(std::vector<std::list<double> > &A){
	for(size_t i=0; i<A.size(); i++){
		for(auto j: A[i]){
			std::cout<<j<<" ";
		}
		std::cout<<std::endl;
	}
}

double max(std::list<double> &vec){
	double val = *(vec.begin());
	for(auto i: vec){
		if(val<i) val=i;
	}
	return val;
}

double min_b(std::vector<std::list<double> > &A){
	double val = *(A[0].begin());
	for(size_t i=1; i<A.size(); i++){
		if(val>*(A[i].begin())) val=*(A[i].begin());
	}
	return val;
}

void add_rows(std::list<double> &row_u, std::list<double> &row_l, std::list<double> &sum){
	auto it_u= row_u.begin(); 
	auto it_l= row_l.begin();
	*(sum.begin())=*(row_u.begin())-*(row_l.begin());
	it_u++;
	it_l++;
	for(auto it=std::next(sum.begin()); it!=sum.end(); it++){
		*(it)=-*(it_u)+*(it_l); //we had x_n<=*(it_u) and *(it_l)<=x_n, so now we want *(it_l)-*(it_u)<=0
		it_u++;
		it_l++;
	}
}

void set_last_variable(std::vector<std::list<double> > &A, double val){ //We are given a system of inequalities a1x1+a2x2+...+a_nx_n<=b and we get rid of the last variable of each inequality by plugging in val into the last variable, i.e. we want a1x1+a2x2+...+a_(n-1)x_(n-1)<=b-a_n*val.
	if (A[0].size()==1) return; //this will be handled in backtracking 
	for(size_t i = 0; i<A.size(); i++){
		*(A[i].begin())-=A[i].back()*val; //A[i][0]=b in the i-th inequality
		A[i].pop_back();
	}
	
}

void set_first_variable(std::vector<std::list<double> > &A, double val){ //We are given a system of inequalities a1x1+a2x2+...+a_nx_n<=b and we get rid of the first variable of each inequality by plugging in val into the first variable, i.e. we want a2x2+...+a_nx_n<=b-a_1*val.
	for(size_t i = 0; i<A.size(); i++){
		auto it = A[i].begin();
		it++; //We go to second element of linked list because the first one is from b...
		*(A[i].begin())-=(*it)*val; //A[i].begin()=b in the i-th inequality
		A[i].erase(it); //This is the exact reason we used linked lists for each inequality
	}
}

void put_variables_to_other_side(std::list<double> &vec){
		
		for(auto it=std::next(vec.begin()); it!=vec.end(); it++){
		*(it)=*(it)*-1;
	}
}

void fouriermotzkin(std::vector<std::list<double> > &mat, size_t num_rows, size_t num_cols){
	std::vector<std::vector<std::list<double> > > A;
	std::vector<double> admissible_solution(mat[0].size()-1, 0); //num_cols=number of variables+1 (because we also store b in a column)
	A.push_back(mat);
	
	while(A.back()[0].size()!=1){	//A.back() is the last matrix that we added. If its first row (and hence all rows) has size 1, it means that there are 0 variables.
		std::vector<std::list<double> > U;
		std::vector<std::list<double> > L;
		std::vector<std::list<double> > N;
			
		//classifying the inequalities
		for(size_t k=0; k<A.back().size(); k++){ //We always look at the last system to classify its inequalities (which we will combine and store onto the next element, aka matrix of A) 
			if(A.back()[k].back()==0){
				N.push_back(A.back()[k]); //these don't contain x_n
				N.back().pop_back();
			} 
			if(A.back()[k].back()>0){
				U.push_back(A.back()[k]); //these are of form x_n<=sth.
				normalize(U.back());
				put_variables_to_other_side(U.back()); //so that x_1,... x_{n-1} go to the right side, we switch their signs. Now we have 0<=l1x1+... etc.

				(U.back()).pop_back();
			}
			if(A.back()[k].back()<0) {
				L.push_back(A.back()[k]); //these are of form x_n>=sth.
				normalize(L.back());
				put_variables_to_other_side(L.back()); //so that x_1,... x_{n-1} go to the right side, we switch their signs. Now we have 0<=l1x1+... etc.
				(L.back()).pop_back();
			}
		}
		//we get the new inequalities
		std::vector<std::list<double> > B;
		if((not U.empty()) and not L.empty()){
			 
			for(size_t k=0; k<U.size(); k++){
				for(size_t j=0; j<L.size(); j++){
					std::list<double> sum(U[k].size(), 0);
					add_rows(U[k], L[j], sum);
					B.push_back(sum);
				}
			}
		}
		
		for(size_t k=0; k<N.size(); k++){
			B.push_back(N[k]);
		}
		if(B.empty()){//in this case we know that the last variable either has upper bounds only or lower bounds only. So we can actually already set a large or small value to it and continue. 
			if(U.empty()){ //then we only have lower bounds for x_n, also in terms of other variables etc. 
				admissible_solution[L[0].size()]=big; //We are at the L[0].size()-1-th variable at this time, again because we store b as a column. However, we just popped back the last variable from L, so size() is one smaller than it should have been normally.
				
				set_last_variable(A.back(),big);
			}
			if(L.empty()){//in this case, x_n only has upper bounds. So we can set a very small number and continue.
				admissible_solution[U[0].size()]=small;
				set_last_variable(A.back(),small) ;
			}
			
			continue; //B is empty, so we don't want to put it as a phase, instead, we work with the same system as in this iteration, except we evaluate at x_n.
		} 
		A.push_back(B);
	}
	
	/*deciding if the original matrix is admissible, i.e. the last system is admissible. The last system consists of 0<=b, where the 
	  b's are components of 1-column matrix in the last system. Hence it suffices to check all b's are larger than 0;*/
	bool admissible = true;
	for(size_t i =0 ; i<A.back().size(); i++){ //iterating over all columns
		if(*(A.back()[i].begin())<0) admissible=false;
	}
	/*backtracking*/
	if(admissible){
		size_t var_index=0; //Throughout commentary, let us call this number n. 
		for(size_t phase=0; phase<A.size()-1; phase++){ //A.size()=number of phases. We set some values without adding a new phase to them, so if we have such a situation, we need to skip that. In that case, we increment var_index and go on in the for-loop. That is, we set the variables that we found not using the variable i (because that only counts the actual phases), but using the variable var_index.
			size_t i = A.size()-2-phase;
			while(admissible_solution[var_index]!=0){ //this means that we already set a variable here without a new phase
				var_index++;
			}
			//Starting from the second-to-last iteration, where we had inequalities of form x_1<=some_upper_bound, x_1>=some_lower_bound, 0<=sth. By the thm from the lecture, we know that the smallest upper bound is still at least the biggest of the lower bounds. So we can directly set either of those as a value.
			double val=A[i][0].front()/A[i][0].back();
			for(size_t j=0; j<A[i].size();j++){ //we only look at inequalities of kind x_n<=sth. A[i][j].back() is the coefficient of the last variable.
				if((val>A[i][j].front()/A[i][j].back() and A[i][j].back()>0) or (val<A[i][j].front()/A[i][j].back() and A[i][j].back()<0)){//in the first case, the 0-term must have been an upper bound and in the second case a lower bound
					val = A[i][j].front()/A[i][j].back();
					
				} 
			}
			admissible_solution[var_index]=val;
			if(i!=0) set_first_variable(A[i-1], val); //In one phase prior, we plug in the value of x_n.
			//printmatrix(A[i-1]);
			var_index++;
		}
		for(size_t i = 0; i<admissible_solution.size(); i++){
			std::cout<<admissible_solution[i]<<" ";
		}
		std::cout<<std::endl;
	}
	else{
		std::cout<<"empty ";
		//keine Ahnung
	}
}

int main(int argc, char* argv[])
{ 
    size_t num_rows=0, num_cols=0;
    std::list<double> c;
    std::vector< std::list<double> > A; //b ist die erste Spalte von A
    readfile(num_rows, num_cols, c, A, argv[1]);
	fouriermotzkin(A, num_rows, num_cols);

    return 0;
}