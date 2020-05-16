#include<iostream>
#include<vector>
#include<list>
#include <fstream> // for file-access
#include <string>
#include<sstream>
using namespace std;

int sign(double x) {
	if (x > 0) {
		return 1;
	}

	if (x < 0) {
		return -1;
	}

	return 0;
}

void printer(list<vector<double>> &A) {
	for (auto it = A.begin(); it != A.end(); it++) {
		for (size_t i = 0; i < it->size(); ++i) {
			cout<< (*it)[i];
			cout<< " ";
		}
		cout << "\n";
	}
}

void scalar_multiplication(double c, vector<double> &a) {
	for (size_t i = 0; i < a.size(); ++i) {
		a[i] *= c;
	}
}

/* Save a + b to c */
void vector_addition(const vector<double> &a,
		     const vector<double> &b,
		           vector<double> &c) {
	for (size_t i = 0; i < c.size(); ++i) {
		c[i] += a[i] + b[i];
	}
}

void readfile(size_t num_rows, size_t num_cols, vector<double> &c, list< vector<double> > &A, char* filename){
    ifstream infile(filename); //open the file
    string line;
	list<vector <double> >::iterator it=A.begin();
	if (infile.is_open() && infile.good()) {
        int counter = 0;
        string line;
        while (getline(infile, line)) {
            istringstream row(line);
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
                    vector<double> column;
                        column.push_back(field);
                        A.push_back(column);
                }
            }
            else{
                advance(it,1);
                while(row>>field){
                     it->push_back(field);		
                }
	    }		
            counter++;
        }
    } else {
        cout << "Failed to open file";
    }
}

void fourier_motzkin(list<vector<double>> &A) {
	list<vector<double>> U;
	list<vector<double>> L;

	/* Classify and normalize the inequalities in A */

	/* while there is more than one variable ... */
	while (A.front().size() != 1) {
		auto it = A.begin();
		while (it != A.end()) {
			auto c = it->back();
			it->pop_back();
			if (sign(c) != 0) {
				scalar_multiplication(1.0/abs(c), *it);
				if (sign(c) > 0) {
					U.splice(U.end(), A, it++);
				}
				if (sign(c) < 0) {
					L.splice(L.end(), A, it++);
				}
			} else {
				it++;
			}
		}

		cout << "\nDecomposition and normalization:\n";
		cout << "N\n";
		printer(A);
		cout << "U\n";
		printer(U);
		cout << "L\n";
		printer(L);

		/** Construct and add new inequalities to A **/
		if (!U.empty() && !L.empty()) {
			for (auto it_u = U.begin(); it_u != U.end(); it_u++) {
				for (auto it_l = L.begin(); it_l != L.end(); it_l++) {
					A.emplace_back(it_u->size());
					vector_addition(*it_u, *it_l, A.back());
				}
			}
			L.erase(L.begin(), L.end());
		        U.erase(U.begin(), U.end());
		} else {
			A.splice(A.end(), L);
			A.splice(A.end(), U);
		}

		cout << "\nReduced system:\n";
		cout << "A\n";
		printer(A);
		cout << "U\n";
		printer(U);
		cout << "L\n";
		printer(L);
	}

	/* Since the first column is the b vector, we must check if it is >= 0 */
	bool admissible = true;
	for (auto it = A.begin(); it != A.end(); it++) {
		admissible = admissible && (it->at(0) >= 0);
	}

	if (admissible) {
		cout << "The system is admissible.\n";
	} else {
		cout << "The system is not admissible.\n";
	}
}

int main() {
	list<vector<double>> A = {{0, 1}, {0, -1}};
	printer(A);
	
	fourier_motzkin(A);
}
