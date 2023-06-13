#include <iostream>
#include <string>
#include <vector>

int pow (int a){
    if (a % 2 == 0)
        return 1;
    else
        return -1;
}

template <class type, class dimension>
class matrix{
public:
    dimension n, m;
    std::vector<std::vector<type>> out;
    std::vector<type> in;

    matrix(const dimension& n, const dimension&m) : n(n), m(m) {};
    matrix(const std::string& diag, const type& a, const dimension& n){
        if (diag == "diag"){
            this->n = n;     this->m = n;
            for (dimension i = 0; i < n; ++i) {
                for (dimension j = 0; j < m; ++j)
                    if (i == j)
                        in.push_back(a);
                    else
                        in.push_back(0);
                out.push_back(in);
                in.clear();
            }
        } else   {std::cerr << "\n               error              \n"; exit(EXIT_FAILURE);}
    }
    matrix(const char& type_of_matrix, const dimension& n){
        if (type_of_matrix == 'E')     *this = matrix ("diag", 1, n);
        else
        if (type_of_matrix == '0') {
            this->n = n;
            this->m = n;
            for (dimension i = 0; i < n; ++i) {
                for (dimension j = 0; j < n; ++j)
                    in.push_back(0);
                out.push_back(in);
                in.clear();
            }
        } else {std::cerr << "\n                       error                         \n"; exit(EXIT_FAILURE);}
    }
    matrix(std::vector<type> mas, dimension n, dimension m){
        this->n = n;  this->m = m;
        if (m * n != mas.size()){
            std::cerr << "\n                    error                             \n"; exit(EXIT_FAILURE);}
        for (dimension i = 0; i < n; ++i){
            for (dimension j = 0; j < m; ++j)
                in.push_back(mas[i * m + j]);
            out.push_back(in);
            in.clear();
        }
    }


    matrix operator* (matrix B){
        matrix<type, dimension> answer (this->n, B.m);
        if (this->m != B.n){
            std::cerr << "\n                   error            \n"; exit(EXIT_FAILURE);}
        type sum = 0;
        for (dimension i = 0; i < n; ++i){
            for (dimension j = 0; j < m; ++j){
                for (dimension k = 0; k < m; ++k)
                    sum += ( this->out[i][k] * B.out[k][j] );
                answer.in.push_back(sum);
                sum = 0;
            }
            answer.out.push_back(answer.in);
            answer.in.clear();
        }
        return answer;
    }

    matrix operator* (type a){
        matrix<type, dimension> answer (n, m);
        for (dimension i = 0; i < n; ++i){
            for (dimension j = 0; j < m; ++j)
                if (this->out[i][j] == 0) answer.in.push_back(0);
                else answer.in.push_back(a * (this->out[i][j]));
            answer.out.push_back(answer.in);
            answer.in.clear();
        }
        return answer;
    }

    matrix operator+ (matrix B){
        matrix<type, dimension> answer (this->n, this->m);
        if (this->m != B.m and this->n != B.n){
            std::cerr << "\n                   error            \n"; exit(EXIT_FAILURE);}
        for (dimension i = 0; i < n; ++i){
            for (dimension j = 0; j < m; ++j)
                answer.in.push_back(this->out[i][j] + B.out[i][j]);
            answer.out.push_back(answer.in);
            answer.in.clear();
        }
        return answer;
    }

    matrix operator- (matrix B){
        matrix<type, dimension> answer (this->n, this->m);
        if (this->m != B.m and this->n != B.n)
            std::cerr << "\n                   error            \n";
        for (dimension i = 0; i < n; ++i){
            for (dimension j = 0; j < m; ++j)
                answer.in.push_back(this->out[i][j] - B.out[i][j]);
            answer.out.push_back(answer.in);
            answer.in.clear();
        }
        return answer;
    }

    bool operator== (matrix B){
        if (this->m != B.m and this->n != B.n){
            std::cerr << "\n                   error            \n"; exit(EXIT_FAILURE);}
        for (dimension i = 0; i < n; ++i){
            for (dimension j = 0; j < m; ++j)
                if (this->out[i][j] != B.out[i][j])
                    return false;
        }
        return true;
    }

    bool operator!= (matrix B){
        dimension flag = 0;
        if (this->m != B.m and this->n != B.n) {
            std::cerr << "\n                   error            \n"; exit(EXIT_FAILURE);}
        for (dimension i = 0; i < n; ++i){
            for (dimension j = 0; j < m; ++j)
                if (this->out[i][j] == B.out[i][j])
                    ++flag;
        }
        if (flag == this->m * this->n)
            return false;
        else
            return true;
    }

    bool operator== (type a){
        matrix<type, dimension> A ("diag", a, this->n);
        if (A == *this)
            return true;
        else
            return false;
    }

    matrix operator=(matrix<type, dimension> A){
        this->m = A.m;
        this->n = A.n;
        this->out.clear();
        for (int i = 0; i < A.n; ++i) {
            for (int j = 0; j < A.m; ++j)
                this->in.push_back(A.out[i][j]);
            this->out.push_back(this->in);
            this->in.clear();
        }
        return *this;
    }

    matrix minor (const dimension& a, const dimension& b) {
        matrix<type, dimension> answer(n - 1, m - 1);
        for (dimension i = 0; i < n; ++i) {
            if (i == a)
                continue;
            else {
                for (dimension j = 0; j < m; ++j)
                    if (j != b)
                        answer.in.push_back(out[i][j]);
                answer.out.push_back(answer.in);
                answer.in.clear();
            }
        }
        return answer;
    }

    type det (){
        type det_sum = 0;
        for (dimension i = 0; i < m; ++i){
            if (this->n == 1)
                return (this->out[0][0]);
            else
                det_sum = det_sum + pow(i) * this->out[0][i] * (this->minor(0, i)).det();
        }
        return det_sum;
    }

    matrix transp(){
        matrix<type, dimension> answer (n, m);
        for (dimension i = 0; i < n; ++i) {
            for (dimension j = 0; j < m; ++j)
                answer.in.push_back(this->out[j][i]);
            answer.out.push_back(answer.in);
            answer.in.clear();
        }
        return answer;
    }

    matrix <double, dimension> operator!(){
        matrix<double, dimension> answer(n, m);
        for (dimension i = 0; i < n; ++i) {
            for (dimension j = 0; j < m; ++j)
                answer.in.push_back(pow(i + j) * this->minor(i, j).det());
            answer.out.push_back(answer.in);
            answer.in.clear();
        }
        return (answer.transp()) * (1 / this->det());
    }
};

template <typename type, typename dimension>
matrix<type, dimension> operator* (const type& a, matrix<type, dimension> A){
    return A * a;
}

template<typename type, typename dimension>
std::ostream& operator<<(std::ostream& stream, const matrix<type, dimension>& A)
{
    for (dimension i = 0; i < A.n; ++i) {
        for (dimension j = 0; j < A.m; ++j) {
            stream << A.out[i][j] << " ";
        }
        stream << "\n";
    }
    return stream;
}

int main() {
    matrix<double, int> A("diag", 2, 5);
    std::cout << A << "\n\n";

    matrix<double, int> E('E', 4);
    std::cout << E << "\n\n";

    matrix<double, int> O('0', 5);
    std::cout << O << "\n\n";

    std::vector<double> mas = {1, 2, 3, 4, 5, 6, 7, 8, 9};

    matrix<double, int> B (mas, 3, 3);
    std::cout << B << "\n\n";

    std::cout << B * B << "\n\n";

    std::cout << B.minor(1, 1) << "\n\n";

    std::cout << B.det() << "\n\n";

    std::cout << 3.0 * B << "\n\n";

    std::cout << (E + E) << "\n\n";

    std::cout << !(A) << "\n\n";

    std::cout << (E == 1);
}
