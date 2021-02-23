#include <cmath>
#include <iostream>
#include <initializer_list>
#include <memory>
#include <map>
#include <stdexcept>
#include <utility>

using namespace std;

template <typename T>
class Vector
{
    // Your implementation of the Vector class starts here
    private:
        T* data;        // dynamical array
        int length; 

    public:
        Vector()
        :length(0), data(nullptr)
        {
        }

        Vector(int n)
        :length(n), data(new T[n])
        {
        }

        int len() const
        {
            return length;
        }

        Vector(const Vector<T>& other) // Copy constructor
        :Vector(other.length)
        {
            for(int i = 0; i<length; i++)
            {
                data[i] = other.data[i];
            }    
        }
        
        Vector(Vector<T>&& other) // Move constructor
        :length(other.len()), data(other.data)
        {
            other.length = 0;
            other.data = nullptr;  
        }


        Vector(const initializer_list<T>& list) // Initializer list constructor
        :Vector((int)list.size())
        {
            uninitialized_copy(list.begin(), list.end(), data);
        }

        ~Vector() {delete[] data;length=0;data=nullptr;} // Destructor

        Vector<T>& operator=(const Vector<T>& other) // Copy operator
        {   
            if (this != &other)
            {
                length = other.length;
                delete[] data;
                data = new T[length];           
                for(int i=0;i<length;i++)
                {
                    data[i] = other.data[i];
                }
            }
            return *this;
        }

        Vector<T>& operator=(Vector<T>&& other) // Move operator
        {   
            if (this != &other)
            {
                delete[] data;
                data = other.data;
                other.data = nullptr;
                length = other.length;
                other.length = 0;
            }
            return *this;
        }

        T& operator[](int i) 
        {   
            if(i>=length)
            { 
                cout<<"Index of vector out of bounds. Exiting..."<<endl;
                abort();
                
            };  
            return data[i]; 
        }

        const T& operator[](int i) const
        {   
            if(i>=length)
            { 
                cout<<"Index of vector out of bounds. Exiting..."<<endl;
                abort();
                
            };
            return data[i];
        }
        
        template<typename S>
        Vector<typename common_type<S,T>::type> operator+(const Vector<S>& other) const
        {
            if (length != other.len())
            {
                throw "Vectors of different length";
            }
            else
            {
                Vector<typename common_type<S,T>::type> result_vector(length);
                for(int i=0;i<length;i++)
                {
                    result_vector.data[i] = data[i] + other.data[i];
                }
                return result_vector;
            }
        }

        template<typename S>
        Vector<typename common_type<S,T>::type> operator-(const Vector<S>& other) const
        {
            if (length != other.len())
            {
                throw "Vectors of different length";
            }
            else
            {
                Vector<typename common_type<S,T>::type> result_vector(length);
                for(int i=0;i<length;i++)
                {
                    result_vector.data[i] = data[i] - other.data[i];
                }
                return result_vector;
            }
        }

        template<typename S>
        Vector<typename common_type<S,T>::type> operator*(const S scalar) const
        {
            Vector<typename common_type<S,T>::type> result_vector(length);
            for(int i=0;i<length;i++)
            {
                result_vector[i] = scalar*data[i];
            }
            return result_vector;     
        }


};


template<typename S, typename T>
Vector<typename common_type<S,T>::type> operator*(const S scalar ,const Vector<T>& other)
{
    Vector<typename common_type<S,T>::type> result_vector(other.len());
    for(int i=0;i<other.len();i++)
    {
        result_vector[i] = scalar*other[i];
    }
    return result_vector;
}

template<typename T, typename U>
typename common_type<T,U>::type 
dot(const Vector<T>& lhs, const Vector<U>& rhs)
{
    if (rhs.len() != lhs.len())
        {
            throw "Vectors of different length";
        }
        else
        {
            typename common_type<T,U>::type result = 0;
            for(int i=0;i<lhs.len();i++)
            {
                result = result + rhs[i] * lhs[i];
            }
            return result;
        }
}

template <typename T>
class Matrix
{
    private:
        int rows;
        int columns;
        map <pair<int, int>, T> matrixMap;
    
    public:
        
        Matrix()
        :rows(0), columns(0), matrixMap()
        {
        }
        Matrix(int n, int m)
        :rows(n), columns(m)
        {
        }
       ~Matrix()
        {
            rows=0;
            columns=0;
            matrixMap.clear();
        }

        T& operator[](const std::pair<int, int>& ij)
        {
            int i = ij.first;
            int j = ij.second;

            if (i >= rows || j >= columns)
            {
                throw "Index out of range";
            }

            else
            {
                for(auto it = matrixMap.begin(); it != matrixMap.end(); ++it)
                {
                    if (i == it->first.first && j == it->first.second)
                    {
                        return it->second;
                    }
                }
                matrixMap.insert( {ij, 0});
                return matrixMap.at(ij);
            }
        }

        const T& operator()(const std::pair<int, int>& ij) const
        {
            int i = ij.first;
            int j = ij.second;

            if (i >= rows || j >= columns)
            {
                throw "Index out of range";
            }
            else
            {
                
                for(auto it = matrixMap.begin(); it != matrixMap.end(); ++it)
                {
                    if (i == it->first.first && j == it->first.second)
                    {
                        return it->second;
                    }
                }
                throw "Value does not exist";
            }
        }

        int col() const
        {
            return columns;
        }

        int row() const
        {
            return rows;
        }

        map <pair<int, int>, T> Map() const   //reference?????
        {
            return matrixMap;
        }

        void print_matrix()
        {
            for(int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                    {
                        cout << matrixMap[{i,j}]<<" ";
                    }
                cout<<""<<endl; 
            }
        }
};

template<typename T, typename U>
Vector<typename std::common_type<T,U>::type> operator*(const Matrix<T>& lhs, const Vector<U>& rhs)
{
    if(lhs.col() != rhs.len())
    {
        throw "Incompatible sizes";
    }
    else
    {
        //map <pair<int, int>, T>  matmap = lhs.Map();
        Vector<typename std::common_type<T,U>::type> result_vector(lhs.row());
        for(int i = 0; i < lhs.row(); i++)      
        {
            result_vector[i] = 0;
        }

        for (auto const& it : lhs.Map()) // This is a c++14 thing that works
        {
            int i = it.first.first;
            int j = it.first.second;
            result_vector[i] += rhs[j]*it.second;
        }
        /*for(auto it = matmap.begin(); it != matmap.end(); ++it)
        {
            int i = it->first.first;
            int j = it->first.second;
            result_vector[i] += rhs[j]*matmap.at({i,j});
        }*/
    return result_vector;
    }
}

template<typename T>
int cg(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x, T tol = (T)1e-8, int maxiter = 100)
{

    // CHECK IF SYMMETRIC?????

    Vector<T> p = b - A*x;
    Vector<T> r = b - A*x;


    T alpha, beta, prev_dot;
    int k=0;
    while (k< maxiter)
    {
        prev_dot = dot(r,r);

        alpha = prev_dot / dot(A*p, p);
        x = x + (alpha*p);
        r = r - alpha*(A*p);
        
        if (dot(r,r) <= tol*tol)
        {
            return k;
        }

        beta = dot(r,r) / prev_dot;
        p = r + beta*p;
        ++k;
    }
    
    return -1;
}

template <int n, typename T>
class Heat
{
    // Your implementation of the heat class starts here
    private:
        T alpha;
        int m;
        T dt;
        Matrix<T> M;
    public:
        
        Heat(T a, int mu ,T timestep) : alpha(a),m(mu),dt(timestep)
        {
            int points = (int) pow(mu,n);
            Matrix<T> Mt( points, points );
            T coeff = alpha*(m+1)*(m+1)*dt;
            
            for(int j=points-1 ; j>=0 ; j--)
            {
                Mt[{j,j}] = 1 + coeff*(2*n);

                for(int k=0 ; k<n ; k++)
                {
                    bool boundary_exception = !(j%(int)pow(m,k+1) < (int)pow(m,k));
                    if(boundary_exception)
                    {
                        for(int i=0 ; i<j ; i++)
                        {
                            bool is_neighbour = (j-i)==(int)pow(m,k);
                            if(is_neighbour)
                            {
                                Mt[{i,j}] = -coeff;
                                Mt[{j,i}] = -coeff;
                            }
                        }
                    }   
                }
            }
            M = Mt;
        }

        void print_heat_matrix()
        {
            M.print_matrix();
        }
        
        Matrix<T> get_matrix() const
        {
            return M;
        }


        Vector<T> exact(T t) const
        {   
            int points = (int)pow(m,n);
            Vector<T> sol(points);
            T dx = (T)1/(m+1); //
            for(int i = 0; i < points; i++)
            {
                sol[i] = exp(-t*n*alpha* pow(M_PI,2) );
                for (int k = 0; k < n; k++)
                {
                    int xk_index = i%(int)pow(m,k+1)/((int)pow(m,k));
                    T xk = xk_index * dx + dx;
                    sol[i] *= sin( M_PI * xk );
                }
            }
            return sol;
        }

        Vector<T> solve(T t) const
        {
            int points = (int)pow(m,n);   //  Mwl+1=wl
            Vector<T> sol(points);
            T dx = (T)1/(m+1);
            
            for(int i = 0; i < points; i++)
            {
                sol[i] = 1;
                for (int k = 0; k < n; k++)
                {
                    int xk_index = i%(int)pow(m,k+1)/((int)pow(m,k));
                    T xk = xk_index * dx + dx;
                    sol[i] *= sin( M_PI * xk );
                }
            }
            int nnc = 0;
            for (T l = 0; l < t ; l+=dt )
            {
                int iter = cg( M , sol, sol);
                if(iter==-1)
                {
                    ++nnc;
                }
            }
            cout<<" Non convergences: "<<nnc<<endl;
            return sol;
        }
        /*
        T getMatrix()
        {
            int const points = (int) pow(m,n);
            T Matr[points][points];
            for(int i = 0; i < points; i++)
            {
                for (int j = 0; j < points; j++)
                    {
                        Matr[i][j] = M[{i,j}];
                    }
            }
            return  Matr;
        }
        */

};


int main(int argc, char* argv[])
{
    const Vector<int> a();
    
    Vector<int> b({1,2});
    Vector<int> c(b);
    Vector<int> f = b + c;
    Vector<int> x({1,1});
    f=f*8;

    Matrix<int> A(2,2);

    A[{0,0}] =2;
    A[{1,1}] =3;
    cout<< (A*b)[1] << endl;

    int S = cg<int>(A, b, x);
    
    
    //cout<<A[{0,1}];
    
    //Vector<int> d = A*b;
    //cout<<d[0];

    Heat<2,double> Aa(0.3125, 99, 0.001);
    
    Vector<double> ex = Aa.exact(1.);
    Vector<double> num = Aa.solve(1.);      

    for(int i = 0; i < ex.len(); i++)
    {
        cout<<ex[i] - num[i]<<" "<<endl;
    }
    
    /*/
    Matrix<double> mat = Aa.get_matrix();
    for (int i = 0; i < mat.row(); i++) 
    { 
        for (int j = 0; j < mat.col(); j++) 
        { 
            cout << mat[{i,j}]<< " "; 
        }
        
    // Newline for new row 
    cout << endl; 
    }   
    */
    

    return 0;
}

/*


int main(int argc, char* argv[])
{
    
    // Verification of 1d system matrix
    Heat<1,double> Heat1d(0.3125, 0.1, 3);
    std::cout << Heat1d.getMatrix() << std::endl;
    // Verification of 2d system matrix
    Heat<2,double> Heat2d(0.3125, 0.1, 3);
    std::cout << Heat2d.getMatrix() << std::endl;
    return 0;
}
*/



           
        
        // Ina's implementation
        /*for (int i = 0; i < points; i++)
        {
            for (int j = 0; j< points; j++)
            {

                if(i == j)
                {
                    M[i,j] = 1 + alpha*(dt/(dx*dx))*(2*n);
                }
                else
                {
                    for (auto k = 0; k<n-1; k++)
                    {   
                        if (i%pow(m,k) != 0 && (i+1)%pow(m,k) != 0 ) // NOT boundary points
                        {
                            M[i,j + pow(m,k)] = 1 + alpha*(dt/(dx*dx));
                            M[i,j - pow(m,k)] = 1 + alpha*(dt/(dx*dx));
                        }
                        else if (i%pow(m,k) = 0)  // no left neighbour in dimension k
                        {
                            M[i,j + pow(m,k)] = 1 + alpha*(dt/(dx*dx)); // only right neighbour in dim k
                        }
                        else if ((i+1)%pow(m,k) = 0)  // no right neighbour in dimension k
                        {
                            M[i,j - pow(m,k)] = 1 + alpha*(dt/(dx*dx)); // only left neighbour in dim k
                        }
                        
                    }
                    
                }      
            } 
        }

            // smito TEST
            int points = (int) pow(mu,n);
            Matrix<T> M( points, points );
            T dx=1/(mu+1);
            
            
            for (int i = 0; i < points; i++)
            {
                M[{i,i}] = 1 + alpha*(dt/(dx*dx))*(2*n);
                
                for (int j = i+1; j < points; j++)
                {      
                    for(int k = 0; k < n; k++)
                    {
                        if((j-i)==pow(m,k) && !(j%pow(m,k+1) < pow(m,k)) )
                        {
                            M[{i,j}] = -alpha*(dt/(dx*dx));
                            M[{j,i}] = -alpha*(dt/(dx*dx));
                        }
                    }                   
                }
            }
            */
            // Loop magic              }
            
            
            // Loop magic