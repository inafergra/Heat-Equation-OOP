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
                matrixMap.insert( {ij  ,0});
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

        map <pair<int, int>, T> Map() const 
        {
            return matrixMap;
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
        map <pair<int, int>, T>  matmap = lhs.Map();
        Vector<typename std::common_type<T,U>::type> result_vector(lhs.row());
        for(int i = 0; i < lhs.row(); i++)      
        {
            result_vector[i] = 0;
        }
        for(auto it = matmap.begin(); it != matmap.end(); ++it)
        {
            int i = it->first.first;
            int j = it->first.second;
            result_vector[i] += rhs[j]*matmap.at({i,j});
        }
    return result_vector;
    }
}



int main()
{
    const Vector<int> a();
    
    Vector<int> b({1,2});
    Vector<int> c(b);
    Vector<int> f = b + c;
    f=f*8;

    Matrix<int> A(10,2);

    A[{0,1}] =2;
    cout<<A[{0,1}];
    
    Vector<int> d = A*b;
    cout<<d[0];
    return 0;
}
/*

template<typename T>
int cg(const Matrix<T>& A, 
       const Vector<T>& b, 
       Vector<T>&       x, 
       T                tol     = (T)1e-8, 
       int              maxiter = 100)
{
    // Your implementation of the cg function starts here
}

template <int n, typename T>
class Heat
{
    // Your implementation of the heat class starts here
};


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