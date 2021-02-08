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

        /*
        Vector(const Vector<T>& other) // Copy constructor
        {
            length = other.length;
            delete[] data;
            data = new T[other.length];
            
            for(int i=0;i<other.length;i++)
            {
                data[i] = other.data[i];
            }
        }
        */

        Vector(const Vector<T>& other) // Copy constructor
        :Vector(other.length)
        {
            for(int i = 0; i<other.length; i++)
            {
                data[i] = other.data[i];
            }    
        }
        
        Vector(Vector<T>&& other) // Move constructor
        :length(other.length), data(other.data)
        {
            other.length = 0;
            other.data = nullptr;  
        }


        Vector(const initializer_list<T>& list)
        :Vector((int)list.size())
        {
            uninitialized_copy(list.begin(), list.end(), data);
        }

        ~Vector() {delete[] data;length=0;data=nullptr;}

        Vector<T>& operator=(const Vector<T>& other) // Copy operator
        {   
            if (this != &other)
            {
                length = other.length;
                delete[] data;
                data = new T[other.length];           
                for(int i=0;i<other.length;i++)
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
            return data[i]; 
        }

        const T& operator[](int i) const
        {   
            return data[i];
        }
        
        
};


int main()
{
    const Vector<int> a();
    
    Vector<int> b({1,2});
    Vector<int> c(b);
    Vector<int> d(move(b));
    return 0;
}

/*
template<typename T, typename U>
typename std::common_type<T,U>::type 
dot(const Vector<T>& lhs, 
    const Vector<U>& rhs)
{
    // Your implementation of the dot function starts here
}

template <typename T>
class Matrix
{
    // Start your implementation of the matrix class here
};

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