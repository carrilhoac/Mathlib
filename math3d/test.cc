
#include "math3d.h"

int main(int argc, char **argv)
{
    math::Vector4d a(0,0,1,0);
    math::Vector4d b(0,0,0,0);

    std::cout << (a > b) << std::endl;
    std::cout << (a >= b) << std::endl;
    //std::cout << (b > a) << std::endl;

#if 0
    math::Matrix4d A;

    A.m[ 0]= 1; A.m[ 1]= 2; A.m[ 2]= 3; 
    A.m[ 4]= 2; A.m[ 5]= 4; A.m[ 6]= 5; 
    A.m[ 8]= 3; A.m[ 9]= 5; A.m[10]= 6; 

    math::Vector4d e;
    A.eigval3(e);

    std::cout << "Matrix:" << std::endl;
    std::cout << A << std::endl;
    std::cout << "eigen-values:" << std::endl;
    std::cout << e << std::endl;
#endif 
#if 0
    math::Matrix4d A;

    A.m[ 0]= 2; A.m[ 1]= 1; A.m[ 2]=-3; A.m[ 3]= 4;
    A.m[ 4]=-1; A.m[ 5]= 0; A.m[ 6]= 2; A.m[ 7]= 5;
    A.m[ 8]= 3; A.m[ 9]= 2; A.m[10]= 1; A.m[11]= 0;
    A.m[12]= 4; A.m[13]=-2; A.m[14]= 3; A.m[15]= 1;

    math::Matrix4d B(A);
    math::Matrix4d C;

    std::cout << A << std::endl;
    B.invert();
    std::cout << B << std::endl;
    A *= B;
    C.identity();
    A *= C;
    std::cout << A << std::endl;
#endif 

#if 0
    math::Matrix4d A;
    A.identity();
    A.m[ 3] = 2;
    A.m[ 7] = 3;
    A.m[11] = 4;

    math::Vector4d a(1,1,1,1);
    math::Vector4d b;

    A.invert();
    A.transform(b, a);

    std::cout << A << std::endl;
    std::cout << a << std::endl;
    std::cout << b << std::endl;
#endif 

#if 0
    math::Vector4d a(1,0,0);
    math::Vector4d b(0,1,0);
    math::Vector4d c;
    math::cross(c, b, a);
    std::cout << a << std::endl;
    std::cout << b << std::endl;
    std::cout << c << std::endl;
#endif 

#if 0
    math::Matrix4f A;
    A.identity();
    A.m[ 3] = 2;
    A.m[ 7] = 3;
    A.m[11] = 4;

    math::Vector4f a(1,1,1,1);
    math::Vector4f b;

    A.invert();
    A.transform(b, a);

    std::cout << A << std::endl;
    std::cout << a << std::endl;
    std::cout << b << std::endl;
#endif 
    return 0;
}
