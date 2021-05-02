#ifndef VECTOR
#define VECTOR

#include "master.cpp"

class Vector {
public:
    explicit Vector(double coords[3]){
        coord[0] = coords[0]; coord[1] = coords[1]; coord[2] = coords[2];
    }
    explicit Vector(double x = 0., double y = 0., double z=0.){
        coord[0] = x; coord[1] = y; coord[2] = z;
    }
    Vector& operator+=(const Vector& b){
        coord[0] += b[0];
        coord[1] += b[1];
        coord[2] += b[2];
        return *this;
    }
    double coord[3];  
    double& operator[](int i) {return coord[i];}  
    const double& operator[](int i) const {return coord[i];}

};

const Vector NULLVEC(0,0,0);

Vector operator+(const Vector& a, const Vector &b){
    return Vector(a[0]+b[0], a[1]+b[1], a[2]+b[2]);
}

Vector operator-(const Vector& a, const Vector &b){
    return Vector(a[0]-b[0], a[1]-b[1], a[2]-b[2]);
}

Vector operator-(const Vector& a, double t){
    return Vector(a[0]-t, a[1]-t, a[2]-t);
}

Vector operator*(double t, const Vector &b){
    return Vector(t*b[0], t*b[1], t*b[2]);
}

Vector operator*(const Vector &b, double t){
    return t*b;
}

Vector operator*(const Vector &a, const Vector &b){
    return Vector(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}

Vector operator/(const Vector &b, double t){
    return Vector(b[0]/t, b[1]/t, b[2]/t);
}

double dot(const Vector& a, const Vector& b){
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

Vector cross(const Vector& a, const Vector& b){
    return Vector(a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]);
}

double norm(const Vector& v){
    return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

double norm2(const Vector& v){
    return v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
}

Vector normalize(const Vector& v){
    if (norm(v)==0)
        print("Normalise null vector (division by 0)");
    return v/norm(v);
}

Vector min(const Vector &a, const Vector &b){
    double x = min(a[0], b[0]);
    double y = min(a[1], b[1]);
    double z = min(a[2], b[2]);
    return Vector(x,y,z);
}

Vector max(const Vector &a, const Vector &b){
    double x = max(a[0], b[0]);
    double y = max(a[1], b[1]);
    double z = max(a[2], b[2]);
    return Vector(x,y,z);
}

void print(const Vector &i){
    std::cout<<'('<<i[0]<<','<<i[1]<<','<<i[2]<<')'<<std::endl;
}

bool operator==(const Vector &a, const Vector &b){
    return a[0]==b[0]&&a[1]==b[1]&&a[2]==b[2];
}

bool operator!=(const Vector &a, const Vector &b){
    return a[0]!=b[0]||a[1]!=b[1]||a[2]!=b[2];
}

#endif