#include "../spheres/vector.cpp"
#include "lbfgs.h"

// Polygon class using Vector
class Polygon{
public:
    std::vector<Vector> p;
    Polygon(){}
    Polygon(std::vector<Vector> points){
        p = points;
    }
    Vector& operator[](int i){return p[i];}  
    const Vector& operator[](int i) const{return p[i];}
    const unsigned int size(){return p.size();}
    void add(Vector tmp){p.push_back(tmp);}
};

std::vector<Polygon> voronoiParallel(std::vector<Vector> points, Polygon boundBox){
    const double e = 10e-6;
    std::vector<Polygon> set(points.size());
    
    // parallel threading
    #pragma omp parallel for
    for(int ic = 0; ic < points.size(); ic++){
        Vector center = points[ic];
        Polygon cell = boundBox;

        for(int ip = 0; ip < points.size(); ip++){
            if(ic == ip) continue;

            Vector N = normalize(points[ip] - center);
            double t = dot(center + ((points[ip] - center) / 2), N);

            int n = cell.size();
            std::vector<double> ts(n);
            for(int i = 0; i < n; i ++) ts[i] = dot(N, cell[i]) - t;

            Polygon ncell;
            for(int i = 0; i < n; i ++){
                int j = (i+1)%n;
                int inter = 0;
                if(ts[i] < e){
                    ncell.add(cell[i]);
                    if(ts[j] > -e) 
                        inter = 1;
                } else if(ts[j] < e) 
                    inter = -1;

                if(inter){
                    double ijN = dot((cell[j] - cell[i]), N);
                    ncell.add((std::abs(ijN) < e)? cell[(inter==1)?j:i] : cell[i]+(cell[j]-cell[i])*(t-dot(cell[i],N))/ijN);
                }
            }
            for(int i = 0; i < ncell.size(); i ++) 
                ncell[i][2] = 0;
            cell = ncell;
            
        }
        set[ic] = cell;
    }
    return set;
}

std::vector<Polygon> voronoiParallel(std::vector<Vector> points){
    Vector boxMin = points[0];
	Vector boxMax = points[0];
    for(int i=1; i<int(points.size()); i++){
        boxMin = min(boxMin, points[i]);
	    boxMax = max(boxMax, points[i]);
    }
    std::vector<Vector> tmp = {boxMin, Vector(boxMin[0], boxMax[1], 0), boxMax, Vector(boxMax[0], boxMin[1], 0)};
    Polygon boundBox(tmp);

    return voronoiParallel(points, boundBox);
}

double area(Polygon p){
    if(p.size()==0) return 0;
    double A = 0;
    for(int i = 1; i < p.size() - 1; i ++)
        A += cross(p[i] - p[0], p[i+1] - p[0])[2];

    return A;
}

double distance(Polygon p, Vector point){
    if(p.size()==0) return 0;
    double A = 0;
    std::vector<Vector> tmp ={p[0], NULLVEC, NULLVEC};
    Polygon c(tmp);
    for(int i = 1; i < p.size() - 1; i ++){
        double x = 0;
        c[1] = p[i];
        c[2] = p[i+1];

        for(int k = 0; k < 3; k ++){
            for(int l = k; l < 3; l ++)
                x += dot(c[k] - point, c[l] - point);
        }
        //cross product from vector.cpp
        A += std::abs(cross(c[0] - c[1], c[2] - c[1])[2]) / 6 * x;
    }
    return A;
}


class Instance{
public:
    std::vector<Vector> p;
    std::vector<double> l;
    Polygon box;
    Instance(std::vector<Vector> points, std::vector<double> lambda, Polygon boundBox){
        p = points;
        l = lambda;
        box = boundBox;
    }
};

lbfgsfloatval_t evaluate(void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step){
    std::vector<Vector> points = ((Instance *)instance)->p;
    std::vector<double> lambda = ((Instance *)instance)->l;
    Polygon boundBox = ((Instance *)instance)->box;

    for(int i = 0; i < n; i++) 
        points[i][2] = x[i];
    
    std::vector<Polygon> voronoisDiag = voronoiParallel(points, boundBox);

    for(int i = 0; i < n; i++) 
        points[i][2] = 0;
    
    double sum = 0;
    lbfgsfloatval_t fx = 0.;
    for(int i = 0; i < n; i++){
        double A = std::abs(area(voronoisDiag[i]));
        sum += A;
        g[i] = A - lambda[i];
        double aT = distance(voronoisDiag[i], points[i]);
        fx += x[i] * g[i] - std::abs(distance(voronoisDiag[i], points[i]));
    }
    return fx;
}