#include "svg.cpp"

double rand01(){
    return double(rand()) / double(RAND_MAX);
}

double gauss(Vector x, Vector m, double s) {
    return  exp(-norm2(m-x)/(2 * s * s)) / (sqrt(2 * M_PI) * s);
}


std::vector<Vector> tutte(std::vector<Vector> points, std::vector<std::vector<int>> adj, std::vector<int> boundaryId) {
    int n = boundaryId.size();
    
    double s = 0;
    for (int i=0; i<n; i++)
        s += norm(points[boundaryId[i]] - points[boundaryId[(i+1)<n?i+1:0]]);
    
    double cs = 0;
    std::vector<Vector> v;

    for (int i=0; i<n; i++){
        double theta = 2*PI*cs/s;
        points[boundaryId[i]] = Vector(cos(theta), sin(theta), 0);
        cs += norm(points[boundaryId[i]] - points[boundaryId[(i+1)<n?i+1:0]]);
    }

    for (int j=0; j<100; j++){
        std::vector<Vector> tmp = points;
        for(int i = 0; i < n; i ++) {
            tmp[i] = NULLVEC;
            for(int k=0; k<adj[i].size(); k++) 
                tmp[i] += points[adj[i][k]];
            tmp[i] = tmp[i] / double(adj[i].size());
        } 
        for(int i=0; i<n; i++) tmp[boundaryId[i]] = points[boundaryId[i]];
        points = tmp;
    }
    return points;
}

int main(int argc, char *argv[]){
    //srand(0);

    int N = 100;
    std::vector<Vector> points(N);
    std::vector<double> lambdas(N);
    Vector middle = Vector(0.5, 0.5, 0);
    double sigma = 0.1;
    for(int i = 0; i < N; i ++) {
        points[i] = Vector(rand01(), rand01(), 0);
        //points[i] = Vector(rand01(), rand01(), rand01());
        lambdas[i] = gauss(points[i], middle, sigma);
        points[i][2] = lambdas[i];
    }
    std::vector<Polygon> voronoisDiag = voronoiParallel(points);
    print("Done");
    save_svg(voronoisDiag, "res100_power.svg");
}