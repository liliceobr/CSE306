#ifndef HELP
#define HELP

#include "objects.cpp"

const double gammaConst = 4; // usual: 2.2
const double rgbCorrection = pow(255., (gammaConst-1.)/gammaConst);

Camera cam;
vector<Sphere> scene;

Vector getColor(Vector p, int si, Light light, int depth=10, Vector previous=cam.p);

double random01(){
    return double(rand()) / double(RAND_MAX);
}

void addSphere(double x, double y, double z, double r, double g, double b, double R, Materials m=opaque){
    Vector tmpPoint(x,y,z);
    Vector tmpColor(r,g,b);
    scene.push_back(Sphere(tmpPoint, R, tmpColor, m));
}

void addHollowSphere(double x, double y, double z, double t, double n1, double n2, double R){
    addSphere(x, y, z, t, n1, n2, R, transparent);
    addSphere(x, y, z, t, n2, n1, R-0.5, transparent);
}

void buildScene(){
    addSphere(0, 0, 1000, 255, 0, 255, 940);
    addSphere(0, -1000, 0, 0, 0, 255, 990);
    addSphere(0, 1000, 0, 255, 0, 0, 940);
    addSphere(0, 0, -1000, 0, 255, 0, 940);
    addSphere(-1000, 0, 0, 0, 255, 255, 940);
    addSphere(1000, 0, 0, 255, 255, 0, 940);
    //addSphere(0, 0, 0, 255, 255, 255, 10);
    //addSphere(-23, 0, 0, 100, 100, 100, 10, mirror);
    //addSphere(23,0,0, 0, 1, 1.5, 10, transparent);
    //addHollowSphere(0, 0, 0, 0, 1, 1.5, 10);

    addSphere(0, 0, 25, 255, 255, 255, 10);
    addSphere(-20, 0, 10, 100, 100, 100, 10, mirror);
    //addHollowSphere(-20, 0, 0, 0, 1, 1.5, 10);
    //for transparents color: (transparency (0 opaque, 1 transparent), n1, n2)
    addSphere(10,0,30, 0, 1, 1.5, 10, transparent);

    //mesh.readOBJ("model_cat/cat.obj");
}

Vector getPixCoord(double i, double j){
    double x = i;
    double y = h-j-1;
    return Vector(cam.p[0]+x+0.5-w/2., cam.p[1]+y+0.5-h/2., cam.p[2]-cam.f);
}

Vector intersect(Sphere s, Ray r){
    Vector tmp = r.p-s.p;
    double t = dot(r.d, tmp);
    double delta = t*t;
    delta -= (norm(tmp)*norm(tmp) - s.R*s.R);
    t = -t;

    if (delta<0 || t<0) return NULLVEC;

    if (delta == 0.){
        return r.p + t*r.d;
    }
    
    delta = sqrt(delta);
    if (t<delta){
        t += delta;
    }
    else{
        t -= delta;
    }

    return r.p + t*r.d;
}

struct sphereIpointP{
    int i = -1;
    Vector inter = NULLVEC;
};

sphereIpointP intersectScene(Ray ray, int skip=-1){
    sphereIpointP res;
    double closest = 0;

    for(int k=0; k<int(scene.size()); k++){
        if (k==skip) continue;
        Sphere& tmp = scene[k];
        Vector inter = intersect(tmp, ray);
        if (inter!=NULLVEC){
            double d = norm(cam.p-inter);
            if(d < closest || res.i==-1){
                res.i = k;
                res.inter = inter;
                closest = d;
            }

        }
    }
    return res;
}

double visibility(Vector p, Vector light, int si){
    Vector tmp = light-p;
    double d = norm(tmp);
    tmp = tmp/d;
    Ray beam(p, tmp); 
    double passThrough = 1;
    for(int i=0; i<int(scene.size()); i++){
        if (i==si) continue;
        Vector inter = intersect(scene[i], beam);
        if (inter!=NULLVEC){
            inter = p-inter;
            if(norm(inter) < d){
                switch (scene[i].m) {
                case transparent:
                    passThrough = scene[i].c[0];
                case opaque:    
                default:
                    passThrough = 0.;
                    break;
                }
            }
        }
    }
    return passThrough;
}

Vector lambertian(Vector p, int si, Light light){
    Vector n = normalize(p - scene[si].p);
    Vector tmp = light.p-p;
    double d = norm(tmp);
    tmp = (light.I/(4.*PI*PI*d*d) * visibility(p, light.p, si) * max(dot(n, tmp/d), 0.)) * scene[si].c;
    return min(tmp, Vector(255,255,255));
}

Vector mirrorSurface(Vector p, int si, Light light, int depth, Vector previous){
    Vector omegaI = normalize(p-previous);
    Vector n = normalize(p - scene[si].p);
    Vector omegaR = omegaI - 2 * dot(omegaI, n) * n;
    omegaR = normalize(omegaR);

    Ray ray(p, omegaR);
    sphereIpointP best = intersectScene(ray);

    if (best.i == -1) return NULLVEC;
    return getColor(best.inter, best.i, light, depth, p);
}

Vector gammaCor(Vector color){
    double x = pow(color[0],1./gammaConst) * rgbCorrection;
    double y = pow(color[1],1./gammaConst) * rgbCorrection;
    double z = pow(color[2],1./gammaConst) * rgbCorrection;
    return Vector(x, y, z);
}

Vector intersectSelf(Sphere s, Ray r){
    Vector tmp = r.p-s.p;
    double t = dot(r.d, tmp);
    double delta = t*t;
    delta -= (norm(tmp)*norm(tmp) - s.R*s.R);
    t = -t;

    if (delta<0 || t<0) return NULLVEC;

    if (delta != 0.)
        delta = sqrt(delta);
    t += delta; 
    if (norm(t*r.d) < s.R/1000.) return NULLVEC;
    return r.p + t*r.d;
}

Vector refract(Vector p, int si, Light light, int depth, Vector previous){
    double n1 = scene[si].c[1];
    double n2 = scene[si].c[2];

    Vector omegaI = normalize(p-previous);
    Vector n = normalize(p - scene[si].p);
    double tmpDot = dot(omegaI, n);

    double R = 1. - 4.*n1*n2/(n1*n1+2.*n1*n2+n2*n2);
    R = R + (1-R)*pow(1.-abs(tmpDot), 5.);
    
    if (random01() < R){
        return mirrorSurface(p, si, light, depth+1, previous);
    }
    
    if (tmpDot>0){
        n = -1*n;
        double a = n1;
        n1 = n2;
        n2 = a;
    }
    tmpDot = dot(omegaI, n);
    double n12 = n1/n2;

    double rad = 1-n12*n12*(1-tmpDot*tmpDot);
    if(rad<0){
        return mirrorSurface(p, si, light, depth+1, previous);
    }

    Vector omegaT = n12*(omegaI-tmpDot*n);
    omegaT = omegaT - n*sqrt(rad);
    omegaT = normalize(omegaT);
    
    Ray ray(p, omegaT);
    sphereIpointP best = intersectScene(ray, si);
    Vector inter = intersectSelf(scene[si], ray);

    if (inter!=NULLVEC){
        if (best.i == -1) return getColor(inter, si, light, depth, p);
        double d1 = norm(p-best.inter);
        double d2 = norm(p-inter);
        if(d2 < d1){
            return getColor(inter, si, light, depth, p);
        }
    }
    
    if (best.i == -1)
        return NULLVEC;
    
    return getColor(best.inter, best.i, light, depth, p);
}

Vector randomVec(const Vector& n){
    const double r1 = random01();
    const double r2 = random01();
    const double x = cos(2.*PI*r1)*sqrt(1.-r2);
    const double y = sin(2.*PI*r1)*sqrt(1.-r2);
    const double z = sqrt(r2);
    Vector T1;
    if(abs(n[0])<abs(n[1])){
        if(abs(n[0])<abs(n[2])){
            T1 = Vector(0, -n[2], n[1]);
        } else{
            T1 = Vector(-n[1], n[0], 0);
        }
    } else{
        if(abs(n[1])<abs(n[2])){
            T1 = Vector(-n[2], 0, n[0]);
        } else{
            T1 = Vector(-n[1], n[0], 0);
        }
    }
    T1 = normalize(T1);
    Vector T2 = normalize(cross(n,T1));
    return normalize(x*T1 + y*T2 + z*n);
}

Vector indirectLight(Vector p, int si, Light light, int depth){
    if (depth<=0) return NULLVEC;
    depth -= 1;
    Vector diffusion(0,0,0);
    Vector n = normalize(p - scene[si].p);

    Vector ranedVec = randomVec(n);
    Ray ray(p, ranedVec);
    sphereIpointP tmpbest = intersectScene(ray);
    if (tmpbest.i != -1)                                
        diffusion += getColor(tmpbest.inter, tmpbest.i, light, depth, p);
    return diffusion;
}

void boxMuller(double& x, double& y, double stdev=0.3){
    const double r1 = random01();
    const double r2 = random01();
    x = x+sqrt(-2*log(r1))*cos(2*PI*r2)*stdev;
    y = y+sqrt(-2*log(r1))*sin(2*PI*r2)*stdev;
}

Vector sphericalLight(Vector p, int si, Light light){
    Vector n = normalize(p - scene[si].p);
    Vector tmp = p-light.p;
    double d = norm(tmp);
    Vector nPrime = randomVec(tmp/d);
    Vector xPrime = light.p + light.R*nPrime;
    Vector omegaI = normalize(xPrime - p);
    double vis = visibility(p, xPrime, si);
    double R = light.R;
    double pdf = dot(nPrime, tmp/d)/PI/R/R;

    if (pdf==0) return NULLVEC; 
    tmp = light.I/(4.*PI*PI*R*R) * scene[si].c/PI * vis * max(dot(n, omegaI), 0.) * max(dot(nPrime, -1*omegaI), 0.) / (pow(norm(xPrime-p), 2)*pdf);
    return tmp;
}

Vector getColor(Vector p, int si, Light light, int depth, Vector previous){
    if (p==previous || depth<=0 || p!=p) return NULLVEC;
    depth -= 1;
    switch (scene[si].m){
    case opaque: // lambertian or sphericalLight
        //return (lambertian(p, si, light) + scene[si].c * indirectLight(p, si, light, depth-1) / 255) / 2;
        return (sphericalLight(p, si, light) + scene[si].c * indirectLight(p, si, light, depth-1) / 255) / 2;
    case mirror:
        return mirrorSurface(p, si, light, depth, previous);
    case transparent:{
        return refract(p, si, light, depth, previous);}
    default:
        return NULLVEC;
    }
}

#endif