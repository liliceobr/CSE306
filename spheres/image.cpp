#include "helpers.cpp"

void writePPM(int w, int h, int* image){
    ofstream output;
    output.open("pic_cool.ppm");
    output << "P3" << "\n";
    output << w << " " << h << "\n";
    output << "255\n";

    for(int i=0; i<h; i++){
        for(int j=0; j<w; j++){
            int tmp = i*w*3 + j*3;
            output << image[tmp] << " " << image[tmp+1] << " " << image[tmp+2] << "\n";
        }
    }

    output.close();
}

void slicedOptimalTransport(int (&input)[], int model[]){
    const double r1 = random01();
    const double r2 = random01();
    const double x = cos(2.*PI*r1)*sqrt(r2-r2*r2);
    const double y = sin(2.*PI*r1)*sqrt(r2-r2*r2);
    const double z = 1.-2.*r2;
    const Vector randVec(x, y, z);

    int wh = sizeList(input)/3;
    if (sizeList(model)/3 != wh){
        print("Incompatible sizes");
        return;
    }

    vector<pair<double, int> > pairInput(wh);
    vector<pair<double, int> > pairModel(wh); //does not need to be a pair, vector<double> should be enough

    #pragma omp parallel for
    for(int i=0; i<wh; i++){
        Vector colorInput(input[i*3+0], input[i*3+1], input[i*3+2]);
        pairInput[i] = make_pair(dot(colorInput, randVec), i);

        Vector colorModel(model[i*3+0], model[i*3+1], model[i*3+2]);
        pairModel[i] = make_pair(dot(colorModel, randVec), i);
    }

    sort(pairInput.begin(), pairInput.end());
    sort(pairModel.begin(), pairModel.end());

    #pragma omp parallel for
    for(int i=0; i<wh; i++){
        input[pairInput[i].second + 0] += (pairModel[i].first - pairInput[i].first) * randVec[0];
        input[pairInput[i].second + 0] = min(255, max(0, int(input[pairInput[i].second + 0])));
        input[pairInput[i].second + 1] += (pairModel[i].first - pairInput[i].first) * randVec[1];
        input[pairInput[i].second + 1] = min(255, max(0, int(input[pairInput[i].second + 1])));
        input[pairInput[i].second + 2] += (pairModel[i].first - pairInput[i].first) * randVec[2];
        input[pairInput[i].second + 2] = min(255, max(0, int(input[pairInput[i].second + 2])));
    }
}

void retargeting(int image[], int (&res)[], int width=w, int height=h){
    int intermediate[width*height];

    // setup of E(x,y)
    for (int i=0; i<width; i++){
        for (int j=0; j<height; j++){
            intermediate[j*width+i]  = std::abs(image[((j+1)*width+i)*3+0]-image[((j-1)*width+i)*3+0] + image[((j+1)*width+i)*3+1]-image[((j-1)*width+i)*3+1] + image[((j+1)*width+i)*3+2]-image[((j-1)*width+i)*3+2]);
            intermediate[j*width+i] += std::abs(image[(j*width+(i+1))*3+0]-image[(j*width+(i-1))*3+0] + image[(j*width+(i+1))*3+1]-image[(j*width+(i-1))*3+1] + image[(j*width+(i+1))*3+2]-image[(j*width+(i-1))*3+2]);
        }
    }

    // calculate C(x,y)
    for (int i=0; i<width; i++){
        for (int j=1; j<height; j++){
            if (i==0)
                intermediate[j*width+i] += min(intermediate[(j-1)*width+i], intermediate[(j-1)*width+i+1]);
            else if (i==width-1)
                intermediate[j*width+i] += min(intermediate[(j-1)*width+i], intermediate[(j-1)*width+i-1]);
            else
                intermediate[j*width+i] += min(min(intermediate[(j-1)*width+i], intermediate[(j-1)*width+i-1]), intermediate[(j-1)*width+i+1]);
        }
    }

    // find the lowest energy
    int minimum = 0;
    int tmpMinC = intermediate[(height-1)*width+0];
    for (int i=1; i<width; i++){
        if(intermediate[(height-1)*width+i]<tmpMinC){
            minimum = i;
            tmpMinC = intermediate[(height-1)*width+i];
        }
    }

    // remove seam at the same time as finding it
    for (int j=height-1; j>=0; j--){
        bool passed = false;
        for (int i=0; i<width; i++){
            if (!passed && i==minimum) {
                if (i==0)
                    minimum = intermediate[(j-1)*width+i]<intermediate[(j-1)*width+i+1]?i:i+1;
                else if (i==width-1)
                    minimum = intermediate[(j-1)*width+i]<intermediate[(j-1)*width+i-1]?i:i-1;
                else{
                    if (intermediate[(j-1)*width+i-1] < intermediate[(j-1)*width+i])
                        minimum = intermediate[(j-1)*width+i-1]<intermediate[(j-1)*width+i+1]?i-1:i+1;
                    else if (intermediate[(j-1)*width+i]<intermediate[(j-1)*width+i+1])
                        minimum = i;
                    else
                        minimum = i+1;
                }
                passed = true;
                continue;
            }
            
            if (!passed){
                res[((j+1)*width+i)*3+0] = image[((j+1)*width+i)*3+0];
                res[((j+1)*width+i)*3+1] = image[((j+1)*width+i)*3+1];
                res[((j+1)*width+i)*3+2] = image[((j+1)*width+i)*3+2];            
            }
            else{
                res[((j+1)*width+i-1)*3+0] = image[((j+1)*width+i)*3+0];
                res[((j+1)*width+i-1)*3+1] = image[((j+1)*width+i)*3+1];
                res[((j+1)*width+i-1)*3+2] = image[((j+1)*width+i)*3+2];  
            }
        }
    }
}