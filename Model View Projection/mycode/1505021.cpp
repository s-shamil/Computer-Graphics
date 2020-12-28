#include <bits/stdc++.h>
using namespace std;
#define pi acos(-1)

struct Matrix{
    double mat[4][4];

    Matrix(){
        for(int i = 0; i<4; i++){
            for(int j = 0; j<4; j++){
                if(i==j) this->mat[i][j] = 1;
                else this->mat[i][j] = 0;
            }
        }
    }
    Matrix(double arr[4][4]){
        for(int i = 0; i<4; i++){
            for(int j = 0; j<4; j++){
                if(i==j) this->mat[i][j] = arr[i][j];
                else this->mat[i][j] = arr[i][j];
            }
        }
    }
    void print(){
        for(int i = 0; i<4; i++){
            for(int j = 0; j<4; j++){
                cout<<this->mat[i][j]<<' ';
            }
            cout<<endl;
        }
    }
    void takeInput(){
        for(int i = 0; i<4; i++){
            for(int j = 0; j<4; j++){
                cin>>this->mat[i][j];
            }
        }
    }
    Matrix operator*(Matrix A) 
    { 
        int N = 4;
        Matrix B;
        for (int i = 0; i < N; i++) 
        { 
            for (int j = 0; j < N; j++) 
            { 
                B.mat[i][j] = 0; 
                for (int k = 0; k < N; k++) 
                { 
                    B.mat[i][j] += (this->mat[i][k]*A.mat[k][j]);
                } 
            } 
        } 
        return B;
    } 
} typedef matrix;

struct Point{
    //double x, y, z;
    double mat[4];
    Point(){
        mat[0] = mat[1] = mat[2] = mat[3] = 1.0;
    }
    Point(double a, double b, double c){
        mat[0] = a;
        mat[1] = b;
        mat[2] = c;
        mat[3] = 1;
    }
    void print(){
        printf("%0.7lf %0.7lf %0.7lf\n", mat[0], mat[1], mat[2]);
    }
    void takeInput(){
        cin>> this->mat[0] >> this->mat[1] >> this->mat[2] ;
    }
} typedef point;

struct Vector{
    double x, y, z;
    Vector(){
        x = y = z = 1;
    }
    Vector(double a, double b, double c){
        x = a;
        y = b;
        z = c;
    }
    void normalize(){
        double val = sqrt(x*x+y*y+z*z);
        x /= val;
        y /= val;
        z /= val;
    }
    Vector operator+(Vector b){
        Vector c;
        c.x = this->x + b.x;
        c.y = this->y + b.y;
        c.z = this->z + b.z;
        return c;
    }
    Vector operator*(double m){
        Vector c;
        c.x = this->x * m;
        c.y = this->y * m;
        c.z = this->z * m;
        return c;
    }
    void print(){
        printf("this vector: %lf %lf %lf\n", x, y, z);
    }
} typedef vector_;

double dot(vector_ a, vector_ b){
    double val = a.x * b.x + a.y * b.y + a.z * b.z;
    return val;
}

vector_ cross(vector_ v1, vector_ v2){
	vector_ tmp;
	tmp.x = v1.y * v2.z - v1.z * v2.y;
	tmp.y = v1.z * v2.x - v1.x * v2.z;
	tmp.z = v1.x * v2.y - v1.y * v2.x;
	return tmp;
}


vector_ rodrigues(vector_ x, vector_ a, double theta){
    //given in degree - convert to radian 
    theta = (theta*pi)/180;
    vector_ v1, v2, v3;
    v1 = x*cos(theta);
    v2 = a*((1-cos(theta))*dot(a,x));
    v3 = cross(a,x)*sin(theta);
    vector_ ret = v1 + v2 + v3;
    return ret;
} 

point transformPoint(matrix t, point p){
    point pnew;
    pnew.mat[0] = p.mat[0] * t.mat[0][0] + p.mat[1] * t.mat[0][1] + p.mat[2] * t.mat[0][2] + p.mat[3] * t.mat[0][3];
    pnew.mat[1] = p.mat[0] * t.mat[1][0] + p.mat[1] * t.mat[1][1] + p.mat[2] * t.mat[1][2] + p.mat[3] * t.mat[1][3];
    pnew.mat[2] = p.mat[0] * t.mat[2][0] + p.mat[1] * t.mat[2][1] + p.mat[2] * t.mat[2][2] + p.mat[3] * t.mat[2][3];
    pnew.mat[3] = p.mat[0] * t.mat[3][0] + p.mat[1] * t.mat[3][1] + p.mat[2] * t.mat[3][2] + p.mat[3] * t.mat[3][3];
    //making w = 1
    pnew.mat[0] /= pnew.mat[3];
    pnew.mat[1] /= pnew.mat[3];
    pnew.mat[2] /= pnew.mat[3];
    pnew.mat[3] /= pnew.mat[3];
    return pnew;
}

int main(){
    freopen("scene.txt", "r", stdin);
    freopen("stage1.txt", "w", stdout);

    stack<int> lastPush;
    /*
        lastPush keeps track of when the push command is found
        when command == push : size of stack<matrix> S is pushed here
        so when we get a pop command reduce our stack S to this size and pop this once
    */
    double eyeX, eyeY, eyeZ;
    cin>>eyeX>>eyeY>>eyeZ;
    double lookX, lookY, lookZ;
    cin>>lookX>>lookY>>lookZ;
    double upX, upY, upZ;
    cin>>upX>>upY>>upZ;
    double fovy, aspectRatio, near, far;
    cin>>fovy>>aspectRatio>>near>>far;
    
    //initialize empty stack S
    stack<matrix> S;
    matrix I;
    S.push(I);

    while(true){
        string command;
        cin>>command;
        if(command == "triangle"){
            for(int i = 0; i<3; i++) {
                point p;
                p.takeInput();
                point newpoint = transformPoint(S.top(), p);
                newpoint.print();
            }
            cout<<endl;
        }
        else if(command == "translate"){
            double tx, ty, tz;
            cin>>tx>>ty>>tz;
            matrix T;
            T.mat[0][3] = tx;
            T.mat[1][3] = ty;
            T.mat[2][3] = tz;
            S.push(S.top()*T);
        }
        else if(command == "scale"){
            double sx, sy, sz;
            cin>>sx>>sy>>sz;
            matrix T;
            T.mat[0][0] = sx;
            T.mat[1][1] = sy;
            T.mat[2][2] = sz;
            S.push(S.top()*T);
        }
        else if(command == "rotate"){
            double angle, ax, ay, az;
            cin>>angle>>ax>>ay>>az;
            vector_ a(ax, ay, az);
            a.normalize();
            vector_ ii(1,0,0), jj(0,1,0), kk(0,0,1);
            vector_ c1 = rodrigues(ii, a, angle);
            vector_ c2 = rodrigues(jj, a, angle);
            vector_ c3 = rodrigues(kk, a, angle);
            matrix T;
            T.mat[0][0] = c1.x ; T.mat[0][1] = c2.x ; T.mat[0][2] = c3.x ;
            T.mat[1][0] = c1.y ; T.mat[1][1] = c2.y ; T.mat[1][2] = c3.y ;
            T.mat[2][0] = c1.z ; T.mat[2][1] = c2.z ; T.mat[2][2] = c3.z ;
            S.push(S.top()*T);
        }

        else if(command == "push"){
            int curr_sz = S.size();
            lastPush.push(curr_sz);
        }
        else if(command == "pop"){
            int curr_sz = S.size();
            if(lastPush.size()>0){
                int past_sz = lastPush.top();
                lastPush.pop();
                int cnt = curr_sz - past_sz;
                while(cnt--){
                    S.pop();
                }
            }
            else{
                printf("Warning: Nothing to pop!!!\n");
            }
        }
        
        else if(command == "end"){
            break;
        }
    }



    //stage 2 begins
    vector_ l(lookX-eyeX, lookY-eyeY, lookZ-eyeZ);
    l.normalize();
    vector_ u(upX, upY, upZ);
    vector_ r = cross(l, u);
    r.normalize();
    u = cross(r, l);

    matrix T; 
    T.mat[0][3] = -eyeX;
    T.mat[1][3] = -eyeY;
    T.mat[2][3] = -eyeZ;
    
    matrix R;
    R.mat[0][0] = r.x ; R.mat[0][1] = r.y ; R.mat[0][2] = r.z ;
    R.mat[1][0] = u.x ; R.mat[1][1] = u.y ; R.mat[1][2] = u.z ;
    R.mat[2][0] = -l.x; R.mat[2][1] = -l.y; R.mat[2][2] = -l.z;
    
    matrix V = R*T;

    ifstream infile;
    infile.open("stage1.txt");
    ofstream outfile;
    outfile.open("stage2.txt");
    double a, b, c;
    int cnt = 0;
    while(infile>>a){
        infile>>b>>c;

        point p1(a, b, c);
        point p2 = transformPoint(V, p1);
        a = p2.mat[0];
        b = p2.mat[1];
        c = p2.mat[2];
        outfile<<setprecision(7)<<fixed<<a<<" "<<b<<" "<<c<<endl;
        cnt++;
        if(cnt==3){
            cnt = 0;
            outfile<<endl;
        }
    }


    //stage 3 begins

    double fovx = fovy * aspectRatio;
    double t_ = near * tan((fovy/2)*(pi/180));
    double r_ = near * tan((fovx/2)*(pi/180));

    matrix P;
    P.mat[0][0] = near/r_;
    P.mat[1][1] = near/t_;
    P.mat[2][2] = -(far+near)/(far-near);
    P.mat[3][3] = 0;
    P.mat[2][3] = -(2*far*near)/(far-near);
    P.mat[3][2] = -1;
    
    ifstream infile_;
    infile_.open("stage2.txt");
    ofstream outfile_;
    outfile_.open("stage3.txt");
    cnt = 0;
    while(infile_>>a){
        infile_>>b>>c;

        point p1(a, b, c);
        point p2 = transformPoint(P, p1);
        a = p2.mat[0];
        b = p2.mat[1];
        c = p2.mat[2];
        outfile_<<setprecision(7)<<fixed<<a<<" "<<b<<" "<<c<<endl;
        cnt++;
        if(cnt==3){
            cnt = 0;
            outfile_<<endl;
        }
    }
    return 0;
}