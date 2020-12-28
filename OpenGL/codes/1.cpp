#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include <bits/stdc++.h>

//#include <windows.h>
#include <GL/glut.h>
#include <GL/glx.h>
#include <GL/gl.h>

using namespace std;

#define pi (2*acos(0.0))
#define theta pi/180 //amount of angle rotation per click
#define delta 0.5 //change in radius of sphere-cube

double cameraHeight;
double cameraAngle;
int drawaxes;
double angle;
double a, r; //for sphere-cube

struct point
{
	double x,y,z;

	point(double a, double b, double c){
		x = a; y = b; z = c;
	}

	point(){
		x = y = z = 0;
	}

	void operator+=(struct point const &p){
		this->x += p.x;
		this->y += p.y;
		this->z += p.z;
	}

	void operator-=(struct point const &p){
		this->x -= p.x;
		this->y -= p.y;
		this->z -= p.z;
	}
	struct point operator+(struct point const &p){
		struct point tmp;
		tmp.x = this->x + p.x;
		tmp.y = this->y + p.y;
		tmp.z = this->z + p.z;
		return tmp;
	}
	struct point operator*(double r){
		struct point tmp;
		tmp.x = this->x * r;
		tmp.y = this->y * r;
		tmp.z = this->z * r;
		return tmp;
	}
} typedef Point;

Point typedef Vector;
Point pos; //camera position
Vector uvec, rvec, lvec; // up, right, lool (not left)

void initCamParams(){
	pos = Point(100, 100, 0);
	double x = cos(pi/4);
	uvec = Vector(0, 0, 1);
	rvec = Vector(-x, x, 0);
	lvec = Vector(-x, -x, 0);
}

Vector crossProduct(Vector v1, Vector v2){
	Vector tmp;
	tmp.x = v1.y * v2.z - v1.z * v2.y;
	tmp.y = v1.z * v2.x - v1.x * v2.z;
	tmp.z = v1.x * v2.y - v1.y * v2.x;
	return tmp;
}

Vector rotateVector(Vector axis, Vector v, double rotAngle, int ccw){
	Vector perp;
	if(ccw==1){
		//counterclock-wise rotation
		perp = crossProduct(axis, v);
	}
	else{
		//clockwise rotation
		perp = crossProduct(v, axis);
	}
	Vector scaled_v = v*cos(rotAngle);
	Vector scaled_perp = perp*sin(rotAngle);
	Vector new_v = scaled_v + scaled_perp;
	return new_v;
}

void drawAxes()
{
	if(drawaxes==1)
	{
		glColor3f(1.0, 1.0, 1.0);
		glBegin(GL_LINES);{
			glVertex3f( 100,0,0);
			glVertex3f(-100,0,0);

			glVertex3f(0,-100,0);
			glVertex3f(0, 100,0);

			glVertex3f(0,0, 100);
			glVertex3f(0,0,-100);
		}glEnd();
	}
}


void drawSquare(double a)
{
    //glColor3f(1.0,0.0,0.0);
	glBegin(GL_QUADS);{
		glVertex3f( a, a,0);
		glVertex3f( a,-a,0);
		glVertex3f(-a,-a,0);
		glVertex3f(-a, a,0);
	}
	glEnd();
}

void drawCylBy4(double radius,double height,int segments)
{
    int i;
    struct point points[100];
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*(pi/2));
        points[i].y=radius*sin(((double)i/(double)segments)*(pi/2));
    }
    //draw quads using generated points
    glColor3f(0,1,0);

	for(i=0;i<segments;i++)
    {
        glBegin(GL_QUADS);
        {
			glVertex3f(points[i].x,points[i].y,height);
			glVertex3f(points[i].x,points[i].y,-height);
			glVertex3f(points[i+1].x,points[i+1].y,-height);
			glVertex3f(points[i+1].x,points[i+1].y,height);
        }
        glEnd();
    }

	glColor3f(1,1,1);

}

void drawSphereBy8(double radius,int slices,int stacks)
{
	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*(pi/2));
			points[i][j].y=r*sin(((double)j/(double)slices)*(pi/2));
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);{
			    //upper hemisphere
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                //lower hemisphere
                /*
                glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
                */
			}glEnd();
		}
	}
}

void draw6Squares(double a, double r){
	double dist = a+r;
	double sidelen = a;
	//cout<<dist<<endl;
	
	glPushMatrix();
	glTranslatef(0,0,dist);
	drawSquare(sidelen);
	glPopMatrix();

	glPushMatrix();
	glTranslatef(0,0,-dist);
	drawSquare(sidelen);
	glPopMatrix();
	
	glPushMatrix();
	glTranslatef(0,dist,0);
	glRotatef(90, 1, 0, 0);
	drawSquare(sidelen);
	glPopMatrix();

	glPushMatrix();
	glTranslatef(0,-dist,0);
	glRotatef(90, 1, 0, 0);
	drawSquare(sidelen);
	glPopMatrix();	

	glPushMatrix();
	glTranslatef(dist,0,0);
	glRotatef(90, 0, 1, 0);
	drawSquare(sidelen);
	glPopMatrix();

	glPushMatrix();
	glTranslatef(-dist,0,0);
	glRotatef(90, 0, 1, 0);
	drawSquare(sidelen);
	glPopMatrix();	
}

void draw4Spheres(double a, double r){
	glPushMatrix();
	glTranslatef(a,a,a);
	drawSphereBy8(r, 30, 30);
	glPopMatrix();

	glPushMatrix();
	glTranslatef(a,a,-a);
	glRotatef(-90,1,0,0);
	drawSphereBy8(r, 30, 30);
	glPopMatrix();

	glPushMatrix();
	glTranslatef(a,-a,a);
	glRotatef(90,1,0,0);
	drawSphereBy8(r, 30, 30);
	glPopMatrix();

	glPushMatrix();
	glTranslatef(a,-a,-a);
	glRotatef(180,1,0,0);
	drawSphereBy8(r, 30, 30);
	glPopMatrix();
}

void draw8Spheres(double a, double r){
	glColor3f(1,0,0);

    draw4Spheres(a, r);

	glPushMatrix();
	{
		glRotatef(180,0,0,1);

        draw4Spheres(a, r);
	}
	glPopMatrix();
	
	glColor3f(1,1,1);
}

void draw4Cylinders(double rad, double hei){
    glPushMatrix();
	glTranslatef(hei, hei, 0);
	drawCylBy4(rad, hei, 40);
	glPopMatrix();

	glPushMatrix();
	glTranslatef(-hei, hei, 0);
	glRotatef(90, 0, 0, 1);
	drawCylBy4(rad, hei, 40);
	glPopMatrix();

	glPushMatrix();
	glTranslatef(-hei, -hei, 0);
	glRotatef(180, 0, 0, 1);
	drawCylBy4(rad, hei, 40);
	glPopMatrix();

	glPushMatrix();
	glTranslatef(hei, -hei, 0);
	glRotatef(-90, 0, 0, 1);
	drawCylBy4(rad, hei, 40);
	glPopMatrix();
}
void draw12Cylinders(double rad, double hei){
    draw4Cylinders(rad, hei);

	glPushMatrix();
	{
		glRotatef(90, 1, 0, 0);
		draw4Cylinders(rad, hei);
	}
	glPopMatrix();

	glPushMatrix();
	{
		glRotatef(90, 0, 1, 0);
        draw4Cylinders(rad, hei);
	}
	glPopMatrix();
}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){

		case '1':
			lvec = rotateVector(uvec, lvec, theta, 1);
			rvec = rotateVector(uvec, rvec, theta, 1);
			break;
		case '2':
			lvec = rotateVector(uvec, lvec, theta, 0);
			rvec = rotateVector(uvec, rvec, theta, 0);
			break;
		case '3':
			lvec = rotateVector(rvec, lvec, theta, 1);
			uvec = rotateVector(rvec, uvec, theta, 1);
			break;
		case '4':
			lvec = rotateVector(rvec, lvec, theta, 0);
			uvec = rotateVector(rvec, uvec, theta, 0);
			break;
		case '5':
			uvec = rotateVector(lvec, uvec, theta, 1);
			rvec = rotateVector(lvec, rvec, theta, 1);
			break;
		case '6':
			uvec = rotateVector(lvec, uvec, theta, 0);
			rvec = rotateVector(lvec, rvec, theta, 0);
			break;

		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			cameraHeight -= 3.0;

			//move backward
			pos -= lvec;
			break;
		case GLUT_KEY_UP:		// up arrow key
			cameraHeight += 3.0;

			//move forward
			pos += lvec;
			break;

		case GLUT_KEY_RIGHT:
			cameraAngle += 0.03;

			//move right
			pos += rvec;
			break;

		case GLUT_KEY_LEFT:
			cameraAngle -= 0.03;

			//move right
			pos -= rvec;
			break;

		case GLUT_KEY_PAGE_UP:

			//move up
			pos += uvec;
			break;

		case GLUT_KEY_PAGE_DOWN:

			//move down
			pos -= uvec;
			break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
			if(a>0){
				a -= delta;
				r += delta;
			}
			break;
		case GLUT_KEY_END:
			if(r>0){
				a += delta;
				r -= delta;
			}
			break;

		default:
			break;
	}
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
				drawaxes=1-drawaxes;
			}
			break;

		case GLUT_RIGHT_BUTTON:
			//........
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}



void display(){

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();

	//now give three info
	//1. where is the camera (viewer)?
	//2. where is the camera looking?
	//3. Which direction is the camera's UP direction?

	//gluLookAt(100,100,100,	0,0,0,	0,0,1);
    //gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
	//gluLookAt(0,0,200,	0,0,0,	0,1,0);
	gluLookAt(pos.x, pos.y, pos.z, pos.x + lvec.x, pos.y + lvec.y, pos.z + lvec.z, uvec.x, uvec.y, uvec.z);
	//cout<<pos.x << ' '<<pos.y<<' '<<pos.z<<endl;
	//cout<<lvec.x << ' '<<lvec.y<<' '<<lvec.z<<endl;

	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

	drawAxes();

    //glColor3f(1,0,0);
    //drawSquare(10);

    //glTranslatef(30, 0, 0);
	//drawSquare(40);


	draw6Squares(a,r);
	draw8Spheres(a,r);
	draw12Cylinders(r, a);

	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}


void animate(){
	angle+=0.05;
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	//codes for initialization
	drawaxes=1;
	cameraHeight=150.0;
	cameraAngle=1.0;
	angle=0;
    initCamParams();
	a = 20.0, r = 5.0; //sphere-cube

	//clear the screen
	glClearColor(0,0,0,0);

	/************************
	/ set-up projection here
	************************/
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	gluPerspective(80,	1,	1,	1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}

int main(int argc, char **argv){
	//initCamParams();
	/*
	//cross product test
	Vector ek(1,-2,3), dui(4,5,-6);
	Vector tin = crossProduct(ek, dui);
	cout<<tin.x<<" "<<tin.y<<" "<<tin.z<<endl;
	*/
	/*
	//rotation test
	Vector xx(1,0,0), yy(0,1,0);
	Vector zz = rotateVector(yy, xx, pi/2, 1);
	printf("%0.2lf %0.2lf %0.2lf \n", zz.x,zz.y,zz.z);
	*/
	glutInit(&argc,argv);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");

	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}
