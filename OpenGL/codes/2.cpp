#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <iostream>
using namespace std;

//#include <windows.h>
#include <GL/glut.h>

#define pi (2*acos(0.0))
#define rotationAngleDelta 5.0

double cameraHeight;
double cameraAngle;
int drawgrid;
double angle;
int wheelRadius, wheelWidth;
double yrotation, zrotation;
/*
    yrotation and zrotation controls the wheel orientation
     zrotation (z axis wrt) -> D (cw), A (ccw)
     yrotation (y axis wrt) -> W (cw), S (ccw)
*/
double linDistPerClick; //the distance of COM of wheel by one click

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

Point centerPos;



void drawGrid()
{
	int i;
	if(drawgrid==1)
	{
		glColor3f(0.6, 0.6, 0.6);	//grey
		glBegin(GL_LINES);{
			for(i=-10;i<=10;i++){

				//if(i==0)
				//	continue;	//SKIP the MAIN axes

				//lines parallel to Y-axis
				glVertex3f(i*10, -120, 0);
				glVertex3f(i*10,  120, 0);

				//lines parallel to X-axis
				glVertex3f(-120, i*10, 0);
				glVertex3f( 120, i*10, 0);
			}
		}glEnd();
	}
}

void drawWheel(double radius,double height,int segments)
{
    int i;
    double shade;
    struct point points[100];
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].z=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw triangles using generated points
    for(i=0;i<segments;i++)
    {
        //create shading effect
        if(i<segments/2)shade=2*(double)i/(double)segments;
        else shade=2*(1.0-(double)i/(double)segments);
        glColor3f(shade,shade,shade);

        glBegin(GL_QUADS);
        {
            glVertex3f(points[i].x, -(height/2),points[i].z);
            glVertex3f(points[i+1].x, -(height/2),points[i+1].z);
            glVertex3f(points[i+1].x, (height/2),points[i+1].z);
            glVertex3f(points[i].x, (height/2),points[i].z);
        }
        glEnd();
    }
    //spikes
    glColor3f(0.5,0.5,0.5);
    glBegin(GL_QUADS);
    {
        glVertex3f( wheelRadius, -(wheelWidth/3), 0);
        glVertex3f( wheelRadius,  (wheelWidth/3), 0);
        glVertex3f(-wheelRadius,  (wheelWidth/3), 0);
        glVertex3f(-wheelRadius, -(wheelWidth/3), 0);
    }
    glEnd();
    glPushMatrix();
    glRotatef(90, 0, 1, 0);
    glBegin(GL_QUADS);
    {
        glVertex3f( wheelRadius, -(wheelWidth/3), 0);
        glVertex3f( wheelRadius,  (wheelWidth/3), 0);
        glVertex3f(-wheelRadius,  (wheelWidth/3), 0);
        glVertex3f(-wheelRadius, -(wheelWidth/3), 0);
    }
    glEnd();
    glPopMatrix();

    glEnd();
}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){
		case 'w':
            //cout<<"BEFORE: "<<centerPos.x<<' '<<centerPos.y<<endl;
			yrotation -= rotationAngleDelta;
            if(yrotation <= -180){
                yrotation =  360 + yrotation;
            }
            centerPos.x = centerPos.x - linDistPerClick*cos(zrotation*pi/180);
            centerPos.y = centerPos.y - linDistPerClick*sin(zrotation*pi/180);
            //cout<<"AFTER: "<<centerPos.x<<' '<<centerPos.y<<endl;
            break;
        case 's':
			yrotation += rotationAngleDelta;
            if(yrotation > 180){
                yrotation = -360 + yrotation;
            }
            centerPos.x = centerPos.x + linDistPerClick*cos(zrotation*pi/180);
            centerPos.y = centerPos.y + linDistPerClick*sin(zrotation*pi/180);
            break;
        case 'a':
            zrotation += rotationAngleDelta;
            if(zrotation > 180){
                zrotation = -360 + zrotation;
            }
            break;
        case 'd':
			zrotation -= rotationAngleDelta;
            if(zrotation <= -180){
                zrotation =  360 + zrotation;
            }
            break;
        

		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			cameraHeight -= 3.0;
			break;
		case GLUT_KEY_UP:		// up arrow key
			cameraHeight += 3.0;
			break;

		case GLUT_KEY_RIGHT:
			cameraAngle += 0.03;
			break;
		case GLUT_KEY_LEFT:
			cameraAngle -= 0.03;
			break;

		case GLUT_KEY_PAGE_UP:
			break;
		case GLUT_KEY_PAGE_DOWN:
			break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
			break;
		case GLUT_KEY_END:
			break;

		default:
			break;
	}
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			//........
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
	gluLookAt(130*cos(cameraAngle), 130*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
	//gluLookAt(0,0,200,	0,0,0,	0,1,0);


	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

	drawGrid();

    //glColor3f(1,0,0);

    glTranslatef(centerPos.x, centerPos.y, centerPos.z);
    glRotatef(zrotation, 0, 0, 1);
    glRotatef(yrotation, 0, 1, 0);
    drawWheel(wheelRadius, wheelWidth, 40);

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
	drawgrid=1;
	cameraHeight=100.0;
	cameraAngle=pi/4;
	angle=0;

    wheelRadius = 25.0;
    wheelWidth = 10.0;
    centerPos = Point(0,0,wheelRadius);
    yrotation = 0.0;
    zrotation = 0.0;
    linDistPerClick = (2*pi*wheelRadius*rotationAngleDelta)/360;

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
