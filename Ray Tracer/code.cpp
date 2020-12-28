#include "bitmap_image.hpp"
#include <bits/stdc++.h>

//#include <windows.h>
#include <GL/glut.h>
#include <GL/glx.h>
#include <GL/gl.h>


using namespace std;

#define pi (2*acos(0.0))
#define theta pi/180 //amount of angle rotation per click
#define MAX_PIXELS 1024
#define inf 1e9
#define eps 1e-6
#define CHECKERBOARD_SIDELEN 30

double cameraHeight;
double cameraAngle;
int drawaxes;
double angle;

double near_dist = 1.0;
double fovX = 90.0, fovY = 90.0;

int level_of_recursion;


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
    struct point operator-(struct point const &p){
		struct point tmp;
		tmp.x = this->x - p.x;
		tmp.y = this->y - p.y;
		tmp.z = this->z - p.z;
		return tmp;
	}
	struct point operator*(double r){
		struct point tmp;
		tmp.x = this->x * r;
		tmp.y = this->y * r;
		tmp.z = this->z * r;
		return tmp;
	}
    void normalize(){
        double val = sqrt(x*x+y*y+z*z);
        x /= val;
        y /= val;
        z /= val;
    }
    void print(){
        cout<<"Point: ("<<this->x<<", "<<this->y<<", "<<this->z<<")"<<endl;
    }
} typedef Point;

typedef point color;

Point typedef Vector;
Point pos; //camera position
Vector uvec, rvec, lvec; // up, right, lool (not left)


void initCamParams(){
	pos = Point(100, 100, 100);
	double x = cos(pi/4);
	uvec = Vector(0, 0, 1);
	rvec = Vector(-x, x, 0);
	lvec = Vector(-x, -x, 0);
}
double dotProduct(Vector v1, Vector v2){
	double tmp;
	tmp = v1.x*v2.x + v1.y*v2.y + v1.z * v2.z;
    return tmp;
}
Vector crossProduct(Vector v1, Vector v2){
	Vector tmp;
	tmp.x = v1.y * v2.z - v1.z * v2.y;
	tmp.y = v1.z * v2.x - v1.x * v2.z;
	tmp.z = v1.x * v2.y - v1.y * v2.x;
	return tmp;
}

double distBetweenPoints(point p1, point p2){
    return sqrt( (p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y) + (p1.z - p2.z)*(p1.z - p2.z) );
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


struct Ray{
    point origin;
    Vector dir;
    Ray(point origin_, Vector dir_){
        origin = origin_;
        dir = dir_;
    }
};


//color finder function 

color getColorValues(point p_hit, Vector norm, color obj_color, Ray ray, double a, double d, double s, double r, double sp_ex, int level);


double sign (point p1, point p2, point p3)
{
    return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
}

struct triangle{
	point pt1, pt2, pt3, clr;

	triangle(point p0, point p1, point p2, color clr_){
		pt1 = p0;
		pt2 = p1;
		pt3 = p2;
        clr = clr_;
	}

	pair<double, color> getIntersectT(Ray ray, double ambient_coeff, double diffuse_coeff, double specular_coeff, double reflection_coeff, double specular_exponent, int level){
		double t = inf;
        color blk(0,0,0);
		Vector side1 = pt2 - pt1;
		Vector side2 = pt3 - pt2;
		Vector side3 = pt1 - pt3;
		// side1.normalize();
		// side2.normalize();
		// side3.normalize();

		Vector norm = crossProduct(side1, side2);
		norm.normalize();


		double dist_of_plane_from_origin = dotProduct(norm, pt1);
		//cout<<dist_of_plane_from_origin<<endl;
		double deno = dotProduct(norm, ray.dir);

		if(abs(deno) <eps){
			return make_pair(inf, blk);
		}
		
		//EDITED THIS FORMULA
		double tmp =  (dotProduct(norm, pt1-ray.origin)) / deno;
		point p_hit = ray.origin + ray.dir * tmp; 
		/*
		Vector op = ray.origin - pt1;
		op.normalize();
		double denominator = dotProduct(norm, ray.dir);
		double numerator = dotProduct(norm, op);

		if(abs(denominator) < eps) {
			return t;
		}
		double t_hit = -numerator/denominator;
		point p_hit = ray.origin + ray.dir * t_hit;
		*/

		if(tmp<0){
			return make_pair(inf, blk);
		}
		//check inside 
		
		Vector v1 = p_hit-pt1;
		Vector v2 = p_hit-pt2;
		Vector v3 = p_hit-pt3;
		if(
			dotProduct(norm, crossProduct(side1, v1)) > (-eps) &&
			dotProduct(norm, crossProduct(side2, v2)) > (-eps) &&
			dotProduct(norm, crossProduct(side3, v3)) > (-eps) 
		){
            point p_intersect = ray.origin + ray.dir * tmp;
            color clr_found = getColorValues(p_intersect, norm, clr, ray, ambient_coeff, diffuse_coeff, specular_coeff, reflection_coeff, specular_exponent, level);
            return make_pair(tmp, clr_found);
		}
		else{
			return make_pair(inf, blk);
		}
		
		

	    /* alternative formula */
		/*
		double d1, d2, d3;
		bool has_neg, has_pos;

		d1 = sign(p_hit, pt1, pt2);
		d2 = sign(p_hit, pt2, pt3);
		d3 = sign(p_hit, pt3, pt1);

		has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
		has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

		if( !(has_neg && has_pos)) return tmp;
		else return inf;
		*/
		

	}

    double getIntersectTonly(Ray ray){
		double t = inf;
		Vector side1 = pt2 - pt1;
		Vector side2 = pt3 - pt2;
		Vector side3 = pt1 - pt3;
		// side1.normalize();
		// side2.normalize();
		// side3.normalize();

		Vector norm = crossProduct(side1, side2);
		norm.normalize();


		double dist_of_plane_from_origin = dotProduct(norm, pt1);
		//cout<<dist_of_plane_from_origin<<endl;
		double deno = dotProduct(norm, ray.dir);

		if(abs(deno) <eps){
			return inf;
		}
		
		//EDITED THIS FORMULA
		double tmp =  (dotProduct(norm, pt1-ray.origin)) / deno;
		point p_hit = ray.origin + ray.dir * tmp; 
		/*
		Vector op = ray.origin - pt1;
		op.normalize();
		double denominator = dotProduct(norm, ray.dir);
		double numerator = dotProduct(norm, op);

		if(abs(denominator) < eps) {
			return t;
		}
		double t_hit = -numerator/denominator;
		point p_hit = ray.origin + ray.dir * t_hit;
		*/

		if(tmp<0){
			return inf;
		}
		//check inside 
		
		Vector v1 = p_hit-pt1;
		Vector v2 = p_hit-pt2;
		Vector v3 = p_hit-pt3;
		if(
			dotProduct(norm, crossProduct(side1, v1)) > (-eps) &&
			dotProduct(norm, crossProduct(side2, v2)) > (-eps) &&
			dotProduct(norm, crossProduct(side3, v3)) > (-eps) 
		){
			return tmp;
		}
		else{
			return inf;
		}
		
		

	    /* alternative formula */
		/*
		double d1, d2, d3;
		bool has_neg, has_pos;

		d1 = sign(p_hit, pt1, pt2);
		d2 = sign(p_hit, pt2, pt3);
		d3 = sign(p_hit, pt3, pt1);

		has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
		has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

		if( !(has_neg && has_pos)) return tmp;
		else return inf;
		*/
		

	}
};


int num_pixels;
int num_objects;
int num_light_src;

struct object{
	int type_;
	//type 1 - sphere
	point center_;
	double radius_;
	//type 2 - pyramid
	point lowest_point_;
	double len_of_base_;
	double height_; 

	color color_;
	double ambient_coeff;
	double diffuse_coeff;
	double specular_coeff;
	double reflection_coeff; 
	double specular_exponent;

	void get_input(int type){
		type_ = type;
		if(type == 1){
			cin>>center_.x>>center_.y>>center_.z;
			cin>>radius_;
		}
		else{
			cin>>lowest_point_.x>>lowest_point_.y>>lowest_point_.z;
			cin>>len_of_base_>>height_;
		}
		cin>>color_.x>>color_.y>>color_.z;
		cin>>ambient_coeff>>diffuse_coeff>>specular_coeff>>reflection_coeff>>specular_exponent;
	}

    pair<double, color> getIntersectionT(Ray ray, int level){
        color blk(0,0,0);
        if(type_ == 1){
            //sphere
			double t = inf;
            point origin = ray.origin; // point
            Vector dir = ray.dir;
            point center_sphere = this->center_;
            double radius = this->radius_;

            Vector o_minus_c = origin - center_sphere;
            double a = 1;
            double b = 2 * dotProduct(dir, o_minus_c);
            double c = dotProduct(o_minus_c, o_minus_c) - radius*radius;

            double discr = b * b - 4 * a * c; 
            //cout<<"here"<<endl;
            if (discr < 0) return make_pair(t, blk); 
			//ADD CODE: Check for eps instead of zero
            //else if (abs(discr) < 0) return  (- 0.5 * b / a); 
            else { 
                double q = (b > 0) ? -0.5 * (b + sqrt(discr)) : -0.5 * (b - sqrt(discr)); 
                double t1 = q / a; 
                double t2 = c / q;
                double t_small = min(t1, t2); 
                double t_big = max(t1, t2);
                
				if(t_small>=0){
                    point p_intersect = ray.origin + ray.dir * t_small;
                    point normal_ = p_intersect - center_;
                    normal_.normalize();
                    color clr_found = getColorValues(p_intersect, normal_, color_, ray, ambient_coeff, diffuse_coeff, specular_coeff, reflection_coeff, specular_exponent, level);
					return make_pair(t_small, clr_found);
				} else if(t_big>=0){
					point p_intersect = ray.origin + ray.dir * t_big;
                    point normal_ = p_intersect - center_;
                    normal_.normalize();
                    color clr_found = getColorValues(p_intersect, normal_, color_, ray, ambient_coeff, diffuse_coeff, specular_coeff, reflection_coeff, specular_exponent, level);
					return make_pair(t_big, clr_found);
				} else{
					return make_pair(inf, blk);
				}
            } 
        }
        else{
			//type 2 - pyramid
			point p1( lowest_point_.x,  lowest_point_.y,  lowest_point_.z);
			point p2( lowest_point_.x+ len_of_base_,  lowest_point_.y,  lowest_point_.z);
			point p3( lowest_point_.x+ len_of_base_,  lowest_point_.y+len_of_base_, lowest_point_.z);
			point p4( lowest_point_.x,  lowest_point_.y+ len_of_base_,  lowest_point_.z);
			point p_top(lowest_point_.x + len_of_base_/2, lowest_point_.y + len_of_base_/2, lowest_point_.z + height_);
			vector<triangle> triangle_array;

			
			triangle_array.push_back(triangle(p_top, p1, p2, color_));
			triangle_array.push_back(triangle(p_top, p2, p3, color_));
		 	triangle_array.push_back(triangle(p_top, p3, p4, color_));
			triangle_array.push_back(triangle(p_top, p4, p1, color_));
			triangle_array.push_back(triangle(p1, p2, p3, color_));
		    triangle_array.push_back(triangle(p3, p4, p1, color_));
			/*p_top.print();
			p3.print();
			p4.print();
			*/
			double t_ = inf;
			int t_idx_ = -1;
            color clr_found(0,0,0);
			int sz = triangle_array.size();
			for(int i = 0; i<sz; i++){
				pair<double, color> tmpNclr = triangle_array[i].getIntersectT(ray, ambient_coeff, diffuse_coeff, specular_coeff, reflection_coeff, specular_exponent, level);
				double tmp = tmpNclr.first;

                if(tmp<t_){
					t_ = tmp;
					t_idx_ = i;
                    clr_found = tmpNclr.second;
				}
			}
			if(t_idx_>=0){
				return make_pair(t_, clr_found);
			}
			else{
				return make_pair(inf, blk);
			}
        }
    }
	

    double getIntersectionTonly(Ray ray){
        if(type_ == 1){
            //sphere
			double t = inf;
            point origin = ray.origin; // point
            Vector dir = ray.dir;
            point center_sphere = this->center_;
            double radius = this->radius_;

            Vector o_minus_c = origin - center_sphere;
            double a = 1;
            double b = 2 * dotProduct(dir, o_minus_c);
            double c = dotProduct(o_minus_c, o_minus_c) - radius*radius;

            double discr = b * b - 4 * a * c; 
            if (discr < 0) return t; 
			//ADD CODE: Check for eps instead of zero
            //else if (abs(discr) < 0) return  (- 0.5 * b / a); 
            else { 
                double q = (b > 0) ? -0.5 * (b + sqrt(discr)) : -0.5 * (b - sqrt(discr)); 
                double t1 = q / a; 
                double t2 = c / q;
                double t_small = min(t1, t2); 
                double t_big = max(t1, t2);
				if(t_small>=0){
					return t_small;
				} else if(t_big>=0){
					return t_big;
				} else{
					return inf;
				}
            } 
        }
        else{
			//type 2 - pyramid
			point p1( lowest_point_.x,  lowest_point_.y,  lowest_point_.z);
			point p2( lowest_point_.x+ len_of_base_,  lowest_point_.y,  lowest_point_.z);
			point p3( lowest_point_.x+ len_of_base_,  lowest_point_.y+len_of_base_,  lowest_point_.z);
			point p4( lowest_point_.x,  lowest_point_.y+ len_of_base_,  lowest_point_.z);
			point p_top(lowest_point_.x + len_of_base_/2, lowest_point_.y + len_of_base_/2, lowest_point_.z + height_);
			vector<triangle> triangle_array;

			
			triangle_array.push_back(triangle(p_top, p1, p2, color_));
		    triangle_array.push_back(triangle(p_top, p2, p3, color_));
			triangle_array.push_back(triangle(p_top, p3, p4, color_));
			triangle_array.push_back(triangle(p_top, p4, p1, color_));
			triangle_array.push_back(triangle(p1, p2, p3, color_));
			triangle_array.push_back(triangle(p3, p4, p1, color_));
			/*p_top.print();
			p3.print();
			p4.print();
			*/
			double t_ = inf;
            int sz = triangle_array.size();
			int t_idx_ = -1;
			for(int i = 0; i<sz; i++){
				double tmp = triangle_array[i].getIntersectTonly(ray);
				if(tmp<t_){
					t_ = tmp;
					t_idx_ = i;
				}
			}
			if(t_idx_>=0){
				return t_;
			}
			else{
				return inf;
			}
        }
    }
};

vector<object> obj_list;
vector<point> light_sources;
color getColorFromRay(Ray ray, int level);

color getColorValues(point p_hit, Vector norm, color obj_color, Ray ray, double a, double d, double s, double r, double sp_ex, int level){
    //normalize
	norm.normalize();
    int ls_cnt = light_sources.size(), obj_cnt = obj_list.size();
    double diffuse=0.0, specular=0.0;
    int effective_src_cnt = 0;


    for(int i = 0; i<ls_cnt; i++){
		Vector direc(light_sources[i] - p_hit);
		direc.normalize();
        Ray back_ray_to_src( p_hit+direc*1.0, direc);
		        
		double t = distBetweenPoints(light_sources[i], p_hit+direc*1.0);
		Vector L = (p_hit - light_sources[i]);
        L.normalize();

		Vector V = ray.origin - p_hit;
		V.normalize();
		Vector R = L - (norm * (2*dotProduct(L, norm)));
		R.normalize();

        int j = 0;
        for(; j<obj_cnt; j++){
            double tmp  = obj_list[j].getIntersectionTonly(back_ray_to_src);
            if(tmp<t){
                break;
            }
        }
        if(j==obj_cnt){
            //no objects in the path
            effective_src_cnt++;
            //norm.normalize();
            diffuse += (d*max(0.0, dotProduct(L*(-1), norm))); //minus L taken as the acute angle is needed
            specular += (s*max(0.0,pow(dotProduct(R, V),sp_ex)));
        }
    }
    color tmpclr(0,0,0);
    if(effective_src_cnt>0){
        //tmpclr = obj_color * (diffuse) + obj_color * (specular);
        tmpclr = obj_color * (diffuse) + color(1,1,1) * (specular);
    }
	Vector reflection_of_eyeRay = ray.dir - norm * (2*dotProduct(ray.dir, norm));
	reflection_of_eyeRay.normalize();
	Ray newray(p_hit + reflection_of_eyeRay*1.0, reflection_of_eyeRay);
	color reflection_color = getColorFromRay(newray, level+1);

    color clr_found = obj_color * a + tmpclr + reflection_color * r;
    return clr_found;
}


point midpoints[MAX_PIXELS][MAX_PIXELS];
color color_array[MAX_PIXELS][MAX_PIXELS];

double checkForCheckerboard(Ray ray){
    point origin = ray.origin;
    Vector dir = ray.dir;
    if(abs(dir.z) < eps ){
        //parallel
        return inf;
    }
    //z must be 0 for the ray
    double t = -origin.z/dir.z;
    if(t>0) return t;
    else return inf;
}

color getColorFromRay(Ray ray, int level){
	if(level == level_of_recursion) {
		color blk(0,0,0);
		return blk;
	}
    color blk(0, 0, 0);
    color wht(1, 1, 1);
    int n_obj = obj_list.size();
    double t = inf;
	int t_idx = -1;
    color clr_found(0,0,0);
    for(int i = 0; i<n_obj; i++){
        pair<double, color> tmpNclr = obj_list[i].getIntersectionT(ray, level);
        //cout<<i<< ' ';
        double tmp = tmpNclr.first;
        if(tmp<t){
            t = tmp;
            t_idx = i;
            clr_found = tmpNclr.second;
        }
    }
    //cout<<endl;

    double tmp = checkForCheckerboard(ray);
    if(tmp<t){
		t_idx = num_objects; // no object at obj_list[num_objects]
        //intersects the CB
        //find x and y
        point intersection_pt = ray.origin + ray.dir * tmp;
        double x_cb = intersection_pt.x;
        double y_cb = intersection_pt.y;
        int xi = abs(x_cb/CHECKERBOARD_SIDELEN);
        int yi = abs(y_cb/CHECKERBOARD_SIDELEN);
		color obj_clr;
        if(x_cb*y_cb >= 0){
            if((xi+yi)%2==0){
                obj_clr = wht;
            }
            else{
                obj_clr = blk;
            }
        }
        else{
            if((xi+yi)%2==0){
                obj_clr = blk;
            }
            else{
                obj_clr = wht;
            }
        }
		Vector norm(0,0,1);
		clr_found = getColorValues(intersection_pt, norm, obj_clr, ray, 0.2, 0.3, 0.2, 0.3, 5, level);
    }
    
    if(t_idx>=0){    
        return clr_found;
    }
    else{
        return blk;
    }
}

void generateImage(){
    cout<<"generating image..."<<endl;
    double img_width = 2 * near_dist * tan((fovX/2)*(pi/180));
    double img_height = 2 * near_dist * tan((fovY/2)*(pi/180));

    double pixel_width = img_width/num_pixels;
    double pixel_height = img_height/num_pixels;

    point central_point = pos + (lvec*near_dist);
    point top_left_midpoint = central_point + (rvec*(-img_width/2)) + (uvec*(img_height/2)) + (rvec*(pixel_width/2)) + (uvec*(-pixel_height/2));
    //generating midpoint array now
    for(int i = 0; i<num_pixels; i++){
        for(int j = 0; j<num_pixels; j++){
            midpoints[i][j] = top_left_midpoint + (uvec * (i * (-pixel_height))) + (rvec * (j * pixel_width)); 
        }
    }
	/*
	midpoints[0][0].print();
	midpoints[0][767].print();
	midpoints[767][0].print();
	midpoints[767][767].print();
	*/

    for(int i = 0; i<num_pixels; i++){
        for(int j = 0; j<num_pixels; j++){
            //midpoints[i][j] = top_left_midpoint + (uvec * (i * (-pixel_height))) + (rvec * (j * pixel_width)); 
            Vector direc = midpoints[i][j]-pos;
            direc.normalize();
            //Ray ray(pos, direc);
            Ray ray(midpoints[i][j], direc);
            color_array[i][j] = getColorFromRay(ray, 0);
            //cout<<i<<' ' <<j<<endl;
            //color_array[i][j].print();
        }
    }



    
    bitmap_image image(num_pixels,num_pixels);
    for (int x = 0; x < num_pixels; x++) {
        for (int y = 0; y < num_pixels; y++) {
            image.set_pixel(y, x, min(255.0, color_array[x][y].x*255), min(255.0, color_array[x][y].y*255), min(255.0, color_array[x][y].z*255));
        }
    }


    image.save_image("out.bmp");


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
	}glEnd();
}



void drawSphere(double radius,int slices,int stacks, color color_)
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
			points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
    glColor3f(color_.x, color_.y, color_.z);
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
                glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
			}glEnd();
		}
	}
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

        case '0':
            generateImage();

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

void drawCheckerBoard(){
	float sidelen = CHECKERBOARD_SIDELEN;
	float how_far = 1000;
	float offset = sidelen * how_far;
	offset /= 2; 
	glPushMatrix();
	{
		for(int i = 0; i<how_far; i++){
			for(int j = 0; j<how_far; j++){
				if((i+j)%2){
				    glColor3f(0,0,0);
				}
				else{
                    glColor3f(1,1,1);
				}
				glBegin(GL_QUADS);{
                    glVertex3f( - offset + sidelen*i, - offset + sidelen*j,0);
                    glVertex3f( - offset + sidelen*(i+1), - offset + sidelen*j,0);
                    glVertex3f( - offset + sidelen*(i+1), - offset + sidelen*(j+1),0);
                    glVertex3f( - offset + sidelen*i, - offset + sidelen*(j+1),0);
                }glEnd();
				glColor3f(1,1,1);
			}
		}
	}
	glPopMatrix();
}


void drawObjects(){
	glPushMatrix();
	{
		int sz = obj_list.size();
		for(int i = 0; i<sz; i++){
            glPushMatrix();
            {
                object obj = obj_list[i];
                if(obj.type_==1){
                    //sphere
                    glTranslatef(obj.center_.x, obj.center_.y, obj.center_.z);
                    drawSphere(obj.radius_, 30, 30, obj.color_);
                }
                else{
                    //pyramid
                    glColor3f(obj.color_.x, obj.color_.y, obj.color_.z);
                    glBegin(GL_QUADS);{
                        glVertex3f( obj.lowest_point_.x, obj.lowest_point_.y, obj.lowest_point_.z);
                        glVertex3f( obj.lowest_point_.x+obj.len_of_base_, obj.lowest_point_.y, obj.lowest_point_.z);
                        glVertex3f( obj.lowest_point_.x+obj.len_of_base_, obj.lowest_point_.y+obj.len_of_base_, obj.lowest_point_.z);
                        glVertex3f( obj.lowest_point_.x, obj.lowest_point_.y+obj.len_of_base_, obj.lowest_point_.z);
                    }glEnd();
                    double a, b, c;
                    a = obj.lowest_point_.x + obj.len_of_base_/2;
                    b = obj.lowest_point_.y + obj.len_of_base_/2;
                    c = obj.lowest_point_.z + obj.height_;
                    glBegin(GL_TRIANGLES);{
                        glVertex3f( obj.lowest_point_.x, obj.lowest_point_.y, obj.lowest_point_.z);
                        glVertex3f( obj.lowest_point_.x+obj.len_of_base_, obj.lowest_point_.y, obj.lowest_point_.z);
                        glVertex3f( a, b, c);
                    }glEnd();

                    glBegin(GL_TRIANGLES);{
                        glVertex3f( obj.lowest_point_.x+obj.len_of_base_, obj.lowest_point_.y, obj.lowest_point_.z);
                        glVertex3f( obj.lowest_point_.x+obj.len_of_base_, obj.lowest_point_.y+obj.len_of_base_, obj.lowest_point_.z);
                        glVertex3f( a, b, c);
                    }glEnd();

                    glBegin(GL_TRIANGLES);{
                        glVertex3f( a, b, c);
                        glVertex3f( obj.lowest_point_.x+obj.len_of_base_, obj.lowest_point_.y+obj.len_of_base_, obj.lowest_point_.z);
                        glVertex3f( obj.lowest_point_.x, obj.lowest_point_.y+obj.len_of_base_, obj.lowest_point_.z);
                    }glEnd();

                    glBegin(GL_TRIANGLES);{
                        glVertex3f( a, b, c);
                        glVertex3f( obj.lowest_point_.x, obj.lowest_point_.y, obj.lowest_point_.z);
                        glVertex3f( obj.lowest_point_.x, obj.lowest_point_.y+obj.len_of_base_, obj.lowest_point_.z);
                    }glEnd();
                    
                }
            }
            glPopMatrix();
		}

        int sz1 = light_sources.size();
		for(int i = 0; i<sz1; i++){
            glPushMatrix();
            {
                point obj = light_sources[i];
                //sphere
                glTranslatef(obj.x, obj.y, obj.z);
                drawSphere(2, 30, 30, color(0.5,0.0,0.5));
            }
            glPopMatrix();
        }
	}
	glPopMatrix();
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
    drawObjects();
	//drawAxes();
	drawCheckerBoard();
   

    //glColor3f(1,0,0);
    //drawSquare(10);

    //drawSS();


	




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
	gluPerspective(90,	1,	1,	1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}

int main(int argc, char **argv){

    //Input Parsing 
	freopen("description.txt", "r", stdin);
	cin>>level_of_recursion>>num_pixels;
	//cout<<"GOT: "<<level_of_recursion<< "  "<< num_pixels<<endl;

	cin>>num_objects;
	int obj_cnt = 0;
	while(obj_cnt<num_objects){
		
		string typeOfObj;
		//getline(cin, typeOfObj);
		cin>>typeOfObj;
		if(typeOfObj.length()>0){
			obj_cnt++;
			object obj;
			if(typeOfObj == "sphere"){
				obj.get_input(1);
			}
			else{
				obj.get_input(2);
			}
			obj_list.push_back(obj);
		}
	}
	
	int light_sources_cnt_;
	cin>>light_sources_cnt_;

	for(int i = 0; i<light_sources_cnt_; i++){
		point p_;
		cin>>p_.x>>p_.y>>p_.z;
		light_sources.push_back(p_);
	}



	initCamParams();
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