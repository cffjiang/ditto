#include <GL/glut.h> 
#include "ALGEBRA.h"
#include "GEOMETRY.h"
#include "DEFORMABLE_OBJECTS.h"
#include "FIXED_COROTATED.h"
#include <fstream>
#include <iostream>

using namespace ALGEBRA;
using namespace GEOMETRY;
using namespace std;

#define PI 3.14159265

int  winWidth, winHeight;
double costheta,sintheta;
MATRIX_3X3<double> ACCUMULATED=ACCUMULATED.Identity();
VECTOR_3D<double> eye(0.0, 0.0, 5.0);   //original eye position
VECTOR_3D<double> ori(0.0, 1.0, 0.0);   //original orientation of camera
VECTOR_3D<double> u(0.0, 0.0, 5.0);		//current position
VECTOR_3D<double> v1(0.0, 0.0, 5.0);	//last position
double one_over_sqrt_3=(double)1/sqrt((double)3);
double rotation_scale=(double)10;

bool trackingMouse=false;
bool trackingMove=false;

//Geometry to render
TETRAHEDRON_MESH* tet_mesh=0;
VECTOR<double>* particles;

//Simulation data
BACKWARD_EULER_TIME_STEPPING_3D<double>* backward_euler;
VECTOR<int>* cons_left;
VECTOR<int>* cons_right;
VECTOR<int>* cons;
VECTOR<double>* cons_pos;
bool constrain_both_sides=false;
double boundary_condition_speed=(double)0;

void processNormalKeys(unsigned char key, int x,int y){
	if(key==27)
		exit(0);
	if(key=='z'||key=='-')
		eye = eye*1.01	;
	if(key=='Z'||key=='+')
		eye = eye*0.99;
}

void processSpecialKeys(int key, int x, int y){
	if(key==GLUT_KEY_UP) eye = eye*0.99;
	if(key==GLUT_KEY_DOWN) eye = eye*1.01;
}

void renderTetMesh(void) {
	glColor3f(0,0,0);
	
	if(tet_mesh){//draw a tet mesh
		if(tet_mesh->boundary_triangle_mesh){
			for(int tri=0;tri<tet_mesh->boundary_triangle_mesh->Size();tri++){
				glBegin(GL_LINE_LOOP);
				int i=(*(tet_mesh->boundary_triangle_mesh))(tri)(0);
				int j=(*(tet_mesh->boundary_triangle_mesh))(tri)(1);
				int k=(*(tet_mesh->boundary_triangle_mesh))(tri)(2);
				glVertex3f((*particles)(3*i),(*particles)(3*i+1),(*particles)(3*i+2));
				glVertex3f((*particles)(3*j),(*particles)(3*j+1),(*particles)(3*j+2));
				glVertex3f((*particles)(3*k),(*particles)(3*k+1),(*particles)(3*k+2));
				glEnd();}}}
	//glutWireCube(2.0);
}

void track_eye(VECTOR_3D<double> u, VECTOR_3D<double> v1)
{
    VECTOR_3D<double> v2,v3;
	v1=ACCUMULATED*v1;
	v1.Normalize();
	u=ACCUMULATED*u;
	u.Normalize();
	v3=VECTOR_3D<double>::Cross_Product(v1,u);
	v3.Normalize();
	v2=u-v1.Dot_Product(v1,u)*v1;
	v2.Normalize();
	
	costheta = v1.Dot_Product(v1, u);   //use dot product, v1 and u have already been normalized
	sintheta = sqrt(1-costheta*costheta);	
	
	MATRIX_3X3<double> Q(v1,v2,v3);
	MATRIX_3X3<double> R(costheta, -sintheta,0,sintheta,costheta,0,0,0,1);
	eye=Q*R*Q.Transposed()*eye;
	ori=Q*R*Q.Transposed()*ori;
	ACCUMULATED=Q*R*Q.Transposed()*ACCUMULATED;
}

void startMotion(int x, int y)
{
    trackingMouse = true;
	trackingMove = true;
	
	v1(0)=rotation_scale*(((double)2*(double)x-(double)winWidth)/(double)winWidth)*one_over_sqrt_3;
	v1(1)=rotation_scale*(((double)winHeight-(double)2*(double)y)/(double)winHeight)*one_over_sqrt_3;
}

void stopMotion(int x, int y)
{
    trackingMouse=false;
	trackingMove=false;
	v1(0)=u(0);
	v1(1)=u(1);
}

void mouseMotion(int x, int y)		
{
	u(0)=rotation_scale*(((double)2*(double)x-(double)winWidth)/(double)winWidth)*one_over_sqrt_3;
	u(1)=rotation_scale*(((double)winHeight-(double)2*(double)y)/(double)winHeight)*one_over_sqrt_3;
	
    if(trackingMouse){
		if(fabs(v1(0)-u(0)) + fabs(v1(1) - u(1)) >0.1) {
			track_eye(u, v1);
			v1(0)= u(0);
			v1(1)= u(1);}} 
    glutPostRedisplay();
	
}

void mouseButton(int button, int state, int x, int y)
{
    if(button!=GLUT_LEFT_BUTTON) return;
    switch(state) 
    {
		case GLUT_DOWN:
			startMotion(x,y);
			break;
		case GLUT_UP:
			stopMotion(x,y);
			break;
    } 
}

void Update_Constrained_Nodes(float time,VECTOR<int>& constrained_left,VECTOR<int>& constrained_right,VECTOR<double>& constrained_locations,const bool left_and_right=true,double speed=(double).25){
	
	float z_right=1;
	z_right+=.25*sin(2*PI*time*speed);
	
	float z_left=0;
	z_left-=.25*sin(2*PI*time*speed);
	
	for(int p=0;p<constrained_right.Size()/3;p++){
		constrained_locations(constrained_left.Size()+3*p+2)=z_right;}
	
	if(left_and_right){
		for(int p=0;p<constrained_left.Size()/3;p++){
			constrained_locations(3*p+2)=z_left;}}
	
}

void Advance_One_Time_Step(){
	Update_Constrained_Nodes(backward_euler->Time(),*cons_left,*cons_right,*cons_pos,constrain_both_sides,boundary_condition_speed);
	backward_euler->Set_Boundary_Conditions(*cons,*cons_pos);
	backward_euler->Advance_One_Time_Step();
}

void display (void) {
	
	glClearColor(1.0f, .84f, .0f, 1.0f); 
	glClear(GL_COLOR_BUFFER_BIT); 
	glLoadIdentity();
	
	gluLookAt(eye(0),eye(1),eye(2), 0.0,0.0,0.0, ori(0),ori(1),ori(2));
	Advance_One_Time_Step();
	renderTetMesh();
	glFlush();
	
}

void reshape (int width, int height) {
	glViewport(0, 0, (GLsizei)width, (GLsizei)height);
	glMatrixMode(GL_PROJECTION); 
	glLoadIdentity(); 
	gluPerspective(60, (GLfloat)width / (GLfloat)height, 1.0, 100.0);
	glMatrixMode(GL_MODELVIEW); 
}

VECTOR_2D<int> Count_Constrained_Nodes(VECTOR<double>& positions,const bool left_and_right=true){
	
	float max_z=-FLT_MAX;
	float min_z=FLT_MAX;
	
	for(int p=0;p<positions.Size()/3;p++){
		if(positions(3*p+2)>max_z) max_z=positions(3*p+2);
		if(positions(3*p+2)<min_z) min_z=positions(3*p+2);}
	
	float epsilon=1e-4;
	
	int num_max=0;
	int num_min=0;
	
	for(int p=0;p<positions.Size()/3;p++){
		if(max_z-positions(3*p+2)<epsilon) num_max++;
		if(positions(3*p+2)-min_z<epsilon) num_min++;}
	
	int num_right;
	if(left_and_right) num_right=3*num_max;
	else num_right=0;
	
	return VECTOR_2D<int>(3*num_min,num_right);
}

void Initialize_Constrained_Nodes(VECTOR<int>& constrained_left,VECTOR<int>& constrained_right,VECTOR<double>& positions,const bool left_and_right=true){
	float max_z=-FLT_MAX;
	float min_z=FLT_MAX;
	
	for(int p=0;p<positions.Size()/3;p++){
		if(positions(3*p+2)>max_z) max_z=positions(3*p+2);
		if(positions(3*p+2)<min_z) min_z=positions(3*p+2);}
	
	float epsilon=1e-4;
	
	int num_max=0;
	int num_min=0;
	
	for(int p=0;p<positions.Size()/3;p++){
		if(max_z-positions(3*p+2)<epsilon && left_and_right){constrained_right(num_max++)=3*p;constrained_right(num_max++)=3*p+1;constrained_right(num_max++)=3*p+2;}
		if(positions(3*p+2)-min_z<epsilon){constrained_left(num_min++)=3*p;constrained_left(num_min++)=3*p+1;constrained_left(num_min++)=3*p+2;}}
	
	assert(num_max==constrained_right.Size());
	assert(num_min==constrained_left.Size());
}

int main(int argc, char **argv)
{
	//set up geometry
	int n=5;//the number of particles per side,n>=2	
	particles=new VECTOR<double>(3*n*n*n);	
	//x_p = (i*d_x, j*d_x, k*d_x)
	//p = i + n*j + n*n*k
	
	//number the particles
	double d_x=(double)1.0/(double)(n-1);	
	GRID_3D<double> grid(n,d_x,(double)0,(double)0,(double)0);
	DEFORMABLE_OBJECT_3D<double> deformable_object(grid);
	tet_mesh=&(deformable_object.Tetrahedron_Mesh());
	tet_mesh->Initialize_Oriented_Boundary_Triangles();
	particles=&(deformable_object.Positions());
	for(int k=0;k<n;k++){
		for(int j=0;j<n;j++){
			for(int i=0;i<n;i++){
				(*particles)(3*(i+n*j+n*n*k))=i*d_x;
				(*particles)(3*(i+n*j+n*n*k)+1)=j*d_x;
				(*particles)(3*(i+n*j+n*n*k)+2)=k*d_x;}}}
	
	//set up physics
	//elastic constitutive model
	//LINEAR_ELASTICITY_3D<double> le((float)10000,(float).3);
	FIXED_COROTATED_3D<double> le((float)100000,(float).3);

	FEM_HYPERELASTICITY_3D<double> fem(deformable_object.Tetrahedron_Mesh(),deformable_object.Positions());
	fem.Set_Constitutive_Model(le);
	fem.Initialize_Undeformed_Configuration();
	
	//Choose time stepping: backward Euler
	double dt=(double)1/(double)300;
	int number_steps=30;
	float start_time=(float)0;
	float end_time=(float)number_steps*dt;
	
	BACKWARD_EULER_TIME_STEPPING_3D<double> be(dt,end_time,start_time,deformable_object);
	backward_euler=&be;
	be.Set_Elastic_Forces(fem);
	VECTOR<double> nodal_volumes(deformable_object.Positions().Size());
	fem.Nodal_Volume_Fractions(nodal_volumes);
	
	be.Initialize_BE_Matrix(nodal_volumes);
	be.Initialize_CG();
	
	//boundary conditions
	VECTOR_2D<int> numbers_constrained=Count_Constrained_Nodes(*particles,constrain_both_sides);
	VECTOR<int> constrained_left(numbers_constrained.x());
	cons_left=&constrained_left;
	VECTOR<int> constrained_right(numbers_constrained.y());
	cons_right=&constrained_right;
	Initialize_Constrained_Nodes(constrained_left,constrained_right,*particles,constrain_both_sides);	
	
	int number_constrained=constrained_left.Size()+constrained_right.Size();
	VECTOR<int> constrained(number_constrained);
	cons=&constrained;
	for(int i=0;i<constrained_left.Size();i++) constrained(i)=constrained_left(i);
	for(int i=0;i<constrained_right.Size();i++) constrained(i+constrained_left.Size())=constrained_right(i);
	VECTOR<double> constrained_locations(number_constrained);
	for(int i=0;i<number_constrained;i++) constrained_locations(i)=(*particles)(constrained(i));
	cons_pos=&constrained_locations;
		
	glutInit(&argc, argv); 
	glutInitDisplayMode(GLUT_SINGLE); 
	glutInitWindowSize(500,500); 
	glutInitWindowPosition(100,100);
	glutCreateWindow("Simulation Visualization"); 
	
	winWidth=glutGet( GLUT_WINDOW_WIDTH );
	winHeight=glutGet( GLUT_WINDOW_HEIGHT);
	
	glutDisplayFunc(display); 
	glutIdleFunc(display); 
	glutMouseFunc(mouseButton);
	glutMotionFunc(mouseMotion);
	glutKeyboardFunc(processNormalKeys);
	glutSpecialFunc(processSpecialKeys);
	glutReshapeFunc(reshape); 	
	glutMainLoop(); 
	
	//if(tet_mesh) delete tet_mesh;
	if(particles) delete particles;
}

